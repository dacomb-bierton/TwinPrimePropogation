#!/usr/bin/env python3

import os
import sys
import time
import math
import pickle
import csv
import numpy as np
from pathlib import Path

OUTPUT_DIR       = Path("true_twin_propagation")
CHECKPOINT_EVERY = 40_000
SEGMENT_SIZE     = 150_000_000
MAX_SUCCESS_ROWS = 4_000_000

RESUME_FROM = 274_877_906_944

OUTPUT_DIR.mkdir(exist_ok=True)

def path(name: str) -> Path:
    return OUTPUT_DIR / f"true_{name}"

try:
    import gmpy2
    def is_prime(n):
        return gmpy2.is_prime(int(n))
    print("Using gmpy2 for primality (fast)")
except ImportError:
    def is_prime(n):
        n = int(n)
        if n < 2:
            return False
        if n < 4:
            return True
        if n % 2 == 0 or n % 3 == 0:
            return False

        d = n - 1
        s = 0
        while d % 2 == 0:
            d //= 2
            s += 1

        for a in (2, 3, 5, 7, 11, 13, 23):
            if a >= n:
                continue
            x = pow(a, d, n)
            if x == 1 or x == n - 1:
                continue
            for _ in range(s - 1):
                x = pow(x, 2, n)
                if x == n - 1:
                    break
            else:
                return False
        return True
    print("Using fast Miller-Rabin for primality")

def is_twin_member(x):
    x = int(x)
    if not is_prime(x):
        return False
    return is_prime(x - 2) or is_prime(x + 2)

def generate_small_primes(limit: int) -> np.ndarray:
    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[0:2] = False
    for i in range(2, int(math.isqrt(limit)) + 1):
        if is_prime[i]:
            is_prime[i*i::i] = False
    return np.nonzero(is_prime)[0]

def sieve_segment(low: int, high: int, small_primes: np.ndarray) -> np.ndarray:
    size = high - low
    is_prime = np.ones(size, dtype=bool)
    if low == 0:
        is_prime[0:2] = False
    for p in small_primes:
        start = max(p * p, ((low + p - 1) // p) * p)
        is_prime[start - low::p] = False
    return is_prime

def save_checkpoint(last_high, last_twin, dyadic, success_count):
    data = {
        "last_high": last_high,
        "last_twin": last_twin,
        "dyadic": dyadic,
        "success_count": success_count,
        "timestamp": time.time()
    }
    tmp = path("checkpoint.tmp")
    final = path("checkpoint.pkl")
    with open(tmp, "wb") as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.flush()
        os.fsync(f.fileno())
    tmp.replace(final)
    print(f"  [checkpoint] {last_high:,}  |  last twin {last_twin:,}")

def load_checkpoint():
    f = path("checkpoint.pkl")
    if not f.exists():
        return None
    try:
        with open(f, "rb") as fh:
            return pickle.load(fh)
    except Exception as e:
        print(f"Checkpoint bad ({e}) – falling back")
        return None

def update_dyadic(dyadic, p, success):
    if p < 16:
        return
    k = p.bit_length() - 1
    key = (1 << k, 1 << (k + 1))
    if key not in dyadic:
        dyadic[key] = {"pairs": 0, "success": 0}
    dyadic[key]["pairs"] += 1
    if success:
        dyadic[key]["success"] += 1

def write_dyadic_csv(dyadic):
    rows = []
    for (lo, hi), v in sorted(dyadic.items()):
        rate = v["success"] / v["pairs"] if v["pairs"] else 0.0
        rows.append([lo, hi, v["pairs"], v["success"], f"{rate:.10f}"])
    with open(path("propagation_dyadic_rates.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["lo", "hi", "pairs", "success", "rate"])
        writer.writerows(rows)

def main():
    print("=" * 60)
    print("True Twin-Propagation – Faster Reliable Version")
    print("=" * 60)

    ckpt = load_checkpoint()

    if ckpt:
        last_high     = ckpt["last_high"]
        last_twin     = ckpt.get("last_twin", 0)
        dyadic        = {tuple(k) if isinstance(k, (list, tuple)) else k: v
                         for k, v in ckpt["dyadic"].items()}
        success_count = ckpt.get("success_count", 0)
        print(f"Resumed from checkpoint → {last_high:,}")
    else:
        last_high     = RESUME_FROM
        last_twin     = 0
        dyadic        = {}
        success_count = 0
        print(f"Starting from {RESUME_FROM:,}")

        rates_file = path("propagation_dyadic_rates.csv")
        if rates_file.exists():
            print("Loading dyadic history...")
            with open(rates_file, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    lo = int(row["lo"])
                    if lo >= 274877906944:
                        continue
                    dyadic[(lo, int(row["hi"]))] = {
                        "pairs": int(row["pairs"]),
                        "success": int(row["success"])
                    }
            print(f"Loaded {len(dyadic)} intervals")

    small_primes = generate_small_primes(3_000_000)

    success_file = path("successful_propagations.csv")
    write_header = not success_file.exists()
    success_f = open(success_file, "a", buffering=1024*1024)
    if write_header:
        success_f.write("p_n,p_n+1,C,D,C_is_twin,D_is_twin\n")

    current_low = last_high
    pairs_since_ckpt = 0
    prev_twin = last_twin

    try:
        while True:
            current_high = current_low + SEGMENT_SIZE
            print(f"\nSegment [{current_low:,} → {current_high:,}]")

            is_prime = sieve_segment(current_low, current_high, small_primes)

            new_twins = []
            for i in range(len(is_prime) - 2):
                if is_prime[i] and is_prime[i+2]:
                    p = current_low + i
                    if p >= 3:
                        new_twins.append(p)

            if not new_twins:
                current_low = current_high
                continue

            pairs = []
            if prev_twin > 0:
                pairs.append((prev_twin, new_twins[0]))
            for i in range(len(new_twins)-1):
                pairs.append((new_twins[i], new_twins[i+1]))

            for p, q in pairs:
                C = p + q + 1
                D = p + q + 3
                c_twin = is_twin_member(C)
                d_twin = is_twin_member(D)
                succ = c_twin or d_twin

                update_dyadic(dyadic, p, succ)

                if succ:
                    success_count += 1
                    if success_count <= MAX_SUCCESS_ROWS:
                        success_f.write(f"{p},{q},{C},{D},{c_twin},{d_twin}\n")

            prev_twin = new_twins[-1]
            current_low = current_high
            pairs_since_ckpt += len(new_twins)

            if pairs_since_ckpt >= CHECKPOINT_EVERY:
                success_f.flush()
                write_dyadic_csv(dyadic)
                save_checkpoint(current_low, prev_twin, dyadic, success_count)
                pairs_since_ckpt = 0

    except KeyboardInterrupt:
        print("\nInterrupted – saving...")
    finally:
        success_f.flush()
        success_f.close()
        write_dyadic_csv(dyadic)
        save_checkpoint(current_low, prev_twin, dyadic, success_count)
        print("Shutdown complete.")
        print(f"Last high: {current_low:,}")
        print(f"Successes: {success_count:,}")

if __name__ == "__main__":
    main()
