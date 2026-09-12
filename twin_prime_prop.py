#!/usr/bin/env python3
"""Computational search for the Twin-Prime Propagation Conjecture.

For consecutive lower twin primes p < q, let C = p + q + 1 and D = p + q + 3.
The pair (p, q) is recorded as *propagating* when C or D is a member of a
twin-prime pair.  The script sieves successive segments, tests every
consecutive pair, records successes and per-dyadic-interval success rates,
and checkpoints so that a long run can be interrupted and resumed.
"""

import argparse
import csv
import math
import os
import pickle
import signal
import sys
import time
from pathlib import Path

import numpy as np

OUTPUT_DIR       = Path("true_twin_propagation")
CHECKPOINT_EVERY = 40_000
SEGMENT_SIZE     = 150_000_000
MAX_SUCCESS_ROWS = 4_000_000
RESUME_FROM      = 274_877_906_944  # 2**38

try:
    import gmpy2
except ImportError:
    gmpy2 = None

# Miller-Rabin with the first twelve primes as witnesses is deterministic for
# every n < 318 665 857 834 031 151 167 461 (~3.2e23), far beyond 64 bits.
_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def _is_prime_miller_rabin(n) -> bool:
    n = int(n)
    if n < 2:
        return False
    for p in _MR_BASES:
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


if gmpy2 is not None:
    PRIMALITY_BACKEND = "gmpy2"

    def is_prime(n) -> bool:
        return bool(gmpy2.is_prime(int(n)))
else:
    PRIMALITY_BACKEND = "deterministic Miller-Rabin"
    is_prime = _is_prime_miller_rabin


def is_twin_member(x) -> bool:
    """True if x is prime and x-2 or x+2 is prime, i.e. x belongs to a twin pair."""
    x = int(x)
    return is_prime(x) and (is_prime(x - 2) or is_prime(x + 2))


def propagation_status(p: int, q: int):
    """Return (C, D, C_is_twin, D_is_twin) for a consecutive lower-twin pair."""
    C = p + q + 1
    D = C + 2
    c_prime = is_prime(C)
    d_prime = is_prime(D)
    c_twin = c_prime and (d_prime or is_prime(C - 2))
    d_twin = d_prime and (c_prime or is_prime(D + 2))
    return C, D, c_twin, d_twin


def generate_small_primes(limit: int) -> np.ndarray:
    flags = np.ones(limit + 1, dtype=bool)
    flags[:2] = False
    for i in range(2, math.isqrt(limit) + 1):
        if flags[i]:
            flags[i * i::i] = False
    return np.nonzero(flags)[0]


def sieve_odd_segment(low: int, high: int, small_primes: np.ndarray):
    """Primality flags for the odd numbers in [low, high).

    Returns (first_odd, flags) where flags[j] refers to n = first_odd + 2*j.
    small_primes must contain every prime <= sqrt(high).
    """
    first_odd = low | 1
    size = max(0, (high - first_odd + 1) // 2)
    flags = np.ones(size, dtype=bool)
    if first_odd == 1 and size:
        flags[0] = False
    for p in small_primes.tolist():
        if p * p >= high:
            break
        if p == 2:
            continue
        start = max(p * p, -(-low // p) * p)
        if start % 2 == 0:
            start += p
        if start < high:
            flags[(start - first_odd) // 2::p] = False
    return first_odd, flags


def lower_twins_in_segment(low: int, high: int, small_primes: np.ndarray) -> np.ndarray:
    """Every lower twin prime p with low <= p < high, in increasing order.

    The sieve is extended by two so that a pair straddling the upper boundary
    (p = high-1, p+2 = high+1) is not lost.
    """
    first_odd, flags = sieve_odd_segment(low, high + 2, small_primes)
    n_in_range = max(0, (high - first_odd + 1) // 2)
    twins = 2 * np.flatnonzero(flags[:n_in_range] & flags[1:n_in_range + 1]) + first_odd
    return twins[twins >= 3]


class SmallPrimeCache:
    """Small-prime table that grows as the search frontier advances."""

    def __init__(self):
        self.limit = 0
        self.primes = np.empty(0, dtype=np.int64)

    def covering(self, high: int) -> np.ndarray:
        need = math.isqrt(high) + 1
        if need > self.limit:
            self.limit = max(need + need // 4, 3_000_000)
            self.primes = generate_small_primes(self.limit)
        return self.primes


def save_checkpoint(out, last_high, last_twin, dyadic, success_count):
    data = {
        "last_high": last_high,
        "last_twin": last_twin,
        "dyadic": dyadic,
        "success_count": success_count,
        "timestamp": time.time(),
    }
    tmp = out("checkpoint.tmp")
    final = out("checkpoint.pkl")
    with open(tmp, "wb") as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.flush()
        os.fsync(f.fileno())
    tmp.replace(final)
    print(f"  [checkpoint] {last_high:,}  |  last twin {last_twin:,}")


def load_checkpoint(out):
    f = out("checkpoint.pkl")
    if not f.exists():
        return None
    try:
        with open(f, "rb") as fh:
            return pickle.load(fh)
    except Exception as e:
        print(f"Checkpoint unreadable ({e}) - ignoring it")
        return None


def update_dyadic(dyadic, p, success):
    if p < 16:
        return
    k = p.bit_length() - 1
    key = (1 << k, 1 << (k + 1))
    entry = dyadic.setdefault(key, {"pairs": 0, "success": 0})
    entry["pairs"] += 1
    if success:
        entry["success"] += 1


def write_dyadic_csv(out, dyadic):
    final = out("propagation_dyadic_rates.csv")
    tmp = out("propagation_dyadic_rates.tmp")
    with open(tmp, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["lo", "hi", "pairs", "success", "rate"])
        for (lo, hi), v in sorted(dyadic.items()):
            rate = v["success"] / v["pairs"] if v["pairs"] else 0.0
            writer.writerow([lo, hi, v["pairs"], v["success"], f"{rate:.10f}"])
        f.flush()
        os.fsync(f.fileno())
    tmp.replace(final)


def load_dyadic_history(out, start):
    """Intervals from an earlier run that lie below the starting point."""
    rates_file = out("propagation_dyadic_rates.csv")
    dyadic = {}
    if not rates_file.exists():
        return dyadic
    with open(rates_file, "r", newline="") as f:
        for row in csv.DictReader(f):
            lo = int(row["lo"])
            if lo >= start:
                continue
            dyadic[(lo, int(row["hi"]))] = {
                "pairs": int(row["pairs"]),
                "success": int(row["success"]),
            }
    print(f"Loaded {len(dyadic)} dyadic intervals from previous results")
    return dyadic


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="Search for propagating pairs of consecutive lower twin primes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--start", type=int, default=RESUME_FROM,
                    help="lower bound of the search when no checkpoint exists")
    ap.add_argument("--stop-at", type=int, default=None,
                    help="stop once the search frontier reaches this value (default: run until interrupted)")
    ap.add_argument("--segment-size", type=int, default=SEGMENT_SIZE,
                    help="numbers sieved per segment")
    ap.add_argument("--checkpoint-every", type=int, default=CHECKPOINT_EVERY,
                    help="checkpoint after this many new twin primes")
    ap.add_argument("--max-success-rows", type=int, default=MAX_SUCCESS_ROWS,
                    help="cap on rows written to the successes CSV")
    ap.add_argument("--output-dir", type=Path, default=OUTPUT_DIR,
                    help="directory for CSV output and checkpoints")
    args = ap.parse_args(argv)
    if args.segment_size < 4:
        ap.error("--segment-size must be at least 4")
    if args.start < 0:
        ap.error("--start must be non-negative")
    return args


def main(argv=None):
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    def out(name: str) -> Path:
        return args.output_dir / f"true_{name}"

    print("=" * 60)
    print("True Twin-Propagation")
    print("=" * 60)
    print(f"Primality backend : {PRIMALITY_BACKEND}")
    print(f"Output directory  : {args.output_dir.resolve()}")

    ckpt = load_checkpoint(out)
    if ckpt:
        last_high     = ckpt["last_high"]
        last_twin     = ckpt.get("last_twin", 0)
        dyadic        = {tuple(k): v for k, v in ckpt["dyadic"].items()}
        success_count = ckpt.get("success_count", 0)
        print(f"Resumed from checkpoint at {last_high:,}")
        if args.start != RESUME_FROM:
            print("  (--start ignored because a checkpoint exists)")
    else:
        last_high     = args.start
        last_twin     = 0
        dyadic        = load_dyadic_history(out, args.start)
        success_count = 0
        print(f"Starting from {args.start:,}")

    small_primes = SmallPrimeCache()

    success_file = out("successful_propagations.csv")
    write_header = not success_file.exists() or success_file.stat().st_size == 0
    success_f = open(success_file, "a", buffering=1024 * 1024)
    if write_header:
        success_f.write("p_n,p_n+1,C,D,C_is_twin,D_is_twin\n")

    current_low = last_high
    prev_twin = last_twin
    twins_since_ckpt = 0
    pending_rows = []
    t_start = time.time()

    # Success rows are buffered and written together with the checkpoint so
    # that the CSV never runs ahead of the checkpointed frontier.
    def persist():
        success_f.writelines(pending_rows)
        pending_rows.clear()
        success_f.flush()
        os.fsync(success_f.fileno())
        write_dyadic_csv(out, dyadic)
        save_checkpoint(out, current_low, prev_twin, dyadic, success_count)

    # The first Ctrl+C (or SIGTERM) lets the current segment finish so the
    # saved state is consistent; a second one aborts as soon as the in-memory
    # state is not mid-update.
    state = {"stop": False, "abort": False, "critical": False}

    def request_stop(signum, frame):
        if state["stop"]:
            state["abort"] = True
            if not state["critical"]:
                raise KeyboardInterrupt
        else:
            state["stop"] = True
            print("\nStop requested - finishing current segment (press Ctrl+C again to abort)")

    signal.signal(signal.SIGINT, request_stop)
    signal.signal(signal.SIGTERM, request_stop)

    try:
        while not state["stop"] and (args.stop_at is None or current_low < args.stop_at):
            seg_low = current_low
            current_high = seg_low + args.segment_size
            if args.stop_at is not None:
                current_high = min(current_high, args.stop_at)
            t0 = time.time()

            twins = lower_twins_in_segment(
                seg_low, current_high, small_primes.covering(current_high + 4)
            ).tolist()

            # Results are collected per segment and applied in one step below,
            # so an abort mid-segment leaves the persisted state consistent.
            seg_dyadic = {}
            seg_rows = []
            seg_pairs = 0
            if twins:
                chain = [prev_twin] + twins if prev_twin > 0 else twins
                for p, q in zip(chain, chain[1:]):
                    C, D, c_twin, d_twin = propagation_status(p, q)
                    succ = c_twin or d_twin
                    update_dyadic(seg_dyadic, p, succ)
                    seg_pairs += 1
                    if succ:
                        seg_rows.append(f"{p},{q},{C},{D},{c_twin},{d_twin}\n")

            state["critical"] = True
            for key, v in seg_dyadic.items():
                entry = dyadic.setdefault(key, {"pairs": 0, "success": 0})
                entry["pairs"] += v["pairs"]
                entry["success"] += v["success"]
            room = max(0, args.max_success_rows - success_count)
            pending_rows.extend(seg_rows[:room])
            success_count += len(seg_rows)
            if twins:
                prev_twin = twins[-1]
                twins_since_ckpt += len(twins)
            current_low = current_high
            if twins_since_ckpt >= args.checkpoint_every:
                persist()
                twins_since_ckpt = 0
            state["critical"] = False
            if state["abort"]:
                raise KeyboardInterrupt

            print(f"Segment [{seg_low:,} -> {current_high:,}]  "
                  f"twins {len(twins):,}  pairs {seg_pairs:,}  "
                  f"propagated {len(seg_rows):,}  total {success_count:,}  "
                  f"({time.time() - t0:.1f}s)")

    except KeyboardInterrupt:
        print("\nAborted - saving last consistent state...")
    finally:
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        persist()
        success_f.close()
        print("Shutdown complete.")
        print(f"Frontier  : {current_low:,}")
        print(f"Successes : {success_count:,}")
        print(f"Elapsed   : {time.time() - t_start:,.1f}s")

    return 130 if state["stop"] else 0


if __name__ == "__main__":
    sys.exit(main())
