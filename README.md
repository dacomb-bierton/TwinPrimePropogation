# True Twin-Propagation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22040224.svg)](https://doi.org/10.5281/zenodo.22040224)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22017371.svg)](https://doi.org/10.5281/zenodo.22017371)

High-performance computational engine for testing the **Twin-Prime Propagation Conjecture**.

Paper: [A Twin-Prime Propagation Conjecture](https://doi.org/10.5281/zenodo.22040224)  
Author: Dacomb Bierton (August 2026)

---

## The Conjecture

Let \( p_n \) denote the \( n \)-th **lower twin prime** (the smaller prime in a twin-prime pair).

For consecutive lower twin primes \( p_n \) and \( p_{n+1} \), define:

\[
C_n = p_n + p_{n+1} + 1, \qquad D_n = p_n + p_{n+1} + 3
\]

The pair \( (p_n, p_{n+1}) \) is said to **propagate** if at least one of \( C_n \) or \( D_n \) is itself a lower twin prime.

**Twin-Prime Propagation Conjecture**  
There are infinitely many indices \( n \) for which the pair \( (p_n, p_{n+1}) \) propagates.

Extensive computation up to \( 10^{13} \) shows that the success rate declines slowly and regularly. The absolute number of successful propagations continues to increase across dyadic intervals, consistent with the conjectured infinitude. A conditional argument under a sufficiently strong form of the Hardy–Littlewood prime tuples conjecture is outlined in the paper.

---

## What This Code Does

This program systematically searches for consecutive lower twin primes and tests the propagation condition:

1. Uses a vectorised **segmented sieve** (odd numbers only) to generate primes efficiently over large intervals. The small-prime table grows automatically with the search frontier, so the sieve stays exact however far the run goes.
2. Identifies lower twin primes (pairs differing by 2), including pairs that straddle a segment boundary.
3. For every pair of consecutive lower twin primes \( (p, q) \), computes \( C = p + q + 1 \) and \( D = p + q + 3 \).
4. Tests whether \( C \) or \( D \) is a member of a twin-prime pair using fast primality tests (gmpy2 when available, otherwise a Miller–Rabin test that is deterministic for all 64-bit inputs).
5. Records successful propagations and tracks empirical success rates in **dyadic intervals**.
6. Supports **safe checkpointing and resumption**. Results are committed one whole segment at a time and success rows are written together with the checkpoint, so an interrupted run resumes without losing or duplicating pairs.

### Output files

All output is written to the `true_twin_propagation/` directory (configurable with `--output-dir`):

| File | Description |
|------|-------------|
| `true_successful_propagations.csv` | Successful propagating pairs (capped for size) |
| `true_propagation_dyadic_rates.csv` | Success rates by dyadic interval |
| `true_checkpoint.pkl` | Resume state |

---

## Requirements

- Python 3.8+
- `numpy`
- Optional but strongly recommended: `gmpy2` (significantly faster primality testing)

```bash
pip install -r requirements.txt
```

---

## Usage

```bash
python twin_prime_prop.py
```

The script will:

- Resume from the last checkpoint if one exists, or start from `--start` (default \( 2^{38} \)).
- Process successive large segments.
- Checkpoint periodically.
- Handle clean interruption (`Ctrl+C`) by finishing the current segment and saving state. A second `Ctrl+C` aborts immediately; the saved state is still consistent.

Every setting can be given on the command line (defaults shown):

```text
--start            274877906944   lower bound when no checkpoint exists
--stop-at          (none)         stop once the frontier reaches this value
--segment-size     150000000      numbers sieved per segment
--checkpoint-every 40000          checkpoint after this many new twin primes
--max-success-rows 4000000        cap on rows written to the successes CSV
--output-dir       true_twin_propagation
```

For example, a short self-contained run from the origin that is easy to verify by hand:

```bash
python twin_prime_prop.py --start 0 --stop-at 10000000 --segment-size 1000000 --output-dir /tmp/tpp-check
```

---

## Citation

If you use this software, please cite it as below.

**APA**

> Bierton, D. (2026). *Twin Prime Propagation* (v1.0). Zenodo. https://doi.org/10.5281/zenodo.22040224

**BibTeX**

```bibtex
@software{Bierton_Twin_Prime_Propagation_2026,
  author       = {Bierton, Dacomb},
  title        = {Twin Prime Propagation},
  month        = aug,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0},
  doi          = {10.5281/zenodo.22040224},
  url          = {https://doi.org/10.5281/zenodo.22040224}
}
```

A `CITATION.cff` file is also included in this repository.

---

## License

Copyright © 2026 Dacomb Bierton

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

You are free to:

- **Share** — copy and redistribute the material in any medium or format
- **Adapt** — remix, transform, and build upon the material for any purpose, even commercially

Under the following terms:

- **Attribution** — You must give appropriate credit, provide a link to the license, and indicate if changes were made.

---

## Related Work

- Full paper and supporting materials: [https://doi.org/10.5281/zenodo.22040224](https://doi.org/10.5281/zenodo.22040224)
