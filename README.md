# True Twin-Propagation

[![DOI](https://zenodo.org/badge/1341432193.svg)](https://doi.org/10.5281/zenodo.22040224)

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

1. Uses a **segmented sieve** to generate primes efficiently over large intervals.
2. Identifies twin primes (pairs differing by 2).
3. For every pair of consecutive lower twin primes \( (p, q) \), computes \( C = p + q + 1 \) and \( D = p + q + 3 \).
4. Tests whether \( C \) or \( D \) is itself a lower twin prime using fast primality tests (gmpy2 when available, otherwise deterministic Miller-Rabin).
5. Records successful propagations and tracks empirical success rates in **dyadic intervals**.
6. Supports **safe checkpointing and resumption**, allowing long-running searches to be interrupted and continued.

### Output files

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
pip install numpy gmpy2
```

---

## Usage

```bash
python true_twin_propagation.py
```

The script will:

- Resume from the last checkpoint if one exists, or start from the configured `RESUME_FROM` value.
- Process successive large segments.
- Checkpoint periodically.
- Handle clean interruption (`Ctrl+C`) by saving state.

Key configuration constants are at the top of the script:

```python
SEGMENT_SIZE     = 150_000_000
CHECKPOINT_EVERY = 40_000
RESUME_FROM      = 274_877_906_944
MAX_SUCCESS_ROWS = 4_000_000
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

A `CITATION.cff` file is also included in this repository for automated citation tools.

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
```
