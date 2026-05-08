# FMM Validation Test Suite

This script validates and benchmarks the Fast Multipole Method (FMM) demagnetisation solver against the CUDA reference implementation in `MagTense`.

The test:

- Generates a granular magnetic structure
- Builds or reloads the exchange matrix (`ExMa`)
- Runs a CUDA baseline simulation
- Runs multiple FMM configurations
- Compares demagnetisation fields and magnetisation trajectories
- Computes error metrics
- Saves plots and simulation data
- Validates results against predefined thresholds

---

# Requirements
The script depends on:

- `numpy`
- `pathlib`
- `argparse`
- `MagTense`
- Local project modules:
  - `grain_generator`
  - `fmm_test_setup`
  - `experiment_io`
  - `auxiliary_functions`

# Running the Test

## Standard Run

```bash
python fmm_test_main.py
```

## Quick Test

Runs only the fastest FMM configurations (`nlmax=2`).

```bash
python fmm_test_main.py --quick-test
```

## Custom Seed

Seed is used both for generating the grain structure, and for assiging random magnetizations.

```bash
python fmm_test_main.py --seed 123
```

Default seed value is 100

## External Field Direction

```bash
python fmm_test_main.py --Hext-dir-type 1
```

| Value | Description |
|---|---|
| `1` | Fixed field along +Z |
| `2` | Random but reproducible direction (based on seed) |
| `3` | Fully random direction |

Default is 2.

---

# FMM Configurations

## Default

| `nlmax` | `fmm_nterms` |
|---|---|
| 2 | 6, 10, 15 |
| 3 | 6, 10, 15 |

## Quick Test

| `nlmax` | `fmm_nterms` |
|---|---|
| 2 | 6, 10, 15 |

---

# Generated Files

## Exchange Matrix  & Demagnetisation Tensor 

Generated once and reused automatically:

```text
seed_<seed>_exMa.npz

seed_<seed>_demag.bin
```

## Timer and Trace Logs

Examples:

```text
cuda_100_timer.log
FMM_L2_N10_100_timer.log

cuda_100_trace.log
FMM_L3_N15_100_trace.log
```

Note that "trace" is deactivated and must be manually enabled in the fmm_test_setup.py file

## Saved Experiment Data

Saved in results/seed_<seed_value>/

Includes:
- simulation results
- error metrics
- comparison plots

---

# Validation Thresholds

The thresholds are arbitrarily set-

## `nlmax = 2`

| Metric | Threshold |
|---|---|
| `Hdem_rel_avg` | `1e-3` |
| `Hdem_rel_max` | `0.05` |
| `Mout_max` | `0.03` |

## `nlmax = 3`

| Metric | Threshold |
|---|---|
| `Hdem_rel_avg` | `1e-2` |
| `Hdem_rel_max` | `0.4` |
| `Mout_max` | `0.25` |

---

# Exit Codes

| Code | Meaning |
|---|---|
| `0` | All tests passed |
| `1` | Validation failed |

Suitable for CI pipelines.

---

# Example

```bash
python fmm_test_main.py \
    --quick-test \
    --seed 42 \
    --Hext-dir-type 2
```