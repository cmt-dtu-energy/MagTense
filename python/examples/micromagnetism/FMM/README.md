# dip-fmm micromagnetism example

`fmm_vs_regular.py` and `fmm_vs_regular.ipynb` evolve the same small uniform
micromagnetic problem for 40 ns twice:

1. with the regular MagTense dense demag calculation; and
2. with the persistent dip-fmm plan enabled by `use_cdfmm=True`.

Build the Python extension from the repository root first:

```bash
module load cuda
LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$CONDA_PREFIX/targets/x86_64-linux/lib:$LD_LIBRARY_PATH" \
  make python USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0 USE_FMM3D=0 USE_CDFMM=1
```

Run the CUDA-full example:

```bash
LD_LIBRARY_PATH="$PWD/dip-fmm/local/lib:$CONDA_PREFIX/lib:$CONDA_PREFIX/targets/x86_64-linux/lib:$LD_LIBRARY_PATH" \
  python python/examples/micromagnetism/FMM/fmm_vs_regular.py
```

The same MagTense `cuda=True` input enables the regular CUDA path and selects
dip-fmm's CUDA-full backend, provided MagTense was compiled with `USE_CUDA=1`.

Or perform a CPU-only smoke test without an NVIDIA device. This sets the normal
MagTense `cuda` input to false; dip-fmm then selects oneMKL when it was compiled
in and otherwise uses the portable CPU implementation:

```bash
MAGTENSE_USE_CUDA=0 \
LD_LIBRARY_PATH="$PWD/dip-fmm/local/lib:$CONDA_PREFIX/lib:$LD_LIBRARY_PATH" \
  python python/examples/micromagnetism/FMM/fmm_vs_regular.py
```

The core solver prints initialization and evaluation times separately for each
run. `Evaluation time` covers the ODE solve and its repeated effective-field
evaluations; it excludes construction of the regular demag tensor or persistent
dip-fmm plan. `Full time` is the sum of those two measured phases. The notebook
contains the same calls split into cells so that the regular and FMM
magnetisation trajectories, timings, and final-state difference can be
inspected interactively.
