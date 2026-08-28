# dip-fmm micromagnetism sweep

The examples use non-periodic spherical dip-fmm plans for all combinations of:

- grids: `15^3`, `20^3`, `25^3`, and `30^3` cells;
- expansion orders: 1 through 10; and
- tree depths: 2 through 5.

Build the Python extension from the repository root first:

```bash
module load cuda
LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$CONDA_PREFIX/targets/x86_64-linux/lib:$LD_LIBRARY_PATH" \
  make python USE_CUDA=1 USE_CVODE=1 USE_MATLAB=0 USE_FMM3D=0 USE_CDFMM=1
```

Prepare all 160 persistent plans without running time evolution:

```bash
LD_LIBRARY_PATH="$PWD/dip-fmm/local/lib:$CONDA_PREFIX/lib:$CONDA_PREFIX/targets/x86_64-linux/lib:$LD_LIBRARY_PATH" \
  python python/examples/micromagnetism/FMM/prepare_fmm_cache.py
```

dip-fmm validates and loads an existing plan instead of rebuilding it. Plans
are backend-independent, so preparation uses the CPU by default and avoids
needless GPU transfers. Set `MAGTENSE_USE_CUDA=1` only if you also want to
exercise the production CUDA setup path while warming the cache.

Run the 40 ns comparison sweep:

```bash
LD_LIBRARY_PATH="$PWD/dip-fmm/local/lib:$CONDA_PREFIX/lib:$CONDA_PREFIX/targets/x86_64-linux/lib:$LD_LIBRARY_PATH" \
  python python/examples/micromagnetism/FMM/fmm_vs_regular.py
```

The regular demagnetisation result is evaluated once per grid size, then reused
for all order/depth comparisons. The core reports evaluation time separately
from initialization time. Set `MAGTENSE_USE_CUDA=0` for CPU execution; dip-fmm
then prefers oneMKL and falls back to its portable CPU backend.

`fmm_vs_regular.ipynb` runs the same sweep interactively and plots final-state
relative RMS error versus order for every grid size and depth.
