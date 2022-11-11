# Python Interface

The Fortran code is compiled and wrapped to a module that can be directly called from Python.
The tool `f2py` of the NumPy package is used to wrap the interface file `lib/FortranToPythonIO.f90`.

## Deployment with Conda

### Requirements

- Python >= 3.9

- [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)

- [Intel® oneAPI Math Kernel Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)

- [Linux] Prepare your terminal, so that ifort and MKL can be found:

    ```bash
    . ~/intel/oneapi/setvars.sh
    ```

- [Windows + MacOS] Installation of Make utility

    ```bash
    conda install -y -c conda-forge make
    ```

- Required python packages

    ```bash
    conda install -y numpy matplotlib
    ```

- Additional python packages to run data creation scripts

    ```bash
    conda install -y h5py tqdm
    ```

### Installation from source

Create an importable Python module from Fortran source code.
Navigate to folder `python/magtense/lib/`, run `make`, and install the package:

```bash
cd /path/to/repo/python/magtense/lib/
make SHELL=cmd
cd /path/to/repo/python/
pip install -e .
```

## [Linux] Set LD_LIBRARY_PATH for specific conda environment only

Adapted from https://stackoverflow.com/questions/46826497/conda-set-ld-library-path-for-env-only.

1. Create these subdirectories and files in the directory of the specific conda environment:
    ```sh
    cd /CONDA_ENV_PATH/
    mkdir -p ./etc/conda/activate.d
    mkdir -p ./etc/conda/deactivate.d
    touch ./etc/conda/activate.d/env_vars.sh
    touch ./etc/conda/deactivate.d/env_vars.sh
    ```
2. Edit `./etc/conda/activate.d/env_vars.sh`
    ```sh
    #!/bin/sh

    export OLD_LD_PRELOAD=${LD_PRELOAD}
    export LD_PRELOAD=/CONDA_ENV_PATH/lib/libmkl_core.so:/CONDA_ENV_PATH/lib/libmkl_intel_lp64.so:/CONDA_ENV_PATH/lib/libmkl_intel_thread.so:/CONDA_ENV_PATH/lib/libiomp5.so:${LD_PRELOAD}
    export OLD_LD_LIBRARY_PATH=${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH=/CONDA_ENV_PATH/lib/:${LD_LIBRARY_PATH}
    ```

3. Edit `./etc/conda/deactivate.d/env_vars.sh`
    ```sh
    #!/bin/sh

    export LD_PRELOAD=${OLD_LD_PRELOAD}
    export LD_LIBRARY_PATH=${OLD_LD_LIBRARY_PATH}
    unset OLD_LD_PRELOAD
    unset OLD_LD_LIBRARY_PATH
    ```

## Read-in customized M-H-curve

This feature is currently only supported for soft magnetic tiles ([type=2](magtense/magtense.py#L49)).

In  [iterate_magnetization()](magtense/magtense.py#L611), an arbitrary number of state functions (M-H-curves) can be defined:

```python
mu_r = 100
datapath = f'./magtense/mat/Fe_mur_{mu_r}_Ms_2_1.csv'

 ...

data_statefcn = numpy.genfromtxt(datapath, delimiter=';')
n_statefcn = 1
```

[Here](magtense/mat), three sample M-H-curves for Fe with different relative permeabilities and a saturation magnetization of 2.1 T are stored as CSV-files. The data format is as follows:

```csv
0; Temp0; Temp1; ...
H0-field; M0@Temp0; M0@Temp1;...
H1-field; M1@Temp0; M1@Temp1;...
.
.
H100-field; M100@Temp0; M100@Temp1; ...
.
```

With only one state function given, the same M-H-curve applies to all tiles of type 2.

When the soft tiles differ in their M-H-curves, multiple state function can be combined. In order to match a specific M-H-curve with the corresponding tile, the variable [stfcn_index](magtense/magtense.py#L54) can be set.
