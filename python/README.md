# Python Interface

The Fortran code is compiled and wrapped to a module that can be directly called from Python.
The tool `f2py` of the NumPy package is used to wrap the interface file `MagTense/python/magtense/lib/FortranToPythonIO.f90`.

## Deployment with Conda

### Requirements

- Python >= 3.9

- [Intel® C++ Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#inpage-nav-6-undefined) and [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)

  - [Linux] Prepare your terminal, so that ifort compiler can be found:

    ```bash
    . /opt/intel/oneapi/setvars.sh
    ```

  - [Windows] Activation via [VS Code extension for Intel® oneAPI Toolkits](https://github.com/intel/vscode-oneapi-environment-configurator)

- Required python packages

    ```bash
    conda install -y -c intel mkl mkl-include mkl-static
    conda install -y -c "nvidia/label/cuda-${CUDA_VERSION}" cuda-nvcc libcusparse-dev libcublas-dev cuda-cudart-dev
    conda install -y numpy matplotlib
    ```

    - [ Linux ]

        ```bash
        conda install -y -c intel intel-fortran-rt
        ```

    - [ Windows / MacOS ] Installation of Make utility

        ```bash
        conda install -y -c conda-forge make
        ```

- Additional python packages to run data creation scripts

    ```bash
    conda install -y h5py tqdm
    ```

### Preparation of cuda objects

Create importable objects from cuda source code.

HINT: Use `nvcc --version` or `nvidia-smi` to detect the correct CUDA version.

```bash
cd MagTense/source/MagTenseFortranCuda/cuda/
nvcc -shared -Xcompiler -fPIC -c MagTenseCudaBlas.cu
icx -fPIC -c MagTenseCudaBlasICLWrapper.cxx
```

### Installation from source

Create an importable Python module from Fortran source code.

Navigate to folder `MagTense/python/magtense/lib/`, run `make`, and install the package:

```bash
cd MagTense/python/magtense/lib/
make SHELL=cmd
cd MagTense/python/
pip install -e .
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
