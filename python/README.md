# Python Interface

The Fortran code is compiled and wrapped to a module that can be directly called from Python. The tool **f2py** of the NumPy package is used to wrap the interface file **lib/FortranToPythonIO.f90**.

## Deployment with Conda

### Step 1

Installation of required python packages

```bash
conda install -y numpy matplotlib
```

### Step 2 - Magnetostatic part only

GFortran compiler and Make utility (Windows + MacOS only)

- Windows

  - Installation in conda environment

    ```bash
    conda install -y -c conda-forge make
    conda install -y -c msys m2w64-gcc-fortran
    ```

  - Installation from binary | [MinGW](https://gcc.gnu.org/wiki/GFortranBinaries#Windows)

- MacOS:
  - Installation from binary | [HPC Mac OS X](http://hpc.sourceforge.net/)
  - Installation with [Homebrew](https://brew.sh/) ( **brew install gcc** )


### Step 2 - Magnetostatic + Micromagnetic part

INFO: Works currently only with Windows and Python >= 3.9 | On Linux, it compiles and runs, but results diverge

- Software
    - [Intel® Fortran Compiler](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)
    - [Intel® oneAPI Math Kernel Library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)

- Linux
    #### Compilation Time

    Prepare your terminal, so that ifort and MKL can be found.

    ```bash
    . ~/intel/oneapi/setvars.sh
    export LDFLAGS=-Wl,-rpath=/opt/intel/oneapi/mkl/latest/lib/intel64/
    ```
    
    #### Runtime

    Environment variables of oneapi have to be set.

    ```bash
    . ~/intel/oneapi/setvars.sh
    ```
    

### Step 3

Creation of an importable Python module from Fortran source code

Navigate to folder **MagTense/python/magtense/lib/**, run **make**, and install the package

```bash
cd /path/to/repo/python/magtense/lib/
make
cd /path/to/repo/python/
pip install -e .
```


## Read-in customized M-H-curve
This feature is currently only supported for soft magnetic tiles ([type=2](magtense/magtense.py#L49)).

In  [iterate_magnetization()](magtense/magtense.py#L611), an arbitrary number of state functions (M-H-curves) can be defined:

```python
mu_r = 100
datapath = f'./magtense/utils/data/Fe_mur_{mu_r}_Ms_2_1.csv'

 ...

data_statefcn = numpy.genfromtxt(datapath, delimiter=';')
n_statefcn = 1
```

[Here](magtense/utils/data), three sample M-H-curves for Fe with different relative permeabilities and a saturation magnetization of 2.1 T are stored as CSV-files. The data format is as follows:

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
