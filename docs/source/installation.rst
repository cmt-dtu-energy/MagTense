Installation
==============================================

==============================================
Matlab
==============================================

----------------------------------------------
Prerequisites
----------------------------------------------
* Matlab :math:`\geq` 2018a

MagTense is directly useable in Matlab by downloading the already compiled
`MEX-files <https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/MEX_files>`_.


==============================================
Python
==============================================

----------------------------------------------
Prerequisites
----------------------------------------------
* Python :math:`\geq` 3.9
* `Intel® Fortran Compiler <https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran>`_
* `Intel® oneAPI Math Kernel Library <https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html>`_

----------------------------------------------
Required Packages
----------------------------------------------
[Info] On Linux, Make utility is pre-installed.

::

    conda install -y numpy matplotlib
    conda install -y -c conda-forge make

----------------------------------------------
Prepare terminal
----------------------------------------------
Setting environment variables to recognize ifort and MKL.
[Info] On Windows, the environment variables might be set already.

::

    . ~/intel/oneapi/setvars.sh

----------------------------------------------
Build Python module
----------------------------------------------
Fortran source code becomes accessible from Python.

::

    cd /path/to/repo/python/src/magtense/lib/
    make SHELL=cmd
    cd /path/to/repo/python/
    pip install -e .
