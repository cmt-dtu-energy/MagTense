Installation
==============================================

==============================================
Matlab
==============================================
MagTense is directly useable in Matlab by downloading the 
already compiled `MEX-files <https://github.com/cmt-dtu-energy/MagTense/tree/master/matlab/MEX_files>`_. The files are directly useable 
with no compilation required, although Matlab 2018a or greater 
is required.

==============================================
Python
==============================================

----------------------------------------------
Prerequisites
----------------------------------------------
* python installation
* numpy package
* matplotlib
* gfortran :math:`\geq` 5

----------------------------------------------
Add scripts to PYTHONPATH
----------------------------------------------
* python/lib_mag/
* python/source/

----------------------------------------------
Build python module from Fortran source code 
----------------------------------------------
::

    cd python/lib_mag/
    make


==============================================
Standalone
==============================================
An EXE-file is available `here <https://github.com/cmt-dtu-energy/MagTense/tree/master/executable>`_.
It can be used with a `Python script <https://github.com/cmt-dtu-energy/MagTense/blob/master/python/source/MagTenseStandalone.py>`_.
