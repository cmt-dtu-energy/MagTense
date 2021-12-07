import setuptools
from numpy.distutils.core import Extension, setup
from numpy.distutils.fcompiler import get_default_fcompiler
from magtense import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

# fortran compiler
compiler = get_default_fcompiler()
# set some fortran compiler-dependent flags
f90flags = []
if compiler == "gnu95":
    #f90flags.append("-fPIC -O3 -fopenmp -fdefault-real-8 -ffree-line-length-512")
    f90flags.append("-fopenmp")
    f90flags.append("-fdefault-real-8")
    f90flags.append("-ffree-line-length-512")
elif compiler == "intel" or compiler == "intelvem":
    f90flags = "/nologo /O3 /assume:nocc_omp /Qopenmp /real-size:64 /fp:precise /libs:static /threads"

ext = Extension(
    name="magtensesource",
    sources=[
        "lib_mag/FortranToPythonIO.f90"
    ],
    library_dirs=["lib_mag/"],
    extra_f90_compile_args=f90flags,
    #extra_objects=["lib_mag/"],
    extra_link_args=["-lsrc"],
)

setup(
    name="magtense",
    version=__version__,
    description="test",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
    ],
    url="https://github.com/zhucaoxiang/CoilPy",
    author="Caoxiang Zhu",
    author_email="caoxiangzhu@gmail.com",
    license="GNU 3.0",
    packages=setuptools.find_packages(),
    package_data={'magtense': ['../util/*']},
    ext_modules=[ext],
)