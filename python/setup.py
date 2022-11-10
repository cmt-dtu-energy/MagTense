from setuptools import setup, find_packages

setup(
    name='magtense',
    version='2.0.1',
    description="MagTense - a micromagnetism and magnetostatic framework",
    long_description=open('README_top.md').read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering",
    ],
    keywords=[
        "Micromagnetism",
        "Magnetostatics",
        "Demagnetization tensor",
    ],
    url="https://www.magtense.org/",
    project_urls={
	    'Source': 'https://github.com/cmt-dtu-energy/MagTense',
	    'Documentation': 'https://cmt-dtu-energy.github.io/MagTense/',
    },
    author="Stefan Pollok",
    author_email="spol@dtu.dk",
    license="GPL 3.0",
    packages=find_packages(include=['magtense']),
    python_requires=">=3.9, <3.11",
    include_package_data=True,
    install_requires=["numpy", "matplotlib", "intel-fortran-rt==2022.1.0"],
)
