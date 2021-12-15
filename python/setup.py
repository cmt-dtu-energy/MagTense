from setuptools import setup, find_packages

with open("../README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='magtense',
    version='0.1.0',
    description="MagTense - a micromagnetism and magnetostatic framework",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
    ],
    url="www.magtense.org",
    author="cmt-dtu-energy",
    license="GPL 3.0",
    packages=find_packages(include=['magtense', 'magtense.*']),
    package_data={'magtense': ['utils/data/*.csv']},
    include_package_data=True,
)
