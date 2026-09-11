#!/usr/bin/env bash
# Script to setup up the MagTense Development Environment
#
# Trying to setup a Miniconda enviroment in a macOS with AppleSilicon
set -e
build_dir=$HOME/build
mkdir -p $build_dir

cwd=$(pwd)
cd "$build_dir"
prefix="$HOME/.local/miniconda"
mkdir -p "$prefix"
curl -o miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
chmod +x ./miniconda.sh
./miniconda.sh -bfp "$prefix"
rm miniconda.sh

export PATH="$prefix/bin:$PATH"
source "$prefix/bin/activate"
conda init bash
conda create --name magtense python=3.11
conda activate magtense
conda install magtense -c cmt-dtu-energy/label/cpu -c https://software.repos.intel.com/python/conda/ -c conda-forge
