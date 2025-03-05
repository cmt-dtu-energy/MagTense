#!/bin/bash

# Navigate to wrapper file for python module
# cd $RECIPE_DIR/../python/src/magtense/lib
# make
# TODO Copy compiled library to python module

# Navigate to distribution script
cd $RECIPE_DIR/../python
$PYTHON -m pip install . -vvv --use-deprecated=legacy-resolver --no-deps --ignore-installed --no-cache-dir

# Clean up
cd $RECIPE_DIR/../python/src/magtense/lib
make clean
rm -r $RECIPE_DIR/../python/src/magtense.egg-info $RECIPE_DIR/../python/build
