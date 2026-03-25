#!/usr/bin/bash


# Remove all build files
rm -rf build/
rm -rf source/_build/

# Rebuild
python3 -m sphinx -E -b html source build/html
