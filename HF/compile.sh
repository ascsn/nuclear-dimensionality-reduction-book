#!/bin/bash

# This script compiles the HF code and the module for the RBM notebook

make

f2py -c skyrme_hpsi.f90 -m skyrme_hpsi

mv skyrme_hpsi.cpython*.so ../RBM/skyrme_hpsi.so
