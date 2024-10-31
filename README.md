# Algebraic Geometry Codes for Secure Distributed Matrix Multiplication

This repository contains code for constructing the PoleGap SDMM scheme from https://arxiv.org/abs/2303.15429.

## Constructing the scheme

The file `construction.sage` contains all the steps to construct the scheme and test it on randomly chosen matrices. The file works at least on `sage` version 9.5.

## Generating the plots

The file `generate_plot.py` can be used to generate the plots in the paper (Fig. 1). Running this requires `matplotlib` and `numpy`.

## Testing where our scheme is superior

The file `count_statistics.py` can be used to compute how many parameter values gives a better construction for our scheme.
