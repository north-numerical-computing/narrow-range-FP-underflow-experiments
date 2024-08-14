# MATLAB experiments with narrow-range floating-point formats
This repository contains the data plotted in [Sec. 5, 1] and the MATLAB code used to generate it.

The scripts rely on the [CPFloat](https://github.com/north-numerical-computing/cpfloat) library for implementing custom-precision arithmetics.

The errors in [data](./data) were produced with MATLAB version R2024a Update 2 and the CPFloat version under the SHA [0ddba87dd19adf5da154e228b68f86e60051f102](https://github.com/north-numerical-computing/cpfloat/tree/0ddba87dd19adf5da154e228b68f86e60051f102).
Run the script [experiments.m](./experiments.m) to regenerate the data files in [data](./data).

### References

 [1] T. Mary and M. Mikaitis. [*Error Analysis of Matrix Multiplication with Narrow Range Floating-Point Arithmetic*](https://). In Prep. Aug. 2024.

### License

This software is distributed under the terms of the 2-clause BSD software license (see [LICENSE](./LICENSE)).

The CPFloat C library is distributed under the [GNU Lesser General Public License, Version 2.1 or later](https://raw.githubusercontent.com/mfasi/cpfloat/master/LICENSES/LGPL-2.1-or-later.txt).
