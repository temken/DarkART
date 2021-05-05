[![Build Status](https://github.com/temken/DarkARCplusplus/workflows/Build%20Status/badge.svg)](https://github.com/temken/DarkARCplusplus/actions)
[![codecov](https://codecov.io/gh/temken/DarkARCplusplus/branch/master/graph/badge.svg)](https://codecov.io/gh/temken/DarkARCplusplus)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# Dark Matter-induced Atomic Response Code in C++ (DarkARC++)

<!-- [![DOI](https://zenodo.org/badge/202155266.svg)](https://zenodo.org/badge/latestdoi/202155266) -->
[![arXiv](https://img.shields.io/badge/arXiv-1912.08204-B31B1B.svg)](https://arxiv.org/abs/1912.08204)

DarkARC++ is an update of a [previous version](https://github.com/temken/DarkARC) written in python, which was simply too slow.

## GENERAL NOTES

- This code computes the four atomic response functions introduced in the paper [[arXiv:1912.08204]](https://arxiv.org/abs/1912.08204).
- The computations are performed in parallel using [*openmp*](https://www.openmp.org/) library.

## DEPENDENCIES

- [libphysica](https://github.com/temken/libphysica)
- boost
- GSL
- arb
- build with [CMake](https://cmake.org/)
- continuous integration with [Github Actions](https://github.com/actions)
- unit testing with [googletest](https://github.com/google/googletest)
- code coverage with [codecov](https://codecov.io/).

## CONTENT

The included folders are:

## CITING THIS CODE

If you decide to use this code, please cite the latest archived version,

<!-- > [[DOI:10.5281/zenodo.3581334]](https://doi.org/10.5281/zenodo.3581334) -->

as well as the original publications,

>Catena, R., Emken, T. , Spaldin, N., and Tarantino, W., **Atomic responses to general dark matter-electron interactions**, [[arXiv:1912.08204]](https://arxiv.org/abs/1912.08204).

## VERSIONS

## AUTHORS & CONTACT

The author of this tool is Timon Emken.

For questions, bug reports or other suggestions please contact [timon.emken@fysik.su.se](mailto:timon.emken@fysik.su.se).

## LICENSE

This project is licensed under the MIT License - see the LICENSE file.