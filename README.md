[![Build Status](https://github.com/temken/DarkART/workflows/Build%20Status/badge.svg)](https://github.com/temken/DarkART/actions)
[![codecov](https://codecov.io/gh/temken/DarkART/branch/main/graph/badge.svg)](https://codecov.io/gh/temken/DarkART)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

# Dark Atomic Response Tabulator (DarkART)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6046224.svg)](https://doi.org/10.5281/zenodo.6046224)
[![PRR](https://img.shields.io/badge/Phys.Rev.Research-2(2020),033195-255773.svg)](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.033195)
[![arXiv](https://img.shields.io/badge/arXiv-1912.08204-B31B1B.svg)](https://arxiv.org/abs/1912.08204)


<img src="https://github.com/temken/DarkART/blob/dev/logo.png?raw=true" width="400">

DarkART is a C++ tool for the computation and tabulation of atomic response functions for direct sub-GeV dark matter (DM) searches.
It supersedes the previous python tool [DarkARC](https://github.com/temken/DarkARC) for better performance and improved software design.

## GENERAL NOTES

- This code computes the four atomic response functions introduced in the paper [[arXiv:1912.08204]](https://arxiv.org/abs/1912.08204). The first response function is also known as the ionization form factor.
- The computations are performed in parallel using [*openmp*](https://www.openmp.org/) library.

<img src="https://user-images.githubusercontent.com/29034913/70995423-d0683c80-20d0-11ea-85bd-fdcb91d972eb.png" width="600">

<details><summary>Repository content</summary>
<p>

The included folders are:

- *bin/*: This folder contains the executable after successful installation together with the configuration files.
- *data/*: Contains files including the RHF coefficients of the initial electron wavefunctions.
- *external/*: This folder will only be created and filled during the build with CMake and will contain the [libphysica](https://github.com/temken/libphysica) library.
- *include/*: All header files of DarkART can be found here.
- *results/*: Each run of DarkART generates result files in a dedicated sub-folder named after the run's ID, which is specified in the configuration file.
- *src/*: Here you find the source code of DarkART.
- *tests/*: All code and executable files of the unit tests are stored here.

</p>
</details>

## DEPENDENCIES

- [arb](https://arblib.org/)
- [boost](https://www.boost.org/)
- [CMake](https://cmake.org/)
- [GSL](https://www.gnu.org/software/gsl/)
- [libconfig](https://github.com/temken/libphysica)
- [libphysica](https://github.com/temken/libphysica)
- [openmp](https://www.openmp.org/)

The *libphysica* library will be downloaded and built automatically using CMake.
The other libraries need to be installed beforehand.

<details><summary>Installation of dependencies</summary>
<p>

On macOS, you can install the dependencies using [homebrew](https://brew.sh/),

```
>brew install arb
>brew install boost
>brew install cmake
>brew install gsl
>brew install libconfig
>brew install libomp
```

or on a Linux machine with APT:

```
>sudo apt-get update
>sudo apt-get install arb
>sudo apt-get install libboost-all-dev
>sudo apt-get install cmake
>sudo apt-get install libgsl-dev
>sudo apt-get install -y libconfig-dev
>sudo apt-get install libomp-dev
```

Alternatively, *libconfig* can be also built from the source files via

```
>wget https://hyperrealm.github.io/libconfig/dist/libconfig-1.7.2.tar.gz
>tar -xvzf libconfig-1.7.2.tar.gz
>pushd libconfig-1.7.2
>./configure
>make
>sudo make install
>popd
```

</p>
</details>

## DOWNLOAD & INSTALLATION

The DarkART source code can be downloaded by cloning this git repository:

```
>git clone https://github.com/temken/DarkART.git 
>cd DarkART
```

The code is compiled and the executable is created using CMake.

```
>cmake -E make_directory build
>cd build
>cmake -DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF ..
>cmake --build . --config Release
>cmake --install .
```

If everything worked well, there should be the executable *DarkART* in the */bin/* folder.

## CITING THIS CODE

If you decide to use this code, please cite the latest archived version,

> [[DOI:10.5281/zenodo.6046225]](https://doi.org/10.5281/zenodo.6046225)

<details><summary>Bibtex entry</summary>
<p>

```
@software{DarkART,
  author = {Emken, Timon},
  title = {{Dark Atomic Response Tabulator (DarkART)[Code, v0.1.0]}},
  year         = {2021},
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {DOI:10.5281/zenodo.6046225},
  url          = {https://doi.org/10.5281/zenodo.6046225},
  howpublished={The code can be found under \url{https://github.com/temken/darkart}.}
}
```

</p>
</details>

as well as the original publications,

>Catena, R., Emken, T. , Spaldin, N., and Tarantino, W., **Atomic responses to general dark matter-electron interactions**, [Phys.Rev.Research 2 (2020) 033195](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.033195), [[arXiv:1912.08204]](https://arxiv.org/abs/1912.08204).



## VERSIONS

- 11.02.2021: Release of version 0.1.0

## AUTHORS & CONTACT

The author of this tool is Timon Emken.

For questions, bug reports or other suggestions please contact [timon.emken@fysik.su.se](mailto:timon.emken@fysik.su.se) or, even better, open an [issue](https://github.com/temken/DarkART/issues).

## LICENSE

This project is licensed under the MIT License - see the LICENSE file.