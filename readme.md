# Semileptonic Data Generation (SemiLep)

This project contains semileptonic data production code for [Mike's PhD](http://lqcd.me/PhD/).

If you wish to perform data analysis only, **don't** build this package. Instead, build the [subpackage](MLU/readme.md): Meson (aka Mike's) Lattice Utilities (`MLU`).

## Folder structure

| Item | Description |
| --- | --- |
| `MLU` | [Subpackage](MLU/readme.md) Meson (aka Mike's) Lattice Utilities (`MLU`) | 
| `SemiLep` | `xml3pt` generates the SemiLep dataset on Tursa (and previously Tesseract). Includes other utilities with dependencies on `Grid` and `Hadrons` |
| `xml` | Sample xml used for data production | 

# Building

## 1. Install dependencies

### You must install these mandatory dependencies

You must install the mandatory dependencies for `MLU`

1. [HDF5]
2. [OpenMP]
3. [GNU Scientific Library][gsl] (`GSL`), version 2.5 or later (I used 2.7 in production) or from [Mike's PhD page][MikeGSL]
4. [Bash] with associative arrays (version >= 4.2, latest preferred)

and additional dependencies for `SemiLep`

5. [c-lime][lime]
6. [Grid] Latest version (with HDF5, OpenMP and c-lime)
7. [Hadrons] commit 49fd078f7b8d43d5c6c2c4814e1a8987afc244f0

This commit of Hadrons requires `-fpermissive` to be included in `CXXFLAGS` in order to build successfully with g++ (`-Xcompiler -fpermissive` for nvcc). It also has a dependency on a Grid deflation header which has moved, requiring (assuming Grid has been installed in `$Prefix`)

    cd ~/$Prefix/include/Grid/algorithms/iterative
    ln -s ../deflation/Deflation.h

[hdf5]: https://www.hdfgroup.org/solutions/hdf5/
[openmp]: https://www.openmp.org
[gsl]: https://www.gnu.org/software/gsl/doc/html/index.html
[MikeGSL]: http://lqcd.me/PhD/tar/gsl-2.7.tar.gz
[bash]: https://www.gnu.org/software/bash/
[lime]: https://usqcd-software.github.io/c-lime/
[grid]: https://github.com/paboyle/Grid
[hadrons]: https://github.com/aportelli/Hadrons

## 2. Bootstrap

Choose one of the following options:

#### 2a. Build both `SemiLep` and `MLU` together

    ./bootstrap.sh

#### 2b. Build `SemiLep` only (if you already have `MLU` installed)

    autoreconf -fvi --no-recursive

## 3. Configure

Use the usual incantation

    mkdir build
    cd build

In `configure`, if you wish to build `MLU` as a subpackage, set `--with-MLU=yes` (or leave it out entirely - `yes` is the default). Otherwise set `--with-MLU=*path*` to link to a pre-installed copy of `MLU` in `*path*`, or '--with-MLU=no' if `MLU` can be found via `CPPFLAGS` etc.

*Optionally* set `with-grid=*path*`, `with-hadrons=*path*`, `--with-hdf5=*path*`, `--with-gsl=*path*` and `--with-lime=*path*`. E.g.

    ../configure --with-grid=*path* --with-hadrons=*path* --with-gsl=*path* --prefix=*InstallPath*

If you need to specify additional paths, use `CPPFLAGS`, `LDFLAGS`, `LIBS` etc.
E.g. to pick up prerequisites installed in a location specified by `$GridPre`: 

    ../configure CPPFLAGS=-I$GridPre/include LDFLAGS=-L$GridPre/lib --with-grid=*path* --with-hadrons=*path* --prefix=*InstallPath*

You can also specify `CXX`, `CC`, `CXXLD`, `CCLD` to choose specific c++ and c compilers and linkers. If not set, `CXX` and `CXXLD` will be taken from 'hadrons-config'. `configure` will also do its best to apportion `hadrons-config --cxxflags` between CPPFLAGS (c & c++ pre-processor flags) and CXXFLAGS (c++ compiler flags). E.g. if you also wish to use Xcode clang and tell it to use OpenMP installed in $GridPkg (download the OpenMP library separately, e.g. with MacPorts):

    ../configure CXX=clang++ CC=clang CPPFLAGS="-I$GridPkg/include -I$GridPkg/include/libomp -Xpreprocessor -fopenmp" LDFLAGS="-L$GridPkg/lib -L$GridPkg/lib/libomp" LIBS=-lomp --with-grid=*path* --with-hadrons=*path* --prefix=*InstallPath*

MacOS example configure scripts can be found in

* `ConfigAll.sh` configure SemiLep and MLU as a single package
* `ConfigSemiLep.sh` configure SemiLep only

If you wish to build `SemiLep` for gpu, it is recommended that `MLU` still be built for cpu so that the `MLU` utilities (`bootstrap`, `MultiFit`, `Continuum` etc) can be run from the shell. If building both at the same time, specify different compiler settings for `MLU` using `--with-cxx`, `--with-cc`, `--with-cxxld`, `--with-ccld`, `--with-cxxflags`, `--with-cppflags`, `--with-ldflags` and `--with-libs`.  

E.g. building for gpu on Tursa (assuming Grid, Hadrons and GSL are all available in $Prefix)

    ../configure CXX= CPPFLAGS="-I$Prefix/include" LDFLAGS="-L$Prefix/lib" --with-grid=$Prefix --with-hadrons=$Prefix --with-cxx=g++ --with-cxxld=g++ --prefix=$Prefix

More help is available

    ../configure --help=recursive

## 4. Build and install

Use the usual incantation:

    make -j 20
    make install

or, using Xcode, open `SemiLep.xcodeproj` and build the scheme `All`.

# Semileptonic data production

This is a very brief introduction to how data for this PhD are produced on the (Tursa) supercomputer.

The project folder for heavy-light semileptonics is:

    /mnt/lustre/tursafs1/home/dp207/dp207/shared/projects/semilep

This contains 3 folders:

1. `runs` contains scripts to produce and post-process data
2. `data` contains the results of those runs
3. `study` contains prior studies

Within `data` and `runs` we have subdirectories for each ensemble, `C1`, `C2`, `F1M`, `M1`, `M2` and `M3`.

`runs` contains a `code` subdirectory containing versions 2, 3 and 4 of the code used for data production (version 1 ran on Tesseract). Each version contains a `sourceme.sh` Bash script (for GPU) or a `sourcemecpu.sh` Bash script (for the CPU nodes).

The `data` directories for each ensemble generally contain 1 or more `run` directories containing the output from `Grid/Hadrons`, as well as an `analyse` directory containing post processing.

Within each `data/run` directory there are subdirectories for each configuration, each of which containing a `log` subdirectory containing the `.sh` script to perform the run, the `.log` from the run, the `.xml` input to the executable run from the script, etc.

The source code for the production runs is `SemiLep/xml3pt.cpp` with a sample `.xml` in `Xml/xml3pt.xml`

NB: The `MakeBS.sh` scripts from each `runs/ensemble` subdirectory are the definitive record of which data comprise the data set.

## Post-processing and analysis

See the [MLU documentation](MLU/readme.md).
