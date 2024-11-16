# Semileptonic Data Generation (SemiLep)

This project contains semileptonic data production code for [Mike's PhD](http://lqcd.me/PhD/).

If you wish to perform data analysis only, **don't** install this package. Instead, install the [subpackage](MLU/readme.md): Meson (aka Mike's) Lattice Utilities (MLU).

## Folder structure

| Item | Description |
| --- | --- |
| `MLU` | [Subpackage](MLU/readme.md) Meson (aka Mike's) Lattice Utilities (MLU) | 
| `SemiLep` | `xml3pt` generates the SemiLep dataset on Tursa (and previously Tesseract). Includes other utilities with dependencies on `Grid` and `Hadrons` |
| `xml` | Sample xml used for data production | 

# Building

## 1. Install dependencies

### You must install these mandatory dependencies

1. [HDF5]
2. [OpenMP]
3. [c-lime][lime]
4. [Grid] Latest version (with HDF5, OpenMP and c-lime)
5. [Hadrons] commit 49fd078f7b8d43d5c6c2c4814e1a8987afc244f0

This commit of Hadrons requires `-fpermissive` to be included in `CXXFLAGS` in order to build successfully. It also has a dependency on a Grid deflation header which has moved, requiring (assuming Grid has been installed in `$Prefix`)

    cd ~/$Prefix/include/Grid/algorithms/iterative
    ln -s ../deflation/Deflation.h

[hdf5]: https://www.hdfgroup.org/solutions/hdf5/
[openmp]: https://www.openmp.org
[lime]: https://usqcd-software.github.io/c-lime/
[grid]: https://github.com/paboyle/Grid
[hadrons]: https://github.com/aportelli/Hadrons

## 2. Bootstrap

Choose one of the following options:

#### 2a. Build both `SemiLep` and `MLU` together

    ./bootstrap.sh

`bootstrap.sh` will download its own version of GSL and Minuit2 into `.Package` subdirectory, then run `autoreconf -fvi`

#### 2b. Build `SemiLep` only (if you already have `MLU` installed)

    autoreconf -fvi --no-recursive

## 3. Configure

If you have built Grid with OpenMP, HDF5 and lime, then use the usual incantation (the arguments are *optional*):

    mkdir build
    cd build
    ../configure --with-grid=*path* --with-hadrons=*path* --with-gsl=*path* --with-minuit2=*path* --prefix=*InstallPath*

If you need to specify additional paths, use `CXXFLAGS`, `LDFLAGS`, etc.
E.g. to pick up prerequisites installed in a location specified by environment variable `GridPre`: 

    ../configure CXXFLAGS=-I$GridPre/include LDFLAGS=-L$GridPre/lib --with-grid=*path* --with-hadrons=*path* --with-gsl=*path* --with-minuit2=*path* --prefix=*InstallPath*

E.g. if you also wish to use Xcode clang and tell it to use OpenMP (download the OpenMP library separately, e.g. with MacPorts):

    ../configure CXX=clang++ CC=clang CXXFLAGS="-I$GridPre/include -I$GridPre/include/libomp -Xpreprocessor -fopenmp" LDFLAGS="-L$GridPre/lib -L$GridPre/lib/libomp" LIBS=-lomp --with-grid=*path* --with-hadrons=*path* --with-gsl=*path* --with-minuit2=*path* --prefix=*InstallPath*
 
MacOS example configure scripts can be found in

* `ConfigAll.sh` configure SemiLep and MLU as a single package
* `ConfigSemiLep.sh` configure SemiLep only

More help is available

    ../configure --help=recursive

## 4. Build and install

Use the usual incantation:

    make -j 20
    make install

# Semileptonic data production

This documentation now explains how data for this PhD are:

1. Produced on the Tursa supercomputer
2. Post-processed on Tursa to produce a manageably small dataset
3. Analysed on a workstation (see [MLU](MLU/readme.md))

## Data production on Tursa

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

## Post-processing on Tursa

This section outlines how the `Grid/Hadrons` output is post-processed on Tursa to produce a manageably small dataset.

Results from all of these steps are contained within the `data/ensemble/analyse` subdirectory.

These scripts are each designed to be run from the `data/ensemble` directory.

**Input:** the `data` directory contains one tarred, gzipped file containing the results of each production data run. These can be copied to another machine and the following steps run on these data. NB: You will need the `MakeBS.sh` scripts from each `runs/ensemble` directory as well. These require about 4TB of storage uncompressed.

**Output:** simply copy the `data/analyse` directories from Tursa instead of performing the post-processing described below.

### 1. Bootstrap random number seeds and environment variables

Post processing utilities rely on the following environment variables (search for `getenv` in `MLU/JackBoot.cpp`):

| Environment variable | Description | Default |
| --- | --- | --- |
| `MLUSeed` | Integer random number seed for bootstrap or "`Jackknife`". `MLUSeed` is used each time the file is loaded to perform the bootstrap (or jackknife). | 1835672416 |
| `MLUHost` | Name of the host the random number seeds were generated on | Local host name. |
| `MLUCache` | Directory of cache containing random number files. Include a trailing '/', otherwise the trailing text becomes a filename prefix. | Cache file names consist of five parts: *MLUHost*.*NumSamples*.`Random`.*MLUSeed*.`h5` |
| `MLUNumBoot` | Number of bootstrap replicas | 10000 |
| `MLUFat` | When defined, causes random numbers to be saved with each sample file | Random numbers stored in cache only. |
| `MLUCovarCentral` | When defined, the central replica is used as the estimate of the mean when computing the covariance matrix. | Mean of the binned data is used as estimate of mean. |

**NB:** for Mike's original analysis, `MLUHost=Tursa` and bootstrap random number seeds for each ensemble were (`Scripts/Plot/EnsembleC1.sh` etc):

| Ensemble | Seed |
| :-: | --: |
| C1 | 2263212701 |
| C2 | 694228835 |
| M1 | 1835672416 |
| M2 | 2475361470 |
| M3 | 3301941204 |
| F1M | 3285139232 |

### 2. `bootstrap` – Collecting observables into a Bootstrap

The `MLU` utility `bootstrap` is used to gather the raw data into a single file for each observable.

Because multiple data production runs were performed, sometimes correcting earlier data production errors, the definitive reference for which data are considered good is contained within the `runs/ensemble/MakeBS.sh` script. This script is unique per ensemble and **not** contained within the git repository for `MLU`.

The output of `MakeBS.sh` is a Bash script which needs to be saved, made executable then run to perform the gather. That script uses the `MLU` utility `bootstrap` to produce bootstrap files.

These correspond to the class `Sample<std::complex<double>>` in `MLU/Sample.hpp'. I.e. each temporal component consists of a complex number.

### 3. `corr_raw` – Choosing real/imaginary and folding (if 2pt) 
 
The `MLU` script `Scripts/MakeCorrRaw.sh` reads the `bootstrap` directory and creates `corr_raw` using the `fold` MLU utility. 

These correspond to the class `Sample<std::complex<double>>` in `MLU/Sample.hpp'. I.e. each temporal component consists of a complex number.

These correspond to the class `Fold<double>` in `MLU/Fold.hpp'. I.e. each temporal component consists of a double, where the appropriate real/imaginary component has been chosen from the `bootstrap` file – and recorded in the HDF5 (.h5) metadata.

### 4. `corr` – Choosing which 2pt data to use

As explained in the thesis, inserting momenta on the quark vs the anti-quark in each meson results in a different smearing.

The `MLU` script `Scripts/MakeCorr.sh` outputs a script which when saved, made executable and run creates the `corr` directory as a set of links into the correct `corr_raw` data.

### 5. `ratioE1ZV1` – Choosing which 2pt data to use

The `MLU` script `Scripts/MakeRatio.sh` produces a script, which when saved, made executable and run produces the ratios `R1`, `R2` and `R3` defined in the thesis.

These ratios are defined to require energies and overlap coefficients extracted from 2pt function fits. However for a first analysis, the energies and overlap coefficients are set to 1 with:

    fixe= fixzv= MakeRatio.sh

This produces the `ratioE1ZV1` directory.

Other `ratio` folders are produced later in the analysis, feeding in chosen 2pt function fits. These need to be manually copied/moved into the appropriate `data/ensemble/analyse/ratio...' folder.
