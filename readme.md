# Meson Lattice Utilities (MLU)
aka
# Mike's Lattice Utilities (MLU)

## Mike's PhD -- Data production and analysis

| Item | Description |
| --- | --- |
| `libMLU.a` | Library containing helper functions (HDF5 wrappers, I/O, filename globbing, ...) | 
| `bootstrap`, `MultiFit`, `Continuum` ... | Generic analysis code |
| `xml3pt` | Production job for generating semilep dataset on Tesseract and Tursa |
| `Scripts` | All analysis jobs to turn correlators into form factors |

## 1. Install Prerequisites

### You must install these dependencies

1. [Grid] Latest version
2. [Hadrons] commit 49fd078f7b8d43d5c6c2c4814e1a8987afc244f0

[grid]: https://github.com/paboyle/Grid
[hadrons]: https://github.com/aportelli/Hadrons

### Install manually -- or let `bootstrap.sh` download

3. [GNU Scientific Library][gsl] (GSL), version 2.5 or later (I'm using 2.7 in production) or from [Mike's PhD page][MikeGSL]
4. ** optional ** [CERN Minuit2][minuit2] standalone version (see paragraph 3, "... can be downloaded [here]") or from [Mike's PhD page][MikeMinuit2]

[gsl]: https://www.gnu.org/software/gsl/
[minuit2]: https://seal.web.cern.ch/seal/MathLibs/Minuit2/html/index.html
[here]: https://seal.web.cern.ch/seal/MathLibs/Minuit2/Minuit2.tar.gz
[MikeMinuit2]: http://lqcd.me/PhD/tar/Minuit2-5.34.14.tar.gz
[MikeGSL]: http://lqcd.me/PhD/tar/gsl-2.7.tar.gz

## 2. Bootstrap

Choose one of the following options:

### Either 2a) Use `bootstrap.sh` to download GSL and Minuit2

    ./bootstrap.sh

`bootstrap.sh` will download its own version of GSL and Minuit2 into `.Package` subdirectory, then run `autoreconf -fvi`

### Or 2b) Shared GSL and optional Minuit2

    autoreconf -fvi

Do this if you:

1. already have GSL installed and wish to use that instance; *or*
2. do not wish to use Minuit2.
  
** NB: This is the only way to build without Minuit2 **

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
 

** Naughty: If either GSL or Minuit2 are not in the path *and* they exist in ../.Package, they will be built and installed by `configure`**

## 4. Build and install

Use the usual incantation:

    make -j 20
    make install

Update your path to include Script subdirectories:

    Utility=*path_to_Utility*
    PATH="$PATH:$Utility/Scripts:$Utility/Scripts/Plot:$Utility/Scripts/Study1Plateau"`

The scripts require a recent version of Bash that supports associative arrays.
** Mac users in particular take note**, make sure the version of Bash you wish to use is earliest in your path.

These scripts require `gnuplot`

    sudo port install gnuplot

The continuum comparison plot scripts use `pdflatex`.
This can be installed on a Mac using [MacTeX][mactex].

[mactex]: https://tug.org/mactex/

## 5. Copy Mike's data

*OPTIONAL* -- if you wish to compare against Mike's data

### 5.1 Download

1. [Mike's desktop analysis dataset][dataset] 18GB. RBC/UKQCD members can get this from [Columbia][dataset_columbia]
2. [Random number cache][cache] 15MB

[dataset]: http://lqcd.me/PhD/Data/NoSync.tar.gz
[dataset_columbia]: https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/marshall/NoSync.tar.gz
[cache]: http://lqcd.me/PhD/Data/MLUCache.tar.gz 

### 5.2 Extract

1. The full dataset to $HOME/NoSync, i.e `cd; tar -xf NoSync.tar.gz`
2. You can extract MLUCache.tgz anywhere you like.
3. `export MLUCache=`*path_to_cache* in .profile. **INCLUDE A TRAILING** `/`
4. Edit the paths in MLUCache/mPi.txt, changing `/Users/mike` to your *home* directory

Apologies in advance for the hardcoded paths.

### 5.3 Test -- continuum fit

Scripts assume NoSync is current directory

    cd ~/NoSync

Backup my version of continuum fits so you can compare against them:

    mv Cont Cont.golden

#### a) Check executables run

    Continuum

Success: Dumps instructions on how to use the continuum fitter.

#### b) Recreate reference continuum fit

    E=3 DisableZ=34 ContFit.sh

Success: `Cont/test/Simul/renormE3-CZ34` should contain fit results (e.g. pdfs).

Compare the fit results against `Ch 10 Reference fit – alternate fit choices for C1` [My continuum fit results][ContFit]

[ContFit]: http://lqcd.me/PhD/Data/Continuum.pdf

Additional success test: compare `F3_K_Ds.corr_f0_fplus.g5P_g5W.model.h5` with version in `Cont.golden/cubic/Simul`, e.g.

    h5diff Cont/test/Simul/renormE3-CZ34/F3_K_Ds.corr_f0_fplus.g5P_g5W.model.h5 Cont.golden/cubic/Simul/renormE3-CZ34/F3_K_Ds.corr_f0_fplus.g5P_g5W.model.h5; echo $?

Should produce no output other than the return code of `0`.

NB: the following two options to h5diff might be useful:

1. `-p 1e-10` relative difference must be greater than 1e-10
2. `-n 10` only show the first 10 differences when comparing files 

#### c) Perform all fit variations

    ContVariations.sh

Success: Cont directory now has many more subdirectories containing alternate fits.

#### d) Plot fit sensitivities

    Make= ContCompare.sh

Success: `Cont` subdirectory contains `Plot` and `Var` subdirectories.
 
Compare the plots against `11 Fit sensitivity / systematic uncertainties` [My continuum fit results][ContFit]

NB: The continuum comparison plot scripts use `pdflatex`. This can be installed on a Mac using [MacTeX][mactex].

## 6. Data production on Tursa

This documentation now shifts gears and explains how data for this PhD are:

1. Produced on the Tursa supercomputer
2. Post-processed on Tursa to produce a manageably small dataset
3. Analysed on a workstation 

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

## 7. Post-processing on Tursa

This section outlines how the `Grid/Hadrons` output is post-processed on Tursa to produce a manageably small dataset.

Results from all of these steps are contained within the `data/ensemble/analyse` subdirectory.

These scripts are each designed to be run from the `data/ensemble` directory.

**Input:** the `data` directory contains one tarred, gzipped file containing the results of each production data run. These can be copied to another machine and the following steps run on these data. NB: You will need the `MakeBS.sh` scripts from each `runs/ensemble` directory as well. These require about 4TB of storage uncompressed.

**Output:** simply copy the `data/analyse` directories from Tursa instead of performing the post-processing described below.

### 7.0 Bootstrap random number seeds and environment variables

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

### 7.1 `bootstrap` – Collecting observables into a Bootstrap

The `MLU` utility `bootstrap` is used to gather the raw data into a single file for each observable.

Because multiple data production runs were performed, sometimes correcting earlier data production errors, the definitive reference for which data are considered good is contained within the `runs/ensemble/MakeBS.sh` script. This script is unique per ensemble and **not** contained within the git repository for `MLU`.

The output of `MakeBS.sh` is a Bash script which needs to be saved, made executable then run to perform the gather. That script uses the `MLU` utility `bootstrap` to produce bootstrap files.

These correspond to the class `Sample<std::complex<double>>` in `MLU/Sample.hpp'. I.e. each temporal component consists of a complex number.

### 7.2 `corr_raw` – Choosing real/imaginary and folding (if 2pt) 
 
The `MLU` script `Scripts/MakeCorrRaw.sh` reads the `bootstrap` directory and creates `corr_raw` using the `fold` MLU utility. 

These correspond to the class `Sample<std::complex<double>>` in `MLU/Sample.hpp'. I.e. each temporal component consists of a complex number.

These correspond to the class `Fold<double>` in `MLU/Fold.hpp'. I.e. each temporal component consists of a double, where the appropriate real/imaginary component has been chosen from the `bootstrap` file – and recorded in the HDF5 (.h5) metadata.

### 7.3 `corr` – Choosing which 2pt data to use

As explained in the thesis, inserting momenta on the quark vs the anti-quark in each meson results in a different smearing.

The `MLU` script `Scripts/MakeCorr.sh` outputs a script which when saved, made executable and run creates the `corr` directory as a set of links into the correct `corr_raw` data.

### 7.4 `ratioE1ZV1` – Choosing which 2pt data to use

The `MLU` script `Scripts/MakeRatio.sh` produces a script, which when saved, made executable and run produces the ratios `R1`, `R2` and `R3` defined in the thesis.

These ratios are defined to require energies and overlap coefficients extracted from 2pt function fits. However for a first analysis, the energies and overlap coefficients are set to 1 with:

    fixe= fixzv= MakeRatio.sh

This produces the `ratioE1ZV1` directory.

Other `ratio` folders are produced later in the analysis, feeding in chosen 2pt function fits. These need to be manually copied/moved into the appropriate `data/ensemble/analyse/ratio...' folder.

## 8 Workstation analysis -- per ensemble

Most workstation analysis scripts source `PlotCommon.sh` which contains common Utilities. In turn, this source `EnsembleX.sh` (where `X` is the actual ensemble script and `Ensemble.sh` is a link to the default ensemble).

Each ensemble script contains hard-coded paths to the location of the post-processed data files (see section 7) on your analysis workstation.   

The workstation analysis scripts take the post processed `ensemble/analyse` data and perform the per-ensemble analysis, i.e. to extract the matrix elements at each kinematic point.

These scripts should work regardless of what the workstation analysis directory is called or where it's located on your workstation – though this is untested. However, my thesis (a separate git repository) assumes it is `$HOME/NoSync'.

These scripts generally consist of an ensemble specific version containing final fit selections for each ensemble. These then call a generic script to do the work.

In the sections that follow, I will describe processing for `C1` ensemble. Repeat the steps below, replacing `C1` with each of the other ensembles `C2`, `F1M`, `M1`, `M2` and `M3`.

### 8.1 Renormalisation

Use `DoAll= Scripts/FitZVC1.sh` to produce `C1/Renorm`

### 8.2 2pt fits

Use `DoAll= DoDisp= Scripts/FitTwoPointC1.sh` to produce `C1/MELFit/2ptp2`

### 8.3 Create ratios using chosen 2pt fits

Use `MakeRatioNoSync.sh` to produce `data/C1/analyse/ratio`

Use `Suffix=AltZV MakeRatioNoSync.sh` to produce `data/C1/analyse/ratioAltZV`

Use `Suffix=Jan24 series=Jan24 MakeRatioNoSync.sh` to produce `data/C1/analyse/ratioJan24`

### 8.4 3pt fits

Use `Scripts/FitMELC1.sh` to produce `C1/MELFit/2ptp2`

### 8.5 Form factors

Use `Scripts/Plot/PlotFormFactorC1.sh` to produce `C1/FormFactor`

## 9 Workstation analysis -- Global chiral continuum fit

### 9.1 Perform all fit variations

    ContVariations.sh

### 9.2 Plot fit sensitivities

    Make= ContCompare.sh

## CBLAS

See `ReadMeCblasVsAccelerate.md` for CBLAS validation and benchmarking.
