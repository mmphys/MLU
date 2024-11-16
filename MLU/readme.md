# Meson (aka Mike's) Lattice Utilities (MLU)

This package contains the MLU Library and Analysis Utilities used for performing the data analysis for [Mike's PhD](http://lqcd.me/PhD/).

Follow these instructions to build a stand-alone version of the MLU Library and Analysis Utilities.

(For the full Semileptonic Data Generation package see [`readme.md`](../readme.md) in the containing directory.)

## Folder structure

| Item | Description |
| --- | --- |
| `MLU` | Library `libMLU.a` supports file formats for bootstrap (and jackknife) samples, fits, models, etc. Also contains helper functions (HDF5 wrappers, I/O, filename globbing, ...) | 
| `Analyse` | `bootstrap`, `CRatio`, `Divide` ... code for performing each step of analysis |
| `Fit` | `MultiFit` allows multiple correlators to be fitted to multiple models simultaneously, optionally constraining values to the results of prior fits. `Continuum` performs the final continuum fit. |
| `Scripts` | All analysis jobs to turn correlators into form factors |

# Building

## 1. Install dependencies

Mandatory

1. [HDF5]
2. [GNU Scientific Library][gsl] (GSL), version 2.5 or later (I used 2.7 in production) or from [Mike's PhD page][MikeGSL]
3. [Bash] with associative arrays (version >= 4.2, latest preferred)

NB: GSL can be installed as a subpackage of MLU (see below).

Optional

1. **highly recommended** [OpenMP] (except for debug builds)
2. ** optional ** [CERN `Minuit2`][minuit2] standalone version (see paragraph 3, "... can be downloaded [here]") or from [Mike's PhD page][MikeMinuit2]

[bash]: https://www.gnu.org/software/bash/
[hdf5]: https://www.hdfgroup.org/solutions/hdf5/
[openmp]: https://www.openmp.org
[gsl]: https://www.gnu.org/software/gsl/doc/html/index.html
[minuit2]: https://seal.web.cern.ch/seal/MathLibs/Minuit2/html/index.html
[here]: https://seal.web.cern.ch/seal/MathLibs/Minuit2/Minuit2.tar.gz
[MikeMinuit2]: http://lqcd.me/PhD/tar/Minuit2-5.34.14.tar.gz
[MikeGSL]: http://lqcd.me/PhD/tar/gsl-2.7.tar.gz

## 2. Bootstrap

Run

    bootstrap.sh

Downloads archives containing the versions of `GSL` and `Minuit2` used in production, but **doesn't** expand them.

### 2a. Optional: Build `GSL` as a subpackage of MLU

**Skip this step if you already have `GSL` installed.**

If GSL exists in the subdirectory `gsl-2.7`, MLU will build it as a subpackage:

    tar -xf gsl-2.7.tar.gz
    cd gsl-2.7
    autoreconf -fvi

NB: You can put any version of GSL in this directory, but it must be named `gsl-2.7`

## 3. Configure

Perform the usual incantation

    mkdir build
    cd build
    ../configure

MacOS example configure script can be found in `Config.sh`

More help is available

    ../configure --help=recursive

## 4. Build and install

Use the usual incantation:

    make -j 20
    make install

Update your path to include Script subdirectories:

    MLU=*path_to_MLU*
    PATH+=":$MLU/Scripts:$MLU/Scripts/Plot:$MLU/Scripts/Study1Plateau"`

The scripts require a recent version of Bash (>= 4.2) that supports associative arrays.
** Mac users in particular take note**, make sure the version of Bash you wish to use is earliest in your path.

### Script dependencies

These scripts require `gnuplot`

    sudo port install gnuplot

The continuum comparison plot scripts use `pdflatex`.
This can be installed on a Mac using [MacTeX][mactex].

[mactex]: https://tug.org/mactex/

# Workstation data analysis

We now explain how the [SemiLep](../readme.md) data are analysed on a workstation to produce a final result:

1. Per-ensemble analysis
2. Global chiral continuum fit

## 1 Workstation analysis -- per ensemble

Most workstation analysis scripts source `PlotCommon.sh` which contains common Utilities. In turn, this source `EnsembleX.sh` (where `X` is the actual ensemble script and `Ensemble.sh` is a link to the default ensemble).

Each ensemble script contains hard-coded paths to the location of the post-processed data files (see section 7) on your analysis workstation.   

The workstation analysis scripts take the post processed `ensemble/analyse` data and perform the per-ensemble analysis, i.e. to extract the matrix elements at each kinematic point.

These scripts should work regardless of what the workstation analysis directory is called or where it's located on your workstation – though this is untested. However, my thesis (a separate git repository) assumes it is `$HOME/NoSync'.

These scripts generally consist of an ensemble specific version containing final fit selections for each ensemble. These then call a generic script to do the work.

In the sections that follow, I will describe processing for `C1` ensemble. Repeat the steps below, replacing `C1` with each of the other ensembles `C2`, `F1M`, `M1`, `M2` and `M3`.

### 1.1 Renormalisation

Use `DoAll= Scripts/FitZVC1.sh` to produce `C1/Renorm`

### 1.2 2pt fits

Use `DoAll= DoDisp= Scripts/FitTwoPointC1.sh` to produce `C1/MELFit/2ptp2`

### 1.3 Create ratios using chosen 2pt fits

Use `MakeRatioNoSync.sh` to produce `data/C1/analyse/ratio`

Use `Suffix=AltZV MakeRatioNoSync.sh` to produce `data/C1/analyse/ratioAltZV`

Use `Suffix=Jan24 series=Jan24 MakeRatioNoSync.sh` to produce `data/C1/analyse/ratioJan24`

### 1.4 3pt fits

Use `Scripts/FitMELC1.sh` to produce `C1/MELFit/2ptp2`

### 1.5 Form factors

Use `Scripts/Plot/PlotFormFactorC1.sh` to produce `C1/FormFactor`

## 2 Workstation analysis -- Global chiral continuum fit

### 2.1 Perform all fit variations

    ContVariations.sh

### 2.2 Plot fit sensitivities

    Make= ContCompare.sh

# Compare against reference data

*OPTIONAL* -- if you wish to compare against Mike's PhD data analysis

### 1 Download

1. [Mike's desktop analysis dataset][dataset] 18GB. RBC/UKQCD members can get this from [Columbia][dataset_columbia]
2. [Random number cache][cache] 15MB

[dataset]: http://lqcd.me/PhD/Data/NoSync.tar.gz
[dataset_columbia]: https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/marshall/NoSync.tar.gz
[cache]: http://lqcd.me/PhD/Data/MLUCache.tar.gz 

### 2 Extract

1. The full dataset to $HOME/NoSync, i.e `cd; tar -xf NoSync.tar.gz`
2. You can extract MLUCache.tgz anywhere you like.
3. `export MLUCache=`*path_to_cache* in .profile. **INCLUDE A TRAILING** `/`
4. Edit the paths in MLUCache/mPi.txt, changing `/Users/mike` to your *home* directory

Apologies in advance for the hardcoded paths.

### 3 Test -- continuum fit

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

# CBLAS

See [`ReadMeCblasVsAccelerate.md`](ReadMeCblasVsAccelerate.md) for CBLAS validation and benchmarking.
