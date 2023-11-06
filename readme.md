# Meson Lattice Utilities (MLU)
aka
# Mike's Lattice Utilities (MLU)

## Mike's PhD -- Data production and analysis

Item | Description
| --- | ---
`libMLU.a` | Library containing helper functions (HDF5 wrappers, I/O, filename globbing, ...) 
`bootstrap`, `MultiFit`, `Continuum` ... | Generic analysis code  
`xml3pt` | Production job for generating semilep dataset on Tesseract and Tursa
`Scripts` | All analysis jobs to turn correlators into form factors

## 1. Install Prerequisites

### You must install these dependencies

1. [Grid] Latest version
2. [Hadrons] commit 49fd078f7b8d43d5c6c2c4814e1a8987afc244f0

[grid]: https://github.com/paboyle/Grid
[hadrons]: https://github.com/aportelli/Hadrons

### Install manually -- or let `bootstrap.sh` download

3. [GNU Scientific Library][gsl] (GSL), version 2.5 or later (I'm using 2.7 in production) or from [Mike's home page][MikeGSL]
4. ** optional ** [CERN Minuit2][minuit2] standalone version (see paragraph 3, "... can be downloaded [here]") or from [Mike's home page][MikeMinuit2]

[gsl]: https://www.gnu.org/software/gsl/
[minuit2]: https://seal.web.cern.ch/seal/MathLibs/Minuit2/html/index.html
[here]: https://seal.web.cern.ch/seal/MathLibs/Minuit2/Minuit2.tar.gz
[MikeMinuit2]: https://www2.ph.ed.ac.uk/~s1786208/Minuit2-5.34.14.tar.gz
[MikeGSL]: https://www2.ph.ed.ac.uk/~s1786208/gsl-2.7.tar.gz

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

1. [Mike's desktop analysis dataset][dataset] 9.1GB
2. [Random number cache][cache] 15MB

[dataset]: https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/marshall/NoSync.tar.gz
[cache]: https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/marshall/MLUCache.tgz 

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

    DisableZ=34 ContFit.sh

Success: `Cont/Simul/renorm-CZ34` should contain fit results (e.g. pdfs).

Compare the fit results against `Ch 10 Reference fit – alternate fit choices for C1` [My continuum fit results][ContFit]

[ContFit]: https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/marshall/SemiLep/Continuum.pdf

#### c) Perform all fit variations

    ContVariations.sh

Success: Cont directory now has many more subdirectories containing alternate fits.

#### d) Plot fit sensitivities

    Make= ContCompare.sh

Success: `Cont` subdirectory contains `Plot` and `Var` subdirectories.
 
Compare the plots against `11 Fit sensitivity / systematic uncertainties` [My continuum fit results][ContFit]

NB: The continuum comparison plot scripts use `pdflatex`. This can be installed on a Mac using [MacTeX][mactex].
