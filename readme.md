# Meson Lattice Utilities (MLU)
aka
# Mike's Lattice Utilities (MLU)

## Mike's PhD -- Data production and analysis

Item | Description
| --- | ---
`libMLU.a` | Library containing helper functions (HDF5 wrappers, I/O, filename globbing, ...) 
`bootstrap`, `MultiFit`, ... | Generic analysis code  
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

### 1. Shared GSL and optional Minuit2

    autoreconf -fvi

Do this if you:

1. already have GSL installed and wish to use that instance; *or*
2. do not wish to use Minuit2.
  
** NB: This is the only way to build without Minuit2 **

### 2. Use `bootstrap.sh` to download GSL and Minuit2

    ./bootstrap.sh

`bootstrap.sh` will download its own version of GSL and Minuit2 into `.Package` subdirectory, then run `autoreconf -fvi`

## 3. Configure

Use the usual incantation (the arguments are *optional*):

    mkdir build
    cd build
    ../configure --with-grid=*path* --with-hadrons=*path* --with-gsl=*path* --with-minuit2=*path* --prefix=*InstallPath*

** Naughty: If either GSL or Minuit2 are not in the path *and* they exist in ../.Package, they will be built and installed by `configure`**

## 4. Build and install

Use the usual incantation:

    make -j 20
    make install
