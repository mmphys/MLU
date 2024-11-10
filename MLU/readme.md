# Meson Lattice Utilities (MLU)

This subpackage of `SemiLep` contains the MLU Library and Utilities.

# Dependencies

Mandatory

1. `HDF5`
2. `GSL`

Optional

1. `Open MP`
2. `Minuit2`

NB: `Open MP` is highly recommended (except for debug builds).

# Building

Follow these instructions to build a stand-alone version of the MLU Library and analysis utilities.

See `readme.md` in the containing directory to build the entire package, i.e. including data generation.

## 1. Bootstrap

If you already have [GSL](https://www.gnu.org/software/gsl/doc/html/index.html) on your system

    mkdir gsl-2.7
    autoreconf -fvi

Otherwise, you can make GSL as a subpackage of MLU

    bootstrap.sh

which will download `GSL` and `Minuit2` and expand `GSL` into the `gsl-2.7` directory.
By default, `GSL` will be made from this source.

If `GSL` is already on your system, you can delete the **contents** of this directory (**the `gsl-2.7` directory must exist**) then run `autoreconf -fvi` 

## 2. Configure and build

Perform the usual incantation `./configure && make && make install`

For MacOS a more complete configure and build can be found in `Config.sh`

More help is available

    ../configure --help=recursive
