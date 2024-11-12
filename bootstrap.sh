#!/usr/bin/env bash

# Fail on any error
set -e

echo "SemiLep (semileptonic) data production and Meson Lattice Utilities (MLU)"

cd MLU
./bootstrap.sh
cd ..

echo "Generating configure for SemiLep (semileptonic) data production"
autoreconf -fvi
