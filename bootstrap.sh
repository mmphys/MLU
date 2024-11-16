#!/usr/bin/env bash

# Fail on any error
set -e

cd "${0%/*}"
echo "Bootstrapping SemiLep (semileptonic data generation) in $PWD"

[ -f MLU/bootstrap.sh ] && MLU/bootstrap.sh

echo "Generating configure for SemiLep (semileptonic data generation)"
autoreconf -fvi
