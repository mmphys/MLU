#!/usr/bin/env bash

# Fail on any error
set -e

echo "Generating configure for Mike's Lattice Utilities (MLU)"
mkdir -p .buildutils/m4
autoreconf -fvi
