#!/usr/bin/env bash

#  Created by Michael Marshall on 02/11/2024.
# Sample configure for Mac OS with MacPorts

cd "${0%/*}"
echo "Configuring MLU in $PWD"

BuildDir=build
Prefix=$HOME/.localTestMLU
Pkg=/opt/local # MacPorts packages

# Clean just my directory
#shopt -s globstar
rm -rf autom4te.cache "$BuildDir"
rm aclocal.m4 ar-lib compile config.{guess,sub} configure depcomp install-sh ltmain.sh missing {,{MLU,Analyse}/}Makefile.in MLUconfig.h.in *~

set -e
set -x

autoreconf -fvi
mkdir "$BuildDir"
cd "$BuildDir"
#../configure --with-hdf5="$Pkg" --prefix=$Prefix
../configure CC=clang CXX=clang++ CPPFLAGS="-I$Pkg/include -I$Pkg/include/libomp -Xpreprocessor -fopenmp" LDFLAGS="-L$Pkg/lib -L$Pkg/lib/libomp" LIBS="-lomp" --with-cblas="-framework Accelerate" --with-minuit2="$HOME/.local" --prefix=$Prefix
make -j 20
make install
otool -L "$Prefix"/{bin/{bo,Mu,Co}*,lib/*.*.dylib}

set +e
$Prefix/bin/bootstrap
$Prefix/bin/MultiFit
$Prefix/bin/Continuum
