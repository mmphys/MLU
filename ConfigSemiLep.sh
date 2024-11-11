#!/usr/bin/env bash

#  Created by Michael Marshall on 02/11/2024.

cd "${0%/*}"
echo "Configuring SemiLep in $PWD"

BuildDir=buildSemiLep
Prefix=$HOME/.localTestSemiLep
Pkg=/opt/local # MacPorts packages
HadronsDir="${GridPre}Release"

# Clean just my directory
#shopt -s globstar
rm -rf autom4te.cache "$BuildDir"
rm aclocal.m4 ar-lib compile config.{guess,sub} configure depcomp install-sh ltmain.sh missing {,SemiLep/}Makefile.in SemiLepConfig.h.in *~

set -e
set -x

autoreconf -fvi --no-recursive
mkdir "$BuildDir"
cd "$BuildDir"
#../configure CC=clang CPPFLAGS="-I$Pkg/include/libomp -Xpreprocessor -fopenmp" LDFLAGS="-L$Pkg/lib/libomp" LIBS="-lomp" --with-MLU="$HOME/.localTestMLU" --with-lime="$HOME/.local"  --with-hdf5="$Pkg" --prefix="$Prefix"
../configure CC=clang --with-MLU="$HOME/.localTestMLU" --prefix="$Prefix"
make -j 20
make install
otool -L "$Prefix"/bin/{li,De,xml3}*

set +e
$Prefix/bin/lime_truncate
$Prefix/bin/Debug
$Prefix/bin/xml3pt
