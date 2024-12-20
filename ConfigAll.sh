#!/usr/bin/env bash

#  Created by Michael Marshall on 02/11/2024.

set -e

cd "${0%/*}"
echo "Configuring SemiLep and MLU in $PWD"

Config=${Config-Debug}
BuildDir=build$Config
Prefix=$HOME/.local$Config
Pkg=/opt/local # MacPorts packages

if [ "${1@L}" = clean ]
then
shopt -s globstar
rm -rf **/autom4te.cache {,MLU/}build* || :
rm **/{aclocal.m4,ar-lib,compile,config.guess,config.sub,configure,depcomp,install-sh,ltmain.sh,missing,Makefile.in,*onfig.h.in,*~} || :

autoreconf -fvi
fi
set -x

mkdir "$BuildDir"
cd "$BuildDir"
#../configure CXX=mpicxx-openmpi-mp CC=mpicc-openmpi-mp CPPFLAGS="-I$Pkg/include -I$Pkg/include/libomp -Xpreprocessor -fopenmp" LDFLAGS="-L$Pkg/lib -L$Pkg/lib/libomp" LIBS="-lomp" --with-cxx=clang++ --with-cc=clang --with-cblas="-framework Accelerate" --with-minuit2="$HOME/.local" --prefix="$Prefix"
../configure CC=clang CPPFLAGS="-I$Pkg/include -I$Pkg/include/libomp -Xpreprocessor -fopenmp" LDFLAGS="-L$Pkg/lib -L$Pkg/lib/libomp" LIBS="-lomp" --with-cblas="-framework Accelerate" --with-minuit2="$HOME/.local" --prefix="$Prefix"
make -j 20
make install
otool -L "$Prefix"/{bin/{Co,De,xml3}*,lib/*.*.dylib}

set +e
$Prefix/bin/Continuum
$Prefix/bin/Debug
$Prefix/bin/xml3pt
