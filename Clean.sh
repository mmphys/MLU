#!/usr/bin/env bash

#  Created by Michael Marshall on 02/11/2024.

cd "${0%/*}"
echo "Cleaning $PWD"

shopt -s globstar
rm -rf **/autom4te.cache {,MLU/}build*
rm **/{aclocal.m4,ar-lib,compile,config.guess,config.sub,configure,depcomp,install-sh,ltmain.sh,missing,Makefile.in,*onfig.h.in,*~}
