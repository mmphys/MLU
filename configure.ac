## Process this file with autoconf to produce a configure script.
# Mike's Lattice Utilities

AC_PREREQ([2.69])
AC_INIT([MLU], [0.1])
AC_COPYRIGHT([Michael Marshall 2021])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR([MLU/CommonGrid.cpp])
#AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([MLUconfig.h])
#AC_CONFIG_HEADERS([MLUconfig.h],[sed -i 's|PACKAGE_|MLU_|' MLUconfig.h])
#AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects 1.13])
AM_INIT_AUTOMAKE([-Wall -Werror foreign 1.13])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

# Optional arguments for Grid / Hadrons.
# If they are installed in same prefix, only one is required (either one)

AC_ARG_WITH([grid],
    [AS_HELP_STRING([--with-grid=<prefix>],
    [Install prefix of Grid])],
    [PATH="$with_grid/bin$PATH_SEPARATOR$PATH"])

AC_ARG_WITH([hadrons],
    [AS_HELP_STRING([--with-hadrons=<prefix>],
    [Install prefix of Hadrons])],
    [if test x"$with_hadrons" != x"$with_grid" ; then]
        [PATH="$with_hadrons/bin$PATH_SEPARATOR$PATH"]
    [fi])

AC_ARG_WITH([gsl],
    [AS_HELP_STRING([--with-gsl=prefix],
    [try this for a non-standard install prefix of the GNU Scientific Library])],
    [AM_CPPFLAGS="-I$with_gsl/include $AM_CPPFLAGS"]
    [AM_LDFLAGS="-L$with_gsl/lib $AM_LDFLAGS"])

AC_ARG_WITH([minuit2],
    [AS_HELP_STRING([--with-minuit2=prefix],
    [try this for a non-standard install prefix of Minuit2])],
    [AM_CPPFLAGS="-I$with_minuit2/include $AM_CPPFLAGS"]
    [AM_LDFLAGS="-L$with_minuit2/lib $AM_LDFLAGS"])

# Make sure Grid is installed

AC_CHECK_PROG([GRIDCONF],[grid-config],[yes])
if test x"$GRIDCONF" != xyes ; then
    AC_MSG_ERROR([grid-config not found])
fi
GRID_PREFIX="`grid-config --prefix`"
if test x"$GRID_PREFIX" == x ; then
    AC_MSG_ERROR([grid-config --prefix returns empty string])
fi
AC_MSG_RESULT([grid-config --prefix=$GRID_PREFIX])

# Make sure Hadrons is installed

AC_CHECK_PROG([HADCONF],[hadrons-config],[yes])
if test x"$HADCONF" != xyes ; then
    AC_MSG_ERROR([hadrons-config not found])
fi
HADRONS_PREFIX="`hadrons-config --prefix`"
if test x"$HADRONS_PREFIX" == x ; then
    AC_MSG_ERROR([hadrons-config --prefix returns empty string])
fi
AC_MSG_RESULT([hadrons-config --prefix=$HADRONS_PREFIX])

# If the user hasn't specified a CXX compiler, get it from Hadrons
if test x"$CXX" == x ; then
    HADRONS_CXX="`hadrons-config --cxx`"
    if test x"$HADRONS_CXX" != x ; then
        AM_CXX="$HADRONS_CXX"
        CXX="$AM_CXX"
        AC_MSG_RESULT([hadrons-config --cxx=$HADRONS_CXX])
    else
        AC_MSG_WARN([hadrons-config --cxx unknown])
    fi
fi

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_AR
AC_LANG([C++])

# Now setup my environment

AM_CPPFLAGS="-I$GRID_PREFIX/include $AM_CPPFLAGS"
AM_LDFLAGS="-L$GRID_PREFIX/lib $AM_LDFLAGS"
if test x"$HADRONS_PREFIX" != x"$GRID_PREFIX" ; then
  AM_CPPFLAGS="$AM_CPPFLAGS -I$HADRONS_PREFIX/include"
  AM_LDFLAGS="$AM_LDFLAGS -L$HADRONS_PREFIX/lib"
fi
AM_CXXFLAGS="`hadrons-config --cxxflags`"
AM_LDFLAGS="$LDFLAGS `hadrons-config --ldflags`"
AM_LIBS="-lHadrons `hadrons-config --libs`"

# TODO: Should I really set these variables here?

AC_MSG_CHECKING([that a minimal Grid and Hadrons program compiles]);
CXX="$AM_CXX"
CPPFLAGS="$AM_CPPFLAGS"
CXXFLAGS="$AM_CXXFLAGS"
LDFLAGS="$AM_LDFLAGS"
LIBS="$AM_LIBS"
AC_LINK_IFELSE(
        [AC_LANG_SOURCE([[
    #include <Grid/Grid.h>
    #include <Hadrons/Application.hpp>

    using namespace Grid;
    using namespace Hadrons;

    int main(int argc, char *argv[])
    {
        Grid_init(&argc, &argv);
        Application application;
        application.run();
        Grid_finalize();
        return 0;
    }

    ]])],
        [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [AC_MSG_ERROR([Could not compile a minimal Grid and Hadrons program])])

# Other pre-requisite libraries

AC_SEARCH_LIBS([gsl_set_error_handler], [gsl], [],
               [AC_MSG_ERROR([GNU Scientific Library not found])])
if test x"$ac_cv_search_gsl_set_error_handler" != "none required" ; then
    AM_LIBS="$AM_LIBS $ac_cv_search_gsl_set_error_handler"
fi

# Not sure how to test for Minuit2 as it's pure C++
AM_LIBS="$AM_LIBS -lMinuit2"
LIBS="$AM_LIBS"

cwd=`pwd -P`; cd ${srcdir}; abs_srcdir=`pwd -P`; cd ${cwd}
AM_CPPFLAGS="-I${abs_srcdir} $AM_CPPFLAGS -I$cwd"
#AM_LDFLAGS="-L${cwd}/MLU $AM_LDFLAGS"

AC_DEFINE_UNQUOTED([MLU_SHA], ["`git rev-parse HEAD`"], [git commit SHA])
AC_DEFINE_UNQUOTED([MLU_SHORT_SHA], ["`git rev-parse --short HEAD`"], [git commit short SHA])
AC_DEFINE_UNQUOTED([MLU_BRANCH], ["`git rev-parse --abbrev-ref HEAD`"], [git branch])

AC_SUBST([MLU_SHORT])
AC_SUBST([MLU_SHA])
AC_SUBST([MLU_BRANCH])
AC_SUBST([AM_CXX])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LIBS])
AC_SUBST([AM_LDFLAGS])
AC_CONFIG_FILES([Makefile MLU/Makefile Analyze/Makefile Study1/Makefile Study2/Makefile])
AC_OUTPUT
