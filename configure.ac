## Process this file with autoconf to produce a configure script.
# Semileptonic data production (SemiLep) and Meson Lattice Utilities (MLU)

AC_INIT([SemiLep], [0.2])
AC_COPYRIGHT([Michael Marshall 2021-2024])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_SRCDIR([SemiLep/xml3pt.cpp])
AM_INIT_AUTOMAKE([foreign])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

# Optional arguments for Grid / Hadrons.
# If they are installed in same prefix, only one is required (either one)

GridConfig="grid-config"
AC_ARG_WITH([grid],
    [AS_HELP_STRING([--with-grid=<prefix>],
    [Install prefix of Grid])],
    [AC_MSG_NOTICE([   --with-grid=$with_grid])]
    [GridConfig="$with_grid/bin/$GridConfig"],
    [AC_CHECK_PROG([GRIDCONF],[$GridConfig],[yes],[no])]
    [if test x"$GRIDCONF" != xyes ; then]
        [AC_MSG_ERROR([$GridConfig not found])]
    [fi]
    [AC_MSG_CHECKING([   Grid --prefix])]
    [with_grid="`$GridConfig --prefix`"]
    [if test x"$with_grid" = x ; then]
        [AC_MSG_ERROR([grid-config --prefix returns empty string])]
    [fi]
    [AC_MSG_RESULT([$with_grid])])

HadronsConfig="hadrons-config"
AC_ARG_WITH([hadrons],
    [AS_HELP_STRING([--with-hadrons=<prefix>],
    [Install prefix of Hadrons])],
    [AC_MSG_NOTICE([--with-hadrons=$with_hadrons])]
    [HadronsConfig="$with_hadrons/bin/$HadronsConfig"],
    [AC_CHECK_PROG([HADCONF],[$HadronsConfig],[yes],[no])]
    [if test x"$HADCONF" != xyes ; then]
        [AC_MSG_ERROR([$HadronsConfig not found])]
    [fi]
    [AC_MSG_CHECKING([Hadrons --prefix])]
    [with_hadrons="`$HadronsConfig --prefix`"]
    [if test x"$with_hadrons" = x ; then]
        [AC_MSG_ERROR([hadrons-config --prefix returns empty string])]
    [fi]
    [AC_MSG_RESULT([$with_hadrons])])

AC_ARG_ENABLE([check], AS_HELP_STRING([--disable-check], [Disable check of Grid and Hadrons build]))

# If the user hasn't specified a CXX compiler, get it from Hadrons

if test x"$CXX" = x ; then
    CXX="`$HadronsConfig --cxx`"
    if test x"$CXX" != x ; then
        AC_MSG_NOTICE([hadrons-config --cxx=$CXX])
    else
        AC_MSG_WARN([hadrons-config --cxx unknown])
    fi
fi
if test x"$CXXLD" = x ; then
    CXXLD="`$HadronsConfig --cxxld`"
    if test x"$CXXLD" != x ; then
        AC_MSG_NOTICE([hadrons-config --cxxld=$CXXLD])
    else
        AC_MSG_WARN([hadrons-config --cxxld unknown])
    fi
fi

# Check for programs
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX(17, [noext], [mandatory])
AC_PROG_CC
AM_PROG_AR
LT_INIT
AC_LANG([C++])

# Decode --with-MLU option
AC_ARG_WITH([MLU],
  [AS_HELP_STRING([--with-MLU],
    [yes to build MLU sub-package, no assumes MLU can be found in environment,
    otherwise path to installation location of MLU @<:@default=yes@:>@])],
  [],
  [: m4_divert_text([DEFAULTS], [with_MLU=yes])])

# --with-MLU[=yes] means we're building MLU subpackage
# Otherwise --with-MLU=optional prefix to pre-installed MLU
AM_CONDITIONAL([WANT_BUILD_MLU], [test "x$with_MLU" = xyes])
AS_IF([test "x$with_MLU" = xyes],
      [AC_CONFIG_SUBDIRS([MLU])],
      [AS_IF([test "x$with_MLU" != xno],
            [CPPFLAGS+=" -I$with_MLU/include"]
            [LDFLAGS+=" -L$with_MLU/lib"])]
)

# Optional path to hdf5_cpp header (we can rely on MLU for link to library)
AC_ARG_WITH([hdf5], [AS_HELP_STRING([--with-hdf5], [path to installation location of hdf5])],
  AS_IF([test "x$with_hdf5" != xyes && test "x$with_hdf5" != xno],
        [CPPFLAGS+=" -I$with_hdf5/include"]))

# Optional path to GSL header (we can rely on MLU for link to libraries)
AC_ARG_WITH([gsl], [AS_HELP_STRING([--with-gsl], [path to installation location of GSL])],
  AS_IF([test "x$with_gsl" != xyes && test "x$with_gsl" != xno],
        [CPPFLAGS+=" -I$with_gsl/include"]))

# Optional path to lime
AC_ARG_WITH([lime], [AS_HELP_STRING([--with-lime], [path to installation location of lime])],
  AS_IF([test "x$with_lime" != xyes && test "x$with_lime" != xno],
      [CPPFLAGS+=" -I$with_lime/include"]
      [LDFLAGS+=" -L$with_lime/lib"]
      [AC_CHECK_LIB([lime],[limeCreateReader],,[AC_MSG_FAILURE([lime library not found])])]
      ))

# Now setup my environment
CPPFLAGS="-I$with_hadrons/include $CPPFLAGS"
LDFLAGS="-L$with_hadrons/lib $LDFLAGS"
if test x"$with_hadrons" != x"$with_grid" ; then
  CPPFLAGS="-I$with_grid/include $CPPFLAGS"
  LDFLAGS="-L$with_grid/lib $LDFLAGS"
fi

CXXFLAGS="`$HadronsConfig --cxxflags` $CXXFLAGS"
LDFLAGS="`$HadronsConfig --ldflags` $LDFLAGS"
LIBS="-lHadrons `$HadronsConfig --libs` $LIBS"

# Grid and Hadrons use CXXFLAGS (compiler) but should have used CPPFLAGS (preprocessor)
# Copy all the -I flags from CXXFLAGS to CPPFLAGS
[CPPFLAGS="`echo $CXXFLAGS | sed -E 's/(^| )([^-]|-[^I])[^ ]*//g'` $CPPFLAGS"]
# Remove all the -I flags from CXXFLAGS
[CXXFLAGS="`echo $CXXFLAGS | sed -E 's/(^| )-I[^ ]*//g'`"]
# If CXXFlags contains -Xpreprocessor -fopenmp, put it in CPPFLAGS
[CPPFLAGS="$CPPFLAGS `echo $CXXFLAGS | sed -E 's/.*(-Xpreprocessor -fopenmp).*/-Xpreprocessor -fopenmp/'`"]
# Remove -Xpreprocessor -fopenmp from CXXFlags
[CXXFLAGS="`echo $CXXFLAGS | sed -E 's/-Xpreprocessor -fopenmp//g'`"]

AC_MSG_NOTICE([CPPFLAGS=$CPPFLAGS])
AC_MSG_NOTICE([LDFLAGS=$LDFLAGS])
AC_MSG_NOTICE([CXXFLAGS=$CXXFLAGS])
AC_MSG_NOTICE([LIBS=$LIBS])
AC_MSG_NOTICE([CXXLD=$CXXLD])

# See whether OpenMP is available
AC_OPENMP

AS_IF([test x"$enable_check" != xno], [
    AC_MSG_CHECKING([that a minimal Grid and Hadrons program compiles]);
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
    ])

#AC_DEFINE_UNQUOTED([SemiLep_SHA], ["`git rev-parse HEAD`"], [git commit SHA])
#AC_DEFINE_UNQUOTED([SemiLep_SHORT_SHA], ["`git rev-parse --short HEAD`"], [git commit short SHA])
#AC_DEFINE_UNQUOTED([SemiLep_BRANCH], ["`git rev-parse --abbrev-ref HEAD`"], [git branch])
AC_DEFINE_UNQUOTED([SemiLep_GIT_SUMMARY], ["MLU `git rev-parse --abbrev-ref HEAD` `git rev-parse HEAD`, Grid `$GridConfig --git`, Hadrons `$HadronsConfig --git`"], [git branch and commit for MLU, Grid and Hadrons])

CCLD="$CXXLD"
AC_SUBST([AM_CPPFLAGS],["$CPPFLAGS"])
AC_SUBST([AM_CXXFLAGS],["$CXXFLAGS"])
AC_SUBST([AM_LDFLAGS],["$LDFLAGS $LIBS"])
AC_SUBST([CXXLD])
AC_SUBST([CCLD])
AC_SUBST([SemiLep_GIT_SUMMARY])
CPPFLAGS=""
CXXFLAGS=""
LDFLAGS=""
LIBS=""
AC_CONFIG_HEADERS([SemiLepConfig.h],
    [$SED_INPLACE -e 's|#define PACKAGE|#define SemiLep_PACKAGE|g' -e 's|#define VERSION|#define SemiLep_VERSION|g' SemiLepConfig.h],
    [if test x"$host_os" = x"${host_os#darwin}" ; then]
        [SED_INPLACE="sed -i"]
    [else]
        [SED_INPLACE="sed -i .bak"]
    [fi])
AC_CONFIG_FILES([Makefile SemiLep/Makefile])
AC_OUTPUT
