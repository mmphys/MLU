# MLU

AC_INIT([MLU], [0.2])
AC_COPYRIGHT([Michael Marshall 2021-2024])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_SRCDIR([MLU/Physics.cpp])
AM_INIT_AUTOMAKE([foreign])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

# Check for programs
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX(17, [noext], [mandatory])
AC_PROG_CC
AM_PROG_AR
LT_INIT
AC_LANG([C++])

# See whether OpenMP is available
AC_OPENMP

# These are the configure options for callers of MLU
AC_SUBST([MLU_CALLER_CPPFLAGS],["$CPPFLAGS $OPENMP_CXXFLAGS"])
AC_SUBST([MLU_CALLER_CXXFLAGS],["$CXXFLAGS $OPENMP_CXXFLAGS"])
AC_SUBST([MLU_CALLER_LDFLAGS],["$LDFLAGS"])
AC_SUBST([MLU_CALLER_LIBS],["$LIBS"])

# Path to hdf5_cpp
AC_ARG_WITH([hdf5],
  [AS_HELP_STRING([--with-hdf5],
    [path to installation location of hdf5])],
  AS_IF([test "x$with_hdf5" != xyes && test "x$with_hdf5" != xno],
          [CPPFLAGS+=" -I$with_hdf5/include"]
          [LDFLAGS+=" -L$with_hdf5/lib"])
  )
AC_CHECK_LIB([hdf5],[H5check_version],[LIBS+=" -lhdf5_cpp -lhdf5"],[AC_MSG_FAILURE([hdf5 library not found])])
AC_MSG_CHECKING([for hdf5_cpp]);
AC_LINK_IFELSE(
        [AC_LANG_SOURCE([[
    #include <stdio.h>
    #include <iostream>
    #include <H5Cpp.h>
    int main(int argc, char *argv[]) {
        H5::DataSpace ds;
        std::cout << "H5::DataSpace ds; ds.isSimple()=" << ( ds.isSimple() ? "true" : "false" ) << std::endl;
        return 0; }
    ]])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [AC_MSG_ERROR([Could not compile a minimal hdf5_cpp program])])

# Check which GSL to use
AC_ARG_WITH([gsl],
  [AS_HELP_STRING([--with-gsl],
    [yes to build GSL sub-package, no assumes GSL can be found in environment,
    otherwise path to installation location of GSL @<:@default=yes@:>@])],
  [],
  [: m4_divert_text([DEFAULTS], [with_gsl=check])]
  )

# Check which GSL to use
AC_ARG_WITH([cblas],
  [AS_HELP_STRING([--with-cblas],
    [When not building GSL, CBLAS libraries to link to, e.g. "-lnvblas -lgslcblas" for NVidia BLAS
    or "-framework Accelerate" on Mac OS @<:@default=-lgslcblas@:>@])],
  [],
  [: m4_divert_text([DEFAULTS], [with_cblas=-lgslcblas])]
  )

# Don't build GSL if not present
AS_IF([test "x$with_gsl" = xcheck || test "x$with_gsl" = xyes],
      [AC_MSG_CHECKING([whether GSL subpackage installed])]
      [AS_IF([test -f "$srcdir/gsl-2.7/configure.ac"],
        [AC_MSG_RESULT([yes])]
        [with_gsl=yes]
        [AC_MSG_NOTICE([Making gsl from subpackage in gsl-2.7 directory])],
        [AC_MSG_RESULT([no])]
        [AS_IF([test "x$with_gsl" = xyes],
          [AC_MSG_ERROR([gsl-2.7/configure.ac missing. Run bootstrap.sh or specify --with-gsl=path])],
          [with_gsl=no]
        )]
      )]
)

AC_MSG_NOTICE([After GSL subpackage check with_gsl=$with_gsl])

# Call GSL's ./configure script recursively.
AM_CONDITIONAL([WANT_BUILD_GSL], [test "x$with_gsl" = xyes])
AM_CONDITIONAL([WANT_LINK_GSL], [test "x$with_gsl" != xyes && test "x$with_gsl" != xno])
AS_IF([test "x$with_gsl" = xyes],
        [AC_CONFIG_SUBDIRS([gsl-2.7])],
      [AS_IF([test "x$with_gsl" != xno],
        [CPPFLAGS+=" -I$with_gsl/include"]
        [LDFLAGS+=" -L$with_gsl/lib"]
      )]
      AC_MSG_NOTICE([Linking to CBLAS libraries with_cblas=$with_cblas])
      AS_IF([test "x$with_cblas" = "x-lgslcblas"],
        [AC_CHECK_LIB([m],[cos],,[AC_MSG_FAILURE([math library (-lm) not found])])]
        [AC_CHECK_LIB([gslcblas],[cblas_dgemm],,[AC_MSG_FAILURE([gslcblas library not found])])]
        [AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_FAILURE([gsl library not found])])],
        LIBS+=" $with_cblas -lgsl"
      )
)

# Check whether to use Minuit2
AC_ARG_WITH([minuit2],
  [AS_HELP_STRING([--with-minuit2],
    [yes assumes Minuit2 can be found in environment, no disable Minuit2,
     otherwise path to installation location of Minuit2 @<:@default=no@:>@])],
  [],
  [: m4_divert_text([DEFAULTS], [with_minuit2=no])])

# If Minuit2 is requested, see whether we can build and link to it (because it's pure C++)
AM_CONDITIONAL([WANT_LINK_MINUIT2], [test "x$with_minuit2" != xno])
AS_IF([test "x$with_minuit2" != xno],
      [AS_IF([test "x$with_minuit2" != xyes],
              [CPPFLAGS+=" -I$with_minuit2/include"]
              [LDFLAGS+=" -L$with_minuit2/lib"])]
      [LIBS_Backup="$LIBS"]
      [LIBS+=" -lMinuit2"]
      [AC_MSG_CHECKING([whether a minimal Minuit2 program compiles])]
      [AC_LINK_IFELSE(
              [AC_LANG_SOURCE([[
          #include <Minuit2/Minuit2Minimizer.h>
          int main(int argc, char *argv[]) {
              ROOT::Minuit2::MnUserParameterState Minuit2State;
              return 0; }
          ]])],
          [AC_MSG_RESULT([yes])]
          [AC_DEFINE([HAVE_MINUIT2], [], [Optional Minuit2 fitter])],
          [AC_MSG_RESULT([no])]
          [AC_MSG_ERROR([Could not compile a minimal Minuit2 program])])]
      [LIBS="$LIBS_Backup"]
      )

#AC_DEFINE_UNQUOTED([MLU_SHA], ["`git rev-parse HEAD`"], [git commit SHA])
#AC_DEFINE_UNQUOTED([MLU_SHORT_SHA], ["`git rev-parse --short HEAD`"], [git commit short SHA])
#AC_DEFINE_UNQUOTED([MLU_BRANCH], ["`git rev-parse --abbrev-ref HEAD`"], [git branch])
AC_DEFINE_UNQUOTED([MLU_GIT_SUMMARY], ["MLU `git rev-parse --abbrev-ref HEAD` `git rev-parse HEAD`"], [git branch and commit for MLU])

# These are the configure options needed to build MLU
AC_SUBST([MLU_CPPFLAGS],["$CPPFLAGS"])
AC_SUBST([MLU_CXXFLAGS],["$CXXFLAGS"])
AC_SUBST([MLU_LDFLAGS],["$LDFLAGS"])
AC_SUBST([MLU_LIBS],["$LIBS"])
AC_SUBST([MLU_GIT_SUMMARY])
AC_CONFIG_HEADERS([MLUconfig.h],
    [$SED_INPLACE -e 's|#define PACKAGE|#define MLU_PACKAGE|g' -e 's|#define VERSION|#define MLU_VERSION|g' MLUconfig.h],
    [if test x"$host_os" = x"${host_os#darwin}" ; then]
        [SED_INPLACE="sed -i"]
    [else]
        [SED_INPLACE="sed -i .bak"]
    [fi])
AC_CONFIG_FILES([Makefile MLU/Makefile Analyse/Makefile Fit/Makefile])
AC_OUTPUT
