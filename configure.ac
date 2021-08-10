## Process this file with autoconf to produce a configure script.
# Meson Lattice Utilities

AC_PREREQ([2.69])
AC_INIT([MLU], [0.1])
AC_COPYRIGHT([Michael Marshall 2021])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR([MLU/CommonGrid.cpp])
#AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([MLUconfig.h],
    [$SED_INPLACE -e 's|#define PACKAGE|#define MLU_PACKAGE|g' -e 's|#define VERSION|#define MLU_VERSION|g' MLUconfig.h],
    [if test x"$host_os" == x"${host_os#darwin}" ; then]
        [SED_INPLACE="sed -i"]
    [else]
        [SED_INPLACE="sed -i .bak"]
    [fi])
AM_INIT_AUTOMAKE([-Wall -Werror foreign 1.13])
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
    [if test x"$with_grid" == x ; then]
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
    [if test x"$with_hadrons" == x ; then]
        [AC_MSG_ERROR([hadrons-config --prefix returns empty string])]
    [fi]
    [AC_MSG_RESULT([$with_hadrons])])

AC_ARG_ENABLE([check], AS_HELP_STRING([--disable-check], [Disable check of Grid and Hadrons build]))

# Optional arguments for other pre-requisite libraries

GSLConfig=gsl-config
GSLInstall=no
AC_ARG_WITH([gsl],
    [AS_HELP_STRING([--with-gsl=prefix],
    [Install path of the GNU Scientific Library])],
    [AC_MSG_NOTICE([--with-gsl=$with_gsl])]
    [GSLConfig="$with_gsl/bin/$GSLConfig"],
    [AC_CHECK_PROG([GSLCONF],[${GSLConfig}],[yes],[no])]
    [if test x"$GSLCONF" != xyes ; then]
        [GSLInstall=yes]
    [fi])

AC_ARG_WITH([minuit2],
    [AS_HELP_STRING([--with-minuit2=prefix],
    [try this for a non-standard install prefix of Minuit2])],
    [CPPFLAGS="-I$with_minuit2/include $CPPFLAGS"]
    [LDFLAGS="-L$with_minuit2/lib $LDFLAGS"])

# If the user hasn't specified a CXX compiler, get it from Hadrons

if test x"$CXX" == x ; then
    CXX="`$HadronsConfig --cxx`"
    if test x"$CXX" != x ; then
        AC_MSG_NOTICE([hadrons-config --cxx=$CXX])
    else
        AC_MSG_WARN([hadrons-config --cxx unknown])
    fi
fi
if test x"$CXXLD" == x ; then
    CXXLD="`$HadronsConfig --cxxld`"
    if test x"$CXXLD" != x ; then
        AC_MSG_NOTICE([hadrons-config --cxxld=$CXXLD])
    else
        AC_MSG_WARN([hadrons-config --cxxld unknown])
    fi
fi

# Now setup my environment
CPPFLAGS="-I$with_grid/include $CPPFLAGS"
LDFLAGS="-L$with_grid/lib $LDFLAGS"
if test x"$with_hadrons" != x"$with_grid" ; then
  CPPFLAGS="-I$with_hadrons/include $CPPFLAGS"
  LDFLAGS="-L$with_hadrons/lib $LDFLAGS"
fi
CXXFLAGS="`$HadronsConfig --cxxflags` $CXXFLAGS"
LDFLAGS="`$HadronsConfig --ldflags` $LDFLAGS"
LIBS="-lHadrons `$HadronsConfig --libs` $LIBS"

# I'd like to use c++17
#CXXFLAGS="${CXXFLAGS//-std=c++* /-std=c++17 }"

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_AR
AC_LANG([C++])

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

AbsMLUBuildDir=`pwd -P`; cd ${srcdir}; abs_srcdir=`pwd -P`; cd ${AbsMLUBuildDir}
MLUBuildDir=${AbsMLUBuildDir##*/}
AC_MSG_NOTICE([AbsMLUBuildDir=$AbsMLUBuildDir])
AC_MSG_NOTICE([MLUBuildDir=$MLUBuildDir])

# Check for Minuit2 by making a small test program (can't check libraries as pure C++)

LIBS+=" -lMinuit2"
GotMinuit2=yes
AC_MSG_CHECKING([that a minimal Minuit2 program compiles]);
AC_LINK_IFELSE(
        [AC_LANG_SOURCE([[
    #include <Minuit2/Minuit2Minimizer.h>
    int main(int argc, char *argv[]) {
        ROOT::Minuit2::MnUserParameterState Minuit2State;
        return 0; }
    ]])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [if test x"$with_minuit2" != x ; then]
        [AC_MSG_ERROR([Could not compile a minimal Minuit2 program])]
    [fi]
    [GotMinuit2=no])
if test x$GotMinuit2 != xyes ; then
    # See whether I've already installed a version in prefix
    if test -d $prefix/include/Minuit2 ; then
        AC_MSG_NOTICE([Using Minuit2 previously installed in --prefix])
    else
        # Install my version in prefix
        Minuit2Dir="$abs_srcdir/.Package/Minuit2-5.34.14"
        Minuit2ConfigureDir="."
        if test x$AbsMLUBuildDir != x$abs_srcdir ; then
            Minuit2Dir="$Minuit2Dir/$MLUBuildDir"
            Minuit2ConfigureDir=".."
        fi
        AC_MSG_NOTICE([Making $Minuit2Dir])
        mkdir -p "$Minuit2Dir"
        cd "$Minuit2Dir"
        OldCXX="$CXX"
        if test x${OldCXX:0:4} == xnvcc ; then
            OldCC="$CC"
            OldCPPFLAGS="$CPPFLAGS"
            OldCXXDEPMODE="$CXXDEPMODE"
            OldCXXFLAGS="$CXXFLAGS"
            OldCXXLD="$CXXLD"
            CXX="g++"
            CC=gcc
            CPPFLAGS=""
            CXXDEPMODE="depmode=gcc3"
            CXXFLAGS=""
            unset CXXLD
            for f in $OldCXXFlags ; do
                if test x${f:0:2} == x-I ; then
                    CXXFLAGS="$CXXFLAGS $f"
                fi
            done
        fi
        $Minuit2ConfigureDir/configure --prefix="$prefix"
        make -j 12
        make install
        if test x${OldCXX:0:4} == xnvcc ; then
            CXX="$OldCXX"
            CC="$OldCC"
            CPPFLAGS="$OldCPPFLAGS"
            CXXDEPMODE="$OldCXXDEPMODE"
            CXXFLAGS="$OldCXXFLAGS"
            CXXLD="$OldCXXLD"
        fi
        cd $AbsMLUBuildDir
        if ! test -d $prefix/include/Minuit2 ; then
            AC_MSG_ERROR([Could not install ${Minuit2Config}])
        fi
    fi
    # Use the version I installed in prefix
    CPPFLAGS="-I$prefix/include $CPPFLAGS"
    LDFLAGS="-L$prefix/lib $LDFLAGS"
    AC_MSG_CHECKING([that Minuit2 installed ok]);
    AC_LINK_IFELSE(
            [AC_LANG_SOURCE([[
        #include <Minuit2/Minuit2Minimizer.h>
        int main(int argc, char *argv[]) {
            ROOT::Minuit2::MnUserParameterState Minuit2State;
            return 0; }
        ]])],
        [AC_MSG_RESULT([yes])]
        [GotMinuit2=yes],
        [AC_MSG_RESULT([no])]
        [AC_MSG_ERROR([Could not install Minuit2])])
fi

# Check for GSL

if test x$GSLInstall == xyes ; then
    GSLConfig="$prefix/bin/$GSLConfig"
    if test -x "$GSLConfig" ; then
        AC_MSG_NOTICE([Using GSL previously installed in --prefix])
    else
        GSLDir="$abs_srcdir/.Package/gsl-2.7"
        GSLConfigureDir="."
        if test x$AbsMLUBuildDir != x$abs_srcdir ; then
            GSLDir="$GSLDir/$MLUBuildDir"
            GSLConfigureDir=".."
        fi
        AC_MSG_NOTICE([Making $GSLDir])
        mkdir -p "$GSLDir"
        cd "$GSLDir"
        $GSLConfigureDir/configure --prefix="$prefix"
        make -j 12
        make install
        cd $AbsMLUBuildDir
        if ! test -x "$GSLConfig" ; then
            AC_MSG_ERROR([Could not install ${GSLConfig}])
        fi
    fi
fi
LDFLAGS+=" `$GSLConfig --libs`"
CPPFLAGS+=" `$GSLConfig --cflags`"

#AC_SEARCH_LIBS([gsl_set_error_handler], [gsl], [],
#               [AC_MSG_ERROR([GNU Scientific Library not found])])
#if test x"$ac_cv_search_gsl_set_error_handler" != x"none required" ; then
#    LIBS+=" $ac_cv_search_gsl_set_error_handler"
#fi

CPPFLAGS="-I${abs_srcdir} $CPPFLAGS -I$AbsMLUBuildDir"

#AC_DEFINE_UNQUOTED([MLU_SHA], ["`git rev-parse HEAD`"], [git commit SHA])
#AC_DEFINE_UNQUOTED([MLU_SHORT_SHA], ["`git rev-parse --short HEAD`"], [git commit short SHA])
#AC_DEFINE_UNQUOTED([MLU_BRANCH], ["`git rev-parse --abbrev-ref HEAD`"], [git branch])
AC_DEFINE_UNQUOTED([MLU_GIT_SUMMARY], ["MLU `git rev-parse --abbrev-ref HEAD` `git rev-parse HEAD`, Grid `$GridConfig --git`, Hadrons `$HadronsConfig --git`"], [git branch and commit for MLU, Grid and Hadrons])

AM_CPPFLAGS="$CPPFLAGS"
AM_CXXFLAGS="$CXXFLAGS"
AM_CFLAGS="$CFLAGS"
AM_LDFLAGS="$LDFLAGS"
MLU_LIBS="$LIBS"

unset CPPFLAGS
unset CXXFLAGS
unset CFLAGS
unset LDFLAGS
unset LIBS

#AC_SUBST([MLU_SHORT])
#AC_SUBST([MLU_SHA])
#AC_SUBST([MLU_BRANCH])
AC_SUBST([MLU_GIT_SUMMARY])
AC_SUBST([CXXLD])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([MLU_LIBS])
AC_CONFIG_FILES([Makefile MLU/Makefile Analyze/Makefile Study/Makefile])
#AC_CONFIG_SUBDIRS([.Package/gsl-2.7 .Package/Minuit2-5.34.14])
AC_OUTPUT
