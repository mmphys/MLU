# Mike's Lattice Utilities - libMLU

AM_CPPFLAGS = -I$(srcdir)/.. @MLU_CPPFLAGS@
AM_CXXFLAGS = @MLU_CXXFLAGS@
AM_LDFLAGS = @MLU_LDFLAGS@ @MLU_LIBS@

if WANT_LINK_MINUIT2
AM_LDFLAGS += -lMinuit2
endif

lib_LTLIBRARIES = libMLU.la

libMLU_la_SOURCES = CorrFile.cpp CorrFile.hpp DataSet.cpp DataSet.hpp \
  DebugInfo.cpp DebugInfo.hpp FitRange.cpp FitRange.hpp FitRangeImp.hpp \
  Fold.cpp Fold.hpp GSLVecMat.cpp GSLVecMat.hpp GSLVecMatImp.hpp \
  HDF5.cpp HDF5.hpp JackBoot.cpp JackBoot.hpp Math.cpp Math.hpp \
  MLU.hpp MLUFirst.hpp Model.cpp Model.hpp Param.cpp Param.hpp \
  Physics.cpp Physics.hpp Posix.cpp Posix.hpp Sample.cpp Sample.hpp \
  String.cpp String.hpp Utility.cpp Utility.hpp

# Don't forget to update nobase_include_HEADERS in ../Makefile.am

if WANT_BUILD_GSL
libMLU_la_LIBADD = ../gsl-2.7/libgsl.la
endif
