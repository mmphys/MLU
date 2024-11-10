# Mike's Lattice Utilities - MultiFit

AM_CPPFLAGS = @MLU_CALLER_CPPFLAGS@ -I$(srcdir)/..
AM_CXXFLAGS = @MLU_CALLER_CXXFLAGS@
AM_LDFLAGS = @MLU_CALLER_LDFLAGS@ @MLU_CALLER_LIBS@

if WANT_LINK_MINUIT2
AM_LDFLAGS += -lMinuit2
endif

# link with the MLU library.
LDADD = ../MLU/libMLU.la

bin_PROGRAMS = Continuum MultiFit

Continuum_SOURCES = Continuum.cpp Continuum.hpp Covar.cpp Covar.hpp \
  FitGSL.cpp FitGSL.hpp Fitter.cpp Fitter.hpp \
  FitterThread.cpp FitterThread.hpp Model.cpp Model.hpp \
  ModelContinuum.cpp ModelContinuum.hpp
MultiFit_SOURCES = Covar.cpp Covar.hpp FitGSL.cpp FitGSL.hpp \
  Fitter.cpp Fitter.hpp FitterThread.cpp FitterThread.hpp \
  Model.cpp Model.hpp Model2pt.cpp Model2pt.hpp \
  Model3pt.cpp Model3pt.hpp ModelCommon.cpp ModelCommon.hpp \
  ModelConstant.cpp ModelConstant.hpp ModelMisc.cpp ModelMisc.hpp \
  ModelRatio.cpp ModelRatio.hpp ModelZV.cpp ModelZV.hpp \
  MultiFit.cpp MultiFit.hpp
if WANT_LINK_MINUIT2
Continuum_SOURCES += FitMinuit2.cpp FitMinuit2.hpp
MultiFit_SOURCES += FitMinuit2.cpp FitMinuit2.hpp
endif
