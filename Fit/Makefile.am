# Mike's Lattice Utilities - MultiFit

LDADD = ../MLU/libMLU.a $(MLU_LIBS)

bin_PROGRAMS = Continuum MultiFit

Continuum_SOURCES = Continuum.cpp Covar.cpp FitGSL.cpp Fitter.cpp FitterThread.cpp Model.cpp ModelContinuum.cpp
MultiFit_SOURCES = Covar.cpp FitGSL.cpp Fitter.cpp FitterThread.cpp Model.cpp Model2pt.cpp Model3pt.cpp ModelCommon.cpp ModelConstant.cpp ModelMisc.cpp ModelRatio.cpp ModelZV.cpp MultiFit.cpp
if MINUIT2
Continuum_SOURCES += FitMinuit2.cpp
MultiFit_SOURCES += FitMinuit2.cpp
endif
