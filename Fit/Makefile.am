# Mike's Lattice Utilities - MultiFit

LDADD = ../MLU/libMLU.a $(MLU_LIBS)

bin_PROGRAMS = MultiFit

MultiFit_SOURCES = Covar.cpp FitGSL.cpp Fitter.cpp FitterThread.cpp ModelBase.cpp MultiFit.cpp Param.cpp
if MINUIT2
MultiFit_SOURCES += FitMinuit2.cpp
endif
