# Mike's Lattice Utilities - MultiFit

LDADD = ../MLU/libMLU.a $(MLU_LIBS)

bin_PROGRAMS = MultiFit

MultiFit_SOURCES = FitGSL.cpp FitModel.cpp MultiFit.cpp
if MINUIT2
MultiFit_SOURCES += FitMinuit2.cpp
endif
