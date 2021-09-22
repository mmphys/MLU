# Mike's Lattice Utilities - Analyze

LDADD = ../MLU/libMLU.a $(MLU_LIBS)

bin_PROGRAMS = bootstrap corr corrstat Debug Divide FitSummary hlxml MultiFit repair

bootstrap_SOURCES = bootstrap.cpp
corr_SOURCES = corr.cpp
corrstat_SOURCES = corrstat.cpp
Debug_SOURCES = Debug.cpp
Divide_SOURCES = Divide.cpp
FitSummary_SOURCES =  FitSummary.cpp
hlxml_SOURCES = hlxml.cpp
MultiFit_SOURCES = FitGSL.cpp FitModel.cpp MultiFit.cpp
if MINUIT2
MultiFit_SOURCES += FitMinuit2.cpp
endif
repair_SOURCES = repair.cpp
