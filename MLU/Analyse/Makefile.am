# Mike's Lattice Utilities - Analyze

AM_CPPFLAGS = -I$(srcdir)/.. @MLU_CALLER_CPPFLAGS@
AM_CXXFLAGS = @MLU_CALLER_CXXFLAGS@
AM_LDFLAGS = @MLU_CALLER_LDFLAGS@ @MLU_CALLER_LIBS@

# link with the MLU library.
LDADD = ../MLU/libMLU.la

bin_PROGRAMS = bootstrap corr corrstat CRatio Divide FitSummary \
  GetColumn ImportCorr MakeEnsembleInfo mixedop repair SynthData

bootstrap_SOURCES = bootstrap.cpp bootstrap.hpp
corr_SOURCES = corr.cpp
corrstat_SOURCES = corrstat.cpp
CRatio_SOURCES = CRatio.cpp CRatio.hpp
Divide_SOURCES = Divide.cpp Divide.hpp
FitSummary_SOURCES =  FitSummary.cpp FitSummary.hpp
GetColumn_SOURCES = GetColumn.cpp
ImportCorr_SOURCES = ImportCorr.cpp ImportCorr.hpp
MakeEnsembleInfo_SOURCES = MakeEnsembleInfo.cpp MakeEnsembleInfo.
mixedop_SOURCES = mixedop.cpp
repair_SOURCES = repair.cpp
SynthData_SOURCES = SynthData.cpp SynthData.hpp
