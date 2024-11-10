# Mike's Lattice Utilities - Analyze

AM_CPPFLAGS = @MLU_CALLER_CPPFLAGS@ -I$(srcdir)/..
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
# Debug_SOURCES = Debug.cpp
Divide_SOURCES = Divide.cpp Divide.hpp
FitSummary_SOURCES =  FitSummary.cpp FitSummary.hpp
#GaugeCmp_SOURCES = GaugeCmp.cpp
GetColumn_SOURCES = GetColumn.cpp
#hlxml_SOURCES = hlxml.cpp
ImportCorr_SOURCES = ImportCorr.cpp ImportCorr.hpp
#lime_truncate = lime_truncate.c
MakeEnsembleInfo_SOURCES = MakeEnsembleInfo.cpp MakeEnsembleInfo.
mixedop_SOURCES = mixedop.cpp
repair_SOURCES = repair.cpp
SynthData_SOURCES = SynthData.cpp SynthData.hpp
