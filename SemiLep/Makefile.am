# SemiLep

# Include MLU when compiling this directory
AM_CPPFLAGS = -I$(srcdir)/../MLU @AM_CPPFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@
if WANT_BUILD_MLU
LDADD = ../MLU/MLU/libMLU.la
else
AM_LDFLAGS += -lMLU
endif

bin_PROGRAMS = Debug GaugeCmp hlxml peramb-compare TensorDemo xml3pt xmlGFWall xmlPeramb xmlZ2
# lime_truncate

Debug_SOURCES = Debug.cpp
GaugeCmp_SOURCES = GaugeCmp.cpp
hlxml_SOURCES = hlxml.cpp
# lime_truncate_SOURCES = lime_truncate.c
peramb_compare_SOURCES = peramb-compare.cpp
TensorDemo_SOURCES = TensorDemo.cpp
xml3pt_SOURCES = xml3pt.cpp HadronsApp.cpp
xmlGFWall_SOURCES = xmlGFWall.cpp
xmlPeramb_SOURCES = xmlPeramb.cpp
xmlZ2_SOURCES = xmlZ2.cpp
