# This will make each .c* file into separate executable
# PRE-REQUISITE ENVIRONMENT VARIABLES:
#   Grid            Path to Grid       source
#   GridBuild       Path to Grid       build directory
#   LatAnalyze      Path to LatAnalyze source
#   LatAnalyzeBuild Path to LatAnalyze build directory

GridConf  := $(GridBuild)/grid-config
LatAnConf := $(LatAnalyzeBuild)/latan-config
CXX       := $(shell $(GridConf) --cxx)
CC        ?= $(CXX)

# De-duplicate header include (-I) paths, followed by other options
CXXFLAGS  := $(shell $(GridConf) --cxxflags) $(shell $(LatAnConf) --cxxflags)
INCLUDE   := $(sort $(patsubst -I%,%,$(filter -I%,$(CXXFLAGS))))
CXXFLAGS  := $(addprefix -I,$(GridBuild)/Grid $(Grid) $(LatAnalyzeBuild) $(LatAnalyze)/lib $(INCLUDE)) $(filter-out -I%,$(CXXFLAGS))
# Turn off compiler optimisation
CXXFLAGS  := $(filter-out -O%,$(CXXFLAGS))

# De-duplicate library (-L) paths and libraries (-l), followed by other options
LDLIBS    := -lHadrons -lGrid -lLatAnalyze -lLexers $(sort $(shell $(GridConf) --libs) $(shell $(LatAnConf) --libs))
LDFLAGS   := $(shell $(GridConf) --ldflags) $(shell $(LatAnConf) --ldflags)
LDPATHS   := $(sort $(patsubst -L%,%,$(filter -L%,$(LDFLAGS))))
LDFLAGS   := $(addprefix -L,$(GridBuild)/Grid $(GridBuild)/Hadrons $(LatAnalyzeBuild)/lib/.libs $(GridPre)/lib $(LDPATHS)) $(filter-out -L%,$(LDFLAGS))

srcdir    := ../Analyze
VPATH     := $(srcdir)
srcfiles  := $(filter-out $(srcdir)/Common.cpp,$(filter-out $(srcdir)/CommonLatAn.cpp,$(filter-out %~,$(wildcard $(srcdir)/*.c*))))
basefiles := $(basename $(notdir $(srcfiles)))
objects   := $(addsuffix .o,$(basefiles)) Common.o
target    := $(basefiles)

all: $(target)

MultiFit: MultiFit.o Common.o CommonLatAn.o

bootstrap: bootstrap.o Common.o

hlxml: hlxml.o Common.o

clean:
	rm -rf $(target) $(objects) $(depfile) *~ $(addsuffix .dSYM,$(basefiles))

debug:
	echo srcfiles="$(srcfiles)"
	echo basefiles="$(basefiles)"
	echo objects="$(objects)"
	echo target="$(target)"
