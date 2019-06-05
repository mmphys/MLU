# This will make each .c* file into separate executable
# PRE-REQUISITE ENVIRONMENT VARIABLES:
#   Grid            Path to Grid       source
#   GridBuild       Path to Grid       build directory
#   LatAnalyze      Path to LatAnalyze source
#   LatAnalyzeBuild Path to LatAnalyze build directory

GridConf  := $(GridBuild)/grid-config
LatAnConf := $(LatAnalyzeBuild)/latan-config
CXX       := $(shell $(GridConf) --cxx)
CC        := $(CXX)
CXXFLAGS  := $(shell $(GridConf) --cxxflags) $(shell $(LatAnConf) --cxxflags)
# Turn off compiler optimisation
CXXFLAGS  := $(filter-out -O%,$(CXXFLAGS))
INCLUDE   := $(sort $(patsubst -I%,%,$(filter -I%,$(CXXFLAGS))))
# Now add deduplicated list of compiler flags
CXXFLAGS  := $(addprefix -I,$(GridBuild)/Grid $(Grid) $(LatAnalyzeBuild) $(LatAnalyze)/lib $(INCLUDE)) $(filter-out -I%,$(CXXFLAGS))

LDLIBS    := $(sort -lGrid -lHadrons -lLatAnalyze -lLexers $(shell $(GridConf) --libs) $(shell $(LatAnConf) --libs))
LDFLAGS   := $(shell $(GridConf) --ldflags) $(shell $(LatAnConf) --ldflags)
LDPATHS   := $(sort $(patsubst -L%,%,$(filter -L%,$(LDFLAGS))))
LDFLAGS   := $(addprefix -L,$(GridBuild)/Grid $(GridBuild)/Hadrons $(LatAnalyzeBuild)/lib/.libs $(GridPre)/lib $(LDPATHS)) $(filter-out -L%,$(LDFLAGS))

srcfiles  := $(filter-out %~,$(wildcard *.c*))
base      := $(basename $(srcfiles))
objects   := $(addsuffix .o,$(base))
#target    := $(addprefix build/,$(base))
target    := $(base)

all: $(target)
	mkdir -p build

clean:
	rm -rf $(target) $(objects) $(depfile) *~ $(addsuffix .dSYM,$(base))
