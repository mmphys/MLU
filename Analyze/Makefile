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
CXXFLAGS  := $(addprefix -I,$(GridBuild)/Grid $(Grid) $(LatAnalyzeBuild) $(LatAnalyze)/lib) $(shell $(GridConf) --cxxflags) $(shell $(LatAnConf) --cxxflags)
# Turn off compiler optimisation
CXXFLAGS  := $(filter-out -O%,$(CXXFLAGS))

LDLIBS    := -lGrid -lHadrons -lLatAnalyze -lLexers $(shell $(GridConf) --libs) $(shell $(LatAnConf) --libs)
LDFLAGS   := $(addprefix -L,$(GridBuild)/Grid $(GridBuild)/Hadrons $(LatAnalyzeBuild)/lib/.libs) $(shell $(GridConf) --ldflags) $(shell $(LatAnConf) --ldflags)

srcfiles  := $(filter-out %~,$(wildcard *.c*))
base      := $(basename $(srcfiles))
objects   := $(addsuffix .o,$(base))
#target    := $(addprefix build/,$(base))
target    := $(base)

all: $(target)
	mkdir -p build

clean:
	rm -rf $(target) $(objects) $(depfile) *~ $(addsuffix .dSYM,$(base))
