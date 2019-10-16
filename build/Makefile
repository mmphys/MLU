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

# Make a separate executable for every source in subdirectories except build*, run* or script*
# Except filenames starting with "Common", for which explicit target rules must be provided
VPATH     := $(shell find .. -maxdepth 1 -type d \! \( -iname 'build*' -or -iname 'run*' -or -iname 'script*' -or -iname '*.*' \))
srcfiles  := $(filter-out %~,$(shell find $(VPATH) -iname '*.c*'))
basefiles := $(basename $(notdir $(srcfiles)))
objects   := $(addsuffix .o,$(basefiles))
target    := $(filter-out Common%,$(basefiles))

all: $(target)

MultiFit: MultiFit.o Common.o CommonLatAn.o

bootstrap: bootstrap.o Common.o CommonLatAn.o CommonGrid.o

hlxml: hlxml.o Common.o

clean:
	rm -rf $(target) $(objects) *~ *.dSYM

debug:
	echo VPATH="$(VPATH)"
	echo srcfiles="$(srcfiles)"
	echo basefiles="$(basefiles)"
	echo objects="$(objects)"
	echo target="$(target)"
