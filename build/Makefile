# This will make each .c* file into separate executable
# PRE-REQUISITE ENVIRONMENT VARIABLES:
#   GridPre         Path to Grid       pre-requisites (custom built)
#   GridPkg         Path to Grid       pre-requisites (from package installs)
#   Grid            Path to Grid       source
#   GridBuild       Path to Grid       build directory
#   Hadrons         Path to Hadrons    source
#   HadronsBuild    Path to Hadrons    build directory
# Not used at present (but could be reinstated)
#   LatAnalyze      Path to LatAnalyze source
#   LatAnalyzeBuild Path to LatAnalyze build directory

GridConf  := $(GridBuild)/grid-config
HadConf   := $(HadronsBuild)/hadrons-config
#LatAnConf:= $(LatAnalyzeBuild)/latan-config
CXX       := $(shell $(GridConf) --cxx)

# De-duplicate header include (-I) paths, followed by other options
CXXFLAGS  := $(shell $(GridConf) --cxxflags) #$(shell $(LatAnConf) --cxxflags)
INCLUDE   := $(sort $(patsubst -I%,%,$(filter -I%,$(CXXFLAGS))))
CXXFLAGS  := $(addprefix -I,$(GridBuild)/Grid $(Grid) $(Hadrons) $(GridPre)/include $(GridPkg)/include $(LatAnalyze) $(Grid)/Eigen $(INCLUDE)) $(filter-out -I%,$(CXXFLAGS))
# Turn off compiler optimisation
CXXFLAGS  := $(filter-out -O%,$(CXXFLAGS))
# Disable warnings for macro redefinition (since I use this in Common.hpp)
ifeq ($(CXX),clang++)
CXXFLAGS  := $(CXXFLAGS) -Wno-macro-redefined
endif

# De-duplicate library (-L) paths and libraries (-l), followed by other options
LDLIBS    := -lHadrons $(shell $(HadConf) --libs) -lMinuit2 -lgsl -lgslcblas
LDFLAGS   := $(shell $(HadConf) --ldflags)
LDPATHS   := $(sort $(patsubst -L%,%,$(filter -L%,$(LDFLAGS))))
LDFLAGS   := $(addprefix -L,$(HadronsBuild)/Hadrons $(GridPre)/lib $(if $(GridPkg),$(GridPkg)/lib,) $(LDPATHS)) $(filter-out -L%,$(LDFLAGS))

# Make a separate executable for every source in subdirectories except build*, run* or script*
# Except filenames starting with "Common", for which explicit target rules must be provided
VPATH     := $(shell find .. -maxdepth 1 -type d \! \( -iname 'build*' -or -iname 'run*' -or -iname 'test*' -or -iname 'script*' -or -iname '.*' -or -iname '*.xcodeproj' \))
srcfiles  := $(filter-out %~,$(shell find $(VPATH) -iname '*.c*'))
basefiles := $(basename $(notdir $(srcfiles)))
objects   := $(addsuffix .o,$(basefiles))
target    := $(filter-out Common%,$(filter-out FitModel,$(basefiles)))

all: $(target)

bootstrap: bootstrap.o Common.o

corrstat: corrstat.o Common.o

Debug: Debug.o Common.o

hlxml: hlxml.o Common.o

FitSummary: FitSummary.o Common.o

fold: fold.o Common.o

hlxml: hlxml.o Common.o

MultiFit: MultiFit.o Common.o FitModel.o

peramb-compare: peramb-compare.o Common.o

repair: repair.o Common.o

mixedop: mixedop.o Common.o

xmlGFWall: xmlGFWall.o Common.o

xmlPeramb: xmlPeramb.o Common.o

xmlZ2: xmlZ2.o Common.o

xml3pt: xml3pt.o Common.o

clean:
	rm -rf $(target) $(objects) *~ *.dSYM

debug:
	echo VPATH="$(VPATH)"
	echo srcfiles="$(srcfiles)"
	echo basefiles="$(basefiles)"
	echo objects="$(objects)"
	echo target="$(target)"
