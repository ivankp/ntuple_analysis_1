SHELL := bash
CXX := g++
CPPFLAGS := -std=c++14 -Iinclude
CXXFLAGS := -Wall -O3 -flto -fmax-errors=3 $(CPPFLAGS)
# CXXFLAGS := -Wall -g -fmax-errors=3 $(CPPFLAGS) -DDEBUG_AT
LDFLAGS :=
LDLIBS :=

BLD := .build
EXT := .cc

.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean)))

include flags.mk

SRCS := $(shell find -L src -type f -name '*$(EXT)')
DEPS := $(patsubst src/%$(EXT),$(BLD)/%.d,$(SRCS))

GREP_EXES := grep -rl '^ *int \+main *(' src --include='*$(EXT)'
EXES := $(patsubst src%$(EXT),bin%, \
  $(shell $(GREP_EXES)) \
  $(wildcard src/analyses/*$(EXT)))

# STUDIES_CMD := find studies -regex '.*/\(.+\)/\1\$(EXT)' | sed 's:\$(EXT)$$::'
# STUDIES := $(shell $(STUDIES_CMD))

all: $(EXES)

-include $(DEPS)

# -------------------------------------------------------------------
C_Higgs2diphoton := $(ROOT_CXXFLAGS)

C_reweighter := $(ROOT_CXXFLAGS) $(LHAPDF_CXXFLAGS)

L_merge_json := -lboost_iostreams -lboost_regex
L_merge := -lboost_iostreams -lboost_regex

C_merge_root := $(ROOT_CXXFLAGS)
L_merge_root := $(ROOT_LDLIBS)

C_check_tree := $(ROOT_CXXFLAGS)
L_check_tree := $(ROOT_LDLIBS) -lTreePlayer

bin/merge bin/merge_json: \
  $(BLD)/ivanp/program_options/program_options.o

bin/merge: \
  $(BLD)/ivanp/io/mem_file.o \
  $(BLD)/ivanp/scribe.o

bin/read_hist bin/list_hists: \
  $(BLD)/ivanp/io/mem_file.o \
  $(BLD)/ivanp/scribe.o

C_analyses/hist_Hjets := -DLOOPSIM
bin/analyses/hist_Hjets: \
  $(BLD)/loopsim/LoopSim.o \
  $(BLD)/loopsim/TreeLevel.o \
  $(BLD)/loopsim/Event.o \
  $(BLD)/loopsim/Flavour.o \
  $(BLD)/loopsim/FlavourPlugin.o
# -------------------------------------------------------------------

$(DEPS): $(BLD)/%.d: src/%$(EXT)
	@mkdir -pv $(dir $@)
	$(CXX) $(CPPFLAGS) $(C_$*) -MM -MT '$(@:.d=.o) $@' $< -MF $@

$(BLD)/%.o:
	@mkdir -pv $(dir $@)
	$(CXX) $(CXXFLAGS) $(C_$*) -c $(filter %$(EXT),$^) -o $@

bin/%: $(BLD)/%.o
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

$(BLD)/analyses/%.o:
	@mkdir -pv $(dir $@)
	$(CXX) $(CXXFLAGS) $(C_analyses/$*) $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS) -c $(filter %$(EXT),$^) -o $@

bin/analyses/%: $(BLD)/analyses/%.o \
  $(BLD)/ivanp/program_options/program_options.o \
  $(BLD)/ivanp/binner/re_axes.o \
  $(BLD)/glob.o \
  $(BLD)/ivanp/io/mem_file.o \
  $(BLD)/ivanp/scribe.o \
  $(BLD)/Higgs2diphoton.o \
  $(BLD)/reweighter.o
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS) $(ROOT_LDLIBS) -lTreePlayer $(FJ_LDLIBS) $(LHAPDF_LDLIBS) -lboost_iostreams $(L_analyses/$*)

endif

clean:
	@rm -rfv $(BLD) bin

