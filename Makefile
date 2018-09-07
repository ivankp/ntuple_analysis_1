SHELL := bash
CXX := g++
CPPFLAGS := -std=c++14 -Iinclude
# CXXFLAGS := -Wall -O3 -flto -fmax-errors=3 $(CPPFLAGS)
CXXFLAGS := -Wall -g -fmax-errors=3 $(CPPFLAGS)
LDFLAGS :=
LDLIBS :=

BLD := .build
EXT := .cc

.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean)))

LDFLAGS += $(shell sed -r 's/([^:]+)(:|$$)/ -L\1/g' <<< "$$LIBRARY_PATH")

ROOT_CXXFLAGS := $(shell root-config --cflags | sed 's/ -std=c++[^ ]\+//')
ROOT_LDFLAGS  := $(shell root-config --ldflags)
ROOT_LDLIBS   := $(shell root-config --libs)

FJ_PREFIX   := $(shell fastjet-config --prefix)
FJ_CXXFLAGS := -I$(FJ_PREFIX)/include
FJ_LDLIBS   := -L$(FJ_PREFIX)/lib -Wl,-rpath=$(FJ_PREFIX)/lib -lfastjet

LHAPDF_PREFIX   := $(shell lhapdf-config --prefix)
LHAPDF_CXXFLAGS := $(shell lhapdf-config --cppflags)
LHAPDF_LDLIBS   := $(shell lhapdf-config --libs) -Wl,-rpath=$(LHAPDF_PREFIX)/lib

# RPATH
rpath_script := ldd $(shell root-config --libdir)/libTreePlayer.so \
  | sed -nr 's|.*=> (.+)/.+\.so[.0-9]* \(.*|\1|p' \
  | sort -u \
  | sed -nr '/^(\/usr)?\/lib/!s/^/-Wl,-rpath=/p'
ROOT_LDLIBS += $(shell $(rpath_script))

SRCS := $(shell find -L src -type f -name '*$(EXT)')
DEPS := $(patsubst src/%$(EXT),$(BLD)/%.d,$(SRCS))

GREP_EXES := grep -rl '^ *int \+main *(' src --include='*$(EXT)'
EXES := $(patsubst src%$(EXT),bin%, \
  $(shell $(GREP_EXES)) \
  $(wildcard src/analyses/*$(EXT)))

all: $(EXES)

-include $(DEPS)

# -------------------------------------------------------------------
C_Higgs2diphoton := $(ROOT_CXXFLAGS)

C_reweighter := $(ROOT_CXXFLAGS) $(LHAPDF_CXXFLAGS)

L_merge := -lboost_iostreams

bin/merge: \
  $(BLD)/ivanp/program_options/program_options.o
# -------------------------------------------------------------------

$(DEPS): $(BLD)/%.d: src/%$(EXT)
	@mkdir -pv $(dir $@)
	$(CXX) $(CPPFLAGS) $(C_$*) -MM -MT '$(@:.d=.o)' $< -MF $@

$(BLD)/%.o:
	@mkdir -pv $(dir $@)
	$(CXX) $(CXXFLAGS) $(C_$*) -c $(filter %$(EXT),$^) -o $@

bin/%: $(BLD)/%.o
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

$(BLD)/analyses/%.o:
	@mkdir -pv $(dir $@)
	$(CXX) $(CXXFLAGS) $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS) -c $(filter %$(EXT),$^) -o $@

bin/analyses/%: $(BLD)/analyses/%.o \
  $(BLD)/ivanp/program_options/program_options.o \
  $(BLD)/ivanp/binner/re_axes.o \
  $(BLD)/glob.o \
  $(BLD)/copy_file.o \
  $(BLD)/Higgs2diphoton.o \
  $(BLD)/reweighter.o
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS) $(ROOT_LDLIBS) -lTreePlayer $(FJ_LDLIBS) $(LHAPDF_LDLIBS) -lboost_iostreams

endif

clean:
	@rm -rfv $(BLD) bin

