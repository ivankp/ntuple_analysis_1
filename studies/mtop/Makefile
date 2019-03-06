SHELL := bash
CXX := g++
CPPFLAGS := -std=c++14 -I. -I../../include -I../../src/analyses
CXXFLAGS := -Wall -O3 -flto -fmax-errors=3 $(CPPFLAGS)
# CXXFLAGS := -Wall -g -fmax-errors=3 $(CPPFLAGS) -DDEBUG_AT
LDFLAGS :=
LDLIBS :=

.PHONY: all clean

EXE := $(shell ls | sed -n 's/\.cc$$//p')
DEPS := $(EXE:%=.%.d)

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean)))

include ../../flags.mk

all: $(EXE)

-include $(DEPS)

# -------------------------------------------------------------------
C_analysis := $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS) $(LDFLAGS)
L_analysis := $(ROOT_LDLIBS) -lTreePlayer $(FJ_LDLIBS) $(LHAPDF_LDLIBS) -lpcre

C_analysis2 := $(ROOT_CXXFLAGS) $(FJ_CXXFLAGS) $(LDFLAGS)
L_analysis2 := $(ROOT_LDLIBS) -lTreePlayer $(FJ_LDLIBS) $(LHAPDF_LDLIBS) -lpcre

VPATH := ../../.build

analysis analysis2: \
  ivanp/program_options/program_options.o \
  ivanp/binner/re_axes.o \
  glob.o \
  Higgs2diphoton.o \
  reweighter.o
# -------------------------------------------------------------------

$(DEPS): .%.d: %.cc
	$(CXX) $(CPPFLAGS) -MM -MT '$(@:.%.d=%) $@' $< -MF $@

$(EXE): %: %.cc
	$(CXX) $(CXXFLAGS) $(C_$*) $(LDFLAGS) $(filter %.cc,$^) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

endif

clean:
	@rm -fv $(EXE) $(EXE:%=.%.d)

