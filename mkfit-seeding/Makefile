# Build seedfind against the standalone build's libMkFitCore.so.
#
#   make                  # -> $(BLD)/test-seedgeom/bin/seedfind
#   make BLD=/path/to/standalone
#   make ARCH=-march=native   # replace the build's -mavx, for an A/B only
#
# The -D set and codegen flags come from the build's `make echo-aclic`, NOT
# from a copy here: standalone/Event.h gains members under MKFIT_TRACE, so a
# mismatch silently misreads every Event member.  The shared build is switched
# between the trace build and a ROOT-off build by other work, so the flags are
# re-read on every make, and the binary depends on libMkFitCore.so: a newer
# library forces a rebuild with whatever flags it was built with.  ROOT is
# linked only if libMkFitCore.so itself links libCore.

BLD   ?= $(abspath $(CURDIR)/../..)
OUT   := $(BLD)/test-seedgeom/bin
SRC   := $(abspath $(BLD)/..)
SADIR := $(shell sed -n 's/^export SADIR *:= *//p' $(BLD)/Makefile)

ROOTLIB := $(shell ldd $(BLD)/libMkFitCore.so | sed -n 's|.*libCore.so => \(.*\)/libCore.so.*|\1|p')
ifneq ($(ROOTLIB),)
  ROOTSYS_B := $(abspath $(ROOTLIB)/..)
  ROOT_ENV  := ROOTSYS=$(ROOTSYS_B) PATH=$(ROOTSYS_B)/bin:$$PATH
  ROOT_CF   := $(shell $(ROOTSYS_B)/bin/root-config --cflags)
  ROOT_LIBS := $(shell $(ROOTSYS_B)/bin/root-config --libs) -Wl,-rpath,$(ROOTLIB)
endif

ECHO := $(shell $(ROOT_ENV) make -s --no-print-directory -C $(BLD)/objs-Core -f $(SADIR)/Makefile \
          SADIR=$(SADIR) SRCDIR=$(SRC) BLDDIR=$(BLD) echo-aclic)
# the first word of each line is glued to its key ("ACLIC_DEFS=-DFOO"), so
# split the keys off before filtering
ECHO_W := $(subst ACLIC_DEFS=,ACLIC_DEFS= ,$(subst ACLIC_OPT=,ACLIC_OPT= ,$(ECHO)))
DEFS := $(filter -D%,$(ECHO_W))
OPT  := $(filter-out -D% ACLIC_DEFS= ACLIC_OPT=,$(ECHO_W))
ifdef ARCH
  OPT := $(filter-out -mavx -mavx2 -march=%,$(OPT)) $(ARCH)
endif

CXXFLAGS := -std=c++20 $(OPT) $(DEFS) $(ROOT_CF) -I$(SRC) -I$(BLD)/mkFit-external -I$(CURDIR)
LDFLAGS  := -L$(BLD) -lMkFitCore -Wl,-rpath,$(BLD) $(ROOT_LIBS) -ltbb

HDRS := $(wildcard $(CURDIR)/*.h)

$(OUT)/seedfind: $(CURDIR)/seedfind.cc $(HDRS) $(BLD)/libMkFitCore.so | $(OUT)
	@echo "[seedfind] DEFS: $(DEFS)"
	@echo "[seedfind] OPT : $(OPT)"
	@echo "[seedfind] ROOT: $(if $(ROOTLIB),$(ROOTLIB),none)"
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

$(OUT):
	mkdir -p $@

.PHONY: flags
flags:
	@echo "DEFS=$(DEFS)"
	@echo "OPT=$(OPT)"
	@echo "ROOTLIB=$(ROOTLIB)"
