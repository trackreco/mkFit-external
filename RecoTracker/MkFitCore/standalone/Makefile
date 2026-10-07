include ${SADIR}/Makefile.config

CPPFLAGS := -I. -I${SRCDIR} -I${BLDDIR}/mkFit-external ${CPPFLAGS}

CORE_DIR := ${SRCDIR}/RecoTracker/MkFitCore

LIB_CORE := ../libMkFitCore.so

TGTS := ${LIB_CORE}

.PHONY: all clean distclean


SRCS := $(wildcard ${CORE_DIR}/src/*.cc) \
		$(wildcard ${CORE_DIR}/src/Matriplex/*.cc) \
		$(wildcard ${SADIR}/*.cc)

ifdef WITH_ROOT
SRCS += $(wildcard ${SADIR}/RntDumper/*.cc)
endif
ifdef MKFIT_TRACE
SRCS += $(wildcard ${SADIR}/RdfTrace/*.cc)
endif

SRCB := $(notdir ${SRCS})
DEPS := $(SRCB:.cc=.d)
OBJS := $(SRCB:.cc=.o)

vpath %.cc ${CORE_DIR}/src ${CORE_DIR}/src/Matriplex ${SADIR}

AUTO_TGTS :=

# Begin ROOT dumpers
ifdef WITH_ROOT
LIB_ROOTDF_NAME := libMkFitRootDataFormats
LIB_ROOTDF := ../${LIB_ROOTDF_NAME}.so

ROOTDF_LINKDEF_H := ${SADIR}/DataFormats/RootDFs_LinkDef.h
ROOTDF_DICT := RootDFs_Dict
ROOTDF_DICT_H := ${SADIR}/DataFormats/${ROOTDF_DICT}.h
ROOTDF_DICT_CC := ${ROOTDF_DICT}.cc
ROOTDF_DICT_O := ${ROOTDF_DICT}.o
ROOTDF_DICT_PCM := ../${ROOTDF_DICT}_rdict.pcm
DEPS += ${ROOTDF_DICT_CC}.d
TGTS += ${LIB_ROOTDF}

ROOTDF_SRCS := $(wildcard ${SADIR}/DataFormats/*.cc)
ROOTDF_SRCB := $(notdir ${ROOTDF_SRCS})
ROOTDF_OBJS := $(ROOTDF_SRCB:.cc=.o) ${ROOTDF_DICT_O}
DEPS += $(ROOTDF_SRCB:.cc=.d)

LIB_CORE_DEPS := ${LIB_ROOTDF}
LIB_CORE_LD_EXTRA := -L.. -lMkFitRootDataFormats

vpath %.cc ${SADIR}/DataFormats ${SADIR}/RntDumper ${SADIR}/RdfTrace
endif
# End ROOT dumpers

all: ${TGTS}

# Begin Matriplex

auto-matriplex:
### Matriplex auto-generated fragments are included in source.
###	${MAKE} -f ${CORE_DIR}/src/Matriplex auto && touch $@
	touch $@

AUTO_TGTS += auto-matriplex

${DEPS}: auto-matriplex

# End Matriplex

ifeq ($(filter clean-local clean distclean, ${MAKECMDGOALS}),)
include ${DEPS}
endif

clean-local:
	-rm -f ${TGTS} *.d *.o *.om *.so
	-rm -f ${ROOTDF_DICT_CC} auto-matriplex main.dSYM
	-rm -f plotting/*.so plotting/*.d plotting/*.pcm

clean: clean-local

distclean: clean-local
	-rm -f ${AUTO_TGTS}
	-rm -f *.optrpt
	-rm -f ${LIB_CORE}
###	cd Matriplex && ${MAKE} distclean

${LIB_CORE}: ${OBJS} ${LIB_CORE_DEPS}
	@mkdir -p $(@D)
	${CXX} ${CXXFLAGS} ${VEC_HOST} ${OBJS} -shared -o $@ ${LDFLAGS_HOST} ${LDFLAGS} ${LIB_CORE_LD_EXTRA}

ifdef WITH_ROOT
ROOTDF_DICT_IN := ${ROOTDF_DICT_H} ${ROOTDF_LINKDEF_H}
${ROOTDF_DICT_CC} ${ROOTDF_DICT_PCM} &: ${ROOTDF_DICT_IN}
	rootcling --umbrellaHeader ${CPPFLAGS} -f ${ROOTDF_DICT_CC} ${ROOTDF_DICT_IN}
	mv ${ROOTDF_DICT}_rdict.pcm ${ROOTDF_DICT_PCM}

${LIB_ROOTDF}: ${ROOTDF_OBJS}
	${CXX} ${CXXFLAGS} ${VEC_HOST} ${ROOTDF_OBJS} -shared -o $@ ${LDFLAGS_HOST} ${LDFLAGS} -lROOTEve
endif

${OBJS} ${ROOTDF_OBJS}: %.o: %.cc %.d
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -c -o $@ $<

%.s: %.cc
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -S -o $@ $<

%.d: %.cc
	${MAKEDEPEND} -o $@ $<

ifdef WITH_ROOT
${ROOTDF_DICT_CC}.d: ${ROOTDF_DICT_H}
	${MAKEDEPEND_BASE} -MQ ${ROOTDF_DICT_CC} -o $@ $<
endif

echo:
	@echo "CXX=${CXX}"
	@echo SRCS = ${SRCS}
	@echo DEPS = ${DEPS}
	@echo OBJS = ${OBJS}
	@echo ROOT_OBJS = ${ROOT_OBJS}
	@echo CPPFLAGS = ${CPPFLAGS}
	@echo LDFLAGS = ${LDFLAGS}
	@echo LDFLAGS_HOST = ${LDFLAGS_HOST}
	@echo LIBCORE_LD_EXTRA = ${LIB_CORE_LD_EXTRA}

echo_cc_defs:
	${CXX} -dM -E -mavx2 - < /dev/null

echo-srcs:
	@echo ${SRCS}

echo-flags:
	@echo "CPPFLAGS=${CPPFLAGS}"

echo-tbb:
	@echo "TBB_GCC=${TBB_GCC}, TBB_PREFIX=${TBB_PREFIX}, TBB_ROOT=${TBB_ROOT}"

# Flags for compiling a study OUTSIDE this build with ACLiC -- see
# mkFit-external/mkfit-standalone-seedgeom/seedgeom-load.C, which reads them
# from here rather than hardcoding a copy.
#
# ACLIC_CPP is ABI-critical, not a convenience: standalone/Event.h gains members
# under MKFIT_TRACE, so an ACLiC build compiled without the identical -D set does
# not fail to link -- it silently reads every member at the wrong offset.
# ACLIC_OPT keeps the numerics identical (-Ofast implies fast-math, so borderline
# float comparisons differ without it) and makes timings comparable.
# NOTE the split: this emits only what the mkFit build owns. ROOT's own -I and
# -std are deliberately NOT included -- ACLiC supplies those from the RUNNING
# ROOT, and root-config here may resolve to a different one (without ROOTSYS it
# falls back to system ROOT and injects -std=c++17).
echo-aclic:
	@echo "ACLIC_DEFS=$(filter -D%,${CPPFLAGS})"
	@echo "ACLIC_OPT=${OPT} ${VEC_HOST} $(filter -f% -m%,${CXXFLAGS})"
