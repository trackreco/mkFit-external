include ${SADIR}/Makefile.config

CPPFLAGS := -I. -I${SRCDIR} -I${BLDDIR}/mkFit-external ${CPPFLAGS}

CMS_DIR := ${SRCDIR}/RecoTracker/MkFitCMS

LIB_CMS := ../libMkFitCMS.so
MAIN    := ../mkFit
MAIN_DEPS := ${LIB_CMS}
MAIN_LIBS := -lMkFitCore -lMkFitCMS
WRMEMF  := ../writeMemoryFile
WMF_DICT_PCM := ../WriteMemFileDict_rdict.pcm
SHELL_DICT_PCM := ../Shell_Dict_rdict.pcm
ROOTOUT := WriteMemFileDict.cc Shell_Dict.cc

TGTS := ${LIB_CMS} ${MAIN}

ifdef WITH_ROOT
MAIN_DEPS += ../libMkFitRootDataFormats.so
MAIN_LIBS += -lMkFitRootDataFormats
TGTS += ${WRMEMF} ${WMF_DICT_PCM} ${SHELL_DICT_PCM}
endif

.PHONY: all clean distclean

all: ${TGTS}

SRCS := $(wildcard ${CMS_DIR}/src/*.cc) $(wildcard ${SACMS}/*.cc)

ifdef WITH_ROOT
SRCS += ${SACMS}/tkNtuple/WriteMemoryFile.cc
WriteMemFileDict.cc ${WMF_DICT_PCM} &: ${SACMS}/tkNtuple/DictsLinkDef.h
	rootcling ${CPPFLAGS} -f WriteMemFileDict.cc $<
	mv WriteMemFileDict_rdict.pcm ${WMF_DICT_PCM}

SRCS += Shell_Dict.cc
SHELL_DICT_IN := ${SACMS}/Shell_Dict.h ${SACMS}/Shell_LinkDef.h
Shell_Dict.cc ${SHELL_DICT_PCM} &: ${SHELL_DICT_IN}
	rootcling --umbrellaHeader ${CPPFLAGS} -f Shell_Dict.cc ${SHELL_DICT_IN}
	mv Shell_Dict_rdict.pcm ${SHELL_DICT_PCM}
else
SRCS := $(filter-out ${SACMS}/Shell.cc ${SACMS}/Shell-LST.cc, ${SRCS})
endif

SRCB := $(notdir ${SRCS})
DEPS := $(SRCB:.cc=.d)
OBJS := $(SRCB:.cc=.o)

CMS_OBJS := $(filter-out mkFit.o WriteMemoryFile.o, ${OBJS})

vpath %.cc ${CMS_DIR}/src ${SACMS} ${SACMS}/tkNtuple

ifdef WITH_ROOT
DEPS += Shell_Dict.cc.d
endif

ifeq ($(filter clean-local clean distclean, ${MAKECMDGOALS}),)
include ${DEPS}
endif

clean-local:
	-rm -f ${TGTS} *.d *.o *.om *.so *.pcm ${ROOTOUT}
	-rm -rf main.dSYM
	-rm -rf plotting/*.so plotting/*.d plotting/*.pcm

clean: clean-local

distclean: clean-local
	-rm -f *.optrpt
	-rm -f ${TGTS}

${LIB_CMS}: ${CMS_OBJS}
	@mkdir -p $(@D)
	${CXX} ${CXXFLAGS} ${VEC_HOST} ${CMS_OBJS} -shared -o $@ ${LDFLAGS_HOST} ${LDFLAGS}

${MAIN}: ${MAIN_DEPS} mkFit.o
	${CXX} ${CXXFLAGS} ${VEC_HOST} ${LDFLAGS} mkFit.o -o $@ ${LDFLAGS_HOST} -ltbb -L.. ${MAIN_LIBS} -Wl,-rpath=.

${WRMEMF}: ${MAIN_DEPS} WriteMemoryFile.o WriteMemFileDict.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} WriteMemoryFile.o WriteMemFileDict.o -o $@ ${LDFLAGS_HOST} -ltbb -L.. ${MAIN_LIBS} -Wl,-rpath=.

${OBJS}: %.o: %.cc %.d
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${VEC_HOST} -c -o $@ $<

%.d: %.cc
	${MAKEDEPEND} -o $@ $<

ifdef WITH_ROOT
Shell_Dict.cc.d: ${SACMS}/Shell_Dict.h
	${MAKEDEPEND_BASE} -MQ Shell_Dict.cc -o $@ $<
endif

echo:
	@echo "CXX=${CXX}"
	@echo SRCS = ${SRCS}
	@echo DEPS = ${DEPS}
	@echo OBJS = ${OBJS}
	@echo MAIN_LIBS = ${MAIN_LIBS}

echo_cc_defs:
	${CXX} -dM -E -mavx2 - < /dev/null

echo-srcs:
	@echo ${SRCS}

echo-flags:
	@echo "CPPFLAGS=${CPPFLAGS}"

echo-tbb:
	@echo "TBB_GCC=${TBB_GCC}, TBB_PREFIX=${TBB_PREFIX}, TBB_ROOT=${TBB_ROOT}"
