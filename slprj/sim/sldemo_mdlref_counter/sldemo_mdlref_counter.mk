# Copyright 1994-2016 The MathWorks, Inc.
#
# File    : ert_unix.tmf   
#
# Abstract:
#       Template makefile for building a UNIX-based stand-alone embedded 
#       real-time version of Simulink model using generated C code.
#
#       This makefile attempts to conform to the guidelines specified in the
#       IEEE Std 1003.2-1992 (POSIX) standard. It is designed to be used
#       with GNU Make which is located in matlabroot/rtw/bin.
#
#       Note that this template is automatically customized by the build 
#       procedure to create "<model>.mk"
#
#       The following defines can be used to modify the behavior of the
#       build:
#         OPT_OPTS       - Optimization options. Default is -O.
#         CPP_OPTS       - C++ compiler options.	
#         OPTS           - User specific compile options.
#         USER_SRCS      - Additional user sources, such as files needed by
#                          S-functions.
#         USER_INCLUDES  - Additional include paths
#                          (i.e. USER_INCLUDES="-Iwhere-ever -Iwhere-ever2")
#
#       To enable debugging:
#         set DEBUG_BUILD = 1 below, which will trigger OPTS=-g and
#          LDFLAGS += -g (may vary with compiler version, see compiler doc) 
#
#       This template makefile is designed to be used with a system target
#       file that contains 'rtwgensettings.BuildDirSuffix' see ert.tlc


#------------------------ Macros read by make_rtw ------------------------------
#
# The following macros are read by the build procedure:
#
#  MAKECMD         - This is the command used to invoke the make utility
#  HOST            - What platform this template makefile is targeted for
#                    (i.e. PC or UNIX)
#  BUILD           - Invoke make from the build procedure (yes/no)?
#  SYS_TARGET_FILE - Name of system target file.

MAKECMD         = /Applications/MATLAB_R2017b.app/bin/maci64/gmake
HOST            = UNIX
BUILD           = yes
SYS_TARGET_FILE = any
COMPILER_TOOL_CHAIN = unix

#---------------------- Tokens expanded by make_rtw ----------------------------
#
# The following tokens, when wrapped with "|>" and "<|" are expanded by the
# build procedure.
#
#  MODEL_NAME          - Name of the Simulink block diagram
#  MODEL_MODULES       - Any additional generated source modules
#  MAKEFILE_NAME       - Name of makefile created from template makefile <model>.mk
#  MATLAB_ROOT         - Path to where MATLAB is installed.
#  S_FUNCTIONS         - List of additional S-function modules.
#  S_FUNCTIONS_LIB     - List of S-functions libraries to link. 
#  NUMST               - Number of sample times
#  NCSTATES            - Number of continuous states
#  COMPUTER            - Computer type. See the MATLAB computer command.
#  BUILDARGS           - Options passed in at the command line.
#  MULTITASKING        - yes (1) or no (0): Is solver mode multitasking
#  INTEGER_CODE        - yes (1) or no (0): Is generated code purely integer
#  MAT_FILE            - yes (1) or no (0): Should mat file logging be done,
#                        if 0, the generated code runs indefinitely
#  EXT_MODE            - yes (1) or no (0): Build for external mode
#  TMW_EXTMODE_TESTING - yes (1) or no (0): Build ext_test.c for external mode
#                        testing.
#  EXTMODE_TRANSPORT   - Index of transport mechanism (e.g. tcpip, serial) for extmode
#  MULTI_INSTANCE_CODE - Is the generated code multi instantiable (1/0)?
#  GEN_MDLREF_SFCN     - (1/0): are we generating model reference wrapper s-function
#  MODELREFS           - List of referenced models
#  EXTMODE_STATIC      - yes (1) or no (0): Use static instead of dynamic mem alloc.
#  EXTMODE_STATIC_SIZE - Size of static memory allocation buffer.
#  TGT_FCN_LIB         - Target Funtion library to use
#  PORTABLE_WORDSIZES  - Is this build intented for SIL simulation with portable word sizes (1/0)?
#  SHRLIBTARGET        - Is this build intended for generation of a shared library instead 
#                        of executable (1/0)?
#  MAKEFILEBUILDER_TGT - Is this build performed by the MakefileBuilder class
#                        e.g. to create a PIL executable?
#  STANDALONE_SUPPRESS_EXE - Build the standalone target but only create object code modules 
#                            and do not build an executable

MODEL                   = sldemo_mdlref_counter
MODULES                 = sldemo_mdlref_counter.c sldemo_mdlref_counter_capi.c
MAKEFILE                = sldemo_mdlref_counter.mk
MATLAB_ROOT             = /Applications/MATLAB_R2017b.app
ALT_MATLAB_ROOT         = /Applications/MATLAB_R2017b.app
MASTER_ANCHOR_DIR       = 
START_DIR               = /Users/ghleila/Downloads/WaveEnergyConversion/wec2/wec
S_FUNCTIONS             = 
S_FUNCTIONS_LIB         = 
NUMST                   = 1
NCSTATES                = 0
COMPUTER                = MACI64
BUILDARGS               =  GENERATE_ASAP2=0 OPTS="-DSLMSG_ALLOW_SYSTEM_ALLOC -DMATLAB_MEX_FILE"
MULTITASKING            = 0
INTEGER_CODE            = 0
MAT_FILE                = 0
ALLOCATIONFCN           = 0
ONESTEPFCN              = 0
TERMFCN                 = 1
B_ERTSFCN               = 0
MEXEXT                  = mexmaci64
ENABLE_SLEXEC_SSBRIDGE  = 11
EXT_MODE                = 0
TMW_EXTMODE_TESTING     = 0
EXTMODE_TRANSPORT       = 0
EXTMODE_STATIC          = 0
EXTMODE_STATIC_SIZE     = 1000000
MULTI_INSTANCE_CODE     = 0
CODE_INTERFACE_PACKAGING = Nonreusable function
CLASSIC_INTERFACE       = 0
TGT_FCN_LIB             = None
MODELREFS               = 
SHARED_SRC              = ../../../slprj/sim/_sharedutils/*.c*
SHARED_SRC_DIR          = ../../../slprj/sim/_sharedutils
SHARED_BIN_DIR          = ../../../slprj/sim/_sharedutils
SHARED_LIB              = ../../../slprj/sim/_sharedutils/rtwshared.a
GEN_SAMPLE_MAIN         = 0
TARGET_LANG_EXT         = c
PORTABLE_WORDSIZES      = 0
SHRLIBTARGET            = 0
MAKEFILEBUILDER_TGT     = 0
STANDALONE_SUPPRESS_EXE = 0
OPTIMIZATION_FLAGS      = -O0 -DNDEBUG
ADDITIONAL_LDFLAGS      = 
DEFINES_CUSTOM          = 
SYSTEM_LIBS             =  -L$(MATLAB_ROOT)/bin/maci64 -lm -lfixedpoint -lmwmathutil -lmwipp -lmwsl_AsyncioQueue -lmwsf_runtime -lmwslexec_simlog -lmwcoder_ToAsyncQueueTgtAppSvc -lmwsl_simtarget_instrumentation -lmwsl_simtarget_core -lmwsimulink

# To enable debugging:
# set DEBUG_BUILD = 1
DEBUG_BUILD             = 0

#--------------------------- Model and reference models -----------------------
MODELLIB                  = sldemo_mdlref_counterlib.a
MODELREF_LINK_LIBS        = 
MODELREF_INC_PATH         = 
RELATIVE_PATH_TO_ANCHOR   = ../../..
# NONE: standalone, SIM: modelref sim, RTW: modelref coder target
MODELREF_TARGET_TYPE       = SIM
MODELREF_SFCN_SUFFIX       = _msf
ISPROTECTINGMODEL          = NOTPROTECTING
PROT_CAPIC_SUFFIX          = _capi.c
PROT_CAPIO_SUFFIX          = _capi_host.o


#-- For gcc, link model reference libraries as a group to resolve circular references ---
#   (NOTE: Clang does not recognize these switches but it already resolves circular references)
MODELREF_LINK_LIBS_GROUP = 
ifeq ($(COMPUTER),GLNXA64)
  MODELREF_LINK_LIBS_GROUP = -Wl,--start-group  -Wl,--end-group
endif

#-- In the case when directory name contains space ---
ifneq ($(MATLAB_ROOT),$(ALT_MATLAB_ROOT))
MATLAB_ROOT := $(ALT_MATLAB_ROOT)
endif

#----------------------------- External mode -----------------------------------
# Uncomment -DVERBOSE to have information printed to stdout
# To add a new transport layer, see the comments in
#   <matlabroot>/toolbox/simulink/simulink/extmode_transports.m
ifeq ($(EXT_MODE),1)
  EXT_CC_OPTS = -DEXT_MODE -D$(COMPUTER) #-DVERBOSE
  EXT_LIB     =
  EXT_SRC     = ext_svr.c updown.c ext_work.c
  LINT_EXT_COMMON_SRC =
  LINT_EXT_TCPIP_SRC  =
  ifeq ($(EXTMODE_TRANSPORT),0) #tcpip
    EXT_SRC +=  rtiostream_interface.c rtiostream_tcpip.c
    LINT_EXT_COMMON_SRC = ext_svr.c updown.c ext_work.c
    LINT_EXT_TCPIP_SRC  = rtiostream_interface.c rtiostream_tcpip.c
  endif
  ifeq ($(EXTMODE_TRANSPORT),1) #serial
    EXT_SRC += ext_svr_serial_transport.c ext_serial_pkt.c
    EXT_SRC += rtiostream_serial_interface.c rtiostream_serial.c
  endif
  ifeq ($(TMW_EXTMODE_TESTING),1)
    EXT_SRC     += ext_test.c
    EXT_CC_OPTS += -DTMW_EXTMODE_TESTING
  endif
  ifeq ($(EXTMODE_STATIC),1)
    EXT_SRC             += mem_mgr.c
    LINT_EXT_COMMON_SRC += mem_mgr.c
    EXT_CC_OPTS         += -DEXTMODE_STATIC -DEXTMODE_STATIC_SIZE=$(EXTMODE_STATIC_SIZE)
  endif
  ifeq ($(COMPUTER),SOL64)
    EXT_LIB = -lsocket -lnsl
  endif
endif

#--------------------------- Tool Specifications -------------------------------

include $(MATLAB_ROOT)/rtw/c/tools/unixtools.mk

# Determine if we are generating an s-function
SFCN = 0
ifeq ($(MODELREF_TARGET_TYPE),SIM)
  SFCN = 1
endif
ifeq ($(B_ERTSFCN),1)
  SFCN = 1
endif

# Use GCC_TEST to do a test compile of the local source (add DO_GCC_TEST=1)
ifeq ($(DO_GCC_TEST), 1)
  GCC_TEST     = echo "Doing gcc test compile"; gcc -c -o /dev/null -Wall 
  GCC_TEST_OUT = 
else
  GCC_TEST     = echo
  GCC_TEST_OUT = > /dev/null
endif

#------------------------------ Include Path -----------------------------------

MATLAB_INCLUDES = \
	-I$(MATLAB_ROOT)/extern/include \
	-I$(MATLAB_ROOT)/simulink/include \
	-I$(MATLAB_ROOT)/rtw/c/src \
	-I$(MATLAB_ROOT)/rtw/c/src/ext_mode/common

# Additional includes 
ADD_INCLUDES = \
	-I$(START_DIR) \
	-I$(START_DIR)/slprj/sim/sldemo_mdlref_counter \
	-I$(MATLAB_ROOT)/toolbox/simulink/simdemos/simfeatures/modelreference \


SHARED_INCLUDES =
ifneq ($(SHARED_SRC_DIR),)
SHARED_INCLUDES = -I$(SHARED_SRC_DIR) 
endif

INCLUDES = -I. -I$(RELATIVE_PATH_TO_ANCHOR) $(MATLAB_INCLUDES) $(ADD_INCLUDES) $(USER_INCLUDES) \
	$(INSTRUMENT_INCLUDES) $(MODELREF_INC_PATH) $(SHARED_INCLUDES)


#-------------------------------- C Flags --------------------------------------

# Optimization Options
ifndef OPT_OPTS
OPT_OPTS = $(DEFAULT_OPT_OPTS)
endif

# General User Options
ifeq ($(DEBUG_BUILD),0)
DBG_FLAG =
else
#   Set OPTS=-g and any additional flags for debugging
DBG_FLAG = -g
LDFLAGS += -g
endif

LDFLAGS += -Wl,-rpath,"$(MATLAB_ROOT)/bin/$(ARCH)" -L"$(MATLAB_ROOT)/bin/$(ARCH)"

# Compiler options, etc: 
CPP_REQ_DEFINES1 = -DMODEL=$(MODEL) -DNUMST=$(NUMST) -DNCSTATES=$(NCSTATES) \
		  -DUNIX -DMAT_FILE=$(MAT_FILE) -DINTEGER_CODE=$(INTEGER_CODE) \
		  -DONESTEPFCN=$(ONESTEPFCN) -DTERMFCN=$(TERMFCN) \
		  -DHAVESTDIO -DMULTI_INSTANCE_CODE=$(MULTI_INSTANCE_CODE) \
		  -DCLASSIC_INTERFACE=$(CLASSIC_INTERFACE) \
		  -DALLOCATIONFCN=$(ALLOCATIONFCN)

ifneq ($(ENABLE_SLEXEC_SSBRIDGE), 0)
    ifeq ($(MODELREF_TARGET_TYPE),SIM)
        CPP_REQ_DEFINES1 += -DENABLE_SLEXEC_SSBRIDGE=$(ENABLE_SLEXEC_SSBRIDGE)
    endif
endif

CPP_REQ_DEFINES = $(CPP_REQ_DEFINES1)
ifeq ($(MODELREF_TARGET_TYPE),SIM)
CPP_REQ_DEFINES += -DMDL_REF_SIM_TGT=1 
else
CPP_REQ_DEFINES += -DMT=$(MULTITASKING)
endif

ifeq ($(PORTABLE_WORDSIZES),1)
CPP_REQ_DEFINES += -DPORTABLE_WORDSIZES
endif

ifneq ($(OPTIMIZATION_FLAGS),)
CC_OPTS = $(OPTS) $(EXT_CC_OPTS) $(OPTIMIZATION_FLAGS)
else
CC_OPTS = $(OPT_OPTS) $(OPTS) $(EXT_CC_OPTS)
endif

CFLAGS = $(ANSI_OPTS) $(DBG_FLAG) $(CC_OPTS) $(DEFINES_CUSTOM) $(CPP_REQ_DEFINES) $(INCLUDES)
CPPFLAGS = $(CPP_ANSI_OPTS) $(DBG_FLAG) $(CPP_OPTS) $(CC_OPTS) $(CPP_REQ_DEFINES) $(INCLUDES)
ifeq ($(SFCN),1)
ifneq ($(OPTIMIZATION_FLAGS),)
CC_OPTS_SFCN = COPTIMFLAGS="$(ANSI_OPTS) $(OPTIMIZATION_FLAGS)"
else
CC_OPTS_SFCN = COPTIMFLAGS="$(OPT_OPTS) $(ANSI_OPTS)"
endif
CFLAGS_SFCN = $(CC_OPTS_SFCN) $(CPP_REQ_DEFINES1) $(INCLUDES)
endif

#-------------------------- Additional Libraries ------------------------------

SYSTEM_LIBS += $(EXT_LIB)
ifeq ($(SFCN),0)
SYSTEM_LIBS += -lm
endif

LIBS =
 
LIBS += $(S_FUNCTIONS_LIB) $(INSTRUMENT_LIBS)

ifeq ($(SFCN),1)	
LIBFIXPT = -L$(MATLAB_ROOT)/bin/$(ARCH) -lfixedpoint
else
LIBFIXPT = 
endif

MODELREF_LIB_DEF = -L$(MATLAB_ROOT)/bin/$(ARCH)

ifneq ($(ENABLE_SLEXEC_SSBRIDGE), 0)
  MODELREF_LIB_DEF += -lmwslexec_simbridge
endif

MODELREF_LIB_DEF += -lmwmathutil
MODELREF_LIB_DEF += -lmwsl_fileio
MODELREF_LIB_DEF += -lmwsl_simtarget_instrumentation
MODELREF_LIB_DEF += -lmwsl_simtarget_core
MODELREF_LIB_DEF += -lmwcoder_ToAsyncQueueTgtAppSvc
MODELREF_LIB_DEF += -lmwsl_AsyncioQueue
MODELREF_LIB_DEF += -lmwsf_runtime
MODELREF_LIB_DEF += -ldl
MODELREF_LIB_DEF += -lmwslio_core
MODELREF_LIB_DEF += -lmwslio_clients
MODELREF_LIB_DEF += -lmwsimulink


ifeq ($(MODELREF_TARGET_TYPE),SIM)
LIBMWIPP = -L$(MATLAB_ROOT)/bin/$(ARCH) -lmwipp
else
LIBMWIPP = 
endif

#----------------------------- Source Files ------------------------------------
ADD_SRCS =

ifneq ($(ADDITIONAL_LDFLAGS),)
MEX_LDFLAGS = CLIBS='$$CLIBS $(ADDITIONAL_LDFLAGS)' LINKLIBS='$$LINKLIBS $(MODELLIB) $(MODELREF_LINK_LIBS_GROUP) $(SHARED_LIB)'
else
MEX_LDFLAGS = LINKLIBS='$(MODELLIB) $(MODELREF_LINK_LIBS_GROUP) $(SHARED_LIB) $$LINKLIBS'
endif

ifeq ($(SFCN),0)
  SRCS  = $(ADD_SRCS) $(S_FUNCTIONS)
  SRC_DEP =
  ifeq ($(MODELREF_TARGET_TYPE), NONE)
    ifeq ($(SHRLIBTARGET), 1)
      # Shared object/dynamic library
      ifeq ($(COMPUTER),MACI64)
        PRODUCT            = $(MODEL).dylib
      else
        PRODUCT            = $(MODEL).so
      endif
      BIN_SETTING        = $(LD) $(SHRLIBLDFLAGS)$(MODEL).def -o $(PRODUCT)
      BUILD_PRODUCT_TYPE = "shared object"
      SRCS               += $(MODULES) $(MODEL).$(TARGET_LANG_EXT) $(EXT_SRC)
    else
      ifeq ($(MAKEFILEBUILDER_TGT), 1)
        # Standalone executable (e.g. for PIL)
        PRODUCT            = $(MODEL)
        BIN_SETTING        = $(LD) $(LDFLAGS) $(ADDITIONAL_LDFLAGS) -o $(PRODUCT)
        BUILD_PRODUCT_TYPE = "executable"
      else
        SRCS               += $(MODULES) $(MODEL).$(TARGET_LANG_EXT) $(EXT_SRC)
        ifeq ($(STANDALONE_SUPPRESS_EXE), 1)
          # Build object code only for top level model
          PRODUCT            = "ObjectModulesOnly"
          BUILD_PRODUCT_TYPE = "object modules"
        else
          # ERT standalone
          PRODUCT            = $(RELATIVE_PATH_TO_ANCHOR)/$(MODEL)
          BIN_SETTING        = $(LD) $(LDFLAGS) $(ADDITIONAL_LDFLAGS) -o $(PRODUCT)
          BUILD_PRODUCT_TYPE = "executable"
          ifeq ($(GEN_SAMPLE_MAIN),0)
             ifeq ($(MULTI_INSTANCE_CODE), 1)
                 ifeq "$(CODE_INTERFACE_PACKAGING)" "C++ class"
                     SRCS += rt_cppclass_main.cpp
                 else
                     SRCS += rt_malloc_main.c
                 endif
             else
                SRCS += rt_main.c
             endif
          else
            SRCS += ert_main.$(TARGET_LANG_EXT)
          endif
        endif
      endif
    endif
  else
    # Model reference coder target
    SRCS               += $(MODULES)
    PRODUCT            = $(MODELLIB)
    BUILD_PRODUCT_TYPE = "library"
  endif
else
  # Model Reference Simulation target, ERT S-function target
  MEX                 = $(MATLAB_ARCH_BIN)/mex
  ifeq ($(MODELREF_TARGET_TYPE), SIM)
  PRODUCT            = $(RELATIVE_PATH_TO_ANCHOR)/$(MODEL)$(MODELREF_SFCN_SUFFIX).$(MEXEXT)
  RTW_SFUN_SRC       = $(MODEL)$(MODELREF_SFCN_SUFFIX).$(TARGET_LANG_EXT)
  SRCS               = $(MODULES) $(ADD_SRCS)
  CFLAGS_SFCN       += -DMDL_REF_SIM_TGT=1
  else
  PRODUCT            = $(RELATIVE_PATH_TO_ANCHOR)/$(MODEL)_sf.$(MEXEXT)
  RTW_SFUN_SRC       = $(MODEL)_sf.$(TARGET_LANG_EXT)
  SRCS               = $(MODULES) $(ADD_SRCS) $(S_FUNCTIONS)

  ifeq ($(MODELREF_TARGET_TYPE),NONE)
  MODELLIB           = 
  endif

  endif
  BIN_SETTING        = $(MEX) -MATLAB_ARCH=$(ARCH) $(CFLAGS_SFCN) $(RTW_SFUN_SRC) $(MEX_LDFLAGS) -outdir $(RELATIVE_PATH_TO_ANCHOR) -silent
  BUILD_PRODUCT_TYPE = "mex file"
  ifeq ($(B_ERTSFCN),1)
    SRCS              += $(MODEL).$(TARGET_LANG_EXT)
  endif
  SRC_DEP            = $(RTW_SFUN_SRC)
endif

USER_SRCS =

USER_OBJS       = $(addsuffix .o, $(basename $(USER_SRCS)))
LOCAL_USER_OBJS = $(notdir $(USER_OBJS))

OBJS      = $(addsuffix .o, $(basename $(SRCS))) $(USER_OBJS)
PROT_CAPIC  = $(addsuffix $(PROT_CAPIC_SUFFIX), $(MODEL))
PROT_CAPIO  = $(addsuffix $(PROT_CAPIO_SUFFIX), $(MODEL))
LINK_OBJS = $(addsuffix .o, $(basename $(SRCS))) $(LOCAL_USER_OBJS)

SHARED_SRC := $(wildcard $(SHARED_SRC))
SHARED_OBJS_NO_PATH = $(addsuffix .o, $(basename $(notdir $(SHARED_SRC))))
SHARED_OBJS = $(addprefix $(join $(SHARED_BIN_DIR), /), $(SHARED_OBJS_NO_PATH))

ADDITIONAL_LDFLAGS += $(ARCH_SPECIFIC_LDFLAGS)

#----------------------------- Lint (sol2 only) --------------------------------

ifeq "$(CODE_INTERFACE_PACKAGING)" "C++ class"
  LINT_SRCS = $(MATLAB_ROOT)/rtw/c/src/common/rt_cppclass_main.cpp
else
  LINT_SRCS = $(MATLAB_ROOT)/rtw/c/src/common/rt_main.c
endif
LINT_SRCS += $(MODEL).$(TARGET_LANG_EXT) $(MODULES) $(USER_SRCS) $(S_FUNCTIONS)
LINTOPTSFILE = $(MODEL).lintopts

ifneq ($(LINT_EXT_COMMON_SRC), )
  LINT_SRCS += $(MATLAB_ROOT)/rtw/c/src/ext_mode/common/$(LINT_EXT_COMMON_SRC)
endif
ifneq ($(LINT_EXT_TCPIP_SRC), )
  LINT_SRCS += $(MATLAB_ROOT)/rtw/c/src/rtiostream/rtiostreamtcpip/$(LINT_EXT_TCPIP_SRC)
endif

LINT_ERROFF1 = E_NAME_DEF_NOT_USED2,E_NAME_DECL_NOT_USED_DEF2
LINT_ERROFF2 = $(LINT_ERROFF1),E_FUNC_ARG_UNUSED
LINT_ERROFF  = $(LINT_ERROFF2),E_INDISTING_FROM_TRUNC2,E_NAME_USED_NOT_DEF2

#--------------------------------- Rules ---------------------------------------
ifeq ($(MODELREF_TARGET_TYPE),NONE)
  ifeq ($(SHRLIBTARGET), 1)
    ifeq ($(GEN_SAMPLE_MAIN), 1)
$(MODEL) : $(PRODUCT) ert_main.o
	$(LD) $(LDFLAGS) $(ADDITIONAL_LDFLAGS) -o $(MODEL) ert_main.o $(PRODUCT) $(SYSTEM_LIBS)
	@mv $(PRODUCT) $(RELATIVE_PATH_TO_ANCHOR)/
	@echo "### Created $(BUILD_PRODUCT_TYPE): $(PRODUCT)"
	@mv $(MODEL) $(RELATIVE_PATH_TO_ANCHOR)/
	@echo "### Created executable: $@"
    else
mvshrlib: $(PRODUCT)
	@mv $(PRODUCT) $(RELATIVE_PATH_TO_ANCHOR)/
    endif
  endif
  ifeq ($(MAKEFILEBUILDER_TGT), 1)
$(PRODUCT) : $(OBJS) $(SHARED_LIB) $(LIBS) $(SRC_DEP) $(MODELLIB) $(MODELREF_LINK_LIBS)
	$(BIN_SETTING) $(LINK_OBJS) $(MODELLIB) $(MODELREF_LINK_LIBS_GROUP) $(SHARED_LIB) $(LIBFIXPT) $(LIBS)  $(SYSTEM_LIBS)
	@echo "### Created $(BUILD_PRODUCT_TYPE): $@"
  else
    ifeq ($(STANDALONE_SUPPRESS_EXE), 1)
.PHONY: $(PRODUCT)
$(PRODUCT) : $(OBJS) $(SHARED_LIB) $(LIBS) $(SRC_DEP) $(MODELREF_LINK_LIBS)
	@echo "### Created $(BUILD_PRODUCT_TYPE): $@"
    else
      ifeq ($(SFCN),0)
$(PRODUCT) : $(OBJS) $(SHARED_LIB) $(LIBS) $(SRC_DEP) $(MODELREF_LINK_LIBS)
	$(BIN_SETTING) $(LINK_OBJS) $(MODELREF_LINK_LIBS_GROUP) $(SHARED_LIB) $(LIBFIXPT) $(LIBS)  $(SYSTEM_LIBS)
	@echo "### Created $(BUILD_PRODUCT_TYPE): $@"
      else
$(PRODUCT) : $(OBJS) $(SHARED_LIB) $(LIBS) $(SRC_DEP) $(MODELREF_LINK_LIBS)
	$(BIN_SETTING) $(LINK_OBJS) $(SHARED_LIB) $(LIBFIXPT) $(LIBS)  $(SYSTEM_LIBS)
	@echo "### Created $(BUILD_PRODUCT_TYPE): $@"
      endif
    endif
  endif
else
 ifeq ($(MODELREF_TARGET_TYPE),SIM)
  ifeq ($(ISPROTECTINGMODEL),PROTECTING)
  all : $(PRODUCT) $(PROT_CAPIO)
  endif
  $(PRODUCT) : $(OBJS) $(SHARED_LIB) $(LIBS) $(SRC_DEP)
	@rm -f $(MODELLIB)
	$(AR) ruvs $(MODELLIB) $(LINK_OBJS)
	@echo "### Created $(MODELLIB)"
	$(BIN_SETTING) $(LIBFIXPT) $(LIBS) $(LIBMWIPP) $(MODELREF_LIB_DEF) $(SYSTEM_LIBS)
	@echo "### Created $(BUILD_PRODUCT_TYPE): $@"
  ifeq ($(ISPROTECTINGMODEL),PROTECTING)
  $(PROT_CAPIO) : $(PROT_CAPIC)
	$(CC) -c -o $(PROT_CAPIO) $(CFLAGS) $(GCC_WALL_FLAG) -DHOST_CAPI_BUILD $(PROT_CAPIC)	
  endif
 else
  $(PRODUCT) : $(OBJS) $(SHARED_LIB) $(SRC_DEP)
	@rm -f $(MODELLIB)
	$(AR) ruvs $(MODELLIB) $(LINK_OBJS)
	@echo "### Created $(MODELLIB)"
	@echo "### Created $(BUILD_PRODUCT_TYPE): $@"
 endif
endif


#--------------- Support for building referenced models -----------------------


#-------------------------- Support for building modules ----------------------

ifneq ($(SHARED_SRC_DIR),)
$(SHARED_BIN_DIR)/%.o : $(SHARED_SRC_DIR)/%.c
	$(CC) -c -o $@ $(CFLAGS) $(GCC_WALL_FLAG_MAX) "$<"

$(SHARED_BIN_DIR)/%.o : $(SHARED_SRC_DIR)/%.cpp
	$(CPP) -c -o $@ $(CPPFLAGS) $(GCC_WALL_FLAG_MAX) "$<"
endif

%.o : %.c
	@$(GCC_TEST) $(CPP_REQ_DEFINES) $(INCLUDES) "$<" $(GCC_TEST_OUT)
	$(CC) -c $(CFLAGS) $(GCC_WALL_FLAG) "$<"

%.o : %.cpp
	@$(GCC_TEST) $(CPP_REQ_DEFINES) $(INCLUDES) "$<" $(GCC_TEST_OUT)
	$(CPP) -c $(CPPFLAGS) $(GCC_WALL_FLAG) "$<"

%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.c
	@$(GCC_TEST) $(CPP_REQ_DEFINES) $(INCLUDES) "$<" $(GCC_TEST_OUT)
	$(CC) -c $(CFLAGS) $(GCC_WALL_FLAG) "$<"

%.o : $(RELATIVE_PATH_TO_ANCHOR)/%.cpp
	@$(GCC_TEST) $(CPP_REQ_DEFINES) $(INCLUDES) "$<" $(GCC_TEST_OUT)
	$(CPP) -c $(CPPFLAGS) $(GCC_WALL_FLAG) "$<"

ifeq ($(GEN_SAMPLE_MAIN),0)
ifeq ($(TARGET_LANG_EXT),cpp)
%.o : $(MATLAB_ROOT)/rtw/c/src/common/%.cpp
	$(CPP) -c $(CPPFLAGS) $(GCC_WALL_FLAG_MAX) "$<"

%.o : $(MATLAB_ROOT)/rtw/c/src/common/%.c
	$(CPP) -c $(CPPFLAGS) $(GCC_WALL_FLAG_MAX) "$<"
else
%.o : $(MATLAB_ROOT)/rtw/c/src/common/%.c
	$(CC) -c $(CFLAGS) $(GCC_WALL_FLAG_MAX) "$<"
endif
endif

%.o : $(MATLAB_ROOT)/rtw/c/src/%.c
	$(CC) -c $(CFLAGS) $(GCC_WALL_FLAG_MAX) "$<"

%.o : $(MATLAB_ROOT)/rtw/c/src/ext_mode/common/%.c
	$(CC) -c $(CFLAGS) $(GCC_WALL_FLAG_MAX) "$<"

%.o : $(MATLAB_ROOT)/rtw/c/src/rtiostream/rtiostreamtcpip/%.c
	$(CC) -c $(CFLAGS) $(GCC_WALL_FLAG_MAX) "$<"

%.o : $(MATLAB_ROOT)/rtw/c/src/ext_mode/custom/%.c
	$(CC) -c $(CFLAGS) $(GCC_WALL_FLAG_MAX) "$<"

%.o : $(MATLAB_ROOT)/rtw/c/src/%.c
	$(CC) -c $(CFLAGS) $(GCC_WALL_FLAG_MAX) "$<"

%.o : $(MATLAB_ROOT)/simulink/src/%.c
	$(CC) -c $(CFLAGS) $(GCC_WALL_FLAG_MAX) "$<"



%.o : $(MATLAB_ROOT)/rtw/c/src/%.cpp
	$(CPP) -c $(CPPFLAGS) $(GCC_WALL_FLAG_MAX) "$<"

%.o : $(MATLAB_ROOT)/simulink/src/%.cpp
	$(CPP) -c $(CPPFLAGS) $(GCC_WALL_FLAG_MAX) "$<"



%.o : $(MATLAB_ROOT)/simulink/src/%.cpp
	$(CPP) -c $(CPPFLAGS) $(GCC_WALL_FLAG_MAX) "$<"

%.o : $(MATLAB_ROOT)/simulink/src/%.c
	$(CC) -c $(CFLAGS) $(GCC_WALL_FLAG_MAX) "$<"

#------------------------------- Libraries -------------------------------------





#----------------------------- Dependencies ------------------------------------

$(OBJS) : $(MAKEFILE) rtw_proj.tmw

$(SHARED_LIB) : $(SHARED_OBJS)
	@echo "### Creating $@ "
	$(AR) ruvs $@ $(SHARED_OBJS)
	@echo "### $@ Created "


#--------- Miscellaneous rules to purge, clean and lint (sol2 only) ------------

purge : clean
	@echo "### Deleting the generated source code for $(MODEL)"
	@\rm -f $(MODEL).c $(MODEL).h $(MODEL)_types.h $(MODEL)_data.c \
	        $(MODEL)_private.h $(MODEL).rtw $(MODULES) rtw_proj.tmw $(MAKEFILE)

clean :
	@echo "### Deleting the objects, $(PRODUCT)"
	@\rm -f $(LINK_OBJS) $(PRODUCT) 

lint  : rtwlib.ln
	@lint -errchk -errhdr=%user -errtags=yes -F -L. -lrtwlib -x -Xc \
	      -erroff=$(LINT_ERROFF) \
	      -D_POSIX_C_SOURCE $(CFLAGS) $(LINT_SRCS)
	@\rm -f $(LINTOPTSFILE)
	@echo
	@echo "### Created lint output only, no executable"
	@echo

rtwlib.ln : $(MAKEFILE) rtw_proj.tmw
	@echo
	@echo "### Linting ..."
	@echo
	@\rm -f llib-lrtwlib.ln $(LINTOPTSFILE)
	@echo "-dirout=. -errchk -errhdr=%user " >> $(LINTOPTSFILE)
	@echo "-errtags -F -ortwlib -x -Xc " >> $(LINTOPTSFILE)
	@echo "-erroff=$(LINT_ERROFF) " >> $(LINTOPTSFILE)
	@echo "-D_POSIX_C_SOURCE $(CFLAGS) " >> $(LINTOPTSFILE)
	@for file in $(MATLAB_ROOT)/rtw/c/libsrc/*.c; do \
	  echo "$$file " >> $(LINTOPTSFILE); \
	done
	lint -flagsrc=$(LINTOPTSFILE)


# EOF: ert_unix.tmf
