################################################################################
# File:     makefile
#
# Purpose:  Makefile for the b64 project, for the Arturius C/C++ Compiler Multiplexer
#           compiler, on Win32
#
# Created:  15th Feburary 2005
# Updated:  4th February 2012
#
################################################################################


################################################################################
# Symbols

# Identification

COMP_TAG            =   ar
PROJ_TAG            =   b64
PROJ_VER_MAJOR      =   1
PROJ_TAG_VER        =   $(PROJ_TAG).$(PROJ_VER_MAJOR)
COMP_ARCH_OS_TAG	=	$(COMP_TAG)$(OS_TAG)$(ARCH_TAG)

XCONTRACT_VER_MAJOR =   0
XCONTRACT_TAG_VER   =   xcontract.$(XCONTRACT_VER_MAJOR)

XTESTS_VER_MAJOR    =   0
XTESTS_TAG_VER      =   xtests.$(XTESTS_VER_MAJOR)

PROJ_MAKEFILE_NAME  =   makefile

########################################
# Testing verbosity

TEST_VERBOSITY		=	2

########################################
# Directories

PROJ_BASE_DIR       =   ../..
PROJ_BIN_DIR        =   $(PROJ_BASE_DIR)/bin
PROJ_INCLUDE_DIR    =   $(PROJ_BASE_DIR)/include
PROJ_LIB_DIR        =   $(PROJ_BASE_DIR)/lib
PROJ_OBJ_DIR        =   .
PROJ_SRC_DIR        =   $(PROJ_BASE_DIR)/src
PROJ_MAPPINGS_DIR   =   $(PROJ_BASE_DIR)/mappings
PROJ_TEST_UNIT_DIR  =   $(PROJ_BASE_DIR)/test/unit
PROJ_EXAMPLES_DIR   =   $(PROJ_BASE_DIR)/examples

XCONTRACT_INCLUDE_DIR	=	$(PROJ_BASE_DIR)/include
XCONTRACT_SRC_DIR		=	$(PROJ_SRC_DIR)/xcontract

XTESTS_INCLUDE_DIR	=	$(PROJ_BASE_DIR)/include
XTESTS_SRC_DIR		=	$(PROJ_SRC_DIR)/xtests

########################################
# tools

CC                  =   arcc.debug
CPPC                =   arcc.debug
AR                  =   lib
LD                  =   arld
RM                  =   del /f
CP                  =   copy
RC                  =   rc


# Includes

INCLUDES_LEADER     =   

PROJ_INCLUDE        =   -I$(PROJ_INCLUDE_DIR)

COMMON_INCLUDES     =

STLSOFT_1_10_INCLUDES	=	


STLSOFT_INCLUDES=\
 $(STLSOFT_1_10_INCLUDES)\
 -I"$(STLSOFT)/include"\
 \

CC_INCLUDES = \
 $(INCLUDES_LEADER)\
 $(PROJ_INCLUDE)\
 $(COMMON_INCLUDES)\
 $(STLSOFT_INCLUDES)\
 \

CPPC_INCLUDES = \
 $(INCLUDES_LEADER)\
 $(PROJ_INCLUDE)\
 $(COMMON_INCLUDES)\
 $(STLSOFT_INCLUDES)\
 \

3PTY_LIBS__D = \
 \

3PTY_LIBS__R = \
 \

3PTY_LIBS__T = \
 \


LIB_PATHS =\
 \

LIBS_LEADER         =   

COMMON_LIB          =

# #defines

CMN_DEF_NO_IMPLICIT_LINK =\
 -DB64_NO_IMPLICIT_LINK\
 \

COMMON_DEFS         =   -DWIN32\
	$(CMN_DEF_NO_IMPLICIT_LINK)\
	-DXCONTRACT_NO_USE_WIN_EVENTLOG

CC_DEFS             =   $(COMMON_DEFS)

CPPC_DEFS           =   $(COMMON_DEFS)

# Flags

CC_F_NOLOGO             =   
CC_F_COMPILE_ONLY       =   -c
CC_F_WARNING_MAX        =   --warning-level=maximum
CC_F_WARNING_ERR        =   --warnings-as-errors=on

CPPC_F_NOLOGO       =   
CPPC_F_COMPILE_ONLY =   -c
CPPC_F_WARNING_MAX  =   --warning-level=maximum
CPPC_F_WARNING_ERR  =   --warnings-as-errors=on
CPPC_F_EXCEPTIONS   =   --exception-handling=on
CPPC_F_RTTI         =   

CC_TOOL_FLAGS       =   --compilers=bc56,dmb-stlport,gcc34,icl8,vc6,vc71,ow12 --announce-tools --active-compiler=vc/6 --announce-tool-commands

CC_FLAGS = \
 $(CC_TOOL_FLAGS)\
 $(CC_F_NOLOGO)\
 $(CC_F_COMPILE_ONLY)\
 $(CC_F_WARNING_MAX)\
 $(CC_F_WARNING_ERR)\
 \

CPPC_TOOL_FLAGS     =   --compilers=bc56,dmb-stlport,gcc34,icl8,vc6,vc71 --announce-tools --active-compiler=vc/6 --announce-tool-commands

CPPC_FLAGS = \
 $(CPPC_TOOL_FLAGS)\
 $(CPPC_F_NOLOGO)\
 $(CPPC_F_COMPILE_ONLY)\
 $(CPPC_F_WARNING_MAX)\
 $(CPPC_F_WARNING_ERR)\
 $(CPPC_F_EXCEPTIONS)\
 $(CPPC_F_RTTI)\
 \

# Arguments (= flags + defines + includes)

CC_ARGS_D = -D_DEBUG -UB64_VARIANT_TEST -UNDEBUG $(CC_FLAGS) $(CC_DEFS) $(CC_INCLUDES)
CC_ARGS_R = -U_DEBUG -UB64_VARIANT_TEST -DNDEBUG $(CC_FLAGS) $(CC_DEFS) $(CC_INCLUDES)
CC_ARGS_T = -U_DEBUG -DB64_VARIANT_TEST -UNDEBUG $(CC_FLAGS) $(CC_DEFS) $(CC_INCLUDES)

CPPC_ARGS_D = -D_DEBUG -UB64_VARIANT_TEST -UNDEBUG $(CPPC_FLAGS) $(CPPC_DEFS) $(CPPC_INCLUDES)
CPPC_ARGS_R = -U_DEBUG -UB64_VARIANT_TEST -DNDEBUG $(CPPC_FLAGS) $(CPPC_DEFS) $(CPPC_INCLUDES)
CPPC_ARGS_T = -U_DEBUG -DB64_VARIANT_TEST -UNDEBUG $(CPPC_FLAGS) $(CPPC_DEFS) $(CPPC_INCLUDES)

# archiver flags and arguments

AR_ARGS = 

# linker flags and arguments

LD_F_NOLOGO = 

LD_ARGS_D =  $(LD_F_NOLOGO) $(LIB_PATHS) $(COMMON_LIB)
LD_ARGS_R =  $(LD_F_NOLOGO) $(LIB_PATHS) $(COMMON_LIB)
LD_ARGS_T =  $(LD_F_NOLOGO) $(LIB_PATHS) $(COMMON_LIB)

################################################################################
# Targets

########################################
# TARGETS_LIBS

TARGET_LIB_D = $(PROJ_LIB_DIR)/$(PROJ_TAG_VER).$(COMP_TAG).debug.lib
TARGET_LIB_UN_D = $(PROJ_LIB_DIR)/$(PROJ_TAG).$(COMP_TAG).debug.lib
LIB_CORE_D = $(PROJ_LIB_DIR)/$(PROJ_TAG_VER).core.$(COMP_TAG).debug.lib
LIB_CORE_UN_D = $(PROJ_LIB_DIR)/$(PROJ_TAG).core.$(COMP_TAG).debug.lib
TARGET_LIB_R = $(PROJ_LIB_DIR)/$(PROJ_TAG_VER).$(COMP_TAG).release.lib
TARGET_LIB_UN_R = $(PROJ_LIB_DIR)/$(PROJ_TAG).$(COMP_TAG).release.lib
LIB_CORE_R = $(PROJ_LIB_DIR)/$(PROJ_TAG_VER).core.$(COMP_TAG).release.lib
LIB_CORE_UN_R = $(PROJ_LIB_DIR)/$(PROJ_TAG).core.$(COMP_TAG).release.lib
TARGET_LIB_T = $(PROJ_LIB_DIR)/$(PROJ_TAG_VER).$(COMP_TAG).test.lib
TARGET_LIB_UN_T = $(PROJ_LIB_DIR)/$(PROJ_TAG).$(COMP_TAG).test.lib
LIB_CORE_T = $(PROJ_LIB_DIR)/$(PROJ_TAG_VER).core.$(COMP_TAG).test.lib
LIB_CORE_UN_T = $(PROJ_LIB_DIR)/$(PROJ_TAG).core.$(COMP_TAG).test.lib

# xcontract - unit/component testing only
LIB_XCONTRACT_D		=	$(PROJ_LIB_DIR)/$(XCONTRACT_TAG_VER).b64.$(COMP_ARCH_OS_TAG).debug.lib
LIB_XCONTRACT_R		=	$(PROJ_LIB_DIR)/$(XCONTRACT_TAG_VER).b64.$(COMP_ARCH_OS_TAG).release.lib
LIB_XCONTRACT_T		=	$(PROJ_LIB_DIR)/$(XCONTRACT_TAG_VER).b64.$(COMP_ARCH_OS_TAG).test.lib

# xtests - unit/component testing only
LIB_XTESTS_D			=	$(PROJ_LIB_DIR)/$(XTESTS_TAG_VER).b64.$(COMP_ARCH_OS_TAG).debug.lib
LIB_XTESTS_R			=	$(PROJ_LIB_DIR)/$(XTESTS_TAG_VER).b64.$(COMP_ARCH_OS_TAG).release.lib
LIB_XTESTS_T			=	$(PROJ_LIB_DIR)/$(XTESTS_TAG_VER).b64.$(COMP_ARCH_OS_TAG).test.lib

########################################
# TARGETS_TESTS

TARG_TEST_UNIT_C_1_D = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).test.c.1.$(COMP_TAG).debug.exe
TARG_TEST_UNIT_CPP_1_D = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).test.cpp.1.$(COMP_TAG).debug.exe
TARG_TEST_UNIT_C_1_R = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).test.c.1.$(COMP_TAG).release.exe
TARG_TEST_UNIT_CPP_1_R = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).test.cpp.1.$(COMP_TAG).release.exe
TARG_TEST_UNIT_C_1_T = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).test.c.1.$(COMP_TAG).test.exe
TARG_TEST_UNIT_CPP_1_T = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).test.cpp.1.$(COMP_TAG).test.exe

TARGETS_TESTS =\
 $(TARG_TEST_UNIT_C_1_D)\
 $(TARG_TEST_UNIT_C_1_R)\
 $(TARG_TEST_UNIT_C_1_T)\
 $(TARG_TEST_UNIT_CPP_1_D)\
 $(TARG_TEST_UNIT_CPP_1_R)\
 $(TARG_TEST_UNIT_CPP_1_T)\
 \

########################################
# TARGETS_EXAMPLES

TARG_EX_C_1_D = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.c.1.$(COMP_TAG).debug.exe
TARG_EX_C_2_D = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.c.2.$(COMP_TAG).debug.exe
TARG_EX_CPP_1_D = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.cpp.1.$(COMP_TAG).debug.exe
TARG_EX_CPP_2_D = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.cpp.2.$(COMP_TAG).debug.exe
TARG_EX_C_1_R = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.c.1.$(COMP_TAG).release.exe
TARG_EX_C_2_R = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.c.2.$(COMP_TAG).release.exe
TARG_EX_CPP_1_R = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.cpp.1.$(COMP_TAG).release.exe
TARG_EX_CPP_2_R = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.cpp.2.$(COMP_TAG).release.exe
TARG_EX_C_1_T = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.c.1.$(COMP_TAG).test.exe
TARG_EX_C_2_T = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.c.2.$(COMP_TAG).test.exe
TARG_EX_CPP_1_T = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.cpp.1.$(COMP_TAG).test.exe
TARG_EX_CPP_2_T = $(PROJ_BIN_DIR)/$(PROJ_TAG_VER).example.cpp.2.$(COMP_TAG).test.exe

TARGETS_EXAMPLES =\
 $(TARG_EX_C_1_D)\
 $(TARG_EX_C_2_D)\
 $(TARG_EX_C_1_R)\
 $(TARG_EX_C_2_R)\
 $(TARG_EX_C_1_T)\
 $(TARG_EX_C_2_T)\
 $(TARG_EX_CPP_1_D)\
 $(TARG_EX_CPP_2_D)\
 $(TARG_EX_CPP_1_R)\
 $(TARG_EX_CPP_2_R)\
 $(TARG_EX_CPP_1_T)\
 $(TARG_EX_CPP_2_T)\
 \

########################################
# Targets

TARGETS_LIB_CORE = \
 $(TARGET_LIB_D)\
 $(TARGET_LIB_UN_D)\
 $(LIB_CORE_D)\
 $(LIB_CORE_UN_D)\
 $(TARGET_LIB_R)\
 $(TARGET_LIB_UN_R)\
 $(LIB_CORE_R)\
 $(LIB_CORE_UN_R)\
 $(TARGET_LIB_T)\
 $(TARGET_LIB_UN_T)\
 $(LIB_CORE_T)\
 $(LIB_CORE_UN_T)\
 \

TARGETS_LIB_XCONTRACT = \
 $(LIB_XCONTRACT_D)\
 $(LIB_XCONTRACT_R)\
 $(LIB_XCONTRACT_T)\
 \

TARGETS_LIB_XTESTS = \
 $(LIB_XTESTS_D)\
 $(LIB_XTESTS_R)\
 $(LIB_XTESTS_T)\
 \

########################################
# OBJS_CORE

OBJ_CORE_b64__D = $(PROJ_OBJ_DIR)/b64.debug.obj

OBJS_CORE_D = \
 $(OBJ_CORE_b64__D)\
 \


OBJ_CORE_b64__R = $(PROJ_OBJ_DIR)/b64.release.obj

OBJS_CORE_R = \
 $(OBJ_CORE_b64__R)\
 \


OBJ_CORE_b64__T = $(PROJ_OBJ_DIR)/b64.test.obj

OBJS_CORE_T = \
 $(OBJ_CORE_b64__T)\
 \


OBJS_CORE =\
 $(OBJS_CORE_D)\
 $(OBJS_CORE_R)\
 $(OBJS_CORE_T)\
 \

# xcontract

OBJ_xct_core_D           = $(PROJ_OBJ_DIR)/xcontract.b64.core.api.debug.obj

OBJ_XCONTRACT_D = \
 $(OBJ_xct_core_D)\
 \

OBJ_xct_core_R           = $(PROJ_OBJ_DIR)/xcontract.b64.core.api.release.obj

OBJ_XCONTRACT_R = \
 $(OBJ_xct_core_R)\
 \

OBJ_xct_core_T           = $(PROJ_OBJ_DIR)/xcontract.b64.core.api.test.obj

OBJ_XCONTRACT_T = \
 $(OBJ_xct_core_T)\
 \


# xtests

OBJ_xts_core_D           = $(PROJ_OBJ_DIR)/xtests.b64.core.api.debug.obj

OBJ_XTESTS_D = \
 $(OBJ_xts_core_D)\
 \

OBJ_xts_core_R           = $(PROJ_OBJ_DIR)/xtests.b64.core.api.release.obj

OBJ_XTESTS_R = \
 $(OBJ_xts_core_R)\
 \

OBJ_xts_core_T           = $(PROJ_OBJ_DIR)/xtests.b64.core.api.test.obj

OBJ_XTESTS_T = \
 $(OBJ_xts_core_T)\
 \


########################################
# OBJS_LIBS

OBJS_LIBS =\
 $(OBJS_CORE)\
 \

########################################
# OBJS_TESTS

OBJ_TEST_UNIT_C_1_D = $(PROJ_OBJ_DIR)/test.c.1.debug.obj
OBJ_TEST_UNIT_CPP_1_D = $(PROJ_OBJ_DIR)/test.cpp.1.debug.obj
OBJ_TEST_UNIT_C_1_R = $(PROJ_OBJ_DIR)/test.c.1.release.obj
OBJ_TEST_UNIT_CPP_1_R = $(PROJ_OBJ_DIR)/test.cpp.1.release.obj
OBJ_TEST_UNIT_C_1_T = $(PROJ_OBJ_DIR)/test.c.1.test.obj
OBJ_TEST_UNIT_CPP_1_T = $(PROJ_OBJ_DIR)/test.cpp.1.test.obj

OBJS_TESTS =\
 $(OBJ_TEST_UNIT_C_1_D)\
 $(OBJ_TEST_UNIT_C_1_R)\
 $(OBJ_TEST_UNIT_C_1_T)\
 $(OBJ_TEST_UNIT_CPP_1_D)\
 $(OBJ_TEST_UNIT_CPP_1_R)\
 $(OBJ_TEST_UNIT_CPP_1_T)\
 \

########################################
# OBJS_EXAMPLES

OBJ_EX_C_1_D = $(PROJ_OBJ_DIR)/example.c.1.debug.obj
OBJ_EX_C_2_D = $(PROJ_OBJ_DIR)/example.c.2.debug.obj
OBJ_EX_CPP_1_D = $(PROJ_OBJ_DIR)/example.cpp.1.debug.obj
OBJ_EX_CPP_2_D = $(PROJ_OBJ_DIR)/example.cpp.2.debug.obj
OBJ_EX_C_1_R = $(PROJ_OBJ_DIR)/example.c.1.release.obj
OBJ_EX_C_2_R = $(PROJ_OBJ_DIR)/example.c.2.release.obj
OBJ_EX_CPP_1_R = $(PROJ_OBJ_DIR)/example.cpp.1.release.obj
OBJ_EX_CPP_2_R = $(PROJ_OBJ_DIR)/example.cpp.2.release.obj
OBJ_EX_C_1_T = $(PROJ_OBJ_DIR)/example.c.1.test.obj
OBJ_EX_C_2_T = $(PROJ_OBJ_DIR)/example.c.2.test.obj
OBJ_EX_CPP_1_T = $(PROJ_OBJ_DIR)/example.cpp.1.test.obj
OBJ_EX_CPP_2_T = $(PROJ_OBJ_DIR)/example.cpp.2.test.obj

OBJS_EXAMPLES =\
 $(OBJ_EX_C_1_D)\
 $(OBJ_EX_C_2_D)\
 $(OBJ_EX_C_1_R)\
 $(OBJ_EX_C_2_R)\
 $(OBJ_EX_C_1_T)\
 $(OBJ_EX_C_2_T)\
 $(OBJ_EX_CPP_1_D)\
 $(OBJ_EX_CPP_2_D)\
 $(OBJ_EX_CPP_1_R)\
 $(OBJ_EX_CPP_2_R)\
 $(OBJ_EX_CPP_1_T)\
 $(OBJ_EX_CPP_2_T)\
 \


########################################
# OBJS

OBJS = \
 $(OBJS_LIBS)\
 $(OBJS_TESTS)\
 $(OBJS_EXAMPLES)\
 \

########################################
# SRC_CORE

LIB_CORE_b64_H = $(PROJ_INCLUDE_DIR)/b64/b64.h
LIB_CORE_implicit_link_H = $(PROJ_INCLUDE_DIR)/b64/implicit_link.h

LIB_CORE_b64_C = $(PROJ_SRC_DIR)/b64.c

########################################
# SRC_HEADERS

HEADERS_LIB	= \
 $(LIB_CORE_b64_H)\
 $(LIB_CORE_implicit_link_H)\
 \

HEADERS_MAPPING_CPP =\
 $(PROJ_INCLUDE_DIR)/b64/b64.hpp\
 \

########################################
# SRC_TESTS

SRC_TEST_C_1 = $(PROJ_TEST_UNIT_DIR)/test.unit.c.1/test.unit.c.1.c
SRC_TEST_CPP_1 = $(PROJ_TEST_UNIT_DIR)/test.unit.cpp.1/test.unit.cpp.1.cpp

########################################
# SRC_EXAMPLES

SRC_EXAMPLE_C_1 = $(PROJ_EXAMPLES_DIR)/c/example.c.1/example.c.1.c
SRC_EXAMPLE_C_2 = $(PROJ_EXAMPLES_DIR)/c/example.c.2/example.c.2.c
SRC_EXAMPLE_CPP_1 = $(PROJ_EXAMPLES_DIR)/cpp/example.cpp.1/example.cpp.1.cpp
SRC_EXAMPLE_CPP_2 = $(PROJ_EXAMPLES_DIR)/cpp/example.cpp.2/example.cpp.2.cpp

################################################################################
# External Targets (invokable by make)

all:                build.all

check.lib.dir:
	@-if not exist $(PROJ_LIB_DIR)\NUL mkdir $(PROJ_LIB_DIR)

check.bin.dir:
	@-if not exist $(PROJ_BIN_DIR)\NUL mkdir $(PROJ_BIN_DIR)

build.libs.core:		check.lib.dir $(TARGETS_LIB_CORE)

build.libs.xcontract:	check.lib.dir $(TARGETS_LIB_XCONTRACT)

build.libs.xtests:		check.lib.dir $(TARGETS_LIB_XTESTS)

build.libs:				build.libs.core build.libs.xcontract build.libs.xtests

build.examples:			check.bin.dir $(TARGETS_EXAMPLES)

build.test:				check.bin.dir $(TARGETS_TESTS)

build.all:				build.libs build.examples build.test

build:					build.all

compile:				$(OBJS)


test_:              makefile
	$(TARG_TEST_UNIT_C_1_D) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_C_1_R) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_C_1_T) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_CPP_1_D) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_CPP_1_R) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_CPP_1_T) --verbosity=$(TEST_VERBOSITY)

test:		build.test test_

test.all:	test

test-all:	test.all	# backwards compatibility


usage:
	@echo "make [usage | clean | compile | all | test | target] - default is 'all'"

clean.libs.core:
	@-$(RM) $(TARGET_LIB_D)
	@-$(RM) $(TARGET_LIB_UN_D)
	@-$(RM) $(LIB_CORE_D)
	@-$(RM) $(LIB_CORE_UN_D)
	@-$(RM) $(TARGET_LIB_R)
	@-$(RM) $(TARGET_LIB_UN_R)
	@-$(RM) $(LIB_CORE_R)
	@-$(RM) $(LIB_CORE_UN_R)
	@-$(RM) $(TARGET_LIB_T)
	@-$(RM) $(TARGET_LIB_UN_T)
	@-$(RM) $(LIB_CORE_T)
	@-$(RM) $(LIB_CORE_UN_T)
	@-$(RM) $(OBJ_CORE_b64__D)
	@-$(RM) $(OBJ_CORE_b64__R)
	@-$(RM) $(OBJ_CORE_b64__T)

clean.libs.xcontract:	makefile

clean.libs.xtests:		makefile

clean.libs:				clean.libs.core clean.libs.xcontract clean.libs.xtests

clean.test:
	@-$(RM) $(TARG_TEST_UNIT_C_1_D)
	@-$(RM) $(TARG_TEST_UNIT_C_1_R)
	@-$(RM) $(TARG_TEST_UNIT_C_1_T)
	@-$(RM) $(TARG_TEST_UNIT_CPP_1_D)
	@-$(RM) $(TARG_TEST_UNIT_CPP_1_R)
	@-$(RM) $(TARG_TEST_UNIT_CPP_1_T)
	@-$(RM) $(OBJ_TEST_UNIT_C_1_D)
	@-$(RM) $(OBJ_TEST_UNIT_C_1_R)
	@-$(RM) $(OBJ_TEST_UNIT_C_1_T)
	@-$(RM) $(OBJ_TEST_UNIT_CPP_1_D)
	@-$(RM) $(OBJ_TEST_UNIT_CPP_1_R)
	@-$(RM) $(OBJ_TEST_UNIT_CPP_1_T)

clean.examples:
	@-$(RM) $(TARG_EX_C_1_D)
	@-$(RM) $(TARG_EX_C_2_D)
	@-$(RM) $(TARG_EX_C_1_R)
	@-$(RM) $(TARG_EX_C_2_R)
	@-$(RM) $(TARG_EX_C_1_T)
	@-$(RM) $(TARG_EX_C_2_T)
	@-$(RM) $(TARG_EX_CPP_1_D)
	@-$(RM) $(TARG_EX_CPP_2_D)
	@-$(RM) $(TARG_EX_CPP_1_R)
	@-$(RM) $(TARG_EX_CPP_2_R)
	@-$(RM) $(TARG_EX_CPP_1_T)
	@-$(RM) $(TARG_EX_CPP_2_T)
	@-$(RM) $(OBJ_EX_C_1_D)
	@-$(RM) $(OBJ_EX_C_2_D)
	@-$(RM) $(OBJ_EX_C_1_R)
	@-$(RM) $(OBJ_EX_C_2_R)
	@-$(RM) $(OBJ_EX_C_1_T)
	@-$(RM) $(OBJ_EX_C_2_T)
	@-$(RM) $(OBJ_EX_CPP_1_D)
	@-$(RM) $(OBJ_EX_CPP_2_D)
	@-$(RM) $(OBJ_EX_CPP_1_R)
	@-$(RM) $(OBJ_EX_CPP_2_R)
	@-$(RM) $(OBJ_EX_CPP_1_T)
	@-$(RM) $(OBJ_EX_CPP_2_T)

clean: clean.libs clean.test clean.examples

################################################################################
# Internal Targets (invokable by make)

########################################
# Main library - 'b64.X.core.a'


$(TARGET_LIB_D):	$(LIB_CORE_D)
	@$(CP) $(LIB_CORE_D) $@

$(TARGET_LIB_UN_D):	$(LIB_CORE_D)
	@$(CP) $(LIB_CORE_D) $@

$(LIB_CORE_UN_D): $(LIB_CORE_D)
	@$(CP) $(LIB_CORE_D) $@

$(LIB_CORE_D): $(OBJS_CORE_D) makefile
	$(AR) $(AR_ARGS) $@ $(OBJS_CORE_D)

$(OBJ_CORE_b64__D): $(PROJ_SRC_DIR)/b64.c $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_D) --output-path=$@ $(PROJ_SRC_DIR)/b64.c


$(TARGET_LIB_R):	$(LIB_CORE_R)
	@$(CP) $(LIB_CORE_R) $@

$(TARGET_LIB_UN_R):	$(LIB_CORE_R)
	@$(CP) $(LIB_CORE_R) $@

$(LIB_CORE_UN_R): $(LIB_CORE_R)
	@$(CP) $(LIB_CORE_R) $@

$(LIB_CORE_R): $(OBJS_CORE_R) makefile
	$(AR) $(AR_ARGS) $@ $(OBJS_CORE_R)

$(OBJ_CORE_b64__R): $(PROJ_SRC_DIR)/b64.c $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_R) --output-path=$@ $(PROJ_SRC_DIR)/b64.c


$(TARGET_LIB_T):	$(LIB_CORE_T)
	@$(CP) $(LIB_CORE_T) $@

$(TARGET_LIB_UN_T):	$(LIB_CORE_T)
	@$(CP) $(LIB_CORE_T) $@

$(LIB_CORE_UN_T): $(LIB_CORE_T)
	@$(CP) $(LIB_CORE_T) $@

$(LIB_CORE_T): $(OBJS_CORE_T) makefile
	$(AR) $(AR_ARGS) $@ $(OBJS_CORE_T)

$(OBJ_CORE_b64__T): $(PROJ_SRC_DIR)/b64.c $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_T) --output-path=$@ $(PROJ_SRC_DIR)/b64.c


########################################
# xContract library - 'xcontract.X.b64.core.a'

$(LIB_XCONTRACT_D): $(OBJ_XCONTRACT_D) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)/NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) $@ $(OBJ_XCONTRACT_D)

$(OBJ_xct_core_D): $(XCONTRACT_SRC_DIR)/xcontract.core.c $(XCONTRACT_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CC) $(CC_ARGS_D) $(CC_XCONTRACT_INCLUDES_D) --output-path=$@ $(XCONTRACT_SRC_DIR)/xcontract.core.c

$(LIB_XCONTRACT_R): $(OBJ_XCONTRACT_R) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)/NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) $@ $(OBJ_XCONTRACT_R)

$(OBJ_xct_core_R): $(XCONTRACT_SRC_DIR)/xcontract.core.c $(XCONTRACT_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CC) $(CC_ARGS_R) $(CC_XCONTRACT_INCLUDES_R) --output-path=$@ $(XCONTRACT_SRC_DIR)/xcontract.core.c

$(LIB_XCONTRACT_T): $(OBJ_XCONTRACT_T) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)/NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) $@ $(OBJ_XCONTRACT_T)

$(OBJ_xct_core_T): $(XCONTRACT_SRC_DIR)/xcontract.core.c $(XCONTRACT_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CC) $(CC_ARGS_T) $(CC_XCONTRACT_INCLUDES_T) --output-path=$@ $(XCONTRACT_SRC_DIR)/xcontract.core.c

########################################
# xTests library - 'xtests.X.b64.core.a'

$(LIB_XTESTS_D): $(OBJ_XTESTS_D) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)/NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) $@ $(OBJ_XTESTS_D)

$(OBJ_xts_core_D): $(XTESTS_SRC_DIR)/xtests.core.cpp $(XTESTS_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CPPC) $(CPPC_ARGS_D) $(CPPC_XTESTS_INCLUDES_D) --output-path=$@ $(XTESTS_SRC_DIR)/xtests.core.cpp

$(LIB_XTESTS_R): $(OBJ_XTESTS_R) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)/NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) $@ $(OBJ_XTESTS_R)

$(OBJ_xts_core_R): $(XTESTS_SRC_DIR)/xtests.core.cpp $(XTESTS_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CPPC) $(CPPC_ARGS_R) $(CPPC_XTESTS_INCLUDES_R) --output-path=$@ $(XTESTS_SRC_DIR)/xtests.core.cpp

$(LIB_XTESTS_T): $(OBJ_XTESTS_T) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)/NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) $@ $(OBJ_XTESTS_T)

$(OBJ_xts_core_T): $(XTESTS_SRC_DIR)/xtests.core.cpp $(XTESTS_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CPPC) $(CPPC_ARGS_T) $(CPPC_XTESTS_INCLUDES_T) --output-path=$@ $(XTESTS_SRC_DIR)/xtests.core.cpp

########################################
# test

# C test client

$(TARG_TEST_UNIT_C_1_D):	\
	$(OBJ_TEST_UNIT_C_1_D)\
	$(TARGET_LIB_D)\
	$(LIB_XCONTRACT_D)\
	$(LIB_XTESTS_D)\
	makefile
	$(LD) $(LD_ARGS_D) $@ $(OBJ_TEST_UNIT_C_1_D) $(TARGET_LIB_D) $(LIB_XTESTS_D) $(LIB_XCONTRACT_D) $(3PTY_LIBS__D)

$(OBJ_TEST_UNIT_C_1_D):	$(SRC_TEST_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_D) --output-path=$@ $(SRC_TEST_C_1)

$(TARG_TEST_UNIT_C_1_R):	\
	$(OBJ_TEST_UNIT_C_1_R)\
	$(TARGET_LIB_R)\
	$(LIB_XCONTRACT_R)\
	$(LIB_XTESTS_R)\
	makefile
	$(LD) $(LD_ARGS_R) $@ $(OBJ_TEST_UNIT_C_1_R) $(TARGET_LIB_R) $(LIB_XTESTS_R) $(LIB_XCONTRACT_R) $(3PTY_LIBS__R)

$(OBJ_TEST_UNIT_C_1_R):	$(SRC_TEST_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_R) --output-path=$@ $(SRC_TEST_C_1)

$(TARG_TEST_UNIT_C_1_T):	\
	$(OBJ_TEST_UNIT_C_1_T)\
	$(TARGET_LIB_T)\
	$(LIB_XCONTRACT_T)\
	$(LIB_XTESTS_T)\
	makefile
	$(LD) $(LD_ARGS_T) $@ $(OBJ_TEST_UNIT_C_1_T) $(TARGET_LIB_T) $(LIB_XTESTS_T) $(LIB_XCONTRACT_T) $(3PTY_LIBS__T)

$(OBJ_TEST_UNIT_C_1_T):	$(SRC_TEST_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_T) --output-path=$@ $(SRC_TEST_C_1)

# C++ test client

$(TARG_TEST_UNIT_CPP_1_D):	\
	$(OBJ_TEST_UNIT_CPP_1_D)\
	$(TARGET_LIB_D)\
	$(LIB_XCONTRACT_D)\
	$(LIB_XTESTS_D)\
	makefile
	$(LD) $(LD_ARGS_D) $@ $(OBJ_TEST_UNIT_CPP_1_D) $(TARGET_LIB_D) $(LIB_XTESTS_D) $(LIB_XCONTRACT_D) $(3PTY_LIBS__D)

$(OBJ_TEST_UNIT_CPP_1_D):	$(SRC_TEST_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_D) --output-path=$@ $(SRC_TEST_CPP_1)

$(TARG_TEST_UNIT_CPP_1_R):	\
	$(OBJ_TEST_UNIT_CPP_1_R)\
	$(TARGET_LIB_R)\
	$(LIB_XCONTRACT_R)\
	$(LIB_XTESTS_R)\
	makefile
	$(LD) $(LD_ARGS_R) $@ $(OBJ_TEST_UNIT_CPP_1_R) $(TARGET_LIB_R) $(LIB_XTESTS_R) $(LIB_XCONTRACT_R) $(3PTY_LIBS__R)

$(OBJ_TEST_UNIT_CPP_1_R):	$(SRC_TEST_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_R) --output-path=$@ $(SRC_TEST_CPP_1)

$(TARG_TEST_UNIT_CPP_1_T):	\
	$(OBJ_TEST_UNIT_CPP_1_T)\
	$(TARGET_LIB_T)\
	$(LIB_XCONTRACT_T)\
	$(LIB_XTESTS_T)\
	makefile
	$(LD) $(LD_ARGS_T) $@ $(OBJ_TEST_UNIT_CPP_1_T) $(TARGET_LIB_T) $(LIB_XTESTS_T) $(LIB_XCONTRACT_T) $(3PTY_LIBS__T)

$(OBJ_TEST_UNIT_CPP_1_T):	$(SRC_TEST_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_T) --output-path=$@ $(SRC_TEST_CPP_1)

########################################
# examples

# C example client

$(TARG_EX_C_1_D):	\
	$(OBJ_EX_C_1_D)\
	$(TARGET_LIB_D)\
	$(LIB_XCONTRACT_D)\
	makefile
	$(LD) $(LD_ARGS_D) $@ $(OBJ_EX_C_1_D) $(TARGET_LIB_D) $(LIB_XCONTRACT_D) $(3PTY_LIBS__D)

$(OBJ_EX_C_1_D):	$(SRC_EXAMPLE_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_D) --output-path=$@ $(SRC_EXAMPLE_C_1)

$(TARG_EX_C_2_D):	\
	$(OBJ_EX_C_2_D)\
	$(TARGET_LIB_D)\
	$(LIB_XCONTRACT_D)\
	makefile
	$(LD) $(LD_ARGS_D) $@ $(OBJ_EX_C_2_D) $(TARGET_LIB_D) $(LIB_XCONTRACT_D) $(3PTY_LIBS__D)

$(OBJ_EX_C_2_D):	$(SRC_EXAMPLE_C_2) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_D) --output-path=$@ $(SRC_EXAMPLE_C_2)

$(TARG_EX_C_1_R):	\
	$(OBJ_EX_C_1_R)\
	$(TARGET_LIB_R)\
	$(LIB_XCONTRACT_R)\
	makefile
	$(LD) $(LD_ARGS_R) $@ $(OBJ_EX_C_1_R) $(TARGET_LIB_R) $(LIB_XCONTRACT_R) $(3PTY_LIBS__R)

$(OBJ_EX_C_1_R):	$(SRC_EXAMPLE_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_R) --output-path=$@ $(SRC_EXAMPLE_C_1)

$(TARG_EX_C_2_R):	\
	$(OBJ_EX_C_2_R)\
	$(TARGET_LIB_R)\
	$(LIB_XCONTRACT_R)\
	makefile
	$(LD) $(LD_ARGS_R) $@ $(OBJ_EX_C_2_R) $(TARGET_LIB_R) $(LIB_XCONTRACT_R) $(3PTY_LIBS__R)

$(OBJ_EX_C_2_R):	$(SRC_EXAMPLE_C_2) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_R) --output-path=$@ $(SRC_EXAMPLE_C_2)

$(TARG_EX_C_1_T):	\
	$(OBJ_EX_C_1_T)\
	$(TARGET_LIB_T)\
	$(LIB_XCONTRACT_T)\
	makefile
	$(LD) $(LD_ARGS_T) $@ $(OBJ_EX_C_1_T) $(TARGET_LIB_T) $(LIB_XCONTRACT_T) $(3PTY_LIBS__T)

$(OBJ_EX_C_1_T):	$(SRC_EXAMPLE_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_T) --output-path=$@ $(SRC_EXAMPLE_C_1)

$(TARG_EX_C_2_T):	\
	$(OBJ_EX_C_2_T)\
	$(TARGET_LIB_T)\
	$(LIB_XCONTRACT_T)\
	makefile
	$(LD) $(LD_ARGS_T) $@ $(OBJ_EX_C_2_T) $(TARGET_LIB_T) $(LIB_XCONTRACT_T) $(3PTY_LIBS__T)

$(OBJ_EX_C_2_T):	$(SRC_EXAMPLE_C_2) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_T) --output-path=$@ $(SRC_EXAMPLE_C_2)

# C++ example client

$(TARG_EX_CPP_1_D):	\
	$(OBJ_EX_CPP_1_D)\
	$(TARGET_LIB_D)\
	$(LIB_XCONTRACT_D)\
	makefile
	$(LD) $(LD_ARGS_D) $@ $(OBJ_EX_CPP_1_D) $(TARGET_LIB_D) $(LIB_XCONTRACT_D) $(3PTY_LIBS__D)

$(OBJ_EX_CPP_1_D):	$(SRC_EXAMPLE_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_D) --output-path=$@ $(SRC_EXAMPLE_CPP_1)

$(TARG_EX_CPP_2_D):	\
	$(OBJ_EX_CPP_2_D)\
	$(TARGET_LIB_D)\
	$(LIB_XCONTRACT_D)\
	makefile
	$(LD) $(LD_ARGS_D) $@ $(OBJ_EX_CPP_2_D) $(TARGET_LIB_D) $(LIB_XCONTRACT_D) $(3PTY_LIBS__D)

$(OBJ_EX_CPP_2_D):	$(SRC_EXAMPLE_CPP_2) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_D) --output-path=$@ $(SRC_EXAMPLE_CPP_2)

$(TARG_EX_CPP_1_R):	\
	$(OBJ_EX_CPP_1_R)\
	$(TARGET_LIB_R)\
	$(LIB_XCONTRACT_R)\
	makefile
	$(LD) $(LD_ARGS_R) $@ $(OBJ_EX_CPP_1_R) $(TARGET_LIB_R) $(LIB_XCONTRACT_R) $(3PTY_LIBS__R)

$(OBJ_EX_CPP_1_R):	$(SRC_EXAMPLE_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_R) --output-path=$@ $(SRC_EXAMPLE_CPP_1)

$(TARG_EX_CPP_2_R):	\
	$(OBJ_EX_CPP_2_R)\
	$(TARGET_LIB_R)\
	$(LIB_XCONTRACT_R)\
	makefile
	$(LD) $(LD_ARGS_R) $@ $(OBJ_EX_CPP_2_R) $(TARGET_LIB_R) $(LIB_XCONTRACT_R) $(3PTY_LIBS__R)

$(OBJ_EX_CPP_2_R):	$(SRC_EXAMPLE_CPP_2) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_R) --output-path=$@ $(SRC_EXAMPLE_CPP_2)

$(TARG_EX_CPP_1_T):	\
	$(OBJ_EX_CPP_1_T)\
	$(TARGET_LIB_T)\
	$(LIB_XCONTRACT_T)\
	makefile
	$(LD) $(LD_ARGS_T) $@ $(OBJ_EX_CPP_1_T) $(TARGET_LIB_T) $(LIB_XCONTRACT_T) $(3PTY_LIBS__T)

$(OBJ_EX_CPP_1_T):	$(SRC_EXAMPLE_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_T) --output-path=$@ $(SRC_EXAMPLE_CPP_1)

$(TARG_EX_CPP_2_T):	\
	$(OBJ_EX_CPP_2_T)\
	$(TARGET_LIB_T)\
	$(LIB_XCONTRACT_T)\
	makefile
	$(LD) $(LD_ARGS_T) $@ $(OBJ_EX_CPP_2_T) $(TARGET_LIB_T) $(LIB_XCONTRACT_T) $(3PTY_LIBS__T)

$(OBJ_EX_CPP_2_T):	$(SRC_EXAMPLE_CPP_2) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_T) --output-path=$@ $(SRC_EXAMPLE_CPP_2)

################################################################################
