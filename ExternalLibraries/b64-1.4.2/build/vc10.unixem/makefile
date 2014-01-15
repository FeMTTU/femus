################################################################################
# File:     makefile
#
# Purpose:  Makefile for the b64 project, for the Visual C++ 10.0 (UNIXEm)
#           compiler, on UNIX
#
# Created:  15th Feburary 2005
# Updated:  4th February 2012
#
################################################################################


################################################################################
# Environment check
#

!ifndef STLSOFT
!error STLSOFT must be defined
!endif # STLSOFT

################################################################################
# Symbols

# Identification

COMP_TAG            =   vc10
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

PROJ_BASE_DIR       =   ..\..
PROJ_BIN_DIR        =   $(PROJ_BASE_DIR)\bin
PROJ_INCLUDE_DIR    =   $(PROJ_BASE_DIR)\include
PROJ_LIB_DIR        =   $(PROJ_BASE_DIR)\lib
PROJ_OBJ_DIR        =   .
PROJ_SRC_DIR        =   $(PROJ_BASE_DIR)\src
PROJ_MAPPINGS_DIR   =   $(PROJ_BASE_DIR)\mappings
PROJ_TEST_UNIT_DIR  =   $(PROJ_BASE_DIR)\test\unit
PROJ_EXAMPLES_DIR   =   $(PROJ_BASE_DIR)\examples

XCONTRACT_INCLUDE_DIR	=	$(PROJ_BASE_DIR)\include
XCONTRACT_SRC_DIR		=	$(PROJ_SRC_DIR)\xcontract

XTESTS_INCLUDE_DIR	=	$(PROJ_BASE_DIR)\include
XTESTS_SRC_DIR		=	$(PROJ_SRC_DIR)\xtests

########################################
# tools

CC                  =   cl
CPPC                =   cl
AR                  =   lib
LD                  =   link
RM                  =   del /f
CP                  =   copy
RC                  =   rc


# Includes

INCLUDES_LEADER     =   

PROJ_INCLUDE        =   -I$(PROJ_INCLUDE_DIR)

COMMON_INCLUDES     =

!ifdef STLSOFT_1_10
STLSOFT_1_10_INCLUDES	=	-I"$(STLSOFT_1_10)\include"
!else # STLSOFT_1_10
STLSOFT_1_10_INCLUDES	=	
!endif # STLSOFT_1_10


STLSOFT_INCLUDES=\
 $(STLSOFT_1_10_INCLUDES)\
 -I"$(STLSOFT)\include"\
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

3PTY_LIBS__DLL_T = \
 \

3PTY_LIBS__MT_T = \
 \

3PTY_LIBS__ST_T = \
 \

3PTY_LIBS__ST_UX_DEBUG = \
 \


LIB_PATHS =\
 \

LIBS_LEADER         =   

COMMON_LIB          =

# #defines

CMN_DEF_NO_IMPLICIT_LINK =\
 -DB64_NO_IMPLICIT_LINK\
 \

COMMON_DEFS         =   -DUNIX\
	$(CMN_DEF_NO_IMPLICIT_LINK)\
	-DXCONTRACT_NO_USE_WIN_EVENTLOG

CC_DEFS             =   $(COMMON_DEFS)

CPPC_DEFS           =   $(COMMON_DEFS)

# Flags

CC_F_NOLOGO             =   -nologo
CC_F_COMPILE_ONLY       =   -c
CC_F_WARNING_MAX        =   -W4                
CC_F_WARNING_ERR        =   -WX

CPPC_F_NOLOGO       =   -nologo
CPPC_F_COMPILE_ONLY =   -c
CPPC_F_WARNING_MAX  =   -W4 -wd4800 -wd4996
CPPC_F_WARNING_ERR  =   -WX
CPPC_F_EXCEPTIONS   =   -GX
CPPC_F_RTTI         =   -GR

CC_TOOL_FLAGS       =

CC_FLAGS = \
 $(CC_TOOL_FLAGS)\
 $(CC_F_NOLOGO)\
 $(CC_F_COMPILE_ONLY)\
 $(CC_F_WARNING_MAX)\
 $(CC_F_WARNING_ERR)\
 \

CPPC_TOOL_FLAGS     =

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

CC_ARGS_DLL_T = -MD -U_DEBUG -DB64_VARIANT_TEST -DNDEBUG $(CC_FLAGS) $(CC_DEFS) $(CC_INCLUDES)
CC_ARGS_MT_T = -MT -U_DEBUG -DB64_VARIANT_TEST -DNDEBUG $(CC_FLAGS) $(CC_DEFS) $(CC_INCLUDES)
CC_ARGS_ST_T = -ML -U_DEBUG -DB64_VARIANT_TEST -DNDEBUG $(CC_FLAGS) $(CC_DEFS) $(CC_INCLUDES)
CC_ARGS_ST_UX_DEBUG = -MLd -D_DEBUG -D_STLSOFT_FORCE_ANY_COMPILER -Zi -Fd$*.pdb $(CC_FLAGS) $(CC_DEFS) $(CC_INCLUDES)

CPPC_ARGS_DLL_T = -MD -U_DEBUG -DB64_VARIANT_TEST -DNDEBUG $(CPPC_FLAGS) $(CPPC_DEFS) $(CPPC_INCLUDES)
CPPC_ARGS_MT_T = -MT -U_DEBUG -DB64_VARIANT_TEST -DNDEBUG $(CPPC_FLAGS) $(CPPC_DEFS) $(CPPC_INCLUDES)
CPPC_ARGS_ST_T = -ML -U_DEBUG -DB64_VARIANT_TEST -DNDEBUG $(CPPC_FLAGS) $(CPPC_DEFS) $(CPPC_INCLUDES)
CPPC_ARGS_ST_UX_DEBUG = -MLd -D_DEBUG -D_STLSOFT_FORCE_ANY_COMPILER -Zi -Fd$*.pdb $(CPPC_FLAGS) $(CPPC_DEFS) $(CPPC_INCLUDES)

# archiver flags and arguments

AR_ARGS = -nologo

# linker flags and arguments

LD_F_NOLOGO = -nologo

LD_ARGS_DLL_T =  $(LD_F_NOLOGO) $(LIB_PATHS) $(COMMON_LIB)
LD_ARGS_MT_T =  $(LD_F_NOLOGO) $(LIB_PATHS) $(COMMON_LIB)
LD_ARGS_ST_T =  $(LD_F_NOLOGO) $(LIB_PATHS) $(COMMON_LIB)
LD_ARGS_ST_UX_DEBUG = libWinSysLog.$(COMP_TAG).debug.lib unixem.1.$(COMP_TAG).debug.lib -debug $(LD_F_NOLOGO) $(LIB_PATHS) $(COMMON_LIB)

################################################################################
# Targets

########################################
# TARGETS_LIBS

TARGET_LIB_DLL_T = $(PROJ_LIB_DIR)\$(PROJ_TAG_VER).$(COMP_TAG).dll.test.lib
TARGET_LIB_UN_DLL_T = $(PROJ_LIB_DIR)\$(PROJ_TAG).$(COMP_TAG).dll.test.lib
LIB_CORE_DLL_T = $(PROJ_LIB_DIR)\$(PROJ_TAG_VER).core.$(COMP_TAG).dll.test.lib
LIB_CORE_UN_DLL_T = $(PROJ_LIB_DIR)\$(PROJ_TAG).core.$(COMP_TAG).dll.test.lib
TARGET_LIB_MT_T = $(PROJ_LIB_DIR)\$(PROJ_TAG_VER).$(COMP_TAG).mt.test.lib
TARGET_LIB_UN_MT_T = $(PROJ_LIB_DIR)\$(PROJ_TAG).$(COMP_TAG).mt.test.lib
LIB_CORE_MT_T = $(PROJ_LIB_DIR)\$(PROJ_TAG_VER).core.$(COMP_TAG).mt.test.lib
LIB_CORE_UN_MT_T = $(PROJ_LIB_DIR)\$(PROJ_TAG).core.$(COMP_TAG).mt.test.lib
TARGET_LIB_ST_T = $(PROJ_LIB_DIR)\$(PROJ_TAG_VER).$(COMP_TAG).test.lib
TARGET_LIB_UN_ST_T = $(PROJ_LIB_DIR)\$(PROJ_TAG).$(COMP_TAG).test.lib
LIB_CORE_ST_T = $(PROJ_LIB_DIR)\$(PROJ_TAG_VER).core.$(COMP_TAG).test.lib
LIB_CORE_UN_ST_T = $(PROJ_LIB_DIR)\$(PROJ_TAG).core.$(COMP_TAG).test.lib
TARGET_LIB_ST_UX_DEBUG = $(PROJ_LIB_DIR)\$(PROJ_TAG_VER).$(COMP_TAG).unix.debug.lib
TARGET_LIB_UN_ST_UX_DEBUG = $(PROJ_LIB_DIR)\$(PROJ_TAG).$(COMP_TAG).unix.debug.lib
LIB_CORE_ST_UX_DEBUG = $(PROJ_LIB_DIR)\$(PROJ_TAG_VER).core.$(COMP_TAG).unix.debug.lib
LIB_CORE_UN_ST_UX_DEBUG = $(PROJ_LIB_DIR)\$(PROJ_TAG).core.$(COMP_TAG).unix.debug.lib
TARGET_PDB_ST_UX_DEBUG = $(PROJ_LIB_DIR)\$(PROJ_TAG_VER).$(COMP_TAG).unix.debug.pdb

# xcontract - unit/component testing only
LIB_XCONTRACT_DLL_T		=	$(PROJ_LIB_DIR)\$(XCONTRACT_TAG_VER).b64.$(COMP_ARCH_OS_TAG).dll.test.lib
LIB_XCONTRACT_MT_T		=	$(PROJ_LIB_DIR)\$(XCONTRACT_TAG_VER).b64.$(COMP_ARCH_OS_TAG).mt.test.lib
LIB_XCONTRACT_ST_T		=	$(PROJ_LIB_DIR)\$(XCONTRACT_TAG_VER).b64.$(COMP_ARCH_OS_TAG).test.lib
LIB_XCONTRACT_ST_UX_DEBUG		=	$(PROJ_LIB_DIR)\$(XCONTRACT_TAG_VER).b64.$(COMP_ARCH_OS_TAG).unix.debug.lib

# xtests - unit/component testing only
LIB_XTESTS_DLL_T			=	$(PROJ_LIB_DIR)\$(XTESTS_TAG_VER).b64.$(COMP_ARCH_OS_TAG).dll.test.lib
LIB_XTESTS_MT_T			=	$(PROJ_LIB_DIR)\$(XTESTS_TAG_VER).b64.$(COMP_ARCH_OS_TAG).mt.test.lib
LIB_XTESTS_ST_T			=	$(PROJ_LIB_DIR)\$(XTESTS_TAG_VER).b64.$(COMP_ARCH_OS_TAG).test.lib
LIB_XTESTS_ST_UX_DEBUG			=	$(PROJ_LIB_DIR)\$(XTESTS_TAG_VER).b64.$(COMP_ARCH_OS_TAG).unix.debug.lib

########################################
# TARGETS_TESTS

TARG_TEST_UNIT_C_1_DLL_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).test.c.1.$(COMP_TAG).dll.test.exe
TARG_TEST_UNIT_CPP_1_DLL_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).test.cpp.1.$(COMP_TAG).dll.test.exe
TARG_TEST_UNIT_C_1_MT_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).test.c.1.$(COMP_TAG).mt.test.exe
TARG_TEST_UNIT_CPP_1_MT_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).test.cpp.1.$(COMP_TAG).mt.test.exe
TARG_TEST_UNIT_C_1_ST_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).test.c.1.$(COMP_TAG).test.exe
TARG_TEST_UNIT_CPP_1_ST_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).test.cpp.1.$(COMP_TAG).test.exe
TARG_TEST_UNIT_C_1_ST_UX_DEBUG = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).test.c.1.$(COMP_TAG).unix.debug.exe
TARG_TEST_UNIT_CPP_1_ST_UX_DEBUG = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).test.cpp.1.$(COMP_TAG).unix.debug.exe

TARGETS_TESTS =\
 $(TARG_TEST_UNIT_C_1_DLL_T)\
 $(TARG_TEST_UNIT_C_1_MT_T)\
 $(TARG_TEST_UNIT_C_1_ST_T)\
 $(TARG_TEST_UNIT_C_1_ST_UX_DEBUG)\
 $(TARG_TEST_UNIT_CPP_1_DLL_T)\
 $(TARG_TEST_UNIT_CPP_1_MT_T)\
 $(TARG_TEST_UNIT_CPP_1_ST_T)\
 $(TARG_TEST_UNIT_CPP_1_ST_UX_DEBUG)\
 \

########################################
# TARGETS_EXAMPLES

TARG_EX_C_1_DLL_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.c.1.$(COMP_TAG).dll.test.exe
TARG_EX_C_2_DLL_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.c.2.$(COMP_TAG).dll.test.exe
TARG_EX_CPP_1_DLL_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.cpp.1.$(COMP_TAG).dll.test.exe
TARG_EX_CPP_2_DLL_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.cpp.2.$(COMP_TAG).dll.test.exe
TARG_EX_C_1_MT_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.c.1.$(COMP_TAG).mt.test.exe
TARG_EX_C_2_MT_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.c.2.$(COMP_TAG).mt.test.exe
TARG_EX_CPP_1_MT_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.cpp.1.$(COMP_TAG).mt.test.exe
TARG_EX_CPP_2_MT_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.cpp.2.$(COMP_TAG).mt.test.exe
TARG_EX_C_1_ST_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.c.1.$(COMP_TAG).test.exe
TARG_EX_C_2_ST_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.c.2.$(COMP_TAG).test.exe
TARG_EX_CPP_1_ST_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.cpp.1.$(COMP_TAG).test.exe
TARG_EX_CPP_2_ST_T = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.cpp.2.$(COMP_TAG).test.exe
TARG_EX_C_1_ST_UX_DEBUG = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.c.1.$(COMP_TAG).unix.debug.exe
TARG_EX_C_2_ST_UX_DEBUG = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.c.2.$(COMP_TAG).unix.debug.exe
TARG_EX_CPP_1_ST_UX_DEBUG = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.cpp.1.$(COMP_TAG).unix.debug.exe
TARG_EX_CPP_2_ST_UX_DEBUG = $(PROJ_BIN_DIR)\$(PROJ_TAG_VER).example.cpp.2.$(COMP_TAG).unix.debug.exe

TARGETS_EXAMPLES =\
 $(TARG_EX_C_1_DLL_T)\
 $(TARG_EX_C_2_DLL_T)\
 $(TARG_EX_C_1_MT_T)\
 $(TARG_EX_C_2_MT_T)\
 $(TARG_EX_C_1_ST_T)\
 $(TARG_EX_C_2_ST_T)\
 $(TARG_EX_C_1_ST_UX_DEBUG)\
 $(TARG_EX_C_2_ST_UX_DEBUG)\
 $(TARG_EX_CPP_1_DLL_T)\
 $(TARG_EX_CPP_2_DLL_T)\
 $(TARG_EX_CPP_1_MT_T)\
 $(TARG_EX_CPP_2_MT_T)\
 $(TARG_EX_CPP_1_ST_T)\
 $(TARG_EX_CPP_2_ST_T)\
 $(TARG_EX_CPP_1_ST_UX_DEBUG)\
 $(TARG_EX_CPP_2_ST_UX_DEBUG)\
 \

########################################
# Targets

TARGETS_LIB_CORE = \
 $(TARGET_LIB_DLL_T)\
 $(TARGET_LIB_UN_DLL_T)\
 $(LIB_CORE_DLL_T)\
 $(LIB_CORE_UN_DLL_T)\
 $(TARGET_LIB_MT_T)\
 $(TARGET_LIB_UN_MT_T)\
 $(LIB_CORE_MT_T)\
 $(LIB_CORE_UN_MT_T)\
 $(TARGET_LIB_ST_T)\
 $(TARGET_LIB_UN_ST_T)\
 $(LIB_CORE_ST_T)\
 $(LIB_CORE_UN_ST_T)\
 $(TARGET_LIB_ST_UX_DEBUG)\
 $(TARGET_LIB_UN_ST_UX_DEBUG)\
 $(LIB_CORE_ST_UX_DEBUG)\
 $(LIB_CORE_UN_ST_UX_DEBUG)\
 $(TARGET_PDB_ST_UX_DEBUG)\
 \

TARGETS_LIB_XCONTRACT = \
 $(LIB_XCONTRACT_DLL_T)\
 $(LIB_XCONTRACT_MT_T)\
 $(LIB_XCONTRACT_ST_T)\
 $(LIB_XCONTRACT_ST_UX_DEBUG)\
 \

TARGETS_LIB_XTESTS = \
 $(LIB_XTESTS_DLL_T)\
 $(LIB_XTESTS_MT_T)\
 $(LIB_XTESTS_ST_T)\
 $(LIB_XTESTS_ST_UX_DEBUG)\
 \

########################################
# OBJS_CORE

OBJ_CORE_b64__DLL_T = $(PROJ_OBJ_DIR)\b64.dll.test.obj

OBJS_CORE_DLL_T = \
 $(OBJ_CORE_b64__DLL_T)\
 \


OBJ_CORE_b64__MT_T = $(PROJ_OBJ_DIR)\b64.mt.test.obj

OBJS_CORE_MT_T = \
 $(OBJ_CORE_b64__MT_T)\
 \


OBJ_CORE_b64__ST_T = $(PROJ_OBJ_DIR)\b64.test.obj

OBJS_CORE_ST_T = \
 $(OBJ_CORE_b64__ST_T)\
 \


OBJ_CORE_b64__ST_UX_DEBUG = $(PROJ_OBJ_DIR)\b64.unix.debug.obj

OBJS_CORE_ST_UX_DEBUG = \
 $(OBJ_CORE_b64__ST_UX_DEBUG)\
 \


OBJS_CORE =\
 $(OBJS_CORE_DLL_T)\
 $(OBJS_CORE_MT_T)\
 $(OBJS_CORE_ST_T)\
 $(OBJS_CORE_ST_UX_DEBUG)\
 \

# xcontract

OBJ_xct_core_DLL_T           = $(PROJ_OBJ_DIR)\xcontract.b64.core.api.dll.test.obj

OBJ_XCONTRACT_DLL_T = \
 $(OBJ_xct_core_DLL_T)\
 \

OBJ_xct_core_MT_T           = $(PROJ_OBJ_DIR)\xcontract.b64.core.api.mt.test.obj

OBJ_XCONTRACT_MT_T = \
 $(OBJ_xct_core_MT_T)\
 \

OBJ_xct_core_ST_T           = $(PROJ_OBJ_DIR)\xcontract.b64.core.api.test.obj

OBJ_XCONTRACT_ST_T = \
 $(OBJ_xct_core_ST_T)\
 \

OBJ_xct_core_ST_UX_DEBUG           = $(PROJ_OBJ_DIR)\xcontract.b64.core.api.unix.debug.obj
PDB_xct_core_ST_UX_DEBUG           = $(PROJ_OBJ_DIR)\xcontract.b64.core.api.unix.debug.pdb

OBJ_XCONTRACT_ST_UX_DEBUG = \
 $(OBJ_xct_core_ST_UX_DEBUG)\
 \

PDB_CORE_ST_UX_DEBUG = \
 $(PDB_xct_core_ST_UX_DEBUG)\
 \


# xtests

OBJ_xts_core_DLL_T           = $(PROJ_OBJ_DIR)\xtests.b64.core.api.dll.test.obj

OBJ_XTESTS_DLL_T = \
 $(OBJ_xts_core_DLL_T)\
 \

OBJ_xts_core_MT_T           = $(PROJ_OBJ_DIR)\xtests.b64.core.api.mt.test.obj

OBJ_XTESTS_MT_T = \
 $(OBJ_xts_core_MT_T)\
 \

OBJ_xts_core_ST_T           = $(PROJ_OBJ_DIR)\xtests.b64.core.api.test.obj

OBJ_XTESTS_ST_T = \
 $(OBJ_xts_core_ST_T)\
 \

OBJ_xts_core_ST_UX_DEBUG           = $(PROJ_OBJ_DIR)\xtests.b64.core.api.unix.debug.obj
PDB_xts_core_ST_UX_DEBUG           = $(PROJ_OBJ_DIR)\xtests.b64.core.api.unix.debug.pdb

OBJ_XTESTS_ST_UX_DEBUG = \
 $(OBJ_xts_core_ST_UX_DEBUG)\
 \

PDB_CORE_ST_UX_DEBUG = \
 $(PDB_xts_core_ST_UX_DEBUG)\
 \


########################################
# OBJS_LIBS

OBJS_LIBS =\
 $(OBJS_CORE)\
 \

########################################
# OBJS_TESTS

OBJ_TEST_UNIT_C_1_DLL_T = $(PROJ_OBJ_DIR)\test.c.1.dll.test.obj
OBJ_TEST_UNIT_CPP_1_DLL_T = $(PROJ_OBJ_DIR)\test.cpp.1.dll.test.obj
OBJ_TEST_UNIT_C_1_MT_T = $(PROJ_OBJ_DIR)\test.c.1.mt.test.obj
OBJ_TEST_UNIT_CPP_1_MT_T = $(PROJ_OBJ_DIR)\test.cpp.1.mt.test.obj
OBJ_TEST_UNIT_C_1_ST_T = $(PROJ_OBJ_DIR)\test.c.1.test.obj
OBJ_TEST_UNIT_CPP_1_ST_T = $(PROJ_OBJ_DIR)\test.cpp.1.test.obj
OBJ_TEST_UNIT_C_1_ST_UX_DEBUG = $(PROJ_OBJ_DIR)\test.c.1.unix.debug.obj
OBJ_TEST_UNIT_CPP_1_ST_UX_DEBUG = $(PROJ_OBJ_DIR)\test.cpp.1.unix.debug.obj

OBJS_TESTS =\
 $(OBJ_TEST_UNIT_C_1_DLL_T)\
 $(OBJ_TEST_UNIT_C_1_MT_T)\
 $(OBJ_TEST_UNIT_C_1_ST_T)\
 $(OBJ_TEST_UNIT_C_1_ST_UX_DEBUG)\
 $(OBJ_TEST_UNIT_CPP_1_DLL_T)\
 $(OBJ_TEST_UNIT_CPP_1_MT_T)\
 $(OBJ_TEST_UNIT_CPP_1_ST_T)\
 $(OBJ_TEST_UNIT_CPP_1_ST_UX_DEBUG)\
 \

########################################
# OBJS_EXAMPLES

OBJ_EX_C_1_DLL_T = $(PROJ_OBJ_DIR)\example.c.1.dll.test.obj
OBJ_EX_C_2_DLL_T = $(PROJ_OBJ_DIR)\example.c.2.dll.test.obj
OBJ_EX_CPP_1_DLL_T = $(PROJ_OBJ_DIR)\example.cpp.1.dll.test.obj
OBJ_EX_CPP_2_DLL_T = $(PROJ_OBJ_DIR)\example.cpp.2.dll.test.obj
OBJ_EX_C_1_MT_T = $(PROJ_OBJ_DIR)\example.c.1.mt.test.obj
OBJ_EX_C_2_MT_T = $(PROJ_OBJ_DIR)\example.c.2.mt.test.obj
OBJ_EX_CPP_1_MT_T = $(PROJ_OBJ_DIR)\example.cpp.1.mt.test.obj
OBJ_EX_CPP_2_MT_T = $(PROJ_OBJ_DIR)\example.cpp.2.mt.test.obj
OBJ_EX_C_1_ST_T = $(PROJ_OBJ_DIR)\example.c.1.test.obj
OBJ_EX_C_2_ST_T = $(PROJ_OBJ_DIR)\example.c.2.test.obj
OBJ_EX_CPP_1_ST_T = $(PROJ_OBJ_DIR)\example.cpp.1.test.obj
OBJ_EX_CPP_2_ST_T = $(PROJ_OBJ_DIR)\example.cpp.2.test.obj
OBJ_EX_C_1_ST_UX_DEBUG = $(PROJ_OBJ_DIR)\example.c.1.unix.debug.obj
OBJ_EX_C_2_ST_UX_DEBUG = $(PROJ_OBJ_DIR)\example.c.2.unix.debug.obj
OBJ_EX_CPP_1_ST_UX_DEBUG = $(PROJ_OBJ_DIR)\example.cpp.1.unix.debug.obj
OBJ_EX_CPP_2_ST_UX_DEBUG = $(PROJ_OBJ_DIR)\example.cpp.2.unix.debug.obj

OBJS_EXAMPLES =\
 $(OBJ_EX_C_1_DLL_T)\
 $(OBJ_EX_C_2_DLL_T)\
 $(OBJ_EX_C_1_MT_T)\
 $(OBJ_EX_C_2_MT_T)\
 $(OBJ_EX_C_1_ST_T)\
 $(OBJ_EX_C_2_ST_T)\
 $(OBJ_EX_C_1_ST_UX_DEBUG)\
 $(OBJ_EX_C_2_ST_UX_DEBUG)\
 $(OBJ_EX_CPP_1_DLL_T)\
 $(OBJ_EX_CPP_2_DLL_T)\
 $(OBJ_EX_CPP_1_MT_T)\
 $(OBJ_EX_CPP_2_MT_T)\
 $(OBJ_EX_CPP_1_ST_T)\
 $(OBJ_EX_CPP_2_ST_T)\
 $(OBJ_EX_CPP_1_ST_UX_DEBUG)\
 $(OBJ_EX_CPP_2_ST_UX_DEBUG)\
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

LIB_CORE_b64_H = $(PROJ_INCLUDE_DIR)\b64\b64.h
LIB_CORE_implicit_link_H = $(PROJ_INCLUDE_DIR)\b64\implicit_link.h

LIB_CORE_b64_C = $(PROJ_SRC_DIR)\b64.c

########################################
# SRC_HEADERS

HEADERS_LIB	= \
 $(LIB_CORE_b64_H)\
 $(LIB_CORE_implicit_link_H)\
 \

HEADERS_MAPPING_CPP =\
 $(PROJ_INCLUDE_DIR)\b64\b64.hpp\
 \

########################################
# SRC_TESTS

SRC_TEST_C_1 = $(PROJ_TEST_UNIT_DIR)\test.unit.c.1\test.unit.c.1.c
SRC_TEST_CPP_1 = $(PROJ_TEST_UNIT_DIR)\test.unit.cpp.1\test.unit.cpp.1.cpp

########################################
# SRC_EXAMPLES

SRC_EXAMPLE_C_1 = $(PROJ_EXAMPLES_DIR)\c\example.c.1\example.c.1.c
SRC_EXAMPLE_C_2 = $(PROJ_EXAMPLES_DIR)\c\example.c.2\example.c.2.c
SRC_EXAMPLE_CPP_1 = $(PROJ_EXAMPLES_DIR)\cpp\example.cpp.1\example.cpp.1.cpp
SRC_EXAMPLE_CPP_2 = $(PROJ_EXAMPLES_DIR)\cpp\example.cpp.2\example.cpp.2.cpp

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
	$(TARG_TEST_UNIT_C_1_DLL_T) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_C_1_MT_T) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_C_1_ST_T) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_C_1_ST_UX_DEBUG) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_CPP_1_DLL_T) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_CPP_1_MT_T) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_CPP_1_ST_T) --verbosity=$(TEST_VERBOSITY)
	$(TARG_TEST_UNIT_CPP_1_ST_UX_DEBUG) --verbosity=$(TEST_VERBOSITY)

test:		build.test test_

test.all:	test

test-all:	test.all	# backwards compatibility


usage:
	@echo "make [usage | clean | compile | all | test | target] - default is 'all'"

clean.libs.core:
	@-$(RM) $(TARGET_LIB_DLL_T)
	@-$(RM) $(TARGET_LIB_UN_DLL_T)
	@-$(RM) $(LIB_CORE_DLL_T)
	@-$(RM) $(LIB_CORE_UN_DLL_T)
	@-$(RM) $(TARGET_LIB_MT_T)
	@-$(RM) $(TARGET_LIB_UN_MT_T)
	@-$(RM) $(LIB_CORE_MT_T)
	@-$(RM) $(LIB_CORE_UN_MT_T)
	@-$(RM) $(TARGET_LIB_ST_T)
	@-$(RM) $(TARGET_LIB_UN_ST_T)
	@-$(RM) $(LIB_CORE_ST_T)
	@-$(RM) $(LIB_CORE_UN_ST_T)
	@-$(RM) $(TARGET_LIB_ST_UX_DEBUG)
	@-$(RM) $(TARGET_LIB_UN_ST_UX_DEBUG)
	@-$(RM) $(LIB_CORE_ST_UX_DEBUG)
	@-$(RM) $(LIB_CORE_UN_ST_UX_DEBUG)
	@-$(RM) $(OBJ_CORE_b64__DLL_T)
	@-$(RM) $(OBJ_CORE_b64__MT_T)
	@-$(RM) $(OBJ_CORE_b64__ST_T)
	@-$(RM) $(OBJ_CORE_b64__ST_UX_DEBUG)

clean.libs.xcontract:	makefile

clean.libs.xtests:		makefile

clean.libs:				clean.libs.core clean.libs.xcontract clean.libs.xtests

clean.test:
	@-$(RM) $(TARG_TEST_UNIT_C_1_DLL_T)
	@-$(RM) $(TARG_TEST_UNIT_C_1_MT_T)
	@-$(RM) $(TARG_TEST_UNIT_C_1_ST_T)
	@-$(RM) $(TARG_TEST_UNIT_C_1_ST_UX_DEBUG)
	@-$(RM) $(TARG_TEST_UNIT_CPP_1_DLL_T)
	@-$(RM) $(TARG_TEST_UNIT_CPP_1_MT_T)
	@-$(RM) $(TARG_TEST_UNIT_CPP_1_ST_T)
	@-$(RM) $(TARG_TEST_UNIT_CPP_1_ST_UX_DEBUG)
	@-$(RM) $(OBJ_TEST_UNIT_C_1_DLL_T)
	@-$(RM) $(OBJ_TEST_UNIT_C_1_MT_T)
	@-$(RM) $(OBJ_TEST_UNIT_C_1_ST_T)
	@-$(RM) $(OBJ_TEST_UNIT_C_1_ST_UX_DEBUG)
	@-$(RM) $(OBJ_TEST_UNIT_CPP_1_DLL_T)
	@-$(RM) $(OBJ_TEST_UNIT_CPP_1_MT_T)
	@-$(RM) $(OBJ_TEST_UNIT_CPP_1_ST_T)
	@-$(RM) $(OBJ_TEST_UNIT_CPP_1_ST_UX_DEBUG)

clean.examples:
	@-$(RM) $(TARG_EX_C_1_DLL_T)
	@-$(RM) $(TARG_EX_C_2_DLL_T)
	@-$(RM) $(TARG_EX_C_1_MT_T)
	@-$(RM) $(TARG_EX_C_2_MT_T)
	@-$(RM) $(TARG_EX_C_1_ST_T)
	@-$(RM) $(TARG_EX_C_2_ST_T)
	@-$(RM) $(TARG_EX_C_1_ST_UX_DEBUG)
	@-$(RM) $(TARG_EX_C_2_ST_UX_DEBUG)
	@-$(RM) $(TARG_EX_CPP_1_DLL_T)
	@-$(RM) $(TARG_EX_CPP_2_DLL_T)
	@-$(RM) $(TARG_EX_CPP_1_MT_T)
	@-$(RM) $(TARG_EX_CPP_2_MT_T)
	@-$(RM) $(TARG_EX_CPP_1_ST_T)
	@-$(RM) $(TARG_EX_CPP_2_ST_T)
	@-$(RM) $(TARG_EX_CPP_1_ST_UX_DEBUG)
	@-$(RM) $(TARG_EX_CPP_2_ST_UX_DEBUG)
	@-$(RM) $(OBJ_EX_C_1_DLL_T)
	@-$(RM) $(OBJ_EX_C_2_DLL_T)
	@-$(RM) $(OBJ_EX_C_1_MT_T)
	@-$(RM) $(OBJ_EX_C_2_MT_T)
	@-$(RM) $(OBJ_EX_C_1_ST_T)
	@-$(RM) $(OBJ_EX_C_2_ST_T)
	@-$(RM) $(OBJ_EX_C_1_ST_UX_DEBUG)
	@-$(RM) $(OBJ_EX_C_2_ST_UX_DEBUG)
	@-$(RM) $(OBJ_EX_CPP_1_DLL_T)
	@-$(RM) $(OBJ_EX_CPP_2_DLL_T)
	@-$(RM) $(OBJ_EX_CPP_1_MT_T)
	@-$(RM) $(OBJ_EX_CPP_2_MT_T)
	@-$(RM) $(OBJ_EX_CPP_1_ST_T)
	@-$(RM) $(OBJ_EX_CPP_2_ST_T)
	@-$(RM) $(OBJ_EX_CPP_1_ST_UX_DEBUG)
	@-$(RM) $(OBJ_EX_CPP_2_ST_UX_DEBUG)

clean: clean.libs clean.test clean.examples

################################################################################
# Internal Targets (invokable by make)

########################################
# Main library - 'b64.X.core.a'


$(TARGET_LIB_DLL_T):	$(LIB_CORE_DLL_T)
	@$(CP) $(LIB_CORE_DLL_T) $@

$(TARGET_LIB_UN_DLL_T):	$(LIB_CORE_DLL_T)
	@$(CP) $(LIB_CORE_DLL_T) $@

$(LIB_CORE_UN_DLL_T): $(LIB_CORE_DLL_T)
	@$(CP) $(LIB_CORE_DLL_T) $@

$(LIB_CORE_DLL_T): $(OBJS_CORE_DLL_T) makefile
	$(AR) $(AR_ARGS) -out:$@ $(OBJS_CORE_DLL_T)

$(OBJ_CORE_b64__DLL_T): $(PROJ_SRC_DIR)\b64.c $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_DLL_T) -Fo$@ $(PROJ_SRC_DIR)\b64.c


$(TARGET_LIB_MT_T):	$(LIB_CORE_MT_T)
	@$(CP) $(LIB_CORE_MT_T) $@

$(TARGET_LIB_UN_MT_T):	$(LIB_CORE_MT_T)
	@$(CP) $(LIB_CORE_MT_T) $@

$(LIB_CORE_UN_MT_T): $(LIB_CORE_MT_T)
	@$(CP) $(LIB_CORE_MT_T) $@

$(LIB_CORE_MT_T): $(OBJS_CORE_MT_T) makefile
	$(AR) $(AR_ARGS) -out:$@ $(OBJS_CORE_MT_T)

$(OBJ_CORE_b64__MT_T): $(PROJ_SRC_DIR)\b64.c $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_MT_T) -Fo$@ $(PROJ_SRC_DIR)\b64.c


$(TARGET_LIB_ST_T):	$(LIB_CORE_ST_T)
	@$(CP) $(LIB_CORE_ST_T) $@

$(TARGET_LIB_UN_ST_T):	$(LIB_CORE_ST_T)
	@$(CP) $(LIB_CORE_ST_T) $@

$(LIB_CORE_UN_ST_T): $(LIB_CORE_ST_T)
	@$(CP) $(LIB_CORE_ST_T) $@

$(LIB_CORE_ST_T): $(OBJS_CORE_ST_T) makefile
	$(AR) $(AR_ARGS) -out:$@ $(OBJS_CORE_ST_T)

$(OBJ_CORE_b64__ST_T): $(PROJ_SRC_DIR)\b64.c $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_ST_T) -Fo$@ $(PROJ_SRC_DIR)\b64.c

$(TARGET_PDB_ST_UX_DEBUG):	$(TARGET_LIB_ST_UX_DEBUG)

$(TARGET_LIB_ST_UX_DEBUG):	$(LIB_CORE_ST_UX_DEBUG)
	@$(CP) $(LIB_CORE_ST_UX_DEBUG) $@

$(TARGET_LIB_UN_ST_UX_DEBUG):	$(LIB_CORE_ST_UX_DEBUG)
	@$(CP) $(LIB_CORE_ST_UX_DEBUG) $@

$(LIB_CORE_UN_ST_UX_DEBUG): $(LIB_CORE_ST_UX_DEBUG)
	@$(CP) $(LIB_CORE_ST_UX_DEBUG) $@

$(LIB_CORE_ST_UX_DEBUG): $(OBJS_CORE_ST_UX_DEBUG) makefile
	$(AR) $(AR_ARGS) -out:$@ $(OBJS_CORE_ST_UX_DEBUG)

$(OBJ_CORE_b64__ST_UX_DEBUG): $(PROJ_SRC_DIR)\b64.c $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_ST_UX_DEBUG) -Fo$@ $(PROJ_SRC_DIR)\b64.c


########################################
# xContract library - 'xcontract.X.b64.core.a'

$(LIB_XCONTRACT_DLL_T): $(OBJ_XCONTRACT_DLL_T) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)\NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) -out:$@ $(OBJ_XCONTRACT_DLL_T)

$(OBJ_xct_core_DLL_T): $(XCONTRACT_SRC_DIR)\xcontract.core.c $(XCONTRACT_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CC) $(CC_ARGS_DLL_T) $(CC_XCONTRACT_INCLUDES_DLL_T) -Fo$@ $(XCONTRACT_SRC_DIR)\xcontract.core.c

$(LIB_XCONTRACT_MT_T): $(OBJ_XCONTRACT_MT_T) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)\NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) -out:$@ $(OBJ_XCONTRACT_MT_T)

$(OBJ_xct_core_MT_T): $(XCONTRACT_SRC_DIR)\xcontract.core.c $(XCONTRACT_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CC) $(CC_ARGS_MT_T) $(CC_XCONTRACT_INCLUDES_MT_T) -Fo$@ $(XCONTRACT_SRC_DIR)\xcontract.core.c

$(LIB_XCONTRACT_ST_T): $(OBJ_XCONTRACT_ST_T) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)\NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) -out:$@ $(OBJ_XCONTRACT_ST_T)

$(OBJ_xct_core_ST_T): $(XCONTRACT_SRC_DIR)\xcontract.core.c $(XCONTRACT_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CC) $(CC_ARGS_ST_T) $(CC_XCONTRACT_INCLUDES_ST_T) -Fo$@ $(XCONTRACT_SRC_DIR)\xcontract.core.c

$(LIB_XCONTRACT_ST_UX_DEBUG): $(OBJ_XCONTRACT_ST_UX_DEBUG) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)\NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) -out:$@ $(OBJ_XCONTRACT_ST_UX_DEBUG)

$(OBJ_xct_core_ST_UX_DEBUG): $(XCONTRACT_SRC_DIR)\xcontract.core.c $(XCONTRACT_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CC) $(CC_ARGS_ST_UX_DEBUG) $(CC_XCONTRACT_INCLUDES_ST_UX_DEBUG) -Fo$@ $(XCONTRACT_SRC_DIR)\xcontract.core.c

########################################
# xTests library - 'xtests.X.b64.core.a'

$(LIB_XTESTS_DLL_T): $(OBJ_XTESTS_DLL_T) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)\NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) -out:$@ $(OBJ_XTESTS_DLL_T)

$(OBJ_xts_core_DLL_T): $(XTESTS_SRC_DIR)\xtests.core.cpp $(XTESTS_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CPPC) $(CPPC_ARGS_DLL_T) $(CPPC_XTESTS_INCLUDES_DLL_T) -Fo$@ $(XTESTS_SRC_DIR)\xtests.core.cpp

$(LIB_XTESTS_MT_T): $(OBJ_XTESTS_MT_T) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)\NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) -out:$@ $(OBJ_XTESTS_MT_T)

$(OBJ_xts_core_MT_T): $(XTESTS_SRC_DIR)\xtests.core.cpp $(XTESTS_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CPPC) $(CPPC_ARGS_MT_T) $(CPPC_XTESTS_INCLUDES_MT_T) -Fo$@ $(XTESTS_SRC_DIR)\xtests.core.cpp

$(LIB_XTESTS_ST_T): $(OBJ_XTESTS_ST_T) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)\NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) -out:$@ $(OBJ_XTESTS_ST_T)

$(OBJ_xts_core_ST_T): $(XTESTS_SRC_DIR)\xtests.core.cpp $(XTESTS_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CPPC) $(CPPC_ARGS_ST_T) $(CPPC_XTESTS_INCLUDES_ST_T) -Fo$@ $(XTESTS_SRC_DIR)\xtests.core.cpp

$(LIB_XTESTS_ST_UX_DEBUG): $(OBJ_XTESTS_ST_UX_DEBUG) $(PROJ_MAKEFILE_NAME)
	@-if not exist $(PROJ_LIB_DIR)\NUL mkdir $(PROJ_LIB_DIR)
	$(AR) $(AR_ARGS) -out:$@ $(OBJ_XTESTS_ST_UX_DEBUG)

$(OBJ_xts_core_ST_UX_DEBUG): $(XTESTS_SRC_DIR)\xtests.core.cpp $(XTESTS_HEADERS) $(PROJ_MAKEFILE_NAME)
	$(CPPC) $(CPPC_ARGS_ST_UX_DEBUG) $(CPPC_XTESTS_INCLUDES_ST_UX_DEBUG) -Fo$@ $(XTESTS_SRC_DIR)\xtests.core.cpp

########################################
# test

# C test client

$(TARG_TEST_UNIT_C_1_DLL_T):	\
	$(OBJ_TEST_UNIT_C_1_DLL_T)\
	$(TARGET_LIB_DLL_T)\
	$(LIB_XCONTRACT_DLL_T)\
	$(LIB_XTESTS_DLL_T)\
	makefile
	$(LD) $(LD_ARGS_DLL_T) -out:$@ $(OBJ_TEST_UNIT_C_1_DLL_T) $(TARGET_LIB_DLL_T) $(LIB_XTESTS_DLL_T) $(LIB_XCONTRACT_DLL_T) $(3PTY_LIBS__DLL_T)

$(OBJ_TEST_UNIT_C_1_DLL_T):	$(SRC_TEST_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_DLL_T) -Fo$@ $(SRC_TEST_C_1)

$(TARG_TEST_UNIT_C_1_MT_T):	\
	$(OBJ_TEST_UNIT_C_1_MT_T)\
	$(TARGET_LIB_MT_T)\
	$(LIB_XCONTRACT_MT_T)\
	$(LIB_XTESTS_MT_T)\
	makefile
	$(LD) $(LD_ARGS_MT_T) -out:$@ $(OBJ_TEST_UNIT_C_1_MT_T) $(TARGET_LIB_MT_T) $(LIB_XTESTS_MT_T) $(LIB_XCONTRACT_MT_T) $(3PTY_LIBS__MT_T)

$(OBJ_TEST_UNIT_C_1_MT_T):	$(SRC_TEST_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_MT_T) -Fo$@ $(SRC_TEST_C_1)

$(TARG_TEST_UNIT_C_1_ST_T):	\
	$(OBJ_TEST_UNIT_C_1_ST_T)\
	$(TARGET_LIB_ST_T)\
	$(LIB_XCONTRACT_ST_T)\
	$(LIB_XTESTS_ST_T)\
	makefile
	$(LD) $(LD_ARGS_ST_T) -out:$@ $(OBJ_TEST_UNIT_C_1_ST_T) $(TARGET_LIB_ST_T) $(LIB_XTESTS_ST_T) $(LIB_XCONTRACT_ST_T) $(3PTY_LIBS__ST_T)

$(OBJ_TEST_UNIT_C_1_ST_T):	$(SRC_TEST_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_ST_T) -Fo$@ $(SRC_TEST_C_1)

$(TARG_TEST_UNIT_C_1_ST_UX_DEBUG):	\
	$(OBJ_TEST_UNIT_C_1_ST_UX_DEBUG)\
	$(TARGET_LIB_ST_UX_DEBUG)\
	$(LIB_XCONTRACT_ST_UX_DEBUG)\
	$(LIB_XTESTS_ST_UX_DEBUG)\
	makefile
	$(LD) $(LD_ARGS_ST_UX_DEBUG) -out:$@ $(OBJ_TEST_UNIT_C_1_ST_UX_DEBUG) $(TARGET_LIB_ST_UX_DEBUG) $(LIB_XTESTS_ST_UX_DEBUG) $(LIB_XCONTRACT_ST_UX_DEBUG) $(3PTY_LIBS__ST_UX_DEBUG)

$(OBJ_TEST_UNIT_C_1_ST_UX_DEBUG):	$(SRC_TEST_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_ST_UX_DEBUG) -Fo$@ $(SRC_TEST_C_1)

# C++ test client

$(TARG_TEST_UNIT_CPP_1_DLL_T):	\
	$(OBJ_TEST_UNIT_CPP_1_DLL_T)\
	$(TARGET_LIB_DLL_T)\
	$(LIB_XCONTRACT_DLL_T)\
	$(LIB_XTESTS_DLL_T)\
	makefile
	$(LD) $(LD_ARGS_DLL_T) -out:$@ $(OBJ_TEST_UNIT_CPP_1_DLL_T) $(TARGET_LIB_DLL_T) $(LIB_XTESTS_DLL_T) $(LIB_XCONTRACT_DLL_T) $(3PTY_LIBS__DLL_T)

$(OBJ_TEST_UNIT_CPP_1_DLL_T):	$(SRC_TEST_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_DLL_T) -Fo$@ $(SRC_TEST_CPP_1)

$(TARG_TEST_UNIT_CPP_1_MT_T):	\
	$(OBJ_TEST_UNIT_CPP_1_MT_T)\
	$(TARGET_LIB_MT_T)\
	$(LIB_XCONTRACT_MT_T)\
	$(LIB_XTESTS_MT_T)\
	makefile
	$(LD) $(LD_ARGS_MT_T) -out:$@ $(OBJ_TEST_UNIT_CPP_1_MT_T) $(TARGET_LIB_MT_T) $(LIB_XTESTS_MT_T) $(LIB_XCONTRACT_MT_T) $(3PTY_LIBS__MT_T)

$(OBJ_TEST_UNIT_CPP_1_MT_T):	$(SRC_TEST_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_MT_T) -Fo$@ $(SRC_TEST_CPP_1)

$(TARG_TEST_UNIT_CPP_1_ST_T):	\
	$(OBJ_TEST_UNIT_CPP_1_ST_T)\
	$(TARGET_LIB_ST_T)\
	$(LIB_XCONTRACT_ST_T)\
	$(LIB_XTESTS_ST_T)\
	makefile
	$(LD) $(LD_ARGS_ST_T) -out:$@ $(OBJ_TEST_UNIT_CPP_1_ST_T) $(TARGET_LIB_ST_T) $(LIB_XTESTS_ST_T) $(LIB_XCONTRACT_ST_T) $(3PTY_LIBS__ST_T)

$(OBJ_TEST_UNIT_CPP_1_ST_T):	$(SRC_TEST_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_ST_T) -Fo$@ $(SRC_TEST_CPP_1)

$(TARG_TEST_UNIT_CPP_1_ST_UX_DEBUG):	\
	$(OBJ_TEST_UNIT_CPP_1_ST_UX_DEBUG)\
	$(TARGET_LIB_ST_UX_DEBUG)\
	$(LIB_XCONTRACT_ST_UX_DEBUG)\
	$(LIB_XTESTS_ST_UX_DEBUG)\
	makefile
	$(LD) $(LD_ARGS_ST_UX_DEBUG) -out:$@ $(OBJ_TEST_UNIT_CPP_1_ST_UX_DEBUG) $(TARGET_LIB_ST_UX_DEBUG) $(LIB_XTESTS_ST_UX_DEBUG) $(LIB_XCONTRACT_ST_UX_DEBUG) $(3PTY_LIBS__ST_UX_DEBUG)

$(OBJ_TEST_UNIT_CPP_1_ST_UX_DEBUG):	$(SRC_TEST_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_ST_UX_DEBUG) -Fo$@ $(SRC_TEST_CPP_1)

########################################
# examples

# C example client

$(TARG_EX_C_1_DLL_T):	\
	$(OBJ_EX_C_1_DLL_T)\
	$(TARGET_LIB_DLL_T)\
	$(LIB_XCONTRACT_DLL_T)\
	makefile
	$(LD) $(LD_ARGS_DLL_T) -out:$@ $(OBJ_EX_C_1_DLL_T) $(TARGET_LIB_DLL_T) $(LIB_XCONTRACT_DLL_T) $(3PTY_LIBS__DLL_T)

$(OBJ_EX_C_1_DLL_T):	$(SRC_EXAMPLE_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_DLL_T) -Fo$@ $(SRC_EXAMPLE_C_1)

$(TARG_EX_C_2_DLL_T):	\
	$(OBJ_EX_C_2_DLL_T)\
	$(TARGET_LIB_DLL_T)\
	$(LIB_XCONTRACT_DLL_T)\
	makefile
	$(LD) $(LD_ARGS_DLL_T) -out:$@ $(OBJ_EX_C_2_DLL_T) $(TARGET_LIB_DLL_T) $(LIB_XCONTRACT_DLL_T) $(3PTY_LIBS__DLL_T)

$(OBJ_EX_C_2_DLL_T):	$(SRC_EXAMPLE_C_2) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_DLL_T) -Fo$@ $(SRC_EXAMPLE_C_2)

$(TARG_EX_C_1_MT_T):	\
	$(OBJ_EX_C_1_MT_T)\
	$(TARGET_LIB_MT_T)\
	$(LIB_XCONTRACT_MT_T)\
	makefile
	$(LD) $(LD_ARGS_MT_T) -out:$@ $(OBJ_EX_C_1_MT_T) $(TARGET_LIB_MT_T) $(LIB_XCONTRACT_MT_T) $(3PTY_LIBS__MT_T)

$(OBJ_EX_C_1_MT_T):	$(SRC_EXAMPLE_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_MT_T) -Fo$@ $(SRC_EXAMPLE_C_1)

$(TARG_EX_C_2_MT_T):	\
	$(OBJ_EX_C_2_MT_T)\
	$(TARGET_LIB_MT_T)\
	$(LIB_XCONTRACT_MT_T)\
	makefile
	$(LD) $(LD_ARGS_MT_T) -out:$@ $(OBJ_EX_C_2_MT_T) $(TARGET_LIB_MT_T) $(LIB_XCONTRACT_MT_T) $(3PTY_LIBS__MT_T)

$(OBJ_EX_C_2_MT_T):	$(SRC_EXAMPLE_C_2) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_MT_T) -Fo$@ $(SRC_EXAMPLE_C_2)

$(TARG_EX_C_1_ST_T):	\
	$(OBJ_EX_C_1_ST_T)\
	$(TARGET_LIB_ST_T)\
	$(LIB_XCONTRACT_ST_T)\
	makefile
	$(LD) $(LD_ARGS_ST_T) -out:$@ $(OBJ_EX_C_1_ST_T) $(TARGET_LIB_ST_T) $(LIB_XCONTRACT_ST_T) $(3PTY_LIBS__ST_T)

$(OBJ_EX_C_1_ST_T):	$(SRC_EXAMPLE_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_ST_T) -Fo$@ $(SRC_EXAMPLE_C_1)

$(TARG_EX_C_2_ST_T):	\
	$(OBJ_EX_C_2_ST_T)\
	$(TARGET_LIB_ST_T)\
	$(LIB_XCONTRACT_ST_T)\
	makefile
	$(LD) $(LD_ARGS_ST_T) -out:$@ $(OBJ_EX_C_2_ST_T) $(TARGET_LIB_ST_T) $(LIB_XCONTRACT_ST_T) $(3PTY_LIBS__ST_T)

$(OBJ_EX_C_2_ST_T):	$(SRC_EXAMPLE_C_2) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_ST_T) -Fo$@ $(SRC_EXAMPLE_C_2)

$(TARG_EX_C_1_ST_UX_DEBUG):	\
	$(OBJ_EX_C_1_ST_UX_DEBUG)\
	$(TARGET_LIB_ST_UX_DEBUG)\
	$(LIB_XCONTRACT_ST_UX_DEBUG)\
	makefile
	$(LD) $(LD_ARGS_ST_UX_DEBUG) -out:$@ $(OBJ_EX_C_1_ST_UX_DEBUG) $(TARGET_LIB_ST_UX_DEBUG) $(LIB_XCONTRACT_ST_UX_DEBUG) $(3PTY_LIBS__ST_UX_DEBUG)

$(OBJ_EX_C_1_ST_UX_DEBUG):	$(SRC_EXAMPLE_C_1) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_ST_UX_DEBUG) -Fo$@ $(SRC_EXAMPLE_C_1)

$(TARG_EX_C_2_ST_UX_DEBUG):	\
	$(OBJ_EX_C_2_ST_UX_DEBUG)\
	$(TARGET_LIB_ST_UX_DEBUG)\
	$(LIB_XCONTRACT_ST_UX_DEBUG)\
	makefile
	$(LD) $(LD_ARGS_ST_UX_DEBUG) -out:$@ $(OBJ_EX_C_2_ST_UX_DEBUG) $(TARGET_LIB_ST_UX_DEBUG) $(LIB_XCONTRACT_ST_UX_DEBUG) $(3PTY_LIBS__ST_UX_DEBUG)

$(OBJ_EX_C_2_ST_UX_DEBUG):	$(SRC_EXAMPLE_C_2) $(HEADERS_LIB) makefile
	$(CC) $(CC_ARGS_ST_UX_DEBUG) -Fo$@ $(SRC_EXAMPLE_C_2)

# C++ example client

$(TARG_EX_CPP_1_DLL_T):	\
	$(OBJ_EX_CPP_1_DLL_T)\
	$(TARGET_LIB_DLL_T)\
	$(LIB_XCONTRACT_DLL_T)\
	makefile
	$(LD) $(LD_ARGS_DLL_T) -out:$@ $(OBJ_EX_CPP_1_DLL_T) $(TARGET_LIB_DLL_T) $(LIB_XCONTRACT_DLL_T) $(3PTY_LIBS__DLL_T)

$(OBJ_EX_CPP_1_DLL_T):	$(SRC_EXAMPLE_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_DLL_T) -Fo$@ $(SRC_EXAMPLE_CPP_1)

$(TARG_EX_CPP_2_DLL_T):	\
	$(OBJ_EX_CPP_2_DLL_T)\
	$(TARGET_LIB_DLL_T)\
	$(LIB_XCONTRACT_DLL_T)\
	makefile
	$(LD) $(LD_ARGS_DLL_T) -out:$@ $(OBJ_EX_CPP_2_DLL_T) $(TARGET_LIB_DLL_T) $(LIB_XCONTRACT_DLL_T) $(3PTY_LIBS__DLL_T)

$(OBJ_EX_CPP_2_DLL_T):	$(SRC_EXAMPLE_CPP_2) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_DLL_T) -Fo$@ $(SRC_EXAMPLE_CPP_2)

$(TARG_EX_CPP_1_MT_T):	\
	$(OBJ_EX_CPP_1_MT_T)\
	$(TARGET_LIB_MT_T)\
	$(LIB_XCONTRACT_MT_T)\
	makefile
	$(LD) $(LD_ARGS_MT_T) -out:$@ $(OBJ_EX_CPP_1_MT_T) $(TARGET_LIB_MT_T) $(LIB_XCONTRACT_MT_T) $(3PTY_LIBS__MT_T)

$(OBJ_EX_CPP_1_MT_T):	$(SRC_EXAMPLE_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_MT_T) -Fo$@ $(SRC_EXAMPLE_CPP_1)

$(TARG_EX_CPP_2_MT_T):	\
	$(OBJ_EX_CPP_2_MT_T)\
	$(TARGET_LIB_MT_T)\
	$(LIB_XCONTRACT_MT_T)\
	makefile
	$(LD) $(LD_ARGS_MT_T) -out:$@ $(OBJ_EX_CPP_2_MT_T) $(TARGET_LIB_MT_T) $(LIB_XCONTRACT_MT_T) $(3PTY_LIBS__MT_T)

$(OBJ_EX_CPP_2_MT_T):	$(SRC_EXAMPLE_CPP_2) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_MT_T) -Fo$@ $(SRC_EXAMPLE_CPP_2)

$(TARG_EX_CPP_1_ST_T):	\
	$(OBJ_EX_CPP_1_ST_T)\
	$(TARGET_LIB_ST_T)\
	$(LIB_XCONTRACT_ST_T)\
	makefile
	$(LD) $(LD_ARGS_ST_T) -out:$@ $(OBJ_EX_CPP_1_ST_T) $(TARGET_LIB_ST_T) $(LIB_XCONTRACT_ST_T) $(3PTY_LIBS__ST_T)

$(OBJ_EX_CPP_1_ST_T):	$(SRC_EXAMPLE_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_ST_T) -Fo$@ $(SRC_EXAMPLE_CPP_1)

$(TARG_EX_CPP_2_ST_T):	\
	$(OBJ_EX_CPP_2_ST_T)\
	$(TARGET_LIB_ST_T)\
	$(LIB_XCONTRACT_ST_T)\
	makefile
	$(LD) $(LD_ARGS_ST_T) -out:$@ $(OBJ_EX_CPP_2_ST_T) $(TARGET_LIB_ST_T) $(LIB_XCONTRACT_ST_T) $(3PTY_LIBS__ST_T)

$(OBJ_EX_CPP_2_ST_T):	$(SRC_EXAMPLE_CPP_2) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_ST_T) -Fo$@ $(SRC_EXAMPLE_CPP_2)

$(TARG_EX_CPP_1_ST_UX_DEBUG):	\
	$(OBJ_EX_CPP_1_ST_UX_DEBUG)\
	$(TARGET_LIB_ST_UX_DEBUG)\
	$(LIB_XCONTRACT_ST_UX_DEBUG)\
	makefile
	$(LD) $(LD_ARGS_ST_UX_DEBUG) -out:$@ $(OBJ_EX_CPP_1_ST_UX_DEBUG) $(TARGET_LIB_ST_UX_DEBUG) $(LIB_XCONTRACT_ST_UX_DEBUG) $(3PTY_LIBS__ST_UX_DEBUG)

$(OBJ_EX_CPP_1_ST_UX_DEBUG):	$(SRC_EXAMPLE_CPP_1) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_ST_UX_DEBUG) -Fo$@ $(SRC_EXAMPLE_CPP_1)

$(TARG_EX_CPP_2_ST_UX_DEBUG):	\
	$(OBJ_EX_CPP_2_ST_UX_DEBUG)\
	$(TARGET_LIB_ST_UX_DEBUG)\
	$(LIB_XCONTRACT_ST_UX_DEBUG)\
	makefile
	$(LD) $(LD_ARGS_ST_UX_DEBUG) -out:$@ $(OBJ_EX_CPP_2_ST_UX_DEBUG) $(TARGET_LIB_ST_UX_DEBUG) $(LIB_XCONTRACT_ST_UX_DEBUG) $(3PTY_LIBS__ST_UX_DEBUG)

$(OBJ_EX_CPP_2_ST_UX_DEBUG):	$(SRC_EXAMPLE_CPP_2) $(HEADERS_LIB) $(HEADERS_MAPPING_CPP) makefile
	$(CPPC) $(CPPC_ARGS_ST_UX_DEBUG) -Fo$@ $(SRC_EXAMPLE_CPP_2)

################################################################################
