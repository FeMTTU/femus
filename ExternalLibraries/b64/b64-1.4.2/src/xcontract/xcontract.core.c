/* /////////////////////////////////////////////////////////////////////////
 * File:        xcontract.core.c (formerly part of Synesis' internal test codebase)
 *
 * Purpose:     Main implementation file for xContract, a contract enforcement
 *              library.
 *
 * Created:     2nd January 2001
 * Updated:     20th January 2010
 *
 * Home:        http://stlsoft.org/
 *
 * Copyright (c) 2001-2010, Matthew Wilson and Synesis Software
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * - Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 * - Neither the name(s) of Matthew Wilson and Synesis Software nor the
 *   names of any contributors may be used to endorse or promote products
 *   derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ////////////////////////////////////////////////////////////////////// */


/* /////////////////////////////////////////////////////////////////////////
 * Includes - 1
 */

#include <xcontract/xcontract.h>
#include <xcontract/internal/safestr.h>

/* STLSoft Header Files */
#include <stlsoft/stlsoft.h>

/* /////////////////////////////////////////////////////////////////////////
 * Compiler compatibility / features - 1
 */

#if _STLSOFT_VER < 0x010959ff
# error xContract requires version 1.9.89 (or later) of STLSoft; download from www.stlsoft.org
#endif /* _STLSOFT_VER */

/* /////////////////////////////////////////////////////////////////////////
 * Includes - 2
 */

/* STLSoft Header Files */
#include <platformstl/platformstl.h>
#include <stlsoft/shims/access/string/std/c_string.h>

#ifdef PLATFORMSTL_OS_IS_UNIX
# include <unistd.h>
#endif /* OS */

/* /////////////////////////////////////////////////////////////////////////
 * Platform/feature discrimination
 */

/* syslog (UNIX) */
#if !defined(XCONTRACT_NO_USE_SYSLOG) && \
    (   defined(PLATFORMSTL_OS_IS_UNIX) || \
        defined(XCONTRACT_USE_SYNESIS_SYSLOG))
# ifndef XCONTRACT_USE_SYSLOG
#  define XCONTRACT_USE_SYSLOG
# endif /* !XCONTRACT_USE_SYSLOG */
#endif /* !PLATFORMSTL_OS_IS_???? */

/* OutputDebugString (Windows) */
#if !defined(XCONTRACT_NO_USE_WIN_OUTPUTDEBUGSTRING) && \
    (   defined(_WIN32) || \
        defined(_WIN64) || \
        defined(PLATFORMSTL_OS_IS_WINDOWS))
# ifndef XCONTRACT_USE_WIN_OUTPUTDEBUGSTRING
#  define XCONTRACT_USE_WIN_OUTPUTDEBUGSTRING
# endif /* !XCONTRACT_USE_WIN_OUTPUTDEBUGSTRING */
#endif /* !PLATFORMSTL_OS_IS_???? */

/* Event Log (Windows) */
#if !defined(XCONTRACT_NO_USE_WIN_EVENTLOG) && \
    (   defined(_WIN32) || \
        defined(_WIN64) || \
        defined(PLATFORMSTL_OS_IS_WINDOWS))
# ifndef XCONTRACT_USE_WIN_EVENTLOG
#  define XCONTRACT_USE_WIN_EVENTLOG
# endif /* !XCONTRACT_USE_WIN_EVENTLOG */

  /* Event Identifier - defaults to 1004, which is the identifier of the
   * default message exported by Pantheios.COM
   */
# ifndef XCONTRACT_WIN_EVENTLOG_EVENTID
#  define XCONTRACT_WIN_EVENTLOG_EVENTID            (1004)
# endif /* !XCONTRACT_WIN_EVENTLOG_EVENTID */

  /* Category - defaults to 14, which is the identifier of the
   * contract category exported by Pantheios.COM
   */
# ifndef XCONTRACT_WIN_EVENTLOG_CATEGORY
#  define XCONTRACT_WIN_EVENTLOG_CATEGORY           (14)
# endif /* !XCONTRACT_WIN_EVENTLOG_CATEGORY */

#endif /* !PLATFORMSTL_OS_IS_???? */

/* /////////////////////////////////////////////////////////////////////////
 * Includes - 3
 */

/* Standard C Header Files */
#include <stdio.h>
#include <stdlib.h>

/* Operating system-/feature-specific Header Files */

#ifdef XCONTRACT_USE_SYSLOG
# include <syslog.h>

 /* Determine the Syslog flags here. Need:
  *
  * - LOG_EMERG - the level is emergency, indicating design violation
  * - LOG_USER - the facility appropriate to application programs
  * - WINSYSLOG_F_EVENTLOG - use the Windows event log, if available (i.e.
  *    if using the Synesis implementation of Syslog)
  */
# ifdef WINSYSLOG_F_EVENTLOG
#  define XCONTRACT_SYSLOG_FLAGS    (LOG_MAKEPRI(LOG_USER, LOG_EMERG) | WINSYSLOG_F_EVENTLOG)
# else /* ? WINSYSLOG_F_EVENTLOG */
#  define XCONTRACT_SYSLOG_FLAGS    LOG_MAKEPRI(LOG_USER, LOG_EMERG)
# endif /* WINSYSLOG_F_EVENTLOG */
#endif /* XCONTRACT_USE_SYSLOG */

#ifdef PLATFORMSTL_OS_IS_UNIX
# include <unistd.h>
#endif /* PLATFORMSTL_OS_IS_UNIX */

#if defined(XCONTRACT_USE_WIN_OUTPUTDEBUGSTRING) || \
    defined(XCONTRACT_USE_WIN_EVENTLOG)
# include <windows.h>
#endif /* XCONTRACT_USE_WIN_OUTPUTDEBUGSTRING || XCONTRACT_USE_WIN_EVENTLOG */

/* /////////////////////////////////////////////////////////////////////////
 * Compiler compatibility / features - 2
 */

#if defined(STLSOFT_COMPILER_IS_MSVC) || \
    (   defined(STLSOFT_COMPILER_IS_INTEL) && \
        defined(PLATFORMSTL_OS_IS_WINDOWS))
# define XCONTRACT_COMPILER_SUPPORTS_SEH
# ifndef EXCEPTION_EXECUTE_HANDLER
#  define EXCEPTION_EXECUTE_HANDLER     (1)
# endif /* !EXCEPTION_EXECUTE_HANDLER */
#endif /* compiler / OS */

#if !defined(XCONTRACT_USING_SAFE_STR_FUNCTIONS)
# if defined(STLSOFT_COMPILER_IS_MSVC) || \
     (   defined(STLSOFT_COMPILER_IS_INTEL) && \
         defined(PLATFORMSTL_OS_IS_WINDOWS))
#  define snprintf   _snprintf
# endif /* compiler / OS */
#endif /* !XCONTRACT_USING_SAFE_STR_FUNCTIONS */

/* /////////////////////////////////////////////////////////////////////////
 * Forward declarations
 */

static char const* xcontract_lookup_violation_type_string_(xContract_violation_type_t type, size_t* len);

#ifdef PLATFORMSTL_OS_IS_WINDOWS
# define xcontract_call_exit_()     ExitProcess(1)
#else /* ? PLATFORMSTL_OS_???? */
# define xcontract_call_exit_()     _exit(1)
#endif /* PLATFORMSTL_OS_???? */

/* /////////////////////////////////////////////////////////////////////////
 * API
 */

XCONTRACT_CALL(void) xContract_violationReport( char const*                 file
                                            ,   unsigned                    line
                                            ,   char const*                 function
                                            ,   char const*                 expression
                                            ,   char const*                 message
                                            ,   char const*                 qualifier
                                            ,   xContract_violation_type_t  type
                                            ,   int                         level)
{
    /* 1. formulate the termination messages */

    char        sz[2001];
    size_t      type_len = 0;
    char const* type_str = xcontract_lookup_violation_type_string_(type, &type_len);
    int         cch;

    STLSOFT_SUPPRESS_UNUSED(level);

#ifdef XCONTRACT_USING_SAFE_STR_FUNCTIONS
    cch = sprintf_s(
#else /* ? XCONTRACT_USING_SAFE_STR_FUNCTIONS */
    cch = snprintf(
#endif /* XCONTRACT_USING_SAFE_STR_FUNCTIONS */
                    &sz[0]
                ,   STLSOFT_NUM_ELEMENTS(sz) - 2
                ,   "Contract violation: %s at %s(%d)%s%s: %s%s%s%s%s"
                /* %s */
                ,   type_str
                /* %s(%d)%s%s */
                ,   file
                ,   line
                ,   (NULL == function) ? "" : " in function "
                ,   stlsoft_ns_qual(c_str_ptr_a)(function)
                /* %s%s%s */
                ,   message
                ,   (NULL == expression) ? "" : " - failed expression: "
                ,   stlsoft_ns_qual(c_str_ptr_a)(expression)
                /* %s%s */
                ,   (NULL == qualifier || '\0' == 0[qualifier]) ? "" : ": "
                ,   stlsoft_ns_qual(c_str_ptr_a)(qualifier)
                );

    if(cch < 0)
    {
        cch = STLSOFT_NUM_ELEMENTS(sz) - 2;
    }

    sz[cch]     = '\n';
    sz[cch + 1] = '\0';

    /* 2. write out the termination messages */

    fputs(sz, stderr);

#ifdef XCONTRACT_USE_WIN_OUTPUTDEBUGSTRING
    OutputDebugStringA(sz);
#endif /* XCONTRACT_USE_WIN_OUTPUTDEBUGSTRING */

    sz[cch] = '\0';

#ifdef XCONTRACT_USE_WIN_EVENTLOG
    {
        HANDLE hEvLog = RegisterEventSource(NULL, "xContract");

        if(NULL != hEvLog)
        {
            WORD        wType       =   EVENTLOG_ERROR_TYPE;
            WORD        category    =   XCONTRACT_WIN_EVENTLOG_CATEGORY;
            DWORD       eventId     =   XCONTRACT_WIN_EVENTLOG_EVENTID;
            PSID        lpUserSid   =   NULL;
            WORD        wNumStrings =   1;
            LPCSTR      entry       =   &sz[0];
            LPCSTR*     lpStrings   =   &entry;
            DWORD       dwDataSize  =   0;
            LPVOID      lpRawData   =   NULL;

            ReportEventA(   hEvLog
                        ,   wType
                        ,   category
                        ,   eventId
                        ,   lpUserSid
                        ,   wNumStrings
                        ,   dwDataSize
                        ,   lpStrings
                        ,   lpRawData);

            DeregisterEventSource(hEvLog);
        }
    }
#endif /* XCONTRACT_USE_WIN_EVENTLOG */

#ifdef XCONTRACT_USE_SYSLOG
    syslog(XCONTRACT_SYSLOG_FLAGS, "%s", sz);
#endif /* XCONTRACT_USE_SYSLOG */
}

XCONTRACT_CALL(void) xContract_contract_violation(  char const*                     file
                                                ,   unsigned                        line
                                                ,   char const*                     expression
                                                ,   char const*                     message
                                                ,   xContract_violation_type_t      type)
{
    xContract_contract_violation2(file, line, NULL, expression, message, NULL, type, 0, xContract_violationReport);
}


XCONTRACT_CALL(void) xContract_contract_violation2( char const*                     file
                                                ,   unsigned                        line
                                                ,   char const*                     function
                                                ,   char const*                     expression
                                                ,   char const*                     message
                                                ,   char const*                     qualifier
                                                ,   xContract_violation_type_t      type
                                                ,   int                             level
                                                ,   xContract_violationReport_fn_t  pfnReport)
{
    if(NULL == pfnReport)
    {
        pfnReport = xContract_violationReport;
    }

#ifdef XCONTRACT_COMPILER_SUPPORTS_SEH
    __try
    {
#endif /* XCONTRACT_COMPILER_SUPPORTS_SEH */

        /* 1. Report */

        (*pfnReport)(file, line, function, expression, message, qualifier, type, level);

#ifdef XCONTRACT_COMPILER_SUPPORTS_SEH
    }
    __finally
    {
#endif /* XCONTRACT_COMPILER_SUPPORTS_SEH */

        /* 2. Respond: terminate the process */

        xcontract_call_exit_();

#ifdef XCONTRACT_COMPILER_SUPPORTS_SEH
    }
#endif /* XCONTRACT_COMPILER_SUPPORTS_SEH */
}

/* /////////////////////////////////////////////////////////////////////////
 * Strings
 */

typedef struct xContract_type_string    xContract_type_string;
struct xContract_type_string
{
    xContract_violation_type_t  type;   /*!< The type code.     */
    char const*                 str;    /*!< The string.        */
    size_t                      len;    /*!< The string length. */
};

#define XCONTRACT_TYPE_STR_DECL(rc, desc)                                                                   \
                                                                                                            \
    static const char                   s_str##rc[] =   desc;                                               \
    static const xContract_type_string  s_rct##rc = { rc, s_str##rc, STLSOFT_NUM_ELEMENTS(s_str##rc) - 1 }

#define XCONTRACT_TYPE_STR_ENTRY(rc)                                                                        \
                                                                                                            \
    &s_rct##rc


static char const* xContract_lookup_type_(xContract_violation_type_t type, xContract_type_string const** mappings, size_t cMappings, size_t* len)
{
    /* Use Null Object (Variable) here for len, so do not need to check
     * elsewhere.
     */
    size_t  len_;

    if(NULL == len)
    {
        len = &len_;
    }

    /* Checked, indexed search. */
    if( type >= 0 &&
        type < xContract_maximum_value)
    {
        if(type == mappings[type]->type)
        {
            return (*len = mappings[type]->len, mappings[type]->str);
        }
    }

    /* Linear search. Should only be needed if order in
     * pantheios_LookupxContract_type_stringA_() messed up.
     */
    { size_t i; for(i = 0; i < cMappings; ++i)
    {
        if(type == mappings[i]->type)
        {
            return (*len = mappings[i]->len, mappings[i]->str);
        }
    }}

    return (*len = 0, "");
}

static char const* xcontract_lookup_violation_type_string_(xContract_violation_type_t type, size_t* len)
{
    XCONTRACT_TYPE_STR_DECL(xContract_unexpectedCondition       ,   "unexpected condition"          );
    XCONTRACT_TYPE_STR_DECL(xContract_precondition_logic        ,   "precondition (logic)"          );
    XCONTRACT_TYPE_STR_DECL(xContract_precondition_parameters   ,   "precondition (parameters)"     );
    XCONTRACT_TYPE_STR_DECL(xContract_postcondition_returnValue ,   "postcondition (return value)"  );
    XCONTRACT_TYPE_STR_DECL(xContract_postcondition_logic       ,   "postcondition (logic)"         );
    XCONTRACT_TYPE_STR_DECL(xContract_postcondition_parameters  ,   "postcondition (parameters)"    );
    XCONTRACT_TYPE_STR_DECL(xContract_invariant_class           ,   "class invariant"               );
    XCONTRACT_TYPE_STR_DECL(xContract_invariant_global          ,   "global invariant"              );
    XCONTRACT_TYPE_STR_DECL(xContract_staticData                ,   "static data"                   );
    XCONTRACT_TYPE_STR_DECL(xContract_intermediateAssumption    ,   "global invariant"              );


    static const xContract_type_string* s_strings[] = 
    {
        XCONTRACT_TYPE_STR_ENTRY(xContract_unexpectedCondition),
        XCONTRACT_TYPE_STR_ENTRY(xContract_precondition_logic),
        XCONTRACT_TYPE_STR_ENTRY(xContract_precondition_parameters),
        XCONTRACT_TYPE_STR_ENTRY(xContract_postcondition_returnValue),
        XCONTRACT_TYPE_STR_ENTRY(xContract_postcondition_logic),
        XCONTRACT_TYPE_STR_ENTRY(xContract_postcondition_parameters),
        XCONTRACT_TYPE_STR_ENTRY(xContract_invariant_class),
        XCONTRACT_TYPE_STR_ENTRY(xContract_invariant_global),
        XCONTRACT_TYPE_STR_ENTRY(xContract_staticData),
        XCONTRACT_TYPE_STR_ENTRY(xContract_intermediateAssumption),
    };

    return xContract_lookup_type_(type, s_strings, STLSOFT_NUM_ELEMENTS(s_strings), len);
}

XCONTRACT_CALL(char const*) xContract_getViolationTypeString(xContract_violation_type_t type)
{
    return xcontract_lookup_violation_type_string_(type, NULL);
}

XCONTRACT_CALL(size_t) xContract_getViolationTypeStringLength(xContract_violation_type_t type)
{
    size_t len;

    return (xcontract_lookup_violation_type_string_(type, &len), len);
}

/* ///////////////////////////// end of file //////////////////////////// */
