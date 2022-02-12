/* /////////////////////////////////////////////////////////////////////////
 * File:        xcontract/xcontract.h (formerly part of Synesis' internal test codebase)
 *
 * Purpose:     Main header file for the xContract library.
 *
 * Created:     2nd January 2001
 * Updated:     4th July 2010
 *
 * Home:        http://xcontract.org/
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


/** \file xcontract/xcontract.h
 *
 * [C, C++] Main header file for the <b>xContract</b> library.
 */

#ifndef XCONTRACT_INCL_XCONTRACT_H_XCONTRACT
#define XCONTRACT_INCL_XCONTRACT_H_XCONTRACT

#ifndef XCONTRACT_DOCUMENTATION_SKIP_SECTION
# define XCONTRACT_VER_XCONTRACT_H_XCONTRACT_MAJOR       2
# define XCONTRACT_VER_XCONTRACT_H_XCONTRACT_MINOR       5
# define XCONTRACT_VER_XCONTRACT_H_XCONTRACT_REVISION    4
# define XCONTRACT_VER_XCONTRACT_H_XCONTRACT_EDIT        181
#endif /* !XCONTRACT_DOCUMENTATION_SKIP_SECTION */

/* /////////////////////////////////////////////////////////////////////////
 * Version information
 */

/**
 * \def XCONTRACT_VER_MAJOR
 * The Major version number of the xContract library
 *
 * \def XCONTRACT_VER_MINOR
 * Minor version number of the xContract library
 *
 * \def XCONTRACT_VER_REVISION
 * The revision number of the xContract library
 *
 * \def XCONTRACT_VER
 * The composite version of the xContract library
 */

#define XCONTRACT_VER_MAJOR         0
#define XCONTRACT_VER_MINOR         3
#define XCONTRACT_VER_REVISION      7

#define XCONTRACT_VER               0x000307ff

/* /////////////////////////////////////////////////////////////////////////
 * Includes
 */

#ifdef __cplusplus
# if defined(UNIX) || \
     defined(unix) || \
     defined(unix__) || \
     defined(__unix) || \
     defined(__unix__)
 /* We include threading.h to prevent the definition of _REENTRANT standard
  * headers on some UNIX operating systems from confusing the feature
  * discrimination in UNIXSTL and having it think that we're multithreading
  * when we're not.
  */
#  include <unixstl/synch/util/features.h>
# endif /* UNIX */
#endif /* __cplusplus */

#include <stddef.h>

/* /////////////////////////////////////////////////////////////////////////
 * Namespace
 */

#if defined(_STLSOFT_NO_NAMESPACE) || \
    !defined(__cplusplus)
# define _XCONTRACT_NO_NAMESPACE
#endif /* _STLSOFT_NO_NAMESPACE */

#ifndef _XCONTRACT_NO_NAMESPACE
/** The xContract namespace
 *
 * All types and functions are defined within this namespace in C++
 */
namespace xcontract
{
#endif /* !_XCONTRACT_NO_NAMESPACE */

/* /////////////////////////////////////////////////////////////////////////
 * Documentation
 */

/** \defgroup group__application_layer Application Layer
 *
 * \remarks The names of these macros are intentionally long-winded in order
 *  that they are totally unambiguous. The intention is that you will define
 *  your own, shorter, macros specific to your library/application
 */

/** \defgroup group__application_defined Application-defined Functions
 *
 * These functions are defined by the application to modify the behaviour
 * of the library in response to contract enforcement violations
 */

/** \defgroup group__api_core Core API
 *
 * The core API functions are used by the \ref group__application_layer
 * macros
 */

/* /////////////////////////////////////////////////////////////////////////
 * Compatibility
 */

#ifndef XCONTRACT_DOCUMENTATION_SKIP_SECTION

 /* static cast */
# ifdef __cplusplus
#  define XCONTRACT_STATIC_CAST_(t, v)			static_cast< t>(v)
# else /* ? __cplusplus */
#  define XCONTRACT_STATIC_CAST_(t, v)			(t)(v)
# endif /* __cplusplus */

 /* __FUNCTION__ support */
# if defined(__BORLANDC__)
# elif defined(__COMO__)
#  define XCONTRACT_PPF_FUNCTION_SYMBOL_SUPPORT
# elif defined(__DMC__)
#  if __DMC__ >= 0x850
#   define XCONTRACT_PPF_FUNCTION_SYMBOL_SUPPORT
#  endif
# elif defined(__INTEL_COMPILER)
#  if __INTEL_COMPILER >= 700
#   define XCONTRACT_PPF_FUNCTION_SYMBOL_SUPPORT
#  endif
# elif defined(__MWERKS__)
#  if (__MWERKS__ & 0xFF00) >= 0x3000
#   define XCONTRACT_PPF_FUNCTION_SYMBOL_SUPPORT
#  endif
# elif defined(__WATCOMC__)
#  if __WATCOMC__ >= 1240
#   define XCONTRACT_PPF_FUNCTION_SYMBOL_SUPPORT
#  endif
# elif defined(__GNUC___)
#  if __GNUC___ >= 3
#   define XCONTRACT_PPF_FUNCTION_SYMBOL_SUPPORT
#  endif
# elif defined(__SUNPRO_C)
#  define XCONTRACT_PPF_func_SYMBOL_SUPPORT
# elif defined(__SUNPRO_CC)
#  define XCONTRACT_PPF_func_SYMBOL_SUPPORT
# elif defined(_MSC_VER)
#  if _MSC_VER >= 1200
#   define XCONTRACT_PPF_DECLSPEC_NORETURN_SYMBOL_SUPPORT
#  endif
#  if _MSC_VER >= 1300
#   define XCONTRACT_PPF_FUNCTION_SYMBOL_SUPPORT
#  endif
# endif

#endif /* !XCONTRACT_DOCUMENTATION_SKIP_SECTION */

/* /////////////////////////////////////////////////////////////////////////
 * Features
 */

#ifndef XCONTRACT_DOCUMENTATION_SKIP_SECTION

# if !defined(XCONTRACT_DECLSPEC)
#  define XCONTRACT_DECLSPEC
# endif /* !XCONTRACT_DECLSPEC */

# if defined(__cplusplus) || \
     defined(XCONTRACT_DOCUMENTATION_SKIP_SECTION)
#  define XCONTRACT_EXTERN_C                    extern "C"
# else /* ? __cplusplus */
#  define XCONTRACT_EXTERN_C                    extern
# endif /* !__cplusplus */

# if !defined(XCONTRACT_CALLCONV)
#  define XCONTRACT_CALLCONV
# endif /* !XCONTRACT_CALLCONV */

# define XCONTRACT_CALL(rt)                     XCONTRACT_DECLSPEC XCONTRACT_EXTERN_C rt XCONTRACT_CALLCONV

# ifdef XCONTRACT_PPF_FUNCTION_SYMBOL_SUPPORT
#  define XCONTRACT_GET_FUNCTION_()             __FUNCTION__
# elif defined(XCONTRACT_PPF_func_SYMBOL_SUPPORT)
#  define XCONTRACT_GET_FUNCTION_()             __func__
# else /* ? XCONTRACT_PPF_FUNC(TION)_SYMBOL_SUPPORT */
#  define XCONTRACT_GET_FUNCTION_()             XCONTRACT_STATIC_CAST_(char const*, 0)
# endif /* XCONTRACT_PPF_FUNC(TION)_SYMBOL_SUPPORT */

# ifndef _XCONTRACT_NO_NAMESPACE
#  define XCONTRACT_NS_QUAL(sym)                xcontract::sym
# else /* ? _XCONTRACT_NO_NAMESPACE */
#  define XCONTRACT_NS_QUAL(sym)                sym
# endif /* _XCONTRACT_NO_NAMESPACE */

# ifdef XCONTRACT_PPF_DECLSPEC_NORETURN_SYMBOL_SUPPORT
#  define XCONTRACT_NORETURN                    __declspec(noreturn)
# else /* ? compiler */
#  define XCONTRACT_NORETURN
# endif /* compiler */

#endif /* !XCONTRACT_DOCUMENTATION_SKIP_SECTION */

/* /////////////////////////////////////////////////////////////////////////
 * Modes:
 *
 *  - Debug:
 *      - pre-condition
 *          - logic             - on
 *          - parameters        - on
 *      - post-condition
 *          - logic             - on
 *          - ret val           - on
 *          - out parameters    - on
 *      - invariant
 *          - class             - on
 *          - global            - on
 *      - unexpected condition  - on
 *
 *  - Test:
 *      - pre-condition
 *          - logic             - on
 *          - parameters        - on
 *      - post-condition
 *          - logic             - on
 *          - ret val           - on
 *          - out parameters    - off
 *      - invariant
 *          - class             - on
 *          - global            - on
 *      - unexpected condition  - on
 *
 *  - Release:
 *      - pre-condition
 *          - logic             - on
 *          - parameters        - on
 *      - post-condition
 *          - logic             - off
 *          - ret val           - off
 *          - out parameters    - off
 *      - invariant
 *          - class             - on
 *          - global            - on
 *      - unexpected condition  - on
 *
 *  - Off:
 *      - pre-condition
 *          - logic             - off
 *          - parameters        - off
 *      - post-condition
 *          - logic             - off
 *          - ret val           - off
 *          - out parameters    - off
 *      - invariant
 *          - class             - off
 *          - global            - off
 *      - unexpected condition  - off
 *
 */

#if !defined(XCONTRACT_NO_STOCK_MODES) && \
    !defined(XCONTRACT_DOCUMENTATION_SKIP_SECTION)

# if !defined(XCONTRACT_MODE_DEBUG) && \
     !defined(XCONTRACT_MODE_TEST) && \
     !defined(XCONTRACT_MODE_RELEASE) && \
     !defined(XCONTRACT_MODE_OFF)
  /* No mode specified, so we deduce from presence/absence of debug
   * symbols
   */
#  if !defined(NDEBUG)
#   define XCONTRACT_MODE_DEBUG
#  else /* ? NDEBUG */
#   define XCONTRACT_MODE_RELEASE
#  endif /* NDEBUG */
# endif /* XCONTRACT_MODE_???? */

# if defined(XCONTRACT_MODE_DEBUG)
   /* Everything is on in Debug */
#  define XCONTRACT_ENFORCING_UNEXPECTED_CONDITION              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_0              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_1              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_2              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_3              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_4              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_0         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_1         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_2         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_3         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_4         (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_RETURN_0            (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_RETURN_1            (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_RETURN_2            (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_RETURN_3            (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_RETURN_4            (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_0             (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_1             (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_2             (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_3             (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_4             (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_0        (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_1        (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_2        (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_3        (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_4        (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_0                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_1                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_2                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_3                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_4                 (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_0                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_1                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_2                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_3                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_4                (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_0                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_1                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_2                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_3                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_4                     (1)
#  define XCONTRACT_ENFORCING_ASSUMPTION_0                      (1)
#  define XCONTRACT_ENFORCING_ASSUMPTION_1                      (1)
#  define XCONTRACT_ENFORCING_ASSUMPTION_2                      (1)
#  define XCONTRACT_ENFORCING_ASSUMPTION_3                      (1)
#  define XCONTRACT_ENFORCING_ASSUMPTION_4                      (1)
# endif /* XCONTRACT_MODE_DEBUG */

# if defined(XCONTRACT_MODE_TEST)
   /* Everything except postcondition out-parameters is on in Test */
#  define XCONTRACT_ENFORCING_UNEXPECTED_CONDITION              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_0              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_1              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_2              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_3              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_4              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_0         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_1         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_2         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_3         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_4         (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_RETURN_0            (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_RETURN_1            (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_RETURN_2            (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_RETURN_3            (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_RETURN_4            (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_0             (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_1             (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_2             (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_3             (1)
#  define XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_4             (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_0                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_1                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_2                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_3                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_4                 (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_0                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_1                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_2                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_3                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_4                (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_0                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_1                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_2                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_3                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_4                     (1)
#  define XCONTRACT_ENFORCING_ASSUMPTION_0                      (1)
#  define XCONTRACT_ENFORCING_ASSUMPTION_1                      (1)
#  define XCONTRACT_ENFORCING_ASSUMPTION_2                      (1)
#  define XCONTRACT_ENFORCING_ASSUMPTION_3                      (1)
#  define XCONTRACT_ENFORCING_ASSUMPTION_4                      (1)
# endif /* XCONTRACT_MODE_TEST */

# if defined(XCONTRACT_MODE_RELEASE)
   /* Everything except postconditions and intermediate assumptions are on in Release */
#  define XCONTRACT_ENFORCING_UNEXPECTED_CONDITION              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_0              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_1              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_2              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_3              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_LOGIC_4              (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_0         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_1         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_2         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_3         (1)
#  define XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_4         (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_0                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_1                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_2                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_3                 (1)
#  define XCONTRACT_ENFORCING_CLASS_INVARIANT_4                 (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_0                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_1                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_2                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_3                (1)
#  define XCONTRACT_ENFORCING_GLOBAL_INVARIANT_4                (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_0                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_1                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_2                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_3                     (1)
#  define XCONTRACT_ENFORCING_STATIC_DATA_4                     (1)
# endif /* XCONTRACT_MODE_RELEASE */

# if defined(XCONTRACT_MODE_OFF)
   /* Nothing is on in Off */
# endif /* XCONTRACT_MODE_OFF */


#endif /* !XCONTRACT_NO_STOCK_MODES */

/* /////////////////////////////////////////////////////////////////////////
 * Types
 */

/** The possible violation types
 *
 * \see xContract_getViolationTypeString, xContract_getViolationTypeStringLength
 */
enum xContract_violation_type_t
{
        xContract_unexpectedCondition       /*!< Indicates an unexpected condition */
    ,   xContract_precondition_logic        /*!< Indicates a precondition logic violation */
    ,   xContract_precondition_parameters   /*!< Indicates a precondition parameter violation */
    ,   xContract_postcondition_returnValue /*!< Indicates a postcondition return-value violation */
    ,   xContract_postcondition_logic       /*!< Indicates a postcondition logic violation */
    ,   xContract_postcondition_parameters  /*!< Indicates a postcondition parameter violation */
    ,   xContract_invariant_class           /*!< Indicates a class invariant violation */
    ,   xContract_invariant_global          /*!< Indicates a global invariant violation */
    ,   xContract_staticData                /*!< Indicates that static data is in violatation */
    ,   xContract_intermediateAssumption    /*!< Indicates an immediate assumption violatation */

#ifndef XCONTRACT_DOCUMENTATION_SKIP_SECTION
    ,   xContract_maximum_value
#endif /* !XCONTRACT_DOCUMENTATION_SKIP_SECTION */
};
#ifndef __cplusplus
typedef enum xContract_violation_type_t         xContract_violation_type_t;
#endif /* !__cplusplus */
#ifndef _XCONTRACT_NO_NAMESPACE
/** Alias for xContract_violation_type_t
 */
typedef xContract_violation_type_t              violation_type_t;
#endif /* !_XCONTRACT_NO_NAMESPACE */

/** Reporting function prototype
 */
typedef void (XCONTRACT_CALLCONV *xContract_violationReport_fn_t)(
	char const*                 file
,   unsigned                    line
,   char const*                 function
,   char const*                 expression
,   char const*                 message
,   char const*                 qualifier
,   xContract_violation_type_t  type
,   int                         level
);

#ifndef _XCONTRACT_NO_NAMESPACE
/** Alias for xContract_violationReport_fn_t
 */
typedef xContract_violationReport_fn_t          violationReport_fn_t;
#endif /* !_XCONTRACT_NO_NAMESPACE */

/* /////////////////////////////////////////////////////////////////////////
 * Macros
 *
 * NOTE: These macro names are deliberately unappealingly large, so that
 * you'll be motivated to define your own, company/project-specific macros,
 * that do exactly what you want, rather than the very general facilities
 * provided here.
 */



#ifndef XCONTRACT_DOCUMENTATION_SKIP_SECTION

# ifdef XCONTRACT_VERIFY_CONDITION_AT_RUNTIME
#  define XCONTRACT_ENFORCE_CONDITION_VERIFIER_(type, level, expr, msg)     XCONTRACT_NS_QUAL(xContract_isConditionVerified)(type, level, __FILE__, __LINE__, XCONTRACT_GET_FUNCTION_(), expr)
# else /* ? XCONTRACT_VERIFY_CONDITION_AT_RUNTIME */
#  define XCONTRACT_ENFORCE_CONDITION_VERIFIER_(type, level, expr, msg)     (1)
# endif /* XCONTRACT_VERIFY_CONDITION_AT_RUNTIME */

# ifndef XCONTRACT_CUSTOM_REPORTER
#  define XCONTRACT_CUSTOM_REPORTER                                         XCONTRACT_NS_QUAL(xContract_violationReport)
# endif /* !XCONTRACT_CUSTOM_REPORTER */

# define XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(type, level, expr, msg)                \
                                                                                        \
    (((0 == XCONTRACT_ENFORCE_CONDITION_VERIFIER_(type, level, #expr, msg)) || (expr))  \
        ? XCONTRACT_STATIC_CAST_(void, 0)                                                  \
        : XCONTRACT_NS_QUAL(xContract_contract_violation2)(__FILE__, __LINE__, XCONTRACT_GET_FUNCTION_(), #expr, msg, NULL, type, level, XCONTRACT_CUSTOM_REPORTER))

#endif /* !XCONTRACT_DOCUMENTATION_SKIP_SECTION */


/* Unexpected condition */

/** \def XCONTRACT_ENFORCE_UNEXPECTED_CONDITION(msg)
 *
 * Expresses an unexpected condition enforcement
 *
 * \ingroup group__application_layer
 *
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_UNEXPECTED_CONDITION
# define XCONTRACT_ENFORCE_UNEXPECTED_CONDITION(msg)                XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_unexpectedCondition), 0, NULL, msg)
#else /* ? XCONTRACT_ENFORCING_UNEXPECTED_CONDITION */
# define XCONTRACT_ENFORCE_UNEXPECTED_CONDITION(msg)                ((void)0)
#endif /* XCONTRACT_ENFORCING_UNEXPECTED_CONDITION */


/* Precondition (logic) */

#if 0
# ifdef XCONTRACT_ENFORCING_PRECONDITION_LOGIC_LEVEL > 0
#  define XCONTRACT_ENFORCE_PRECONDITION_LOGIC(expr, msg, level)    XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_precondition_logic), 0, expr, msg)
# else /* ? XCONTRACT_ENFORCING_PRECONDITION_LOGIC_LEVEL */
#  define XCONTRACT_ENFORCE_PRECONDITION_LOGIC(expr, msg, level)    ((void)0)
# endif /* XCONTRACT_ENFORCING_PRECONDITION_LOGIC_LEVEL */
#endif /* 0 */

/** \def XCONTRACT_ENFORCE_PRECONDITION_LOGIC_0(expr, msg)
 *
 * Expresses a precondition (logic) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_PRECONDITION_LOGIC_0
# define XCONTRACT_ENFORCE_PRECONDITION_LOGIC_0(expr, msg)          XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_precondition_logic), 0, expr, msg)
#else /* ? XCONTRACT_ENFORCING_PRECONDITION_LOGIC_0 */
# define XCONTRACT_ENFORCE_PRECONDITION_LOGIC_0(expr, msg)          ((void)0)
#endif /* XCONTRACT_ENFORCING_PRECONDITION_LOGIC_0 */

/** \def XCONTRACT_ENFORCE_PRECONDITION_LOGIC_1(expr, msg)
 *
 * Expresses a precondition (logic) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_PRECONDITION_LOGIC_1
# define XCONTRACT_ENFORCE_PRECONDITION_LOGIC_1(expr, msg)          XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_precondition_logic), 1, expr, msg)
#else /* ? XCONTRACT_ENFORCING_PRECONDITION_LOGIC_1 */
# define XCONTRACT_ENFORCE_PRECONDITION_LOGIC_1(expr, msg)          ((void)0)
#endif /* XCONTRACT_ENFORCING_PRECONDITION_LOGIC_1 */

/** \def XCONTRACT_ENFORCE_PRECONDITION_LOGIC_2(expr, msg)
 *
 * Expresses a precondition (logic) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_PRECONDITION_LOGIC_2
# define XCONTRACT_ENFORCE_PRECONDITION_LOGIC_2(expr, msg)          XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_precondition_logic), 2, expr, msg)
#else /* ? XCONTRACT_ENFORCING_PRECONDITION_LOGIC_2 */
# define XCONTRACT_ENFORCE_PRECONDITION_LOGIC_2(expr, msg)          ((void)0)
#endif /* XCONTRACT_ENFORCING_PRECONDITION_LOGIC_2 */

/** \def XCONTRACT_ENFORCE_PRECONDITION_LOGIC_3(expr, msg)
 *
 * Expresses a precondition (logic) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_PRECONDITION_LOGIC_3
# define XCONTRACT_ENFORCE_PRECONDITION_LOGIC_3(expr, msg)          XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_precondition_logic), 3, expr, msg)
#else /* ? XCONTRACT_ENFORCING_PRECONDITION_LOGIC_3 */
# define XCONTRACT_ENFORCE_PRECONDITION_LOGIC_3(expr, msg)          ((void)0)
#endif /* XCONTRACT_ENFORCING_PRECONDITION_LOGIC_3 */

/** \def XCONTRACT_ENFORCE_PRECONDITION_LOGIC_4(expr, msg)
 *
 * Expresses a precondition (logic) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_PRECONDITION_LOGIC_4
# define XCONTRACT_ENFORCE_PRECONDITION_LOGIC_4(expr, msg)          XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_precondition_logic), 4, expr, msg)
#else /* ? XCONTRACT_ENFORCING_PRECONDITION_LOGIC_4 */
# define XCONTRACT_ENFORCE_PRECONDITION_LOGIC_4(expr, msg)          ((void)0)
#endif /* XCONTRACT_ENFORCING_PRECONDITION_LOGIC_4 */


/* Precondition (parameters) */

/** \def XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_0(expr, msg)
 *
 * Expresses a precondition (parameters) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_0
# define XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_0(expr, msg)     XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_precondition_parameters), 0, expr, msg)
#else /* ? XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_0 */
# define XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_0(expr, msg)     ((void)0)
#endif /* XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_0 */

/** \def XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_1(expr, msg)
 *
 * Expresses a precondition (parameters) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_1
# define XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_1(expr, msg)     XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_precondition_parameters), 1, expr, msg)
#else /* ? XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_1 */
# define XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_1(expr, msg)     ((void)0)
#endif /* XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_1 */

/** \def XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_2(expr, msg)
 *
 * Expresses a precondition (parameters) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_2
# define XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_2(expr, msg)     XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_precondition_parameters), 2, expr, msg)
#else /* ? XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_2 */
# define XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_2(expr, msg)     ((void)0)
#endif /* XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_2 */

/** \def XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_3(expr, msg)
 *
 * Expresses a precondition (parameters) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_3
# define XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_3(expr, msg)     XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_precondition_parameters), 3, expr, msg)
#else /* ? XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_3 */
# define XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_3(expr, msg)     ((void)0)
#endif /* XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_3 */

/** \def XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_4(expr, msg)
 *
 * Expresses a precondition (parameters) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_4
# define XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_4(expr, msg)     XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_precondition_parameters), 4, expr, msg)
#else /* ? XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_4 */
# define XCONTRACT_ENFORCE_PRECONDITION_PARAMETERS_4(expr, msg)     ((void)0)
#endif /* XCONTRACT_ENFORCING_PRECONDITION_PARAMETERS_4 */


/* Postcondition (return value) */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_RETURN_0(expr, msg)
 *
 * Expresses a postcondition (return value) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_RETURN_0
# define XCONTRACT_ENFORCE_POSTCONDITION_RETURN_0(expr, msg)         XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_returnValue), 0, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_RETURN_0 */
# define XCONTRACT_ENFORCE_POSTCONDITION_RETURN_0(expr, msg)         ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_RETURN_0 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_RETURN_1(expr, msg)
 *
 * Expresses a postcondition (return value) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_RETURN_1
# define XCONTRACT_ENFORCE_POSTCONDITION_RETURN_1(expr, msg)         XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_returnValue), 1, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_RETURN_1 */
# define XCONTRACT_ENFORCE_POSTCONDITION_RETURN_1(expr, msg)         ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_RETURN_1 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_RETURN_2(expr, msg)
 *
 * Expresses a postcondition (return value) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_RETURN_2
# define XCONTRACT_ENFORCE_POSTCONDITION_RETURN_2(expr, msg)         XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_returnValue), 2, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_RETURN_2 */
# define XCONTRACT_ENFORCE_POSTCONDITION_RETURN_2(expr, msg)         ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_RETURN_2 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_RETURN_3(expr, msg)
 *
 * Expresses a postcondition (return value) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_RETURN_3
# define XCONTRACT_ENFORCE_POSTCONDITION_RETURN_3(expr, msg)         XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_returnValue), 3, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_RETURN_3 */
# define XCONTRACT_ENFORCE_POSTCONDITION_RETURN_3(expr, msg)         ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_RETURN_3 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_RETURN_4(expr, msg)
 *
 * Expresses a postcondition (return value) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_RETURN_4
# define XCONTRACT_ENFORCE_POSTCONDITION_RETURN_4(expr, msg)         XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_returnValue), 4, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_RETURN_4 */
# define XCONTRACT_ENFORCE_POSTCONDITION_RETURN_4(expr, msg)         ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_RETURN_4 */


/* Postcondition (logic) */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_0(expr, msg)
 *
 * Expresses a postcondition (logic) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_0
# define XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_0(expr, msg)         XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_logic), 0, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_0 */
# define XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_0(expr, msg)         ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_0 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_1(expr, msg)
 *
 * Expresses a postcondition (logic) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_1
# define XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_1(expr, msg)         XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_logic), 1, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_1 */
# define XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_1(expr, msg)         ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_1 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_2(expr, msg)
 *
 * Expresses a postcondition (logic) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_2
# define XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_2(expr, msg)         XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_logic), 2, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_2 */
# define XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_2(expr, msg)         ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_2 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_3(expr, msg)
 *
 * Expresses a postcondition (logic) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_3
# define XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_3(expr, msg)         XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_logic), 3, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_3 */
# define XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_3(expr, msg)         ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_3 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_4(expr, msg)
 *
 * Expresses a postcondition (logic) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_4
# define XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_4(expr, msg)         XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_logic), 4, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_4 */
# define XCONTRACT_ENFORCE_POSTCONDITION_LOGIC_4(expr, msg)         ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_LOGIC_4 */


/* Postcondition (parameters) */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_0(expr, msg)
 *
 * Expresses a postcondition (parameters) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_0
# define XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_0(expr, msg)    XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_parameters), 0, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_0 */
# define XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_0(expr, msg)    ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_0 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_1(expr, msg)
 *
 * Expresses a postcondition (parameters) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_1
# define XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_1(expr, msg)    XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_parameters), 1, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_1 */
# define XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_1(expr, msg)    ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_1 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_2(expr, msg)
 *
 * Expresses a postcondition (parameters) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_2
# define XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_2(expr, msg)    XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_parameters), 2, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_2 */
# define XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_2(expr, msg)    ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_2 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_3(expr, msg)
 *
 * Expresses a postcondition (parameters) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_3
# define XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_3(expr, msg)    XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_parameters), 3, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_3 */
# define XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_3(expr, msg)    ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_3 */

/** \def XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_4(expr, msg)
 *
 * Expresses a postcondition (parameters) enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_4
# define XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_4(expr, msg)    XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_postcondition_parameters), 4, expr, msg)
#else /* ? XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_4 */
# define XCONTRACT_ENFORCE_POSTCONDITION_PARAMETERS_4(expr, msg)    ((void)0)
#endif /* XCONTRACT_ENFORCING_POSTCONDITION_PARAMETERS_4 */


/* Invariant (class) */

/** \def XCONTRACT_ENFORCE_CLASS_INVARIANT_0(expr, msg)
 *
 * Expresses class invariant enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_CLASS_INVARIANT_0
# define XCONTRACT_ENFORCE_CLASS_INVARIANT_0(expr, msg)             XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_invariant_class), 0, expr, msg)
#else /* ? XCONTRACT_ENFORCING_CLASS_INVARIANT_0 */
# define XCONTRACT_ENFORCE_CLASS_INVARIANT_0(expr, msg)             ((void)0)
#endif /* XCONTRACT_ENFORCING_CLASS_INVARIANT_0 */

/** \def XCONTRACT_ENFORCE_CLASS_INVARIANT_1(expr, msg)
 *
 * Expresses class invariant enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_CLASS_INVARIANT_1
# define XCONTRACT_ENFORCE_CLASS_INVARIANT_1(expr, msg)             XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_invariant_class), 1, expr, msg)
#else /* ? XCONTRACT_ENFORCING_CLASS_INVARIANT_1 */
# define XCONTRACT_ENFORCE_CLASS_INVARIANT_1(expr, msg)             ((void)0)
#endif /* XCONTRACT_ENFORCING_CLASS_INVARIANT_1 */

/** \def XCONTRACT_ENFORCE_CLASS_INVARIANT_2(expr, msg)
 *
 * Expresses class invariant enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_CLASS_INVARIANT_2
# define XCONTRACT_ENFORCE_CLASS_INVARIANT_2(expr, msg)             XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_invariant_class), 2, expr, msg)
#else /* ? XCONTRACT_ENFORCING_CLASS_INVARIANT_2 */
# define XCONTRACT_ENFORCE_CLASS_INVARIANT_2(expr, msg)             ((void)0)
#endif /* XCONTRACT_ENFORCING_CLASS_INVARIANT_2 */

/** \def XCONTRACT_ENFORCE_CLASS_INVARIANT_3(expr, msg)
 *
 * Expresses class invariant enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_CLASS_INVARIANT_3
# define XCONTRACT_ENFORCE_CLASS_INVARIANT_3(expr, msg)             XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_invariant_class), 3, expr, msg)
#else /* ? XCONTRACT_ENFORCING_CLASS_INVARIANT_3 */
# define XCONTRACT_ENFORCE_CLASS_INVARIANT_3(expr, msg)             ((void)0)
#endif /* XCONTRACT_ENFORCING_CLASS_INVARIANT_3 */

/** \def XCONTRACT_ENFORCE_CLASS_INVARIANT_4(expr, msg)
 *
 * Expresses class invariant enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_CLASS_INVARIANT_4
# define XCONTRACT_ENFORCE_CLASS_INVARIANT_4(expr, msg)             XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_invariant_class), 4, expr, msg)
#else /* ? XCONTRACT_ENFORCING_CLASS_INVARIANT_4 */
# define XCONTRACT_ENFORCE_CLASS_INVARIANT_4(expr, msg)             ((void)0)
#endif /* XCONTRACT_ENFORCING_CLASS_INVARIANT_4 */


/* Invariant (global) */

/** \def XCONTRACT_ENFORCE_GLOBAL_INVARIANT_0(expr, msg)
 *
 * Expresses global invariant enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_GLOBAL_INVARIANT_0
# define XCONTRACT_ENFORCE_GLOBAL_INVARIANT_0(expr, msg)            XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_invariant_global), 0, expr, msg)
#else /* ? XCONTRACT_ENFORCING_GLOBAL_INVARIANT_0 */
# define XCONTRACT_ENFORCE_GLOBAL_INVARIANT_0(expr, msg)            ((void)0)
#endif /* XCONTRACT_ENFORCING_GLOBAL_INVARIANT_0 */

/** \def XCONTRACT_ENFORCE_GLOBAL_INVARIANT_1(expr, msg)
 *
 * Expresses global invariant enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_GLOBAL_INVARIANT_1
# define XCONTRACT_ENFORCE_GLOBAL_INVARIANT_1(expr, msg)            XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_invariant_global), 1, expr, msg)
#else /* ? XCONTRACT_ENFORCING_GLOBAL_INVARIANT_1 */
# define XCONTRACT_ENFORCE_GLOBAL_INVARIANT_1(expr, msg)            ((void)0)
#endif /* XCONTRACT_ENFORCING_GLOBAL_INVARIANT_1 */

/** \def XCONTRACT_ENFORCE_GLOBAL_INVARIANT_2(expr, msg)
 *
 * Expresses global invariant enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_GLOBAL_INVARIANT_2
# define XCONTRACT_ENFORCE_GLOBAL_INVARIANT_2(expr, msg)            XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_invariant_global), 2, expr, msg)
#else /* ? XCONTRACT_ENFORCING_GLOBAL_INVARIANT_2 */
# define XCONTRACT_ENFORCE_GLOBAL_INVARIANT_2(expr, msg)            ((void)0)
#endif /* XCONTRACT_ENFORCING_GLOBAL_INVARIANT_2 */

/** \def XCONTRACT_ENFORCE_GLOBAL_INVARIANT_3(expr, msg)
 *
 * Expresses global invariant enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_GLOBAL_INVARIANT_3
# define XCONTRACT_ENFORCE_GLOBAL_INVARIANT_3(expr, msg)            XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_invariant_global), 3, expr, msg)
#else /* ? XCONTRACT_ENFORCING_GLOBAL_INVARIANT_3 */
# define XCONTRACT_ENFORCE_GLOBAL_INVARIANT_3(expr, msg)            ((void)0)
#endif /* XCONTRACT_ENFORCING_GLOBAL_INVARIANT_3 */

/** \def XCONTRACT_ENFORCE_GLOBAL_INVARIANT_4(expr, msg)
 *
 * Expresses global invariant enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_GLOBAL_INVARIANT_4
# define XCONTRACT_ENFORCE_GLOBAL_INVARIANT_4(expr, msg)            XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_invariant_global), 4, expr, msg)
#else /* ? XCONTRACT_ENFORCING_GLOBAL_INVARIANT_4 */
# define XCONTRACT_ENFORCE_GLOBAL_INVARIANT_4(expr, msg)            ((void)0)
#endif /* XCONTRACT_ENFORCING_GLOBAL_INVARIANT_4 */


/* Static data */

/** \def XCONTRACT_ENFORCE_STATIC_DATA_0(expr, msg)
 *
 * Expresses static data enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_STATIC_DATA_0
# define XCONTRACT_ENFORCE_STATIC_DATA_0(expr, msg)                 XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_staticData), 0, expr, msg)
#else /* ? XCONTRACT_ENFORCING_STATIC_DATA_0 */
# define XCONTRACT_ENFORCE_STATIC_DATA_0(expr, msg)                 ((void)0)
#endif /* XCONTRACT_ENFORCING_STATIC_DATA_0 */

/** \def XCONTRACT_ENFORCE_STATIC_DATA_1(expr, msg)
 *
 * Expresses static data enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_STATIC_DATA_1
# define XCONTRACT_ENFORCE_STATIC_DATA_1(expr, msg)                 XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_staticData), 1, expr, msg)
#else /* ? XCONTRACT_ENFORCING_STATIC_DATA_1 */
# define XCONTRACT_ENFORCE_STATIC_DATA_1(expr, msg)                 ((void)0)
#endif /* XCONTRACT_ENFORCING_STATIC_DATA_1 */

/** \def XCONTRACT_ENFORCE_STATIC_DATA_2(expr, msg)
 *
 * Expresses static data enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_STATIC_DATA_2
# define XCONTRACT_ENFORCE_STATIC_DATA_2(expr, msg)                 XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_staticData), 2, expr, msg)
#else /* ? XCONTRACT_ENFORCING_STATIC_DATA_2 */
# define XCONTRACT_ENFORCE_STATIC_DATA_2(expr, msg)                 ((void)0)
#endif /* XCONTRACT_ENFORCING_STATIC_DATA_2 */

/** \def XCONTRACT_ENFORCE_STATIC_DATA_3(expr, msg)
 *
 * Expresses static data enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_STATIC_DATA_3
# define XCONTRACT_ENFORCE_STATIC_DATA_3(expr, msg)                 XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_staticData), 3, expr, msg)
#else /* ? XCONTRACT_ENFORCING_STATIC_DATA_3 */
# define XCONTRACT_ENFORCE_STATIC_DATA_3(expr, msg)                 ((void)0)
#endif /* XCONTRACT_ENFORCING_STATIC_DATA_3 */

/** \def XCONTRACT_ENFORCE_STATIC_DATA_4(expr, msg)
 *
 * Expresses static data enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_STATIC_DATA_4
# define XCONTRACT_ENFORCE_STATIC_DATA_4(expr, msg)                 XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_staticData), 4, expr, msg)
#else /* ? XCONTRACT_ENFORCING_STATIC_DATA_4 */
# define XCONTRACT_ENFORCE_STATIC_DATA_4(expr, msg)                 ((void)0)
#endif /* XCONTRACT_ENFORCING_STATIC_DATA_4 */


/* Intermediate assumption */

/** \def XCONTRACT_ENFORCE_ASSUMPTION_0(expr, msg)
 *
 * Expresses an intermediate assumption enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_ASSUMPTION_0
# define XCONTRACT_ENFORCE_ASSUMPTION_0(expr, msg)                  XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_intermediateAssumption), 0, expr, msg)
#else /* ? XCONTRACT_ENFORCING_ASSUMPTION_0 */
# define XCONTRACT_ENFORCE_ASSUMPTION_0(expr, msg)                  ((void)0)
#endif /* XCONTRACT_ENFORCING_ASSUMPTION_0 */

/** \def XCONTRACT_ENFORCE_ASSUMPTION_1(expr, msg)
 *
 * Expresses an intermediate assumption enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_ASSUMPTION_1
# define XCONTRACT_ENFORCE_ASSUMPTION_1(expr, msg)                  XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_intermediateAssumption), 1, expr, msg)
#else /* ? XCONTRACT_ENFORCING_ASSUMPTION_1 */
# define XCONTRACT_ENFORCE_ASSUMPTION_1(expr, msg)                  ((void)0)
#endif /* XCONTRACT_ENFORCING_ASSUMPTION_1 */

/** \def XCONTRACT_ENFORCE_ASSUMPTION_2(expr, msg)
 *
 * Expresses an intermediate assumption enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_ASSUMPTION_2
# define XCONTRACT_ENFORCE_ASSUMPTION_2(expr, msg)                  XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_intermediateAssumption), 2, expr, msg)
#else /* ? XCONTRACT_ENFORCING_ASSUMPTION_2 */
# define XCONTRACT_ENFORCE_ASSUMPTION_2(expr, msg)                  ((void)0)
#endif /* XCONTRACT_ENFORCING_ASSUMPTION_2 */

/** \def XCONTRACT_ENFORCE_ASSUMPTION_3(expr, msg)
 *
 * Expresses an intermediate assumption enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_ASSUMPTION_3
# define XCONTRACT_ENFORCE_ASSUMPTION_3(expr, msg)                  XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_intermediateAssumption), 3, expr, msg)
#else /* ? XCONTRACT_ENFORCING_ASSUMPTION_3 */
# define XCONTRACT_ENFORCE_ASSUMPTION_3(expr, msg)                  ((void)0)
#endif /* XCONTRACT_ENFORCING_ASSUMPTION_3 */

/** \def XCONTRACT_ENFORCE_ASSUMPTION_4(expr, msg)
 *
 * Expresses an intermediate assumption enforcement
 *
 * \ingroup group__application_layer
 *
 * \param expr The expression whose truth is enforced
 * \param msg The message associated with the enforcement
 */
#ifdef XCONTRACT_ENFORCING_ASSUMPTION_4
# define XCONTRACT_ENFORCE_ASSUMPTION_4(expr, msg)                  XCONTRACT_ENFORCE_CONDITION_TYPE_LEVEL_(XCONTRACT_NS_QUAL(xContract_intermediateAssumption), 4, expr, msg)
#else /* ? XCONTRACT_ENFORCING_ASSUMPTION_4 */
# define XCONTRACT_ENFORCE_ASSUMPTION_4(expr, msg)                  ((void)0)
#endif /* XCONTRACT_ENFORCING_ASSUMPTION_4 */

/* /////////////////////////////////////////////////////////////////////////
 * Application-defined functions
 */

/** Application-defined function that determines whether a test for a given
 * violation type at a given level should be verified.
 *
 * \ingroup group__application_defined
 *
 * \param type The violation type
 * \param level The level
 * \param file The file in which the enforcement is expressed
 * \param line The file on which the enforcement is expressed
 * \param function The file within which the enforcement is expressed, May
 *   be NULL.
 * \param expression The enforcement expression
 *
 * \return A value that controls whether this instance
 * \retval 0 The condition will not be verified
 * \retval !=0 The condition will be verified
 */
XCONTRACT_CALL(int) xContract_isConditionVerified(
    xContract_violation_type_t  type
,   int                         level
,   char const*                 file
,   unsigned                    line
,   char const*                 function
,   char const*                 expression
);

/* /////////////////////////////////////////////////////////////////////////
 * API functions
 */

#ifndef XCONTRACT_DOCUMENTATION_SKIP_SECTION
XCONTRACT_CALL(int) xContract_doWhileCondition(char const* );
#endif /* !XCONTRACT_DOCUMENTATION_SKIP_SECTION */

/** Stock reporting function, invoked when a custom function is not
 *   specified
 *
 * \ingroup group__api_core
 *
 * \param file The file in which the violation occurred
 * \param line The line at which the violation occurred
 * \param function The function within which the violation occurred. Will be
 *   <code>NULL</code> for compilers that do not support the standard
 *   <code>__FUNCTION__</code> symbol
 * \param expression The expression that caused the violation
 * \param message The message associated with the enforcement that
 *   experienced the violation
 * \param qualifier A qualifier to the message. May be <code>NULL</code>
 * \param type The violation type (one of xContract_violation_type_t)
 * \param level The abstraction level of the enforcement
 *
 * \remarks This function writes a detailed informative statement about the
 *   contract violation to the standard error stream, and to appropriate
 *   platform-specific outputs. On UNIX it writes to <code>syslog()</code>,
 *   and on Windows it writes to the system debugger (via
 *   <code>OutputDebugString()</code>) and to the Event Log.
 */
XCONTRACT_CALL(void) xContract_violationReport(
    char const                      *file
,   unsigned                        line
,   char const                      *function
,   char const                      *expression
,   char const                      *message
,   char const                      *qualifier
,   xContract_violation_type_t      type
,   int                             level
);

/** Invoked to report and respond to a contract violation
 *
 * \ingroup group__api_core
 *
 * \param file The file in which the violation occurred
 * \param line The line at which the violation occurred
 * \param expression The expression that caused the violation
 * \param message The message associated with the enforcement that
 *   experienced the violation
 * \param type The violation type (one of xContract_violation_type_t)
 */
XCONTRACT_NORETURN
XCONTRACT_CALL(void) xContract_contract_violation(
    char const*                     file
,   unsigned                        line
,   char const*                     expression
,   char const*                     message
,   xContract_violation_type_t      type
);

/** Invoked to report and respond to a contract violation
 *
 * \ingroup group__api_core
 *
 * \param file The file in which the violation occurred
 * \param line The line at which the violation occurred
 * \param function The function within which the violation occurred. Will be
 *   <code>NULL</code> for compilers that do not support the standard
 *   <code>__FUNCTION__</code> symbol
 * \param expression The expression that caused the violation
 * \param message The message associated with the enforcement that
 *   experienced the violation
 * \param qualifier A qualifier to the message. May be <code>NULL</code>
 * \param type The violation type (one of xContract_violation_type_t)
 * \param level The abstraction level of the enforcement
 * \param pfnReport Reporting function
 */
XCONTRACT_NORETURN
XCONTRACT_CALL(void) xContract_contract_violation2(
    char const*                     file
,   unsigned                        line
,   char const*                     function
,   char const*                     expression
,   char const*                     message
,   char const*                     qualifier
,   xContract_violation_type_t      type
,   int                             level
,   xContract_violationReport_fn_t  pfnReport
);

/** C-style string corresponding to the given violation type
 *
 * \ingroup group__api_core
 *
 * \param type The violation type
 */
XCONTRACT_CALL(char const*) xContract_getViolationTypeString(xContract_violation_type_t type);

/** Length of the C-style string corresponding to the given violation type
 *
 * \ingroup group__api_core
 *
 * \param type The violation type
 */
XCONTRACT_CALL(size_t) xContract_getViolationTypeStringLength(xContract_violation_type_t type);

/* /////////////////////////////////////////////////////////////////////////
 * String access shims
 */

#ifdef __cplusplus
inline char const* c_str_data_a(xContract_violation_type_t type)
{
	return xContract_getViolationTypeString(type);
}
inline size_t c_str_len_a(xContract_violation_type_t type)
{
	return xContract_getViolationTypeStringLength(type);
}
inline char const* c_str_ptr_a(xContract_violation_type_t type)
{
	return xContract_getViolationTypeString(type);
}
#endif /* __cplusplus */

/* /////////////////////////////////////////////////////////////////////////
 * Namespace
 */

#ifndef _XCONTRACT_NO_NAMESPACE
} /* namespace xcontract */
#endif /* !_XCONTRACT_NO_NAMESPACE */

#if !defined(_STLSOFT_NO_NAMESPACE) && \
    defined(__cplusplus)
namespace stlsoft
{

#  ifndef _XCONTRACT_NO_NAMESPACE
    using ::xcontract::c_str_data_a;
    using ::xcontract::c_str_len_a;
    using ::xcontract::c_str_ptr_a;
#  else /* ? !_XCONTRACT_NO_NAMESPACE */
    using ::c_str_data_a;
    using ::c_str_len_a;
    using ::c_str_ptr_a;
#  endif /* !_XCONTRACT_NO_NAMESPACE */

} /* namespace stlsoft */
#endif /* _STLSOFT_NO_NAMESPACE */

/* ////////////////////////////////////////////////////////////////////// */

#endif /* !XCONTRACT_INCL_XCONTRACT_H_XCONTRACT */

/* ///////////////////////////// end of file //////////////////////////// */
