/* /////////////////////////////////////////////////////////////////////////
 * File:        xcontract/implicit_link.h
 *
 * Purpose:     Implicit linking for the xContract API
 *
 * Created:     3rd March 2003
 * Updated:     1st April 2008
 *
 * Home:        http://xcontract.org/
 *
 * Copyright (c) 2003-2008, Matthew Wilson and Synesis Software
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


/** \file xcontract/implicit_link.h
 *
 * [C, C++] Implicit linking for the xContract library.
 *
 * This file may be <code>\#include</code>d in any compilation unit within a
 * given link unit to implicitly link in the appropriate <b>xContract</b>
 * archive (library).
 *
 * \note For compilers that do not support implicit linking, inclusion of
 *  the file has no effect.
 */

#ifndef XCONTRACT_INCL_XCONTRACT_H_IMPLICIT_LINK
#define XCONTRACT_INCL_XCONTRACT_H_IMPLICIT_LINK

/* /////////////////////////////////////////////////////////////////////////
 * Version information
 */

#ifndef XCONTRACT_DOCUMENTATION_SKIP_SECTION
# define XCONTRACT_VER_XCONTRACT_H_IMPLICIT_LINK_MAJOR      1
# define XCONTRACT_VER_XCONTRACT_H_IMPLICIT_LINK_MINOR      2
# define XCONTRACT_VER_XCONTRACT_H_IMPLICIT_LINK_REVISION   1
# define XCONTRACT_VER_XCONTRACT_H_IMPLICIT_LINK_EDIT       22
#endif /* !XCONTRACT_DOCUMENTATION_SKIP_SECTION */

/* /////////////////////////////////////////////////////////////////////////
 * Includes
 */

#ifndef XCONTRACT_INCL_XCONTRACT_H_XCONTRACT
# include <xcontract/xcontract.h>
#endif /* !XCONTRACT_INCL_XCONTRACT_H_XCONTRACT */

/* /////////////////////////////////////////////////////////////////////////
 * Implicit linking
 */

#if defined(WIN32) || \
    defined(_WIN32)

# if defined(__BORLANDC__) || \
     /* defined(__DMC__) || */ \
     defined(__INTEL_COMPILER) || \
     defined(__MWERKS__) || \
     defined(_MSC_VER)
#  define XCONTRACT_IMPLICIT_LINK_SUPPORT
# endif /* compiler */

# if defined(XCONTRACT_IMPLICIT_LINK_SUPPORT) && \
     defined(XCONTRACT_NO_IMPLICIT_LINK)
#  undef XCONTRACT_IMPLICIT_LINK_SUPPORT
# endif /* XCONTRACT_IMPLICIT_LINK_SUPPORT && XCONTRACT_NO_IMPLICIT_LINK */

# if defined(XCONTRACT_IMPLICIT_LINK_SUPPORT)

#  if defined(__BORLANDC__)
#   if __BORLANDC__ == 0x0550
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "bc55"
#   elif (__BORLANDC__ == 0x0551)
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "bc551"
#   elif (__BORLANDC__ == 0x0560)
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "bc56"
#   elif (__BORLANDC__ == 0x0564)
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "bc564"
#   elif (__BORLANDC__ == 0x0582)
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "bc582"
#   else /* ? __BORLANDC__ */
#    error Unrecognised value of __BORLANDC__
#   endif /* __BORLANDC__ */

/*
#  elif defined(__DMC__)
#   define XCONTRACT_IMPL_LINK_COMPILER_NAME                "dm"
 */

#  elif defined(__INTEL_COMPILER)
#   if __INTEL_COMPILER == 600
#    define   XCONTRACT_IMPL_LINK_COMPILER_NAME             "icl6"
#   elif __INTEL_COMPILER == 700
#    define   XCONTRACT_IMPL_LINK_COMPILER_NAME             "icl7"
#   elif __INTEL_COMPILER == 800
#    define   XCONTRACT_IMPL_LINK_COMPILER_NAME             "icl8"
#   elif __INTEL_COMPILER == 900
#    define   XCONTRACT_IMPL_LINK_COMPILER_NAME             "icl9"
#   else
#    error Intel C/C++ version not supported
#   endif /* _MSC_VER */

#  elif defined(__MWERKS__)
#   if ((__MWERKS__ & 0xFF00) == 0x2400)
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "cw7"
#   elif ((__MWERKS__ & 0xFF00) == 0x3000)
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "cw8"
#   elif ((__MWERKS__ & 0xFF00) == 0x3200)
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "cw9"
#   else /* ? __MWERKS__ */
#    error Unrecognised value of __MWERKS__
#   endif /* __MWERKS__ */

#  elif defined(_MSC_VER)
#   if _MSC_VER == 1000
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "vc4"
#   elif _MSC_VER == 1020
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "vc42"
#   elif _MSC_VER == 1100
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "vc5"
#   elif _MSC_VER == 1200
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "vc6"
#   elif _MSC_VER == 1300
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "vc7"
#   elif _MSC_VER == 1310
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "vc71"
#   elif _MSC_VER == 1400
#    define XCONTRACT_IMPL_LINK_COMPILER_NAME               "vc8"
#   else /* ? _MSC_VER */
#    error Visual C++ version not supported
#   endif /* _MSC_VER */

#  else /* ? compiler */
#   error Unrecognised compiler
#  endif /* compiler */


#  if defined(__MT__) || \
      defined(_REENTRANT) || \
      defined(_MT)
#   if defined(_DLL) || \
       defined(__DLL)
#    define XCONTRACT_IMPL_LINK_THREADING_TYPE              ".dll"
#   else /* ? dll */
#    define XCONTRACT_IMPL_LINK_THREADING_TYPE              ".mt"
#   endif /* dll */
#  else /* ? mt */
#    define XCONTRACT_IMPL_LINK_THREADING_TYPE              ""
#  endif /* mt */


#  if !defined(NDEBUG) && \
      defined(_DEBUG)
#   define XCONTRACT_IMPL_LINK_DEBUG_TYPE                   ".debug"
#  else /* ? _DEBUG */
#   define XCONTRACT_IMPL_LINK_DEBUG_TYPE                   ""
#  endif /* _DEBUG */

# define XCONTRACT_IMPL_LINK_LIB_PREFIX

# define XCONTRACT_IMPL_LINK_LIB_EXTENSION                  ".lib"


#  define XCONTRACT_IMPL_LINK_LIBRARY_BASENAME_s_(x)        #x
#  define XCONTRACT_IMPL_LINK_LIBRARY_BASENAME_s(x)         XCONTRACT_IMPL_LINK_LIBRARY_BASENAME_s_(x)
#  define XCONTRACT_IMPL_LINK_LIBRARY_BASENAME              "xcontract." XCONTRACT_IMPL_LINK_LIBRARY_BASENAME_s(_XCONTRACT_VER_MAJOR) ".core."

#  define XCONTRACT_IMPL_LINK_LIBRARY_NAME                  XCONTRACT_IMPL_LINK_LIB_PREFIX XCONTRACT_IMPL_LINK_LIBRARY_BASENAME XCONTRACT_IMPL_LINK_COMPILER_NAME XCONTRACT_IMPL_LINK_THREADING_TYPE XCONTRACT_IMPL_LINK_DEBUG_TYPE XCONTRACT_IMPL_LINK_LIB_EXTENSION

#  pragma message("lib: " XCONTRACT_IMPL_LINK_LIBRARY_NAME)

#  pragma comment(lib, XCONTRACT_IMPL_LINK_LIBRARY_NAME)

# endif /* XCONTRACT_IMPLICIT_LINK_SUPPORT */

#endif /* Win32 */

/* ////////////////////////////////////////////////////////////////////// */

#endif /* !XCONTRACT_INCL_XCONTRACT_H_IMPLICIT_LINK */

/* ////////////////////////////////////////////////////////////////////// */
