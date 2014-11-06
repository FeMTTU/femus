// ad_library_switch.h - Enable advection_schemes.cpp to be compiled multiple times with different AD libraries

#define preallocate_statements(n)
#define preallocate_operations(n)

#ifdef ADEPT
#undef preallocate_statements
#undef preallocate_operations

#include "adept.h"
using namespace adept;
#include "advection_schemes.h"


/* *************************************** */
#else
#ifdef CPPAD 

#include "cppad/cppad.hpp"
using CppAD::AD;
#define adouble AD<double>
#include "advection_schemes.h"

template<class A>
static
inline
A value(const AD<A>& x) { return CppAD::Value(Var2Par(x)); }


/* *************************************** */
#else
#ifdef ADOLC

#include "adolc/adouble.h"
#define abs fabs
static
inline
double value(const adouble& x) { return x.getValue(); }
#include "advection_schemes.h"

/* *************************************** */
#else
#ifdef SACADO

#include "Sacado.hpp"
#define adouble Sacado::Rad::ADvar<double>

static
inline
double value(const adouble& x) { return x.val(); }

#include "advection_schemes.h"

/* *************************************** */
#else
#ifdef SACADO_FAD

#include "Sacado.hpp"
#define adouble Sacado::ELRFad::DFad<double>

static
inline
double value(const adouble& x) { return x.val(); }

#include "advection_schemes.h"

/* *************************************** */
#else

#define adouble double

static
inline
const double& value(const double& val) { return val; }
#include "advection_schemes.h"
#define abs fabs

#endif
#endif
#endif
#endif
#endif
