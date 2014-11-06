#ifdef ADEPT

#include "adept.h"
using namespace adept;
#define areal aReal
#include "ms_ad.h"

/* *************************************** */
#else
#ifdef CPPAD 

#include "cppad/cppad.hpp"
using CppAD::AD;
#define areal AD<ms_real>
#include "ms_ad.h"

template<class A>
static
inline
A value(const AD<A>& x) { return CppAD::Value(Var2Par(x)); }


/* *************************************** */
#else
#ifdef ADOLC

#include "adolc/adouble.h"
//#include "adolc/interfaces.h"
//#include "adolc/taping.h"
#define areal adouble
#define abs fabs
static
inline
double value(const adouble& x) { return x.getValue(); }
#include "ms_ad.h"

/* *************************************** */
#else
#ifdef SACADO

#include "Sacado.hpp"
#define areal Sacado::Rad::ADvar<double>

static
inline
double value(const areal& x) { return x.val(); }

#include "ms_ad.h"
/* *************************************** */
#else
#ifdef SACADO_FAD

#include "Sacado.hpp"
#define areal Sacado::ELRFad::DFad<double>

static
inline
double value(const areal& x) { return x.val(); }

#include "ms_ad.h"

/* *************************************** */
#else

#define areal ms_real

static
inline
const ms_real& value(const ms_real& val) { return val; }

#endif
#endif
#endif
#endif
#endif
