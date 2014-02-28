#ifndef __domain_h__
#define __domain_h__

#include <string>
#include "Typedefs_conf.hpp"

class Utils;

//----------------------------- 
//------------ BASE DOMAIN class -----------  
//------------------------------ 

//if we make the constructor protected, then this class cannot be directly instantiated
// in fact this is only a Base class
//Base for Equations: EqnBase
//Base for Quantities: Quantity
//Base for Domains: Domain

 class Domain {


 protected:

    Domain(Utils& utils_in);
 
 public: 

   ~Domain();
   
   Utils &    _utils;
   std::string _name;
   double      _Lref;
   uint    _spacedim;
   ///transformation to reference Domain frame
   virtual void TransformPointToRef(const double* x_in, double* x_out) const = 0;
   
  };


#endif