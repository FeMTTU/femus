#ifndef __domain_h__
#define __domain_h__

#include <string>
#include "Typedefs.hpp"
#include "RunTimeMap.hpp"


namespace femus {



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

    Domain(const uint spacedim_in, const RunTimeMap<double> & map_in);
 
 public: 

   ~Domain();
   
   const RunTimeMap<double> _domain_rtmap;

   std::string _name;
   double      _Lref;
   uint    _spacedim;
   ///transformation to reference Domain frame
   virtual void TransformPointToRef(const double* x_in, double* x_out) const = 0;
   
  };



} //end namespace femus



#endif