#ifndef __domain_h__
#define __domain_h__

#include <string>
#include "Typedefs.hpp"
#include "FemusInputParser.hpp"


namespace femus {



//----------------------------- 
//------------ BASE DOMAIN class -----------  
//------------------------------ 

//if we make the constructor protected, then this class cannot be directly instantiated
// in fact this is only a Base class
//Base for Equations: SystemTwo
//Base for Quantities: Quantity
//Base for Domains: Domain

 class Domain {


 protected:

    Domain(const uint spacedim_in, const FemusInputParser<double> & map_in);
 
 public: 

   ~Domain();
   
   const FemusInputParser<double> & _domain_rtmap;  //TODO if you dont put the reference it crashes in debug mode... I think it is a const-related thing

   std::string _name;
   double      _Lref;
   uint    _spacedim;
   ///transformation to reference Domain frame
   virtual void TransformPointToRef(const double* x_in, double* x_out) const = 0;
   
   virtual const uint GetDomainFlag()  const = 0;
   
  };



} //end namespace femus



#endif