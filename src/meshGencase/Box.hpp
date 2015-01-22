#ifndef __box_h__
#define __box_h__


#include "Domain.hpp"
#include "FemusInputParser.hpp"


namespace femus {




//------------ BOX -----------  
//------------ BOX ----------- 
//------------ BOX ----------- 
 
 class Box : public Domain {
 
 public:
    
   const uint GetDomainFlag() const { return 0;}
   
    double* _lb;
    double* _le;
 
   Box(const uint spacedim_in, FemusInputParser<double> & map_in);
   ~Box();
   
   void InitAndNondimensionalize(double Lref_in);

   void TransformPointToRef(const double* x_in,double* x_out) const; 
   
};

//if a function has a const and it is the correspondent of a virtual function that does not have const,
//THE TYPES are DIFFERENT!! The const related to the class data is an attribute that changes the things!
// Here, the function ElFlagControl is not in the correct place.
//The idea in fact is that things related to Control should stay in a separate class
//dedicated to the Control.
//Clearly this function is dependent upon the domain, 
//like the BCIC is dependent upon the domain.
//But, i do not want to have an equation dependent on the domain,
//because the equation is just a bunch of operators
//therefore, i could define the control flag as another quantity,
//associated to the Zero order finite elements
//But then the equation is still specialized by another Quantity.

//So the equation is dependent on the Domain in any way, either 
//because of the Physical Quantities or because of the Control Flag

//Now, let us consider using the Specific Domains for Gencase.
//It may be useful, we could tune them with the libmesh calls
//now, the point is that in Gencase we DO NOT HAVE to nondimensionalize the domain!
//So, we still need a box but without Lref! This is good, because otherwise
//we should read the Lref parameter which is now associated to the Physics.
//So, we'll do another constructor that has no need for Lref.
//Then, since the Domain is needed, we'll have to do another constructor
//for the Domain




} //end namespace femus



#endif