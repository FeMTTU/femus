/*=========================================================================

 Program: FEMUS
 Module: Box
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_meshGencase_Box_hpp__
#define __femus_meshGencase_Box_hpp__

#include <vector>

#include "Domain.hpp"
#include "FemusInputParser.hpp"


namespace femus {


 class Box : public Domain {
 
 public:
    
   const uint GetDomainFlag() const { return 0; }
   
    std::vector<double> _lb;
    std::vector<double> _le;
 
    Box(const uint spacedim_in, FemusInputParser<double> & map_in);
   ~Box();
   
   void InitAndNondimensionalize(double Lref_in);

   void TransformPointToRef(const double* x_in,double* x_out) const; 
   
};

} //end namespace femus



#endif