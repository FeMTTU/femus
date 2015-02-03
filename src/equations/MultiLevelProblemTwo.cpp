/*=========================================================================

 Program: FEMUS
 Module: MultiLevelProblemTwo
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "MultiLevelProblemTwo.hpp"

// C++
#include <iomanip>
#include <sstream>

// local includes 
#include "FemusDefault.hpp"
#include "XDMFWriter.hpp"
#include "Quantity.hpp"
#include "MultiLevelMeshTwo.hpp"

#include "paral.hpp"


namespace femus {



// ====================================================
/// This function constructs the equation map

MultiLevelProblemTwo::MultiLevelProblemTwo(const MultiLevelMeshTwo& mesh_in,
			   const std::string quadr_order_in):
        _mesh(mesh_in),
        MultiLevelProblem()  {
	  
// ======  QRule ================================
	  
  _qrule.reserve(_mesh.get_dim());
  for (int idim=0;idim < _mesh.get_dim(); idim++) { 
          Gauss qrule_temp(_mesh._geomelem_id[idim].c_str(),quadr_order_in.c_str());
         _qrule.push_back(qrule_temp);
           }  

  // =======FEElems =====  //remember to delete the FE at the end
  const std::string  FEFamily[QL] = {"biquadratic","linear","constant"}; 
  _elem_type.resize(_mesh.get_dim());
  for (int idim=0;idim < _mesh.get_dim(); idim++)   _elem_type[idim].resize(QL);
  
  for (int idim=0;idim < _mesh.get_dim(); idim++) { 
    for (int fe=0; fe<QL; fe++) {
       _elem_type[idim][fe] = elem_type::build(_mesh._geomelem_id[idim].c_str(),fe,quadr_order_in.c_str());
       _elem_type[idim][fe]->EvaluateShapeAtQP(_mesh._geomelem_id[idim].c_str(),fe);
     }
    }  
           
	  
  }


// ====================================================
/// This function destroys the equations
void MultiLevelProblemTwo::clean() {
    for (MultiLevelProblemTwo::iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++) {
        delete eqn->second;
    }
}


} //end namespace femus
