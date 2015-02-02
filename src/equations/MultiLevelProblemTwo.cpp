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

#include "Files.hpp"
#include "XDMFWriter.hpp"
#include "XDMFWriter.hpp"
#include "Quantity.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "TimeLoop.hpp"

#include "paral.hpp"


namespace femus {



// ====================================================
/// This function constructs the equation map

MultiLevelProblemTwo::MultiLevelProblemTwo(const FemusInputParser<double> & phys_in,
                           const QuantityMap& qtymap_in,
                           const MultiLevelMeshTwo& mesh_in,
                           const std::vector< std::vector<elem_type*> >  & elem_type_in,
			   const std::vector<Gauss>   qrule_in ):
        _phys(phys_in),
        _qtymap(qtymap_in),
        _mesh(mesh_in),
        _elem_type(elem_type_in),
        _qrule(qrule_in),
        MultiLevelProblem()  { }


// ====================================================
/// This function destroys the equations
void MultiLevelProblemTwo::clean() {
    for (MultiLevelProblemTwo::iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++) {
        delete eqn->second;
    }
}


} //end namespace femus
