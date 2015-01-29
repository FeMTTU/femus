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

MultiLevelProblemTwo::MultiLevelProblemTwo(Files& files_in,
                           FemusInputParser<double> & phys_in,
                           QuantityMap& qtymap_in,
                           MultiLevelMeshTwo& mesh_in,
                           std::vector< std::vector<elem_type*> >  & elem_type_in,
			   std::vector<Gauss>   qrule_in ):
        _files(files_in),
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

// ====================================================
/// This sets dof initial and boundary conditions and sets the operators
//inside this functions there are a lot of new of class elements
//so, i must do a corresponding clean NOT in the DESTRUCTOR
//what is new'd in the constructor is deleted in the destructor
//what is new'd with this function is deleted with a corresponding clear function
//all these things are related to the Dofs, so the first thing is to settle them
// ===============================================================
/// This  function reads all the Operators from files
///  initialization of all levels: dofs and matrices;
/// initialize dof map
/// initialize BC
/// initialize MG Operators (TODO the values of the Restrictor Operators DEPEND on the boundary conditions)
/// initialize vectors (could do this even before the BC)
/// The INITIAL conditions can be done only after initVectors();

void  MultiLevelProblemTwo::setDofBcOpIc() {

    for (iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++) {
        SystemTwo* mgsol = eqn->second;
        
#ifdef DEFAULT_PRINT_INFO
    std::cout << "\n Reading "  <<  mgsol -> _eqname << " Dof, Bc, Op, Ic \n";
#endif

//=====================
    mgsol -> _dofmap.ComputeMeshToDof();
//=====================
    mgsol -> GenerateBdc();
    mgsol -> GenerateBdcElem();
//=====================
    mgsol -> ReadMGOps();
//=====================
    mgsol -> initVectors();     //TODO can I do it earlier than this position?
//=====================
    mgsol -> Initialize();              // initial solution

#ifdef DEFAULT_PRINT_INFO
    std::cout << " Dof, Bc, Op, Ic settled for"  <<  mgsol -> _eqname <<  "\n";
#endif
    }
    
    return;
}

} //end namespace femus
