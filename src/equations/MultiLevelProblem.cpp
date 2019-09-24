/*=========================================================================

 Program: FEMuS
 Module: MultiLevelProblem
 Authors: Eugenio Aulisa, Simone Bnà, Giorgio Bornia

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "ImplicitRungeKuttaSystem.hpp"
#include "TransientSystem.hpp"
#include "FemusConfig.hpp"
#include "Parameter.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "GeomElTypeEnum.hpp"
#include <iostream>

namespace femus {

using std::cout;
using std::endl;

bool (* Mesh::_SetRefinementFlag)(const std::vector < double >& x,
				  const int &ElemGroupNumber,const int &level) = NULL;


MultiLevelProblem::MultiLevelProblem( MultiLevelSolution *ml_sol):
				      _ml_sol(ml_sol),
				      _ml_msh(ml_sol->_mlMesh),
				      _gridn(_ml_msh->GetNumberOfLevels()),
				      _files(NULL)
{

}



MultiLevelProblem::~MultiLevelProblem(){
  for( system_iterator iterator = _systems.begin(); iterator != _systems.end(); iterator++) {
    delete iterator->second;
  }
}

System & MultiLevelProblem::add_system (const std::string& sys_type,
				      const std::string& name)
{
  // If the user already built a system with this name, we'll
  // trust them and we'll use it.  That way they can pre-add
  // non-standard derived system classes, and if their restart file
  // has some non-standard sys_type we won't throw an error.
  if (_systems.count(name))
  {
      return this->get_system(name);
  }
  // Build a basic System
  else if (sys_type == "Basic")
    this->add_system<System> (name);

  // Build a Newmark system
//   else if (sys_type == "Newmark")
//     this->add_system<NewmarkSystem> (name);

  // Build an Explicit system
  else if ((sys_type == "Explicit"))
    this->add_system<ExplicitSystem> (name);

  // Build an Implicit system
  else if ((sys_type == "Implicit") ||
	   (sys_type == "Steady"  ))
    this->add_system<ImplicitSystem> (name);

  // build a transient implicit linear system
  else if ((sys_type == "Transient") ||
	   (sys_type == "TransientImplicit") ||
	   (sys_type == "TransientLinearImplicit"))
    this->add_system<TransientLinearImplicitSystem> (name);

  // build a transient implicit nonlinear system
  else if (sys_type == "TransientNonlinearImplicit")
    this->add_system<TransientNonlinearImplicitSystem> (name);

  // build a transient explicit system
  else if (sys_type == "TransientExplicit")
    this->add_system<TransientExplicitSystem> (name);

  // build a linear implicit system
  else if (sys_type == "LinearImplicit")
    this->add_system<LinearImplicitSystem> (name);

  // build a nonlinear implicit system
  else if (sys_type == "NonlinearImplicit")
    this->add_system<NonLinearImplicitSystem> (name);

  else
  {
    std::cerr << "ERROR: Unknown system type: " << sys_type << std::endl;
  }

  // Return a reference to the new system
  return this->get_system(name);
}




void MultiLevelProblem::clear ()
{
  // Clear any additional parameters
  parameters.clear();

  // clear the systems.  We must delete them
  // since we newed them!
  while (!_systems.empty())
  {
    system_iterator pos = _systems.begin();

    System *sys = pos->second;
    delete sys;
    sys = NULL;

    _systems.erase (pos);
  }
}


// void MultiLevelProblem::init()
// {
//   const unsigned int n_sys = this->n_systems();
//
//   assert(n_sys != 0);
//
//   for (unsigned int i=0; i != this->n_systems(); ++i)
//     this->get_system(i).init();
// }


 void MultiLevelProblem::SetQruleAndElemType(const std::string quadr_order_in) {


  _qrule.reserve(_ml_msh->GetLevel(LEV_PICK)->GetDimension());
  for (int idim=0;idim < _ml_msh->GetLevel(LEV_PICK)->GetDimension(); idim++) {
          Gauss qrule_temp(_mesh->_geomelem_id[idim].c_str(),quadr_order_in.c_str());
         _qrule.push_back(qrule_temp);
           }

  // =======FEElems =====  //remember to delete the FE at the end
  const std::string  FEFamily[QL] = {"biquadratic","linear","constant"};
  _elem_type.resize(_ml_msh->GetLevel(LEV_PICK)->GetDimension());
  for (int idim=0;idim < _ml_msh->GetLevel(LEV_PICK)->GetDimension(); idim++)   _elem_type[idim].resize(QL);

  for (int idim=0;idim < _ml_msh->GetLevel(LEV_PICK)->GetDimension(); idim++) {
    for (int fe=0; fe<QL; fe++) {
       _elem_type[idim][fe] = _ml_msh->_finiteElement[ _mesh->_geomelem_flag[idim] ][ elem_type::_fe_old_to_new[fe] ];
     }
    }

   return;
}

} //end namespace femus


