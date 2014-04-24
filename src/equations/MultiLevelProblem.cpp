/*=========================================================================

 Program: FEMUS
 Module: MultiLevelProblem
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelProblem.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "FEMTTUConfig.h"
#include "Parameter.hpp"
#include <iostream>

namespace femus {

using std::cout;
using std::endl;

bool (* mesh::_SetRefinementFlag)(const double &x, const double &y, const double &z, 
				  const int &ElemGroupNumber,const int &level) = NULL;

//---------------------------------------------------------------------------------------------------
MultiLevelProblem::MultiLevelProblem( MultiLevelMesh *ml_msh, MultiLevelSolution *ml_sol):
				      _gridn(ml_msh->GetNumberOfGrid()),
				      _gridr(ml_msh->GetNumberOfGridTotallyRefined()),
				      _ml_msh(ml_msh),
				      _ml_sol(ml_sol)
{
 
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


template <typename T_sys>
inline
const T_sys & MultiLevelProblem::get_system (const unsigned int num) const
{
  assert(num < this->n_systems());

  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    if (pos->second->number() == num)
      break;

  // Check for errors
  if (pos == end)
  {
    std::cerr << "ERROR: no system number " << num << " found!" << std::endl;
  }

  // Attempt dynamic cast
  return *static_cast<T_sys*>(pos->second);
}

template <typename T_sys>
inline
T_sys & MultiLevelProblem::get_system (const unsigned int num)
{
  assert(num < this->n_systems());

  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    if (pos->second->number() == num)
      break;

  // Check for errors
  if (pos == end)
  {
    std::cerr << "ERROR: no system number " << num << " found!" << std::endl;
  }

  // Attempt dynamic cast
  return *static_cast<T_sys*>(pos->second);
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

//---------------------------------------------------------------------------------------------------
int MultiLevelProblem::ComputeBdIntegral(const char pdename[],const char var_name[], const unsigned & kel, const unsigned & jface, unsigned level, unsigned dir) {
  
  
//   unsigned ipde=GetPdeIndex(pdename);
//   //   PetscVector* RESp=static_cast<PetscVector*> (_LinSolver[ipde][level]->_RES);  //TODO
//   //   Vec RES=RESp->vec(); //TODO
//   
    int ierr;
//   double tau;
//   double vx[3][27];
//   double phi[27],gradphi[27][3],Weight;
//   double normal[3];
//   int node[27];
//      
//   unsigned indexvar = GetSolPdeIndex(pdename,var_name);
//   short unsigned kelt = _msh[level]->el->GetElementType(kel);
//   unsigned order_ind = SolType[GetIndex(var_name)];
//   unsigned indX=GetIndex("X");
//   unsigned indY=GetIndex("Y");
//   unsigned indZ=GetIndex("Z");
// 		
//   bool test = _SetBoundaryConditionFunction(0.,0.,0.,var_name,tau,-(_msh[level]->el->GetFaceElementIndex(kel,jface)+1),_time);
//   if(!test) {
// 		
//     const short unsigned NV1[6][2]={{9,9},{6,6},{9,6},{3,3},{3,3},{1,1}};
//     unsigned nve=NV1[kelt][jface<_msh[level]->el->GetElementFaceNumber(kel,0)];
//     const unsigned FELT[6][2]={{3,3},{4,4},{3,4},{5,5},{5,5}};
//     unsigned felt = FELT[kelt][jface<_msh[level]->el->GetElementFaceNumber(kel,0)];
// 
//     //		cout << "felt : "  << felt << endl;
//     // 		cout << "node: " << nve << endl;
// 		
//     for(unsigned i=0;i<nve;i++) {
//       unsigned inode=_msh[level]->el->GetFaceVertexIndex(kel,jface,i)-1u;
//       node[i] = inode + _LinSolver[ipde][level]->KKIndex[indexvar];
//       unsigned inode_Metis=_msh[level]->GetMetisDof(inode,2);
//       
//       vx[0][i]=(*_solution[level]->_Sol[indX])(inode_Metis);  
//       vx[1][i]=(*_solution[level]->_Sol[indY])(inode_Metis);
//       vx[2][i]=(*_solution[level]->_Sol[indZ])(inode_Metis);
//       
//       
//     }
// 
//     for(unsigned igs=0;igs < type_elem[felt][order_ind]->GetGaussPointNumber(); igs++) {
//       (type_elem[felt][order_ind]->*type_elem[felt][order_ind]->Jacobian_sur_ptr)(vx,igs,Weight,phi,gradphi,normal);
//       //   		  cout << Weight << endl;
//       //  		  cout << "normal_x: " << normal[0] << endl;
//       //  		  cout << "normal_y: " << normal[1] << endl;
//       // 		  cout << "normal_z: " << normal[2] << endl;
//       //  		  cout << type_elem[felt][order_ind]->GetGaussPointNumber() << endl;
//       //  		  cout << "  " << endl;
//       for(unsigned i=0;i<nve;i++) {
// 	//     	          cout << phi[i] << "  " << nve << endl;
// 	_SetBoundaryConditionFunction(vx[0][i],vx[1][i],vx[2][i],var_name,tau,-(_msh[level]->el->GetFaceElementIndex(kel,jface)+1),_time);
//  	    
// 	//Manca la moltiplicazione x la normale che e'ancora da fare
// 	PetscScalar value = -phi[i]*tau*normal[dir]*Weight*_dt;
// 		    
// 	// Non voglio chiamare Vecsetvalue ma aggiungere il valore direttamente a F
// 	// per fare questo mi serve la relazione tra i(node locale di surface) e il nodo locale di volume
// 	_LinSolver[ipde][level]->_RES->add(node[i],value);
//       }
//     }
//   }
  return ierr;
}



} //end namespace femus


