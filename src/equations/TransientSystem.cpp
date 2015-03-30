/*=========================================================================

 Program: FEMUS
 Module: TransientSystem
 Authors: Simone Bnà
 
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
#include "TransientSystem.hpp"
#include "ExplicitSystem.hpp"
#include "LinearImplicitSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "NumericVector.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"

namespace femus {



// ------------------------------------------------------------
// TransientSystem implementation
template <class Base>
TransientSystem<Base>::TransientSystem (MultiLevelProblem& ml_probl,
					const std::string& name_in,
					const unsigned int number_in,const MgSmoother & smoother_type) :

  Base (ml_probl, name_in, number_in,smoother_type),
  _is_selective_timestep(false),
  _time(0.),
  _time_step(0),
  _dt(0.1)
{

}

template <class Base>
TransientSystem<Base>::~TransientSystem ()
{
  this->clear();
}

template <class Base>
void TransientSystem<Base>::clear ()
{
  // clear the parent data
  Base::clear(); 

}

template <class Base>
void TransientSystem<Base>::UpdateSolution() {
 
  for (int ig=0; ig< this->_gridn; ig++) {   
    this->_solution[ig]->UpdateSolution();
  }
  
}

template <class Base>
void TransientSystem<Base>::solve() {
 
  if (_is_selective_timestep) {
    _dt = _get_time_interval_function(_time);
  }

  //update time
  _time += _dt;
   
  //update time step
  _time_step++;
  
  std::cout << " Time: " << _time << "   TimeStep: " << _time_step << std::endl;
    
   //update boundary condition
  this->_ml_sol->UpdateBdc(_time);
  
  // call the parent solver
  Base::solve();
  
}


//---------------------------------------------------------------------------------------------------------

template <class Base>
void TransientSystem<Base>::NewmarkAccUpdate() {

  const double gamma = 0.5;
  const double a5    = -1.*(1. - gamma)/gamma;
  const double a1    = 1./(gamma*_dt);
  const double a2    = -1./(gamma*_dt);
  
  unsigned dim = this->_msh[0]->GetDimension();
 
  unsigned axyz[3];
  unsigned vxyz[3];
  const char accname[3][3] = {"AX","AY","AZ"};
  const char velname[3][2] = {"U","V","W"};
   
  for(unsigned i=0; i<dim; i++) {
     axyz[i] = this->_ml_sol->GetIndex(&accname[i][0]);
     vxyz[i] = this->_ml_sol->GetIndex(&velname[i][0]);
  }
   
  for (int ig=0;ig< this->_gridn;ig++) {
    for(unsigned i=0; i<dim; i++) {
      this->_solution[ig]->_Sol[axyz[i]]->scale(a5);
      this->_solution[ig]->_Sol[axyz[i]]->add(a1,*(this->_solution[ig]->_Sol[vxyz[i]]));
      this->_solution[ig]->_Sol[axyz[i]]->add(a2,*(this->_solution[ig]->_SolOld[vxyz[i]]));
    }
  }
}


//------------------------------------------------------------------------------------------------------

// template <class Base>
// int TransientSystem<Base>::SaveData() const {
//   char *filename = new char[80];
//   PetscVector* petsc_vec_sol;
// 
//   PetscErrorCode ierr;
//   PetscViewer bin_viewer;
//   for (unsigned ig=0; ig<_gridn; ig++) {
//     for (unsigned i=0; i<SolType.size(); i++) {
//       sprintf(filename,"./save/save.time_%d.level_%d.%s.bin",_time_step,ig,SolName[i]);
//       ierr = PetscViewerBinaryOpen(MPI_COMM_WORLD,filename,FILE_MODE_WRITE,&bin_viewer);
//       CHKERRQ(ierr);
//       //petsc_vec_sol  = static_cast<PetscVector*>(_LinSolver[0][ig]->Sol_[i]);
//       petsc_vec_sol  = static_cast<PetscVector*>(_solution[ig]->_Sol[i]);
//       ierr = VecView(petsc_vec_sol->vec(),bin_viewer);
//       CHKERRQ(ierr);
//       ierr = PetscViewerDestroy(&bin_viewer);
//     }
//   }
//   delete [] filename;
//   return ierr;
// }
// 
// //------------------------------------------------------------------------------------------------------

// template <class Base>
// int TransientSystem<Base>::InitializeFromRestart(unsigned restart_time_step) {
// 
//   // Set the restart time step
//   SetInitTimeStep(restart_time_step+1);
//   
//   //Set the restart time
//   if(_ats_flag==0) {
//     _time = restart_time_step*_dt;
//   } 
//   else {
//     _time = 0;
//     for(unsigned i=0; i<restart_time_step; i++) {
//       _dt = _set_time_step_function(_time);
//       _time += _dt;
//     }
//   }
//    
//   //Set the restart time step
//   _time_step = restart_time_step+1;
//      
//   char *filename = new char[80];
//   PetscVector* petsc_vec_sol;
// 
//   PetscErrorCode ierr;
//   PetscViewer bin_viewer;
//   for(unsigned ig=0;ig<_gridn;ig++){
//     for(unsigned i=0;i<SolType.size();i++){
//       sprintf(filename,"./save/save.time_%d.level_%d.%s.bin",restart_time_step,ig,SolName[i]);
//       ierr = PetscViewerBinaryOpen(MPI_COMM_WORLD,filename,FILE_MODE_READ,&bin_viewer); CHKERRQ(ierr);
//       //petsc_vec_sol  = static_cast<PetscVector*>(_LinSolver[0][ig]->Sol_[i]);
//       petsc_vec_sol  = static_cast<PetscVector*>(_solution[ig]->_Sol[i]);
//       ierr = VecLoad(petsc_vec_sol->vec(),bin_viewer); CHKERRQ(ierr);
//       ierr = PetscViewerDestroy(&bin_viewer);
//     }
//     _solution[ig]->UpdateSolution();
//   }
//   
//   delete [] filename;
//   return ierr;
// }



template <class Base>
void TransientSystem<Base>::AttachGetTimeIntervalFunction (double (* get_time_interval_function)(const double time)) {
  _get_time_interval_function = get_time_interval_function;
  _is_selective_timestep = true;
}


// ------------------------------------------------------------
// TransientSystem instantiations
template class TransientSystem<LinearImplicitSystem>;
template class TransientSystem<NonLinearImplicitSystem>;
template class TransientSystem<MonolithicFSINonLinearImplicitSystem>;
template class TransientSystem<ExplicitSystem>;
template class TransientSystem<System>;



} //end namespace femus


