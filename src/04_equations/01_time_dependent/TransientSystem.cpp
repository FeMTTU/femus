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
#include "MultiLevelSolution.hpp"
#include "TransientSystem.hpp"
#include "NumericVector.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"

namespace femus {

// ------------------------------------------------------------
// TransientSystem implementation
template <class Base>
TransientSystem<Base>::TransientSystem (
                    MultiLevelProblem& ml_probl,
					const std::string& name_in,
					const unsigned int number_in,
                    const LinearEquationSolverType & smoother_type) :
  Base (ml_probl, name_in, number_in,smoother_type),
  _is_selective_timestep(false),
  _time(0.),
  _time_step(0),
  _dt(0.1),
  _assembleCounter(0)
{

}

// ------------------------------------------------------------
template <class Base>
TransientSystem<Base>::~TransientSystem ()
{
   // clear the parent data
  _is_selective_timestep = false;
  _time = 0.;
  _time_step = 0;
  _dt = 0.1;
  _assembleCounter= 0;   
}

// ------------------------------------------------------------
template <class Base>
void TransientSystem<Base>::CopySolutionToOldSolution() {

  for (int ig=0; ig< this->_gridn; ig++) {
    this->_solution[ig]->CopySolutionToOldSolution();
  }

}

// ------------------------------------------------------------
template <class Base>
void TransientSystem<Base>::SetUpForSolve(){
  double dtOld = _dt;

  if (_is_selective_timestep) {
    _dt = _get_time_interval_function(_time);
  }

  if(_assembleCounter % 1 == 0 || _dt != dtOld){
    std::cout<<"Assemble Matrix\n";
    Base::_buildSolver = true;
    if( _dt != dtOld )
      _assembleCounter = 0;
  }
  else{
    std::cout<<"Do not Assemble Matrix";
    Base::_buildSolver = false;
  }
  std::cout<<"assemble counter = "<<_assembleCounter<<std::endl;
  _assembleCounter++;


  //update time
  _time += _dt;

  //update time step
  _time_step++;

  std::cout << " Simulation Time: " << _time << " TimeStep: " << _dt << " Iteration: " << _time_step << std::endl;

   //update boundary condition
  this->_ml_sol->UpdateBdc(_time); 
}

// ------------------------------------------------------------
template <class Base>
void TransientSystem<Base>::MGsolve( const MgSmootherType& mgSmootherType ) {

  SetUpForSolve();  

  Base::MGsolve( mgSmootherType );

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

// ------------------------------------------------------------
template <class Base>
void TransientSystem<Base>::AttachGetTimeIntervalFunction (double (* get_time_interval_function)(const double time)) {
  _get_time_interval_function = get_time_interval_function;
  _is_selective_timestep = true;
}

// ------------------------------------------------------------
// TransientSystem forward instantiations
template class TransientSystem<LinearImplicitSystem>;
template class TransientSystem<NonLinearImplicitSystem>;
template class TransientSystem<MonolithicFSINonLinearImplicitSystem>;
template class TransientSystem<ExplicitSystem>;
template class TransientSystem<System>;

} //end namespace femus


