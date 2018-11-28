/*=========================================================================

 Program: FEMUS
 Module: ImplicitRungeKutta
 Authors: Kara Erdi
 
 /*=========================================================================

 Program: FEMUS
 Module: ImplicitRungeKuttaSystem
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
#include "ImplicitRungeKuttaSystem.hpp"
#include "ExplicitSystem.hpp"
#include "LinearImplicitSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "NumericVector.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"

namespace femus {



// ------------------------------------------------------------
// ImplicitRungeKuttaSystem implementation
template <class Base>
ImplicitRungeKuttaSystem<Base>::ImplicitRungeKuttaSystem (MultiLevelProblem& ml_probl,
					const std::string& name_in,
					const unsigned int number_in,const MgSmoother & smoother_type) :

  Base (ml_probl, name_in, number_in,smoother_type),
  _is_selective_timestep(false),
  _time(0.),
  _time_step(0),
  _dt(0.1),
  _assembleCounter(0),
  _RK(1u)
{

}

template <class Base>
ImplicitRungeKuttaSystem<Base>::~ImplicitRungeKuttaSystem ()
{
  this->clear();
}

template <class Base>
void ImplicitRungeKuttaSystem<Base>::clear ()
{
  // clear the parent data
  Base::clear();

}

template <class Base>
void ImplicitRungeKuttaSystem<Base>::AddSolutionToSystemPDE(const char solname[]){
    
  for(unsigned i = 0; i < _RK; i++){
    std::ostringstream solnameki;
    solnameki << solname << "k" << i+1;
    this->Base::AddSolutionToSystemPDE(solnameki.str().c_str());
  }
  
}


template <class Base>
void ImplicitRungeKuttaSystem<Base>::CopySolutionToOldSolution() {

  for (int ig=0; ig< this->_gridn; ig++) {
    this->_solution[ig]->CopySolutionToOldSolution();
  }

}

template <class Base>
void ImplicitRungeKuttaSystem<Base>::MLsolve() {

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

  std::cout << " Simulation Time: " << _time << "   TimeStep: " << _time_step << std::endl;

   //update boundary condition
  this->_ml_sol->UpdateBdc(_time);

  // call the parent solver
  Base::_MLsolver = true;
  Base::_MGsolver = false;

  Base::solve();

}

template <class Base>
void ImplicitRungeKuttaSystem<Base>::MGsolve( const MgSmootherType& mgSmootherType ) {

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

  std::cout << " Simulation Time:  " << _time << "   TimeStep: " << _time_step << std::endl;

   //update boundary condition
  this->_ml_sol->UpdateBdc(_time);

  // call the parent solver
  Base::_MLsolver = false;
  Base::_MGsolver = true;

  Base::solve( mgSmootherType );

}


template <class Base>
void ImplicitRungeKuttaSystem<Base>::AttachGetTimeIntervalFunction (double (* get_time_interval_function)(const double time)) {
  _get_time_interval_function = get_time_interval_function;
  _is_selective_timestep = true;
}


// ------------------------------------------------------------
// ImplicitRungeKuttaSystem instantiations
template class ImplicitRungeKuttaSystem<LinearImplicitSystem>;
template class ImplicitRungeKuttaSystem<NonLinearImplicitSystem>;
template class ImplicitRungeKuttaSystem<MonolithicFSINonLinearImplicitSystem>;
template class ImplicitRungeKuttaSystem<ExplicitSystem>;
template class ImplicitRungeKuttaSystem<System>;



} //end namespace femus


