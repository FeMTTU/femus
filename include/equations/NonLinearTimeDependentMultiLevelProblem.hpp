/*=========================================================================

 Program: FEMUS
 Module: NonLinearTimeDependentMultiLevelProblem
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __NonLinearTimeDependentMultiLevelProblem_hpp__
#define __NonLinearTimeDependentMultiLevelProblem_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelProblem.hpp"

/**
* This class is a black box container to handle time-dependent multilevel problems
* This class have to be deleted as soon as possible
*/

class NonLinearTimeDependentMultiLevelProblem : public MultiLevelProblem {

private:
  unsigned _save_step;

  bool _ats_flag;      //adaptive time step flag

  //pointer function to the set time step function
  double (* _set_time_step_function)(const double time);

  void SetInitTimeStep(const unsigned time_step0);
  
  

public:

  ///Constructor
  NonLinearTimeDependentMultiLevelProblem(const unsigned short &igridn, const unsigned short &igridr,
                    const char mesh_file[], const char GaussOrder[], const double Lref=1., bool (* SetRefinementFlag)(const double &x, const double &y, const double &z, 
                     const int &ElemGroupNumber,const int &level)=NULL);

  ///Destructor
  ~NonLinearTimeDependentMultiLevelProblem();
 
  void SetTimeStep(const double dt);

  void SetPrintTimeStep(const unsigned print_step);

  void SetSaveTimeStep(const unsigned save_step);
  
  unsigned GetInitTimeStep() const;

  unsigned GetNumTimeSteps() const;
  
  unsigned GetPrintTimeStep() const;

  unsigned GetSaveTimeStep() const;

  double GetTimeStep();

  double GetTime();
  
  int SaveData() const;
  
  void _NewmarkAccUpdate();

  void _UpdateSolution();
    
  void AttachSetTimeStepFunction (double (* set_time_step_function)(const double time));

  int InitializeFromRestart(unsigned restart_time_step);

  void SetNumTimeSteps(const double ntimesteps);

  void Solve(const char pdename[], unsigned const &ncycle,  unsigned const &npost, unsigned const &npre, 
		    const char mg_type[]="F-Cycle", const bool &test_linear=0);
  


};

#endif