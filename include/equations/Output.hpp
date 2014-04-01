/*=========================================================================

 Program: FEMUS
 Module: Output
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __output_h_
#define __output_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <vector>
#include <string>
//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class MultiLevelProblem;
class SparseMatrix;
class Vector;


class Output {

public:

  /** Constructor. */
  Output(MultiLevelProblem& ml_probl);

  /** Destructor */
  virtual ~Output();
  
  /** write output function */
  virtual void write_system_solutions() = 0;
  
  
protected:
  
  /** Build the prolongation matrices */ 
  void BuildProlongationMatrices();
  
  /** a set of matrix for the projection of the solution */ 
  static std::vector<SparseMatrix*> _ProlQitoQj_[3][3];
  
  /** a flag to move the output mesh */
  int _moving_mesh;
  
  /** the displacement variables for mesh moving */
  std::vector<std::string> _moving_vars;  
  
  /** the multilevelproblem reference */
  MultiLevelProblem& _ml_probl;

  
private:  
  
};
  
  
#endif