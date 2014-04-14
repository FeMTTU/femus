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
#include "ParallelObject.hpp"


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class MultiLevelSolution;
class SparseMatrix;
class Vector;


class Output : public ParallelObject {

public:

  /** Constructor. */
  Output(MultiLevelSolution& ml_probl);

  /** Destructor */
  virtual ~Output();
  
  /** write output function */
  virtual void write_system_solutions(const char order[], std::vector<std::string>& vars, const unsigned time_step = 0) = 0;
  
  /** set moving mesh */
  void SetMovingMesh(std::vector<std::string>& movvars_in);
  
protected:
  
  /** a set of matrix for the projection of the solution */ 
  static std::vector<SparseMatrix*> _ProlQitoQj[3][3];
  
  /** a flag to move the output mesh */
  int _moving_mesh;
  
  /** the displacement variables for mesh moving */
  std::vector<std::string> _moving_vars;  
  
  /** the multilevelproblem reference */
  MultiLevelSolution& _ml_sol;
  
  int _gridn;
  
  int _gridr;
  
  
private:  
  
  /** This routine generates the matrices for the projection of the solutions from different FE spaces */
  void BuildProlongatorMatrices();
  
};
  


} //end namespace femus


  
#endif