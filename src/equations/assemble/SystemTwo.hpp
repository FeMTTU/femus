/*=========================================================================

 Program: FEMUS
 Module: SystemTwo
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_SystemTwo_hpp__
#define __femus_equations_SystemTwo_hpp__

//C++ includes
#include <vector>
#include <string>
#include <map>

//FEMuS includes 
#include "FemusDefault.hpp"
#include "FemusConfig.hpp"
#include "Typedefs.hpp"
#include "FETypeEnum.hpp"
#include "VBTypeEnum.hpp"
#include "DofMap.hpp"
#include "BoundaryConditions.hpp"
#include "Quantity.hpp"

#include "NonLinearImplicitSystem.hpp"


namespace femus {



class MultiLevelProblem;
class MultiLevelMeshTwo;

class SparseMatrix;
class NumericVector;
class LinearEquationSolver;


class SystemTwo : public NonLinearImplicitSystem {

public:

//=======================================================================
// CONSTRUCTOR / DESTRUCTOR
//=======================================================================
  SystemTwo(MultiLevelProblem & equations_map, const std::string & eq_name_in, const unsigned int number, const LinearEquationSolverType & smoother_type);   //System//
  
  DofMap  _dofmap;  //// LinearEquation (each level)
  
  BoundaryConditions _bcond;
  
  void  initVectors();  ///initialize vectors       //System//

//=======================================================================
//======= Quantities =========
//=======================================================================
      inline const std::vector<Quantity*> & GetUnknownQuantitiesVector() const { 	return _UnknownQuantitiesVector;   }//MultilevelSolution//
      
     
      void AddUnknownToSystemPDE( Quantity* qty_in) {                  //System//
	  unsigned n = _UnknownQuantitiesVector.size();

	_UnknownQuantitiesVector.resize(n+1);
	_UnknownQuantitiesVector[n] = qty_in;
	qty_in->set_eqn(this);
	qty_in->SetPosInAssocEqn(n);
	
	return;
      }
      

	    std::vector<std::string> _var_names;                   //MultilevelSolution//
	  void initVarNames();                                     //MultilevelSolution//

	  std::vector<double>      _refvalue;        //MultilevelSolution//
          void initRefValues();                      //MultilevelSolution//

	  void init_unknown_vars();     //System//  //MultilevelSolution//

// ============ INITIAL CONDITIONS of the equation ====== (procs,levels) ==    
          void    Initialize();           //MultilevelSolution  //this uses x and fills in x_old at all levels
          
protected:
  
  std::vector<Quantity*>          _UnknownQuantitiesVector;  //MultilevelSolution//

};



} //end namespace femus



#endif