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

#ifndef __mgsolbase_hpp_
#define __mgsolbase_hpp_

//C++ includes
#include <vector>
#include <string>
#include <map>

//FEMuS includes 
#include "FemusDefault.hpp"
#include "FEMTTUConfig.h"
#include "Typedefs.hpp"
#include "FETypeEnum.hpp"
#include "VBTypeEnum.hpp"
#include "DofMap.hpp"
#include "BoundaryConditions.hpp"
#include "Quantity.hpp"

#include "LinearImplicitSystem.hpp"


namespace femus {



class MultiLevelProblem;
class MultiLevelMeshTwo;

class SparseMatrix;
class NumericVector;
class LinearEquationSolver;


class SystemTwo : public LinearImplicitSystem {

public:

//=======================================================================
// CONSTRUCTOR / DESTRUCTOR
//=======================================================================
  SystemTwo(MultiLevelProblem & equations_map, const std::string & eq_name_in, const unsigned int number, const MgSmoother & smoother_type);   //System//
  
  virtual ~SystemTwo();                    //System//

  
  DofMap  _dofmap;  //// LinearEquation (each level)
  
  BoundaryConditions _bcond;
  
//=======================================================================
//======== MG Ops ============ (procs,levels) ====
//=======================================================================
  std::vector<SparseMatrix *> _A;  // LinearEquation (each level)
  std::vector<SparseMatrix *> _Rst; // LinearEquation (each level)
  std::vector<SparseMatrix *> _Prl; // LinearEquation (each level)
  
    void ReadMGOps(const std::string output_path); // LinearEquation  (each level)
    void ReadMatrix(const std::string& name); // LinearEquation  (each level)
    void ReadProl(const std::string& name);   // LinearEquation  (each level)
    void ReadRest(const std::string& name);   // LinearEquation  (each level)
    void ComputeMatrix();                     // LinearEquation  (each level)
    void ComputeProl();                       // LinearEquation  (each level)
    void ComputeRest();                       // LinearEquation  (each level)
  
//=======================================================================
//======== Vectors =============== (procs,levels) ==
//=======================================================================

  std::vector<NumericVector *> _b;   //// LinearEquation (each level)
  std::vector<NumericVector *> _x;   //// LinearEquation (each level)
  std::vector<NumericVector *> _res; //// LinearEquation (each level)

  std::vector<NumericVector *> _x_old; //// LinearEquation (each level)
  
  std::vector<NumericVector *> _x_oold;    //this is used by MGTimeStep and also by the OptLoop
  std::vector<NumericVector *> _x_tmp;     //this is used by MGTimeStep and also by the OptLoop
  
          void  initVectors();  ///initialize vectors       //System//

//===================
  void MGSolve(double Eps,int MaxIter, const uint Gamma=DEFAULT_MG_GAMMA, const uint Nc_pre=DEFAULT_NC_PRE,const uint Nc_coarse=DEFAULT_NC_COARSE,const uint Nc_post=DEFAULT_NC_POST);  //LinearImplicitSystem//
  double MGStep(int Level,double Eps1,int MaxIter, const uint Gamma, const uint Nc_pre,const uint Nc_coarse,const uint Nc_post);                                                          //LinearImplicitSystem//

//=======================================================================
//======= Quantities =========
//=======================================================================
      inline const std::vector<Quantity*> & GetUnknownQuantitiesVector() const { //MultilevelSolution//
	return _UnknownQuantitiesVector;
      }
      
     
      void AddUnknownToSystemPDE( Quantity* qty_in) { //System//
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

	  void init_sys();     //System//  //MultilevelSolution//

// ============ INITIAL CONDITIONS of the equation ====== (procs,levels) ==    
          void    Initialize();           //MultilevelSolution  //this uses x and fills in x_old at all levels
          
protected:
  
  std::vector<Quantity*>          _UnknownQuantitiesVector;  //MultilevelSolution//

};



} //end namespace femus



#endif