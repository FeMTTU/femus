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

#include "System.hpp"


namespace femus {



class Physics      ;
class MultiLevelProblem;
class MultiLevelMeshTwo      ;
class FEElemBase   ;
class Quantity     ;
class CurrentQuantity;

class SparseMatrix;
class NumericVector;
class LinearEquationSolver;



class SystemTwo : public System {

public:

//=======================================================================
// CONSTRUCTOR / DESTRUCTOR
//=======================================================================
  SystemTwo(MultiLevelProblem & equations_map, const std::string & eq_name_in, const unsigned int number, const MgSmoother & smoother_type);   //System//
  
  virtual ~SystemTwo();                    //System//
  
//=======================================================================
//======== MG Ops ============ (procs,levels) ====
//=======================================================================
  std::vector<SparseMatrix  *> _A;  // LinearEquation (each level)
  std::vector<SparseMatrix *> _Rst; // LinearEquation (each level)
  std::vector<SparseMatrix *> _Prl; // LinearEquation (each level)
  
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

//======= Linear Solvers for every Level ============
  LinearEquationSolver **_solver;     ///(each level)

//=======================================================================
//======= Quantities =========
//=======================================================================
      inline const std::vector<Quantity*> & GetQtyIntVector() const { 
	return _QtyInternalVector;
      }
      
      void SetQtyIntVector(const std::vector<Quantity*> & vect_in) { 
	_QtyInternalVector = vect_in;
	init_sys();
      }

	    std::vector<std::string> _var_names;                   //MultilevelSolution//
	  void initVarNames();                                     //MultilevelSolution//

	  std::vector<double>      _refvalue;        //MultilevelSolution//
          void initRefValues();                      //MultilevelSolution//

	  void init_sys();     //System//

//=======================================================================
//====== Attributes of the equation ====
//=======================================================================
  const uint      _NoLevels;   ///< level number      //System//

//=======================================================================
//========= MULTIGRID FUNCTIONS (Vectors + A,R,P) ======== (procs,levels) 
//=======================================================================
    void ReadMGOps(const std::string output_path); // LinearEquation  (each level)
    void ReadMatrix(const std::string& name); // LinearEquation  (each level)
    void ReadProl(const std::string& name);   // LinearEquation  (each level)
    void ReadRest(const std::string& name);   // LinearEquation  (each level)
    void ComputeMatrix();                     // LinearEquation  (each level)
    void ComputeProl();                       // LinearEquation  (each level)
    void ComputeRest();                       // LinearEquation  (each level)

 virtual  void GenMatRhs(const uint Level) = 0;  //System//
          void MGSolve(double Eps,int MaxIter, const uint Gamma=DEFAULT_MG_GAMMA, const uint Nc_pre=DEFAULT_NC_PRE,const uint Nc_coarse=DEFAULT_NC_COARSE,const uint Nc_post=DEFAULT_NC_POST);  //LinearImplicitSystem//
        double MGStep(int Level,double Eps1,int MaxIter, const uint Gamma, const uint Nc_pre,const uint Nc_coarse,const uint Nc_post);                                                          //LinearImplicitSystem//
          void MGCheck(int Level) const;

  DofMap  _dofmap;  //// LinearEquation (each level)
  
//=======================================================================
// ============ INITIAL CONDITIONS of the equation ====== (procs,levels) ==    
// ========================================================
          void    Initialize();           //MultilevelSolution  //this uses x and fills in x_old at all levels
          
//=======================================================================
//==== BOUNDARY CONDITIONS of the equation ========= (procs,levels) ==     //MultilevelSolution
//=======================================================================
    int   *_bc;         //==== NODAL DIRICHLET ======== ///< boundary conditions map (top level)  // POINTWISE(NODAL) FLAG for the BOUNDARY DOFS = FLAG for the tEST FUNCTIONS //TODO this should be PrintNumericVector of the equation, integer instead of double! do it when you make it parallel especially! //Later on I will do a bc for every level, considering the ELEMENT DOFS
    int  **_bc_fe_kk;   //==== FE KK DIRICHLET ========
          void    GenerateBdc();          //MultilevelSolution
 //====PENALTY DIRICHLET ======Elem BC=====================
   uint  _Dir_pen_fl;         ///flag for penalty with Dirichlet (0=no penalty, 1=yes penalty) //this penalty is for ALL the QUADRATIC variables //could we do a penalty only for ux and not for uy and uz?
    int  ***_elem_bc;        ///[LEVELS][IPROC][2xELEMENTSxLEV&PROC]
 double  ***_elem_val_norm;  ///[LEVELS][IPROC][1xELEMENTSxLEV&PROC]
 double  ***_elem_val_tg;    ///[LEVELS][IPROC][(1or3)xELEMENTSxLEV&PROC]
 static const int _number_tang_comps[3];  //  {0,1,3} Number of tangential components in 1D,2D,3D (see BC in general domains) : use it only for 2D and 3D
          void    GenerateBdcElem();
          void  clearElBc();
  virtual void elem_bc_read(const double * el_xm, int& surf_id, double * value,int * el_flag) const = 0;
          void Bc_GetElFlagValLevSubd(const uint Level,const uint isubd,const uint iel,int* el_flag,double* el_value ) const;
          void Bc_ConvertToDirichletPenalty(const uint elem_dim, const uint ql, uint* bc) const;
          void Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(const uint *bc_eldofs,const CurrentQuantity &Velold_in,const CurrentQuantity& press_in,uint &press_fl) const;
//========= treating NumericVectors, related to Dirichlet Boundary Conditions! =======
          void Bc_ScaleDofVec(NumericVector * myvec,  double ScaleFac);
          void Bc_AddDofVec(NumericVector* myvec, NumericVector* myvec2 );
          void Bc_AddScaleDofVec(NumericVector* vec_in,NumericVector* vec_out,const double ScaleFac );

protected:
  
  const FemusInputParser<double>  & _phys;   //passed from MultilevelProblem
  const MultiLevelMeshTwo         & _mesh;   //passed from MultilevelProblem

  std::vector<Quantity*>          _QtyInternalVector;  //System//

};



} //end namespace femus



#endif