#ifndef __mgsolbase__
#define __mgsolbase__

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


namespace femus {




class Files        ;
class Physics      ;
class EquationsMap ;
class MeshTwo      ;
class FEElemBase   ;
class Quantity     ;
class QuantityLocal;

class SparseMatrix;
class NumericVector;
class LinearEquationSolver;



class EqnBase  {

public:

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
  std::vector<NumericVector *> _x_tmp;
  
          void  initVectors();                                      ///initialize vectors                                                               //System//

	  void PrintVector(std::string namefile);   ///prints on a "Quadratic-Linearized" Mesh //TODO this should be PrintNumericVector of the equation //Writer//
          void  ReadVector(std::string namefile);                       ///read from a "Quadratic-Linearized" Mesh                                      //Writer/Reader// 

//======= Linear Solvers for every Level ============
  LinearEquationSolver **_solver;     ///(each level)

//=======================================================================
//======= Quantities =========
//=======================================================================
      std::vector<Quantity*>          _QtyInternalVector;  //System//
  
//=======================================================================
  //====== Attributes of the equation ====
//=======================================================================
  const std::string _eqname;   ///< equation name     //System//
  uint               _iproc;   ///< processor rank    //ParallelObject//
  const uint      _NoLevels;   ///< level number      //System//

//=======================================================================
// CONSTRUCTOR / DESTRUCTOR
//=======================================================================
  EqnBase(std::vector<Quantity*> int_map_in,
	  EquationsMap& equations_map,
          std::string eq_name_in="Base",
	  std::string varname_in="u");   //System//
  
  virtual ~EqnBase();                    //System//
  
//=======================================================================
//==== DOF MAP of the equation ============ (procs,levels) ==============   //// LinearEquation (each level)
//=======================================================================
//====== data =======
  int    **  _node_dof; ///< dof map
  uint  *    _Dim;            //number of dofs per level
  uint **    _DofNumLevFE;     
  uint **    _DofOffLevFE;     
  uint ***   _DofLocLevProcFE;     
  uint       _nvars[QL];      ///  number of SCALAR variables          
  uint       _VarOff[QL];
  uint       _n_vars;           ///< number of SCALAR variables
  std::string *_var_names;     /// variable names of every SCALAR variable
//====== functions =======
          void initNVars();
	  void initVarNames(std::string varname_in);
          void ComputeMeshToDof();
          void PrintMeshToDof() const;

  double      *_refvalue;          ///reference values of every SCALAR variable
          void initRefValues();

//=======================================================================
//========= MULTIGRID FUNCTIONS (Vectors + A,R,P) ======== (procs,levels) 
//=======================================================================
    void ReadMGOps();                         // LinearEquation  (each level)
    void ReadMatrix(const std::string& name); // LinearEquation  (each level)
    void ReadProl(const std::string& name);   // LinearEquation  (each level)
    void ReadRest(const std::string& name);   // LinearEquation  (each level)
    void ComputeMatrix();                     // LinearEquation  (each level)
    void ComputeProl();                       // LinearEquation  (each level)
    void ComputeRest();                       // LinearEquation  (each level)

 virtual  void GenMatRhs(const uint Level) = 0;  //System//
          void MGSolve(double Eps,int MaxIter, const uint Gamma=DEFAULT_MG_GAMMA, const uint Nc_pre=DEFAULT_NC_PRE,const uint Nc_coarse=DEFAULT_NC_COARSE,const uint Nc_post=DEFAULT_NC_POST);  //PetscLinearSolverM//
        double MGStep(int Level,double Eps1,int MaxIter, const uint Gamma, const uint Nc_pre,const uint Nc_coarse,const uint Nc_post);                                                          //PetscLinearSolverM//
          void MGCheck(int Level) const;


//=======================================================================
// ============ INITIAL CONDITIONS of the equation ====== (procs,levels) ==
// ========================================================
          void    GenIc();           //MultilevelSolution, Initialize function
  virtual void  ic_read(const double * xp, double * ic,const double * el_xm) const = 0; //TODO see what parameters can be made constant
          
//=======================================================================
//==== BOUNDARY CONDITIONS of the equation ========= (procs,levels) ==     //MultilevelSolution
//=======================================================================
    int   *_bc;         //==== NODAL DIRICHLET ======== ///< boundary conditions map (top level)  // POINTWISE(NODAL) FLAG for the BOUNDARY DOFS = FLAG for the tEST FUNCTIONS //TODO this should be PrintNumericVector of the equation, integer instead of double! do it when you make it parallel especially! //Later on I will do a bc for every level, considering the ELEMENT DOFS
    int  **_bc_fe_kk;   //==== FE KK DIRICHLET ========
          void    GenBc();
  virtual void  bc_read(const double * xp,const double * normal, int * bc) const = 0;
          void  PrintBc(std::string namefile);      
 //====PENALTY DIRICHLET ======Elem BC=====================
   uint  _Dir_pen_fl;         ///flag for penalty with Dirichlet (0=no penalty, 1=yes penalty) //this penalty is for ALL the QUADRATIC variables //could we do a penalty only for ux and not for uy and uz?
    int  ***_elem_bc;        ///[LEVELS][IPROC][2xELEMENTSxLEV&PROC]
 double  ***_elem_val_norm;  ///[LEVELS][IPROC][1xELEMENTSxLEV&PROC]
 double  ***_elem_val_tg;    ///[LEVELS][IPROC][(1or3)xELEMENTSxLEV&PROC]
 int _number_tang_comps[3];  //  {0,1,3} Number of tangential components in 1D,2D,3D (see BC in general domains) : use it only for 2D and 3D
          void    GenElBc();
          void  clearElBc();
  virtual void elem_bc_read(const double * el_xm, int& surf_id, double * value,int * el_flag) const = 0;
          void Bc_GetElFlagValLevSubd(const uint Level,const uint isubd,const uint iel,int* el_flag,double* el_value ) const;
          void Bc_ConvertToDirichletPenalty(const uint elem_dim, const uint ql, uint* bc) const;
          void Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(const uint *bc_eldofs,const QuantityLocal &Velold_in,const QuantityLocal& press_in,uint &press_fl) const;
//========= treating NumericVectors, related to Dirichlet Boundary Conditions! =======
          void Bc_ScaleDofVec(NumericVector * myvec,  double ScaleFac);
          void Bc_AddDofVec(NumericVector* myvec, NumericVector* myvec2 );
          void Bc_AddScaleDofVec(NumericVector* vec_in,NumericVector* vec_out,const double ScaleFac );

protected:
  
//=======================================================================
// ====== data pointer ==========
//=======================================================================
  Files                     & _files;
  Physics                   & _phys;           //passed from MultilevelProblem
  MeshTwo                   & _mesh;
  std::vector<FEElemBase*>  &  _AbstractFE;
  EquationsMap              & _eqnmap;


};



} //end namespace femus



#endif