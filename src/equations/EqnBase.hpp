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
class LinearSolverM;



class EqnBase  {

public:

//=======================================================================
//======== Vectors =============== (procs,levels) ==
//=======================================================================

  std::vector<NumericVector *> _b;   ///< rhs b
  std::vector<NumericVector *> _x;   ///< solution x
  std::vector<NumericVector *> _res; ///< residual

  std::vector<NumericVector *> _x_old;
  std::vector<NumericVector *> _x_oold;
  std::vector<NumericVector *> _x_tmp;

          void  initVectors();                                      ///initialize vectors
          void PrintVector(std::string namefile);   ///prints on a "Quadratic-Linearized" Mesh //TODO this should be PrintNumericVector of the equation 
          void  ReadVector(std::string namefile);                       ///read from a "Quadratic-Linearized" Mesh 

//=======================================================================
//======= Quantities =========
//=======================================================================
      std::vector<Quantity*>          _QtyInternalVector;
  
//=======================================================================
  //====== Attributes of the equation ====
//=======================================================================
  const std::string _eqname;   ///< equation name
  uint               _iproc;   ///< processor rank
  const uint      _NoLevels;   ///< level number

//=======================================================================
// CONSTRUCTOR / DESTRUCTOR
//=======================================================================
  EqnBase(std::vector<Quantity*> int_map_in,
	  EquationsMap& equations_map,
          std::string eq_name_in="Base",
	  std::string varname_in="u");
  
  virtual ~EqnBase();
  
//=======================================================================
//==== DOF MAP of the equation ============ (procs,levels) ==============
//=======================================================================
//====== data =======
  int    **  _node_dof; ///< dof map
  uint  *    _Dim;            ///< dimension //number of dofs per level
  uint **    _DofNumLevFE;     
  uint **    _DofOffLevFE;     
  uint ***   _DofLocLevProcFE;     
  uint       _nvars[QL];      ///  number of SCALAR variables          
  uint       _VarOff[QL];
  uint       _n_vars;           ///< number of SCALAR variables
  std::string *_var_names;     /// variable names of every SCALAR variable
  double      *_refvalue;          ///reference values of every SCALAR variable
//====== functions =======
          void initNVars();
	  void initVarNamesRefValues(std::string varname_in);
          void ComputeMeshToDof();     ///  distributing dof function
          void PrintMeshToDof() const;


//=======================================================================
//========= MULTIGRID FUNCTIONS (Vectors + A,R,P) ======== (procs,levels) 
//=======================================================================
    void initMGOps();                                           ///initialize A, R, P and read
    void ReadMatrix(const std::string& name);  ///< Reading matrix A
    void ReadProl(const std::string& name);    ///< Reading Prolongation
    void ReadRest(const std::string& name);    ///< Reading Restriction
    void ComputeMatrix();
    void ComputeProl();
    void ComputeRest();
    void PrintOneVarMatrixHDF5(std::string name, std::string groupname, uint** n_nodes_all, int count,int* Mat,int* len,int* len_off,int type1, int type2, int* FELevel ) const;
    void PrintOneVarMGOperatorHDF5(std::string filename, std::string groupname, uint* n_dofs_lev, int count,int* Rest,double* values,int* len,int* len_off, int FELevel, int FELevel2, int fe) const;

 virtual  void GenMatRhs(const uint Level) = 0;  //every equation must have one, explicitly!
        double MGTimeStep(const uint iter);                    ///< MG time step solver (backward Euler)
          void MGSolve(double Eps,int MaxIter, const uint Gamma=DEFAULT_MG_GAMMA, const uint Nc_pre=DEFAULT_NC_PRE,const uint Nc_coarse=DEFAULT_NC_COARSE,const uint Nc_post=DEFAULT_NC_POST);    ///< MultiGrid Solver
        double MGStep(int Level,double Eps1,int MaxIter, const uint Gamma, const uint Nc_pre,const uint Nc_coarse,const uint Nc_post);    ///< MultiGrid Step
          void MGCheck(int Level) const;                                           ///< Check Operators


//=======================================================================
// ============ INITIAL CONDITIONS of the equation ====== (procs,levels) ==
// ========================================================
          void    GenIc();
  virtual void  ic_read(const double * xp, double * ic,const double * el_xm) const = 0; //TODO see what parameters can be made constant
          
//=======================================================================
//==== BOUNDARY CONDITIONS of the equation ========= (procs,levels) ==
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

//=======================================================================
//======= MG: Linear Solvers for every Level ============
//=======================================================================
  LinearSolverM **_solver;     ///< linear system solver type (each level)

protected:
  
//=======================================================================
// ====== data pointer ==========
//=======================================================================
  Files                     & _files;        ///<  file class pointer
  Physics                   & _phys;         ///<  parameter class pointer
  MeshTwo                   & _mesh;         ///<  mesh pointer
  std::vector<FEElemBase*>  &  _AbstractFE;  ///<  FE
  EquationsMap              & _eqnmap;       ///<  equation map  pointer

//=======================================================================
//======== MG Ops ============ (procs,levels) ====
//=======================================================================
  std::vector<SparseMatrix  *> _A;    ///< Matrix A
  std::vector<SparseMatrix *> _Rst; ///< Restrictor
  std::vector<SparseMatrix *> _Prl; ///< Prolongation

};



} //end namespace femus



#endif