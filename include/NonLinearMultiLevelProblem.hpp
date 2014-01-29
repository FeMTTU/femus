#ifndef __NonLinearMultiLevelProblem_hpp__
#define __NonLinearMultiLevelProblem_hpp__

#include "LinSysPde.hpp"
#include "Solution.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include <b64/b64.h>

#include <vector>
#include <map>
using std::map;

class elem_type;
class LinearSolverM;
//********************  class solver**************************
//************************************************************
class NonLinearMultiLevelProblem {

private:
  vector < SparseRectangularMatrix* > ProlQitoQj_[3][3];
  vector < map <unsigned,bool> > index;
  vector <unsigned> elr_old;
  int _moving_mesh;
  std::vector<std::string> _moving_vars;
  
  unsigned _non_linear_algorithm;
  bool _is_nonlinear;
  double _non_linear_toll;
  
  bool _init_func_set;
  bool _bdc_func_set;
  
 protected:
  int _nprocs;
  int _iproc;
  unsigned short gridn, gridr;
  bool _test_time;
  unsigned _time_step;
  unsigned _ntime_steps;
  double _dt;
  double _time;
  unsigned _time_step0;
  bool _VankaSchur;
  bool _Schur;
  bool _VankaIsSet;
  short unsigned _NSchurVar;

 
  vector <char*> SolName;
  vector <int> SolTmorder;
  vector <char*> BdcType;
  vector <bool> TestIfPressure;
  
 
  
  
  

  //pointer function to the assemble function
  int (* _assemble_function)(NonLinearMultiLevelProblem &mg, unsigned level, const unsigned &gridn);
  
  bool (* _SetBoundaryConditionFunction) (const double &x, const double &y, const double &z,const char name[], 
                            double &value, const int FaceName, const double time);
  
  double (* _InitVariableFunction)(const double &x, const double &y, const double &z,const char name[]);

 public:
   
  vector <char*> _PDEName;
  vector <unsigned> _PDE_MGIndex;
  void AddPDE(const char name[]);
  unsigned GetPDEIndex(const char name[]) const; 
   
   
   
  vector <int> SolType;
  vector< vector <unsigned> > MGIndex;
  vector <unsigned> VankaIndex;
  const elem_type *type_elem[6][5]; 

  
  ///Array of linear solver
  vector<vector <LinearSolverM*> > _LinSolver;
  vector <Solution*>  _solution;
  vector <mesh*> _msh;
  
  
  ///Physical Object
  Fluid* _fluid;
  Solid* _solid;

  ///Constructor
  NonLinearMultiLevelProblem(const unsigned short &igridn, const unsigned short &igridr, const char mesh_file[],
            const char GaussOrder[], const double Lref=1., bool (* SetRefinementFlag)(const double &x, const double &y, const double &z, 
                                   const int &ElemGroupNumber,const int &level)=NULL );

  ~NonLinearMultiLevelProblem();

  //Attaching Functions
  ///Provides a method for filling the Matrix and the Residual vector
  void AttachAssembleFunction ( int (*function)(NonLinearMultiLevelProblem &mg, unsigned level, const unsigned &gridn));
  
  void AttachSetBoundaryConditionFunction ( bool (* SetBoundaryConditionFunction) (const double &x, const double &y, const double &z,const char name[], 
                            double &value, const int FaceName, const double time) );
  
  
  void AttachInitVariableFunction ( double (* _InitVariableFunction)(const double &x, const double &y, 
                                                                     const double &z,const char name[]) );
  
  //utilities
  void MarkStructureNode();
  double ComputeL2norm();
  bool GetConvergence(const char pdename[],const unsigned gridn);
  int ComputeBdStress(int bd, double cforce[3]);
  int ComputeBdIntegral(const char pdename[],const char var_name[], const unsigned & kel, 
                         const unsigned & jface, unsigned level, unsigned dir);
  unsigned GetNumberOfGrid();
  
  // Config
  void SetMatrixProperties(const char pdename[], const char property[]);
  void AddStabilization(const char pdename[], const bool stab=false, const double compressibility=0.);
  void SetNonLinearAlgorithm(bool isnonlinear=false, const char nnlinalg[]="Linear", const double nl_toll=1.e-03);
  unsigned GetNonLinearAlgorithm();
  bool GetNonLinearCase();

  //Boundary conditions
  void GenerateBdc(const char name[], const char bdc_type[]="Steady");

  // Solver  
  int FullMultiGrid(const char pdename[], unsigned const &ncycle,  unsigned const &npost, unsigned const &npre, 
                    const char mg_type[]="F-Cycle");
  int FreeMultigrid();
  void SetSmoother(const char smoothername[]);
  void SetVankaSchurOptions(bool VankaSchur=1, bool Schur=0, short unsigned NSchurVar=1);
  void SetTolerances(const char pdename[],const double rtol,const double atol,const double divtol, const unsigned maxits);
  void SetSchurTolerances(const char pdename[], const double rtol,const double atol,const double divtol, const unsigned maxits);
  void SetSolverFineGrids(const char pdename[], const char solvertype[] = "GMRES");
  void SetPreconditionerFineGrids(const char pdename[], const char preconditionertype[] = "LU");
  void SetDimVankaBlock(const char pdename[], unsigned const dim_vanka_block);
  void SetDimVankaBlock(const char pdename[], const char dim_vanka_block[] = "All");
  
  // Fluid and Solid object
  void Add_Fluid(Fluid *fluid);
  void Add_Solid(Solid *solid);

  // Vector handling functions
  void AddSolutionVector(const char name[], const char type[],unsigned tmorder=0, const bool &PDE_type=1);
  void AssociatePropertyToSolution(const char solution_name[], const char solution_property[]);
  void ResizeSolutionVector( const char name[]);
  void CheckVectorSize(const unsigned &i);
  void Initialize(const char name[]);
  unsigned GetIndex(const char name[]) const;
  unsigned GetTmOrder(const unsigned i);

  //Index
  void CreateMGStruct();
  void DeleteMGStruct();
  void ClearMGIndex();
  void AddToMGIndex(const char pdename[], const char solname[]);
  unsigned GetMGIndex(const char pdename[], const char name[]);
  unsigned GetSolType(const char name[]) ;
  //char* GetMGIndexName(unsigned index);
  void ClearVankaIndex();
  void AddToVankaIndex( const char pdename[], const char solname[]);

  //MultiGrid tools
  int Restrictor(const unsigned &gridf, const unsigned &ipde);
  int Prolungator(const unsigned &gridf, const unsigned &ipde);
  void ProlungatorSol(const char pdename[], unsigned gridf);
  int BuildProlungatorMatrix(unsigned gridf,const char pdename[]);
  void BuildProlungatorMatrix(unsigned gridf, unsigned SolIndex);
  void BuildProlungatorMatrices();

  //printing and reading solution Functions 
  int printsol_gmv_binary(const char name[]="linear",unsigned igridn=0, bool debug=0) const;
  void printsol_vtu_inline(const char name[],std::vector<std::string>& vars) const;
  void printsol_xdmf_hdf5(const char name[],std::vector<std::string>& vars)const;
//  hid_t print_Dhdf5(hid_t file,const std::string & name, hsize_t dimsf[],double data[]);
  void SetMovingMesh(std::vector<std::string>& myss);
  
};

#endif
// kate: indent-mode cstyle; space-indent on; indent-width 2;
