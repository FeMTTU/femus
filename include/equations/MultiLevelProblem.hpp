/*=========================================================================

 Program: FEMUS
 Module: MultiLevelProblem
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __MultiLevelProblem_hpp__
#define __MultiLevelProblem_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "LinearEquation.hpp"
#include "Solution.hpp"
#include "Parameters.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include <b64/b64.h>

#include <vector>
#include <map>
using std::map;

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class elem_type;
class LinearEquationSolver;
class System;


typedef double (*initfunc) (const double &x, const double &y, const double &z);

/**
* This class is a black box container to handle multilevel problems.
*/

class MultiLevelProblem {

private:
  vector < SparseMatrix* > ProlQitoQj_[3][3];
  vector < map <unsigned,bool> > index;
  vector <unsigned> elr_old;
  int _moving_mesh;
  std::vector<std::string> _moving_vars;
  
//   unsigned _non_linear_algorithm;
//   bool _is_nonlinear;
//   double _non_linear_toll;
  
  bool _init_func_set;
  bool _bdc_func_set;
  
  void GenerateBdc(const unsigned int k, const double time);
  
    
protected:
  
  int _nprocs;
  int _iproc;
  unsigned short _gridn, _gridr;
//   bool _test_time;
//   unsigned _time_step;
//   unsigned _ntime_steps;
//   double _dt;
//   double _time;
//   unsigned _time_step0;
//   bool _Schur;
//   bool _VankaIsSet;
//   short unsigned _NSchurVar;
  
//   unsigned int _print_step;
 
  vector <int> SolTmorder;
  vector <char*> BdcType;
  vector <bool> _TestIfPressure;
//   vector <bool> _TestIfPdeHasDisplacement;
  
  
  //pointer function to the assemble function
//   void (*_assemble_function)(MultiLevelProblem &mg, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix);
  
  bool (*_SetBoundaryConditionFunction) (const double &x, const double &y, const double &z,const char name[], 
                            double &value, const int FaceName, const double time);
  
 public:
  vector <bool> _TestIfDisplacement;
   
   /** Data structure holding the systems. */
  std::map<std::string, System*> _systems;

  /** Typedef for system iterators */
  typedef std::map<std::string, System*>::iterator       system_iterator;

  /** Typedef for constatnt system iterators */
  typedef std::map<std::string, System*>::const_iterator const_system_iterator; 
   
//   vector <char*> _PdeName;
//   vector <unsigned> _PdeIndex;
//   void AddPde(const char name[]);
//   unsigned GetPdeIndex(const char name[]) const; 
  
  /** Add the system of type \p system_type named \p name to the systems array. */
  virtual System & add_system (const std::string& system_type, const std::string& name);

  /** Add the system named \p name to the systems array. */
  template <typename T_sys> T_sys & add_system (const std::string& name);
  
  /**
   * @returns a constant reference to the system named \p name.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> ("sys");
   * is an example of how the method might be used
   */
  template <typename T_sys>
  const T_sys & get_system (const std::string& name) const;

  /**
   * @returns a writeable referene to the system named \p name.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> ("sys");
   * is an example of how the method might be used
   */
  template <typename T_sys>
  T_sys & get_system (const std::string& name);

  /**
   * @returns a constant reference to system number \p num.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> (0);
   * is an example of how the method might be used
   */
  template <typename T_sys>
  const T_sys & get_system (const unsigned int num) const;

  /**
   * @returns a writeable referene to the system number \p num.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> (0);
   * is an example of how the method might be used
   */
  template <typename T_sys>
  T_sys & get_system (const unsigned int num);

  /**
   * @returns a constant reference to the system named \p name.
   */
  const System & get_system (const std::string& name) const;

  /**
   * @returns a writeable referene to the system named \p name.
   */
  System & get_system (const std::string& name);

  /**
   * @returns a constant reference to system number \p num.
   */
  const System & get_system (const unsigned int num) const;

  /**
   * @returns a writeable referene to the system number \p num.
   */
  System & get_system (const unsigned int num);
  
  
  /** @returns the number of equation systems. */
  unsigned int n_systems() const;
  
  /** Clear all the Sytems PDE structures */
  void clear();
  
  /** init the system pde structures */
  void init();

  vector <char*> SolName;
  vector <int> SolType;
//   vector< vector <unsigned> > _SolPdeIndex;
//   vector <unsigned> VankaIndex;
  const elem_type *type_elem[6][5]; 

  /** Array of linear solver */
  vector<vector <LinearEquationSolver*> > _LinSolver;
  
  /** Array of solution */
  vector <Solution*>  _solution;
  
  /** Array of mesh */
  vector <mesh*> _msh;
  
  /** Data structure holding arbitrary parameters. */
  Parameters parameters;

  /** Constructor */
  MultiLevelProblem(const unsigned short &igridn, const unsigned short &igridr, const char mesh_file[],
            const char GaussOrder[], const double Lref=1., bool (* SetRefinementFlag)(const double &x, const double &y, const double &z, 
                                   const int &ElemGroupNumber,const int &level)=NULL );

  /** Destructor */
  ~MultiLevelProblem();

  //Attaching Functions
  /** Provides a method for filling the Matrix and the Residual vector */
//   void AttachAssembleFunction ( void (*function)(MultiLevelProblem &mg, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix));
  
  void AttachSetBoundaryConditionFunction ( bool (* SetBoundaryConditionFunction) (const double &x, const double &y, const double &z,const char name[], 
                            double &value, const int FaceName, const double time) );
  
  //utilities
  void MarkStructureNode();
//   bool GetConvergence(const char pdename[],const unsigned gridn);

  int ComputeBdIntegral(const char pdename[],const char var_name[], const unsigned & kel, 
                         const unsigned & jface, unsigned level, unsigned dir);
  unsigned GetNumberOfGrid();
  unsigned GetNumberOfGridTotallyRefined();

  
  // Config
//   void SetMatrixProperties(const char pdename[], const char property[]);
//   void AddStabilization(const char pdename[], const bool stab=false, const double compressibility=0.);
//   void SetNonLinearAlgorithm(bool isnonlinear=false, const char nnlinalg[]="Linear", const double nl_toll=1.e-03);
//   unsigned GetNonLinearAlgorithm();
//   bool GetNonLinearCase();

  //Boundary conditions
  void GenerateBdc(const char name[], const char bdc_type[]="Steady");
  void UpdateBdc(const double time);
//   void SetDirichletBCsHandling(const char pdename[],const char DirichletMode[]);

  //* Multigrid Solver  */
//   void Solve(const char pdename[], unsigned const &ncycle,  unsigned const &npost, unsigned const &npre, 
//              const char mg_type[]="F-Cycle", const bool &linear=0);
//   void FreeMultigrid();
//   void SetSmoother(const char smoothername[]);
//   void SetVankaSchurOptions(bool Schur=0, short unsigned NSchurVar=0);
//   void SetTolerances(const char pdename[],const double rtol,const double atol,const double divtol, const unsigned maxits);
//   void SetSchurTolerances(const char pdename[], const double rtol,const double atol,const double divtol, const unsigned maxits);
//   void SetSolverFineGrids(const char pdename[], const char solvertype[] = "GMRES");
//   void SetPreconditionerFineGrids(const char pdename[], const char preconditionertype[] = "LU");
//   void SetDimVankaBlock(const char pdename[], unsigned const dim_vanka_block);
//   void SetDimVankaBlock(const char pdename[], const char dim_vanka_block[] = "All");
  
  // Vector handling functions
  void AddSolution(const char name[], const char type[],unsigned tmorder=0, const bool &Pde_type=1);
  void AssociatePropertyToSolution(const char solution_name[], const char solution_property[]);
  void ResizeSolutionVector( const char name[]);
//   void CheckVectorSize(const unsigned &i);
  void Initialize(const char name[], initfunc func = NULL);
  unsigned GetIndex(const char name[]) const;
//   unsigned GetTmOrder(const unsigned i);

  //Index
//   void CreatePdeStructure();
//   void DeletePdeStructure();
//   void ClearSolPdeIndex();
//   void AddSolutionToSolPdeIndex(const char pdename[], const char solname[]);
//   unsigned GetSolPdeIndex(const char pdename[], const char name[]);
  unsigned GetSolType(const char name[]) ;
//   void ClearVankaIndex();
//   void AddToVankaIndex( const char pdename[], const char solname[]);

  //MultiGrid tools
//   void Restrictor(const unsigned &ig, const unsigned &ipde, const unsigned &igridn, 
//                   const unsigned &nlcycle,const unsigned &lcycle, const bool &flagmc);
//   void Prolongator(const unsigned &gridf, const unsigned &ipde);
//   void ProlongatorSol(const char pdename[], unsigned gridf);
//   void BuildProlongatorMatrix(unsigned gridf,const char pdename[]);
  void BuildProlongatorMatrix(unsigned gridf, unsigned SolIndex);
  void BuildProlongatorMatrices();

  //printing and reading solution Functions 
  void printsol_gmv_binary(const char name[]="linear",unsigned igridn=0, bool debug=0) const;
  void printsol_vtu_inline(const char name[], std::vector<std::string>& vars, const unsigned time_step=0) const;
  void printsol_xdmf_hdf5(const char name[],std::vector<std::string>& vars)const;
  void SetMovingMesh(std::vector<std::string>& myss);
  void printsol_xdmf_archive(const char type[]) const;
 
//   char* GetThisPdeName(const unsigned &ipde){ return _PdeName[ipde];}
};

template <typename T_sys>
inline
T_sys & MultiLevelProblem::add_system (const std::string& name)
{
  T_sys* ptr = NULL;

  if (!_systems.count(name))
    {
      ptr = new T_sys(*this, name, this->n_systems());

      _systems.insert (std::make_pair(name, ptr));

    }
  else
    {
      // We now allow redundant add_system calls, to make it
      // easier to load data from files for user-derived system
      // subclasses
      std::cerr << "ERROR: There was already a system"
              << " named " << name
              << std::endl;

      ptr = &(this->get_system<T_sys>(name));
    }

  // Return a dynamically casted reference to the newly added System.
  return *ptr;
}


template <typename T_sys>
inline
const T_sys & MultiLevelProblem::get_system (const std::string& name) const
{
  const_system_iterator pos = _systems.find(name);

  // Check for errors
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: no system named \"" << name << "\" found!"
                    << std::endl;
    }

  // Attempt dynamic cast
  return *static_cast<T_sys*>(pos->second);
}

template <typename T_sys>
inline
T_sys & MultiLevelProblem::get_system (const std::string& name)
{
  system_iterator pos = _systems.find(name);

  // Check for errors
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: no system named " << name << " found!"
                    << std::endl;
    }

  // Attempt dynamic cast
  return *static_cast<T_sys*>(pos->second);
}

inline
const System & MultiLevelProblem::get_system (const std::string& name) const
{
  return this->get_system<System>(name);
}


inline
System & MultiLevelProblem::get_system (const std::string& name)
{
  return this->get_system<System>(name);
}


inline
const System & MultiLevelProblem::get_system (const unsigned int num) const
{
  return this->get_system<System>(num);
}


inline
System & MultiLevelProblem::get_system (const unsigned int num)
{
  return this->get_system<System>(num);
}

inline
unsigned int MultiLevelProblem::n_systems () const
{
  return static_cast<unsigned int>(_systems.size());  //libmesh static_cast_int
}

#endif

