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

#ifndef __MultiLevelSolution_hpp__
#define __MultiLevelSolution_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelMesh.hpp"
#include "Solution.hpp"
#include "ParallelObject.hpp"
#include "FElemTypeEnum.hpp"
#include <vector>


namespace femus {



typedef double (*initfunc) (const double &x, const double &y, const double &z);

/**
 * This class is a black box container to handle multilevel solutions.
 */

class MultiLevelSolution : public ParallelObject {
  
private:
 
  bool _bdc_func_set;
  unsigned short  _gridn;
  vector <int>    _SolType;
  vector <FEFamily> _family;
  vector <FEOrder> _order; 
  vector <char*>  _SolName;
  vector <char*>  _BdcType;
  vector <int>    _SolTmorder;
  vector <bool>   _TestIfPressure;
  vector <bool>   _TestIfDisplacement;
  
  /** Array of solution */
  vector <Solution*>  _solution;

public:
 
  Solution* GetSolutionLevel(const unsigned i) {return _solution[i];};
   
  /** Multilevel mesh */
  MultiLevelMesh* _ml_msh;
  
  /** Constructor */
  MultiLevelSolution(MultiLevelMesh *ml_msh);

  /** Destructor */
  ~MultiLevelSolution();
     
  // Vector handling functions
  void AddSolution(const char name[], const FEFamily fefamily, const FEOrder order, unsigned tmorder=0, const bool &Pde_type=1);
  
  void AssociatePropertyToSolution(const char solution_name[], const char solution_property[]);

  void ResizeSolutionVector( const char name[]);

  void Initialize(const char name[], initfunc func = NULL);

  unsigned GetIndex(const char name[]) const;

  unsigned GetSolType(const char name[]);

  unsigned GetSolutionSize(){ return _SolType.size();};

  vector <char*>  GetSolName(){return _SolName;};

  vector <int>  GetSolType(){return _SolType;};
  
  void BuildProlongatorMatrix(unsigned gridf, unsigned SolIndex);
  
  // boundary condition function pointer
  bool (*_SetBoundaryConditionFunction) (const double &x, const double &y, const double &z,const char name[], 
                                         double &value, const int FaceName, const double time);
  
  void AttachSetBoundaryConditionFunction ( bool (* SetBoundaryConditionFunction) (const double &x, const double &y, const double &z,const char name[], 
										   double &value, const int FaceName, const double time) );
  
  void GenerateBdc(const char name[], const char bdc_type[]="Steady");
  void UpdateBdc(const double time);
  void GenerateBdc(const unsigned int k, const double time);
   
  char* GetSolutionName(unsigned i){return _SolName[i];};
  int   GetSolutionType(unsigned i){return _SolType[i];};
  unsigned GetSolutionType(const char name[]);
  char* GetBdcType(unsigned i){return _BdcType[i];};
  int   GetSolutionTimeOrder(unsigned i){return _SolTmorder[i];};
  bool  TestIfSolutionIsPressure(unsigned i){return _TestIfPressure[i];};
  bool  TestIfSolutionIsDisplacemenet(unsigned i){return _TestIfDisplacement[i];};
  
};


} //end namespace femus



#endif

