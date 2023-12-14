/*=========================================================================

Program: FEMUS
Module: MultiLevelSolution
Authors: Eugenio Aulisa, Simone Bnà, Giorgio Bornia

Copyright (c) FEMTTU
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelSolution.hpp"
#include "ElemType.hpp"
#include "NumericVector.hpp"
#include "FemusConfig.hpp"
#include "Solution_config.hpp"

#include "ParsedFunction.hpp"



//C++ include
#include <iostream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>


namespace femus {



//---------------------------------------------------------------------------------------------------
  MultiLevelSolution::~MultiLevelSolution()  {
    clear();
  }

//---------------------------------------------------------------------------------------------------
  MultiLevelSolution::MultiLevelSolution(/*const*/ MultiLevelMesh* ml_msh) :
    _gridn(ml_msh->GetNumberOfLevels()),
    _mlMesh(ml_msh) {
      
      
    _solution.resize(_gridn);

    for(unsigned i = 0; i < _gridn; i++) {
      _solution[i] = new Solution(GetMLMesh()->GetLevel(i));
    }

    _bdcFuncSet = false;
    _bdcFuncSetMLProb = false;
    _useParsedBCFunction = false;

    _mlBCProblem = NULL;

    _FSI = false;

    _writer = NULL;

  }


//---------------------------------------------------------------------------------------------------
// this is the destructor that can be called explicitly, instead of the automatic destructor
  void MultiLevelSolution::clear() {

    for(unsigned i = 0; i < _gridn; i++) {
      _solution[i]->FreeSolutionVectors();
      delete _solution[i];
    }

    for(unsigned i = 0; i < _solName.size(); i++) delete [] _solName[i];

    for(unsigned i = 0; i < _solName.size(); i++) delete [] _bdcType[i];

    if(_writer != NULL) delete _writer;



  }






  void MultiLevelSolution::AddSolutionLevel() {
    
    // add level solution
    _solution.resize(_gridn + 1);
    
    _solution[_gridn] = new Solution(GetMLMesh()->GetLevel(_gridn));

    // add all current solutions and initialize to zero
    for(unsigned i = 0; i < _solName.size(); i++) {
      _solution[_gridn]->AddSolution(_solName[i], _family[i], _order[i], _solTimeOrder[i], _pdeType[i]);
    }

    /*const*/ Mesh * mesh_current = GetMLMesh()->GetLevel(_gridn);
    
    for(unsigned k = 0; k < _solName.size(); k++) {
      _solution[_gridn]->ResizeSolutionVector(_solName[k]);
      _solution[_gridn]->_Sol[k]->matrix_mult(*_solution[_gridn - 1]->_Sol[k],
                                              * mesh_current->get_prol_matrices().GetCoarseToFineProjection(_solType[k], * mesh_current ));
      _solution[_gridn]->_Sol[k]->close();
      if(_solTimeOrder[k] == TIME_DEPENDENT) {
        _solution[_gridn]->_SolOld[k]->matrix_mult(*_solution[_gridn - 1]->_SolOld[k],
                                                   * mesh_current->get_prol_matrices().GetCoarseToFineProjection(_solType[k], * mesh_current));
        _solution[_gridn]->_SolOld[k]->close();
      }
    }

    _gridn++;

    for(int k = 0; k < _solName.size(); k++) {
      GenerateBdc(k, _gridn - 1, 0.);

    }

  }

//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::AddSolution(const char name[], const FEFamily fefamily, const FEOrder order,
                                       unsigned tmorder, const bool& PdeType) {

    unsigned n = _solType.size();
    _solType.resize(n + 1u);
    _family.resize(n + 1u);
    _order.resize(n + 1u);
    _solName.resize(n + 1u);
    _bdcType.resize(n + 1u);
    _solTimeOrder.resize(n + 1u);
    _pdeType.resize(n + 1u);
    _testIfPressure.resize(n + 1u);
    _addAMRPressureStability.resize(n + 1u);
    _fixSolutionAtOnePoint.resize(n + 1u);
    _solPairIndex.resize(n + 1u);
    _solPairInverseIndex.resize(n + 1u);


    _testIfPressure[n] = 0;
    _addAMRPressureStability[n] = false;
    _fixSolutionAtOnePoint[n] = false;
    _family[n] = fefamily;
    _order[n] = order;
    _solType[n] = order - ((fefamily == LAGRANGE) ? 1 : 0) + fefamily * 3;
    _solName[n]  = new char [SOL_MAXIMUM_NCHARS];
    _bdcType[n]  = new char [20];
    sprintf(_bdcType[n], "undefined");
    strcpy(_solName[n], name);
    _solTimeOrder[n] = tmorder;
    _pdeType[n] = PdeType;
    _solPairIndex[n] = n;
    _solPairInverseIndex[n] = n;


    std::cout << " Add variable " << std::setw(3) << _solName[n] << 
    " discretized with FE Family "  <<   fe_families[ fefamily ] << 
    " and FE order " <<  std::setw(12) << order << 
    " and time discretization order " << tmorder << std::endl;

    for(unsigned ig = 0; ig < _gridn; ig++) {
      _solution[ig]->AddSolution(_solName[n], _family[n], _order[n], _solTimeOrder[n], _pdeType[n]);
    }
  }

 
 
 ///  Weak Galerkin Add Solution -------------------------------------------------------------------------------------
  void MultiLevelSolution::AddSolution(const char name[], const FEFamily fefamily, const FEOrder order_v, const FEOrder order_b,
                                       unsigned tmorder, const bool& PdeType) {


    unsigned n = _solType.size();
    
// ---------------------
//------ first resize -----------
// ---------------------
  

// ID related---
    _solName.resize(n + 1u);
    
// is_an_unknown_of_a_pde---
    _pdeType.resize(n + 1u);
    
//     For pressure variables---
    _testIfPressure.resize(n + 1u);
    _addAMRPressureStability.resize(n + 1u);
//     For pressure variables, fix it at one point for PCFieldSplit, in Boundary Conditions ---
    _fixSolutionAtOnePoint.resize(n + 1u);
    
//     For FSI, pairing velocities and displacements ---
    _solPairIndex.resize(n + 1u);
    _solPairInverseIndex.resize(n + 1u);

//     Time discretization---
    _solTimeOrder.resize(n + 1u);
    
// These are the things that should change for different FE... - BEGIN ------    
//   FE related---  
    _family.resize(n + 1u);
    _order.resize(n + 1u);
    _solType.resize(n + 1u);
    
//     FE related and Time related: Boundary Conditions--- 
    _bdcType.resize(n + 1u);
// These are the things that should change for different FE... - END ------    
    

// ---------------------
//------ then fill -----------
// ---------------------
    
// ID related---
    _solName[n]  = new char [SOL_MAXIMUM_NCHARS];
    strcpy(_solName[n], name);
    
// is_an_unknown_of_a_pde---
    _pdeType[n] = PdeType;
    
//   For pressure variables---
    _testIfPressure[n] = 0;
    _addAMRPressureStability[n] = false;
//   For pressure variables, if solution is fixed at one point (then null space must be removed) ---
    _fixSolutionAtOnePoint[n] = false;
        
//   For FSI, pairing velocities and displacements ---
    _solPairIndex[n] = n;
    _solPairInverseIndex[n] = n;
    
//     Time discretization---
    _solTimeOrder[n] = tmorder;
    
// These are the things that should change for different FE...------    
//   FE related---  
    _family[n] = fefamily;
    _order[n] = order_v;
//     _order_b[n] = order_b;
    _solType[n] = Solution::compute_fe_sol_type(fefamily, order_v, order_b); ///@todo: this normally goes from 0 to 4 for now, I have to understand how to use it
    
//     FE related and Time related: Boundary Conditions
    _bdcType[n]  = new char [20];
    sprintf(_bdcType[n], "undefined");
// These are the things that should change for different FE... - end------    


// ---------------------
    std::cout << " Add variable " << std::setw(3) << _solName[n] <<
    " discretized with FE Family "  <<   fe_families[ fefamily ] << 
    " and FE order " << std::setw(12) << order_v << 
    " and time discretization order " << tmorder << std::endl;    

    for(unsigned ig = 0; ig < _gridn; ig++) {
      _solution[ig]->AddSolution(_solName[n], _family[n], _order[n],  order_b, _solTimeOrder[n], _pdeType[n]);
    }
  }

  
  
  
  void MultiLevelSolution::AddSolutionVector(const unsigned n_components, const std::string name, const FEFamily fefamily, const FEOrder order, unsigned tmorder, const bool& Pde_type) {

    for(unsigned i = 0; i < n_components; i++) {
      std::ostringstream name_cmp;
      name_cmp << name << i;
      AddSolution(name_cmp.str().c_str(), fefamily, order, tmorder, Pde_type);
    }

    return;
  }


  void MultiLevelSolution::ResizeSolution_par(const unsigned new_size)  {

    for(unsigned ig = 0; ig < _solution.size(); ig++) {

      for(unsigned s = 0; s < _solType.size(); s++) {

        _solution[ig]->ResizeSolution_par(new_size);

      }
    }


  }


//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::AssociatePropertyToSolution(const char solution_name[], const char solution_property[],
                                                       const bool& bool_property) {
    unsigned index = GetIndex(solution_name);

    if(!strcmp(solution_property, "pressure") || !strcmp(solution_property, "Pressure")) {
      _testIfPressure[index] = 1;
      _addAMRPressureStability[index] = bool_property;
    }
    else if(!strcmp(solution_property, "default") || !strcmp(solution_property, "Default")) {
      _testIfPressure[index] = 0;
      _addAMRPressureStability[index] == false;
    }
    else {
      std::cout << "Error invalid property in function MultiLevelProblem::AssociatePropertyToSolution" << std::endl;
      exit(0);
    }
  }

// *******************************************************


  void MultiLevelSolution::PairSolution(const char solution_name[], const char solution_pair[]) {
    unsigned index = GetIndex(solution_name);
    unsigned indexPair = GetIndex(solution_pair);
    _solPairIndex[index] = indexPair;
    _solPairInverseIndex[indexPair] = index;
  }

// *******************************************************
  void MultiLevelSolution::Initialize(const char * name, InitFunc func) {
    Initialize(name, func, NULL, NULL);
  }

  void MultiLevelSolution::Initialize(const char * name, InitFuncMLProb func, const MultiLevelProblem* ml_prob) {
    Initialize(name, NULL, func, ml_prob);
  }
  
  
  
  /** A Solution is by default initialized to zero, or by a provided function     */
  void MultiLevelSolution::Initialize(const char * name, InitFunc func, InitFuncMLProb funcMLProb, const MultiLevelProblem* ml_prob) {

    

   std::vector< unsigned > sol_start_end =  solution_start_and_end(std::string (name));
      
  for (unsigned ig = 0; ig < _gridn; ig++) {
      
        const Mesh* msh   = GetMLMesh()->GetLevel(ig);
        const Mesh* msh_2 = _solution[ig]->GetMesh();
        
        if( msh_2 != msh) abort();
        
    for (unsigned i = sol_start_end[0]; i < sol_start_end[1]; i++) {
      const unsigned sol_type = _solType[i];

        
              // initialize vectors - BEGIN -----
        _solution[ig]->ResizeSolutionVector(_solName[i]);
        _solution[ig]->_Sol[i]->zero();
              // initialize vectors - END -----

        if(func || funcMLProb) {


          if(sol_type < NFE_FAMS_C_ZERO_LAGRANGE) {

            for(int isdom = _iproc; isdom < _iproc + 1; isdom++) {

              for(int iel = msh->_elementOffset[isdom];
                  iel < msh->_elementOffset[isdom + 1]; iel++) {
                
                const unsigned nloc_dof = msh->GetElementDofNumber(iel, sol_type);

                for(int j = 0; j < nloc_dof; j++) {

                  const std::vector< double > xx = evaluate_coord_of_dof_carrier_Lagrange(msh, iel, j);

                  const double value = evaluate_at_dof_carrier(name, func, funcMLProb, NULL, ml_prob, xx);
                  
                  const unsigned inode_Metis = msh->GetSolutionDof(j, iel, sol_type);
                  
                  _solution[ig]->_Sol[i]->set(inode_Metis, value);

                  if(_solTimeOrder[i] == TIME_DEPENDENT) {
                    _solution[ig]->_SolOld[i]->set(inode_Metis, value);
                  }
                  
                }
                
              }
              
            }

          }
          else if(sol_type < NFE_FAMS) {

            for(int isdom = _iproc; isdom < _iproc + 1; isdom++) {

              for(int iel = msh->_elementOffset[isdom];
                  iel < msh->_elementOffset[isdom + 1]; iel++) {
                  
                const std::vector < double > xx = evaluate_coord_of_dof_carrier_Discontinuous(msh, iel);

                const double value = evaluate_at_dof_carrier(name, func, funcMLProb, NULL, ml_prob, xx);
                
                const unsigned placeholder_index = 0; ///@todo see how to interpolate disc linear with first derivatives

                const unsigned solDof = msh->GetSolutionDof(placeholder_index, iel, sol_type);

                _solution[ig]->_Sol[i]->set(solDof, value);

                if(_solTimeOrder[i] == TIME_DEPENDENT) {
                  _solution[ig]->_SolOld[i]->set(solDof, value);
                }
                
              }
            }
            
          }

              // close vectors - BEGIN -----
          _solution[ig]->_Sol[i]->close();

          if(_solTimeOrder[i] == TIME_DEPENDENT) {
            _solution[ig]->_SolOld[i]->close();
          }
              // close vectors - END -----

        }
        
      } //end solutions
    
  } //end levels

    return;
  }



  double MultiLevelSolution::evaluate_at_dof_carrier(const char * name,
                                                     InitFunc func, 
                                                     InitFuncMLProb funcMLProb, 
                                                     const Math::Function< double > * func_in,
                                                     const MultiLevelProblem * ml_prob,
                                                     const std::vector<double> xx) const {

                  const double value = (func) ? func(xx) : funcMLProb(ml_prob, xx, name);
                  
                  return value;
  }


  
  std::vector< double > MultiLevelSolution::evaluate_coord_of_dof_carrier_Lagrange(const Mesh * msh, const unsigned iel, const unsigned inode) const {
    
                  const unsigned icoord_Metis = msh->GetSolutionDof(inode, iel, CONTINUOUS_BIQUADRATIC);
                  
                  std::vector < double > xx(3);

                  xx[0] = (*msh->GetTopology()->_Sol[0])(icoord_Metis);
                  xx[1] = (*msh->GetTopology()->_Sol[1])(icoord_Metis);
                  xx[2] = (*msh->GetTopology()->_Sol[2])(icoord_Metis);

                  return xx;
    
  }
  
  
  std::vector< double > MultiLevelSolution::evaluate_coord_of_dof_carrier_Discontinuous(const Mesh * msh, const unsigned iel) const {

                const unsigned nloc_dof = msh->GetElementDofNumber(iel, CONTINUOUS_BIQUADRATIC);
              
                std::vector < double > xx(3, 0.);

                for(int j = 0; j < nloc_dof; j++) {
                  const unsigned icoord_Metis = msh->GetSolutionDof(j, iel, CONTINUOUS_BIQUADRATIC);
                  xx[0] += (*msh->GetTopology()->_Sol[0])(icoord_Metis);
                  xx[1] += (*msh->GetTopology()->_Sol[1])(icoord_Metis);
                  xx[2] += (*msh->GetTopology()->_Sol[2])(icoord_Metis);
                }

                xx[0] /= nloc_dof;
                xx[1] /= nloc_dof;
                xx[2] /= nloc_dof;

                  return xx;
    
  }
  
  
  std::vector< unsigned >  MultiLevelSolution::solution_start_and_end(const std::string name) const {

      std::vector< unsigned > sol_start_end(2, 0);

    if(!strcmp(name.c_str(), "All") || !strcmp(name.c_str(), "all") || !strcmp(name.c_str(), "ALL")) {
      sol_start_end[0] = 0;
      sol_start_end[1] = _solType.size();
    }
    else {
      sol_start_end[0] = GetIndex(name.c_str());
      sol_start_end[1] = sol_start_end[0] + 1u;
    }
    
    return sol_start_end;
    
  }
  
  
//---------------------------------------------------------------------------------------------------
  unsigned MultiLevelSolution::GetIndex(const char name[]) const {
    
    return _solution[_sol_level_to_pick_from]->GetIndex(name);
    
  }

// *******************************************************
  unsigned MultiLevelSolution::GetSolutionType(const char name[]) const {

    return _solution[_sol_level_to_pick_from]->GetSolutionType(name);
  }



//---------------------------------------------------------------------------------------------------
  const FEFamily MultiLevelSolution::GetSolutionFamily(const std::string & sol_name) const {

    const char * name = sol_name.c_str();

    const unsigned int index = GetIndex(name);

    return _family[index];

  }

//---------------------------------------------------------------------------------------------------
  const FEOrder MultiLevelSolution::GetSolutionOrder(const std::string & sol_name) const {

    const char * name = sol_name.c_str();

    const unsigned index = GetIndex(name);

    return _order[index];

  }


//---------------------------------------------------------------------------------------------------
  const int MultiLevelSolution::GetSolutionTimeOrder(const std::string & sol_name) const {

    const char * name = sol_name.c_str();

    const unsigned index = GetIndex(name);

    return _solTimeOrder[index];

  }


//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::AttachSetBoundaryConditionFunction(BoundaryFuncMLProb SetBoundaryConditionFunction_in) {
    _bdcFuncSetMLProb = true;
    _bdcFuncSet = false;
    _SetBoundaryConditionFunctionMLProb = SetBoundaryConditionFunction_in;
    return;
  }

//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::AttachSetBoundaryConditionFunction(BoundaryFunc SetBoundaryConditionFunction) {
    _bdcFuncSet = true;
    _bdcFuncSetMLProb = false;
    _SetBoundaryConditionFunction = SetBoundaryConditionFunction;
    return;
  }

//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::InitializeBdc() {
    _useParsedBCFunction = true;

    int nvars = _solType.size();
    int nfaces = GetMLMesh()->GetLevel(0)->GetBoundaryInfo().size();
    _boundaryConditions.resize(nvars);
    _isHomogeneous.resize(nvars);
    _nonHomogeneousBCFunction.resize(nvars);

    for(int i = 0; i < nvars; i++) {
      _boundaryConditions[i].resize(nfaces);
      _isHomogeneous[i].resize(nfaces);
      _nonHomogeneousBCFunction[i].resize(nfaces);

      for(int j = 0; j < nfaces; j++) {
        _boundaryConditions[i][j] = DIRICHLET;
        _isHomogeneous[i][j] = true;
        _nonHomogeneousBCFunction[i][j] = NULL;
      }
    }
  }

//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::SetBoundaryCondition_new(const std::string name, const std::string facename,
                                                    const BDCType bdctype, const bool istimedependent, FunctionBase* func) {



    unsigned int ivar = GetIndex(name.c_str());
    unsigned int iface = 0;
    bool ishomogeneous = true;

    if(func != NULL) {
      ishomogeneous = false;
    }

    std::map<unsigned int, std::string>::iterator iter;
    iter = GetMLMesh()->GetLevel(0)->GetBoundaryInfo().begin();

    for(iter = GetMLMesh()->GetLevel(0)->GetBoundaryInfo().begin(); iter != GetMLMesh()->GetLevel(0)->GetBoundaryInfo().end(); ++iter) {
      if(iter->second.compare(facename) == 0) {
        iface = iter->first;
        break;
      }
    }

    if(iter == GetMLMesh()->GetLevel(0)->GetBoundaryInfo().end()) {
      std::cout << " Error: the facename " << facename << " does not exist!" << std::endl;
      exit(1);
    }

    _boundaryConditions[ivar][iface]       = bdctype;
    _isHomogeneous[ivar][iface]            = ishomogeneous;
    _nonHomogeneousBCFunction[ivar][iface] = func;

  }







//---------------------------------------------------------------------------------------------------
 /** bdc_type can be "Steady" or "Time_dependent" */
void MultiLevelSolution::GenerateBdc(const char* name, const char* bdc_type, const MultiLevelProblem* ml_prob) {

    if(_useParsedBCFunction == false && _bdcFuncSet == false && _bdcFuncSetMLProb == false) {
      std::cout << "Error: The boundary condition user-function is not set! Please call the AttachSetBoundaryConditionFunction routine"
           << std::endl;
      abort();
    }

    if(_bdcFuncSetMLProb == true && ml_prob != NULL) {
      _mlBCProblem = ml_prob;
    }
    else if(_bdcFuncSetMLProb == true && ml_prob == NULL) {
      std::cout << "Warning, no MultilevelProblem pointer has been passed in MultiLevelSolution::GenerateBdc(...}!!!" << std::endl;
      std::cout << "The MultiLevelProblem pointer in the boundary function will be unavailable, and if used will cause segmentation-fault." << std::endl;
      std::cout << "It is highly recommended either to pass the MultilevelProblem pointer in MultiLevelSolution::GenerateBdc(...} or to use" << std::endl;
      std::cout << "MultiLevelSolution::AttachSetBoundaryConditionFunction( (*BoundaryFunc) ) rather than" << std::endl;
      std::cout << "MultiLevelSolution::AttachSetBoundaryConditionFunction( (*BoundaryFuncMLProb) )" << std::endl;
    }

    unsigned i_start;
    unsigned i_end;

    if(!strcmp(name, "All")  || !strcmp(name, "all") || !strcmp(name, "ALL")) {

      i_start = 0;
      i_end = _solType.size();

      for(unsigned k = i_start; k < i_end; k++) {
        if(_solution[0]->is_unknown_of_system(k) ) {
          sprintf(_bdcType[k], "Steady");        /// @todo since bdc_type is a single string, why do we assume that All variables must be steady? They could also be all Time dependent...
          std::cout << " Set " << std::setw(15) << _bdcType[k] << " Boundary_condition"
               << " for variable " << std::setw(3) << _solName[k] << std::endl;
        }
        else {
          sprintf(_bdcType[k], "Not-available");
        }

      }
    }
    else {
      i_start = GetIndex(name);
      i_end = i_start + 1u;

      if(_solution[0]->is_unknown_of_system(i_start) ) {
        if(!strcmp(bdc_type, "Steady")) {
          strcpy(_bdcType[i_start], bdc_type);
        }
        else if(!strcmp(bdc_type, "Time_dependent")) {
          strcpy(_bdcType[i_start], bdc_type);
        }
        else {
          std::cout << "Error! Invalid boundary condition specified for " << _solName[i_start]
               << " in GenerateBdc function" << std::endl;
          exit(1);
        }

        std::cout << " Set " << std::setw(15) << _bdcType[i_start] << " Boundary_condition"
             << " for variable " << std::setw(3) << _solName[i_start] << std::endl;
      }
      else {
        sprintf(_bdcType[i_start], "Not-available");
      }
    }

    for(unsigned i = i_start; i < i_end; i++) {
      GenerateBdc(i, 0, 0.);
    }
  }

//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::UpdateBdc(const double time) {

    for(int k = 0; k < _solName.size(); k++) {
      if(!strcmp(_bdcType[k], "Time_dependent")) {
        GenerateBdc(k, 0, time);
      }
    }
  }


//---------------------------------------------------------------------------------------------------
 /**  2 Default Neumann
   *  1 AMR artificial Dirichlet = 0 BC
   *  0 Dirichlet
  */  
  void MultiLevelSolution::GenerateBdc(const unsigned int k, const unsigned int grid0, const double time) {

    for(unsigned igridn = grid0; igridn < _gridn; igridn++) {
      
      if(_solution[igridn]->is_unknown_of_system(k) ) {
        
        const Mesh* msh = GetMLMesh()->GetLevel(igridn);

        const std::vector < std::map < unsigned,  std::map < unsigned, double  > > > & amrRestriction = msh->GetAmrRestrictionMap();

        // default Neumann
        for(unsigned j = msh->dofmap_get_dof_offset(_solType[k], _iproc); j < msh->dofmap_get_dof_offset(_solType[k], _iproc + 1); j++) {
          _solution[igridn]->_Bdc[k]->set(j, 2.);
        }

        if(_solType[k] < NFE_FAMS_C_ZERO_LAGRANGE) {  // boundary condition for lagrangian elements


          //AMR - related - BEGIN
          for(int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {
            for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
              if(msh->el->GetBoundaryIndex(iel, jface) == 0) {   // interior boundary (AMR) u = 0
                short unsigned ielt = msh->GetElementType(iel);
                unsigned nv1 = msh->GetElementFaceDofNumber(iel, jface, _solType[k]);  // only the face dofs
                for(unsigned iv = 0; iv < nv1; iv++) {
                  unsigned i = msh->GetLocalFaceVertexIndex(iel, jface, iv);
                  unsigned idof = msh->GetSolutionDof(i, iel, _solType[k]);
                  if(amrRestriction[_solType[k]].find(idof) != amrRestriction[_solType[k]].end() &&
                      amrRestriction[_solType[k]].at(idof).at(idof) == 0) {
                    _solution[igridn]->_Bdc[k]->set(idof, 1.);
                  }
                }
              }
            }
          }
          //AMR - related - END
          
          for(int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {

            for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {

              const int boundary_index = msh->el->GetBoundaryIndex(iel, jface);

              if(boundary_index > 0) {   // exterior boundary u = value
                short unsigned ielt = msh->GetElementType(iel);
                unsigned n_face_dofs = msh->GetElementFaceDofNumber(iel, jface, _solType[k]);

                for(unsigned iv = 0; iv < n_face_dofs; iv++) {

                  unsigned i = msh->GetLocalFaceVertexIndex(iel, jface, iv);

                  
                  if(_useParsedBCFunction) {
                    unsigned int faceIndex = boundary_index;

                    if(GetBoundaryCondition(k, faceIndex - 1u) == DIRICHLET) {
                      
                      const unsigned inode_Metis = msh->GetSolutionDof(i, iel, _solType[k]);
                      _solution[igridn]->_Bdc[k]->set(inode_Metis, 0.);
                      double value = 0.;

                      if(!Ishomogeneous(k, faceIndex - 1u)) {
                        
                        ParsedFunction* bdcfunc = (ParsedFunction*)(GetBdcFunction(k, faceIndex - 1u));

                       const std::vector< double > xx = evaluate_coord_of_dof_carrier_Lagrange(msh, iel, i);

                        double xyzt[4];
                        xyzt[0] = xx[0];
                        xyzt[1] = xx[1];
                        xyzt[2] = xx[2];
                        xyzt[3] = time;
                        
                        value = (*bdcfunc)(xyzt);
                        
                      }
                      _solution[igridn]->_Sol[k]->set(inode_Metis, value);
                    }
                  }
                  else {
                    double value;
                    
                  const std::vector< double > xx = evaluate_coord_of_dof_carrier_Lagrange(msh, iel, i);

                  bool test = (_bdcFuncSetMLProb) ?

                                _SetBoundaryConditionFunctionMLProb(_mlBCProblem, xx, _solName[k], value, boundary_index, time) :
                                _SetBoundaryConditionFunction(xx, _solName[k], value, boundary_index, time);

                    if(test) {
                      unsigned idof = msh->GetSolutionDof(i, iel, _solType[k]);
                      _solution[igridn]->_Bdc[k]->set(idof, 0.);
                      _solution[igridn]->_Sol[k]->set(idof, value);
                    }
                  }
                }

              } //end boundary faces
            }  //end faces
          }  //end element
        }  //end Lagrangian


        if(_fixSolutionAtOnePoint[k] == true  && igridn == 0 && _iproc == 0) {
          _solution[igridn]->_Bdc[k]->set(0, 0.);
          _solution[igridn]->_Sol[k]->set(0, 0.);


        }
        _solution[igridn]->_Sol[k]->close();
        _solution[igridn]->_Bdc[k]->close();
      }
    }


  }


  void MultiLevelSolution::GenerateBdcOnVolumeConstraint(const std::vector<unsigned> &volumeConstraintFlags, const unsigned &solIndex, const unsigned &grid0) {

    unsigned solType = GetSolutionType(solIndex);

    for(unsigned igridn = grid0; igridn < _gridn; igridn++) {
      Mesh* msh = GetMLMesh()->GetLevel(igridn);
      for(int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {

        short unsigned ielGroup = msh->GetElementGroup(iel);

        bool ielIsInVolumeConstraint = false;

        for(unsigned i = 0; i < volumeConstraintFlags.size(); i++) {
          if(volumeConstraintFlags[i] == ielGroup) {
            ielIsInVolumeConstraint = true;
            break;
          }
        }

        if(ielIsInVolumeConstraint) {

          unsigned nDofu  = msh->GetElementDofNumber(iel, solType);

          for (unsigned i = 0; i < nDofu; i++) {
            unsigned solDof = msh->GetSolutionDof (i, iel, solType);
            _solution[igridn]->_Bdc[solIndex]->set (solDof, 0.);


            double value;
            
            const std::vector< double > xCord = evaluate_coord_of_dof_carrier_Lagrange(msh, iel, i);
                   
            bool test = (_bdcFuncSetMLProb) ?
                        _SetBoundaryConditionFunctionMLProb (_mlBCProblem, xCord, _solName[solIndex], value, UINT_MAX, 0.) :
                        _SetBoundaryConditionFunction (xCord, _solName[solIndex], value, UINT_MAX, 0.);

            if (test) {
              _solution[igridn]->_Sol[solIndex]->set (solDof, value);
            }

          }

        }
      }

      _solution[igridn]->_Bdc[solIndex]->close();
      _solution[igridn]->_Sol[solIndex]->close();

    }

  }
  
  
  
  

  bool PRINT = false;
  void MultiLevelSolution::GenerateRKBdc(
    const unsigned int &solIndex, const std::vector<unsigned> &solKiIndex,
    const unsigned int &grid0, const std::vector < double> & itime, const double &time0, const double &dt, const double *AI) {

    for(unsigned k = 0; k < solKiIndex.size(); k++) {
      strcpy(_bdcType[solKiIndex[k]], "RKbdc");
    }
    // 2 Default Neumann
    // 1 AMR artificial Dirichlet = 0 BC
    // 0 Dirichlet
    for(unsigned igridn = grid0; igridn < _gridn; igridn++) {
      if(_solution[igridn]->is_unknown_of_system(solIndex) ) {
        const Mesh* msh = GetMLMesh()->GetLevel(igridn);

        const std::vector < std::map < unsigned,  std::map < unsigned, double  > > > & amrRestriction = msh->GetAmrRestrictionMap();

        // default Neumann
        for(unsigned j = msh->dofmap_get_dof_offset(_solType[solIndex], _iproc); j < msh->dofmap_get_dof_offset(_solType[solIndex], _iproc + 1); j++) {
          for(unsigned k = 0; k < solKiIndex.size(); k++) {
            _solution[igridn]->_Bdc[solKiIndex[k]]->set(j, 2.);
          }
        }
        

        if(_solType[solIndex] < NFE_FAMS_C_ZERO_LAGRANGE) {  // boundary condition for lagrangian elements
          
          
          for(int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {
            for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
              if(msh->el->GetBoundaryIndex(iel, jface) == 0) {   // interior boundary (AMR) u = 0
                short unsigned ielt = msh->GetElementType(iel);
                unsigned nv1 = msh->GetElementFaceDofNumber(iel, jface, _solType[solIndex]);  // only the face dofs
                for(unsigned iv = 0; iv < nv1; iv++) {
                  unsigned i = msh->GetLocalFaceVertexIndex(iel, jface, iv);
                  unsigned idof = msh->GetSolutionDof(i, iel, _solType[solIndex]);
                  if(amrRestriction[_solType[solIndex]].find(idof) != amrRestriction[_solType[solIndex]].end() &&
                      amrRestriction[_solType[solIndex]].at(idof).at(idof) == 0) {
                    for(unsigned k = 0; k < solKiIndex.size(); k++) {
                      _solution[igridn]->_Bdc[solKiIndex[k]]->set(idof, 1.);
                    }
                  }
                }
              }
            }
          }

          for(int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {
            for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
              if(msh->el->GetBoundaryIndex(iel, jface) > 0) {   // exterior boundary u = value
                short unsigned ielt = msh->GetElementType(iel);
                unsigned nv1 = msh->GetElementFaceDofNumber(iel, jface, _solType[solIndex]);

                for(unsigned iv = 0; iv < nv1; iv++) {
                  const unsigned i = msh->GetLocalFaceVertexIndex(iel, jface, iv);
                  const unsigned inode_coord_Metis = msh->GetSolutionDof(i, iel, CONTINUOUS_BIQUADRATIC);
                  
                  if(_useParsedBCFunction) {
                    unsigned int faceIndex = msh->el->GetBoundaryIndex(iel, jface);

                    if(GetBoundaryCondition(solIndex, faceIndex - 1u) == DIRICHLET) {
                      const unsigned inode_Metis = msh->GetSolutionDof(i, iel, _solType[solIndex]);
                      for(unsigned k = 0; k < solKiIndex.size(); k++) {
                        _solution[igridn]->_Bdc[solKiIndex[k]]->set(inode_Metis, 0.);
                      }
                      std::vector < double > ivalue(solKiIndex.size(), 0.);

                      if(!Ishomogeneous(solIndex, faceIndex - 1u)) {

                        ParsedFunction* bdcfunc = (ParsedFunction*)(GetBdcFunction(solIndex, faceIndex - 1u));

                       const std::vector< double > xx = evaluate_coord_of_dof_carrier_Lagrange(msh, iel, i);

                       double xyzt[4];
                        xyzt[0] = xx[0];
                        xyzt[1] = xx[1];
                        xyzt[2] = xx[2];
                        xyzt[3] = time0;
                        
                        double value0 = (*bdcfunc)(xyzt);

                        for(unsigned k = 0; k < solKiIndex.size(); k++) {
                          xyzt[3] = itime[k];
                          ivalue[k] = (*bdcfunc)(xyzt) - value0;
                        }
                        if(_solTimeOrder[solIndex] == TIME_DEPENDENT) {
                          _solution[igridn]->_SolOld[solIndex]->set(inode_Metis, value0);
                        }
                        for(unsigned k1 = 0; k1 < solKiIndex.size(); k1++) {
                          double value = 0.;
                          for(unsigned k2 = 0; k2 < solKiIndex.size(); k2++) {
                            value += AI[k1 * solKiIndex.size() + k2] * ivalue[k2];
                          }
                          _solution[igridn]->_Sol[solKiIndex[k1]]->set(inode_Metis, value / dt);
                        }
                      }
                    }
                  }
                  else {
                    double value0;
                    
                    const std::vector< double > xx = evaluate_coord_of_dof_carrier_Lagrange(msh, iel, i);
                  
                    bool test = (_bdcFuncSetMLProb) ?
                                _SetBoundaryConditionFunctionMLProb(_mlBCProblem, xx, _solName[solIndex], value0, msh->el->GetBoundaryIndex(iel, jface), time0) :
                                _SetBoundaryConditionFunction(xx, _solName[solIndex], value0, msh->el->GetBoundaryIndex(iel, jface), time0);

                    if(test) {

                      if(PRINT) std::cout << inode_coord_Metis << " " << value0 << "\n" ;

                      std::vector < double > ivalue(solKiIndex.size(), 0.);
                      for(unsigned k = 0; k < solKiIndex.size(); k++) {
                        if(_bdcFuncSetMLProb) {
                          _SetBoundaryConditionFunctionMLProb(_mlBCProblem, xx, _solName[solIndex], ivalue[k], msh->el->GetBoundaryIndex(iel, jface), itime[k]);
                        }
                        else {
                          _SetBoundaryConditionFunction(xx, _solName[solIndex], ivalue[k], msh->el->GetBoundaryIndex(iel, jface), itime[k]);
                        }
                        if(PRINT) std::cout << inode_coord_Metis << " " << ivalue[k] << std::endl;
                        ivalue[k] -= value0;
                      }

                      unsigned idof = msh->GetSolutionDof(i, iel, _solType[solIndex]);

                      if(_solTimeOrder[solIndex] == TIME_DEPENDENT) {
                        _solution[igridn]->_SolOld[solIndex]->set(idof, value0);
                      }
                      for(unsigned k1 = 0; k1 < solKiIndex.size(); k1++) {
                        double value = 0.;
                        for(unsigned k2 = 0; k2 < solKiIndex.size(); k2++) {
                          value += AI[k1 * solKiIndex.size() + k2] * ivalue[k2];
                        }
                        _solution[igridn]->_Bdc[solKiIndex[k1]]->set(idof, 0.);
                        _solution[igridn]->_Sol[solKiIndex[k1]]->set(idof, value / dt);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        
        if(_fixSolutionAtOnePoint[solIndex] == true  && igridn == 0 && _iproc == 0) {
          for(unsigned k = 0; k < solKiIndex.size(); k++) {

            _solution[igridn]->_Bdc[solKiIndex[k]]->set(0, 0.);
            _solution[igridn]->_Sol[solKiIndex[k]]->set(0, 0.);

          }
        }
        
        
        if(_solTimeOrder[solIndex] == TIME_DEPENDENT) {
          _solution[igridn]->_SolOld[solIndex]->close();
        }
        
        for(unsigned k = 0; k < solKiIndex.size(); k++) {
          _solution[igridn]->_Sol[solKiIndex[k]]->close();
          _solution[igridn]->_Bdc[solKiIndex[k]]->close();

        }
      }
    }

    PRINT = false;
  }



  void MultiLevelSolution::SaveSolution(const char* filename, const double &time) {

    char composedFileName[100];

    for(int i = 0; i < _solName.size(); i++) {
      sprintf(composedFileName, "./save/%s_time%f_sol%s_level%d", filename, time, _solName[i], _gridn);
      _solution[_gridn - 1]->_Sol[i]->BinaryPrint(composedFileName);
    }
  }

  void MultiLevelSolution::SaveSolution(const char* filename, const unsigned &iteration) {

    char composedFileName[100];

    for(int i = 0; i < _solName.size(); i++) {
      sprintf(composedFileName, "./save/%s_iteration%d_sol%s_level%d", filename, iteration, _solName[i], _gridn);
      _solution[_gridn - 1]->_Sol[i]->BinaryPrint(composedFileName);
    }
  }

  void MultiLevelSolution::LoadSolution(const char* filename) {
    LoadSolution(_gridn, filename);
  }

  void MultiLevelSolution::LoadSolution(const unsigned &level, const char* filename) {

    if(level > _gridn) {
      std::cout << "Error in MultiLevelSolution::LoadSolution function:" << std::endl;
      std::cout << "the solution level = " << level << " is not available in this MultilevelSolution" << std::endl;
      abort();
    }

    char composedFileName[200];
    for(int i = 0; i < _solName.size(); i++) {
      sprintf(composedFileName, "%s_sol%s_level%d", filename, _solName[i], level);
      // check if the file really exists
      if(strncmp(filename, "http://", 7) && strncmp(filename, "ftp://", 6)) {
        struct stat buffer;
        if(stat(composedFileName, &buffer) != 0) {
          std::cerr << "Error: cannot locate file " << composedFileName << std::endl;
          abort();
        }
      }
      _solution[level - 1]->_Sol[i]->BinaryLoad(composedFileName);
    }

    for(int gridf = level; gridf < _gridn; gridf++) {
      Mesh * meshf = GetMLMesh()->GetLevel(gridf);
      for(unsigned i = 0; i < _solName.size(); i++) {
        _solution[gridf]->_Sol[i]->matrix_mult(*_solution[gridf - 1]->_Sol[i],
                                               * meshf->get_prol_matrices().GetCoarseToFineProjection( _solType[i], * meshf) );
        _solution[gridf]->_Sol[i]->close();
      }
    }


  }




  /** Refine the solution at (gridf) level from (gridf - 1) */
  void MultiLevelSolution::RefineSolution(const unsigned &gridf) {

    Mesh *msh = GetMLMesh()->GetLevel(gridf);

    for(unsigned k = 0; k < _solType.size(); k++) {

      unsigned solType = _solType[k];
      _solution[gridf]->_Sol[k]->matrix_mult(*_solution[gridf - 1]->_Sol[k],
                                             * msh->get_prol_matrices().GetCoarseToFineProjection(solType, *msh));
      _solution[gridf]->_Sol[k]->close();
    }
  }




  void MultiLevelSolution::CoarsenSolutionByOneLevel_wrong(const unsigned &grid_fine)  {


    Mesh *msh = GetMLMesh()->GetLevel(grid_fine);
    const unsigned grid_coarse = grid_fine - 1;

    for(unsigned k = 0; k < _solType.size(); k++) {

      unsigned solType = _solType[k];
      _solution[grid_coarse]->_Sol[k]->matrix_mult_transpose(*(_solution[grid_fine]->_Sol[k]), *(msh->get_prol_matrices().GetCoarseToFineProjectionRestrictionOnCoarse(solType, *msh)));
      _solution[grid_coarse]->_Sol[k]->close();
    }

  }

  void MultiLevelSolution::CoarsenSolutionByOneLevel(const unsigned &grid_fine)  {

    const unsigned grid_coarse = grid_fine - 1;
    Mesh *msh = GetMLMesh()->GetLevel(grid_coarse);

    //loop over the coarse elements
    //loop over the dofs of each coarse element
    //loop over the children elements of that coarse element
    //find the child element to which the coarse dof belongs, and find the child dof
    //set the _Sol in the coarse dofs to the fine dofs

  }

  /** Copies from another MLSol object from a given level to some other level.
      One should also check that they belong to the same underlying mesh structure, have the same list of variables, and in the same order */
  void MultiLevelSolution::fill_at_level_from_level(const unsigned lev_out, const unsigned lev_in, const MultiLevelSolution & ml_sol_in)  {

    if (_solType.size() != ml_sol_in.GetSolutionSize()) { std::cout << "Different Solutions" << std::endl; abort(); }

    for(unsigned k = 0; k < _solType.size(); k++) {
      *(_solution[lev_out]->_Sol[k]) = *(ml_sol_in.GetSolutionLevel(lev_in)->_Sol[k]);
    }

  }


  void MultiLevelSolution::add_solution(const unsigned index_read, const unsigned index_write) {

    for(unsigned short i = 0; i < _gridn; i++) {
      _solution[i]->add_solution(index_read, index_write);
    }
  }


  void MultiLevelSolution::CopySolutionToOldSolution() {

    for(unsigned short i = 0; i < _gridn; i++) {
      _solution[i]->CopySolutionToOldSolution();
    }
  }

  void MultiLevelSolution::UpdateSolution(const char * name, InitFunc func, const double& time) {

    unsigned i = GetIndex(name);

    unsigned sol_type = _solType[i];

    
    for(unsigned ig = 0; ig < _gridn; ig++) {
      
        const Mesh* msh = GetMLMesh()->GetLevel(ig);

      if(sol_type < NFE_FAMS_C_ZERO_LAGRANGE) {

        for(int isdom = _iproc; isdom < _iproc + 1; isdom++) {
          for(int iel = msh->_elementOffset[isdom];
              iel < msh->_elementOffset[isdom + 1]; iel++) {
            
            unsigned nloc_dof = msh->GetElementDofNumber(iel, sol_type);

            for(int j = 0; j < nloc_dof; j++) {
              
              const std::vector< double > xx = evaluate_coord_of_dof_carrier_Lagrange(msh, iel, j);

              std::vector < double > xyzt(4);
              
              xyzt[0] = xx[0];
              xyzt[1] = xx[1];
              xyzt[2] = xx[2];
              xyzt[3] = time;
              
              const double value = func(xyzt);
              
              const unsigned inode_Metis = msh->GetSolutionDof(j, iel, sol_type);
              
              _solution[ig]->_Sol[i]->set(inode_Metis, value);
              
            }
          }
        }
      }
      
      else if(sol_type < NFE_FAMS) {
        
        for(int isdom = _iproc; isdom < _iproc + 1; isdom++) {
          for(int iel = msh->_elementOffset[isdom];
              iel < msh->_elementOffset[isdom + 1]; iel++) {
          
            const std::vector < double > xx = evaluate_coord_of_dof_carrier_Discontinuous(msh, iel);

            std::vector < double > xyzt(4, 0.);
          
            xyzt[0] = xx[0];
            xyzt[1] = xx[1];
            xyzt[2] = xx[2];
            xyzt[3] = time;

            const double value =  func(xyzt);
            
            const unsigned placeholder_index = 0; ///@todo see how to interpolate disc linear with first derivatives
            
            const unsigned solDof = msh->GetSolutionDof(placeholder_index, iel, sol_type);
            _solution[ig]->_Sol[i]->set(solDof, value);
          }
        }
      }
      _solution[ig]->_Sol[i]->close();
    }
    return;
  }


  

  

} //end namespace femus




