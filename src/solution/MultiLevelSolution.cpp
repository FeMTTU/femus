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
#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "FemusConfig.hpp"
#include "FemusDefault.hpp"
#include "ParsedFunction.hpp"



//C++ include
#include <iostream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>

namespace femus {


  using std::cout;
  using std::endl;

//---------------------------------------------------------------------------------------------------
  MultiLevelSolution::~MultiLevelSolution()  {

    clear();

  }

//---------------------------------------------------------------------------------------------------
  MultiLevelSolution::MultiLevelSolution (MultiLevelMesh* ml_msh) :
    _gridn (ml_msh->GetNumberOfLevels()),
    _mlMesh (ml_msh) {
    _solution.resize (_gridn);

    for (unsigned i = 0; i < _gridn; i++) {
      _solution[i] = new Solution (_mlMesh->GetLevel (i));
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
    _solution.resize (_gridn + 1);
    _solution[_gridn] = new Solution (_mlMesh->GetLevel (_gridn));

    // add all current solutions and initialize to zero
    for (unsigned i = 0; i < _solName.size(); i++) {
      _solution[_gridn]->AddSolution (_solName[i], _family[i], _order[i], _solTimeOrder[i], _pdeType[i]);
    }

    for (unsigned k = 0; k < _solName.size(); k++) {
      _solution[_gridn]->ResizeSolutionVector (_solName[k]);
      _solution[_gridn]->_Sol[k]->matrix_mult (*_solution[_gridn - 1]->_Sol[k],
                                               *_mlMesh->GetLevel (_gridn)->GetCoarseToFineProjection (_solType[k]));
      _solution[_gridn]->_Sol[k]->close();
      if (_solTimeOrder[k] == 2) {
        _solution[_gridn]->_SolOld[k]->matrix_mult (*_solution[_gridn - 1]->_SolOld[k],
                                                    *_mlMesh->GetLevel (_gridn)->GetCoarseToFineProjection (_solType[k]));
        _solution[_gridn]->_SolOld[k]->close();
      }
    }

    _gridn++;

    for (int k = 0; k < _solName.size(); k++) {
      GenerateBdc (k, _gridn - 1, 0.);

    }

  }

//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::AddSolution (const char name[], const FEFamily fefamily, const FEOrder order,
                                        unsigned tmorder, const bool& PdeType) {

    unsigned n = _solType.size();
    _solType.resize (n + 1u);
    _family.resize (n + 1u);
    _order.resize (n + 1u);
    _solName.resize (n + 1u);
    _bdcType.resize (n + 1u);
    _solTimeOrder.resize (n + 1u);
    _pdeType.resize (n + 1u);
    _testIfPressure.resize (n + 1u);
    _addAMRPressureStability.resize (n + 1u);
    _fixSolutionAtOnePoint.resize (n + 1u);
    _solPairIndex.resize (n + 1u);
    _solPairInverseIndex.resize (n + 1u);


    _testIfPressure[n] = 0;
    _addAMRPressureStability[n] = false;
    _fixSolutionAtOnePoint[n] = false;
    _family[n] = fefamily;
    _order[n] = order;
    _solType[n] = order - ( (fefamily == LAGRANGE) ? 1 : 0) + fefamily * 3;
    _solName[n]  = new char [DEFAULT_SOL_NCHARS];
    _bdcType[n]  = new char [20];
    sprintf (_bdcType[n], "undefined");
    strcpy (_solName[n], name);
    _solTimeOrder[n] = tmorder;
    _pdeType[n] = PdeType;
    _solPairIndex[n] = n;
    _solPairInverseIndex[n] = n;


    cout << " Add variable " << std::setw (3) << _solName[n] << " discretized with FE type "
         << std::setw (12) << order << " and time discretzation order " << tmorder << endl;

    for (unsigned ig = 0; ig < _gridn; ig++) {
      _solution[ig]->AddSolution (_solName[n], _family[n], _order[n], _solTimeOrder[n], _pdeType[n]);
    }
  }

  void MultiLevelSolution::AddSolutionVector (const unsigned n_components, const std::string name, const FEFamily fefamily, const FEOrder order, unsigned tmorder, const bool& Pde_type) {

    for (unsigned i = 0; i < n_components; i++) {
      std::ostringstream name_cmp;
      name_cmp << name << i;
      AddSolution (name_cmp.str().c_str(), fefamily, order, tmorder, Pde_type);
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
  void MultiLevelSolution::AssociatePropertyToSolution (const char solution_name[], const char solution_property[],
                                                        const bool& bool_property) {
    unsigned index = GetIndex (solution_name);

    if (!strcmp (solution_property, "pressure") || !strcmp (solution_property, "Pressure")) {
      _testIfPressure[index] = 1;
      _addAMRPressureStability[index] = bool_property;
    }
    else if (!strcmp (solution_property, "default") || !strcmp (solution_property, "Default")) {
      _testIfPressure[index] = 0;
      _addAMRPressureStability[index] == false;
    }
    else {
      cout << "Error invalid property in function MultiLevelProblem::AssociatePropertyToSolution" << endl;
      exit (0);
    }
  }

// *******************************************************


  void MultiLevelSolution::PairSolution (const char solution_name[], const char solution_pair[]) {
    unsigned index = GetIndex (solution_name);
    unsigned indexPair = GetIndex (solution_pair);
    _solPairIndex[index] = indexPair;
    _solPairInverseIndex[indexPair] = index;
  }

// *******************************************************
  void MultiLevelSolution::Initialize (const char name[], InitFunc func) {
    Initialize (name, func, NULL, NULL);
  }

  void MultiLevelSolution::Initialize (const char name[], InitFuncMLProb func, const MultiLevelProblem* ml_prob) {
    Initialize (name, NULL, func, ml_prob);
  }

  void MultiLevelSolution::Initialize (const char name[], InitFunc func, InitFuncMLProb funcMLProb, const MultiLevelProblem* ml_prob) {

    unsigned i_start;
    unsigned i_end;

    if (!strcmp (name, "All") || !strcmp (name, "all") || !strcmp (name, "ALL")) {
      i_start = 0;
      i_end = _solType.size();
    }
    else {
      i_start = GetIndex (name);
      i_end = i_start + 1u;
    }

    for (unsigned i = i_start; i < i_end; i++) {
      unsigned sol_type = _solType[i];

      for (unsigned ig = 0; ig < _gridn; ig++) {
        unsigned num_el = _mlMesh->GetLevel (ig)->GetNumberOfElements();
        _solution[ig]->ResizeSolutionVector (_solName[i]);
        _solution[ig]->_Sol[i]->zero();

        if (func || funcMLProb) {
          double value;

          if (sol_type < 3) {
            for (int isdom = _iproc; isdom < _iproc + 1; isdom++) {
              for (int iel = _mlMesh->GetLevel (ig)->_elementOffset[isdom];
                   iel < _mlMesh->GetLevel (ig)->_elementOffset[isdom + 1]; iel++) {
                unsigned nloc_dof = _mlMesh->GetLevel (ig)->GetElementDofNumber (iel, sol_type);

                for (int j = 0; j < nloc_dof; j++) {
                  unsigned inode_Metis = _mlMesh->GetLevel (ig)->GetSolutionDof (j, iel, sol_type);
                  unsigned icoord_Metis = _mlMesh->GetLevel (ig)->GetSolutionDof (j, iel, 2);
                  std::vector < double > xx (3);
                  xx[0] = (*_mlMesh->GetLevel (ig)->_topology->_Sol[0]) (icoord_Metis);
                  xx[1] = (*_mlMesh->GetLevel (ig)->_topology->_Sol[1]) (icoord_Metis);
                  xx[2] = (*_mlMesh->GetLevel (ig)->_topology->_Sol[2]) (icoord_Metis);

                  value = (func) ? func (xx) : funcMLProb (ml_prob, xx, name);

                  _solution[ig]->_Sol[i]->set (inode_Metis, value);

                  if (_solTimeOrder[i] == 2) {
                    _solution[ig]->_SolOld[i]->set (inode_Metis, value);
                  }
                }
              }
            }
          }
          else if (sol_type < 5) {
            for (int isdom = _iproc; isdom < _iproc + 1; isdom++) {
              for (int iel = _mlMesh->GetLevel (ig)->_elementOffset[isdom];
                   iel < _mlMesh->GetLevel (ig)->_elementOffset[isdom + 1]; iel++) {
                unsigned nloc_dof = _mlMesh->GetLevel (ig)->GetElementDofNumber (iel, 0);
                std::vector < double > xx (3, 0.);

                for (int j = 0; j < nloc_dof; j++) {
                  unsigned icoord_Metis = _mlMesh->GetLevel (ig)->GetSolutionDof (j, iel, 2);
                  xx[0] += (*_mlMesh->GetLevel (ig)->_topology->_Sol[0]) (icoord_Metis);
                  xx[1] += (*_mlMesh->GetLevel (ig)->_topology->_Sol[1]) (icoord_Metis);
                  xx[2] += (*_mlMesh->GetLevel (ig)->_topology->_Sol[2]) (icoord_Metis);
                }

                xx[0] /= nloc_dof;
                xx[1] /= nloc_dof;
                xx[2] /= nloc_dof;

                value = (func) ? func (xx) : funcMLProb (ml_prob, xx, name);

                unsigned solDof = _mlMesh->GetLevel (ig)->GetSolutionDof (2, iel, sol_type);

                _solution[ig]->_Sol[i]->set (solDof, value);

                if (_solTimeOrder[i] == 2) {
                  _solution[ig]->_SolOld[i]->set (solDof, value);
                }
              }
            }
          }

          _solution[ig]->_Sol[i]->close();

          if (_solTimeOrder[i] == 2) {
            _solution[ig]->_SolOld[i]->close();
          }
        }
      }
    }

    return;
  }


//---------------------------------------------------------------------------------------------------
  unsigned MultiLevelSolution::GetIndex (const char name[]) const {
    unsigned index = 0;

    while (strcmp (_solName[index], name)) {
      index++;

      if (index == _solType.size()) {
        cout << "error! invalid solution name: " << name << " in entry GetIndex(...)" << endl;
        abort();
      }
    }

    return index;
  }

// *******************************************************
  unsigned MultiLevelSolution::GetSolutionType (const char name[]) {
    unsigned index = 0;

    while (strcmp (_solName[index], name)) {
      index++;

      if (index == _solType.size()) {
        cout << "error! invalid name entry GetSolType(...)" << endl;
        abort();
      }
    }

    return _solType[index];
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
  void MultiLevelSolution::AttachSetBoundaryConditionFunction (BoundaryFuncMLProb SetBoundaryConditionFunction_in) {
    _bdcFuncSetMLProb = true;
    _bdcFuncSet = false;
    _SetBoundaryConditionFunctionMLProb = SetBoundaryConditionFunction_in;
    return;
  }

//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::AttachSetBoundaryConditionFunction (BoundaryFunc SetBoundaryConditionFunction) {
    _bdcFuncSet = true;
    _bdcFuncSetMLProb = false;
    _SetBoundaryConditionFunction = SetBoundaryConditionFunction;
    return;
  }

//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::InitializeBdc() {
    _useParsedBCFunction = true;

    int nvars = _solType.size();
    int nfaces = _mlMesh->GetLevel (0)->_boundaryinfo.size();
    _boundaryConditions.resize (nvars);
    _isHomogeneous.resize (nvars);
    _nonHomogeneousBCFunction.resize (nvars);

    for (int i = 0; i < nvars; i++) {
      _boundaryConditions[i].resize (nfaces);
      _isHomogeneous[i].resize (nfaces);
      _nonHomogeneousBCFunction[i].resize (nfaces);

      for (int j = 0; j < nfaces; j++) {
        _boundaryConditions[i][j] = DIRICHLET;
        _isHomogeneous[i][j] = true;
        _nonHomogeneousBCFunction[i][j] = NULL;
      }
    }
  }

//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::SetBoundaryCondition_new (const std::string name, const std::string facename,
                                                     const BDCType bdctype, const bool istimedependent, FunctionBase* func) {



    unsigned int ivar = GetIndex (name.c_str());
    unsigned int iface = 0;
    bool ishomogeneous = true;

    if (func != NULL) {
      ishomogeneous = false;
    }

    std::map<unsigned int, std::string>::iterator iter;
    iter = _mlMesh->GetLevel (0)->_boundaryinfo.begin();

    for (iter = _mlMesh->GetLevel (0)->_boundaryinfo.begin(); iter != _mlMesh->GetLevel (0)->_boundaryinfo.end(); ++iter) {
      if (iter->second.compare (facename) == 0) {
        iface = iter->first;
        break;
      }
    }

    if (iter == _mlMesh->GetLevel (0)->_boundaryinfo.end()) {
      std::cout << " Error: the facename " << facename << " does not exist!" << std::endl;
      exit (1);
    }

    _boundaryConditions[ivar][iface]       = bdctype;
    _isHomogeneous[ivar][iface]            = ishomogeneous;
    _nonHomogeneousBCFunction[ivar][iface] = func;

  }







//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::GenerateBdc (const char* name, const char* bdc_type, const MultiLevelProblem* ml_prob) {

    if (_useParsedBCFunction == false && _bdcFuncSet == false && _bdcFuncSetMLProb == false) {
      cout << "Error: The boundary condition user-function is not set! Please call the AttachSetBoundaryConditionFunction routine"
           << endl;

      abort();

    }

    if (_bdcFuncSetMLProb == true && ml_prob != NULL) {
      _mlBCProblem = ml_prob;
    }
    else if (_bdcFuncSetMLProb == true && ml_prob == NULL) {
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

      for (unsigned k = i_start; k < i_end; k++) {
        if (_solution[0]->_ResEpsBdcFlag[k]) {
          sprintf (_bdcType[k], "Steady");
          cout << " Set " << std::setw (15) << _bdcType[k] << " Boundary_condition"
               << " for variable " << std::setw (3) << _solName[k] << endl;
        }
        else {
          sprintf (_bdcType[k], "Not-available");
        }

      }
    }
    else {
      i_start = GetIndex (name);
      i_end = i_start + 1u;

      if (_solution[0]->_ResEpsBdcFlag[i_start]) {
        if (!strcmp (bdc_type, "Steady")) {
          strcpy (_bdcType[i_start], bdc_type);
        }
        else if (!strcmp (bdc_type, "Time_dependent")) {
          strcpy (_bdcType[i_start], bdc_type);
        }
        else {
          cout << "Error! Invalid boundary condition specified for " << _solName[i_start]
               << " in GenerateBdc function" << endl;
          exit (1);
        }

        cout << " Set " << std::setw (14) << _bdcType[i_start] << " Boundary_condition"
             << " for variable " << std::setw (3) << _solName[i_start] << endl;
      }
      else {
        sprintf (_bdcType[i_start], "Not-available");
      }
    }

    for (unsigned i = i_start; i < i_end; i++) {
      GenerateBdc (i, 0, 0.);
    }
  }

//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::UpdateBdc (const double time) {

    for (int k = 0; k < _solName.size(); k++) {
      if (!strcmp (_bdcType[k], "Time_dependent")) {
        GenerateBdc (k, 0, time);
      }
    }
  }


//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::GenerateBdc (const unsigned int k, const unsigned int grid0, const double time) {

    // 2 Default Neumann
    // 1 AMR artificial Dirichlet = 0 BC
    // 0 Dirichlet
    for (unsigned igridn = grid0; igridn < _gridn; igridn++) {
      if (_solution[igridn]->_ResEpsBdcFlag[k]) {
        Mesh* msh = _mlMesh->GetLevel (igridn);

        std::vector < std::map < unsigned,  std::map < unsigned, double  > > > &amrRestriction = msh->GetAmrRestrictionMap();

        // default Neumann
        for (unsigned j = msh->_dofOffset[_solType[k]][_iproc]; j < msh->_dofOffset[_solType[k]][_iproc + 1]; j++) {
          _solution[igridn]->_Bdc[k]->set (j, 2.);
        }

        if(_solType[k] < 3) {  // boundary condition for lagrangian elements

       //AMR - related            
          for(int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {
            for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
              if(msh->el->GetBoundaryIndex(iel, jface) == 0) {   // interior boundary (AMR) u = 0
                short unsigned ielt = msh->GetElementType(iel);
                unsigned nv1 = msh->GetElementFaceDofNumber(iel, jface, _solType[k]);  // only the face dofs
                for(unsigned iv = 0; iv < nv1; iv++) {
                  unsigned i = msh->GetLocalFaceVertexIndex(iel, jface, iv);
                  unsigned idof = msh->GetSolutionDof(i, iel, _solType[k]);
                  if(amrRestriction[_solType[k]].find(idof) != amrRestriction[_solType[k]].end() &&
                      amrRestriction[_solType[k]][idof][idof] == 0) {
                    _solution[igridn]->_Bdc[k]->set (idof, 1.);
                  }
                }
              }
            }
          }
       //AMR - related - end
       

          for(int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {
              
            for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
                
                const int boundary_index = msh->el->GetBoundaryIndex(iel, jface);
                
              if(boundary_index > 0) {   // exterior boundary u = value
                short unsigned ielt = msh->GetElementType(iel);
                unsigned n_face_dofs = msh->GetElementFaceDofNumber(iel, jface, _solType[k]);

                for(unsigned iv = 0; iv < n_face_dofs; iv++) {
                  unsigned i = msh->GetLocalFaceVertexIndex(iel, jface, iv);
                  unsigned inode_coord_Metis = msh->GetSolutionDof(i, iel, 2);

                  if(_useParsedBCFunction) {
                    unsigned int faceIndex = boundary_index;

                    if (GetBoundaryCondition (k, faceIndex - 1u) == DIRICHLET) {
                      unsigned inode_Metis = msh->GetSolutionDof (i, iel, _solType[k]);
                      _solution[igridn]->_Bdc[k]->set (inode_Metis, 0.);
                      double value = 0.;

                      if (!Ishomogeneous (k, faceIndex - 1u)) {
                        ParsedFunction* bdcfunc = (ParsedFunction*) (GetBdcFunction (k, faceIndex - 1u));
                        double xyzt[4];
                        xyzt[0] = (*msh->_topology->_Sol[0]) (inode_coord_Metis);
                        xyzt[1] = (*msh->_topology->_Sol[1]) (inode_coord_Metis);
                        xyzt[2] = (*msh->_topology->_Sol[2]) (inode_coord_Metis);
                        xyzt[3] = time;
                        value = (*bdcfunc) (xyzt);
                      }

                      _solution[igridn]->_Sol[k]->set (inode_Metis, value);
                    }
                  }
                  else {
                    double value;
                    std::vector < double > xx (3);
                    xx[0] = (*msh->_topology->_Sol[0]) (inode_coord_Metis);
                    xx[1] = (*msh->_topology->_Sol[1]) (inode_coord_Metis);
                    xx[2] = (*msh->_topology->_Sol[2]) (inode_coord_Metis);
                    bool test = (_bdcFuncSetMLProb) ?

                                _SetBoundaryConditionFunctionMLProb(_mlBCProblem, xx, _solName[k], value, boundary_index, time) :
                                _SetBoundaryConditionFunction(xx, _solName[k], value, boundary_index, time);


                    if (test) {
                      unsigned idof = msh->GetSolutionDof (i, iel, _solType[k]);
                      _solution[igridn]->_Bdc[k]->set (idof, 0.);
                      _solution[igridn]->_Sol[k]->set (idof, value);
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


  void MultiLevelSolution::GenerateBdcOnVolumeConstraint (const std::vector<unsigned> &volumeConstraintFlags, const unsigned &solIndex, const unsigned &grid0) {

    unsigned solType = GetSolutionType (solIndex);

    for (unsigned igridn = grid0; igridn < _gridn; igridn++) {
      Mesh* msh = _mlMesh->GetLevel (igridn);
      for (int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {

        short unsigned ielGroup = msh->GetElementGroup (iel);

        bool ielIsInVolumeConstraint = false;

        for (unsigned i = 0; i < volumeConstraintFlags.size(); i++) {
          if (volumeConstraintFlags[i] == ielGroup) {
            ielIsInVolumeConstraint = true;
            break;
          }
        }

        if (ielIsInVolumeConstraint) {

          unsigned nDofu  = msh->GetElementDofNumber (iel, solType);

          for (unsigned i = 0; i < nDofu; i++) {
            unsigned solDof = msh->GetSolutionDof (i, iel, solType);
            _solution[igridn]->_Bdc[solIndex]->set (solDof, 0.);
          

          double value;
          std::vector < double > xCord (3);
          unsigned inode_coord_Metis = msh->GetSolutionDof (i, iel, 2);
          xCord[0] = (*msh->_topology->_Sol[0]) (inode_coord_Metis);
          xCord[1] = (*msh->_topology->_Sol[1]) (inode_coord_Metis);
          xCord[2] = (*msh->_topology->_Sol[2]) (inode_coord_Metis);
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
  void MultiLevelSolution::GenerateRKBdc (
    const unsigned int &solIndex, const std::vector<unsigned> &solKiIndex,
    const unsigned int &grid0, const std::vector < double> & itime, const double &time0, const double &dt, const double *AI) {

    for (unsigned k = 0; k < solKiIndex.size(); k++) {
      strcpy (_bdcType[solKiIndex[k]], "RKbdc");
    }
    // 2 Default Neumann
    // 1 AMR artificial Dirichlet = 0 BC
    // 0 Dirichlet
    for (unsigned igridn = grid0; igridn < _gridn; igridn++) {
      if (_solution[igridn]->_ResEpsBdcFlag[solIndex]) {
        Mesh* msh = _mlMesh->GetLevel (igridn);

        std::vector < std::map < unsigned,  std::map < unsigned, double  > > > &amrRestriction = msh->GetAmrRestrictionMap();

        // default Neumann
        for (unsigned j = msh->_dofOffset[_solType[solIndex]][_iproc]; j < msh->_dofOffset[_solType[solIndex]][_iproc + 1]; j++) {
          for (unsigned k = 0; k < solKiIndex.size(); k++) {
            _solution[igridn]->_Bdc[solKiIndex[k]]->set (j, 2.);
          }
        }

        if (_solType[solIndex] < 3) { // boundary condition for lagrangian elements
          for (int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {
            for (unsigned jface = 0; jface < msh->GetElementFaceNumber (iel); jface++) {
              if (msh->el->GetBoundaryIndex (iel, jface) == 0) { // interior boundary (AMR) u = 0
                short unsigned ielt = msh->GetElementType (iel);
                unsigned nv1 = msh->GetElementFaceDofNumber (iel, jface, _solType[solIndex]); // only the face dofs
                for (unsigned iv = 0; iv < nv1; iv++) {
                  unsigned i = msh->GetLocalFaceVertexIndex (iel, jface, iv);
                  unsigned idof = msh->GetSolutionDof (i, iel, _solType[solIndex]);
                  if (amrRestriction[_solType[solIndex]].find (idof) != amrRestriction[_solType[solIndex]].end() &&
                      amrRestriction[_solType[solIndex]][idof][idof] == 0) {
                    for (unsigned k = 0; k < solKiIndex.size(); k++) {
                      _solution[igridn]->_Bdc[solKiIndex[k]]->set (idof, 1.);
                    }
                  }
                }
              }
            }
          }

          for (int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {
            for (unsigned jface = 0; jface < msh->GetElementFaceNumber (iel); jface++) {
              if (msh->el->GetBoundaryIndex (iel, jface) > 0) { // exterior boundary u = value
                short unsigned ielt = msh->GetElementType (iel);
                unsigned nv1 = msh->GetElementFaceDofNumber (iel, jface, _solType[solIndex]);

                for (unsigned iv = 0; iv < nv1; iv++) {
                  unsigned i = msh->GetLocalFaceVertexIndex (iel, jface, iv);
                  unsigned inode_coord_Metis = msh->GetSolutionDof (i, iel, 2);

                  if (_useParsedBCFunction) {
                    unsigned int faceIndex = msh->el->GetBoundaryIndex (iel, jface);

                    if (GetBoundaryCondition (solIndex, faceIndex - 1u) == DIRICHLET) {
                      unsigned inode_Metis = msh->GetSolutionDof (i, iel, _solType[solIndex]);
                      for (unsigned k = 0; k < solKiIndex.size(); k++) {
                        _solution[igridn]->_Bdc[solKiIndex[k]]->set (inode_Metis, 0.);
                      }
                      std::vector < double > ivalue (solKiIndex.size(), 0.);

                      if (!Ishomogeneous (solIndex, faceIndex - 1u)) {
                        ParsedFunction* bdcfunc = (ParsedFunction*) (GetBdcFunction (solIndex, faceIndex - 1u));
                        double xyzt[4];
                        xyzt[0] = (*msh->_topology->_Sol[0]) (inode_coord_Metis);
                        xyzt[1] = (*msh->_topology->_Sol[1]) (inode_coord_Metis);
                        xyzt[2] = (*msh->_topology->_Sol[2]) (inode_coord_Metis);
                        xyzt[3] = time0;
                        double value0 = (*bdcfunc) (xyzt);

                        for (unsigned k = 0; k < solKiIndex.size(); k++) {
                          xyzt[3] = itime[k];
                          ivalue[k] = (*bdcfunc) (xyzt) - value0;
                        }
                        if (_solTimeOrder[solIndex] == 2) {
                          _solution[igridn]->_SolOld[solIndex]->set (inode_Metis, value0);
                        }
                        for (unsigned k1 = 0; k1 < solKiIndex.size(); k1++) {
                          double value = 0.;
                          for (unsigned k2 = 0; k2 < solKiIndex.size(); k2++) {
                            value += AI[k1 * solKiIndex.size() + k2] * ivalue[k2];
                          }
                          _solution[igridn]->_Sol[solKiIndex[k1]]->set (inode_Metis, value / dt);
                        }
                      }
                    }
                  }
                  else {
                    double value0;
                    std::vector < double > xx (3);
                    xx[0] = (*msh->_topology->_Sol[0]) (inode_coord_Metis);
                    xx[1] = (*msh->_topology->_Sol[1]) (inode_coord_Metis);
                    xx[2] = (*msh->_topology->_Sol[2]) (inode_coord_Metis);
                    bool test = (_bdcFuncSetMLProb) ?
                                _SetBoundaryConditionFunctionMLProb (_mlBCProblem, xx, _solName[solIndex], value0, msh->el->GetBoundaryIndex (iel, jface), time0) :
                                _SetBoundaryConditionFunction (xx, _solName[solIndex], value0, msh->el->GetBoundaryIndex (iel, jface), time0);

                    if (test) {

                      if (PRINT) std::cout << inode_coord_Metis << " " << value0 << "\n" ;

                      std::vector < double > ivalue (solKiIndex.size(), 0.);
                      for (unsigned k = 0; k < solKiIndex.size(); k++) {
                        if (_bdcFuncSetMLProb) {
                          _SetBoundaryConditionFunctionMLProb (_mlBCProblem, xx, _solName[solIndex], ivalue[k], msh->el->GetBoundaryIndex (iel, jface), itime[k]);
                        }
                        else {
                          _SetBoundaryConditionFunction (xx, _solName[solIndex], ivalue[k], msh->el->GetBoundaryIndex (iel, jface), itime[k]);
                        }
                        if (PRINT) std::cout << inode_coord_Metis << " " << ivalue[k] << std::endl;
                        ivalue[k] -= value0;
                      }

                      unsigned idof = msh->GetSolutionDof (i, iel, _solType[solIndex]);

                      if (_solTimeOrder[solIndex] == 2) {
                        _solution[igridn]->_SolOld[solIndex]->set (idof, value0);
                      }
                      for (unsigned k1 = 0; k1 < solKiIndex.size(); k1++) {
                        double value = 0.;
                        for (unsigned k2 = 0; k2 < solKiIndex.size(); k2++) {
                          value += AI[k1 * solKiIndex.size() + k2] * ivalue[k2];
                        }
                        _solution[igridn]->_Bdc[solKiIndex[k1]]->set (idof, 0.);
                        _solution[igridn]->_Sol[solKiIndex[k1]]->set (idof, value / dt);
                      }
                    }
                  }
                }
              }
            }
          }
        }
        if (_fixSolutionAtOnePoint[solIndex] == true  && igridn == 0 && _iproc == 0) {
          for (unsigned k = 0; k < solKiIndex.size(); k++) {

            _solution[igridn]->_Bdc[solKiIndex[k]]->set (0, 0.);
            _solution[igridn]->_Sol[solKiIndex[k]]->set (0, 0.);

          }
        }
        if (_solTimeOrder[solIndex] == 2) {
          _solution[igridn]->_SolOld[solIndex]->close();
        }
        for (unsigned k = 0; k < solKiIndex.size(); k++) {
          _solution[igridn]->_Sol[solKiIndex[k]]->close();
          _solution[igridn]->_Bdc[solKiIndex[k]]->close();

        }
      }
    }

    PRINT = false;
  }



  void MultiLevelSolution::SaveSolution (const char* filename, const double &time) {

    char composedFileName[100];

    for (int i = 0; i < _solName.size(); i++) {
      sprintf (composedFileName, "./save/%s_time%f_sol%s_level%d", filename, time, _solName[i], _gridn);
      _solution[_gridn - 1]->_Sol[i]->BinaryPrint (composedFileName);
    }
  }

  void MultiLevelSolution::SaveSolution (const char* filename, const unsigned &iteration) {

    char composedFileName[100];

    for (int i = 0; i < _solName.size(); i++) {
      sprintf (composedFileName, "./save/%s_iteration%d_sol%s_level%d", filename, iteration, _solName[i], _gridn);
      _solution[_gridn - 1]->_Sol[i]->BinaryPrint (composedFileName);
    }
  }

  void MultiLevelSolution::LoadSolution (const char* filename) {
    LoadSolution (_gridn, filename);
  }

  void MultiLevelSolution::LoadSolution (const unsigned &level, const char* filename) {

    if (level > _gridn) {
      std::cout << "Error in MultiLevelSolution::LoadSolution function:" << std::endl;
      std::cout << "the solution level = " << level << " is not available in this MultilevelSolution" << std::endl;
      abort();
    }

    char composedFileName[200];
    for (int i = 0; i < _solName.size(); i++) {
      sprintf (composedFileName, "%s_sol%s_level%d", filename, _solName[i], level);
      // check if the file really exists
      if (strncmp (filename, "http://", 7) && strncmp (filename, "ftp://", 6)) {
        struct stat buffer;
        if (stat (composedFileName, &buffer) != 0) {
          std::cerr << "Error: cannot locate file " << composedFileName << std::endl;
          abort();
        }
      }
      _solution[level - 1]->_Sol[i]->BinaryLoad (composedFileName);
    }

    for (int gridf = level; gridf < _gridn; gridf++) {
      for (unsigned i = 0; i < _solName.size(); i++) {
        _solution[gridf]->_Sol[i]->matrix_mult (*_solution[gridf - 1]->_Sol[i],
                                                *_mlMesh->GetLevel (gridf)->GetCoarseToFineProjection (_solType[i]));
        _solution[gridf]->_Sol[i]->close();
      }
    }


  }





  void MultiLevelSolution::RefineSolution (const unsigned &gridf) {

    Mesh *msh = _mlMesh->GetLevel (gridf);

    for (unsigned k = 0; k < _solType.size(); k++) {

      unsigned solType = _solType[k];
      _solution[gridf]->_Sol[k]->matrix_mult (*_solution[gridf - 1]->_Sol[k],
                                              *msh->GetCoarseToFineProjection (solType));
      _solution[gridf]->_Sol[k]->close();
    }
  }



    void MultiLevelSolution::CoarsenSolutionByOneLevel_wrong(const unsigned &grid_fine)  {

    Mesh *msh = _mlMesh->GetLevel(grid_fine);
    const unsigned grid_coarse = grid_fine - 1;
    
    for(unsigned k = 0; k < _solType.size(); k++) {

      unsigned solType = _solType[k];
      _solution[grid_coarse]->_Sol[k]->matrix_mult_transpose( *(_solution[grid_fine]->_Sol[k]), *(msh->GetCoarseToFineProjectionRestrictionOnCoarse(solType)) );
      _solution[grid_coarse]->_Sol[k]->close();
    }
    
  }

    void MultiLevelSolution::CoarsenSolutionByOneLevel(const unsigned &grid_fine)  {
        
     const unsigned grid_coarse = grid_fine - 1;
     Mesh *msh = _mlMesh->GetLevel(grid_coarse);
     
        //loop over the coarse elements
        //loop over the dofs of each coarse element
        //loop over the children elements of that coarse element
        //find the child element to which the coarse dof belongs, and find the child dof
        //set the _Sol in the coarse dofs to the fine dofs
        
    }
  
  /** Copies from another MLSol object from a given level to some other level.
      One should also check that they belong to the same underlying mesh structure */
  void MultiLevelSolution::fill_at_level_from_level(const unsigned lev_out, const unsigned lev_in, const MultiLevelSolution & ml_sol_in)  {
      
      assert(_solType.size() == ml_sol_in.GetSolutionSize());
      
          for(unsigned k = 0; k < _solType.size(); k++) {
              *(_solution[lev_out]->_Sol[k]) = *(ml_sol_in.GetSolutionLevel(lev_in)->_Sol[k]);
          }
          
  }

  
  
  void MultiLevelSolution::CopySolutionToOldSolution()
  {
    for(unsigned short i = 0; i < _gridn; i++) {
      _solution[i]->CopySolutionToOldSolution();
    }
  }

  void MultiLevelSolution::UpdateSolution (const char name[], InitFunc func, const double& time) {
    unsigned i = GetIndex (name);

    unsigned sol_type = _solType[i];

    for (unsigned ig = 0; ig < _gridn; ig++) {
      if (sol_type < 3) {
        for (int isdom = _iproc; isdom < _iproc + 1; isdom++) {
          for (int iel = _mlMesh->GetLevel (ig)->_elementOffset[isdom];
               iel < _mlMesh->GetLevel (ig)->_elementOffset[isdom + 1]; iel++) {
            unsigned nloc_dof = _mlMesh->GetLevel (ig)->GetElementDofNumber (iel, sol_type);
            for (int j = 0; j < nloc_dof; j++) {
              unsigned inode_Metis = _mlMesh->GetLevel (ig)->GetSolutionDof (j, iel, sol_type);
              unsigned icoord_Metis = _mlMesh->GetLevel (ig)->GetSolutionDof (j, iel, 2);
              std::vector < double > xx (4);
              xx[0] = (*_mlMesh->GetLevel (ig)->_topology->_Sol[0]) (icoord_Metis);
              xx[1] = (*_mlMesh->GetLevel (ig)->_topology->_Sol[1]) (icoord_Metis);
              xx[2] = (*_mlMesh->GetLevel (ig)->_topology->_Sol[2]) (icoord_Metis);
              xx[3] = time;
              double value = func (xx);
              _solution[ig]->_Sol[i]->set (inode_Metis, value);
            }
          }
        }
      }
      else if (sol_type < 5) {
        for (int isdom = _iproc; isdom < _iproc + 1; isdom++) {
          for (int iel = _mlMesh->GetLevel (ig)->_elementOffset[isdom];
               iel < _mlMesh->GetLevel (ig)->_elementOffset[isdom + 1]; iel++) {
            unsigned nloc_dof = _mlMesh->GetLevel (ig)->GetElementDofNumber (iel, 0);
            std::vector < double > xx (4, 0.);

            for (int j = 0; j < nloc_dof; j++) {
              unsigned icoord_Metis = _mlMesh->GetLevel (ig)->GetSolutionDof (j, iel, 2);
              xx[0] += (*_mlMesh->GetLevel (ig)->_topology->_Sol[0]) (icoord_Metis);
              xx[1] += (*_mlMesh->GetLevel (ig)->_topology->_Sol[1]) (icoord_Metis);
              xx[2] += (*_mlMesh->GetLevel (ig)->_topology->_Sol[2]) (icoord_Metis);
            }
            xx[0] /= nloc_dof;
            xx[1] /= nloc_dof;
            xx[2] /= nloc_dof;
            xx[3] = time;

            double value =  func (xx);
            unsigned solDof = _mlMesh->GetLevel (ig)->GetSolutionDof (2, iel, sol_type);
            _solution[ig]->_Sol[i]->set (solDof, value);
          }
        }
      }
      _solution[ig]->_Sol[i]->close();
    }
    return;
  }



} //end namespace femus



