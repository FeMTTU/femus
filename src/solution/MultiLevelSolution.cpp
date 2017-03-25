/*=========================================================================

Program: FEMUS
Module: MultiLevelProblem
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
  MultiLevelSolution::~MultiLevelSolution() {

    for(unsigned i = 0; i < _gridn; i++) {
      _solution[i]->FreeSolutionVectors();
      delete _solution[i];
    }

    for(unsigned i = 0; i < _solName.size(); i++) delete [] _solName[i];

    for(unsigned i = 0; i < _solName.size(); i++) delete [] _bdcType[i];


  };

//---------------------------------------------------------------------------------------------------
  MultiLevelSolution::MultiLevelSolution(MultiLevelMesh* ml_msh) :
    _gridn(ml_msh->GetNumberOfLevels()),
    _mlMesh(ml_msh) {
    _solution.resize(_gridn);

    for(unsigned i = 0; i < _gridn; i++) {
      _solution[i] = new Solution(_mlMesh->GetLevel(i));
    }

    _bdcFuncSet = false;
    _bdcFuncSetMLProb = false;
    _useParsedBCFunction = false;

    _mlBCProblem = NULL;

  }

  void MultiLevelSolution::AddSolutionLevel() {
    // add level solution
    _solution.resize(_gridn + 1);
    _solution[_gridn] = new Solution(_mlMesh->GetLevel(_gridn));

    // add all current solutions and initialize to zero
    for(unsigned i = 0; i < _solName.size(); i++) {
      _solution[_gridn]->AddSolution(_solName[i], _family[i], _order[i], _solTimeOrder[i], _pdeType[i]);
    }

    for(unsigned i = 0; i < _solName.size(); i++) {
      _solution[_gridn]->ResizeSolutionVector(_solName[i]);

      _solution[_gridn]->_Sol[i]->zero();

      if(_solTimeOrder[i] == 2) {
        _solution[_gridn]->_SolOld[i]->zero();
      }
    }

    _gridn++;
    unsigned  grid0 = _gridn - 1;

    for(int k = 0; k < _solName.size(); k++) {
      GenerateBdc(k, grid0, 0.);
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


    _testIfPressure[n] = 0;
    _addAMRPressureStability[n] = false;
    _fixSolutionAtOnePoint[n] = false;
    _family[n] = fefamily;
    _order[n] = order;
    _solType[n] = order - ((fefamily == LAGRANGE) ? 1 : 0) + fefamily * 3;
    _solName[n]  = new char [DEFAULT_SOL_NCHARS];
    _bdcType[n]  = new char [20];
    strcpy(_solName[n], name);
    _solTimeOrder[n] = tmorder;
    _pdeType[n] = PdeType;
    _solPairIndex[n] = n;

    cout << " Add variable " << std::setw(3) << _solName[n] << " discretized with FE type "
         << std::setw(12) << order << " and time discretzation order " << tmorder << endl;

    for(unsigned ig = 0; ig < _gridn; ig++) {
      _solution[ig]->AddSolution(_solName[n], _family[n], _order[n], _solTimeOrder[n], _pdeType[n]);
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
      cout << "Error invalid property in function MultiLevelProblem::AssociatePropertyToSolution" << endl;
      exit(0);
    }
  }

// *******************************************************


  void MultiLevelSolution::PairSolution(const char solution_name[], const char solution_pair[]) {
    unsigned index = GetIndex(solution_name);
    unsigned indexPair = GetIndex(solution_pair);
    _solPairIndex[index] = indexPair;
  }

// *******************************************************
  void MultiLevelSolution::Initialize(const char name[], InitFunc func) {
    Initialize(name, func, NULL, NULL);
  }

  void MultiLevelSolution::Initialize(const char name[], InitFuncMLProb func, const MultiLevelProblem* ml_prob) {
    Initialize(name, NULL, func, ml_prob);
  }

  void MultiLevelSolution::Initialize(const char name[], InitFunc func, InitFuncMLProb funcMLProb, const MultiLevelProblem* ml_prob) {

    unsigned i_start;
    unsigned i_end;

    if(!strcmp(name, "All") || !strcmp(name, "all") || !strcmp(name, "ALL")) {
      i_start = 0;
      i_end = _solType.size();
    }
    else {
      i_start = GetIndex(name);
      i_end = i_start + 1u;
    }

    for(unsigned i = i_start; i < i_end; i++) {
      unsigned sol_type = _solType[i];

      for(unsigned ig = 0; ig < _gridn; ig++) {
        unsigned num_el = _mlMesh->GetLevel(ig)->GetNumberOfElements();
        _solution[ig]->ResizeSolutionVector(_solName[i]);
        _solution[ig]->_Sol[i]->zero();

        if(func || funcMLProb) {
          double value;

          if(sol_type < 3) {
            for(int isdom = _iproc; isdom < _iproc + 1; isdom++) {
              for(int iel = _mlMesh->GetLevel(ig)->_elementOffset[isdom];
                  iel < _mlMesh->GetLevel(ig)->_elementOffset[isdom + 1]; iel++) {
                unsigned nloc_dof = _mlMesh->GetLevel(ig)->GetElementDofNumber(iel, sol_type);

                for(int j = 0; j < nloc_dof; j++) {
                  unsigned inode_Metis = _mlMesh->GetLevel(ig)->GetSolutionDof(j, iel, sol_type);
                  unsigned icoord_Metis = _mlMesh->GetLevel(ig)->GetSolutionDof(j, iel, 2);
                  std::vector < double > xx(3);
                  xx[0] = (*_mlMesh->GetLevel(ig)->_topology->_Sol[0])(icoord_Metis);
                  xx[1] = (*_mlMesh->GetLevel(ig)->_topology->_Sol[1])(icoord_Metis);
                  xx[2] = (*_mlMesh->GetLevel(ig)->_topology->_Sol[2])(icoord_Metis);

                  value = (func) ? func(xx) : funcMLProb(ml_prob, xx, name);

                  _solution[ig]->_Sol[i]->set(inode_Metis, value);

                  if(_solTimeOrder[i] == 2) {
                    _solution[ig]->_SolOld[i]->set(inode_Metis, value);
                  }
                }
              }
            }
          }
          else if(sol_type < 5) {
            for(int isdom = _iproc; isdom < _iproc + 1; isdom++) {
              for(int iel = _mlMesh->GetLevel(ig)->_elementOffset[isdom];
                  iel < _mlMesh->GetLevel(ig)->_elementOffset[isdom + 1]; iel++) {
                unsigned nloc_dof = _mlMesh->GetLevel(ig)->GetElementDofNumber(iel, 0);
                std::vector < double > xx(3, 0.);

                for(int j = 0; j < nloc_dof; j++) {
                  unsigned icoord_Metis = _mlMesh->GetLevel(ig)->GetSolutionDof(j, iel, 2);
                  xx[0] += (*_mlMesh->GetLevel(ig)->_topology->_Sol[0])(icoord_Metis);
                  xx[1] += (*_mlMesh->GetLevel(ig)->_topology->_Sol[1])(icoord_Metis);
                  xx[2] += (*_mlMesh->GetLevel(ig)->_topology->_Sol[2])(icoord_Metis);
                }

                xx[0] /= nloc_dof;
                xx[1] /= nloc_dof;
                xx[2] /= nloc_dof;

                value = (func) ? func(xx) : funcMLProb(ml_prob, xx, name);

                _solution[ig]->_Sol[i]->set(iel, value);

                if(_solTimeOrder[i] == 2) {
                  _solution[ig]->_SolOld[i]->set(iel, value);
                }
              }
            }
          }

          _solution[ig]->_Sol[i]->close();

          if(_solTimeOrder[i] == 2) {
            _solution[ig]->_SolOld[i]->close();
          }
        }
      }
    }

    return;
  }


//---------------------------------------------------------------------------------------------------
  unsigned MultiLevelSolution::GetIndex(const char name[]) const {
    unsigned index = 0;

    while(strcmp(_solName[index], name)) {
      index++;

      if(index == _solType.size()) {
        cout << "error! invalid solution name " << name << "entry GetIndex(...)" << endl;
        abort();
      }
    }

    return index;
  }

// *******************************************************
  unsigned MultiLevelSolution::GetSolutionType(const char name[]) {
    unsigned index = 0;

    while(strcmp(_solName[index], name)) {
      index++;

      if(index == _solType.size()) {
        cout << "error! invalid name entry GetSolType(...)" << endl;
        abort();
      }
    }

    return _solType[index];
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
    int nfaces = _mlMesh->GetLevel(0)->_boundaryinfo.size();
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
    iter = _mlMesh->GetLevel(0)->_boundaryinfo.begin();

    for(iter = _mlMesh->GetLevel(0)->_boundaryinfo.begin(); iter != _mlMesh->GetLevel(0)->_boundaryinfo.end(); ++iter) {
      if(iter->second.compare(facename) == 0) {
        iface = iter->first;
        break;
      }
    }

    if(iter == _mlMesh->GetLevel(0)->_boundaryinfo.end()) {
      std::cout << " Error: the facename " << facename << " does not exist!" << std::endl;
      exit(1);
    }

    _boundaryConditions[ivar][iface]       = bdctype;
    _isHomogeneous[ivar][iface]            = ishomogeneous;
    _nonHomogeneousBCFunction[ivar][iface] = func;

  }







//---------------------------------------------------------------------------------------------------
  void MultiLevelSolution::GenerateBdc(const char* name, const char* bdc_type, const MultiLevelProblem* ml_prob) {

    if(_useParsedBCFunction == false && _bdcFuncSet == false && _bdcFuncSetMLProb == false) {
      cout << "Error: The boundary condition user-function is not set! Please call the AttachSetBoundaryConditionFunction routine"
           << endl;

      exit(1);
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

    if(!strcmp(name, "All")) {
      i_start = 0;
      i_end = _solType.size();

      for(unsigned k = i_start; k < i_end; k++) {
        if(_solution[0]->_ResEpsBdcFlag[k]) {
          sprintf(_bdcType[k], "Steady");
          cout << " Set " << std::setw(15) << _bdcType[k] << " Boundary_condition"
               << " for variable " << std::setw(3) << _solName[k] << endl;
        }
        else {
          sprintf(_bdcType[k], "Not-available");
        }

      }
    }
    else {
      i_start = GetIndex(name);
      i_end = i_start + 1u;

      if(_solution[0]->_ResEpsBdcFlag[i_start]) {
        if(!strcmp(bdc_type, "Steady")) {
          strcpy(_bdcType[i_start], bdc_type);
        }
        else if(!strcmp(bdc_type, "Time_dependent")) {
          strcpy(_bdcType[i_start], bdc_type);
        }
        else {
          cout << "Error! Invalid boundary condition specified for " << _solName[i_start]
               << " in GenerateBdc function" << endl;
          exit(1);
        }

        cout << " Set " << std::setw(14) << _bdcType[i_start] << " Boundary_condition"
             << " for variable " << std::setw(3) << _solName[i_start] << endl;
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
  void MultiLevelSolution::GenerateBdc(const unsigned int k, const unsigned int grid0, const double time) {

    // 2 Default Neumann
    // 1 AMR artificial Dirichlet = 0 BC
    // 0 Dirichlet
    for(unsigned igridn = grid0; igridn < _gridn; igridn++) {
      if(_solution[igridn]->_ResEpsBdcFlag[k]) {
        Mesh* msh = _mlMesh->GetLevel(igridn);

        // default Neumann
        for(unsigned j = msh->_dofOffset[_solType[k]][_iproc]; j < msh->_dofOffset[_solType[k]][_iproc + 1]; j++) {
          _solution[igridn]->_Bdc[k]->set(j, 2.);
        }

        if(_solType[k] < 3) {  // boundary condition for lagrangian elements
          for(int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {
            for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
              if(msh->el->GetBoundaryIndex(iel, jface) == 0) {   // interior boundary (AMR) u = 0
                short unsigned ielt = msh->GetElementType(iel);
                unsigned nv1 = (_addAMRPressureStability[k] == true) ?
                               msh->GetElementDofNumber(iel, _solType[k]) :  //all the dofs in the element
                               msh->GetElementFaceDofNumber(iel, jface, _solType[k]);  // only the face dofs

                for(unsigned iv = 0; iv < nv1; iv++) {
                  unsigned i = (_addAMRPressureStability[k] == true) ? iv : msh->GetLocalFaceVertexIndex(iel, jface, iv);
                  unsigned idof = msh->GetSolutionDof(i, iel, _solType[k]);
                  _solution[igridn]->_Bdc[k]->set(idof, 1.);
                }
              }
            }
          }

          for(int iel = msh->_elementOffset[_iproc]; iel < msh->_elementOffset[_iproc + 1]; iel++) {
            for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
              if(msh->el->GetBoundaryIndex(iel, jface) > 0) {   // exterior boundary u = value
                short unsigned ielt = msh->GetElementType(iel);
                unsigned nv1 = msh->GetElementFaceDofNumber(iel, jface, _solType[k]);

                for(unsigned iv = 0; iv < nv1; iv++) {
                  unsigned i = msh->GetLocalFaceVertexIndex(iel, jface, iv);
                  unsigned inode_coord_Metis = msh->GetSolutionDof(i, iel, 2);

                  if(_useParsedBCFunction) {
                    unsigned int faceIndex = msh->el->GetBoundaryIndex(iel, jface);

                    if(GetBoundaryCondition(k, faceIndex - 1u) == DIRICHLET) {
                      unsigned inode_Metis = msh->GetSolutionDof(i, iel, _solType[k]);
                      _solution[igridn]->_Bdc[k]->set(inode_Metis, 0.);
                      double value = 0.;

                      if(!Ishomogeneous(k, faceIndex - 1u)) {
                        ParsedFunction* bdcfunc = (ParsedFunction*)(GetBdcFunction(k, faceIndex - 1u));
                        double xyzt[4];
                        xyzt[0] = (*msh->_topology->_Sol[0])(inode_coord_Metis);
                        xyzt[1] = (*msh->_topology->_Sol[1])(inode_coord_Metis);
                        xyzt[2] = (*msh->_topology->_Sol[2])(inode_coord_Metis);
                        xyzt[3] = time;
                        value = (*bdcfunc)(xyzt);
                      }

                      _solution[igridn]->_Sol[k]->set(inode_Metis, value);
                    }
                  }
                  else {
                    double value;
                    std::vector < double > xx(3);
                    xx[0] = (*msh->_topology->_Sol[0])(inode_coord_Metis);
                    xx[1] = (*msh->_topology->_Sol[1])(inode_coord_Metis);
                    xx[2] = (*msh->_topology->_Sol[2])(inode_coord_Metis);
                    bool test = (_bdcFuncSetMLProb) ?
                                _SetBoundaryConditionFunctionMLProb(_mlBCProblem, xx, _solName[k], value, msh->el->GetBoundaryIndex(iel, jface), time) :
                                _SetBoundaryConditionFunction(xx, _solName[k], value, msh->el->GetBoundaryIndex(iel, jface), time);

                    if(test) {
                      unsigned idof = msh->GetSolutionDof(i, iel, _solType[k]);
                      _solution[igridn]->_Bdc[k]->set(idof, 0.);
                      _solution[igridn]->_Sol[k]->set(idof, value);
                    }
                  }
                }
              }
            }
          }
        }
        else if(_addAMRPressureStability[k]) {  // interior boundary (AMR) for discontinuous elements u = 0
          unsigned offset = msh->_elementOffset[_iproc];
          unsigned offsetp1 = msh->_elementOffset[_iproc + 1];
          unsigned owned = offsetp1 - offset;
          std::vector < short unsigned > markedElement(owned, 0);
          // Add all interior boundary elements
          for(unsigned iel = offset; iel < offsetp1; iel++) {
            short unsigned ielt = msh->GetElementType(iel);
            for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
              if(msh->el->GetBoundaryIndex(iel, jface) == 0) {
                markedElement[iel - offset] = 1;
              }
            }
          }
          //remove adjacent interior boundary elements
          for(unsigned i = 0; i < owned; i++) {
            if( markedElement[i] == 1){
              markedElement[i] = 2;
              std::vector < unsigned > seed(1, offset + i);
              while(seed.size() != 0){
                bool testSeed = true;
                unsigned iel = seed[ seed.size() - 1u];
                short unsigned ielt = msh->GetElementType(iel);
                for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
                  int jel = msh->el->GetFaceElementIndex(seed[seed.size()-1], jface) - 1;
                  if( jel >= offset && jel < offsetp1 && markedElement[jel - offset] == 1 ){
                    markedElement[jel - offset] = 0;
                    seed.resize(seed.size() + 1);
                    seed[seed.size() - 1] = jel;
                    testSeed = false;
                  }
                }
                if(testSeed) seed.resize(seed.size() - 1);
              }
            }
          }

          for(unsigned i = 0; i < owned; i++) {
            if( markedElement[i] > 0){
              unsigned iel = offset + i;
              unsigned idof = msh->GetSolutionDof(0, iel, _solType[k]);
              _solution[igridn]->_Bdc[k]->set(idof, 1.);
            }
          }
        }
        if( _fixSolutionAtOnePoint[k] == true  && igridn == 0 && _iproc == 0 ){
          _solution[igridn]->_Bdc[k]->set(0, 0.);
          _solution[igridn]->_Sol[k]->set(0, 0.);
        }
        _solution[igridn]->_Sol[k]->close();
        _solution[igridn]->_Bdc[k]->close();
      }
    }


  }

  void MultiLevelSolution::SaveSolution(const char* filename, const double time) {

    char composedFileName[100];

    for(int i = 0; i < _solName.size(); i++) {
      sprintf(composedFileName, "./save/%s_time%f_sol%s_level%d", filename, time, _solName[i], _gridn);
      _solution[_gridn - 1]->_Sol[i]->BinaryPrint(composedFileName);
    }
  }

  void MultiLevelSolution::LoadSolution(const char* filename) {
    LoadSolution(_gridn, filename);
  }

  void MultiLevelSolution::LoadSolution(const unsigned &level, const char* filename) {

    if(level > _gridn){
      std::cout<< "Error in MultiLevelSolution::LoadSolution function:"<<std::endl;
      std::cout<< "the solution level = "<<level<<" is not available in this MultilevelSolution"<<std::endl;
      abort();
    }

    char composedFileName[200];
    for(int i = 0; i < _solName.size(); i++) {
      sprintf(composedFileName, "%s_sol%s_level%d", filename, _solName[i], level);
      // check if the file really exists
      if ( strncmp(filename, "http://", 7) && strncmp(filename, "ftp://", 6) ){
        struct stat buffer;
        if(stat (composedFileName, &buffer) != 0) {
          std::cerr << "Error: cannot locate file " << composedFileName << std::endl;
          abort();
        }
      }
      _solution[level - 1]->_Sol[i]->BinaryLoad(composedFileName);
    }

    for(int gridf = level; gridf < _gridn; gridf++) {
      for(unsigned i = 0; i < _solName.size(); i++) {
        _solution[gridf]->_Sol[i]->matrix_mult(*_solution[gridf - 1]->_Sol[i],
                                               *_mlMesh->GetLevel(gridf)->GetCoarseToFineProjection(_solType[i]));
        _solution[gridf]->_Sol[i]->close();
      }
    }


  }


} //end namespace femus
