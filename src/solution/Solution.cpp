/*=========================================================================

 Program: FEMUS
 Module: Solution
 Authors: Eugenio Aulisa

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

#include <ctime>
#include <fstream>
#include <algorithm>
#include "Solution.hpp"
#include "FemusDefault.hpp"
#include "ElemType.hpp"
#include "ParalleltypeEnum.hpp"
#include "NumericVector.hpp"


namespace femus {

  using std::cout;
  using std::endl;

  /**
   *  Constructor
   **/
// ------------------------------------------------------------------
  Solution::Solution(Mesh *other_msh) {
    _msh = other_msh;

    for(int i = 0; i < 5; i++) {
      _GradMat[i].resize(_msh->GetDimension());
      _AMR_flag = 0;
    }
    _FSI = false;
  }

  /**
   * Destructor
   **/
// ------------------------------------------------------------------
  Solution::~Solution() {
    for(unsigned i = 0; i < _SolName.size(); i++) {
      delete [] _SolName[i];
    }
  }

  /**
   * Add a new variable called 'name'
   */
  void Solution::AddSolution(const char name[], const FEFamily fefamily, const FEOrder order,
                             const unsigned& tmorder, const bool &Pde_type) {

    unsigned n = _Sol.size();

    _SolType.resize(n + 1u);
    _SolName.resize(n + 1u);
    _SolTmOrder.resize(n + 1u);
    _family.resize(n + 1u);
    _order.resize(n + 1u);

    _Sol.resize(n + 1u);
    _Sol[n] = NULL;

    _Res.resize(n + 1u);
    _Res[n] = NULL;

    _Eps.resize(n + 1u);
    _Eps[n] = NULL;

    _GradVec.resize(n + 1u);
    _GradVec[n].resize(_msh->GetDimension());

    for(int i = 0; i < _msh->GetDimension(); i++) {
      _GradVec[n][i] = NULL;
    }

    _Bdc.resize(n + 1u);
    _Bdc[n] = NULL;
    _ResEpsBdcFlag.resize(n + 1u);
    _ResEpsBdcFlag[n] = Pde_type;

    _family[n] = fefamily;
    _order[n] = order;
    _SolType[n] = order - ((fefamily == LAGRANGE) ? 1 : 0) + fefamily * 3;
    _SolTmOrder[n] = tmorder;
    _SolOld.resize(n + 1u);
    _SolOld[n] = NULL;
    _SolName[n] = new char [DEFAULT_SOL_NCHARS];

    _removeNullSpace.resize(n + 1u);
    _removeNullSpace[n] = false;

    strcpy(_SolName[n], name);

  }


  
  void Solution::AddSolution_par(const int n_sols, const char name[], const FEFamily fefamily, const FEOrder order,
                             const unsigned& tmorder, const bool &Pde_type) {

   const unsigned old_size = _Sol.size();
   const int      new_size = old_size + n_sols;
    
    ResizeSolution_par(new_size);
    
    //initialize
    if ( n_sols > 0 ) {
        
   for(int s = 0; s < n_sols; s++) {
       
       _Sol[old_size + s] = NULL;
    _SolOld[old_size + s] = NULL;
       _Res[old_size + s] = NULL;
       _Eps[old_size + s] = NULL; 
       _Bdc[old_size + s] = NULL;       
       
     _GradVec[old_size + s].resize(_msh->GetDimension());
         for(int i = 0; i < _msh->GetDimension(); i++) {     _GradVec[old_size + s][i] = NULL; }  
         
  _ResEpsBdcFlag[old_size + s] = Pde_type;       
         _family[old_size + s] = fefamily;
          _order[old_size + s] = order;
        _SolType[old_size + s] = order - ((fefamily == LAGRANGE) ? 1 : 0) + fefamily * 3;
     _SolTmOrder[old_size + s] = tmorder;
        _SolName[old_size + s] = new char [DEFAULT_SOL_NCHARS];
 strcpy(_SolName[old_size + s], name);
_removeNullSpace[old_size + s] = false;
    
      }

    }   
   
   
  }  
  


  void Solution::ResizeSolution_par(const int new_size) {
      
// it seems that the destructors of the NumericVectors are not called after a shrink...      
      
         _SolType.resize(new_size);
         _SolName.resize(new_size);
      _SolTmOrder.resize(new_size);
          _family.resize(new_size);
           _order.resize(new_size);
   _ResEpsBdcFlag.resize(new_size);        
 _removeNullSpace.resize(new_size);
           
// NumericVector objects
             _Sol.resize(new_size); 
          _SolOld.resize(new_size); 
         _GradVec.resize(new_size); 
             _Res.resize(new_size); 
             _Eps.resize(new_size); 
             _Bdc.resize(new_size);

      
  }
  

  
  /**
   * Get the solution index for the variable called name
   **/
//-------------------------------------------------------------------
  unsigned Solution::GetIndex(const char name[]) const {
    unsigned index = 0;

    while(strcmp(_SolName[index], name)) {
      index++;

      if(index == _Res.size()) {
        cout << "error! invalid name entry GetIndex(...)" << endl;
        exit(0);
      }
    }

    return index;
  }

  /**
   * Allocate memory for the variable called name
   **/
// ------------------------------------------------------------------
  void Solution::ResizeSolutionVector(const char name[]) {

    unsigned i = GetIndex(name);

    if(_Sol[i])  delete _Sol[i];

    if(_ResEpsBdcFlag[i]) {
      if(_Res[i]) delete _Res[i];

      if(_Eps[i]) delete _Eps[i];

      if(_Bdc[i]) delete _Bdc[i];
    }

    if(_SolTmOrder[i] == 2) {
      if(_SolOld[i]) delete _SolOld[i];
    }

    _Sol[i] = NumericVector::build().release();

    if(n_processors() == 1) {  // IF SERIAL
      _Sol[i]->init(_msh->_dofOffset[_SolType[i]][n_processors()], _msh->_ownSize[_SolType[i]][processor_id()], false, SERIAL);
    }
    else { // IF PARALLEL
      if(_SolType[i] < 3) {
        if(_msh->_ghostDofs[_SolType[i]][processor_id()].size() != 0) {
          _Sol[i]->init(_msh->_dofOffset[_SolType[i]][n_processors()], _msh->_ownSize[_SolType[i]][processor_id()],
                        _msh->_ghostDofs[_SolType[i]][processor_id()], false, GHOSTED);
        }
        else {
          std::vector <int> fake_ghost(1, _msh->_ownSize[_SolType[i]][processor_id()]);
          _Sol[i]->init(_msh->_dofOffset[_SolType[i]][n_processors()], _msh->_ownSize[_SolType[i]][processor_id()],
                        fake_ghost, false, GHOSTED);
        }
      }
      else { //discontinuous pressure has no ghost nodes
        _Sol[i]->init(_msh->_dofOffset[_SolType[i]][n_processors()], _msh->_ownSize[_SolType[i]][processor_id()], false, PARALLEL);
      }
    }

    if(_SolTmOrder[i] == 2) {  // only if the variable is time dependent
      _SolOld[i] = NumericVector::build().release();
      _SolOld[i]->init(*_Sol[i]);
    }

    if(_ResEpsBdcFlag[i]) {  //only if the variable is a Pde type

      _Res[i] = NumericVector::build().release();
      _Res[i]->init(*_Sol[i]);

      _Eps[i] = NumericVector::build().release();
      _Eps[i]->init(*_Sol[i]);

      _Bdc[i] = NumericVector::build().release();
      _Bdc[i]->init(*_Sol[i]);
    }
  }


  /** Init and set to zero The AMR Eps vector */
  void Solution::InitAMREps() {
    _AMR_flag = 1;
    _AMREps.resize(_Sol.size());

    for(int i = 0; i < _Sol.size(); i++) {
      _AMREps[i] = NumericVector::build().release();
      _AMREps[i]->init(*_Sol[i]);
      _AMREps[i]->zero();
    }

  }



  /**
   * Deallocate memory for all variables
   **/
// ------------------------------------------------------------------
  void Solution::FreeSolutionVectors() {
    for(unsigned i = 0; i < _Sol.size(); i++) {
      if(_Sol[i]) delete _Sol[i];

      _Sol[i] = NULL;

      if(_ResEpsBdcFlag[i]) {
        if(_Res[i]) delete _Res[i];

        _Res[i] = NULL;

        if(_Eps[i]) delete _Eps[i];

        _Eps[i] = NULL;

        if(_Bdc[i]) delete _Bdc[i];

        _Bdc[i] = NULL;

      }

      if(_SolTmOrder[i] == 2) {
        if(_SolOld[i]) delete _SolOld[i];

        _SolOld[i] = NULL;
      }

      for(int j = 0; j < _msh->GetDimension(); j++) {
        if(_GradVec[i][j]) {
          if(_GradVec[i][j]) delete _GradVec[i][j];

          _GradVec[i][j] = NULL;
        }
      }

      if(_AMR_flag) {
        if(_AMREps[i]) delete _AMREps[i];

        _AMREps[i] = NULL;
      }
    }

    for(unsigned i = 0; i < 5; i++) {
      for(int j = 0; j < _msh->GetDimension(); j++) {
        if(_GradMat[i][j]) {
          delete _GradMat[i][j];
        }
      }
    }
  }

  /**
   * Update _Sol, _Res and _Eps based on EPS and RES
   **/
// //--------------------------------------------------------------------------------
//   void Solution::UpdateSolAndRes(const vector <unsigned> &_SolPdeIndex,  NumericVector* _EPS,  NumericVector* _RES,
//                                  const vector <vector <unsigned> > &KKoffset) {
// 
//     PetscScalar zero = 0.;
// 
//     for(unsigned k = 0; k < _SolPdeIndex.size(); k++) {
//       unsigned indexSol = _SolPdeIndex[k];
//       unsigned soltype =  _SolType[indexSol];
// 
//       int loc_offset_EPS = KKoffset[k][processor_id()];
// 
//       int glob_offset_eps = _msh->_dofOffset[soltype][processor_id()];
// 
//       vector <int> index(_msh->_ownSize[soltype][processor_id()]);
// 
//       for(int i = 0; i < _msh->_ownSize[soltype][processor_id()]; i++) {
//         index[i] = loc_offset_EPS + i;
//       }
// 
//       vector <double> valueEPS(_msh->_ownSize[soltype][processor_id()]);
//       _EPS->get(index, valueEPS);
//       vector <double> valueRES(_msh->_ownSize[soltype][processor_id()]);
//       _RES->get(index, valueRES);
// 
//       for(int i = 0; i < _msh->_ownSize[soltype][processor_id()]; i++) {
//         _Eps[indexSol]->set(i + glob_offset_eps, valueEPS[i]);
// 
//         if((*_Bdc[indexSol])(i + glob_offset_eps) > 1.1) _Res[indexSol]->set(i + glob_offset_eps, valueRES[i]);
//         else _Res[indexSol]->set(i + glob_offset_eps, zero);
//       }
// 
//       _Res[indexSol]->close();
//       _Eps[indexSol]->close();
//     }
// 
//     for(unsigned k = 0; k < _SolPdeIndex.size(); k++) {
//       unsigned indexSol = _SolPdeIndex[k];
//       _Sol[indexSol]->add(*_Eps[indexSol]);
//       _Sol[indexSol]->close();
// 
//       if(_AMR_flag) {
//         _AMREps[indexSol]->add(*_Eps[indexSol]);
//         _AMREps[indexSol]->close();
//       }
// 
//     }
// 
//   }

  /**
   * Update _Sol
   **/

  void Solution::UpdateSol(const vector <unsigned> &_SolPdeIndex,  NumericVector* _EPS, const vector <vector <unsigned> > &KKoffset) {

    PetscScalar zero = 0.;

    for(unsigned k = 0; k < _SolPdeIndex.size(); k++) {
      unsigned indexSol = _SolPdeIndex[k];
      unsigned soltype =  _SolType[indexSol];

      int loc_offset_EPS = KKoffset[k][processor_id()];

      int glob_offset_eps = _msh->_dofOffset[soltype][processor_id()];

      vector <int> index(_msh->_ownSize[soltype][processor_id()]);

      for(int i = 0; i < _msh->_ownSize[soltype][processor_id()]; i++) {
        index[i] = loc_offset_EPS + i;
      }

      vector <double> valueEPS(_msh->_ownSize[soltype][processor_id()]);
      _EPS->get(index, valueEPS);
      //vector <double> valueRES(_msh->_ownSize[soltype][processor_id()]);
      //_RES->get(index,valueRES);

      for(int i = 0; i < _msh->_ownSize[soltype][processor_id()]; i++) {
        _Eps[indexSol]->set(i + glob_offset_eps, valueEPS[i]);
        //if ((*_Bdc[indexSol])(i+glob_offset_eps)>1.1) _Res[indexSol]->set(i+glob_offset_eps,valueRES[i]);
        //else _Res[indexSol]->set(i+glob_offset_eps,zero);
      }

      //_Res[indexSol]->close();
      _Eps[indexSol]->close();
    }

    for(unsigned k = 0; k < _SolPdeIndex.size(); k++) {
      unsigned indexSol = _SolPdeIndex[k];
      _Sol[indexSol]->add(*_Eps[indexSol]);
      _Sol[indexSol]->close();

      if(_AMR_flag) {
        _AMREps[indexSol]->add(*_Eps[indexSol]);
        _AMREps[indexSol]->close();
      }

    }

  }


  /**
   * Update _Res
   **/
//--------------------------------------------------------------------------------
  void Solution::UpdateRes(const vector <unsigned> &_SolPdeIndex, NumericVector* _RES, const vector <vector <unsigned> > &KKoffset) {

    PetscScalar zero = 0.;

    for(unsigned k = 0; k < _SolPdeIndex.size(); k++) {
      unsigned indexSol = _SolPdeIndex[k];
      unsigned soltype =  _SolType[indexSol];

      int loc_offset_RES = KKoffset[k][processor_id()];

      int glob_offset_res = _msh->_dofOffset[soltype][processor_id()];

      vector <int> index(_msh->_ownSize[soltype][processor_id()]);

      for(int i = 0; i < _msh->_ownSize[soltype][processor_id()]; i++) {
        index[i] = loc_offset_RES + i;
      }

      vector <double> valueRES(_msh->_ownSize[soltype][processor_id()]);
      _RES->get(index, valueRES);

      for(int i = 0; i < _msh->_ownSize[soltype][processor_id()]; i++) {
        if((*_Bdc[indexSol])(i + glob_offset_res) > 1.1) {
          _Res[indexSol]->set(i + glob_offset_res, valueRES[i]);
        }
        else {
          _Res[indexSol]->set(i + glob_offset_res, zero);
        }
      }

      _Res[indexSol]->close();
    }

  }

  bool Solution::FlagAMRRegionBasedOnErroNorm(const vector <unsigned> &solIndex, std::vector <double> &AMRthreshold, const unsigned& normType) {

    unsigned    iproc = _msh->processor_id(); // get the process_id (for parallel computation)
    const unsigned  dim = _msh->GetDimension();

    Solution* AMR = _msh->_topology;
    unsigned  AMRIndex = AMR->GetIndex("AMR");
    AMR->_Sol[AMRIndex]->zero();

    NumericVector *counter_vec;
    counter_vec = NumericVector::build().release();
    counter_vec->init(_msh->n_processors(), 1 , false, AUTOMATIC);
    counter_vec->zero();

    if(AMRthreshold.size() != solIndex.size()) {
      double value = AMRthreshold[0];
      AMRthreshold.assign(solIndex.size(), value);
    }

    for(unsigned k = 0; k < solIndex.size(); k++) {

      vector < double >  sol; // local solution
      unsigned solType = _SolType[solIndex[k]];;    // get the finite element type for "u"

      vector < vector < double > > x(dim);    // local coordinates
      unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

      vector <double> phi;  // local test function
      vector <double> phi_x; // local test function first order partial derivatives
      vector <double> phi_xx; // local test function second order partial derivatives
      double weight; // gauss point weight

      // reserve memory for the local standar vectors
      const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
      sol.reserve(maxSize);

      for(unsigned i = 0; i < dim; i++) {
        x[i].reserve(maxSize);
      }

      phi.reserve(maxSize);
      phi_x.reserve(maxSize * dim);
      unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
      phi_xx.reserve(maxSize * dim2);

      double solNorm2 = 0.;
      double volume = 0.;

      for(int iel = _msh->_elementOffset[iproc]; iel < _msh->_elementOffset[iproc + 1]; iel++) {

        short unsigned ielGeom = _msh->GetElementType(iel);
        unsigned solDofs  = _msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
        unsigned xDofs  = _msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

        // resize local arrays
        sol.resize(solDofs);

        for(int i = 0; i < dim; i++) {
          x[i].resize(xDofs);
        }

        // local storage of global mapping and solution
        for(unsigned i = 0; i < solDofs; i++) {
          unsigned iDof = _msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
          sol[i] = (*_Sol[ solIndex[k]])(iDof);      // global extraction and local storage for the solution
        }

        // local storage of coordinates
        for(unsigned i = 0; i < xDofs; i++) {
          unsigned iDof  = _msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

          for(unsigned j = 0; j < dim; j++) {
            x[j][i] = (*_msh->_topology->_Sol[j])(iDof);      // global extraction and local storage for the element coordinates
          }
        }

        // *** Gauss point loop ***
        for(unsigned ig = 0; ig < _msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
          // *** get gauss point weight, test function and test function partial derivatives ***
          _msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

          double solig = 0.;
          std::vector < double > solGradig(dim, 0.);

          for(unsigned i = 0; i < solDofs; i++) {
            solig += phi[i] * sol[i];

            if(normType > 0) {
              for(int j = 0; j < dim; j++) {
                solGradig[j] += sol[i] * phi_x[i * dim + j];
              }
            }
          }

          solNorm2 += solig * solig * weight;

          if(normType > 0) {
            for(int j = 0; j < dim; j++) {
              solNorm2 += solGradig[j] * solGradig[j] * weight;
            }
          }

          volume += weight;
        }
      }

      NumericVector* parallelVec;
      parallelVec = NumericVector::build().release();
      parallelVec->init(_msh->n_processors(), 1 , false, AUTOMATIC);

      parallelVec->set(iproc, solNorm2);
      parallelVec->close();
      solNorm2 = parallelVec->l1_norm();

      parallelVec->set(iproc, volume);
      parallelVec->close();
      volume = parallelVec->l1_norm();

      for(int iel = _msh->_elementOffset[iproc]; iel < _msh->_elementOffset[iproc + 1]; iel++) {
        if(_msh->el->GetIfElementCanBeRefined(iel) && (*AMR->_Sol[AMRIndex])(iel) == 0.) {

          double ielErrNorm2 = 0.;
          double ielVolume = 0.;

          short unsigned ielGeom = _msh->GetElementType(iel);
          unsigned solDofs  = _msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
          unsigned xDofs  = _msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

          // resize local arrays
          sol.resize(solDofs);

          for(int i = 0; i < dim; i++) {
            x[i].resize(xDofs);
          }

          // local storage of global mapping and solution
          for(unsigned i = 0; i < solDofs; i++) {
            unsigned iDof = _msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
            sol[i] = (*_AMREps[ solIndex[k]])(iDof);      // global extraction and local storage for the solution
          }

          // local storage of coordinates
          for(unsigned i = 0; i < xDofs; i++) {
            unsigned iDof  = _msh->GetSolutionDof(i, iel, xType); // global to global mapping between coordinates node and coordinate dof

            for(unsigned j = 0; j < dim; j++) {
              x[j][i] = (*_msh->_topology->_Sol[j])(iDof); // global extraction and local storage for the element coordinates
            }
          }

          for(unsigned ig = 0; ig < _msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
            _msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

            double solig = 0.;

            std::vector < double > solGradig(dim, 0.);

            for(unsigned i = 0; i < solDofs; i++) {
              solig += phi[i] * sol[i];

              if(normType) {
                for(int j = 0; j < dim; j++) {
                  solGradig[j] += sol[i] * phi_x[i * dim + j];
                }
              }
            }

            ielErrNorm2 += solig * solig * weight;

            if(normType > 0) {
              for(int j = 0; j < dim; j++) {
                ielErrNorm2 += solGradig[j] * solGradig[j] * weight;
              }
            }

            ielVolume += weight;
          }

          if(ielErrNorm2 > AMRthreshold[k] * AMRthreshold[k] * solNorm2 * ielVolume / volume) {
            AMR->_Sol[AMRIndex]->set(iel, 1.);
            counter_vec->add(_iproc, 1.);

// 	    for(unsigned i = 0; i < _msh->GetElementDofNumber(iel, 0); i++) { //loop on the element vertices
// 	      unsigned inode = _msh->el->GetElementDofIndex(iel, i);
// 	      const std::vector < unsigned > & localElementNearVertexNumber = _msh->el->GetLocalElementNearVertex(inode);
// 	      unsigned nve = localElementNearVertexNumber.size();
// 	      for(unsigned j = 0; j < nve; j++) {
// 		unsigned jel = localElementNearVertexNumber[j];
// 		if(jel >= _msh->_elementOffset[iproc] && jel < _msh->_elementOffset[iproc + 1]) {
// 		  if(_msh->el->GetIfElementCanBeRefined(jel) && (*AMR->_Sol[AMRIndex])(jel) == 0.) {
// 		    AMR->_Sol[AMRIndex]->set(jel, 1.);
// 		    counter_vec->add(_iproc, 1.);
// 		  }
// 		}
//               }
//             }
          }
        }
      }
    }

    AMR->_Sol[AMRIndex]->close();

    counter_vec->close();
    double counter = counter_vec->l1_norm();
    bool test = (counter <= _nprocs) ? 1 : 0;

    delete counter_vec;
    return test;
  }



  bool Solution::FlagAMRRegionBasedOnErroNormAdaptive(const vector <unsigned> &solIndex, std::vector <double> &AMRthreshold, 
						      const unsigned& normType, const double &neighborThresholdValue) {

    const double scale2[3][2] = {{0.111111, 1.}, {0.0204081632653, 0.111111}, {0.0204081632653, 0.111111} };
    //const double scale2[3][2] = {{1., 1.}, {1., 1.}, {1., 1.} };


    unsigned    iproc = _msh->processor_id(); // get the process_id (for parallel computation)
    const unsigned  dim = _msh->GetDimension();

    Solution* AMR = _msh->_topology;
    unsigned  AMRIndex = AMR->GetIndex("AMR");
    AMR->_Sol[AMRIndex]->zero();

    if(AMRthreshold.size() != solIndex.size()) {
      double value = AMRthreshold[0];
      AMRthreshold.assign(solIndex.size(), value);
    }

    for(unsigned k = 0; k < solIndex.size(); k++) {

      vector < double >  sol; // local solution
      unsigned solType = _SolType[solIndex[k]];;    // get the finite element type for "u"

      vector < vector < double > > x(dim);    // local coordinates
      unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

      vector <double> phi;  // local test function
      vector <double> phi_x; // local test function first order partial derivatives
      vector <double> phi_xx; // local test function second order partial derivatives
      double weight; // gauss point weight

      // reserve memory for the local standar vectors
      const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
      sol.reserve(maxSize);

      for(unsigned i = 0; i < dim; i++) {
        x[i].reserve(maxSize);
      }

      phi.reserve(maxSize);
      phi_x.reserve(maxSize * dim);
      unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
      phi_xx.reserve(maxSize * dim2);

      double solNorm2 = 0.;
      double volumeRefined = 0.;
      double volume = 0.;

      for(int iel = _msh->_elementOffset[iproc]; iel < _msh->_elementOffset[iproc + 1]; iel++) {

        short unsigned ielGeom = _msh->GetElementType(iel);
        unsigned solDofs  = _msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
        unsigned xDofs  = _msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

        // resize local arrays
        sol.resize(solDofs);

        for(int i = 0; i < dim; i++) {
          x[i].resize(xDofs);
        }

        // local storage of global mapping and solution
        for(unsigned i = 0; i < solDofs; i++) {
          unsigned iDof = _msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
          sol[i] = (*_Sol[ solIndex[k]])(iDof);      // global extraction and local storage for the solution
        }

        // local storage of coordinates
        for(unsigned i = 0; i < xDofs; i++) {
          unsigned iDof  = _msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

          for(unsigned j = 0; j < dim; j++) {
            x[j][i] = (*_msh->_topology->_Sol[j])(iDof);      // global extraction and local storage for the element coordinates
          }
        }

        // *** Gauss point loop ***
        for(unsigned ig = 0; ig < _msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
          // *** get gauss point weight, test function and test function partial derivatives ***
          _msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

          double solig = 0.;
          std::vector < double > solGradig(dim, 0.);

          for(unsigned i = 0; i < solDofs; i++) {
            solig += phi[i] * sol[i];

            if(normType > 0) {
              for(int j = 0; j < dim; j++) {
                solGradig[j] += sol[i] * phi_x[i * dim + j];
              }
            }
          }

          solNorm2 += solig * solig * weight;

          if(normType > 0) {
            for(int j = 0; j < dim; j++) {
              solNorm2 += solGradig[j] * solGradig[j] * weight;
            }
          }

          volume += weight;

          if(_msh->el->GetIfElementCanBeRefined(iel)) {
            volumeRefined += weight;
          }
        }
      }

      NumericVector* parallelVec;
      parallelVec = NumericVector::build().release();
      parallelVec->init(_msh->n_processors(), 1 , false, AUTOMATIC);

      parallelVec->set(iproc, solNorm2);
      parallelVec->close();
      solNorm2 = parallelVec->l1_norm();

      parallelVec->set(iproc, volume);
      parallelVec->close();
      volume = parallelVec->l1_norm();

      parallelVec->set(iproc, volumeRefined);
      parallelVec->close();
      volumeRefined = parallelVec->l1_norm();

      double  volumeTestFalse = 0.;
      double errTestTrue2 = 0.;

      double eps2 = AMRthreshold[k] * AMRthreshold[k] * solNorm2  / volume;

      unsigned offset = _msh->_elementOffset[iproc];
      std::vector < double > ielVolume(_msh->_elementOffset[iproc + 1] - _msh->_elementOffset[iproc],0);
      std::vector < double > ielErrNorm2(_msh->_elementOffset[iproc + 1] - _msh->_elementOffset[iproc],0);
      
      
      for(int iel = _msh->_elementOffset[iproc]; iel < _msh->_elementOffset[iproc + 1]; iel++) {
        if(_msh->el->GetIfElementCanBeRefined(iel)) {

//           double ielErrNorm2 = 0.;
//           double ielVolume = 0.;

          short unsigned ielGeom = _msh->GetElementType(iel);
          unsigned solDofs  = _msh->GetElementDofNumber(iel, solType);    // number of solution element dofs
          unsigned xDofs  = _msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

          // resize local arrays
          sol.resize(solDofs);

          for(int i = 0; i < dim; i++) {
            x[i].resize(xDofs);
          }

          // local storage of global mapping and solution
          for(unsigned i = 0; i < solDofs; i++) {
            unsigned iDof = _msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
            sol[i] = (*_AMREps[ solIndex[k]])(iDof);      // global extraction and local storage for the solution
          }

          // local storage of coordinates
          for(unsigned i = 0; i < xDofs; i++) {
            unsigned iDof  = _msh->GetSolutionDof(i, iel, xType); // global to global mapping between coordinates node and coordinate dof

            for(unsigned j = 0; j < dim; j++) {
              x[j][i] = (*_msh->_topology->_Sol[j])(iDof); // global extraction and local storage for the element coordinates
            }
          }

          for(unsigned ig = 0; ig < _msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++) {
            _msh->_finiteElement[ielGeom][solType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

            double errig = 0.;

            std::vector < double > errGradig(dim, 0.);

            for(unsigned i = 0; i < solDofs; i++) {
              errig += phi[i] * sol[i];

              if(normType) {
                for(int j = 0; j < dim; j++) {
                  errGradig[j] += sol[i] * phi_x[i * dim + j];
                }
              }
            }

            ielErrNorm2[iel-offset] += scale2[solType][normType] * errig * errig * weight;

            if(normType > 0) {
              for(int j = 0; j < dim; j++) {
                ielErrNorm2[iel-offset] += scale2[solType][normType] * errGradig[j] * errGradig[j] * weight;
              }
            }

            ielVolume[iel-offset] += weight;
          }

          if(ielErrNorm2[iel-offset] > eps2 * ielVolume[iel-offset]  || 
	    ( (*AMR->_Sol[AMRIndex])(iel) == 2. && ielErrNorm2[iel-offset] > neighborThresholdValue * eps2 * ielVolume[iel-offset] ) ) {
            AMR->_Sol[AMRIndex]->set(iel, 1.);
            volumeTestFalse += ielVolume[iel-offset];

            if( ielErrNorm2[iel-offset] > eps2 * ielVolume[iel-offset] ) {
              for(unsigned j = 1; j < _msh->el->GetElementNearElementSize(iel,1);j++){
		unsigned jel = _msh->el->GetElementNearElement(iel,j);
		if(jel >= _msh->_elementOffset[iproc] && jel<_msh->_elementOffset[iproc + 1] ){
		  if(_msh->el->GetIfElementCanBeRefined(jel)){
		    if(jel > iel) {
                      AMR->_Sol[AMRIndex]->set(jel, 2.);
                    }
                    else if( (*AMR->_Sol[AMRIndex])(jel) == 0. &&  ielErrNorm2[jel-offset] > neighborThresholdValue * eps2 * ielVolume[jel-offset] ) {
		      errTestTrue2 -= ielErrNorm2[jel-offset];
		      AMR->_Sol[AMRIndex]->set(jel, 1.);
		      volumeTestFalse += ielVolume[jel-offset];
		    }
                  }
                }
              }
            }

          }
          else {
	    AMR->_Sol[AMRIndex]->set(iel, 0.);
            errTestTrue2 += ielErrNorm2[iel-offset];
          }
        }
      }

      parallelVec->set(iproc, volumeTestFalse);
      parallelVec->close();
      volumeTestFalse = parallelVec->l1_norm();

      parallelVec->set(iproc, errTestTrue2);
      parallelVec->close();
      errTestTrue2 = parallelVec->l1_norm();

      if(volumeTestFalse != 0) {
	cout.precision(24);
	printf("%e %e %e %e %e \n",errTestTrue2, solNorm2, volume, volumeRefined, volumeTestFalse);
        //std::cout  << errTestTrue2 << " " << solNorm2 << " " << volume << " " << volumeRefined << " " << volumeTestFalse << std::endl;
	AMRthreshold[k] = sqrt(AMRthreshold[k] * AMRthreshold[k] * volumeRefined / volumeTestFalse - errTestTrue2 / solNorm2 * volume / volumeTestFalse);
        std::cout << AMRthreshold[k] << std::endl;  
      }
      else {
        AMRthreshold[k] = 1.;
      }

      std::cout << "\nNew AMR threshold value =" << AMRthreshold[k] << std::endl;
      
      delete parallelVec;
    }

    AMR->_Sol[AMRIndex]->close();

    double counter = AMR->_Sol[AMRIndex]->l1_norm();
    bool test = (counter * _msh->GetRefIndex() <= _nprocs) ? 1 : 0;

    return test;
  }

// bool Solution::FlagAMRRegionBasedOnl2(const vector <unsigned> &SolIndex,const double &AMRthreshold){
//
//   vector <double> SolMax(SolIndex.size());
//   vector <unsigned> SolType(SolIndex.size());
//   vector <unsigned> SolEndInd(SolIndex.size());
//
//   unsigned END_IND[5]= {0,1,1,4,5};
//
//   for (unsigned k=0; k<SolIndex.size(); k++) {
//     double EPSMAX = _AMREps[SolIndex[k]]->linfty_norm ();
//     double SOLMAX = _Sol[SolIndex[k]]->linfty_norm ();
//     cout << std::endl << "Current maximum relative change = " <<EPSMAX/SOLMAX << endl << endl;
//     SolMax[k] = AMRthreshold * SOLMAX;
//     SolType[k] = _SolType[SolIndex[k]];
//     SolEndInd[k]   = END_IND[SolType[k]];
//   }
//
//   Solution* AMR = _msh->_topology;
//   unsigned  AMRIndex= AMR->GetIndex("AMR");
//   AMR->_Sol[AMRIndex]->zero();
//
//   unsigned nel= _msh->GetNumberOfElements();
//
//   NumericVector *counter_vec;
//   counter_vec = NumericVector::build().release();
//
//   if(_nprocs==1) {
//     counter_vec->init(_nprocs,1,false,SERIAL);
//   }
//   else {
//     counter_vec->init(_nprocs,1,false,PARALLEL);
//   }
//   counter_vec->zero();
//
//   for (int kel = _msh->_elementOffset[_iproc]; kel < _msh->_elementOffset[_iproc+1]; kel++) {
//     short unsigned kelt=_msh->GetElementType(kel);
//     for (unsigned k=0; k<SolIndex.size(); k++) {
//       if(SolType[k]<3){
//         unsigned nve=_msh->GetElementDofNumber(kel,SolEndInd[k]);
// 	for(unsigned i=0; i<nve; i++) {
// 	  unsigned inode_metis=_msh->GetSolutionDof(i,kel,SolType[k]);
// 	  double value = (*_AMREps[SolIndex[k]])(inode_metis);
// 	  if(fabs(value)>SolMax[k]){
// 	    counter_vec->add(_iproc,1.);
// 	    AMR->_Sol[AMRIndex]->set(kel, 1.);
// 	    k=SolIndex.size();
// 	    i=nve;
// 	  }
// 	}
//       }
//     }
//   }
//   AMR->_Sol[AMRIndex]->close();
//
//   counter_vec->close();
//   double counter=counter_vec->l1_norm();
//   bool test=(counter<=_nprocs)?1:0;
//
//   delete counter_vec;
//   return test;
//
// }



//   bool Solution::FlagAMRRegionBasedOnSemiNorm(const vector <unsigned> &SolIndex, const double & AMRthreshold) {
//
//     vector <double> GradSolMax(SolIndex.size());
//     vector <unsigned> SolType(SolIndex.size());
//     vector <unsigned> SolEndInd(SolIndex.size());
//     unsigned dim = _msh->GetDimension();
//
//     unsigned nel = _msh->GetNumberOfElements();
//
//     for(unsigned k = 0; k < SolIndex.size(); k++) {
//
//       if(SolType[k] < 3) {
//
//         SolType[k] = _SolType[SolIndex[k]];
//
//         BuildGradMatrixStructure(SolType[k]);
//
//         GradSolMax[k] = 0.;
//
//         for(int i = 0; i < dim; i++) {
//
//           if(_GradVec[SolIndex[k]][i] == 0) {
//             _GradVec[SolIndex[k]][i] = NumericVector::build().release();
//
//             if(n_processors() == 1) {  // IF SERIAL
//               _GradVec[SolIndex[k]][i]->init(_msh->_dofOffset[3][n_processors()], _msh->_ownSize[3][processor_id()], false, SERIAL);
//             }
//             else { //discontinuous pressure has no ghost nodes
//               _GradVec[SolIndex[k]][i]->init(_msh->_dofOffset[3][n_processors()], _msh->_ownSize[3][processor_id()], false, PARALLEL);
//
//             }
//           }
//
//           _GradVec[SolIndex[k]][i]->matrix_mult(*_Sol[SolIndex[k]], *_GradMat[SolType[k]][i]);
//           double GradSolMaxi = _GradVec[SolIndex[k]][i]->linfty_norm();
//           GradSolMax[k] += GradSolMaxi * GradSolMaxi;
//           _GradVec[SolIndex[k]][i]->close();
//           _GradVec[SolIndex[k]][i]->matrix_mult(*_Eps[SolIndex[k]], *_GradMat[SolType[k]][i]);
//         }
//
//         GradSolMax[k] = AMRthreshold * sqrt(GradSolMax[k]);
//       }
//     }
//
//     Solution* AMR = _msh->_topology;
//     unsigned  AMRIndex = AMR->GetIndex("AMR");
//
//     AMR->_Sol[AMRIndex]->zero();
//
//     NumericVector *counter_vec;
//     counter_vec = NumericVector::build().release();
//
//     if(_nprocs == 1) {
//       counter_vec->init(_nprocs, 1, false, SERIAL);
//     }
//     else {
//       counter_vec->init(_nprocs, 1, false, PARALLEL);
//     }
//
//     counter_vec->zero();
//
//     for(int iel_metis = _msh->_elementOffset[_iproc]; iel_metis < _msh->_elementOffset[_iproc + 1]; iel_metis++) {
//
//       for(unsigned k = 0; k < SolIndex.size(); k++) {
//
//         if(SolType[k] < 3) {
//           double value = 0.;
//
//           for(int i = 0; i < dim; i++) {
//             double valuei = (*_GradVec[SolIndex[k]][i])(iel_metis);
//             value += valuei * valuei;
//           }
//
//           value = sqrt(value);
//
//           if(fabs(value) > GradSolMax[k]) {
//             counter_vec->add(_iproc, 1.);
//             AMR->_Sol[AMRIndex]->set(iel_metis, 1.);
//             k = SolIndex.size();
//           }
//         }
//       }
//     }
//
//     AMR->_Sol[AMRIndex]->close();
//
//     counter_vec->close();
//     double counter = counter_vec->l1_norm();
//     bool test = (counter <= _nprocs) ? 1 : 0;
//
//     return test;
//   }


  void Solution::BuildGradMatrixStructure(unsigned SolType) {

    if(SolType < 3 && _GradMat[SolType][0] == 0) {

      unsigned dim = _msh->GetDimension();

      int nr     = _msh->_dofOffset[3][_nprocs];
      int nc     = _msh->_dofOffset[SolType][_nprocs];
      int nr_loc = _msh->_ownSize[3][_iproc];
      int nc_loc = _msh->_ownSize[SolType][_iproc];

      for(int i = 0; i < dim; i++) {
        _GradMat[SolType][i] = SparseMatrix::build().release();
        _GradMat[SolType][i]->init(nr, nc, nr_loc, nc_loc, 27, 27);
      }

      vector< vector < double> > coordinates(dim);
      vector< int > column_dofs;
      vector< int > row_dof;
      vector <double> phi;
      vector <double> gradphi;
      vector <double> nablaphi;
      double weight;
      vector< vector< double> > B(dim);


      const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

      for(int i = 0; i < dim; i++)
        coordinates[i].reserve(max_size);

      row_dof.reserve(1);
      column_dofs.reserve(max_size);

      phi.reserve(max_size);
      gradphi.reserve(max_size * dim);
      nablaphi.reserve(max_size * (3 * (dim - 1) + !(dim - 1)));

      for(int i = 0; i < dim; i++) {
        B[i].reserve(max_size);
      }

      // Set to zeto all the entries of the Global Matrix
      for(int i = 0; i < dim; i++) {
        _GradMat[SolType][i]->zero();
      }

      unsigned nel = _msh->GetNumberOfElements();

      for(int iel = _msh->_elementOffset[_iproc]; iel < _msh->_elementOffset[_iproc + 1]; iel++) {

        row_dof.resize(1);
        row_dof[0] = iel;

        short unsigned ielt = _msh->GetElementType(iel);

        unsigned nve = _msh->GetElementDofNumber(iel, SolType);

        // resize
        column_dofs.resize(nve);
        phi.resize(nve);
        gradphi.resize(nve * dim);
        nablaphi.resize(nve * (3 * (dim - 1) + !(dim - 1)));

        for(int i = 0; i < dim; i++) {
          coordinates[i].resize(nve);
        }

        // set to zero all the entries of the FE matrices
        for(int i = 0; i < dim; i++) {
          B[i].resize(nve);
          memset(&B[i][0], 0, nve * sizeof(double));
        }

        for(unsigned i = 0; i < nve; i++) {
          unsigned inode_coord_metis = _msh->GetSolutionDof(i, iel, 2);
          column_dofs[i] = _msh->GetSolutionDof(i, iel, SolType);

          for(unsigned ivar = 0; ivar < dim; ivar++) {
            coordinates[ivar][i] = (*_msh->_topology->_Sol[ivar])(inode_coord_metis);
          }
        }

        _msh->_finiteElement[ielt][SolType]->Jacobian(coordinates, 0, weight, phi, gradphi, nablaphi);


        for(int i = 0; i < nve; i++) {
          for(int j = 0; j < dim; j++) {
            B[j][i] = gradphi[i * dim + j];
          }
        }

        for(int i = 0; i < dim; i++) {
          _GradMat[SolType][i]->add_matrix_blocked(B[i], row_dof, column_dofs);
        }
      }


      // End build elem type structure
      for(int i = 0; i < dim; i++) {
        _GradMat[SolType][i]->close();
      }
    }
  }


// ------------------------------------------------------------------
  void Solution::CopySolutionToOldSolution() {
    for(unsigned i = 0; i < _Sol.size(); i++) {
      // Copy the old vector
      if(_SolTmOrder[i] == 2) {
        *(_SolOld[i]) = *(_Sol[i]);
      }
    }
  }
  
  void Solution::ResetSolutionToOldSolution() {
    for(unsigned i = 0; i < _Sol.size(); i++) {
      // Copy the old vector
      if(_SolTmOrder[i] == 2) {
        *(_Sol[i]) = *(_SolOld[i]);
      }
    }
  }


} //end namespace femus




