/*=========================================================================

 Program: FEMUS
 Module: LinearEquation
 Authors: Eugenio Aulisa, Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

#include "LinearEquation.hpp"
#include "Mesh.hpp"
#include "ParalleltypeEnum.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"


#include <cstring>
#include <iostream>
#include <fstream>
#include <ctime>
#include <fstream>
#include <algorithm>



namespace femus {




//--------------------------------------------------------------------------------

  LinearEquation::LinearEquation(Solution *other_solution) :
     _msh ( other_solution->GetMesh() )  
     {
       
    _solution = other_solution;
    
    _EPS = NULL;
    _EPSC = NULL;
    _RES = NULL;
    _RESC = NULL;
    _KK = NULL;
    _KKamr = NULL;
    _numberOfGlobalVariables = 0u;
  }

//--------------------------------------------------------------------------------
  LinearEquation::~LinearEquation() { }


//--------------------------------------------------------------------------------
  unsigned LinearEquation::GetIndex(const char name[]) {
    unsigned index = 0;
    while(strcmp(_SolName[index], name)) {
      index++;
      if(index == _SolType.size()) {
        std::cout << "error! invalid name entry GetIndex(...)" << std::endl;
        exit(0);
      }
    }
    return index;
  }


  unsigned LinearEquation::GetSystemDof(const unsigned &index_sol, const unsigned &kkindex_sol,
                                        const unsigned &i, const unsigned &iel) const {

    unsigned soltype =  _SolType[index_sol];
    unsigned idof = _msh->GetSolutionDof(i, iel, soltype);

    unsigned isubdom = _msh->BisectionSearch_find_processor_of_dof(idof, soltype);
    return KKoffset[kkindex_sol][isubdom] + idof - _msh->_dofOffset[soltype][isubdom];
    
  }

  unsigned LinearEquation::GetSystemDof(const unsigned &soltype, const unsigned &kkindex_sol,
                                        const unsigned &i, const unsigned &iel, const std::vector < std::vector <unsigned> > &otherKKoffset) const {

    //unsigned soltype =  _SolType[index_sol];
    unsigned idof = _msh->GetSolutionDof(i, iel, soltype);

    unsigned isubdom = _msh->BisectionSearch_find_processor_of_dof(idof, soltype);
    return otherKKoffset[kkindex_sol][isubdom] + idof - _msh->_dofOffset[soltype][isubdom];
    
  }


  unsigned LinearEquation::GetSystemDof(const unsigned &index_sol, const unsigned &kkindex_sol,
                                        const unsigned &ielc, const unsigned &i0, const unsigned &i1,
                                        const Mesh* mshc) const {
                                            
    unsigned soltype =  _SolType[index_sol];
    unsigned idof = _msh->GetSolutionDof(ielc, i0, i1, soltype, mshc);

    unsigned isubdom = _msh->BisectionSearch_find_processor_of_dof(idof, soltype);
    return KKoffset[kkindex_sol][isubdom] + idof - _msh->_dofOffset[soltype][isubdom];
    
  }
  

  /** Print numeric vector with structure */
  void LinearEquation::print_residual_with_structure_matlab_friendly(const unsigned iproc, const std::string filename, NumericVector * num_vec_in) const {
      
//       if (iproc == 0) {
      
      std::ofstream file_pr;
      file_pr.open (filename);
      
      const unsigned global_size_for_all_vars = KKIndex[KKIndex.size() - 1];
      
      file_pr << "% " << "global size for all vars " << global_size_for_all_vars << std::endl;
      
    for(int ip = 0; ip < _nprocs; ip++) {
                 const unsigned local_size_for_all_vars = KKoffset[KKIndex.size() - 1][ip] - KKoffset[0][ip];
                  file_pr << "% " << "in proc " << ip << " local size for all vars "  << local_size_for_all_vars << std::endl;
         }
         
      std::vector< double > v_global(global_size_for_all_vars);   
         
      num_vec_in->localize(v_global);
//       Since I want to print the NumericVectors, I need to read from them.
// The problem is that they don't have the structure of each var
// Since these vectors are filled by looping over elements for each proc, that's what we should do
// The problem is that we would repeat the visits of nodes, so we don't want to do that
      
      unsigned running_index = 1;  //to be consistent with potential future Matlab reading
     for(int ip = 0; ip < _nprocs; ip++) {
          file_pr << "% " << "--------- processor " << ip << " ---------" << std::endl;
      for(int ivar = 0; ivar < _SolPdeIndex.size(); ivar++) {
                  file_pr  << "% " << "++++++ " << ivar << " " << _SolName[ivar] << ", ";
                  file_pr << " FE type " << _SolType[ivar] << ", ";

                      const unsigned local_size_for_var = KKoffset[ivar + 1][ip] - KKoffset[ivar][ip];
                  file_pr << " local size "  << local_size_for_var << std::endl;
                  
      for(int loc_ind = 0; loc_ind < local_size_for_var; loc_ind++) {
                  file_pr << running_index << " " << v_global[running_index-1] << std::endl;
                 running_index ++;
            }
         }
      }
      
      file_pr << std::endl;
      

      file_pr.close();          
      
//       }
      
  }

  

  void LinearEquation::sparsity_pattern_print_nonzeros(const std::string filename_base, const std::string on_or_off) {

         
      if (on_or_off != "on" && on_or_off != "off")
      { std::cout << "Must be either on or off diagonal" << std::endl; abort(); }     
          
   std::string filename =  filename_base + on_or_off + ".txt";
      std::ofstream file_pr;
      file_pr.open (filename);
      
         std::vector<int> * nonzeros; 
      if (on_or_off == "on") {
          nonzeros = &d_nnz;
      }
      else if (on_or_off == "off") {
          nonzeros = &o_nnz;
      }
      
      
      for(unsigned j = 0; j < nonzeros->size(); j++) {
                  file_pr << j << " " << (*nonzeros)[j] << std::endl;
                }
                
      
   file_pr.close();          
  
  
  }
  
  
//--------------------------------------------------------------------------------
  void LinearEquation::InitPde(const std::vector <unsigned> &SolPdeIndex_other, 
                               const std::vector <unsigned> &SolType_other,
                               const std::vector <char*> &SolName_other,
                               std::vector <NumericVector*> *Bdc_other,
                               const unsigned &other_gridn, 
                               std::vector <bool> &SparsityPattern_other) {
      
    _SolPdeIndex = SolPdeIndex_other;
    _gridn = other_gridn;

    _SolType = SolType_other;
    _SolName = SolName_other;
    _Bdc = Bdc_other;

    _SparsityPattern = SparsityPattern_other;

    //--- Matrix and vectors offsets - BEGIN ---------------------------------------------------------------------------------------------
    KKIndex.resize(_SolPdeIndex.size() + 1u);
    KKIndex[0] = 0;
    for(unsigned i = 1; i < KKIndex.size(); i++)
      KKIndex[i] = KKIndex[i - 1] + _msh->_dofOffset[_SolType[_SolPdeIndex[i - 1]]][n_processors()];

    KKoffset.resize(_SolPdeIndex.size() + 1);
    for(int i = 0; i < KKoffset.size(); i++) {
      KKoffset[i].resize(_nprocs);
    }

    //KKoffset for each var, at processor 0
    KKoffset[0][0] = 0;
    for(int j = 1; j < KKoffset.size(); j++) {
      unsigned indexSol = _SolPdeIndex[j - 1];
      KKoffset[j][0] = KKoffset[j - 1][0] + (_msh->_dofOffset[_SolType[indexSol]][1] - _msh->_dofOffset[_SolType[indexSol]][0]);
    }

    //KKoffset for each var, for all other processors
    for(int i = 1; i < _nprocs; i++) {
      KKoffset[0][i] = KKoffset[_SolPdeIndex.size()][i - 1];
      for(int j = 1; j < KKoffset.size(); j++) {
        unsigned indexSol = _SolPdeIndex[j - 1];
        KKoffset[j][i] = KKoffset[j - 1][i] + (_msh->_dofOffset[_SolType[indexSol]][i + 1] - _msh->_dofOffset[_SolType[indexSol]][i]);
      }
    }

    //ghost size
    KKghostsize.resize(_nprocs, 0);
    for(int i = 0; i < _nprocs; i++) {
      for(int j = 0; j < _SolPdeIndex.size(); j++) {
        unsigned indexSol = _SolPdeIndex[j];
        KKghostsize[i] += _msh->dofmap_get_ghost_dofs(_SolType[indexSol], i).size();
      }
    }

    //ghost nodes
    KKghost_nd.resize(_nprocs);
    for(int i = 1; i < _nprocs; i++) {
      KKghost_nd[i].resize(KKghostsize[i]);
    }


    for(int i = 0; i < _nprocs; i++) {
      unsigned counter = 0;
      for(int j = 0; j < _SolPdeIndex.size(); j++) {
        unsigned indexSol = _SolPdeIndex[j];
        for(int k = 0; k < _msh->dofmap_get_ghost_dofs(_SolType[indexSol], i).size(); k++) {
          // ghost node
          unsigned idof_metis = _msh->dofmap_get_ghost_dofs(_SolType[indexSol], i)[k];
          unsigned isubdom = _msh->BisectionSearch_find_processor_of_dof(idof_metis, _SolType[indexSol]);
          KKghost_nd[i][counter] = KKoffset[j][isubdom] + idof_metis - _msh->_dofOffset[_SolType[indexSol]][isubdom];
          counter++;
        }
      }
    }
    //--- Matrix and vectors offsets - END ---------------------------------------------------------------------------------------------
    

    //--- Error and residual: build and init - BEGIN --------------------------------------------------------------------------------------------
    int EPSsize = KKIndex[KKIndex.size() - 1];
    _EPS = NumericVector::build().release();
    if(n_processors() == 1) {  // IF SERIAL
      _EPS->init(EPSsize, EPSsize, false, SERIAL);
    }
    else { // IF PARALLEL
      int EPS_local_size = KKoffset[KKIndex.size() - 1][processor_id()] - KKoffset[0][processor_id()];
      _EPS->init(EPSsize, EPS_local_size, KKghost_nd[processor_id()], false, GHOSTED);
    }

    _RES = NumericVector::build().release();
    _RES->init(*_EPS);

    _EPSC = NumericVector::build().release();
    _EPSC->init(*_EPS);

    _RESC = NumericVector::build().release();
    _RESC->init(*_EPS);
    //--- Error and residual: build and init - END --------------------------------------------------------------------------------------------


    //--- Matrix: Sparsity pattern: fill it as if all variables were sparse - BEGIN  -----------------------------------------------------------------------------------------------
    GetSparsityPatternSize();
    //--- Matrix: Sparsity pattern: fill it as if all variables were sparse - END  -----------------------------------------------------------------------------------------------

    
    //--- Matrix: Sparsity pattern: adjust it for dense variables - BEGIN -----------------------------------------------------------------------------------------------
//     const unsigned dim = _msh->GetDimension();
//     int KK_UNIT_SIZE_ = pow(5, dim);
    int KK_size = KKIndex[KKIndex.size() - 1u];
    int KK_local_size = KKoffset[KKIndex.size() - 1][processor_id()] - KKoffset[0][processor_id()]; //sum of all vars

    if(_sparsityPatternMinimumSize.size() > 0) {
        
      for(unsigned i = 0; i < _sparsityPatternVariableIndex.size(); i++) {
                 
        unsigned maxDiagSize    = (_sparsityPatternMinimumSize[i] < KK_local_size) ? _sparsityPatternMinimumSize[i] : KK_local_size; //minimum between SPMinSize and the sum of all vars: I cannot clearly have a number of columns greater than the sum of all vars
        unsigned maxOffDiagSize = (_sparsityPatternMinimumSize[i] < KK_size - KK_local_size) ? _sparsityPatternMinimumSize[i] : KK_size - KK_local_size;  
          
        unsigned idx = _sparsityPatternVariableIndex[i];   
        unsigned jstart = KKoffset[idx][processor_id()] - KKoffset[0][processor_id()];
        unsigned jend = KKoffset[idx + 1][processor_id()] - KKoffset[0][processor_id()];
        
        for(unsigned j = jstart; j < jend; j++) {
          d_nnz[j] = (d_nnz[j] > maxDiagSize) ? d_nnz[j] : maxDiagSize;
          o_nnz[j] = (o_nnz[j] > maxOffDiagSize) ? o_nnz[j] : maxOffDiagSize;
        }
      }
      
    }

    if(_iproc == 0) {
      unsigned globalVariableStart = KKoffset[KKIndex.size() - 2][0];
      for(unsigned k = globalVariableStart; k < globalVariableStart + _numberOfGlobalVariables; k++) {
        d_nnz[k] = KK_local_size;
        o_nnz[k] = KK_size - KK_local_size;
      }
    }
    //--- Matrix: Sparsity pattern: adjust it for dense variables - END -----------------------------------------------------------------------------------------------
    

    //--- Matrix: build and init - BEGIN --------------------------------------------------------------------------------------------
    _KK = SparseMatrix::build().release();
    _KK->init(KK_size, KK_size, KK_local_size, KK_local_size, d_nnz, o_nnz);
    //--- Matrix: build and init - END  --------------------------------------------------------------------------------------------

    //--- Matrix AMR: build - BEGIN --------------------------------------------------------------------------------------------
    _KKamr = SparseMatrix::build().release();
    //--- Matrix AMR: build - END --------------------------------------------------------------------------------------------

  }

  void LinearEquation::SetSparsityPatternMinimumSize(const std::vector < unsigned> &minimumSize, const std::vector < unsigned > &variableIndex) {
    _sparsityPatternMinimumSize = minimumSize;
    _sparsityPatternVariableIndex = variableIndex;

  }

//--------------------------------------------------------------------------------
  void LinearEquation::AddLevel() {
    _gridn++;
  }

//--------------------------------------------------------------------------------
  void LinearEquation::SetResZero() {
    _RES->zero();
  }

//--------------------------------------------------------------------------------
  void LinearEquation::SetEpsZero() {
    _EPS->zero();
    _EPSC->zero();
  }

//--------------------------------------------------------------------------------
  void LinearEquation::SumEpsCToEps() {
    *_EPS += *_EPSC;
  }

//--------------------------------------------------------------------------------
  void LinearEquation::UpdateResidual() {
    _RESC->matrix_mult(*_EPSC, *_KK);
    *_RES -= *_RESC;
  }

//-------------------------------------------------------------------------------------------
  void LinearEquation::DeletePde() {

    if(_KK)
      delete _KK;

    if(_KKamr)
      delete _KKamr;

    if(_EPS)
      delete _EPS;

    if(_EPSC)
      delete _EPSC;

    if(_RES)
      delete _RES;

    if(_RESC)
      delete _RESC;

    _EPS = NULL;
    _EPSC = NULL;
    _RES = NULL;
    _RESC = NULL;
    _KK = NULL;
    _KKamr = NULL;

  }

  void LinearEquation::GetSparsityPatternSize() {

    unsigned SolPdeSize = _SolPdeIndex.size();
    if(_SparsityPattern.size() == 0) {
      _SparsityPattern.resize(SolPdeSize * SolPdeSize);
      for(int i = 0; i < SolPdeSize; i++) {
        for(int j = 0; j < SolPdeSize; j++) {
          _SparsityPattern[SolPdeSize * i + j] = 1;
        }
      }
    }
    else if(_SparsityPattern.size() != SolPdeSize * SolPdeSize) {
      std::cout << "Sparsity Pattern size ( " << _SparsityPattern.size() << " ) does not match system PDE size" << std::endl;
      exit(0);
    }

    const int dim = _msh->GetDimension();

    const int max_size = static_cast< int >(ceil(pow(3, dim)));

    std::vector < std::vector < int > > dofsVAR(_SolPdeIndex.size());

    for(int i = 0; i < _SolPdeIndex.size(); i++) {
      dofsVAR[i].reserve(max_size);
    }

    // mesh and procs
    int nel    = _msh->GetNumberOfElements();
    int igrid  = _msh->GetLevel();
    int this_proc  = _msh->processor_id();
    int nprocs =       _msh->n_processors();

    // *** element loop ***

    int IndexStart = KKoffset[0][this_proc];
    int IndexEnd  = KKoffset[KKIndex.size() - 1][this_proc];
    int owned_dofs    = IndexEnd - IndexStart;  //all dofs of all vars of this_proc

    std::vector < std::map < int, bool > > BlgToMe_d(owned_dofs); //belong to me
    std::vector < std::map < int, bool > > BlgToMe_o(owned_dofs); //belong to me
    std::map < int, std::map <int, bool > > DnBlgToMe_o;     //don't belong to me 
    std::map < int, std::map <int, bool > > DnBlgToMe_d;     //don't belong to me

    for(int kel = _msh->GetElementOffset( this_proc ); kel < _msh->GetElementOffset(this_proc + 1); kel++) {

      short int kelt = _msh->GetElementType(kel);
      std::vector < int > nve(_SolPdeIndex.size());
      
      for(int i = 0; i < _SolPdeIndex.size(); i++) {
        nve[i] = _msh->GetElementDofNumber(kel, _SolType[_SolPdeIndex[i]]);
      }
      for(int i = 0; i < _SolPdeIndex.size(); i++) {
        dofsVAR[i].resize(nve[i]);
      }

      for(int i = 0; i < _SolPdeIndex.size(); i++) {
        int ThisSolType = _SolType[_SolPdeIndex[i]];
        for(int j = 0; j < nve[i]; j++) {
          int inode = _msh->GetSolutionDof(j, kel, ThisSolType);
          dofsVAR[i][j] = GetSystemDof(_SolPdeIndex[i], i, j, kel);
        }
      }
      
      for(int i = 0; i < _SolPdeIndex.size(); i++) {
          
        for(int inode = 0; inode < nve[i]; inode++) {
          int idof_local = dofsVAR[i][inode] - IndexStart;
          
          for(int j = 0; j < _SolPdeIndex.size(); j++) {
              
            if(_SparsityPattern[_SolPdeIndex.size() *i + j]) {
                
              for(int jnode = 0; jnode < nve[j]; jnode++) {
                int jdof_local = dofsVAR[j][jnode] - IndexStart;

                if(idof_local >= 0) {  // i-row belongs to this proc
                  if(jdof_local >= 0) {  // j-row belongs to this proc (diagonal)
                    BlgToMe_d[ idof_local ][ jdof_local ] = 1;
                  }
                  else { // j-row does not belong to this proc (off-diagonal)
                    BlgToMe_o[ idof_local ][ dofsVAR[j][jnode] ] = 1;
                  }
                }
                else { // i-row does not belong to this proc
                  // identify the process the i-row belongs to
                  int iproc = 0;
                  while(dofsVAR[i][inode] >= KKoffset[KKIndex.size() - 1][iproc]) iproc++;

                  // identify the process the j-column belongs to
                  int jproc = 0;
                  while(dofsVAR[j][jnode] >= KKoffset[KKIndex.size() - 1][jproc]) jproc++;
                  if(iproc != jproc) {  // if diagonal
                    DnBlgToMe_o[ dofsVAR[i][inode] ][ dofsVAR[j][jnode] ] = 1;
                  }
                  else if(iproc == jproc) {  // if off-diagonal
                    DnBlgToMe_d[ dofsVAR[i][inode] ][ dofsVAR[j][jnode] ] = 1;
                  }
                }
                
              }  //ndofs jvar
              
            }  //test if there is coupling between ivar and jvar
          } //jvar
        }  //ndofs ivar
      }  //ivar
      
    } //end el loop

    NumericVector  *sizeDnBM_o = NumericVector::build().release();
    sizeDnBM_o->init(*_EPS);
    sizeDnBM_o->zero();
    for(std::map < int, std::map <int, bool > >::iterator it = DnBlgToMe_o.begin(); it != DnBlgToMe_o.end(); ++it) {
      sizeDnBM_o->add(it->first, it->second.size());
    }
    sizeDnBM_o->close();


    NumericVector  *sizeDnBM_d = NumericVector::build().release();
    sizeDnBM_d->init(*_EPS);
    sizeDnBM_d->zero();
    for(std::map < int, std::map <int, bool > >::iterator it = DnBlgToMe_d.begin(); it != DnBlgToMe_d.end(); ++it) {
      sizeDnBM_d->add(it->first, it->second.size());
    }
    sizeDnBM_d->close();


    d_nnz.resize(owned_dofs);
    o_nnz.resize(owned_dofs);

    int d_max = owned_dofs + _numberOfGlobalVariables;
    int o_max = KKIndex[KKIndex.size() - 1u] - owned_dofs + _numberOfGlobalVariables;

    for(int i = 0; i < owned_dofs; i++) {
      d_nnz[i] = static_cast <int>((*sizeDnBM_d)(IndexStart + i)) + BlgToMe_d[i].size() + _numberOfGlobalVariables;
      if(d_nnz[i] > d_max) d_nnz[i] = d_max;
      o_nnz[i] = static_cast <int>((*sizeDnBM_o)(IndexStart + i)) + BlgToMe_o[i].size() + _numberOfGlobalVariables;
      if(o_nnz[i] > o_max) o_nnz[i] = o_max;
    }

    delete sizeDnBM_o;
    delete sizeDnBM_d;
  }
}






