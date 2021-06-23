/*=========================================================================

 Program: FEMUS
 Module: NonLinearImplicitSystem
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"
#include "ElemType.hpp"
#include "MultiLevelSolution.hpp"


namespace femus {





// ------------------------------------------------------------
// NonLinearImplicitSystem implementation
  MonolithicFSINonLinearImplicitSystem::MonolithicFSINonLinearImplicitSystem(MultiLevelProblem& ml_probl,
      const std::string& name_in,
      const unsigned int number_in, const LinearEquationSolverType & smoother_type) :
    NonLinearImplicitSystem(ml_probl, name_in, number_in, smoother_type) {
  }

  MonolithicFSINonLinearImplicitSystem::~MonolithicFSINonLinearImplicitSystem() {
  }

  void MonolithicFSINonLinearImplicitSystem::init() {
    Parent::init();
  }

//---------------------------------------------------------------------------------------------------
// This routine generates the matrix for the projection of the FE matrix to finer grids
//---------------------------------------------------------------------------------------------------

  void MonolithicFSINonLinearImplicitSystem::BuildProlongatorMatrix(unsigned gridf) {

    if(gridf < 1) {
      std::cout << "Error! In function \"BuildProlongatorMatrix\" argument less then 1" << std::endl;
      exit(0);
    }

    int iproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

    LinearEquationSolver* LinSolf = _LinSolver[gridf];
    LinearEquationSolver* LinSolc = _LinSolver[gridf - 1];
    Mesh* mshc = _msh[gridf - 1];
    int nf = LinSolf->KKIndex[LinSolf->KKIndex.size() - 1u];
    int nc = LinSolc->KKIndex[LinSolc->KKIndex.size() - 1u];
    int nf_loc = LinSolf->KKoffset[LinSolf->KKIndex.size() - 1][iproc] - LinSolf->KKoffset[0][iproc];
    int nc_loc = LinSolc->KKoffset[LinSolc->KKIndex.size() - 1][iproc] - LinSolc->KKoffset[0][iproc];

    NumericVector *NNZ_d = NumericVector::build().release();
    NNZ_d->init(*LinSolf->_EPS);
    NNZ_d->zero();

    NumericVector *NNZ_o = NumericVector::build().release();
    NNZ_o->init(*LinSolf->_EPS);
    NNZ_o->zero();

    for(unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
      unsigned SolIndex = _SolSystemPdeIndex[k];
      unsigned  SolType = _ml_sol->GetSolutionType(SolIndex);
      // loop on the coarse grid
      for(int isdom = iproc; isdom < iproc + 1; isdom++) {
        for(int iel = mshc->_elementOffset[isdom]; iel < mshc->_elementOffset[isdom + 1]; iel++) {
          short unsigned ielt = mshc->GetElementType(iel);
          mshc->_finiteElement[ielt][SolType]->GetSparsityPatternSize(*LinSolf, *LinSolc, iel, NNZ_d, NNZ_o, SolIndex, k);
        }
      }
    }

    NNZ_d->close();
    NNZ_o->close();

    unsigned offset = LinSolf->KKoffset[0][iproc];
    vector <int> nnz_d(nf_loc);
    vector <int> nnz_o(nf_loc);
    for(int i = 0; i < nf_loc; i++) {
      nnz_d[i] = static_cast <int>((*NNZ_d)(offset + i));
      nnz_o[i] = static_cast <int>((*NNZ_o)(offset + i));
    }

    delete NNZ_d;
    delete NNZ_o;

    _PP[gridf] = SparseMatrix::build().release();
    _PP[gridf]->init(nf, nc, nf_loc, nc_loc, nnz_d, nnz_o);

    SparseMatrix *RRt;
    RRt = SparseMatrix::build().release();
    RRt->init(nf, nc, nf_loc, nc_loc, nnz_d, nnz_o);

    for(unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
      unsigned SolIndex = _SolSystemPdeIndex[k];
      unsigned solPairIndex = _ml_sol->GetSolutionPairIndex(k);
      unsigned SolType = _ml_sol->GetSolutionType(SolIndex);
      unsigned solPairPdeIndex = GetSolPdeIndex(_ml_sol->GetSolutionName(solPairIndex));

      bool testIfPressure = 0;
      if(_ml_sol->TestIfSolutionIsPressure(SolIndex))   testIfPressure = 1;

      // loop on the coarse grid
      for(int isdom = iproc; isdom < iproc + 1; isdom++) {
        for(int iel = mshc->_elementOffset[isdom]; iel < mshc->_elementOffset[isdom + 1]; iel++) {
          short unsigned ielt = mshc->GetElementType(iel);
          if(!testIfPressure) {
            mshc->_finiteElement[ielt][SolType]->BuildRestrictionTranspose(*LinSolf, *LinSolc, iel, RRt, SolIndex, k, solPairIndex, solPairPdeIndex);
          }
          else {
            mshc->_finiteElement[ielt][SolType]->BuildProlongation(*LinSolf, *LinSolc, iel, RRt, SolIndex, k);
          }
          mshc->_finiteElement[ielt][SolType]->BuildProlongation(*LinSolf, *LinSolc, iel, _PP[gridf], SolIndex, k);
        }
      }
    }

    _PP[gridf]->close();
    RRt->close();

    _RR[gridf] = SparseMatrix::build().release();
    RRt->get_transpose(*_RR[gridf]);
    delete RRt;

  }
  
  void MonolithicFSINonLinearImplicitSystem::BuildAmrProlongatorMatrix(unsigned level) {

    int iproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

    LinearEquationSolver* LinSol = _LinSolver[level];

    Mesh* mesh = _msh[level];
    std::vector < std::map < unsigned,  std::map < unsigned, double  > > > &amrRestriction = mesh->GetAmrRestrictionMap();

    std::vector < std::map < unsigned, bool > > & amrSolidMark = mesh->GetAmrSolidMark();
    
    int n = LinSol->KKIndex[LinSol->KKIndex.size() - 1u];
    int n_loc = LinSol->KKoffset[LinSol->KKIndex.size() - 1][iproc] - LinSol->KKoffset[0][iproc];

    NumericVector* NNZ_d = NumericVector::build().release();
    NNZ_d->init(*LinSol->_EPS);
    NNZ_d->zero();

    NumericVector* NNZ_o = NumericVector::build().release();
    NNZ_o->init(*LinSol->_EPS);
    NNZ_o->zero();

    for(unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {

      unsigned solIndex = _SolSystemPdeIndex[k];
      unsigned solType = _ml_sol->GetSolutionType(solIndex);
      unsigned kOffset = LinSol->KKoffset[k][iproc];

      unsigned solOffset = mesh->_dofOffset[solType][iproc];
      unsigned solOffsetp1 = mesh->_dofOffset[solType][iproc + 1];
      for(unsigned i = solOffset; i < solOffsetp1; i++) {
        if(solType > 2 || amrRestriction[solType].find(i) == amrRestriction[solType].end()) {
          NNZ_d->set(kOffset + (i - solOffset), 1);
        }
        else {
          double cnt_d = 0;
          double cnt_o = 0;
          for(std::map<unsigned, double> ::iterator it = amrRestriction[solType][i].begin(); it != amrRestriction[solType][i].end(); it++) {
            if(it->first >= solOffset && it->first < solOffsetp1) {
              cnt_d++;
            }
            else {
              cnt_o++;
            }
          }
          NNZ_d->set(kOffset + (i - solOffset), cnt_d);
          NNZ_o->set(kOffset + (i - solOffset), cnt_o);
        }
      }
    }

    NNZ_d->close();
    NNZ_o->close();

    unsigned offset = LinSol->KKoffset[0][iproc];
    vector <int> nnz_d(n_loc);
    vector <int> nnz_o(n_loc);

    for(int i = 0; i < n_loc; i++) {
      nnz_d[i] = static_cast <int>((*NNZ_d)(offset + i));
      nnz_o[i] = static_cast <int>((*NNZ_o)(offset + i));
    }

    delete NNZ_d;
    delete NNZ_o;

    _PPamr[level] = SparseMatrix::build().release();
    _PPamr[level]->init(n, n, n_loc, n_loc, nnz_d, nnz_o);
    
    _RRamr[level] = SparseMatrix::build().release();
    _RRamr[level]->init(n, n, n_loc, n_loc, nnz_d, nnz_o);

    for(unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
     
      unsigned solIndex = _SolSystemPdeIndex[k];
      unsigned solPairIndex = _ml_sol->GetSolutionPairInverseIndex(k);
      
      unsigned kPair = GetSolPdeIndex(_ml_sol->GetSolutionName(solPairIndex));
            
      unsigned kOffset = LinSol->KKoffset[k][iproc];
      unsigned kPairOffset = LinSol->KKoffset[kPair][iproc];
      
      unsigned solType = _ml_sol->GetSolutionType(solIndex);

      unsigned solOffset = mesh->_dofOffset[solType][iproc];
      unsigned solOffsetp1 = mesh->_dofOffset[solType][iproc + 1];

      for(unsigned i = solOffset; i < solOffsetp1; i++) {
	unsigned irow = kOffset + (i - solOffset);
        if(solType > 2 || amrRestriction[solType].find(i) == amrRestriction[solType].end()) {
          std::vector <int> col(1, irow);
          double value = 1.;
          _PPamr[level]->insert_row(irow, 1, col, &value);
	  _RRamr[level]->insert_row(irow, 1, col, &value);
        }
        else {
	  bool solidMarki = amrSolidMark[solType][i];
	  int ncols = amrRestriction[solType][i].size();
          std::vector <int> colPP(ncols);
	  std::vector <int> colRR(ncols);
          std::vector <double> valuePP(ncols);
	  std::vector <double> valueRR(ncols);
          unsigned j = 0;
          for(std::map<unsigned, double> ::iterator it = amrRestriction[solType][i].begin(); it != amrRestriction[solType][i].end(); it++) {
	    bool solidMarkj = amrSolidMark[solType][it->first];
            if(it->first >= solOffset && it->first < solOffsetp1) {
              colPP[j] = kOffset + (it->first - solOffset);
	      if(solidMarki == true && solidMarkj == false && k != kPair){
		colRR[j] = kPairOffset + (it->first - solOffset);
	      }
	      else{
		colRR[j] = kOffset + (it->first - solOffset);
	      }
            }
            else {
              unsigned jproc = _msh[level]->IsdomBisectionSearch(it->first, solType);
              colPP[j] = LinSol->KKoffset[k][jproc] + (it->first - mesh->_dofOffset[solType][jproc]);
	      if(solidMarki == true && solidMarkj == false && k != kPair){
		colRR[j] = LinSol->KKoffset[kPair][jproc] + (it->first - mesh->_dofOffset[solType][jproc]);
	      }
	      else{
		colRR[j] = LinSol->KKoffset[k][jproc] + (it->first - mesh->_dofOffset[solType][jproc]);
	      }
            }
            valuePP[j] = it->second;
	    if( solidMarki == true && solidMarkj == false  && k == kPair){
	      valueRR[j] = 0.;   
	    }
	    else{
	      valueRR[j] = it->second;
	    }
	    j++;
          }
          _PPamr[level]->insert_row(irow, ncols, colPP, &valuePP[0]);
	  _RRamr[level]->insert_row(irow, ncols, colRR, &valueRR[0]);
        }
      }
    }
    _PPamr[level]->close();
    _RRamr[level]->close();
    
    _PPamr[level]->get_transpose(*_PPamr[level]);
  }
  
  

  void MonolithicFSINonLinearImplicitSystem::SetElementBlockFluidAll() {
    for(unsigned i = 1; i < _gridn; i++) {
      _LinSolver[i]->SetElementBlockNumberFluid(_msh[i]->GetNumberOfElements(), 0);
    }
  }

  void MonolithicFSINonLinearImplicitSystem::SetElementBlockSolidAll() {
    for(unsigned i = 1; i < _gridn; i++) {
      _LinSolver[i]->SetElementBlockNumberSolid(_msh[i]->GetNumberOfElements(), 0);
    }
  }
  
  void MonolithicFSINonLinearImplicitSystem::SetElementBlockPorousAll() {
    for(unsigned i = 1; i < _gridn; i++) {
      _LinSolver[i]->SetElementBlockNumberPorous(_msh[i]->GetNumberOfElements(), 0);
    }
  }


  void MonolithicFSINonLinearImplicitSystem::SetElementBlockNumberFluid(unsigned const &dim_block, unsigned const &overlap) {
    _numblock_test = 1;
    const unsigned dim = _msh[0]->GetDimension();
    const unsigned base = pow(2, dim);
    _num_block = pow(base, dim_block);

    for(unsigned i = 1; i < _gridn; i++) {
      unsigned num_block2 = std::min(_num_block, _msh[i]->GetNumberOfElements());
      _LinSolver[i]->SetElementBlockNumberFluid(num_block2, overlap);
    }
  }

  void MonolithicFSINonLinearImplicitSystem::SetElementBlockNumberSolid(unsigned const &dim_block, unsigned const &overlap) {
    _numblock_test = 1;
    const unsigned dim = _msh[0]->GetDimension();
    const unsigned base = pow(2, dim);
    _num_block = pow(base, dim_block);

    for(unsigned i = 1; i < _gridn; i++) {
      unsigned num_block2 = std::min(_num_block, _msh[i]->GetNumberOfElements());
      _LinSolver[i]->SetElementBlockNumberSolid(num_block2, overlap);
    }
  }
  
  void MonolithicFSINonLinearImplicitSystem::SetElementBlockNumberPorous(unsigned const &dim_block, unsigned const &overlap) {
    _numblock_test = 1;
    const unsigned dim = _msh[0]->GetDimension();
    const unsigned base = pow(2, dim);
    _num_block = pow(base, dim_block);

    for(unsigned i = 1; i < _gridn; i++) {
      unsigned num_block2 = std::min(_num_block, _msh[i]->GetNumberOfElements());
      _LinSolver[i]->SetElementBlockNumberPorous(num_block2, overlap);
    }
  }


} //end namespace femus

