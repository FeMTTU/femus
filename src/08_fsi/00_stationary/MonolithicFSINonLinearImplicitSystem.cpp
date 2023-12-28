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

  
  

  void MonolithicFSINonLinearImplicitSystem::Build_RestrictionTranspose_OneElement_OneFEFamily_With_Pair_In_System(const LinearEquation& lspdef,
                                            const LinearEquation& lspdec, 
                                            const int& ielc, 
                                            SparseMatrix* Projmat,
                                            const unsigned& index_sol, 
                                            const unsigned& kkindex_sol,
                                            const unsigned& index_pair_sol,
                                            const unsigned& kkindex_pair_sol,
                                            const elem_type * elem_type_in) const
  {
    
    const unsigned soltype_in = elem_type_in->GetSolType();
    const unsigned      ndofs = elem_type_in->GetNDofs();
    const unsigned ndofs_fine = elem_type_in->GetNDofsFine();
    

    if(lspdec.GetMeshFromLinEq()->GetRefinedElementIndex(ielc)) {  // coarse2fine prolongation

      //BEGIN project nodeSolidMark
      std::vector < double > fineNodeSolidMark(ndofs_fine, 0);
      std::vector < bool > coarseNodeSolidMark(ndofs, 0);

      if( soltype_in == 2 ) {
        for(unsigned j = 0; j < ndofs; j++) {
          int jadd = lspdec.GetMeshFromLinEq()->GetSolutionDof(j, ielc, soltype_in);
          coarseNodeSolidMark[j] = lspdec.GetMeshFromLinEq()->GetSolidMark(jadd);
        }

        for(unsigned i = 0; i < ndofs_fine; i++) {
          int ncols = elem_type_in->Get_Prolongator_Num_Columns(i);

          for(int k = 0; k < ncols; k++) {
            int j = elem_type_in->Get_Prolongator_Index(i, k);
            fineNodeSolidMark[i] += elem_type_in->Get_Prolongator_Value(i, k) * coarseNodeSolidMark[j];
          }
        }
      }

      //END project nodeSolidMark

      std::vector <int> cols( ndofs );
      std::vector <double> copy_prol_val;
      copy_prol_val.reserve( ndofs );

      for(int i = 0; i < ndofs_fine; i++) {
        
        const std::pair<int, int> id_0_1 = elem_type_in->GetKVERT_IND(i);

        int irow = lspdef.GetSystemDof(index_sol, kkindex_sol, ielc, id_0_1.first, id_0_1.second, lspdec.GetMeshFromLinEq());
        
        const bool isolidmark = (fineNodeSolidMark[i] > 0.99 && fineNodeSolidMark[i] < 1.01) ? true : false;

        int ncols = elem_type_in->Get_Prolongator_Num_Columns(i);
        cols.assign(ncols, 0);
        copy_prol_val.assign(ncols, 0);

        for(int k = 0; k < ncols; k++) {
          int j = elem_type_in->Get_Prolongator_Index(i, k);
          bool jsolidmark = coarseNodeSolidMark[j];

          if(isolidmark == jsolidmark) {
            int jcolumn = lspdec.GetSystemDof(index_sol, kkindex_sol, j, ielc);
            cols[k] = jcolumn;
            copy_prol_val[k] = elem_type_in->Get_Prolongator_Value(i, k);
          }
          else {
            int jcolumn = lspdec.GetSystemDof(index_pair_sol, kkindex_pair_sol, j, ielc);
            cols[k] = jcolumn;
            copy_prol_val[k] = (index_sol != index_pair_sol) ? elem_type_in->Get_Prolongator_Value(i, k) : 0.;
          }
        }

        Projmat->insert_row(irow, ncols, cols, &copy_prol_val[0]);
      }
      
    }
    else {
      
      std::vector <int> jcol(1);
      double one = 1.;

      for(int i = 0; i < ndofs; i++) {
        
        int irow = lspdef.GetSystemDof(index_sol, kkindex_sol, ielc, 0, i, lspdec.GetMeshFromLinEq());
        jcol[0] = lspdec.GetSystemDof(index_sol, kkindex_sol, i, ielc);
        Projmat->insert_row(irow, 1, jcol, &one);
        
      }
    }
    
    
  }

  
  
//---------------------------------------------------------------------------------------------------
// This routine generates the matrix for the projection of the FE matrix to finer grids
//---------------------------------------------------------------------------------------------------

  void MonolithicFSINonLinearImplicitSystem::BuildProlongatorMatrix(unsigned gridf) {

    // ------------------- Sparsity pattern size - BEGIN
    
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
        for(int iel = mshc->GetElementOffset(isdom); iel < mshc->GetElementOffset(isdom + 1); iel++) {
          short unsigned ielt = mshc->GetElementType(iel);
          Get_Prolongation_SparsityPatternSize_OneElement_OneFEFamily_In_System(*LinSolf, *LinSolc, iel, NNZ_d, NNZ_o, SolIndex, k,
                                                                                                                     mshc->GetFiniteElement(ielt, SolType) );
        }
      }
    }

    NNZ_d->close();
    NNZ_o->close();

    unsigned offset = LinSolf->KKoffset[0][iproc];
    std::vector <int> nnz_d(nf_loc);
    std::vector <int> nnz_o(nf_loc);
    for(int i = 0; i < nf_loc; i++) {
      nnz_d[i] = static_cast <int>((*NNZ_d)(offset + i));
      nnz_o[i] = static_cast <int>((*NNZ_o)(offset + i));
    }

    delete NNZ_d;
    delete NNZ_o;
    // ------------------- Sparsity pattern size - END
    
    // ------------------- Prolongator - BEGIN
    

    _PP[gridf] = SparseMatrix::build().release();
    _PP[gridf]->init(nf, nc, nf_loc, nc_loc, nnz_d, nnz_o);

    SparseMatrix *RRt;
    RRt = SparseMatrix::build().release();
    RRt->init(nf, nc, nf_loc, nc_loc, nnz_d, nnz_o);

    for(unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
      unsigned SolIndex = _SolSystemPdeIndex[k];
      unsigned solPairIndex = _ml_sol->GetSolutionPairIndex(k);
      unsigned SolType = _ml_sol->GetSolutionType(SolIndex);
      unsigned solPairPdeIndex = GetSolPdeIndex(_ml_sol->GetSolName_from_index(solPairIndex));

      bool testIfPressure = 0;
      if(_ml_sol->TestIfSolutionIsPressure(SolIndex))   testIfPressure = 1;

      // loop on the coarse grid
      for(int isdom = iproc; isdom < iproc + 1; isdom++) {
        for(int iel = mshc->GetElementOffset(isdom); iel < mshc->GetElementOffset(isdom + 1); iel++) {
          short unsigned ielt = mshc->GetElementType(iel);
          
          if(!testIfPressure) {
            Build_RestrictionTranspose_OneElement_OneFEFamily_With_Pair_In_System(*LinSolf, *LinSolc, iel, RRt, SolIndex, k, solPairIndex, solPairPdeIndex,  mshc->GetFiniteElement(ielt, SolType));
          }
          else {
            Build_Prolongation_OneElement_OneFEFamily_In_System(*LinSolf, *LinSolc, iel, RRt, SolIndex, k, 
                                                                                                  mshc->GetFiniteElement(ielt, SolType) );
          }
          
            Build_Prolongation_OneElement_OneFEFamily_In_System(*LinSolf, *LinSolc, iel, _PP[gridf], SolIndex, k, 
                                                                                                  mshc->GetFiniteElement(ielt, SolType) );
          
        }
      }
    }

    _PP[gridf]->close();
    RRt->close();

    _RR[gridf] = SparseMatrix::build().release();
    RRt->get_transpose(*_RR[gridf]);
    delete RRt;
    
    // ------------------- Prolongator - END

    
  }
  
  void MonolithicFSINonLinearImplicitSystem::BuildAmrProlongatorMatrix(unsigned level) {

    // ------------------- Sparsity pattern size - BEGIN

    int iproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

    LinearEquationSolver* LinSol = _LinSolver[level];

    const Mesh* mesh = _msh[level];
    const std::vector < std::map < unsigned,  std::map < unsigned, double  > > > & amrRestriction = mesh->GetAmrRestrictionMap();

    const std::vector < std::map < unsigned, bool > > & amrSolidMark = mesh->GetAmrSolidMark();
    
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
          for(std::map<unsigned, double> ::const_iterator it = amrRestriction[solType].at(i).begin(); it != amrRestriction[solType].at(i).end(); it++) {
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
    std::vector <int> nnz_d(n_loc);
    std::vector <int> nnz_o(n_loc);

    for(int i = 0; i < n_loc; i++) {
      nnz_d[i] = static_cast <int>((*NNZ_d)(offset + i));
      nnz_o[i] = static_cast <int>((*NNZ_o)(offset + i));
    }

    delete NNZ_d;
    delete NNZ_o;
    // ------------------- Sparsity pattern size - END
    
    
    // ------------------- Prolongator - BEGIN

    _PPamr[level] = SparseMatrix::build().release();
    _PPamr[level]->init(n, n, n_loc, n_loc, nnz_d, nnz_o);
    
    _RRamr[level] = SparseMatrix::build().release();
    _RRamr[level]->init(n, n, n_loc, n_loc, nnz_d, nnz_o);

    for(unsigned k = 0; k < _SolSystemPdeIndex.size(); k++) {
     
      unsigned solIndex = _SolSystemPdeIndex[k];
      unsigned solPairIndex = _ml_sol->GetSolutionPairInverseIndex(k);
      
      unsigned kPair = GetSolPdeIndex(_ml_sol->GetSolName_from_index(solPairIndex));
            
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
	  bool solidMarki = amrSolidMark[solType].at(i);
	  const int ncols = amrRestriction[solType].at(i).size();
          std::vector <int> colPP(ncols);
	  std::vector <int> colRR(ncols);
          std::vector <double> valuePP(ncols);
	  std::vector <double> valueRR(ncols);
          unsigned j = 0;
          for(std::map<unsigned, double> ::const_iterator it = amrRestriction[solType].at(i).begin(); it != amrRestriction[solType].at(i).end(); it++) {
	    bool solidMarkj = amrSolidMark[solType].at(it->first);
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
              unsigned jproc = _msh[level]->BisectionSearch_find_processor_of_dof(it->first, solType);
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
    
    // ------------------- Prolongator - END
    
    
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

