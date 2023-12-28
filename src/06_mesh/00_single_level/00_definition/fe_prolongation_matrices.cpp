
#include "fe_prolongation_matrices.hpp"

#include "Mesh.hpp"
#include "ElemType.hpp"

#include "NumericVector.hpp"
#include "SparseMatrix.hpp"

#include <cstring>


namespace femus {


    FE_Prolongation_Matrices::FE_Prolongation_Matrices() { 


    _coarseMsh = NULL;

    for(int i = 0; i < NFE_FAMS; i++) {
      _ProjCoarseToFine[i] = NULL;
    }
    
    }
    
    
    
     FE_Prolongation_Matrices::~FE_Prolongation_Matrices() { 
   
       for(unsigned i = 0; i < NFE_FAMS; i++) {
      if(_ProjCoarseToFine[i]) {
        delete _ProjCoarseToFine[i];
        _ProjCoarseToFine[i] = NULL;
      }
    }

         
     }
    
 
 
 
   
  SparseMatrix* FE_Prolongation_Matrices::GetCoarseToFineProjectionRestrictionOnCoarse(const unsigned& solType, const Mesh & mesh_in) {

    if(solType >= NFE_FAMS) {
      std::cout << "Wrong argument range in function \"GetCoarseToFineProjection\": "
                << "solType is greater then SolTypeMax" << std::endl;
      abort();
    }

    if(_ProjCoarseToFine[solType])
      BuildCoarseToFineProjection(solType, "coarse", mesh_in);

    return _ProjCoarseToFine[solType];
  }


  SparseMatrix* FE_Prolongation_Matrices::GetCoarseToFineProjection(const unsigned& solType, const Mesh & mesh_in) {

    if(solType >= NFE_FAMS) {
      std::cout << "Wrong argument range in function \"GetCoarseToFineProjection\": "
                << "solType is greater then SolTypeMax" << std::endl;
      abort();
    }

    if(!_ProjCoarseToFine[solType])
      BuildCoarseToFineProjection(solType, "fine", mesh_in);

    return _ProjCoarseToFine[solType];
  }



  void FE_Prolongation_Matrices::BuildCoarseToFineProjection(const unsigned& solType, const char el_dofs[], const Mesh & mesh_in) {

    if(!_coarseMsh) {
      std::cout << "Error! In function \"BuildCoarseToFineProjection\": the coarse mesh has not been set" << std::endl;
      abort();
    }

    if(!_ProjCoarseToFine[solType]) {

    // ------------------- Sparsity pattern size - BEGIN
      const int nf     = mesh_in.dofmap_get_dof_offset(solType, _nprocs);
      const int nc     = _coarseMsh->dofmap_get_dof_offset(solType, _nprocs);
      const int nf_loc = mesh_in.dofmap_get_own_size(solType, _iproc);
      const int nc_loc = _coarseMsh->dofmap_get_own_size(solType, _iproc);

      //build matrix sparsity pattern size
      NumericVector* NNZ_d = NumericVector::build().release();

      if(n_processors() == 1) {  // IF SERIAL
        NNZ_d->init(nf, nf_loc, false, SERIAL);
      }
      else { // IF PARALLEL
        if(solType < NFE_FAMS_C_ZERO_LAGRANGE) {  // GHOST nodes only for Lagrange FE families
          NNZ_d->init(nf, nf_loc, mesh_in.dofmap_get_ghost_dofs(solType, processor_id() ), false, GHOSTED);
        }
        else { //piecewise discontinuous variables have no ghost nodes
          NNZ_d->init(nf, nf_loc, false, PARALLEL);
        }
      }

      NNZ_d->zero();

      NumericVector* NNZ_o = NumericVector::build().release();
      NNZ_o->init(*NNZ_d);
      NNZ_o->zero();

      for(int isdom = _iproc; isdom < _iproc + 1; isdom++) {
        for(int iel = _coarseMsh->GetElementOffset(isdom); iel < _coarseMsh->GetElementOffset(isdom + 1); iel++) {
          const short unsigned ielt = _coarseMsh->GetElementType(iel);
            Get_Prolongation_SparsityPatternSize_OneElement_OneFEFamily( mesh_in, *_coarseMsh, iel, NNZ_d, NNZ_o, el_dofs, mesh_in.GetFiniteElement(ielt, solType) );
            ///@todo this GetFiniteElement should come from the Coarse one, it is more consistent, but it should be abstract so it should be the same, try to compare both
        }
      }

      NNZ_d->close();
      NNZ_o->close();

      const unsigned offset = mesh_in.dofmap_get_dof_offset(solType, _iproc);
      std::vector <int> nnz_d(nf_loc);
      std::vector <int> nnz_o(nf_loc);

      for(int i = 0; i < nf_loc; i++) {
        nnz_d[i] = static_cast <int>(floor((*NNZ_d)(offset + i) + 0.5));
        nnz_o[i] = static_cast <int>(floor((*NNZ_o)(offset + i) + 0.5));
      }

      delete NNZ_d;
      delete NNZ_o;
    // ------------------- Sparsity pattern size - END


      
    // ------------------- Prolongator - BEGIN
      //build matrix
      _ProjCoarseToFine[solType] = SparseMatrix::build().release();
      _ProjCoarseToFine[solType]->init(nf, nc, nf_loc, nc_loc, nnz_d, nnz_o);

      // loop on the coarse grid
      for(int isdom = _iproc; isdom < _iproc + 1; isdom++) {
        for(int iel = _coarseMsh->GetElementOffset(isdom); iel < _coarseMsh->GetElementOffset(isdom + 1); iel++) {
          const short unsigned ielt = _coarseMsh->GetElementType(iel);
            Build_Prolongation_OneElement_OneFEFamily(mesh_in, *_coarseMsh, iel, _ProjCoarseToFine[solType], el_dofs, mesh_in.GetFiniteElement(ielt, solType) );         }
      }

      _ProjCoarseToFine[solType]->close();
    // ------------------- Prolongator - END
      
    }

  }

  
  
//----------------------------------------------------------------------------------------------------
//BEGIN  build matrix sparsity pattern size and build prolongator matrix for single solution
//-----------------------------------------------------------------------------------------------------

  void FE_Prolongation_Matrices::Get_Prolongation_SparsityPatternSize_OneElement_OneFEFamily(const Mesh& meshf,
                                         const Mesh& meshc,
                                         const int& ielc,
                                         NumericVector* NNZ_d,
                                         NumericVector* NNZ_o,
                                         const char is_fine_or_coarse[],
                                         const elem_type * elem_type_in) const
  {

    const unsigned soltype_in = elem_type_in->GetSolType();
    const unsigned      ndofs = elem_type_in->GetNDofs();
    const unsigned      ndofs_fine = elem_type_in->GetNDofsFine();
    
      unsigned n_elemdofs = 0;
      if ( !strcmp(is_fine_or_coarse, "fine") )        n_elemdofs = ndofs_fine;
      else if ( !strcmp(is_fine_or_coarse, "coarse") ) n_elemdofs = ndofs;
      
    if(meshc.GetRefinedElementIndex(ielc)) {  // coarse2fine prolongation
      
      for(int i = 0; i < n_elemdofs ; i++) {
        
        const std::pair<int, int> id_0_1 = elem_type_in->GetKVERT_IND(i);
        
        const int irow = meshf.GetSolutionDof(ielc, id_0_1.first, id_0_1.second, soltype_in, &meshc);

        const int iproc = meshf.BisectionSearch_find_processor_of_dof(irow, soltype_in);
        
        const int ncols = elem_type_in->Get_Prolongator_Num_Columns(i);
        
        unsigned counter_off_diag = 0;

        for(int k = 0; k < ncols; k++) {
          int j = elem_type_in->Get_Prolongator_Index(i, k);
          int jcolumn = meshc.GetSolutionDof(j, ielc, soltype_in);

          if(jcolumn <  meshc.dofmap_get_dof_offset(soltype_in, iproc)     ||
             jcolumn >= meshc.dofmap_get_dof_offset(soltype_in, iproc + 1)   ) counter_off_diag++;
        }

        NNZ_d->set(irow, ncols - counter_off_diag);
        NNZ_o->set(irow, counter_off_diag);
      }
      
    }
    else { // coarse2coarse prolongation
      
      for(int i = 0; i < ndofs; i++) {
        
        int irow = meshf.GetSolutionDof(ielc, 0, i , soltype_in, &meshc);

        int iproc = meshf.BisectionSearch_find_processor_of_dof(irow, soltype_in);
        int jcolumn = meshc.GetSolutionDof(i, ielc, soltype_in);

        if( jcolumn <  meshc.dofmap_get_dof_offset(soltype_in, iproc)      || 
            jcolumn >= meshc.dofmap_get_dof_offset(soltype_in, iproc + 1)     ) {
          NNZ_o->set(irow, 1);
        }
        else {
          NNZ_d->set(irow, 1);
        }
      }
      
    }
    
    
  }
  
  

  void FE_Prolongation_Matrices::Build_Prolongation_OneElement_OneFEFamily(const Mesh& meshf,
                                    const Mesh& meshc, 
                                    const int& ielc,
                                    SparseMatrix* Projmat, 
                                    const char is_fine_or_coarse[],
                                  const elem_type * elem_type_in) const
  {
 
      const unsigned soltype_in = elem_type_in->GetSolType();
      const unsigned      ndofs = elem_type_in->GetNDofs();
      const unsigned      ndofs_fine = elem_type_in->GetNDofsFine();
    
      unsigned n_elemdofs = 0;
      if ( !strcmp(is_fine_or_coarse, "fine") )        n_elemdofs = ndofs_fine;
      else if ( !strcmp(is_fine_or_coarse, "coarse") ) n_elemdofs = ndofs;
      
      if(meshc.GetRefinedElementIndex(ielc)) {  // coarse2fine prolongation

      std::vector<int> jcols( ndofs );

      for(int i = 0; i < n_elemdofs /*_nf*/; i++) {
        
        const std::pair<int, int> id_0_1 = elem_type_in->GetKVERT_IND(i);
                
        const int irow = meshf.GetSolutionDof(ielc, id_0_1.first, id_0_1.second, soltype_in, &meshc);
        const int ncols  =  elem_type_in->Get_Prolongator_Num_Columns(i);
        
        jcols.assign(ncols, 0);

        for(int k = 0; k < ncols; k++) {
          int j = elem_type_in->Get_Prolongator_Index(i, k);
          int jcolumn = meshc.GetSolutionDof(j, ielc, soltype_in);
          jcols[k] = jcolumn;
        }

        Projmat->insert_row(irow, ncols, jcols, elem_type_in->Get_Prolongator_Values_Row(i) );
      }
      
    }
    else { // coarse2coarse prolongation
      
      std::vector <int> jcol(1);
      double one = 1.;

      for(int i = 0; i < ndofs; i++) {
        
        const int irow = meshf.GetSolutionDof(ielc, 0, i , soltype_in, &meshc);
        jcol[0] = meshc.GetSolutionDof(i, ielc, soltype_in);
        Projmat->insert_row(irow, 1, jcol, &one);
        
      }
      
    }
    
    
  }

//----------------------------------------------------------------------------------------------------
//END  build matrix sparsity pattern size and build prolongator matrix for single solution
//-----------------------------------------------------------------------------------------------------
 
 
 
 
 
    
}
