

#include "fe_projection_matrices_Lagrange_continuous.hpp"
#include "Mesh.hpp"
#include "ElemType.hpp"

#include "NumericVector.hpp"
#include "SparseMatrix.hpp"


namespace femus {

    
    FE_Proj_Matrices::FE_Proj_Matrices() { 
        
    for(int itype = 0; itype < NFE_FAMS_C_ZERO_LAGRANGE; itype++) {
      for(int jtype = 0; jtype < NFE_FAMS_C_ZERO_LAGRANGE; jtype++) {
        _ProjQitoQj[itype][jtype] = NULL;
      }
    }

        
        
    }
    
    FE_Proj_Matrices::~FE_Proj_Matrices() {
        
    for(int itype = 0; itype < NFE_FAMS_C_ZERO_LAGRANGE; itype++) {
      for(int jtype = 0; jtype < NFE_FAMS_C_ZERO_LAGRANGE; jtype++) {
        if(_ProjQitoQj[itype][jtype]) {
          delete _ProjQitoQj[itype][jtype];
          _ProjQitoQj[itype][jtype] = NULL;
        }
      }
    }

        
    }
    
    
 
  SparseMatrix* FE_Proj_Matrices::GetQitoQjProjection(const unsigned& itype, const unsigned& jtype, const Mesh & mesh_in) {
      
    if(itype < NFE_FAMS_C_ZERO_LAGRANGE && jtype < NFE_FAMS_C_ZERO_LAGRANGE) {
      
      if(!_ProjQitoQj[itype][jtype]) {
        BuildQitoQjProjection(itype, jtype, mesh_in);
      }
      
    }
    else {
      std::cout << "Wrong argument range in function"
                << "FE_Proj_Matrices::GetLagrangeProjectionMatrix(const unsigned& itype, const unsigned& jtype)" << std::endl;
      abort();
    }

    return _ProjQitoQj[itype][jtype];
  }
  
  

  void FE_Proj_Matrices::BuildQitoQjProjection(const unsigned& itype, const unsigned& jtype, const Mesh & mesh_in) {

    assert(itype < NFE_FAMS_C_ZERO_LAGRANGE);
    assert(jtype < NFE_FAMS_C_ZERO_LAGRANGE);

    // ------------------- Sparsity pattern size - BEGIN
    const unsigned ni = mesh_in.dofmap_get_dof_offset(itype, _nprocs);
    const unsigned ni_loc = mesh_in.dofmap_get_own_size(itype, _iproc);

    const unsigned nj = mesh_in.dofmap_get_dof_offset(jtype, _nprocs);
    const unsigned nj_loc = mesh_in.dofmap_get_own_size(jtype, _iproc);

    NumericVector* NNZ_d = NumericVector::build().release();

    if(1 == _nprocs) {  // IF SERIAL
      NNZ_d->init(ni, ni_loc, false, SERIAL);
    }
    else {
      NNZ_d->init(ni, ni_loc, mesh_in.dofmap_get_ghost_dofs( itype, processor_id() ), false, GHOSTED);
    }

    NNZ_d->zero();

    NumericVector* NNZ_o = NumericVector::build().release();
    NNZ_o->init(*NNZ_d);
    NNZ_o->zero();

    
    
    for(unsigned isdom = _iproc; isdom < _iproc + 1; isdom++) {
      for(unsigned iel = mesh_in.GetElementOffset(isdom); iel < mesh_in.GetElementOffset(isdom + 1); iel++) {
        const short unsigned ielt = mesh_in.GetElementType(iel);
        const unsigned ndofs_itype_in = mesh_in.GetFiniteElement(ielt, itype)->GetBasis()->n_dofs();
            Get_QitoQjProjection_SparsityPatternSize_OneElement_OneFEFamily_Lagrange_Continuous(mesh_in, iel, NNZ_d, NNZ_o, itype, mesh_in.GetFiniteElement(ielt, jtype), ndofs_itype_in );
      }
    }

    NNZ_d->close();
    NNZ_o->close();

    const unsigned offset = mesh_in.dofmap_get_dof_offset(itype, _iproc);

    std::vector < int > nnz_d(ni_loc);
    std::vector < int > nnz_o(ni_loc);

    for(unsigned i = 0; i < ni_loc; i++) {
      nnz_d[i] = static_cast < int >((*NNZ_d)(offset + i));
      nnz_o[i] = static_cast < int >((*NNZ_o)(offset + i));
    }
    // ------------------- Sparsity pattern size - END

    
    // ------------------- Projection - BEGIN
    
    _ProjQitoQj[itype][jtype] = SparseMatrix::build().release();
    _ProjQitoQj[itype][jtype]->init(ni, nj, ni_loc, nj_loc, nnz_d, nnz_o);

    for(unsigned isdom = _iproc; isdom < _iproc + 1; isdom++) {
      for(unsigned iel = mesh_in.GetElementOffset(isdom); iel < mesh_in.GetElementOffset(isdom + 1); iel++) {
        const short unsigned ielt = mesh_in.GetElementType(iel);
        const unsigned ndofs_itype_in = mesh_in.GetFiniteElement(ielt, itype)->GetBasis()->n_dofs();
          Build_QitoQjProjection_OneElement_OneFEFamily_Lagrange_Continuous(mesh_in, iel, _ProjQitoQj[itype][jtype], NNZ_d, NNZ_o, itype, mesh_in.GetFiniteElement(ielt, jtype), ndofs_itype_in );
      }
    }

    _ProjQitoQj[itype][jtype]->close();
    // ------------------- Projection - END

    
    delete NNZ_d;
    delete NNZ_o;
    
    
  }



//----------------------------------------------------------------------------------------------------
//BEGIN build matrix sparsity pattern size and build prolongator matrix for solution printing
//----------------------------------------------------------------------------------------------------

  void FE_Proj_Matrices::Get_QitoQjProjection_SparsityPatternSize_OneElement_OneFEFamily_Lagrange_Continuous(const Mesh& mesh,
                                         const int& iel,
                                         NumericVector* NNZ_d,
                                         NumericVector* NNZ_o,
                                         const unsigned& itype,
                                         const elem_type * elem_type_in_jtype,
                                         const unsigned ndofs_itype_in) const
  {
    
        assert(itype < NFE_FAMS_C_ZERO_LAGRANGE);

      
    const unsigned soltype_in = elem_type_in_jtype->GetSolType();
    const basis * pt_basis_in = elem_type_in_jtype->GetBasis();
    const unsigned      ndofs = elem_type_in_jtype->GetNDofs();
    
    const unsigned      ndofs_Lagrange = ndofs_itype_in;
    
    const bool identity = ( ndofs_Lagrange <= ndofs ) ? true : false;
    
    for(int i = 0; i < ndofs_Lagrange; i++) {
      
      const int irow = mesh.GetSolutionDof(i, iel, itype);
      const int iproc = mesh.BisectionSearch_find_processor_of_dof(irow, itype);
      const int ncols = (identity) ? 1 : ndofs;
      
      unsigned counter_off_diag = 0;
      unsigned counter_all_diag = 0;
      
      for(int k = 0; k < ncols; k++) {
        
        const double phi = (identity) ? 1. : pt_basis_in->eval_phi( pt_basis_in->GetIND(k), pt_basis_in->GetXcoarse(i) );
        
        if(fabs(phi) > 1.0e-14) {
          counter_all_diag ++;
          int kcolumn = (identity) ? mesh.GetSolutionDof(i, iel, soltype_in) : mesh.GetSolutionDof(k, iel, soltype_in);
          if( kcolumn <  mesh.dofmap_get_dof_offset(soltype_in, iproc)       || 
              kcolumn >= mesh.dofmap_get_dof_offset(soltype_in, iproc + 1) ) counter_off_diag++;
        }
        
      }
      
      NNZ_d->set(irow, counter_all_diag - counter_off_diag);
      NNZ_o->set(irow, counter_off_diag);
    }
    
    
  }
  

  void FE_Proj_Matrices::Build_QitoQjProjection_OneElement_OneFEFamily_Lagrange_Continuous(const Mesh& mesh,
                                    const int& iel,
                                    SparseMatrix* Projmat, 
                                    NumericVector* NNZ_d,
                                    NumericVector* NNZ_o,
                                    const unsigned& itype,
                                  const elem_type * elem_type_in_jtype,
                                         const unsigned ndofs_itype_in) const
  {
    
        assert(itype < NFE_FAMS_C_ZERO_LAGRANGE);
    
    const unsigned soltype_in = elem_type_in_jtype->GetSolType();
    const basis * pt_basis_in = elem_type_in_jtype->GetBasis();
    const unsigned      ndofs = elem_type_in_jtype->GetNDofs();
    
    const unsigned      ndofs_Lagrange = ndofs_itype_in;

    bool identity = ( ndofs_Lagrange <= ndofs ) ? true : false;
    
    std::vector<int> cols( ndofs );
    std::vector<double> value( ndofs );
    
    for(int i = 0; i < ndofs_Lagrange; i++) {
      int irow = mesh.GetSolutionDof(i, iel, itype);
      int ncols = (identity) ? 1 : ndofs;
      unsigned counter = 0;
      cols.resize( ndofs );
      for(int k = 0; k < ncols; k++) {
        
        const double phi = (identity) ? 1. : pt_basis_in->eval_phi(pt_basis_in->GetIND(k), pt_basis_in->GetXcoarse(i));
        
        if(fabs(phi) > 1.0e-14) {
          cols[counter]  = (identity) ? mesh.GetSolutionDof(i, iel, soltype_in) : mesh.GetSolutionDof(k, iel, soltype_in);
          value[counter] = phi;
          counter++;
        }
      }
      
      cols.resize(counter);
      int ncols_stored = static_cast <int>(floor((*NNZ_d)(irow) + (*NNZ_o)(irow) + 0.5));
      if(counter == ncols_stored) {
        Projmat->insert_row(irow, counter, cols, &value[0]);
      }
      
    }
    
    
  }

//----------------------------------------------------------------------------------------------------
//END build matrix sparsity pattern size and build prolongator matrix for solution printing
//----------------------------------------------------------------------------------------------------
   
    
}
