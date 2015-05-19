/*=========================================================================

 Program: FEMUS
 Module: ElemType
 Authors: Eugenio Aulisa, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

#include "GaussPoints.hpp"
#include "ElemType.hpp"
#include "FETypeEnum.hpp"
#include "Elem.hpp"
#include "NumericVector.hpp"

using std::cout;
using std::endl;


namespace femus {

 const unsigned elem_type::_fe_old_to_new[QL] = {2,0,3};
 const unsigned elem_type::_fe_new_to_old[NFE_FAMS] = {1,-7,0,2,-7};

unsigned elem_type::_refindex=1;

//   Constructor
  elem_type::elem_type(const char *geom_elem, const char *order_gauss) : _gauss(geom_elem, order_gauss) {
    isMpGDAllocated = false;
  }


elem_type::~elem_type() {
  delete [] _X;
  delete [] _KVERT_IND;
  delete [] _IND;

  delete [] _prol_val;
  delete [] _prol_ind;
  delete [] _mem_prol_val;
  delete [] _mem_prol_ind;

  delete _pt_basis;

  if(isMpGDAllocated){
    for (int g = 0; g < GetGaussRule().GetGaussPointsNumber(); g++) {
      delete [] _phi_mapGD[g];
      delete [] _dphidxez_mapGD[g];
    }
    delete [] _phi_mapGD;
    delete [] _dphidxez_mapGD;
  }

}


//----------------------------------------------------------------------------------------------------
// evaluate shape functions at all quadrature points  TODO DEALLOCATE at destructor TODO FEFamilies TODO change HEX27 connectivity
//-----------------------------------------------------------------------------------------------------

void elem_type::EvaluateShapeAtQP(const std::string geomel_id_in,const std::string fe_in) {

// if (  (!strcmp(fe_in.c_str(),"disc_linear"))  || (!strcmp(fe_in.c_str(),"quadratic")) ) {  std::cout << "BEWARE, family not supported yet" << std::endl; return; }

// ============== allocate canonical shape ==================================================================
         _phi_mapGD = new double*[GetGaussRule().GetGaussPointsNumber()];// TODO valgrind, remember to DEALLOCATE THESE, e.g. with smart pointers
    _dphidxez_mapGD = new double*[GetGaussRule().GetGaussPointsNumber()];

    for (int g = 0; g < GetGaussRule().GetGaussPointsNumber(); g++) {
           _phi_mapGD[g] = new double[GetNDofs()];
      _dphidxez_mapGD[g] = new double[GetNDofs()*GetDim()];
    }

    isMpGDAllocated = true;
// ============== allocate canonical shape ==================================================================


// HEX 27 CASE ==========================================
// HEX 27 CASE ==========================================
// HEX 27 CASE ==========================================
// from eu connectivity to my (=libmesh) connectivity
const unsigned from_femus_to_libmesh[27] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,24,20,21,22,23,25,26};
//                                          0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26
// from libmesh to eu connectivity
const unsigned from_libmesh_to_femus[27] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,21,22,23,24,20,25,26};

if ( (!strcmp(fe_in.c_str(),"biquadratic")) && GetDim() == 3  && (!strcmp(geomel_id_in.c_str(),"hex")) ) {
//             std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << "REMEMBER THAT ONLY HEX27 HAS A DIFFERENT CONNECTIVITY MAP"  << std::endl;

      for (int ig = 0; ig < GetGaussRule().GetGaussPointsNumber(); ig++) {

      for (int idof=0; idof < GetNDofs(); idof++) {
//                 std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << vb << " " << ig << " " << idof << std::endl;
        _phi_mapGD[ig][idof] = GetPhi(ig)[ from_femus_to_libmesh[idof] ];
// 	std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << vb << " " << ig << " " << dof << " phi " << _phi_mapGD[vb][ig][dof] << std::endl;

// derivatives in canonical element
        for (uint idim = 0; idim < GetDim(); idim++) {
          double* dphi_g =   ( this->*(_DPhiXiEtaZetaPtr[idim]) )(ig);  //how to access a pointer to member function
          _dphidxez_mapGD[ig][ idof + idim*GetNDofs()] =  dphi_g[ from_femus_to_libmesh[idof] ];
//           std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << " " << ig << " " << idof << " " << idim << " dphi         " << _dphidxez_mapGD[ig][ idof + idim*GetNDofs()]  << "                                      "  << std::endl;

        }

      }

    }  // end gauss



}
// HEX 27 CASE ==========================================
// HEX 27 CASE ==========================================
// HEX 27 CASE ==========================================

// ALL THE OTHERS ==========================================
// ALL THE OTHERS ==========================================
// ALL THE OTHERS ==========================================

   else {

    for (int ig = 0; ig < GetGaussRule().GetGaussPointsNumber(); ig++) {

      for (int idof=0; idof < GetNDofs(); idof++) {
        _phi_mapGD[ig][idof] = GetPhi(ig)[idof];

// derivatives in canonical element
        for (uint idim = 0; idim < GetDim(); idim++) {
          double* dphi_g = (this->*(_DPhiXiEtaZetaPtr[idim]) )(ig);  //how to access a pointer to member function
          _dphidxez_mapGD[ig][ idof + idim*GetNDofs()] =  dphi_g[idof];

        }

      }

    }  // end gauss


  } //else HEX27



  return;
}


//----------------------------------------------------------------------------------------------------
// build matrix sparsity pattern size and build prolungator matrix for the LsysPde  Matrix
//-----------------------------------------------------------------------------------------------------

void elem_type::BuildProlongation(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc, SparseMatrix* Projmat,
				  const unsigned &index_sol, const unsigned &kkindex_sol) const {
  vector<int> cols(27);

  for (int i=0; i<_nf; i++) {
    int i0=_KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=lspdec._msh->el->GetChildElement(ielc,i0);
    int i1=_KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=lspdef._msh->el->GetMeshDof(ielf,i1,_SolType);
    int irow=lspdef.GetKKDof(index_sol,kkindex_sol,iadd);  //  local-id to dof
    int ncols=_prol_ind[i+1]-_prol_ind[i];
    cols.assign(ncols,0);
    for (int k=0; k<ncols; k++) {
      int j=_prol_ind[i][k];
      int jadd=lspdec._msh->el->GetMeshDof(ielc,j,_SolType);
      int jj=lspdec.GetKKDof(index_sol,kkindex_sol,jadd);
      cols[k]=jj;
    }
    Projmat->insert_row(irow,ncols,cols,_prol_val[i]);
  }
}


void elem_type::GetSparsityPatternSize(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc,
				       NumericVector* NNZ_d, NumericVector* NNZ_o,
				       const unsigned &index_sol, const unsigned &kkindex_sol) const {

  for (int i=0; i<_nf; i++) {
    int i0=_KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=lspdec._msh->el->GetChildElement(ielc,i0);
    int i1=_KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=lspdef._msh->el->GetMeshDof(ielf,i1,_SolType);
    int irow=lspdef.GetKKDof(index_sol,kkindex_sol,iadd);  //  local-id to dof

    int iproc=0;
    while (irow < lspdef.KKoffset[0][iproc] || irow >= lspdef.KKoffset[lspdef.KKIndex.size()-1][iproc] ) iproc++;
    int ncols=_prol_ind[i+1]-_prol_ind[i];

    int counter_o=0;
    for (int k=0; k<ncols; k++) {
      int j=_prol_ind[i][k];
      int jadd=lspdec._msh->el->GetMeshDof(ielc,j,_SolType);
      int jcolumn=lspdec.GetKKDof(index_sol,kkindex_sol,jadd);
      if(jcolumn < lspdec.KKoffset[0][iproc] || jcolumn >= lspdec.KKoffset[lspdef.KKIndex.size()-1][iproc] ) counter_o++;
    }

    NNZ_d->set(irow,ncols-counter_o);
    NNZ_o->set(irow,counter_o);

  }
}

void elem_type::BuildRestrictionTranspose(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc, SparseMatrix* Projmat,
					  const unsigned &index_sol, const unsigned &kkindex_sol,
					  const unsigned &index_pair_sol, const unsigned &kkindex_pair_sol) const{

  vector<int> cols(27);
  bool fluid_region = (2==lspdec._msh->el->GetElementMaterial(ielc))?1:0;

  vector <double> copy_prol_val;
  copy_prol_val.reserve(27);
  for (int i=0; i<_nf; i++) {
    int i0=_KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=lspdec._msh->el->GetChildElement(ielc,i0);
    int i1=_KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=lspdef._msh->el->GetMeshDof(ielf,i1,_SolType);

    int irow=lspdef.GetKKDof(index_sol,kkindex_sol,iadd);  //  local-id to dof
    int ncols=_prol_ind[i+1]-_prol_ind[i];

    bool isolidmark=lspdef._msh->el->GetNodeRegion(iadd);

    cols.assign(ncols,0);
    copy_prol_val.resize(ncols);
    for (int k=0; k<ncols; k++) {
      int j=_prol_ind[i][k];
      int jadd=lspdec._msh->el->GetMeshDof(ielc,j,_SolType);
      bool jsolidmark=lspdef._msh->el->GetNodeRegion(jadd);
      if(isolidmark == jsolidmark){
	int jcolumn=lspdec.GetKKDof(index_sol,kkindex_sol,jadd);
	cols[k]=jcolumn;
	copy_prol_val[k]=_prol_val[i][k];
      }
      else {
	int jcolumn=lspdec.GetKKDof(index_pair_sol,kkindex_pair_sol,jadd);
	cols[k]=jcolumn;
	copy_prol_val[k]=(index_sol != index_pair_sol) ? _prol_val[i][k]:0.;
      }


      //copy_prol_val[k]=(!fluid_region || isolidmark==jsolidmark)?_prol_val[i][k]:0.;
    }
    Projmat->insert_row(irow,ncols,cols,&copy_prol_val[0]);
  }
}



//----------------------------------------------------------------------------------------------------
// build matrix sparsity pattern size and build prolungator matrix for single solution
//-----------------------------------------------------------------------------------------------------

void elem_type::GetSparsityPatternSize(const Mesh &meshf,const Mesh &meshc, const int& ielc, NumericVector* NNZ_d, NumericVector* NNZ_o) const {

  for (int i=0; i<_nf; i++) {
    int i0=_KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=meshc.el->GetChildElement(ielc,i0);
    int i1=_KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=meshf.el->GetMeshDof(ielf,i1,_SolType);
    int irow=meshf.GetMetisDof(iadd,_SolType);  //  local-id to dof
    int iproc=0;
    while (irow < meshf.MetisOffset[_SolType][iproc] || irow >= meshf.MetisOffset[_SolType][iproc+1] ) iproc++;

    int ncols=_prol_ind[i+1]-_prol_ind[i];
    unsigned counter_o=0;
    for (int k=0; k<ncols; k++) {
      int j=_prol_ind[i][k];
      int jadd=meshc.el->GetMeshDof(ielc,j,_SolType);
      int jcolumn=meshc.GetMetisDof(jadd,_SolType);
      if(jcolumn < meshc.MetisOffset[_SolType][iproc] || jcolumn >= meshc.MetisOffset[_SolType][iproc+1] ) counter_o++;
    }

    NNZ_d->set(irow,ncols-counter_o);
    NNZ_o->set(irow,counter_o);

  }
}

void elem_type::BuildProlongation(const Mesh &meshf,const Mesh &meshc, const int& ielc,
				  SparseMatrix* Projmat) const {
  vector<int> cols(27);
  for (int i=0; i<_nf; i++) {
    int i0=_KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=meshc.el->GetChildElement(ielc,i0);
    int i1=_KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=meshf.el->GetMeshDof(ielf,i1,_SolType);
    int irow=meshf.GetMetisDof(iadd,_SolType);  //  local-id to dof
    int ncols=_prol_ind[i+1]-_prol_ind[i];
    cols.assign(ncols,0);
    for (int k=0; k<ncols; k++) {
      int j=_prol_ind[i][k];
      int jadd=meshc.el->GetMeshDof(ielc,j,_SolType);
      int jcolumn=meshc.GetMetisDof(jadd,_SolType);
      cols[k]=jcolumn;
    }

    Projmat->insert_row(irow,ncols,cols,_prol_val[i]);
  }
}

//----------------------------------------------------------------------------------------------------
// prolungator for solution printing
//----------------------------------------------------------------------------------------------------

void elem_type::GetSparsityPatternSize(const Mesh& mesh,const int& iel, NumericVector* NNZ_d, NumericVector* NNZ_o, const unsigned &itype) const{
  for (int i=0; i<_nlag[itype]; i++) {
    int inode=mesh.el->GetMeshDof(iel,i,_SolType);
    int irow=mesh.GetMetisDof(inode,itype);
    int iproc=0;
    while (irow < mesh.MetisOffset[itype][iproc] || irow >= mesh.MetisOffset[itype][iproc+1] ) iproc++;
    int ncols=_prol_ind[i+1]-_prol_ind[i];
    unsigned counter_o=0;
    for (int k=0; k<ncols; k++) {
      int jj=_prol_ind[i][k];
      int jnode   = mesh.el->GetMeshDof(iel,jj,_SolType);
      int jcolumn = mesh.GetMetisDof(jnode,_SolType);
      if(jcolumn < mesh.MetisOffset[_SolType][iproc] || jcolumn >= mesh.MetisOffset[_SolType][iproc+1] ) counter_o++;
    }
    NNZ_d->set(irow,ncols-counter_o);
    NNZ_o->set(irow,counter_o);
  }
}


void elem_type::BuildProlongation(const Mesh& mesh,const int& iel, SparseMatrix* Projmat,const unsigned &itype) const{
  vector<int> cols(27);
  for (int i=0; i<_nlag[itype]; i++) {
    int inode=mesh.el->GetMeshDof(iel,i,_SolType);
    int irow=mesh.GetMetisDof(inode,itype);
    int ncols=_prol_ind[i+1]-_prol_ind[i];
    cols.assign(ncols,0);
    for (int k=0; k<ncols; k++) {
      int jj=_prol_ind[i][k];
      int jnode=mesh.el->GetMeshDof(iel,jj,_SolType);
      int jcolumn=mesh.GetMetisDof(jnode,_SolType);
      cols[k]=jcolumn;
    }
    Projmat->insert_row(irow,ncols,cols,_prol_val[i]);
  }
}

elem_type_1D::elem_type_1D(const char *geom_elem, const char *order, const char *order_gauss) :
	      elem_type(geom_elem,order_gauss) {

  _dim = 1;
  _DPhiXiEtaZetaPtr.resize(_dim);
  _DPhiXiEtaZetaPtr[0] = &elem_type::GetDPhiDXi;

  //************ BEGIN FE and MG SETUP ******************
  if (!strcmp(order,"linear")) 		 _SolType=0;
  else if (!strcmp(order,"quadratic")) 	 _SolType=1;
  else if (!strcmp(order,"biquadratic")) _SolType=2;
  else if (!strcmp(order,"constant"))    _SolType=3;
  else if (!strcmp(order,"disc_linear")) _SolType=4;
  else {
    cout << order << " is not a valid option for " << geom_elem << endl;
    exit(0);
  }

  if (!strcmp(geom_elem,"line")) { //line
    if 	    (_SolType == 0) _pt_basis = new line1;
    else if (_SolType == 1) _pt_basis = new line2;
    else if (_SolType == 2) _pt_basis = new line2;
    else if (_SolType == 3) _pt_basis = new line0;
    else if (_SolType == 4) _pt_basis = new linepwl;
    else {
      cout << order << " is not a valid option for " << geom_elem << endl;
      exit(0);
    }
  }
  else {
    cout << geom_elem << " is not a valid option" << endl;
    exit(0);
  }

  // get data from basis object
  _nc 	   = _pt_basis->_nc;
  _nf 	   = _pt_basis->_nf;
  _nlag[0] = _pt_basis->_nlag0;
  _nlag[1] = _pt_basis->_nlag1;
  _nlag[2] = _pt_basis->_nlag2;


  _IND=new const int * [_nc];
  for (int i=0; i<_nc; i++){
    _IND[i]=_pt_basis->getIND(i);
  }
  _KVERT_IND=new const int * [_nf];
  _X=new const double * [_nf];
  for (int i=0; i<_nf; i++) {
    _KVERT_IND[i] = _pt_basis->getKVERT_IND(i);
    _X[i] = _pt_basis->getX(i);
  }

  // local projection matrix evaluation
  int counter = 0;
  for (int i = 0; i < _nf; i++) {
    for (int j = 0; j < _nc; j++) {
      double phi = _pt_basis->eval_phi(_IND[j],_X[i]);
      if ( _SolType == 4 ) { //if piece_wise_linear
        if (i/2 == 1)  phi = _pt_basis->eval_dphidx(_IND[j],_X[i]);
      }
      if ( fabs( phi ) >= 1.0e-14 ){
        counter++;
      }
    }
  }
  double *pt_d;
  int *pt_i;

  _prol_val = new double * [_nf+1];
  _prol_ind = new int * [_nf+1];
  _mem_prol_val = new double [counter];
  _mem_prol_ind = new int [counter];

  pt_d= _mem_prol_val;
  pt_i= _mem_prol_ind;
  for (int i=0; i<_nf; i++) {
    _prol_val[i] = pt_d;
    _prol_ind[i] = pt_i;
    for (int j=0; j<_nc; j++) {
      double phi = _pt_basis->eval_phi(_IND[j],_X[i]);
      if ( _SolType==4 ) { //if piece_wise_linear
        if (i/2 == 1) phi = _pt_basis->eval_dphidx(_IND[j],_X[i])/2.;
      }
      if ( fabs(phi) >= 1.0e-14 ){
        *(pt_d++) = phi;
        *(pt_i++) = j;
      }
    }
  }

  _prol_val[_nf] = pt_d;
  _prol_ind[_nf] = pt_i;

  // shape function and its derivatives evaluated at Gauss'points
   int n_gauss = _gauss.GetGaussPointsNumber();

  _phi= new double*[n_gauss];
  _dphidxi  = new double*[n_gauss];
  _d2phidxi2  = new double*[n_gauss];

  _phi_memory=new double [n_gauss*_nc];
  _dphidxi_memory  =new double [n_gauss*_nc];
  _d2phidxi2_memory  =new double [n_gauss*_nc];

  for (unsigned i=0; i<n_gauss; i++) {
    _phi[i]=&_phi_memory[i*_nc];
    _dphidxi[i]  =&_dphidxi_memory[i*_nc];
    _d2phidxi2[i]  =&_d2phidxi2_memory[i*_nc];
  }

 const double *ptx[1]={_gauss.GetGaussWeightsPointer() + n_gauss};
  for (unsigned i=0; i<n_gauss; i++){
    double x[1];
    for (unsigned j=0; j<1;j++) {
      x[j] = *ptx[j];
      ptx[j]++;
    }

    for (int j=0; j<_nc; j++) {
      _phi[i][j] = _pt_basis->eval_phi(_IND[j],x);
      _dphidxi[i][j] = _pt_basis->eval_dphidx(_IND[j],x);
      _d2phidxi2[i][j] = _pt_basis->eval_d2phidx2(_IND[j],x);
    }
  }

 //=====================
  EvaluateShapeAtQP(geom_elem,order);


}


elem_type_2D::elem_type_2D(const char *geom_elem, const char *order, const char *order_gauss):
	      elem_type(geom_elem,order_gauss)
	       {

  _dim = 2;
  _DPhiXiEtaZetaPtr.resize(_dim);
  _DPhiXiEtaZetaPtr[0] = &elem_type::GetDPhiDXi;
  _DPhiXiEtaZetaPtr[1] = &elem_type::GetDPhiDEta;

  //************ BEGIN FE and MG SETUP ******************
  if 	  (!strcmp(order,"linear")) 	 _SolType=0;
  else if (!strcmp(order,"quadratic")) 	 _SolType=1;
  else if (!strcmp(order,"biquadratic")) _SolType=2;
  else if (!strcmp(order,"constant"))    _SolType=3;
  else if (!strcmp(order,"disc_linear")) _SolType=4;
  else {
    cout << order << " is not a valid option for " << geom_elem << endl;
    abort();
  }

  if (!strcmp(geom_elem,"quad")) { //QUAD
    if 	    (_SolType == 0) _pt_basis = new quad1;
    else if (_SolType == 1) _pt_basis = new quadth;
    else if (_SolType == 2) _pt_basis = new quad2;
    else if (_SolType == 3) _pt_basis = new quad0;
    else if (_SolType == 4) _pt_basis = new quadpwl;
    else {
      cout << order << " is not a valid option for " << geom_elem << endl;
      abort();
    }
  }
  else if (!strcmp(geom_elem,"tri")) { //TRIANGLE

    if 	    (_SolType == 0) _pt_basis = new tri1;
    else if (_SolType == 1) _pt_basis = new tri2;
    else if (_SolType == 2) _pt_basis = new tri2;
    else if (_SolType == 3) _pt_basis = new tri0;
    else if (_SolType == 4) _pt_basis = new tripwl;
    else {
    cout << order << " is not a valid option for " << geom_elem << endl;
    abort();
    }
  }
  else {
    cout << geom_elem << " is not a valid option" << endl;
    abort();
  }

  // get data from basis object
  _nc 	   = _pt_basis->_nc;
  _nf 	   = _pt_basis->_nf;
  _nlag[0] = _pt_basis->_nlag0;
  _nlag[1] = _pt_basis->_nlag1;
  _nlag[2] = _pt_basis->_nlag2;

  _IND=new const int * [_nc];
  for (int i=0; i<_nc; i++){
    _IND[i]=_pt_basis->getIND(i);
  }
  _KVERT_IND=new const int * [_nf];
  _X=new const double * [_nf];
  for (int i=0; i<_nf; i++) {
    _KVERT_IND[i]=_pt_basis->getKVERT_IND(i);
    _X[i]=_pt_basis->getX(i);
  }

  // local projection matrix evaluation
  int counter=0;
  for (int i=0; i<_nf; i++) {
    for (int j=0; j<_nc; j++) {
      double phi=_pt_basis->eval_phi(_IND[j],_X[i]);
      if ( _SolType==4 ) { //if piece_wise_linear
        if 	(i/4==1) phi=_pt_basis->eval_dphidx(_IND[j],_X[i]);
        else if (i/4==2) phi=_pt_basis->eval_dphidy(_IND[j],_X[i]);
      }
      if ( fabs(phi) >= 1.0e-14 ){
        counter++;
      }
    }
  }
  double *pt_d;
  int *pt_i;

  _prol_val=new double * [_nf+1];
  _prol_ind=new int * [_nf+1];
  _mem_prol_val=new double [counter];
  _mem_prol_ind=new int [counter];

  pt_d=_mem_prol_val;
  pt_i=_mem_prol_ind;
  for (int i=0; i<_nf; i++) {
    _prol_val[i]=pt_d;
    _prol_ind[i]=pt_i;
    for (int j=0; j<_nc; j++) {
      double phi=_pt_basis->eval_phi(_IND[j],_X[i]);
      if ( _SolType==4 ) { //if piece_wise_linear
        if 	(i/4==1)  phi=_pt_basis->eval_dphidx(_IND[j],_X[i])/2.;
        else if (i/4==2)  phi=_pt_basis->eval_dphidy(_IND[j],_X[i])/2.;
      }
      if ( fabs(phi) >= 1.0e-14 ){
        *(pt_d++)=phi;
        *(pt_i++)=j;
      }
    }
  }

  _prol_val[_nf]=pt_d;
  _prol_ind[_nf]=pt_i;

  // shape function and its derivatives evaluated at Gauss'points
   int n_gauss = _gauss.GetGaussPointsNumber();

  _phi= new double*[n_gauss];
  _dphidxi  = new double*[n_gauss];
  _dphideta = new double*[n_gauss];

  _d2phidxi2  = new double*[n_gauss];
  _d2phideta2 = new double*[n_gauss];

  _d2phidxideta  = new double*[n_gauss];

  _phi_memory=new double [n_gauss*_nc];
  _dphidxi_memory  =new double [n_gauss*_nc];
  _dphideta_memory =new double [n_gauss*_nc];

  _d2phidxi2_memory  =new double [n_gauss*_nc];
  _d2phideta2_memory =new double [n_gauss*_nc];

  _d2phidxideta_memory  =new double [n_gauss*_nc];

  for (unsigned i=0; i<n_gauss; i++) {
    _phi[i]=&_phi_memory[i*_nc];
    _dphidxi[i]  =&_dphidxi_memory[i*_nc];
    _dphideta[i] =&_dphideta_memory[i*_nc];

    _d2phidxi2[i]  =&_d2phidxi2_memory[i*_nc];
    _d2phideta2[i] =&_d2phideta2_memory[i*_nc];

    _d2phidxideta[i]  = &_d2phidxideta_memory[i*_nc];

  }

  const double *ptx[2]={_gauss.GetGaussWeightsPointer() + n_gauss, _gauss.GetGaussWeightsPointer() + 2*n_gauss};
  for (unsigned i=0; i<n_gauss; i++){
    double x[2];
    for (unsigned j=0; j<2;j++) {
      x[j] = *ptx[j];
      ptx[j]++;
    }
    for (int j=0; j<_nc; j++) {
      _phi[i][j] = _pt_basis->eval_phi(_IND[j],x);
      _dphidxi[i][j] = _pt_basis->eval_dphidx(_IND[j],x);
      _dphideta[i][j] = _pt_basis->eval_dphidy(_IND[j],x);
      _d2phidxi2[i][j] = _pt_basis->eval_d2phidx2(_IND[j],x);
      _d2phideta2[i][j] = _pt_basis->eval_d2phidy2(_IND[j],x);
      _d2phidxideta[i][j] = _pt_basis->eval_d2phidxdy(_IND[j],x);
    }
  }

 //=====================
  EvaluateShapeAtQP(geom_elem,order);


}

elem_type_3D::elem_type_3D(const char *geom_elem, const char *order, const char *order_gauss) :
	      elem_type(geom_elem,order_gauss)  {

  _dim = 3;
  _DPhiXiEtaZetaPtr.resize(_dim);
  _DPhiXiEtaZetaPtr[0] = &elem_type::GetDPhiDXi;
  _DPhiXiEtaZetaPtr[1] = &elem_type::GetDPhiDEta;
  _DPhiXiEtaZetaPtr[2] = &elem_type::GetDPhiDZeta;


  //************ BEGIN FE and MG SETUP ******************
  if 	  (!strcmp(order,"linear")) 	 _SolType=0;
  else if (!strcmp(order,"quadratic")) 	 _SolType=1;
  else if (!strcmp(order,"biquadratic")) _SolType=2;
  else if (!strcmp(order,"constant"))    _SolType=3;
  else if (!strcmp(order,"disc_linear")) _SolType=4;
  else {
    cout << order << " is not a valid option for " << geom_elem << endl;
    exit(0);
  }

  if (!strcmp(geom_elem,"hex")) {//HEX

    if 	    (_SolType == 0) _pt_basis = new hex1;
    else if (_SolType == 1) _pt_basis = new hexth;
    else if (_SolType == 2) _pt_basis = new hex2;
    else if (_SolType == 3) _pt_basis = new hex0;
    else if (_SolType == 4) _pt_basis = new hexpwl;
    else {
      cout << order << " is not a valid option for " << geom_elem << endl;
      exit(0);
    }
  }
  else if (!strcmp(geom_elem,"wedge")) { //WEDGE
    if 	    (_SolType == 0) _pt_basis = new wedge1;
    else if (_SolType == 1) _pt_basis = new wedgeth;
    else if (_SolType == 2) _pt_basis = new wedge2;
    else if (_SolType == 3) _pt_basis = new wedge0;
    else if (_SolType == 4) _pt_basis = new wedgepwl;
    else {
      cout << order << " is not a valid option for " << geom_elem << endl;
      exit(0);
    }
  }
  else if (!strcmp(geom_elem,"tet")) { //TETRAHEDRA
    if 	    (_SolType == 0) _pt_basis = new tet1;
    else if (_SolType == 1) _pt_basis = new tet2;
    else if (_SolType == 2) _pt_basis = new tet2;
    else if (_SolType == 3) _pt_basis = new tet0;
    else if (_SolType == 4) _pt_basis = new tetpwl;
    else {
      cout << order << " is not a valid option for " << geom_elem << endl;
      exit(0);
    }
  }
  else {
    cout << geom_elem << " is not a valid option" << endl;
    exit(0);
  }

  // get data from basis object
  _nc 	   = _pt_basis->_nc;
  _nf 	   = _pt_basis->_nf;
  _nlag[0] = _pt_basis->_nlag0;
  _nlag[1] = _pt_basis->_nlag1;
  _nlag[2] = _pt_basis->_nlag2;

  _IND=new const int * [_nc];
  for (int i=0; i<_nc; i++){
    _IND[i]=_pt_basis->getIND(i);
  }
  _KVERT_IND=new const int * [_nf];
  _X=new const double * [_nf];
  for (int i=0; i<_nf; i++) {
    _KVERT_IND[i]=_pt_basis->getKVERT_IND(i);
    _X[i]=_pt_basis->getX(i);
  }

  // local projection matrix evaluation
  int counter=0;
  for (int i=0; i<_nf; i++) {
    for (int j=0; j<_nc; j++) {
      double phi=_pt_basis->eval_phi(_IND[j],_X[i]);
      if ( _SolType==4 ) { //if piece_wise_linear
        if 	( i/8 == 1 ) phi=_pt_basis->eval_dphidx(_IND[j],_X[i]);
        else if ( i/8 == 2 ) phi=_pt_basis->eval_dphidy(_IND[j],_X[i]);
        else if ( i/8 == 3 ) phi=_pt_basis->eval_dphidz(_IND[j],_X[i]);
      }
      if ( fabs(phi) >= 1.0e-14 ){
        counter++;
      }
    }
  }
  double *pt_d;
  int *pt_i;

  _prol_val=new double * [_nf+1];
  _prol_ind=new int * [_nf+1];
  _mem_prol_val=new double [counter];
  _mem_prol_ind=new int [counter];

  pt_d=_mem_prol_val;
  pt_i=_mem_prol_ind;
  for (int i=0; i<_nf; i++) {
    _prol_val[i] = pt_d;
    _prol_ind[i] = pt_i;
    for (int j=0; j<_nc; j++) {
      double phi=_pt_basis->eval_phi(_IND[j],_X[i]);
      if ( _SolType==4 ) { //if piece_wise_linear
        if 	( i/8 == 1 ) phi = _pt_basis->eval_dphidx(_IND[j],_X[i])/2.;
        else if ( i/8 == 2 ) phi = _pt_basis->eval_dphidy(_IND[j],_X[i])/2.;
        else if ( i/8 == 3 ) phi = _pt_basis->eval_dphidz(_IND[j],_X[i])/2.;
      }
      if ( fabs(phi) >= 1.0e-14 ){
        *(pt_d++)=phi;
        *(pt_i++)=j;
      }
    }
  }

  _prol_val[_nf]=pt_d;
  _prol_ind[_nf]=pt_i;

  // shape function and its derivatives evaluated at Gauss'points
   int n_gauss = _gauss.GetGaussPointsNumber();

  _phi= new double*[n_gauss];
  _dphidxi  = new double*[n_gauss];
  _dphideta = new double*[n_gauss];
  _dphidzeta= new double*[n_gauss];

  _d2phidxi2  = new double*[n_gauss];
  _d2phideta2 = new double*[n_gauss];
  _d2phidzeta2= new double*[n_gauss];

  _d2phidxideta  = new double*[n_gauss];
  _d2phidetadzeta = new double*[n_gauss];
  _d2phidzetadxi= new double*[n_gauss];

  _phi_memory=new double [n_gauss*_nc];
  _dphidxi_memory  =new double [n_gauss*_nc];
  _dphideta_memory =new double [n_gauss*_nc];
  _dphidzeta_memory=new double [n_gauss*_nc];

  _d2phidxi2_memory  =new double [n_gauss*_nc];
  _d2phideta2_memory =new double [n_gauss*_nc];
  _d2phidzeta2_memory=new double [n_gauss*_nc];

  _d2phidxideta_memory  =new double [n_gauss*_nc];
  _d2phidetadzeta_memory =new double [n_gauss*_nc];
  _d2phidzetadxi_memory=new double [n_gauss*_nc];

  for (unsigned i=0; i< n_gauss; i++) {
    _phi[i]=&_phi_memory[i*_nc];
    _dphidxi[i]  =&_dphidxi_memory[i*_nc];
    _dphideta[i] =&_dphideta_memory[i*_nc];
    _dphidzeta[i]=&_dphidzeta_memory[i*_nc];

    _d2phidxi2[i]  =&_d2phidxi2_memory[i*_nc];
    _d2phideta2[i] =&_d2phideta2_memory[i*_nc];
    _d2phidzeta2[i]=&_d2phidzeta2_memory[i*_nc];

    _d2phidxideta[i]  = &_d2phidxideta_memory[i*_nc];
    _d2phidetadzeta[i]= &_d2phidetadzeta_memory[i*_nc];
    _d2phidzetadxi[i] = &_d2phidzetadxi_memory[i*_nc];

  }

  const double *ptx[3]={_gauss.GetGaussWeightsPointer() +   n_gauss,
                        _gauss.GetGaussWeightsPointer() + 2*n_gauss,
                        _gauss.GetGaussWeightsPointer() + 3*n_gauss};
  for (unsigned i=0; i<n_gauss; i++){
    double x[3];
    for (unsigned j=0; j<3;j++) {
      x[j] = *ptx[j];
      ptx[j]++;
    }

    for (int j=0; j<_nc; j++) {
      _phi[i][j] = _pt_basis->eval_phi(_IND[j],x);
      _dphidxi[i][j] = _pt_basis->eval_dphidx(_IND[j],x);
      _dphideta[i][j] = _pt_basis->eval_dphidy(_IND[j],x);
      _dphidzeta[i][j] = _pt_basis->eval_dphidz(_IND[j],x);
      _d2phidxi2[i][j] = _pt_basis->eval_d2phidx2(_IND[j],x);
      _d2phideta2[i][j] = _pt_basis->eval_d2phidy2(_IND[j],x);
      _d2phidzeta2[i][j] = _pt_basis->eval_d2phidz2(_IND[j],x);
      _d2phidxideta[i][j] = _pt_basis->eval_d2phidxdy(_IND[j],x);
      _d2phidetadzeta[i][j] = _pt_basis->eval_d2phidydz(_IND[j],x);
      _d2phidzetadxi[i][j] = _pt_basis->eval_d2phidzdx(_IND[j],x);
    }
  }


 //=====================
  EvaluateShapeAtQP(geom_elem,order);

}

//---------------------------------------------------------------------------------------------------------

template <class type>
void elem_type_1D::Jacobian_type(const vector < vector < type > > &vt,const unsigned &ig, type &Weight,
		   vector < double > &phi, vector < type > &gradphi, vector < type > &nablaphi) const{

  phi.resize(_nc);
  gradphi.resize(_nc*1);
  nablaphi.resize(_nc*1);

  type Jac=0.;
  type JacI;

  const double *dxi=_dphidxi[ig];

  for (int inode=0; inode<_nc; inode++,dxi++) {
    Jac+=(*dxi)*vt[0][inode];
  }

  Weight=Jac*_gauss.GetGaussWeightsPointer()[ig];

  JacI=1/Jac;

  dxi = _dphidxi[ig];
  const double *dxi2 = _d2phidxi2[ig];

  for (int inode=0; inode<_nc; inode++,dxi++, dxi2++) {
    phi[inode]=_phi[ig][inode];
    gradphi[inode]=(*dxi)*JacI;
    nablaphi[inode] = (*dxi2)*JacI*JacI;
  }

}

//---------------------------------------------------------------------------------------------------------

template <class type>
void elem_type_1D::JacobianSur_type(const vector < vector < type > > &vt, const unsigned &ig, type &Weight,
				  vector < double > &phi, vector < type > &gradphi, vector < type > &normal) const {

  phi.resize(_nc);
  normal.resize(2);

  type Jac[2][2]={{0.,0.},{0.,0.}};
  type JacI[2][2];

  const double *dfeta=_dphidxi[ig];

  for (int inode=0; inode<_nc; inode++,dfeta++) {
    Jac[0][0] += (*dfeta)*vt[0][inode];
    Jac[1][0] += (*dfeta)*vt[1][inode];
  }

//   normal module
  type modn = sqrt(Jac[0][0]*Jac[0][0] + Jac[1][0]*Jac[1][0]);

  normal[0] =  Jac[1][0]/modn;
  normal[1] = -Jac[0][0]/modn;

  //The derivative of x with respect to eta (dx/deta) has the opposite sign with respect to the normal
  //obtained as cross product between (dx/deta , dy/deta, 0) x (0,0,1)
  //The Jacobian has the structure
  // |dx/deta  -nx|
  // |dy/deta  -ny|
  Jac[0][1] = -normal[0];
  Jac[1][1] = -normal[1];

  //The determinant of that matrix is the area
  type det= (Jac[0][0]*Jac[1][1]-Jac[0][1]*Jac[1][0]);

  JacI[0][0] =  Jac[1][1]/det;
  JacI[0][1] = -Jac[0][1]/det;
  JacI[1][0] = -Jac[1][0]/det;
  JacI[1][1] =  Jac[0][0]/det;

  Weight = det*_gauss.GetGaussWeightsPointer()[ig];

  for(int inode=0;inode<_nc;inode++){
    phi[inode]=_phi[ig][inode];
  }

}

//---------------------------------------------------------------------------------------------------------

template <class type>
void elem_type_2D::Jacobian_type(const vector < vector < type > > &vt,const unsigned &ig, type &Weight,
				 vector < double > &phi, vector < type > &gradphi, vector < type > &nablaphi) const{

  phi.resize(_nc);
  gradphi.resize(_nc*2);
  nablaphi.resize(_nc*3);

  type Jac[2][2]={{0,0},{0,0}};
  type JacI[2][2];
  const double *dxi=_dphidxi[ig];
  const double *deta=_dphideta[ig];
  for (int inode=0; inode<_nc; inode++,dxi++,deta++){
    Jac[0][0] += (*dxi)*vt[0][inode];
    Jac[0][1] += (*dxi)*vt[1][inode];
    Jac[1][0] += (*deta)*vt[0][inode];
    Jac[1][1] += (*deta)*vt[1][inode];
  }
  type det=(Jac[0][0]*Jac[1][1]-Jac[0][1]*Jac[1][0]);

  JacI[0][0]= Jac[1][1]/det;
  JacI[0][1]=-Jac[0][1]/det;
  JacI[1][0]=-Jac[1][0]/det;
  JacI[1][1]= Jac[0][0]/det;

  Weight=det*_gauss.GetGaussWeightsPointer()[ig];

  dxi=_dphidxi[ig];
  deta=_dphideta[ig];

  const double *dxi2=_d2phidxi2[ig];
  const double *deta2=_d2phideta2[ig];
  const double *dxideta=_d2phidxideta[ig];

  for (int inode=0; inode<_nc; inode++, dxi++, deta++, dxi2++, deta2++, dxideta++) {

    phi[inode]=_phi[ig][inode];

    gradphi[2*inode+0]=(*dxi)*JacI[0][0] + (*deta)*JacI[0][1];
    gradphi[2*inode+1]=(*dxi)*JacI[1][0] + (*deta)*JacI[1][1];

    nablaphi[3*inode+0]=
      ( (*dxi2)   *JacI[0][0] + (*dxideta)*JacI[0][1] ) * JacI[0][0] +
      ( (*dxideta)*JacI[0][0] + (*deta2)  *JacI[0][1] ) * JacI[0][1];
    nablaphi[3*inode+1]=
      ( (*dxi2)   *JacI[1][0] + (*dxideta)*JacI[1][1] ) * JacI[1][0] +
      ( (*dxideta)*JacI[1][0] + (*deta2)  *JacI[1][1] ) * JacI[1][1];
    nablaphi[3*inode+2]=
      ( (*dxi2)   *JacI[0][0] + (*dxideta)*JacI[0][1] ) * JacI[1][0] +
      ( (*dxideta)*JacI[0][0] + (*deta2)  *JacI[0][1] ) * JacI[1][1];

  }
}

//---------------------------------------------------------------------------------------------------------

template <class type>
void elem_type_2D::JacobianSur_type(const vector < vector < type > > &vt, const unsigned &ig, type &Weight,
				    vector < double > &phi, vector < type > &gradphi, vector < type > &normal) const {
  phi.resize(_nc);
  normal.resize(3);

  type Jac[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};

  const double *dfx=_dphidxi[ig];
  const double *dfy=_dphideta[ig];

  for(int inode=0; inode<_nc; inode++,dfx++,dfy++){
    Jac[0][0] += (*dfx)*vt[0][inode];
    Jac[1][0] += (*dfx)*vt[1][inode];
    Jac[2][0] += (*dfx)*vt[2][inode];

    Jac[0][1] += (*dfy)*vt[0][inode];
    Jac[1][1] += (*dfy)*vt[1][inode];
    Jac[2][1] += (*dfy)*vt[2][inode];
  }

    //   normal module
  type nx = Jac[1][0]*Jac[2][1] - Jac[1][1]*Jac[2][0];
  type ny = Jac[0][1]*Jac[2][0] - Jac[2][1]*Jac[0][0];
  type nz = Jac[0][0]*Jac[1][1] - Jac[0][1]*Jac[1][0];
  type invModn = 1./sqrt(nx*nx + ny*ny + nz*nz);

  normal[0] =  (nx)*invModn;
  normal[1] =  (ny)*invModn;
  normal[2] =  (nz)*invModn;

  Jac[0][2] = normal[0];
  Jac[1][2] = normal[1];
  Jac[2][2] = normal[2];

  //the determinant of the matrix is the area
  type det=(Jac[0][0]*(Jac[1][1]*Jac[2][2]-Jac[1][2]*Jac[2][1])+
		      Jac[0][1]*(Jac[1][2]*Jac[2][0]-Jac[1][0]*Jac[2][2])+
		      Jac[0][2]*(Jac[1][0]*Jac[2][1]-Jac[1][1]*Jac[2][0]));

  Weight=det*_gauss.GetGaussWeightsPointer()[ig];

  for(int inode=0;inode<_nc;inode++){
    phi[inode]=_phi[ig][inode];
  }

  //TODO warning the surface gradient is missing!!!!!!!!!!!!!!!
}


//---------------------------------------------------------------------------------------------------------
template <class type>
void elem_type_3D::Jacobian_type(const vector < vector < type > > &vt,const unsigned &ig, type &Weight,
		   vector < double > &phi, vector < type > &gradphi, vector < type > &nablaphi) const{

  phi.resize(_nc);
  gradphi.resize(_nc*3);
  nablaphi.resize(_nc*6);


  type Jac[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
  type JacI[3][3];

  const double *dxi=_dphidxi[ig];
  const double *deta=_dphideta[ig];
  const double *dzeta=_dphidzeta[ig];

  for (int inode=0; inode<_nc; inode++,dxi++,deta++,dzeta++) {
    Jac[0][0]+=(*dxi)*vt[0][inode];
    Jac[0][1]+=(*dxi)*vt[1][inode];
    Jac[0][2]+=(*dxi)*vt[2][inode];
    Jac[1][0]+=(*deta)*vt[0][inode];
    Jac[1][1]+=(*deta)*vt[1][inode];
    Jac[1][2]+=(*deta)*vt[2][inode];
    Jac[2][0]+=(*dzeta)*vt[0][inode];
    Jac[2][1]+=(*dzeta)*vt[1][inode];
    Jac[2][2]+=(*dzeta)*vt[2][inode];
  }
  type det=(Jac[0][0]*(Jac[1][1]*Jac[2][2]-Jac[1][2]*Jac[2][1])+
		      Jac[0][1]*(Jac[1][2]*Jac[2][0]-Jac[1][0]*Jac[2][2])+
		      Jac[0][2]*(Jac[1][0]*Jac[2][1]-Jac[1][1]*Jac[2][0]));

  JacI[0][0]= (-Jac[1][2]*Jac[2][1] + Jac[1][1]*Jac[2][2])/det;
  JacI[0][1]= ( Jac[0][2]*Jac[2][1] - Jac[0][1]*Jac[2][2])/det;
  JacI[0][2]= (-Jac[0][2]*Jac[1][1] + Jac[0][1]*Jac[1][2])/det;
  JacI[1][0]= ( Jac[1][2]*Jac[2][0] - Jac[1][0]*Jac[2][2])/det;
  JacI[1][1]= (-Jac[0][2]*Jac[2][0] + Jac[0][0]*Jac[2][2])/det;
  JacI[1][2]= ( Jac[0][2]*Jac[1][0] - Jac[0][0]*Jac[1][2])/det;
  JacI[2][0]= (-Jac[1][1]*Jac[2][0] + Jac[1][0]*Jac[2][1])/det;
  JacI[2][1]= ( Jac[0][1]*Jac[2][0] - Jac[0][0]*Jac[2][1])/det;
  JacI[2][2]= (-Jac[0][1]*Jac[1][0] + Jac[0][0]*Jac[1][1])/det;

  Weight=det*_gauss.GetGaussWeightsPointer()[ig];

  dxi=_dphidxi[ig];
  deta=_dphideta[ig];
  dzeta=_dphidzeta[ig];

  const double *dxi2=_d2phidxi2[ig];
  const double *deta2=_d2phideta2[ig];
  const double *dzeta2=_d2phidzeta2[ig];
  const double *dxideta=_d2phidxideta[ig];
  const double *detadzeta=_d2phidetadzeta[ig];
  const double *dzetadxi=_d2phidzetadxi[ig];

  for (int inode=0; inode<_nc; inode++, dxi++,deta++,dzeta++,dxi2++,deta2++,dzeta2++,dxideta++,detadzeta++,dzetadxi++) {

    phi[inode]=_phi[ig][inode];

    gradphi[3*inode+0]=(*dxi)*JacI[0][0] + (*deta)*JacI[0][1] + (*dzeta)*JacI[0][2];
    gradphi[3*inode+1]=(*dxi)*JacI[1][0] + (*deta)*JacI[1][1] + (*dzeta)*JacI[1][2];
    gradphi[3*inode+2]=(*dxi)*JacI[2][0] + (*deta)*JacI[2][1] + (*dzeta)*JacI[2][2];

    nablaphi[6*inode+0]=
      ( (*dxi2)    *JacI[0][0] + (*dxideta)  *JacI[0][1] + (*dzetadxi) *JacI[0][2] )*JacI[0][0]+
      ( (*dxideta) *JacI[0][0] + (*deta2)    *JacI[0][1] + (*detadzeta)*JacI[0][2] )*JacI[0][1]+
      ( (*dzetadxi)*JacI[0][0] + (*detadzeta)*JacI[0][1] + (*dzeta2)   *JacI[0][2] )*JacI[0][2];
    nablaphi[6*inode+1]=
      ( (*dxi2)    *JacI[1][0] + (*dxideta)  *JacI[1][1] + (*dzetadxi) *JacI[1][2] )*JacI[1][0]+
      ( (*dxideta) *JacI[1][0] + (*deta2)    *JacI[1][1] + (*detadzeta)*JacI[1][2] )*JacI[1][1]+
      ( (*dzetadxi)*JacI[1][0] + (*detadzeta)*JacI[1][1] + (*dzeta2)   *JacI[1][2] )*JacI[1][2];
    nablaphi[6*inode+2]=
      ( (*dxi2)    *JacI[2][0] + (*dxideta)  *JacI[2][1] + (*dzetadxi) *JacI[2][2] )*JacI[2][0]+
      ( (*dxideta) *JacI[2][0] + (*deta2)    *JacI[2][1] + (*detadzeta)*JacI[2][2] )*JacI[2][1]+
      ( (*dzetadxi)*JacI[2][0] + (*detadzeta)*JacI[2][1] + (*dzeta2)   *JacI[2][2] )*JacI[2][2];
    nablaphi[6*inode+3]=
      ( (*dxi2)    *JacI[0][0] + (*dxideta)  *JacI[0][1] + (*dzetadxi) *JacI[0][2] )*JacI[1][0]+
      ( (*dxideta) *JacI[0][0] + (*deta2)    *JacI[0][1] + (*detadzeta)*JacI[0][2] )*JacI[1][1]+
      ( (*dzetadxi)*JacI[0][0] + (*detadzeta)*JacI[0][1] + (*dzeta2)   *JacI[0][2] )*JacI[1][2];
    nablaphi[6*inode+4]=
      ( (*dxi2)    *JacI[1][0] + (*dxideta)  *JacI[1][1] + (*dzetadxi) *JacI[1][2] )*JacI[2][0]+
      ( (*dxideta) *JacI[1][0] + (*deta2)    *JacI[1][1] + (*detadzeta)*JacI[1][2] )*JacI[2][1]+
      ( (*dzetadxi)*JacI[1][0] + (*detadzeta)*JacI[1][1] + (*dzeta2)   *JacI[1][2] )*JacI[2][2];
    nablaphi[6*inode+5]=
      ( (*dxi2)    *JacI[2][0] + (*dxideta)  *JacI[2][1] + (*dzetadxi) *JacI[2][2] )*JacI[0][0]+
      ( (*dxideta) *JacI[2][0] + (*deta2)    *JacI[2][1] + (*detadzeta)*JacI[2][2] )*JacI[0][1]+
      ( (*dzetadxi)*JacI[2][0] + (*detadzeta)*JacI[2][1] + (*dzeta2)   *JacI[2][2] )*JacI[0][2];
  }

}

//---------------------------------------------------------------------------------------------------------

} //end namespace femus


