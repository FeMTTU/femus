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
#include "Elem.hpp"
#include "NumericVector.hpp"

using std::cout;
using std::endl;


namespace femus {

unsigned elem_type::_refindex=1;

elem_type::~elem_type() {
  delete [] X;
  delete [] KVERT_IND;
  delete [] IND;

  delete [] prol_val;
  delete [] prol_ind;
  delete [] mem_prol_val;
  delete [] mem_prol_ind;
  
  delete pt_basis;
  
};

//----------------------------------------------------------------------------------------------------
// build matrix sparsity pattern size and build prolungator matrix for the LsysPde  Matrix
//-----------------------------------------------------------------------------------------------------

void elem_type::BuildProlongation(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc, SparseMatrix* Projmat, 
				  const unsigned &index_sol, const unsigned &kkindex_sol) const {
  vector<int> cols(27);
   
  for (int i=0; i<nf_; i++) {
    int i0=KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=lspdec._msh->el->GetChildElement(ielc,i0);
    int i1=KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=lspdef._msh->el->GetDof(ielf,i1,type_);
    int irow=lspdef.GetKKDof(index_sol,kkindex_sol,iadd);  //  local-id to dof 
    int ncols=prol_ind[i+1]-prol_ind[i];
    cols.assign(ncols,0);
    for (int k=0; k<ncols; k++) {
      int j=prol_ind[i][k]; 
      int jadd=lspdec._msh->el->GetDof(ielc,j,type_);
      int jj=lspdec.GetKKDof(index_sol,kkindex_sol,jadd); 
      cols[k]=jj;
    }
    Projmat->insert_row(irow,ncols,cols,prol_val[i]);
  }
}


void elem_type::GetSparsityPatternSize(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc,  
				       NumericVector* NNZ_d, NumericVector* NNZ_o,
				       const unsigned &index_sol, const unsigned &kkindex_sol) const {
     
  for (int i=0; i<nf_; i++) {
    int i0=KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=lspdec._msh->el->GetChildElement(ielc,i0);
    int i1=KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=lspdef._msh->el->GetDof(ielf,i1,type_);
    int irow=lspdef.GetKKDof(index_sol,kkindex_sol,iadd);  //  local-id to dof 
    
    int iproc=0;
    while (irow < lspdef.KKoffset[0][iproc] || irow >= lspdef.KKoffset[lspdef.KKIndex.size()-1][iproc] ) iproc++;
    int ncols=prol_ind[i+1]-prol_ind[i];
    
    int counter_o=0;
    for (int k=0; k<ncols; k++) {
      int j=prol_ind[i][k]; 
      int jadd=lspdec._msh->el->GetDof(ielc,j,type_);
      int jcolumn=lspdec.GetKKDof(index_sol,kkindex_sol,jadd); 
      if(jcolumn < lspdec.KKoffset[0][iproc] || jcolumn >= lspdec.KKoffset[lspdef.KKIndex.size()-1][iproc] ) counter_o++;
    }
       
    NNZ_d->set(irow,ncols-counter_o);
    NNZ_o->set(irow,counter_o);
    
  }
}

void elem_type::BuildRestrictionTranspose(const LinearEquation &lspdef,const LinearEquation &lspdec, const int& ielc, SparseMatrix* Projmat, 
					  const unsigned &index_sol, const unsigned &kkindex_sol, const bool &TestDisp) const {
  vector<int> cols(27);
  bool fluid_region = (2==lspdec._msh->el->GetElementMaterial(ielc))?1:0;
  
  vector <double> copy_prol_val;
  copy_prol_val.reserve(27); 
  for (int i=0; i<nf_; i++) {
    int i0=KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=lspdec._msh->el->GetChildElement(ielc,i0);
    int i1=KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=lspdef._msh->el->GetDof(ielf,i1,type_);
        
    int irow=lspdef.GetKKDof(index_sol,kkindex_sol,iadd);  //  local-id to dof 
    int ncols=prol_ind[i+1]-prol_ind[i];
    
    bool isolidmark=lspdef._msh->el->GetNodeRegion(iadd);
    
    cols.assign(ncols,0);
    copy_prol_val.resize(ncols);
    for (int k=0; k<ncols; k++) {
      int j=prol_ind[i][k]; 
      int jadd=lspdec._msh->el->GetDof(ielc,j,type_);
      int jcolumn=lspdec.GetKKDof(index_sol,kkindex_sol,jadd); 
      cols[k]=jcolumn;
      
      bool jsolidmark=lspdef._msh->el->GetNodeRegion(jadd); 
      
      copy_prol_val[k]=(!TestDisp || !fluid_region || isolidmark==jsolidmark)?prol_val[i][k]:0.;
    }    
    Projmat->insert_row(irow,ncols,cols,&copy_prol_val[0]);
  }
}



//----------------------------------------------------------------------------------------------------
// build matrix sparsity pattern size and build prolungator matrix for single solution
//-----------------------------------------------------------------------------------------------------

void elem_type::GetSparsityPatternSize(const mesh &meshf,const mesh &meshc, const int& ielc, NumericVector* NNZ_d, NumericVector* NNZ_o) const {
   
  for (int i=0; i<nf_; i++) {
    int i0=KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=meshc.el->GetChildElement(ielc,i0);
    int i1=KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=meshf.el->GetDof(ielf,i1,type_);
    int irow=meshf.GetMetisDof(iadd,_SolType);  //  local-id to dof
    int iproc=0;
    while (irow < meshf.MetisOffset[_SolType][iproc] || irow >= meshf.MetisOffset[_SolType][iproc+1] ) iproc++;
    
    int ncols=prol_ind[i+1]-prol_ind[i];
    unsigned counter_o=0;
    for (int k=0; k<ncols; k++) {
      int j=prol_ind[i][k]; 
      int jadd=meshc.el->GetDof(ielc,j,type_);
      int jcolumn=meshc.GetMetisDof(jadd,_SolType);
      if(jcolumn < meshc.MetisOffset[_SolType][iproc] || jcolumn >= meshc.MetisOffset[_SolType][iproc+1] ) counter_o++;
    }
       
    NNZ_d->set(irow,ncols-counter_o);
    NNZ_o->set(irow,counter_o);
    
  }
}

void elem_type::BuildProlongation(const mesh &meshf,const mesh &meshc, const int& ielc,
				  SparseMatrix* Projmat) const {
  vector<int> cols(27);
  for (int i=0; i<nf_; i++) {
    int i0=KVERT_IND[i][0]; //id of the subdivision of the fine element
    int ielf=meshc.el->GetChildElement(ielc,i0);
    int i1=KVERT_IND[i][1]; //local id node on the subdivision of the fine element
    int iadd=meshf.el->GetDof(ielf,i1,type_);
    int irow=meshf.GetMetisDof(iadd,_SolType);  //  local-id to dof 
    int ncols=prol_ind[i+1]-prol_ind[i];
    cols.assign(ncols,0);
    for (int k=0; k<ncols; k++) {
      int j=prol_ind[i][k]; 
      int jadd=meshc.el->GetDof(ielc,j,type_);
      int jcolumn=meshc.GetMetisDof(jadd,_SolType);
      cols[k]=jcolumn;
    }

    Projmat->insert_row(irow,ncols,cols,prol_val[i]);
  }
}

//----------------------------------------------------------------------------------------------------
// prolungator for solution printing
//----------------------------------------------------------------------------------------------------

void elem_type::GetSparsityPatternSize(const mesh& mesh,const int& iel, NumericVector* NNZ_d, NumericVector* NNZ_o, const unsigned &itype) const{
  for (int i=0; i<ncf_[itype]; i++) {
    int inode=mesh.el->GetDof(iel,i,type_);
    int irow=mesh.GetMetisDof(inode,itype);
    int iproc=0;
    while (irow < mesh.MetisOffset[itype][iproc] || irow >= mesh.MetisOffset[itype][iproc+1] ) iproc++;
    int ncols=prol_ind[i+1]-prol_ind[i];
    unsigned counter_o=0;
    for (int k=0; k<ncols; k++) {
      int jj=prol_ind[i][k];
      int jnode   = mesh.el->GetDof(iel,jj,type_);
      int jcolumn = mesh.GetMetisDof(jnode,_SolType);
      if(jcolumn < mesh.MetisOffset[_SolType][iproc] || jcolumn >= mesh.MetisOffset[_SolType][iproc+1] ) counter_o++;
    }
    NNZ_d->set(irow,ncols-counter_o);
    NNZ_o->set(irow,counter_o);
  }
}


void elem_type::BuildProlongation(const mesh& mesh,const int& iel, SparseMatrix* Projmat,const unsigned &itype) const{
  vector<int> cols(27);
  for (int i=0; i<ncf_[itype]; i++) {
    int inode=mesh.el->GetDof(iel,i,type_);
    int irow=mesh.GetMetisDof(inode,itype);
    int ncols=prol_ind[i+1]-prol_ind[i];
    cols.assign(ncols,0);
    for (int k=0; k<ncols; k++) {
      int jj=prol_ind[i][k];
      int jnode=mesh.el->GetDof(iel,jj,type_);
      int jcolumn=mesh.GetMetisDof(jnode,_SolType);
      cols[k]=jcolumn;
    }
    Projmat->insert_row(irow,ncols,cols,prol_val[i]);
  }
}

elem_type_1D::elem_type_1D(const char *solid, const char *order, const char *order_gauss) :
	      elem_type(){

  //************ BEGIN GAUSS SETUP ******************	
  Gauss gauss(solid, order_gauss);
  GaussWeight = gauss.GaussWeight;
  GaussPoints = gauss.GaussPoints;
  //************ END GAUSS SETUP ******************
	
  //************ BEGIN FE and MG SETUP ******************	
  if (!strcmp(order,"linear")) 		 _SolType=0;   
  else if (!strcmp(order,"quadratic")) 	 _SolType=1;   
  else if (!strcmp(order,"biquadratic")) _SolType=2;   
  else if (!strcmp(order,"constant"))    _SolType=3;   
  else if (!strcmp(order,"disc_linear")) _SolType=4;   
   
  if (!strcmp(solid,"line")) { //line
    if 	    (_SolType == 0) pt_basis = new line1;
    else if (_SolType == 1) pt_basis = new line2;
    else if (_SolType == 2) pt_basis = new line2;
    else if (_SolType == 3) pt_basis = new line0;
    else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    }
  } 
  else {
    cout<<solid<<" is not a valid option"<<endl;
    exit(0);
  }
  
  // get data from basis object
  type_	  = pt_basis->_type; 
  nc_ 	  = pt_basis->_nc;
  nf_ 	  = pt_basis->_nf;
  ncf_[0] = pt_basis->_ncf0;
  ncf_[1] = pt_basis->_ncf1;
  ncf_[2] = pt_basis->_ncf2;
  
  
  IND=new const int * [nc_];
  for (int i=0; i<nc_; i++){
    IND[i]=pt_basis->getIND(i);
  }
  KVERT_IND=new const int * [nf_];
  X=new const double * [nf_];
  for (int i=0; i<nf_; i++) {
    KVERT_IND[i]=pt_basis->getKVERT_IND(i);
    X[i]=pt_basis->getX(i);
  }

  // local projection matrix evaluation
  int counter=0;
  for (int i=0; i<nf_; i++) {
    for (int j=0; j<nc_; j++) {
      double phi=pt_basis->eval_phi(IND[j],X[i]);
      if ( fabs(phi) >= 1.0e-14 ){
        counter++;
      }
    }
  }
  double *pt_d;
  int *pt_i;

  prol_val=new double * [nf_+1];
  prol_ind=new int * [nf_+1];
  mem_prol_val=new double [counter];
  mem_prol_ind=new int [counter];

  pt_d=mem_prol_val;
  pt_i= mem_prol_ind;
  for (int i=0; i<nf_; i++) {
    prol_val[i]=pt_d;
    prol_ind[i]=pt_i;
    for (int j=0; j<nc_; j++) {
      double phi=pt_basis->eval_phi(IND[j],X[i]);
      if ( fabs(phi) >= 1.0e-14 ){
        *(pt_d++)=phi;
        *(pt_i++)=j;
      }
    }
  }
  
  prol_val[nf_]=pt_d;
  prol_ind[nf_]=pt_i;
    
  // shape function and its derivatives evaluated at Gauss'points
  phi= new double*[GaussPoints];
  dphidxi  = new double*[GaussPoints];
  d2phidxi2  = new double*[GaussPoints];

  phi_memory=new double [GaussPoints*nc_];
  dphidxi_memory  =new double [GaussPoints*nc_];
  d2phidxi2_memory  =new double [GaussPoints*nc_];

  for (unsigned i=0; i<GaussPoints; i++) {
    phi[i]=&phi_memory[i*nc_];
    dphidxi[i]  =&dphidxi_memory[i*nc_];
    d2phidxi2[i]  =&d2phidxi2_memory[i*nc_];
  }

 const double *ptx[1]={GaussWeight+GaussPoints};
  for (unsigned i=0; i<GaussPoints; i++){
    double x[1];
    for (unsigned j=0; j<1;j++) {
      x[j] = *ptx[j];
      ptx[j]++;
    }  
    
    for (int j=0; j<nc_; j++) {
      phi[i][j] = pt_basis->eval_phi(IND[j],x); 	
      dphidxi[i][j] = pt_basis->eval_dphidx(IND[j],x); 		
      d2phidxi2[i][j] = pt_basis->eval_d2phidx2(IND[j],x);
    }
  }
}


elem_type_2D::elem_type_2D(const char *solid, const char *order, const char *order_gauss):
	      elem_type(){

  //************ BEGIN GAUSS SETUP ******************	
  Gauss gauss(solid,order_gauss);
  GaussWeight = gauss.GaussWeight;
  GaussPoints = gauss.GaussPoints;
  //************ END GAUSS SETUP ******************
	
  //************ BEGIN FE and MG SETUP ******************	
  if 	  (!strcmp(order,"linear")) 	 _SolType=0;   
  else if (!strcmp(order,"quadratic")) 	 _SolType=1;   
  else if (!strcmp(order,"biquadratic")) _SolType=2;   
  else if (!strcmp(order,"constant"))    _SolType=3;   
  else if (!strcmp(order,"disc_linear")) _SolType=4;   
  
  if (!strcmp(solid,"quad")) { //QUAD
    if 	    (_SolType == 0) pt_basis = new quad1;
    else if (_SolType == 1) pt_basis = new quadth;
    else if (_SolType == 2) pt_basis = new quad2;
    else if (_SolType == 3) pt_basis = new quad0;
    else if (_SolType == 4) pt_basis = new quadpwl;
    else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    }
  }      
  else if (!strcmp(solid,"tri")) { //TRIANGLE
    
    if 	    (_SolType == 0) pt_basis = new tri1;
    else if (_SolType == 1) pt_basis = new tri2;
    else if (_SolType == 2) pt_basis = new tri2;
    else if (_SolType == 3) pt_basis = new tri0;
    else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    } 
  } 
  else {
    cout<<solid<<" is not a valid option"<<endl;
    exit(0);
  }
  
  // get data from basis object
  type_	  = pt_basis->_type; 
  nc_ 	  = pt_basis->_nc;
  nf_ 	  = pt_basis->_nf;
  ncf_[0] = pt_basis->_ncf0;
  ncf_[1] = pt_basis->_ncf1;
  ncf_[2] = pt_basis->_ncf2;

  IND=new const int * [nc_];
  for (int i=0; i<nc_; i++){
    IND[i]=pt_basis->getIND(i);
  }
  KVERT_IND=new const int * [nf_];
  X=new const double * [nf_];
  for (int i=0; i<nf_; i++) {
    KVERT_IND[i]=pt_basis->getKVERT_IND(i);
    X[i]=pt_basis->getX(i);
  }
  
  // local projection matrix evaluation   
  int counter=0;
  for (int i=0; i<nf_; i++) {
    for (int j=0; j<nc_; j++) {
      double phi=pt_basis->eval_phi(IND[j],X[i]);
      if (type_==16) {
        if (i/4==1) phi=pt_basis->eval_dphidx(IND[j],X[i]);
        else if (i/4==2) phi=pt_basis->eval_dphidy(IND[j],X[i]);
      }
      if ( fabs(phi) >= 1.0e-14 ){
        counter++;
      }
    }
  }
  double *pt_d;
  int *pt_i;

  prol_val=new double * [nf_+1];
  prol_ind=new int * [nf_+1];
  mem_prol_val=new double [counter];
  mem_prol_ind=new int [counter];

  pt_d=mem_prol_val;
  pt_i= mem_prol_ind;
  for (int i=0; i<nf_; i++) {
    prol_val[i]=pt_d;
    prol_ind[i]=pt_i;
    for (int j=0; j<nc_; j++) {
      double phi=pt_basis->eval_phi(IND[j],X[i]);
      if (type_==16) {
        if (i/4==1)
          phi=pt_basis->eval_dphidx(IND[j],X[i])/2.;
        else if (i/4==2)
          phi=pt_basis->eval_dphidy(IND[j],X[i])/2.;
      } 
      if ( fabs(phi) >= 1.0e-14 ){
        *(pt_d++)=phi;
        *(pt_i++)=j;
      }
    }
  }
  
  prol_val[nf_]=pt_d;
  prol_ind[nf_]=pt_i;

  // shape function and its derivatives evaluated at Gauss'points
  phi= new double*[GaussPoints];
  dphidxi  = new double*[GaussPoints];
  dphideta = new double*[GaussPoints];
  
  d2phidxi2  = new double*[GaussPoints];
  d2phideta2 = new double*[GaussPoints];
  
  d2phidxideta  = new double*[GaussPoints];

  phi_memory=new double [GaussPoints*nc_];
  dphidxi_memory  =new double [GaussPoints*nc_];
  dphideta_memory =new double [GaussPoints*nc_];
  
  d2phidxi2_memory  =new double [GaussPoints*nc_];
  d2phideta2_memory =new double [GaussPoints*nc_];
  
  d2phidxideta_memory  =new double [GaussPoints*nc_];

  for (unsigned i=0; i<GaussPoints; i++) {
    phi[i]=&phi_memory[i*nc_];
    dphidxi[i]  =&dphidxi_memory[i*nc_];
    dphideta[i] =&dphideta_memory[i*nc_];
    
    d2phidxi2[i]  =&d2phidxi2_memory[i*nc_];
    d2phideta2[i] =&d2phideta2_memory[i*nc_];
    
    d2phidxideta[i]  = &d2phidxideta_memory[i*nc_];
    
  }

  const double *ptx[2]={GaussWeight+GaussPoints, GaussWeight+2*GaussPoints};
  for (unsigned i=0; i<GaussPoints; i++){
    double x[2];
    for (unsigned j=0; j<2;j++) {
      x[j] = *ptx[j];
      ptx[j]++;
    }
    for (int j=0; j<nc_; j++) {
      phi[i][j] = pt_basis->eval_phi(IND[j],x); 	
      dphidxi[i][j] = pt_basis->eval_dphidx(IND[j],x); 	
      dphideta[i][j] = pt_basis->eval_dphidy(IND[j],x); 	
      d2phidxi2[i][j] = pt_basis->eval_d2phidx2(IND[j],x);
      d2phideta2[i][j] = pt_basis->eval_d2phidy2(IND[j],x);
      d2phidxideta[i][j] = pt_basis->eval_d2phidxdy(IND[j],x);
    }
  }
}

elem_type_3D::elem_type_3D(const char *solid, const char *order, const char *order_gauss) :
	      elem_type(){

  //************ BEGIN GAUSS SETUP ******************	
  Gauss gauss(solid,order_gauss);
  GaussWeight = gauss.GaussWeight;
  GaussPoints = gauss.GaussPoints;
  //************ END GAUSS SETUP ******************
	
  //************ BEGIN FE and MG SETUP ******************	
  if 	  (!strcmp(order,"linear")) 	 _SolType=0;   
  else if (!strcmp(order,"quadratic")) 	 _SolType=1;   
  else if (!strcmp(order,"biquadratic")) _SolType=2;   
  else if (!strcmp(order,"constant"))    _SolType=3;   
  else if (!strcmp(order,"disc_linear")) _SolType=4;  
  
  if (!strcmp(solid,"hex")) {//HEX
    
    if 	    (_SolType == 0) pt_basis = new hex1;
    else if (_SolType == 1) pt_basis = new hexth;
    else if (_SolType == 2) pt_basis = new hex2;
    else if (_SolType == 3) pt_basis = new hex0;
    else if (_SolType == 4) pt_basis = new hexpwl;
    else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    }
  }     
  else if (!strcmp(solid,"wedge")) { //WEDGE
    if 	    (_SolType == 0) pt_basis = new wedge1;
    else if (_SolType == 1) pt_basis = new wedgeth;
    else if (_SolType == 2) pt_basis = new wedge2;
    else if (_SolType == 3) pt_basis = new wedge0;     
    else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    }
  } 
  else if (!strcmp(solid,"tet")) { //TETRAHEDRA
    if 	    (_SolType == 0) pt_basis = new tet1;
    else if (_SolType == 1) pt_basis = new tet2;
    else if (_SolType == 2) pt_basis = new tet2;
    else if (_SolType == 3) pt_basis = new tet0;     
    else {
      cout<<order<<" is not a valid option for "<<solid<<endl;
      exit(0);
    }
  } 
  else {
    cout<<solid<<" is not a valid option"<<endl;
    exit(0);
  }

  // get data from basis object
  type_	  = pt_basis->_type; 
  nc_ 	  = pt_basis->_nc;
  nf_ 	  = pt_basis->_nf;
  ncf_[0] = pt_basis->_ncf0;
  ncf_[1] = pt_basis->_ncf1;
  ncf_[2] = pt_basis->_ncf2;
    
  IND=new const int * [nc_];
  for (int i=0; i<nc_; i++){
    IND[i]=pt_basis->getIND(i);
  }
  KVERT_IND=new const int * [nf_];
  X=new const double * [nf_];
  for (int i=0; i<nf_; i++) {
    KVERT_IND[i]=pt_basis->getKVERT_IND(i);
    X[i]=pt_basis->getX(i);
  }
    
  // local projection matrix evaluation
  int counter=0;
  for (int i=0; i<nf_; i++) {
    for (int j=0; j<nc_; j++) {
      double phi=pt_basis->eval_phi(IND[j],X[i]);
      if ( type_ == 18 ) {
        if ( i/8 == 1 )
          phi=pt_basis->eval_dphidx(IND[j],X[i])/2.;
        else if ( i/8 == 2 )
          phi=pt_basis->eval_dphidy(IND[j],X[i])/2.;
        else if ( i/8 == 3 )
          phi=pt_basis->eval_dphidz(IND[j],X[i])/2.;
      }
      if ( fabs(phi) >= 1.0e-14 ){
        counter++;
      }
    }
  }
  double *pt_d;
  int *pt_i;

  prol_val=new double * [nf_+1];
  prol_ind=new int * [nf_+1];
  mem_prol_val=new double [counter];
  mem_prol_ind=new int [counter];

  pt_d=mem_prol_val;
  pt_i= mem_prol_ind;
  for (int i=0; i<nf_; i++) {
    prol_val[i] = pt_d;
    prol_ind[i] = pt_i;
    for (int j=0; j<nc_; j++) {
      double phi=pt_basis->eval_phi(IND[j],X[i]);
      if ( type_== 18 ) {
        if ( i/8 == 1 )
          phi = pt_basis->eval_dphidx(IND[j],X[i])/2.;
        else if ( i/8 == 2 )
          phi = pt_basis->eval_dphidy(IND[j],X[i])/2.;
        else if ( i/8 == 3 )
          phi = pt_basis->eval_dphidz(IND[j],X[i])/2.;
      }
      if ( fabs(phi) >= 1.0e-14 ){
        *(pt_d++)=phi;
        *(pt_i++)=j;
      }
    }
  }
  
  prol_val[nf_]=pt_d;
  prol_ind[nf_]=pt_i;

  // shape function and its derivatives evaluated at Gauss'points

  phi= new double*[GaussPoints];
  dphidxi  = new double*[GaussPoints];
  dphideta = new double*[GaussPoints];
  dphidzeta= new double*[GaussPoints];
  
  d2phidxi2  = new double*[GaussPoints];
  d2phideta2 = new double*[GaussPoints];
  d2phidzeta2= new double*[GaussPoints];
  
  d2phidxideta  = new double*[GaussPoints];
  d2phidetadzeta = new double*[GaussPoints];
  d2phidzetadxi= new double*[GaussPoints];

  phi_memory=new double [GaussPoints*nc_];
  dphidxi_memory  =new double [GaussPoints*nc_];
  dphideta_memory =new double [GaussPoints*nc_];
  dphidzeta_memory=new double [GaussPoints*nc_];
  
  d2phidxi2_memory  =new double [GaussPoints*nc_];
  d2phideta2_memory =new double [GaussPoints*nc_];
  d2phidzeta2_memory=new double [GaussPoints*nc_];
  
  d2phidxideta_memory  =new double [GaussPoints*nc_];
  d2phidetadzeta_memory =new double [GaussPoints*nc_];
  d2phidzetadxi_memory=new double [GaussPoints*nc_];

  for (unsigned i=0; i<GaussPoints; i++) {
    phi[i]=&phi_memory[i*nc_];
    dphidxi[i]  =&dphidxi_memory[i*nc_];
    dphideta[i] =&dphideta_memory[i*nc_];
    dphidzeta[i]=&dphidzeta_memory[i*nc_];
    
    d2phidxi2[i]  =&d2phidxi2_memory[i*nc_];
    d2phideta2[i] =&d2phideta2_memory[i*nc_];
    d2phidzeta2[i]=&d2phidzeta2_memory[i*nc_];
    
    d2phidxideta[i]  = &d2phidxideta_memory[i*nc_];
    d2phidetadzeta[i]= &d2phidetadzeta_memory[i*nc_];
    d2phidzetadxi[i] = &d2phidzetadxi_memory[i*nc_];
    
  }

  const double *ptx[3]={GaussWeight+GaussPoints, GaussWeight+2*GaussPoints, GaussWeight+3*GaussPoints};
  for (unsigned i=0; i<GaussPoints; i++){
    double x[3];
    for (unsigned j=0; j<3;j++) {
      x[j] = *ptx[j];
      ptx[j]++;
    }  
    
    for (int j=0; j<nc_; j++) {
      phi[i][j] = pt_basis->eval_phi(IND[j],x); 	
      dphidxi[i][j] = pt_basis->eval_dphidx(IND[j],x); 	
      dphideta[i][j] = pt_basis->eval_dphidy(IND[j],x); 	
      dphidzeta[i][j] = pt_basis->eval_dphidz(IND[j],x); 	
      d2phidxi2[i][j] = pt_basis->eval_d2phidx2(IND[j],x);
      d2phideta2[i][j] = pt_basis->eval_d2phidy2(IND[j],x);
      d2phidzeta2[i][j] = pt_basis->eval_d2phidz2(IND[j],x);
      d2phidxideta[i][j] = pt_basis->eval_d2phidxdy(IND[j],x);
      d2phidetadzeta[i][j] = pt_basis->eval_d2phidydz(IND[j],x);
      d2phidzetadxi[i][j] = pt_basis->eval_d2phidzdx(IND[j],x);
    }
  }
}

//---------------------------------------------------------------------------------------------------------

void elem_type_1D::Jacobian_AD(const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight, 
			      vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const{
  other_phi.resize(nc_);
  gradphi.resize(nc_*1);
  nablaphi.resize(nc_*1);				
  
  adept::adouble Jac=0.;
  adept::adouble JacI;
  
  const double *dxi=dphidxi[ig];
  
  for (int inode=0; inode<nc_; inode++,dxi++) {
    Jac+=(*dxi)*vt[0][inode];
  }
    
  Weight=Jac*GaussWeight[ig];

  JacI=1/Jac;
  
  dxi = dphidxi[ig];
  const double *dxi2 = d2phidxi2[ig];
  
  for (int inode=0; inode<nc_; inode++,dxi++, dxi2++) {
    other_phi[inode]=phi[ig][inode];
    gradphi[inode]=(*dxi)*JacI;
    nablaphi[inode] = (*dxi2)*JacI*JacI;
  }
  
}

void elem_type_1D::Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight, 
    			    vector < double > &other_phi, vector < double > &gradphi, vector < double > &nablaphi) const{
  other_phi.resize(nc_);
  gradphi.resize(nc_*1);
  nablaphi.resize(nc_*1);				
  
  double Jac=0.;
  double JacI;
  
  const double *dxi=dphidxi[ig];
  
  for (int inode=0; inode<nc_; inode++,dxi++) {
    Jac+=(*dxi)*vt[0][inode];
  }
    
  Weight=Jac*GaussWeight[ig];

  JacI=1/Jac;
  
  dxi = dphidxi[ig];
  const double *dxi2 = d2phidxi2[ig];
  
  for (int inode=0; inode<nc_; inode++,dxi++, dxi2++) {
    other_phi[inode]=phi[ig][inode];
    gradphi[inode]=(*dxi)*JacI;
    nablaphi[inode] = (*dxi2)*JacI*JacI;
  }
  
}

void elem_type_1D::JacobianSur_AD(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight, 
				  vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const {

  other_phi.resize(nc_);
  normal.resize(2);				
			
  adept::adouble Jac[2][2]={{0.,0.},{0.,0.}};
  adept::adouble JacI[2][2];

  const double *dfeta=dphidxi[ig];

  for (int inode=0; inode<nc_; inode++,dfeta++) {
    Jac[0][0] += (*dfeta)*vt[0][inode];
    Jac[1][0] += (*dfeta)*vt[1][inode];
  }

//   normal module
  adept::adouble modn = sqrt(Jac[0][0]*Jac[0][0] + Jac[1][0]*Jac[1][0]);

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
  adept::adouble det= (Jac[0][0]*Jac[1][1]-Jac[0][1]*Jac[1][0]);

  JacI[0][0] =  Jac[1][1]/det;
  JacI[0][1] = -Jac[0][1]/det;
  JacI[1][0] = -Jac[1][0]/det;
  JacI[1][1] =  Jac[0][0]/det;

  Weight = det*GaussWeight[ig];
  
  for(int inode=0;inode<nc_;inode++){
    other_phi[inode]=phi[ig][inode];
  }

}

void elem_type_1D::JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight, 
			       vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const {

  other_phi.resize(nc_);
  normal.resize(2);				
			
  double Jac[2][2]={{0.,0.},{0.,0.}};
  double JacI[2][2];

  const double *dfeta=dphidxi[ig];

  for (int inode=0; inode<nc_; inode++,dfeta++) {
    Jac[0][0] += (*dfeta)*vt[0][inode];
    Jac[1][0] += (*dfeta)*vt[1][inode];
  }

  //   normal module
  double modn = sqrt(Jac[0][0]*Jac[0][0] + Jac[1][0]*Jac[1][0]);

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
  double det= (Jac[0][0]*Jac[1][1]-Jac[0][1]*Jac[1][0]);

  JacI[0][0] =  Jac[1][1]/det;
  JacI[0][1] = -Jac[0][1]/det;
  JacI[1][0] = -Jac[1][0]/det;
  JacI[1][1] =  Jac[0][0]/det;

  Weight = det*GaussWeight[ig];
  
  for(int inode=0;inode<nc_;inode++){
    other_phi[inode]=phi[ig][inode];
  }

}

//---------------------------------------------------------------------------------------------------------
void elem_type_2D::Jacobian_AD(const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight, 
			      vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const{
	
  other_phi.resize(nc_);
  gradphi.resize(nc_*2);
  nablaphi.resize(nc_*3);
				
  adept::adouble Jac[2][2]={{0,0},{0,0}};
  adept::adouble JacI[2][2];
  const double *dxi=dphidxi[ig];
  const double *deta=dphideta[ig];
  for (int inode=0; inode<nc_; inode++,dxi++,deta++){
    Jac[0][0]+=(*dxi)*vt[0][inode];
    Jac[0][1]+=(*dxi)*vt[1][inode];
    Jac[1][0]+=(*deta)*vt[0][inode];
    Jac[1][1]+=(*deta)*vt[1][inode];
  }
  adept::adouble det=(Jac[0][0]*Jac[1][1]-Jac[0][1]*Jac[1][0]);

  JacI[0][0]= Jac[1][1]/det;
  JacI[0][1]=-Jac[0][1]/det;
  JacI[1][0]=-Jac[1][0]/det;
  JacI[1][1]= Jac[0][0]/det;

  Weight=det*GaussWeight[ig];
     
  dxi=dphidxi[ig];
  deta=dphideta[ig];
  
  const double *dxi2=d2phidxi2[ig];
  const double *deta2=d2phideta2[ig];
  const double *dxideta=d2phidxideta[ig];
    
  for (int inode=0; inode<nc_; inode++, dxi++, deta++, dxi2++, deta2++, dxideta++) {
    
    other_phi[inode]=phi[ig][inode];
     
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

void elem_type_2D::Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight, 
			    vector < double > &other_phi, vector < double > &gradphi, vector < double > &nablaphi) const{
	
  other_phi.resize(nc_);
  gradphi.resize(nc_*2);
  nablaphi.resize(nc_*3);
				
  double Jac[2][2]={{0,0},{0,0}};
  double JacI[2][2];
  const double *dxi=dphidxi[ig];
  const double *deta=dphideta[ig];
  for (int inode=0; inode<nc_; inode++,dxi++,deta++){
    Jac[0][0]+=(*dxi)*vt[0][inode];
    Jac[0][1]+=(*dxi)*vt[1][inode];
    Jac[1][0]+=(*deta)*vt[0][inode];
    Jac[1][1]+=(*deta)*vt[1][inode];
  }
  double det=(Jac[0][0]*Jac[1][1]-Jac[0][1]*Jac[1][0]);

  JacI[0][0]= Jac[1][1]/det;
  JacI[0][1]=-Jac[0][1]/det;
  JacI[1][0]=-Jac[1][0]/det;
  JacI[1][1]= Jac[0][0]/det;

  Weight=det*GaussWeight[ig];
     
  dxi=dphidxi[ig];
  deta=dphideta[ig];
  
  const double *dxi2=d2phidxi2[ig];
  const double *deta2=d2phideta2[ig];
  const double *dxideta=d2phidxideta[ig];
    
  for (int inode=0; inode<nc_; inode++, dxi++, deta++, dxi2++, deta2++, dxideta++) {
    
    other_phi[inode]=phi[ig][inode];
     
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


void elem_type_2D::JacobianSur_AD(const vector < vector < adept::adouble > > &vt, const unsigned &ig, adept::adouble &Weight, 
				  vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &normal) const {
  other_phi.resize(nc_);
  normal.resize(3);
  
  adept::adouble Jac[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
 
  const double *dfx=dphidxi[ig];
  const double *dfy=dphideta[ig];
  
  for(int inode=0; inode<nc_; inode++,dfx++,dfy++){       
    Jac[0][0] += (*dfx)*vt[0][inode];
    Jac[1][0] += (*dfx)*vt[1][inode];
    Jac[2][0] += (*dfx)*vt[2][inode];
    
    Jac[0][1] += (*dfy)*vt[0][inode];
    Jac[1][1] += (*dfy)*vt[1][inode];
    Jac[2][1] += (*dfy)*vt[2][inode];
  }
    
    //   normal module
  adept::adouble nx = Jac[1][0]*Jac[2][1] - Jac[1][1]*Jac[2][0];
  adept::adouble ny = Jac[0][1]*Jac[2][0] - Jac[2][1]*Jac[0][0]; 
  adept::adouble nz = Jac[0][0]*Jac[1][1] - Jac[0][1]*Jac[1][0];
  adept::adouble invModn = 1./sqrt(nx*nx + ny*ny + nz*nz);

  normal[0] =  (nx)*invModn;
  normal[1] =  (ny)*invModn;
  normal[2] =  (nz)*invModn;
    
  Jac[0][2] = normal[0];
  Jac[1][2] = normal[1];
  Jac[2][2] = normal[2];
    
  //the determinant of the matrix is the area 
  adept::adouble det=(Jac[0][0]*(Jac[1][1]*Jac[2][2]-Jac[1][2]*Jac[2][1])+
		      Jac[0][1]*(Jac[1][2]*Jac[2][0]-Jac[1][0]*Jac[2][2])+
		      Jac[0][2]*(Jac[1][0]*Jac[2][1]-Jac[1][1]*Jac[2][0]));

  Weight=det*GaussWeight[ig];

  for(int inode=0;inode<nc_;inode++){
    other_phi[inode]=phi[ig][inode];
  }
  
  //TODO warning the surface gradient is missing!!!!!!!!!!!!!!!
}

void elem_type_2D::JacobianSur(const vector < vector < double > > &vt, const unsigned &ig, double &Weight, 
			       vector < double > &other_phi, vector < double > &gradphi, vector < double > &normal) const {
  other_phi.resize(nc_);
  normal.resize(3);
  
  double Jac[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
 
  const double *dfx=dphidxi[ig];
  const double *dfy=dphideta[ig];
  
  for(int inode=0; inode<nc_; inode++,dfx++,dfy++){       
    Jac[0][0] += (*dfx)*vt[0][inode];
    Jac[1][0] += (*dfx)*vt[1][inode];
    Jac[2][0] += (*dfx)*vt[2][inode];
    
    Jac[0][1] += (*dfy)*vt[0][inode];
    Jac[1][1] += (*dfy)*vt[1][inode];
    Jac[2][1] += (*dfy)*vt[2][inode];
  }
    
    //   normal module
  double nx = Jac[1][0]*Jac[2][1] - Jac[1][1]*Jac[2][0];
  double ny = Jac[0][1]*Jac[2][0] - Jac[2][1]*Jac[0][0]; 
  double nz = Jac[0][0]*Jac[1][1] - Jac[0][1]*Jac[1][0];
  double invModn = 1./sqrt(nx*nx + ny*ny + nz*nz);

  normal[0] =  (nx)*invModn;
  normal[1] =  (ny)*invModn;
  normal[2] =  (nz)*invModn;
    
  Jac[0][2] = normal[0];
  Jac[1][2] = normal[1];
  Jac[2][2] = normal[2];
    
  //the determinant of the matrix is the area 
  double det=(Jac[0][0]*(Jac[1][1]*Jac[2][2]-Jac[1][2]*Jac[2][1])+
		      Jac[0][1]*(Jac[1][2]*Jac[2][0]-Jac[1][0]*Jac[2][2])+
		      Jac[0][2]*(Jac[1][0]*Jac[2][1]-Jac[1][1]*Jac[2][0]));

  Weight=det*GaussWeight[ig];

  for(int inode=0;inode<nc_;inode++){
    other_phi[inode]=phi[ig][inode];
  }
  
  //TODO warning the surface gradient is missing!!!!!!!!!!!!!!!
}

//---------------------------------------------------------------------------------------------------------
void elem_type_3D::Jacobian(const vector < vector < double > > &vt,const unsigned &ig, double &Weight, 
			      vector < double > &other_phi, vector < double > &gradphi, vector < double > &nablaphi) const{
  
  other_phi.resize(nc_);
  gradphi.resize(nc_*3);
  nablaphi.resize(nc_*6);
				
				
  double Jac[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
  double JacI[3][3];
  
  const double *dxi=dphidxi[ig];
  const double *deta=dphideta[ig];
  const double *dzeta=dphidzeta[ig];

  for (int inode=0; inode<nc_; inode++,dxi++,deta++,dzeta++) {
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
  double det=(Jac[0][0]*(Jac[1][1]*Jac[2][2]-Jac[1][2]*Jac[2][1])+
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

  Weight=det*GaussWeight[ig];
 
  dxi=dphidxi[ig];
  deta=dphideta[ig];
  dzeta=dphidzeta[ig];
  
  const double *dxi2=d2phidxi2[ig];
  const double *deta2=d2phideta2[ig];
  const double *dzeta2=d2phidzeta2[ig];
  const double *dxideta=d2phidxideta[ig];
  const double *detadzeta=d2phidetadzeta[ig];
  const double *dzetadxi=d2phidzetadxi[ig];
  
  for (int inode=0; inode<nc_; inode++, dxi++,deta++,dzeta++,dxi2++,deta2++,dzeta2++,dxideta++,detadzeta++,dzetadxi++) {
  
    other_phi[inode]=phi[ig][inode];
    
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


void elem_type_3D::Jacobian_AD(const vector < vector < adept::adouble > > &vt,const unsigned &ig, adept::adouble &Weight, 
			      vector < double > &other_phi, vector < adept::adouble > &gradphi, vector < adept::adouble > &nablaphi) const{
  
  other_phi.resize(nc_);
  gradphi.resize(nc_*3);
  nablaphi.resize(nc_*6);
				
				
  adept::adouble Jac[3][3]={{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
  adept::adouble JacI[3][3];
  
  const double *dxi=dphidxi[ig];
  const double *deta=dphideta[ig];
  const double *dzeta=dphidzeta[ig];

  for (int inode=0; inode<nc_; inode++,dxi++,deta++,dzeta++) {
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
  adept::adouble det=(Jac[0][0]*(Jac[1][1]*Jac[2][2]-Jac[1][2]*Jac[2][1])+
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

  Weight=det*GaussWeight[ig];
 
  dxi=dphidxi[ig];
  deta=dphideta[ig];
  dzeta=dphidzeta[ig];
  
  const double *dxi2=d2phidxi2[ig];
  const double *deta2=d2phideta2[ig];
  const double *dzeta2=d2phidzeta2[ig];
  const double *dxideta=d2phidxideta[ig];
  const double *detadzeta=d2phidetadzeta[ig];
  const double *dzetadxi=d2phidzetadxi[ig];
  
  for (int inode=0; inode<nc_; inode++, dxi++,deta++,dzeta++,dxi2++,deta2++,dzeta2++,dxideta++,detadzeta++,dzetadxi++) {
  
    other_phi[inode]=phi[ig][inode];
    
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










} //end namespace femus


