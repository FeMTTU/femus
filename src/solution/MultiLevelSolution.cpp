/*=========================================================================

Program: FEMUS
Module: MultiLevelProblem
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
#include "MultiLevelSolution.hpp"
#include "ElemType.hpp"
#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "FEMTTUConfig.h"
#include "ParsedFunction.hpp"


//C++ include
#include <iostream>
#include <iomanip>


namespace femus {


using std::cout;
using std::endl;

//---------------------------------------------------------------------------------------------------
MultiLevelSolution::~MultiLevelSolution() {

  for (unsigned i=0; i<_gridn; i++) {
    _solution[i]->FreeSolutionVectors();
    delete _solution[i];
  }

  for (unsigned i=0; i<_SolName.size(); i++) delete [] _SolName[i];
  for (unsigned i=0; i<_SolName.size(); i++) delete [] _BdcType[i];


};

//---------------------------------------------------------------------------------------------------
MultiLevelSolution::MultiLevelSolution( MultiLevelMesh *ml_msh):
  _gridn(ml_msh->GetNumberOfLevels()),
  _ml_msh(ml_msh)
{
  _solution.resize(_gridn);
  
  for (unsigned i=0; i<_gridn; i++) {
    _solution[i]=new Solution(_ml_msh->GetLevel(i));
  }

  _bdc_func_set=false;

  
}

void MultiLevelSolution::AddSolutionLevel(){
  // add level solution
  _solution.resize(_gridn+1);
  _solution[_gridn]=new Solution(_ml_msh->GetLevel(_gridn));
  // add all current solutions and initialize to zero
  for(unsigned i=0;i<_SolName.size();i++){
    _solution[_gridn]->AddSolution(_SolName[i],_family[i],_order[i],_SolTmorder[i],_PdeType[i]);
  }
  for(unsigned i=0;i<_SolName.size();i++){
    _solution[_gridn]->ResizeSolutionVector(_SolName[i]);
    BuildProlongatorMatrix(_gridn,i);   
        
    _solution[_gridn]->_Sol[i]->zero();
    if (_SolTmorder[i]==2) {
      _solution[_gridn]->_SolOld[i]->zero();
    }
  }
  _gridn++;
  unsigned  grid0=_gridn-1;
  GenerateBdc(grid0);
}

//---------------------------------------------------------------------------------------------------
void MultiLevelSolution::AddSolution(const char name[], const FEFamily fefamily, const FEOrder order,
				     unsigned tmorder, const bool &PdeType) {
  
  unsigned n=_SolType.size();
  _SolType.resize(n+1u);
  _family.resize(n+1u);
  _order.resize(n+1u);
  _SolName.resize(n+1u);
  _BdcType.resize(n+1u);
  _SolTmorder.resize(n+1u);
  _PdeType.resize(n+1u);
  _TestIfPressure.resize(n+1u);
  _TestIfDisplacement.resize(n+1u);

  _TestIfDisplacement[n]=0;
  _TestIfPressure[n]=0;
  _family[n] = fefamily;
  _order[n] = order;
  _SolType[n] = order - ((fefamily==LAGRANGE)?1:0) + fefamily*3;     
  _SolName[n]  = new char [8];
  _BdcType[n]  = new char [20];
  strcpy(_SolName[n],name);
  _SolTmorder[n]=tmorder;
  _PdeType[n]=PdeType;
 
  cout << " Add variable " << std::setw(3) << _SolName[n] << " discretized with FE type "
       << std::setw(12) << order << " and time discretzation order " << tmorder << endl;

  for (unsigned ig=0; ig<_gridn; ig++) {
    _solution[ig]->AddSolution(_SolName[n],_family[n],_order[n],_SolTmorder[n],_PdeType[n]);
    //_solution[ig]->AddSolution(name,fefamily,order,tmorder,Pde_type);
  }
}

//---------------------------------------------------------------------------------------------------
void MultiLevelSolution::AssociatePropertyToSolution(const char solution_name[], const char solution_property[]){
  unsigned index=GetIndex(solution_name);
  if( !strcmp(solution_property,"pressure") || !strcmp(solution_property,"Pressure") ) _TestIfPressure[index]=1;
  else if( !strcmp(solution_property,"displacement") || !strcmp(solution_property,"Displacement") ) _TestIfDisplacement[index]=1;
  else if( !strcmp(solution_property,"default") || !strcmp(solution_property,"Default") ) {
    _TestIfPressure[index]=0;
    _TestIfDisplacement[index]=0;
  }
  else {
    cout<<"Error invalid property in function MultiLevelProblem::AssociatePropertyToSolution"<<endl;
    exit(0);
  }
}

// *******************************************************
void MultiLevelSolution::Initialize(const char name[], initfunc func) {
  
  unsigned i_start;
  unsigned i_end;
  if (!strcmp(name,"All")) {
    i_start=0;
    i_end=_SolType.size();
  } else {
    i_start=GetIndex(name);
    i_end=i_start+1u;
  }
  
  double value;
  for (unsigned i=i_start; i<i_end; i++) {
    //CheckVectorSize(i);
    unsigned sol_type = _SolType[i];
    for (unsigned ig=0; ig<_gridn; ig++) {
      unsigned num_el = _ml_msh->GetLevel(ig)->GetElementNumber();
      _solution[ig]->ResizeSolutionVector(_SolName[i]);
      if (ig>0) BuildProlongatorMatrix(ig,i);     
      //for parallel
      for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
        for (int iel=_ml_msh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[isdom]; 
	     iel < _ml_msh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
	  unsigned kel_gmt = _ml_msh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];    
	  unsigned sol_ord = _ml_msh->GetLevel(ig)->GetEndIndex(_SolType[i]);
	  unsigned nloc_dof= _ml_msh->GetLevel(ig)->el->GetElementDofNumber(kel_gmt,sol_ord);
	  if(sol_type<3) {
            for(int j=0; j<nloc_dof; j++) {
	      unsigned inode=(sol_type<3)?(_ml_msh->GetLevel(ig)->el->GetElementVertexIndex(kel_gmt,j)-1u):(kel_gmt+j*num_el);
	      int inode_Metis=_ml_msh->GetLevel(ig)->GetMetisDof(inode,sol_type);
	      unsigned icoord_Metis=_ml_msh->GetLevel(ig)->GetMetisDof(inode,2);
	      double xx=(*_ml_msh->GetLevel(ig)->_coordinate->_Sol[0])(icoord_Metis);  
	      double yy=(*_ml_msh->GetLevel(ig)->_coordinate->_Sol[1])(icoord_Metis);
	      double zz=(*_ml_msh->GetLevel(ig)->_coordinate->_Sol[2])(icoord_Metis);
	      
	      if(func) {
		value = (sol_type<3)?func(xx,yy,zz):0;
	      }
	      else {
		value = 0.;
	      }
	      _solution[ig]->_Sol[i]->set(inode_Metis,value);
	      if (_SolTmorder[i]==2) {
		_solution[ig]->_SolOld[i]->set(inode_Metis,value);
	      }
	    }
	  }
	}
      }
      _solution[ig]->_Sol[i]->close();
      if (_SolTmorder[i]==2) {
	_solution[ig]->_SolOld[i]->close();
      }
    }
  }
}
 
//---------------------------------------------------------------------------------------------------
unsigned MultiLevelSolution::GetIndex(const char name[]) const {
  unsigned index=0;
  while (strcmp(_SolName[index],name)) {
    index++;
    if (index==_SolType.size()) {
      cout<<"error! invalid solution name "<< name <<"entry GetIndex(...)"<<endl;
      exit(0);
    }
  }
  return index;
}

// *******************************************************
unsigned MultiLevelSolution::GetSolutionType(const char name[]) {
  unsigned index=0;
  while (strcmp(_SolName[index],name)) {
    index++;
    if (index==_SolType.size()) {
      cout<<"error! invalid name entry GetSolType(...)"<<endl;
      exit(0);
    }
  }
  return _SolType[index];
}

//---------------------------------------------------------------------------------------------------
// This routine generates the matrix for the projection of the solution to finer grids  
//---------------------------------------------------------------------------------------------------

void MultiLevelSolution::BuildProlongatorMatrix(unsigned gridf, unsigned SolIndex) {
  
  if (gridf<1) {
    cout<<"Error! In function \"BuildProlongatorMatrix\" argument less then 1"<<endl;
    exit(0);
  }
  
  unsigned TypeIndex=_SolType[SolIndex];
    
  if(_solution[gridf]->_ProjMatFlag[TypeIndex]==0){
    _solution[gridf]->_ProjMatFlag[TypeIndex]=1;

    int nf     = _ml_msh->GetLevel(gridf)->MetisOffset[_SolType[SolIndex]][_nprocs];
    int nc     = _ml_msh->GetLevel(gridf-1)->MetisOffset[_SolType[SolIndex]][_nprocs];
    int nf_loc = _ml_msh->GetLevel(gridf)->own_size[_SolType[SolIndex]][_iproc];
    int nc_loc = _ml_msh->GetLevel(gridf-1)->own_size[_SolType[SolIndex]][_iproc]; 

    _solution[gridf]->_ProjMat[TypeIndex] = SparseMatrix::build().release();
    _solution[gridf]->_ProjMat[TypeIndex]->init(nf,nc,nf_loc,nc_loc,27,27);
 
    // loop on the coarse grid 
    for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
      for (int iel_mts=_ml_msh->GetLevel(gridf-1)->IS_Mts2Gmt_elem_offset[isdom]; 
	   iel_mts < _ml_msh->GetLevel(gridf-1)->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = _ml_msh->GetLevel(gridf-1)->IS_Mts2Gmt_elem[iel_mts];
	if(_ml_msh->GetLevel(gridf-1)->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
	  short unsigned ielt=_ml_msh->GetLevel(gridf-1)->el->GetElementType(iel);
	  _ml_msh->_type_elem[ielt][_SolType[SolIndex]]->prolongation(*_ml_msh->GetLevel(gridf),*_ml_msh->GetLevel(gridf-1),iel,
								      _solution[gridf]->_ProjMat[TypeIndex]); 
	}
      }
    }
    _solution[gridf]->_ProjMat[TypeIndex]->close();
  }
}

//---------------------------------------------------------------------------------------------------
void MultiLevelSolution::AttachSetBoundaryConditionFunction ( bool (* SetBoundaryConditionFunction) (const double &x, const double &y, const double &z,const char name[], 
												     double &value, const int FaceName, const double time) ) {
  _bdc_func_set = true;
  _SetBoundaryConditionFunction = SetBoundaryConditionFunction;
  return;
}

//---------------------------------------------------------------------------------------------------
void MultiLevelSolution::InitializeBdc() {
  int nvars = _SolType.size();
  int nfaces = _ml_msh->GetLevel(0)->_boundaryinfo.size();
  _boundaryconditions.resize(nvars);
  _ishomogeneous.resize(nvars);
  _nonhomogeneousbcfunction.resize(nvars);
  for(int i=0; i<nvars; i++) {
    _boundaryconditions[i].resize(nfaces);
    _ishomogeneous[i].resize(nfaces);
    _nonhomogeneousbcfunction[i].resize(nfaces);
    for(int j=0; j<nfaces; j++) {
       _boundaryconditions[i][j] = DIRICHLET;
       _ishomogeneous[i][j] = true;
       _nonhomogeneousbcfunction[i][j] = NULL;
    }
  }
}

//---------------------------------------------------------------------------------------------------
void MultiLevelSolution::SetBoundaryCondition(const std::string name, const std::string facename, 
					      const BDCType bdctype, const bool ishomogeneous, 
					      const bool istimedependent, FunctionBase* func) {
  
  if(ishomogeneous == false && func == NULL) {
    std::cout << " Error: you have set a non-homogeneous boundary condition withtout passing a valid function pointer" << std::endl;
    exit(1);  
  }
  
  unsigned int ivar = GetIndex(name.c_str());
  unsigned int iface = 0;
  
  std::map<unsigned int, std::string>::iterator iter;
  iter = _ml_msh->GetLevel(0)->_boundaryinfo.begin();
  
  for (iter = _ml_msh->GetLevel(0)->_boundaryinfo.begin(); iter != _ml_msh->GetLevel(0)->_boundaryinfo.end(); ++iter) {
    if( iter->second.compare(facename) == 0) {
      iface = iter->first; 
      break;
    }
  }
  if(iter == _ml_msh->GetLevel(0)->_boundaryinfo.end()) {
    std::cout << " Error: the facename " << facename << " does not exist!" << std::endl;
    exit(1);    
  }
  
  _boundaryconditions[ivar][iface]       = bdctype;
  _ishomogeneous[ivar][iface]            = ishomogeneous;
  _nonhomogeneousbcfunction[ivar][iface] = func;
   
}




void MultiLevelSolution::GenerateBdc(const unsigned &grid0) {
  
  const short unsigned NV1[6][2]= {{9,9},{6,6},{9,6},{3,3},{3,3},{1,1}};
  
  int nvars = _SolType.size();
  
  double time = 0;
  
  // 2 Default Neumann
  // 1 DD Dirichlet
  // 0 Dirichlet
  for (unsigned k=0; k<nvars; k++) {
  for (unsigned igridn=grid0; igridn<_gridn; igridn++) {
    if(_solution[igridn]->_ResEpsBdcFlag[k]){
      for (unsigned j=_ml_msh->GetLevel(igridn)->MetisOffset[_SolType[k]][_iproc]; j<_ml_msh->GetLevel(igridn)->MetisOffset[_SolType[k]][_iproc+1]; j++) {
	_solution[igridn]->_Bdc[k]->set(j,2.);
      }
      if (_SolType[k]<3) {  
	for(int isdom=_iproc; isdom<_iproc+1; isdom++) {   // 1 DD Dirichlet
	  for (int iel=_ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom]; 
	       iel < _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
	    unsigned kel_gmt = _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem[iel];
	    for (unsigned jface=0; jface<_ml_msh->GetLevel(igridn)->el->GetElementFaceNumber(kel_gmt); jface++) {
	      if (_ml_msh->GetLevel(igridn)->el->GetFaceElementIndex(kel_gmt,jface)==0) { //Domain Decomposition Dirichlet
		short unsigned ielt=_ml_msh->GetLevel(igridn)->el->GetElementType(kel_gmt);
		unsigned nv1=(!_TestIfPressure[k])?
		  NV1[ielt][jface<_ml_msh->GetLevel(igridn)->el->GetElementFaceNumber(kel_gmt,0)]:
		  _ml_msh->GetLevel(igridn)->el->GetElementDofNumber(iel,_ml_msh->GetLevel(igridn)->GetEndIndex(_SolType[k]));
		for (unsigned iv=0; iv<nv1; iv++) {
		  unsigned inode=(!_TestIfPressure[k])? 
		    _ml_msh->GetLevel(igridn)->el->GetFaceVertexIndex(kel_gmt,jface,iv)-1u:
		    _ml_msh->GetLevel(igridn)->el->GetElementVertexIndex(kel_gmt,iv)-1u;
		  unsigned inode_Metis=_ml_msh->GetLevel(igridn)->GetMetisDof(inode,_SolType[k]);
		  _solution[igridn]->_Bdc[k]->set(inode_Metis,1.);
		}
	      }
	    }
	  }
	}
	for(int isdom=_iproc; isdom<_iproc+1; isdom++) {  // 0 Dirichlet
	  for (int iel=_ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom]; 
	       iel < _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
	    unsigned kel_gmt = _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem[iel];
	    for (unsigned jface=0; jface<_ml_msh->GetLevel(igridn)->el->GetElementFaceNumber(kel_gmt); jface++) {
	      if (_ml_msh->GetLevel(igridn)->el->GetFaceElementIndex(kel_gmt,jface)<0) { 
		short unsigned ielt=_ml_msh->GetLevel(igridn)->el->GetElementType(kel_gmt);
		// unsigned nv1=NV1[ielt][jface<_ml_msh->GetLevel(igridn)->el->GetElementFaceNumber(kel_gmt,0)];
		unsigned nv1 = _ml_msh->GetLevel(igridn)->el->GetElementFaceDofNumber(kel_gmt,jface,_SolType[k]); 
		for (unsigned iv=0; iv<nv1; iv++) {
		  unsigned inode=_ml_msh->GetLevel(igridn)->el->GetFaceVertexIndex(kel_gmt,jface,iv)-1u;
		  unsigned inode_coord_Metis=_ml_msh->GetLevel(igridn)->GetMetisDof(inode,2);
		  double value = 0.;
		  double xx=(*_ml_msh->GetLevel(igridn)->_coordinate->_Sol[0])(inode_coord_Metis);  
		  double yy=(*_ml_msh->GetLevel(igridn)->_coordinate->_Sol[1])(inode_coord_Metis);
		  double zz=(*_ml_msh->GetLevel(igridn)->_coordinate->_Sol[2])(inode_coord_Metis);
		  unsigned int face = -(_ml_msh->GetLevel(igridn)->el->GetFaceElementIndex(kel_gmt,jface)+1)-1 ;
		  if(GetBoundaryCondition(k,face) == DIRICHLET) {
		     unsigned inode_Metis=_ml_msh->GetLevel(igridn)->GetMetisDof(inode,_SolType[k]);
		    _solution[igridn]->_Bdc[k]->set(inode_Metis,0.);   
		    if(!Ishomogeneous(k,face)) {
		      ParsedFunction* bdcfunc = (ParsedFunction*)(GetBdcFunction(k,face));
		      double xyzt[4];
		      xyzt[0] = xx;
		      xyzt[1] = yy;
		      xyzt[2] = zz;
		      xyzt[3] = time;
		      value = (*bdcfunc)(xyzt);
		      _solution[igridn]->_Sol[k]->set(inode_Metis,value); 
		    }
		    else {
		      _solution[igridn]->_Sol[k]->set(inode_Metis,0.);  
		    }    
		  }
		}
	      }
	    }
	  }
	}
      }
      else if(_TestIfPressure[k]){
	for(int isdom=_iproc; isdom<_iproc+1; isdom++) {   // 1 DD Dirichlet for pressure variable only
	  unsigned nel=_ml_msh->GetLevel(igridn)->GetElementNumber();
	  for (int iel=_ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom]; 
	       iel < _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
	    unsigned kel_gmt = _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem[iel];
	    for (unsigned jface=0; jface<_ml_msh->GetLevel(igridn)->el->GetElementFaceNumber(kel_gmt); jface++) {
	      if (_ml_msh->GetLevel(igridn)->el->GetFaceElementIndex(kel_gmt,jface)==0) { //Domain Decomposition Dirichlet
		short unsigned ielt=_ml_msh->GetLevel(igridn)->el->GetElementType(kel_gmt);
		unsigned nv1=_ml_msh->GetLevel(igridn)->el->GetElementDofNumber(kel_gmt,_ml_msh->GetLevel(igridn)->GetEndIndex(_SolType[k]));
		for (unsigned iv=0; iv<nv1; iv++) {
		  unsigned inode=(kel_gmt+iv*nel);
		  unsigned inode_Metis=_ml_msh->GetLevel(igridn)->GetMetisDof(inode,_SolType[k]);
		  _solution[igridn]->_Bdc[k]->set(inode_Metis,1.);
		}
	      }
	    }
	  }
	}
      }	
      _solution[igridn]->_Sol[k]->close();
      _solution[igridn]->_Bdc[k]->close();
    }
   }
  }
  
  
  
}




//---------------------------------------------------------------------------------------------------
void MultiLevelSolution::GenerateBdc(const char name[], const char bdc_type[]) {

  if(_bdc_func_set==false) {
    cout << "Error: The boundary condition user-function is not set! Please call the AttachSetBoundaryConditionFunction routine" 
         << endl;
    exit(1); 
  }

  unsigned i_start;
  unsigned i_end;
  if (!strcmp(name,"All")) {
    i_start=0;
    i_end=_SolType.size();
    for (unsigned k=i_start; k<i_end; k++) {
      if(_solution[0]->_ResEpsBdcFlag[k]){
	sprintf(_BdcType[k],"Steady");
	cout << " Set " << std::setw(15) << _BdcType[k] << " Boundary_condition"
	     << " for variable " << std::setw(3) << _SolName[k] << endl;
      }
      else {
	sprintf(_BdcType[k],"Not-available");
      }
	
    }
  } 
  else {
    i_start=GetIndex(name);
    i_end=i_start+1u;
    if(_solution[0]->_ResEpsBdcFlag[i_start]){
      if (!strcmp(bdc_type,"Steady")) {
	strcpy(_BdcType[i_start],bdc_type);
      } else if (!strcmp(bdc_type,"Time_dependent")) {
	strcpy(_BdcType[i_start],bdc_type);
      } else {
	cout << "Error! Invalid boundary condition specified for " << _SolName[i_start]
	     << " in GenerateBdc function" << endl;
	exit(1);
      }
      cout << " Set " << std::setw(14) <<_BdcType[i_start] << " Boundary_condition"
	   << " for variable " << std::setw(3) << _SolName[i_start] << endl;
    }
    else {
      sprintf(_BdcType[i_start],"Not-available");
    }
  }
  for (unsigned i=i_start; i<i_end; i++) {
    GenerateBdc(i,0.);
  }
}

//---------------------------------------------------------------------------------------------------
void MultiLevelSolution::UpdateBdc(const double time) {

  for (int k=0; k<_SolName.size(); k++) {
    if (!strcmp(_BdcType[k],"Time_dependent")) {
      GenerateBdc(k,time);

    }
  }
}

//---------------------------------------------------------------------------------------------------
void MultiLevelSolution::GenerateBdc(const unsigned int k, const double time) {
  
  const short unsigned NV1[6][2]= {{9,9},{6,6},{9,6},{3,3},{3,3},{1,1}};
  
       
  // 2 Default Neumann
  // 1 DD Dirichlet
  // 0 Dirichlet
  for (unsigned igridn=0; igridn<_gridn; igridn++) {
    if(_solution[igridn]->_ResEpsBdcFlag[k]){
      for (unsigned j=_ml_msh->GetLevel(igridn)->MetisOffset[_SolType[k]][_iproc]; j<_ml_msh->GetLevel(igridn)->MetisOffset[_SolType[k]][_iproc+1]; j++) {
	_solution[igridn]->_Bdc[k]->set(j,2.);
      }
      if (_SolType[k]<3) {  
	for(int isdom=_iproc; isdom<_iproc+1; isdom++) {   // 1 DD Dirichlet
	  for (int iel=_ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom]; 
	       iel < _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
	    unsigned kel_gmt = _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem[iel];
	    for (unsigned jface=0; jface<_ml_msh->GetLevel(igridn)->el->GetElementFaceNumber(kel_gmt); jface++) {
	      if (_ml_msh->GetLevel(igridn)->el->GetFaceElementIndex(kel_gmt,jface)==0) { //Domain Decomposition Dirichlet
		short unsigned ielt=_ml_msh->GetLevel(igridn)->el->GetElementType(kel_gmt);
		unsigned nv1=(!_TestIfPressure[k])?
		  //NV1[ielt][jface<_ml_msh->GetLevel(igridn)->el->GetElementFaceNumber(kel_gmt,0)]:
		  _ml_msh->GetLevel(igridn)->el->GetElementFaceDofNumber(kel_gmt,jface,_SolType[k]):
		  _ml_msh->GetLevel(igridn)->el->GetElementDofNumber(iel,_ml_msh->GetLevel(igridn)->GetEndIndex(_SolType[k]));
		for (unsigned iv=0; iv<nv1; iv++) {
		  unsigned inode=(!_TestIfPressure[k])? 
		    _ml_msh->GetLevel(igridn)->el->GetFaceVertexIndex(kel_gmt,jface,iv)-1u:
		    _ml_msh->GetLevel(igridn)->el->GetElementVertexIndex(kel_gmt,iv)-1u;
		  unsigned inode_Metis=_ml_msh->GetLevel(igridn)->GetMetisDof(inode,_SolType[k]);
		  _solution[igridn]->_Bdc[k]->set(inode_Metis,1.);
		}
	      }
	    }
	  }
	}
	for(int isdom=_iproc; isdom<_iproc+1; isdom++) {  // 0 Dirichlet
	  for (int iel=_ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom]; 
	       iel < _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
	    unsigned kel_gmt = _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem[iel];
	    for (unsigned jface=0; jface<_ml_msh->GetLevel(igridn)->el->GetElementFaceNumber(kel_gmt); jface++) {
	      if (_ml_msh->GetLevel(igridn)->el->GetFaceElementIndex(kel_gmt,jface)<0) { //Dirichlet
		short unsigned ielt=_ml_msh->GetLevel(igridn)->el->GetElementType(kel_gmt);
		// unsigned nv1=NV1[ielt][jface<_ml_msh->GetLevel(igridn)->el->GetElementFaceNumber(kel_gmt,0)];
		unsigned nv1 = _ml_msh->GetLevel(igridn)->el->GetElementFaceDofNumber(kel_gmt,jface,_SolType[k]); 
		for (unsigned iv=0; iv<nv1; iv++) {
		  unsigned inode=_ml_msh->GetLevel(igridn)->el->GetFaceVertexIndex(kel_gmt,jface,iv)-1u;
		  unsigned inode_coord_Metis=_ml_msh->GetLevel(igridn)->GetMetisDof(inode,2);
		  double value;
		  double xx=(*_ml_msh->GetLevel(igridn)->_coordinate->_Sol[0])(inode_coord_Metis);  
		  double yy=(*_ml_msh->GetLevel(igridn)->_coordinate->_Sol[1])(inode_coord_Metis);
		  double zz=(*_ml_msh->GetLevel(igridn)->_coordinate->_Sol[2])(inode_coord_Metis);
		  bool test=_SetBoundaryConditionFunction(xx,yy,zz,_SolName[k],value,-(_ml_msh->GetLevel(igridn)->el->GetFaceElementIndex(kel_gmt,jface)+1),time);
		  if (test) {
		    unsigned inode_Metis=_ml_msh->GetLevel(igridn)->GetMetisDof(inode,_SolType[k]);
		    _solution[igridn]->_Bdc[k]->set(inode_Metis,0.);
		    _solution[igridn]->_Sol[k]->set(inode_Metis,value);
		  }
		}
	      }
	    }
	  }
	}
      }
      else if(_TestIfPressure[k]){
	for(int isdom=_iproc; isdom<_iproc+1; isdom++) {   // 1 DD Dirichlet for pressure variable only
	  unsigned nel=_ml_msh->GetLevel(igridn)->GetElementNumber();
	  for (int iel=_ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom]; 
	       iel < _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
	    unsigned kel_gmt = _ml_msh->GetLevel(igridn)->IS_Mts2Gmt_elem[iel];
	    for (unsigned jface=0; jface<_ml_msh->GetLevel(igridn)->el->GetElementFaceNumber(kel_gmt); jface++) {
	      if (_ml_msh->GetLevel(igridn)->el->GetFaceElementIndex(kel_gmt,jface)==0) { //Domain Decomposition Dirichlet
		short unsigned ielt=_ml_msh->GetLevel(igridn)->el->GetElementType(kel_gmt);
		unsigned nv1=_ml_msh->GetLevel(igridn)->el->GetElementDofNumber(kel_gmt,_ml_msh->GetLevel(igridn)->GetEndIndex(_SolType[k]));
		for (unsigned iv=0; iv<nv1; iv++) {
		  unsigned inode=(kel_gmt+iv*nel);
		  unsigned inode_Metis=_ml_msh->GetLevel(igridn)->GetMetisDof(inode,_SolType[k]);
		  _solution[igridn]->_Bdc[k]->set(inode_Metis,1.);
		}
	      }
	    }
	  }
	}
      }	
      _solution[igridn]->_Sol[k]->close();
      _solution[igridn]->_Bdc[k]->close();
    }
  }
  
}



} //end namespace femus


