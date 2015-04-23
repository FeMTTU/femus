/*=========================================================================

 Program: FEMUS
 Module: GMVWriter
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
#include "GMVWriter.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>   
#include <algorithm>  
#include <cstring> 
#include "Files.hpp"


namespace femus {



GMVWriter::GMVWriter(MultiLevelSolution * ml_sol): Writer(ml_sol)
{
  _debugOutput = false;
}

GMVWriter::GMVWriter(MultiLevelMesh * ml_mesh): Writer(ml_mesh)
{
  _debugOutput = false;
}

GMVWriter::~GMVWriter()
{
  
}

void GMVWriter::write(const std::string output_path, const char order[], const std::vector<std::string>& vars, const unsigned time_step) const { 
  
  unsigned igridn = _gridn; // aggiunta da me
      
  if (igridn==0) igridn=_gridn;
  
  unsigned igridr=(_gridr <= igridn)?_gridr:igridn;

  // ********** linear -> index==0 *** quadratic -> index==1 **********
  unsigned index=(strcmp(order,"linear"))?1:0;

  std::string filename_prefix;
  if( _ml_sol != NULL ) filename_prefix = "sol";
  else filename_prefix = "mesh";
  
  std::ostringstream filename;
    filename << output_path << "/" << filename_prefix << ".level" << _gridn << "." << time_step << "." << order << ".gmv"; 
    
  std::ofstream fout;
  
  if(_iproc!=0) {
    fout.rdbuf();   //redirect to dev_null
  }
  else {
    fout.open(filename.str().c_str());
    if (fout.is_open()) {
      std::cout << std::endl << " The output is printed to file " << filename.str() << " in GMV format" << std::endl; 
    }
    else {
      std::cout << std::endl << " The output file "<< filename.str() <<" cannot be opened.\n";
      abort();
    }
  }  
  
  unsigned nvt=0;
  unsigned nvt_max=0;
  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
    unsigned nvt_ig=_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
    nvt_max=(nvt_max>nvt_ig)?nvt_max:nvt_ig;
    nvt+=nvt_ig;
  }
  
  double *var_nd=new double [nvt_max+1]; //TO FIX Valgrind complaints! In reality it should be only nvt
  vector <NumericVector*> Mysol(igridn);
  for(unsigned ig=igridr-1u; ig<_gridn; ig++) {
    Mysol[ig] = NumericVector::build().release();
    Mysol[ig]->init(_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs],_ml_mesh->GetLevel(ig)->own_size[index][_iproc],true,AUTOMATIC);
  }
     
  // ********** Header **********
  char *det= new char[10];
  sprintf(det,"%s","gmvinput");
  fout.write((char *)det,sizeof(char)*8);
  sprintf(det,"%s","ieeei4r8");
  fout.write((char *)det,sizeof(char)*8);

  // ********** Start printing node coordinates  **********
  sprintf(det,"%s","nodes");
  fout.write((char *)det,sizeof(char)*8);
  fout.write((char *)&nvt,sizeof(unsigned));
    
  for (int i=0; i<3; i++) {
    for (unsigned ig=igridr-1u; ig<igridn; ig++) {
      Mysol[ig]->matrix_mult(*_ml_mesh->GetLevel(ig)->_coordinate->_Sol[i],
			     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,2) );
      vector <double> v_local;
      Mysol[ig]->localize_to_one(v_local,0);
      unsigned nvt_ig=_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];      
      if(_iproc==0){ 
	for (unsigned ii=0; ii<nvt_ig; ii++) 
	  var_nd[ii]= v_local[ii];
      }
      if (_ml_sol != NULL && _moving_mesh  && _ml_mesh->GetLevel(0)->GetDimension() > i)  {
	unsigned indDXDYDZ=_ml_sol->GetIndex(_moving_vars[i].c_str());
	Mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Sol[indDXDYDZ],
			       *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index, _ml_sol->GetSolutionType(indDXDYDZ)) );
	Mysol[ig]->localize_to_one(v_local,0);
	unsigned nvt_ig=_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];      
	if(_iproc==0){ 
	  for (unsigned ii=0; ii<nvt_ig; ii++) 
	    var_nd[ii]+= v_local[ii];
	}
      }
      fout.write((char *)&var_nd[0],nvt_ig*sizeof(double)); 
    }
  }
  // ********** End printing node coordinates  **********

  // ********** Start printing cell connectivity  **********
  const int eltp[2][6]= {{8,4,6,4,3,2},{20,10,15,8,6,3}};
  sprintf(det,"%s","cells");
  fout.write((char *)det,sizeof(char)*8);

  unsigned nel=0;
  for (unsigned ig=igridr-1u; ig<igridn-1u; ig++)
    nel+=( _ml_mesh->GetLevel(ig)->GetNumberOfElements() - _ml_mesh->GetLevel(ig)->el->GetRefinedElementNumber());
  nel+=_ml_mesh->GetLevel(igridn-1u)->GetNumberOfElements();
  fout.write((char *)&nel,sizeof(unsigned));

  unsigned topology[27];
  unsigned offset=1;
  
  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
    for (unsigned ii=0; ii<_ml_mesh->GetLevel(ig)->GetNumberOfElements(); ii++) {
      if ( ig == igridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
        short unsigned ielt=_ml_mesh->GetLevel(ig)->el->GetElementType(ii);
        if (ielt==0) sprintf(det,"phex%d",eltp[index][0]);
        else if (ielt==1) sprintf(det,"ptet%d",eltp[index][1]);
        else if (ielt==2) sprintf(det,"pprism%d",eltp[index][2]);
        else if (ielt==3) {
          if (eltp[index][3]==8) sprintf(det,"%dquad",eltp[index][3]);
          else sprintf(det,"quad");
        } else if (ielt==4) {
          if (eltp[index][4]==6) sprintf(det,"%dtri",eltp[index][4]);
          else sprintf(det,"tri");
        } else if (ielt==5) {
          if (eltp[index][5]==3) sprintf(det,"%dline",eltp[index][5]);
          else sprintf(det,"line");
        }
        fout.write((char *)det,sizeof(char)*8);
        fout.write((char *)&NVE[ielt][index],sizeof(unsigned));
	for(unsigned j=0;j<NVE[ielt][index];j++){
	  
	  unsigned jnode=_ml_mesh->GetLevel(ig)->el->GetElementVertexIndex(ii,j)-1u;
	  unsigned jnode_Metis = _ml_mesh->GetLevel(ig)->GetMetisDof(jnode,index);
	  	  
	  topology[j]=jnode_Metis+offset;
	}
	fout.write((char *)topology,sizeof(unsigned)*NVE[ielt][index]);
      }
    }
    offset+=_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
  }
  // ********** End printing cell connectivity  **********
  
  double *var_el=new double [nel+1]; //TO FIX Valgrind complaints! In reality it should be only nel
  
  // ********** Start printing Variables **********
  const unsigned zero=0u;
  const unsigned one=1u;
  sprintf(det,"%s","variable");
  fout.write((char *)det,sizeof(char)*8);

  // ********** Start printing Regions **********
  strcpy(det,"Regions");
  fout.write((char *)det,sizeof(char)*8);
  fout.write((char *)&zero,sizeof(unsigned));

  int icount=0;
  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
    for (unsigned ii=0; ii<_ml_mesh->GetLevel(ig)->GetNumberOfElements(); ii++) {
      if ( ig==igridn-1u || 0==_ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
	var_el[icount]=_ml_mesh->GetLevel(ig)->el->GetElementGroup(ii);
        icount++;
      }
    }
  }
  fout.write((char *)&var_el[0],nel*sizeof(double));
  
  if(_nprocs>=1){
    strcpy(det,"Reg_proc");
    fout.write((char *)det,sizeof(char)*8);
    fout.write((char *)&zero,sizeof(unsigned));

    int icount=0;
    for (unsigned ig=igridr-1u; ig<igridn; ig++) {
      for (unsigned ii=0; ii<_ml_mesh->GetLevel(ig)->GetNumberOfElements(); ii++) {
	if ( ig==igridn-1u || 0==_ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
	  var_el[icount]=_ml_mesh->GetLevel(ig)->epart[ii];
	  icount++;
	}
      }
    }
    fout.write((char *)&var_el[0],nel*sizeof(double));
  }
  
  // ********** End printing Regions **********
  
  // ********** Start printing Solution **********
  if (_ml_sol != NULL)  {
  
  bool printAll = 0;
  for (unsigned ivar=0; ivar < vars.size(); ivar++){
    printAll += !(vars[ivar].compare("All")) + !(vars[ivar].compare("all")) + !(vars[ivar].compare("ALL"));
  }
   
  for (unsigned ivar=0; ivar< !printAll*vars.size() + printAll*_ml_sol->GetSolutionSize(); ivar++) {
    unsigned i = ( printAll == 0 ) ? _ml_sol->GetIndex( vars[ivar].c_str()) : ivar;
  
    for(int name=0;name<4;name++){
      if (name==0){
	sprintf(det,"%s", _ml_sol->GetSolutionName(i));
      }
      else if (name==1){
	sprintf(det,"%s %s","Bdc",_ml_sol->GetSolutionName(i));
      }
      else if (name==2){
	sprintf(det,"%s %s","Res",_ml_sol->GetSolutionName(i));
      }
      else{
	sprintf(det,"%s %s","Eps",_ml_sol->GetSolutionName(i));
      }
      if(name==0 || ( _debugOutput  && _ml_sol->GetSolutionLevel(igridn-1u)->_ResEpsBdcFlag[i])){
	if (_ml_sol->GetSolutionType(i)<3) {  // **********  on the nodes **********
	  fout.write((char *)det,sizeof(char)*8);
	  fout.write((char *)&one,sizeof(unsigned));
	  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
	    if (name==0){
	      Mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Sol[i],
				     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index, _ml_sol->GetSolutionType(i)) );
	    }
	    else if (name==1){
	      Mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Bdc[i],
				     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index, _ml_sol->GetSolutionType(i)) );
	    }
	    else if (name==2){
	      Mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Res[i],
				     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index, _ml_sol->GetSolutionType(i)) );
	    }
	    else{
	      Mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Eps[i],
				     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index, _ml_sol->GetSolutionType(i)) );
	    }
	    std::vector<double> v_local;
	    Mysol[ig]->localize_to_one(v_local,0);
	    fout.write((char *)&v_local[0],v_local.size()*sizeof(double));
	  }
	}
	else { // ********** on the elements **********
	  fout.write((char *)det,sizeof(char)*8);
	  fout.write((char *)&zero,sizeof(unsigned));
	  int icount=0;
	  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
	    std::vector<double> v_local;
	    if (name==0){
	      _ml_sol->GetSolutionLevel(ig)->_Sol[i]->localize_to_one(v_local,0); 
	    }
	    else if (name==1){
	      _ml_sol->GetSolutionLevel(ig)->_Bdc[i]->localize_to_one(v_local,0); 
	    }
	    else if (name==2){
	      _ml_sol->GetSolutionLevel(ig)->_Res[i]->localize_to_one(v_local,0);
	    }
	    else{
	      _ml_sol->GetSolutionLevel(ig)->_Eps[i]->localize_to_one(v_local,0);
	    }
	    for (unsigned ii=0; ii<_ml_mesh->GetLevel(ig)->GetNumberOfElements(); ii++) {
	      if ( ig==igridn-1u || 0==_ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
		unsigned iel_Metis = _ml_mesh->GetLevel(ig)->GetMetisDof(ii,_ml_sol->GetSolutionType(i));
		var_el[icount]=v_local[iel_Metis];
		icount++;
	      }
	    }
	  }
	  fout.write((char *)&var_el[0],nel*sizeof(double));
	}	
      }
    }
  }
  
  } //end _ml_sol
  // ********** End printing Solution **********
  sprintf(det,"%s","endvars");
  fout.write((char *)det,sizeof(char)*8);
  
  // ********** End printing Variables **********
  sprintf(det,"%s","endgmv");
  fout.write((char *)det,sizeof(char)*8);
  fout.close();
  // ********** End printing file **********
  
  // Free memory
  delete [] var_el;
  delete [] var_nd;
  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
    delete Mysol[ig];
  }
  delete [] det;
  
  return;   
}


void GMVWriter::Pwrite(const std::string output_path, const char order[], const std::vector<std::string>& vars, const unsigned time_step) const { 
  
  if(_nprocs == 1){
    write( output_path, order,  vars, time_step);
    return;
  }
  
  unsigned gridn = _gridn; // aggiunta da me
      
  if ( gridn == 0 ) gridn = _gridn;
  
  unsigned igridr=( _gridr <= gridn ) ? _gridr : gridn;

  // ********** linear -> index==0 *** quadratic -> index==1 **********
  unsigned index=( strcmp(order, "linear") ) ? 1 : 0;

  std::string dirnamePGMV = "GMVParallelFiles/";
  Files files; 
  files.CheckDir(output_path,dirnamePGMV);
  
  std::string filename_prefix;
  if( _ml_sol != NULL ) filename_prefix = "sol";
  else filename_prefix = "mesh";
  
  std::ostringstream filename;
  filename << output_path << "/" << dirnamePGMV << filename_prefix << ".level" << _gridn << "." <<_iproc<<"."<< time_step << "." << order << ".gmv"; 
    
  std::ofstream fout;
  
  fout.open(filename.str().c_str());
  if (!fout.is_open()) {
    std::cout << std::endl << " The output file "<< filename.str() <<" cannot be opened.\n";
    abort();
  }
  
  //count the own node dofs on all levels
  unsigned nvt = 0;
  unsigned nvt_max = 0;
  for (unsigned ig=igridr-1u; ig<gridn; ig++) {
    unsigned nvt_ig = _ml_mesh->GetLevel(ig)->own_size[index][_iproc];
    nvt_max=( nvt_max > nvt_ig ) ? nvt_max : nvt_ig;
    nvt += nvt_ig;
  }
  
   map < unsigned, unsigned > ghostMap; 
   
  // count the ghost node dofs and the own element dofs element on all levels
  unsigned ghostMapCounter = 0;
  unsigned gridOffset = 0;
  
  // count the ghost node dofs and the own element dofs element on all levels
  //unsigned ghost_counter = 0;
  unsigned nel=0;
  for (unsigned ig = igridr-1u; ig<gridn; ig++) {
    unsigned offset_iprc = _ml_mesh->GetLevel(ig)->MetisOffset[index][_iproc];	
    for (int iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
      unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
      if ( ig == gridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
	nel++;
	short unsigned ielt=_ml_mesh->GetLevel(ig)->el->GetElementType(kel);
        for(unsigned j=0; j<NVE[ielt][index]; j++){
	  unsigned jnode =_ml_mesh->GetLevel(ig)->el->GetMeshDof(kel, j, index); 
	  unsigned jnodeMetis = _ml_mesh->GetLevel(ig)->GetMetisDof(jnode, index);
	  if( jnodeMetis < offset_iprc ){ //Is this a ghost node?
	    if( ghostMap.find( gridOffset + jnodeMetis) == ghostMap.end()){
	      ghostMap[ gridOffset + jnodeMetis] = ghostMapCounter;
	      ghostMapCounter++;
	    }
	  }
	}
      }
    }
    gridOffset += _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
  }
  
  nvt_max= ( nvt_max > ghostMap.size() )? nvt_max : ghostMap.size();
  unsigned nvt0 = nvt;
  nvt += ghostMap.size(); // total node dofs (own + ghost)
  
  vector < double > var_nd;
  var_nd.resize(nvt_max+1); //TO FIX Valgrind complaints! In reality it should be only nvt
  vector < double > var_el;
  var_el.resize(nel+1); //TO FIX Valgrind complaints! In reality it should be only nel
  vector < NumericVector* > Mysol(gridn);
  for(unsigned ig=igridr-1u; ig<_gridn; ig++) {
    Mysol[ig] = NumericVector::build().release();
    
    if(n_processors()==1) { // IF SERIAL
      Mysol[ig]->init(_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs],
		      _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs],false,SERIAL);
    } 
    else{ // IF PARALLEL
      Mysol[ig]->init(_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs], 
		      _ml_mesh->GetLevel(ig)->own_size[index][_iproc],
		      _ml_mesh->GetLevel(ig)->ghost_nd_mts[index][_iproc],false,GHOSTED );
    }
  }
     
  // ********** Header **********
  char *det= new char[10];
  sprintf(det,"%s","gmvinput");
  fout.write((char *)det,sizeof(char)*8);
  sprintf(det,"%s","ieeei4r8");
  fout.write((char *)det,sizeof(char)*8);

  // ********** Start printing node coordinates  **********
  sprintf(det,"%s","nodes");
  fout.write((char *)det,sizeof(char)*8);
  fout.write((char *)&nvt,sizeof(unsigned));
    
  for (int i=0; i<3; i++) {
    for (unsigned ig=igridr-1u; ig<gridn; ig++) {
      Mysol[ig]->matrix_mult(*_ml_mesh->GetLevel(ig)->_coordinate->_Sol[i],
			     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,2) );
      
      unsigned offset_iprc = _ml_mesh->GetLevel(ig)->MetisOffset[index][_iproc];
      unsigned nvt_ig= _ml_mesh->GetLevel(ig)->own_size[index][_iproc];
      for (unsigned ii=0; ii<nvt_ig; ii++) 
	var_nd[ii]= (*Mysol[ig])(ii + offset_iprc);
      if (_ml_sol != NULL && _moving_mesh  && _ml_mesh->GetLevel(0)->GetDimension() > i)  { // if moving mesh
	unsigned indDXDYDZ=_ml_sol->GetIndex(_moving_vars[i].c_str());
	Mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Sol[indDXDYDZ],
			       *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,_ml_sol->GetSolutionType(indDXDYDZ)) );
 	for (unsigned ii=0; ii<nvt_ig; ii++) 
 	  var_nd[ii]+= (*Mysol[ig])(ii + offset_iprc);
      }
      fout.write( (char *)&var_nd[0], nvt_ig*sizeof(double) ); 
    }
    //print ghost coordinates
    
    gridOffset = 0;
    unsigned ig = igridr-1u;
    for (std::map <unsigned, unsigned>::iterator it=ghostMap.begin(); it!=ghostMap.end(); ++it){
      while( it->first >= gridOffset + _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs] ) {
	gridOffset += _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
	ig++;
      }
      var_nd[ it->second ] = (*Mysol[ig])( it->first - gridOffset);
    }
    if (_ml_sol != NULL && _moving_mesh  && _ml_mesh->GetLevel(0)->GetDimension() > i) { // if moving mesh  
      for (unsigned ig=igridr-1u; ig<gridn; ig++) {
	Mysol[ig]->matrix_mult(*_ml_mesh->GetLevel(ig)->_coordinate->_Sol[i],
			       *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,2) );
      }
      gridOffset = 0;
      unsigned ig = igridr-1u;
      for (std::map <unsigned, unsigned>::iterator it=ghostMap.begin(); it!=ghostMap.end(); ++it){
	while( it->first >= gridOffset + _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs] ) {
	  gridOffset += _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
	  ig++;
	}
	var_nd[ it->second ] += (*Mysol[ig])( it->first - gridOffset);
      }
    }
    fout.write( (char *)&var_nd[0], ghostMap.size()*sizeof(double) );    
  }
  // ********** End printing node coordinates  **********

  // ********** Start printing cell connectivity  **********
  const int eltp[2][6]= {{8,4,6,4,3,2},{20,10,15,8,6,3}};
  sprintf(det,"%s","cells");
  fout.write((char *)det,sizeof(char)*8);
  fout.write((char *)&nel,sizeof(unsigned));
  
  unsigned topology[27];
  unsigned offset=1;
  
  gridOffset = 0;
  for (unsigned ig=igridr-1u; ig<gridn; ig++) {
    unsigned offset_iprc = _ml_mesh->GetLevel(ig)->MetisOffset[index][_iproc];	
    unsigned nvt_ig= _ml_mesh->GetLevel(ig)->own_size[index][_iproc];
    for (int iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
      unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
      if ( ig == gridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
        short unsigned ielt=_ml_mesh->GetLevel(ig)->el->GetElementType(kel);
        if (ielt==0) sprintf(det,"phex%d",eltp[index][0]);
        else if (ielt==1) sprintf(det,"ptet%d",eltp[index][1]);
        else if (ielt==2) sprintf(det,"pprism%d",eltp[index][2]);
        else if (ielt==3) {
          if (eltp[index][3]==8) sprintf(det,"%dquad",eltp[index][3]);
          else sprintf(det,"quad");
        } 
        else if (ielt==4) {
          if (eltp[index][4]==6) sprintf(det,"%dtri",eltp[index][4]);
          else sprintf(det,"tri");
        } 
        else if (ielt==5) {
          if (eltp[index][5]==3) sprintf(det,"%dline",eltp[index][5]);
          else sprintf(det,"line");
        }
        fout.write((char *)det,sizeof(char)*8);
        fout.write((char *)&NVE[ielt][index],sizeof(unsigned));
	for(unsigned j=0;j<NVE[ielt][index];j++){
	  
	  unsigned jnode = _ml_mesh->GetLevel(ig)->el->GetMeshDof(kel, j, index);
	  unsigned jnodeMetis = _ml_mesh->GetLevel(ig)->GetMetisDof(jnode,index);
	  topology[j]=(jnodeMetis >= offset_iprc )? jnodeMetis - offset_iprc + offset : 
						     nvt0 + ghostMap[gridOffset+jnodeMetis] + 1u;
	}
	fout.write((char *)topology,sizeof(unsigned)*NVE[ielt][index]);
      }
    }
    offset += nvt_ig;
    gridOffset += _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
  }
  // ********** End printing cell connectivity  **********
     
  // ********** Start printing Variables **********
  const unsigned zero=0u;
  const unsigned one=1u;
  sprintf(det,"%s","variable");
  fout.write((char *)det,sizeof(char)*8);

  // ********** Start printing Regions **********
  strcpy(det,"Regions");
  fout.write((char *)det,sizeof(char)*8);
  fout.write((char *)&zero,sizeof(unsigned));

  int icount=0;
  for (unsigned ig=igridr-1u; ig<gridn; ig++) {
    for (int iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
      unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
      if ( ig==gridn-1u || 0==_ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
	var_el[icount]=_ml_mesh->GetLevel(ig)->el->GetElementGroup(kel);
        icount++;
      }
    }
  }
  fout.write((char *)&var_el[0],nel*sizeof(double));
  
  if(_nprocs>=1){
    strcpy(det,"Reg_proc");
    fout.write((char *)det,sizeof(char)*8);
    fout.write((char *)&zero,sizeof(unsigned));

    int icount=0;
    for (unsigned ig=igridr-1u; ig<gridn; ig++) {
      for (int iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
      unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
      if ( ig==gridn-1u || 0==_ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
	  var_el[icount]=_ml_mesh->GetLevel(ig)->epart[kel];
	  icount++;
	}
      }
    }
    fout.write((char *)&var_el[0],nel*sizeof(double));
  }
  
  // ********** End printing Regions **********
  
  // ********** Start printing Solution **********
  if (_ml_sol != NULL)  {
  
  bool printAll = 0;
  for (unsigned ivar=0; ivar < vars.size(); ivar++){
    printAll += !(vars[ivar].compare("All")) + !(vars[ivar].compare("all")) + !(vars[ivar].compare("ALL"));
  }
   
  for (unsigned ivar=0; ivar< !printAll*vars.size() + printAll*_ml_sol->GetSolutionSize(); ivar++) {
    unsigned i = ( printAll == 0 ) ? _ml_sol->GetIndex( vars[ivar].c_str()) : ivar;
  
    for(int name=0;name<4;name++){
      if (name==0){
	sprintf(det,"%s", _ml_sol->GetSolutionName(i));
      }
      else if (name==1){
	sprintf(det,"%s %s","Bdc",_ml_sol->GetSolutionName(i));
      }
      else if (name==2){
	sprintf(det,"%s %s","Res",_ml_sol->GetSolutionName(i));
      }
      else{
	sprintf(det,"%s %s","Eps",_ml_sol->GetSolutionName(i));
      }
      if(name==0 || ( _debugOutput  && _ml_sol->GetSolutionLevel(gridn-1u)->_ResEpsBdcFlag[i])){
	if (_ml_sol->GetSolutionType(i)<3) {  // **********  on the nodes **********
	  fout.write((char *)det,sizeof(char)*8);
	  fout.write((char *)&one,sizeof(unsigned));
	  for (unsigned ig=igridr-1u; ig<gridn; ig++) {
	    if (name==0){
	      Mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Sol[i],
				     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index, _ml_sol->GetSolutionType(i)) );
	    }
	    else if (name==1){
	      Mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Bdc[i],
				     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index, _ml_sol->GetSolutionType(i)) );
	    }
	    else if (name==2){
	      Mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Res[i],
				     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index, _ml_sol->GetSolutionType(i)) );
	    }
	    else{
	      Mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Eps[i],
				     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index, _ml_sol->GetSolutionType(i)) );
	    }
	    unsigned offset_iprc = _ml_mesh->GetLevel(ig)->MetisOffset[index][_iproc];
	    unsigned nvt_ig = _ml_mesh->GetLevel(ig)->own_size[index][_iproc];
	    for (unsigned ii=0; ii<nvt_ig; ii++) 
	      var_nd[ii]= (*Mysol[ig])(ii + offset_iprc);
	    fout.write((char *)&var_nd[0],nvt_ig*sizeof(double));
	  }
	  //print ghost dofs
	  gridOffset = 0;
	  unsigned ig = igridr-1u;
	  for (std::map <unsigned, unsigned>::iterator it=ghostMap.begin(); it!=ghostMap.end(); ++it){
	    while( it->first >= gridOffset + _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs] ) {
	      gridOffset += _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
	      ig++;
	    }
	    var_nd[ it->second ] = (*Mysol[ig])( it->first - gridOffset);
	  }
	  fout.write( (char *)&var_nd[0], ghostMap.size()*sizeof(double) );   	  
	}
	else { // ********** on the elements **********
	  fout.write((char *)det,sizeof(char)*8);
	  fout.write((char *)&zero,sizeof(unsigned));
	  int icount=0;
	  for (unsigned ig=igridr-1u; ig<gridn; ig++) {
	    for (unsigned iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
	      unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
		if ( ig == _gridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
		unsigned iel_Metis = _ml_mesh->GetLevel(ig)->GetMetisDof(kel,_ml_sol->GetSolutionType(i));
		if ( ig==gridn-1u || 0==_ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
		  if (name==0){
		    var_el[icount] = (*_ml_sol->GetSolutionLevel(ig)->_Sol[i])(iel_Metis);
		  }
		  else if (name==1){
		    var_el[icount] = (*_ml_sol->GetSolutionLevel(ig)->_Bdc[i])(iel_Metis);
		  }
		  else if (name==2){
		    var_el[icount] = (*_ml_sol->GetSolutionLevel(ig)->_Res[i])(iel_Metis);
		  }
		  else{
		    var_el[icount] = (*_ml_sol->GetSolutionLevel(ig)->_Eps[i])(iel_Metis);
		  }
		  icount++;
		}
	      }
	    }
	  }
	  fout.write((char *)&var_el[0],nel*sizeof(double));
	}	
      }
    }
  }
  
  } //end _ml_sol
  // ********** End printing Solution **********
  sprintf(det,"%s","endvars");
  fout.write((char *)det,sizeof(char)*8);
  
  // ********** End printing Variables **********
  sprintf(det,"%s","endgmv");
  fout.write((char *)det,sizeof(char)*8);
  fout.close();
  // ********** End printing file **********
  
  // Free memory
  for (unsigned ig=igridr-1u; ig<gridn; ig++) {
    delete Mysol[ig];
  }
  delete [] det;
  
  // print xml wrapper
  
  // *********** open pvtu file *************
  std::ofstream Pfout;
  if(_iproc!=0) {
    Pfout.rdbuf();   //redirect to dev_null
  }
  else {
    std::ostringstream Pfilename;
    Pfilename << output_path << "/" << filename_prefix << "_level" << _gridn <<"_"<< time_step << "_" << order << "_gmv.xml";
    Pfout.open(Pfilename.str().c_str());
    if (Pfout.is_open()) {
      std::cout << std::endl << " The output is printed to file " << Pfilename.str() << " in parallel GMV-XML (64-based) format" << std::endl; 
      std::cout << " Laod it in ParaView -> Tools -> Play Test ... \n";
    }
    else {
      std::cout << std::endl << " The output file "<< Pfilename.str() <<" cannot be opened.\n";
      abort();
    }
  }
  
  // *********** write pvtu header ***********
  Pfout<< "<?xml version=\"1.0\"?>" << std::endl;
  Pfout<< "<pqevents>\n";
  for(int jproc=0;jproc<_nprocs;jproc++){
    Pfout<< "  <pqevent object=\"pqClientMainWindow/MainControlsToolbar/1QToolButton0\" command=\"activate\" arguments=\"\" />\n";
    Pfout<< "  <pqevent object=\"pqClientMainWindow/FileOpenDialog\" command=\"filesSelected\" arguments=\"";
    Pfout<< dirnamePGMV << filename_prefix << ".level" << _gridn << "." <<jproc<<"."<< time_step << "." << order << ".gmv\" />\n";
  }
  
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mousePress\"   arguments=\"1,1,0,50,5,/0:0\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mouseRelease\" arguments=\"1,0,0,50,5,/0:0\" />\n";
  
  Pfout<< "  <pqevent object=\"pqClientMainWindow/propertiesDock/propertiesPanel/Accept\" command=\"activate\"   arguments=\"\" />\n";
    
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mousePress\"   arguments=\"1,1,33554432,50,5,/0:0/"<<_nprocs-1<<":0\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mouseRelease\" arguments=\"1,0,33554432,50,5,/0:0/"<<_nprocs-1<<":0\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"keyEvent\" 	 arguments=\"7,16777248,0,,0,1\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/menubar\" command=\"activate\" arguments=\"menuFilters\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/menubar/menuFilters/Alphabetical\" command=\"activate\"        arguments=\"GroupDataSets\" />\n";
  
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mousePress\"   arguments=\"1,1,0,10,5,/0:0/"<<_nprocs<<":1\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mouseRelease\" arguments=\"1,0,0,10,5,/0:0/"<<_nprocs<<":1\" />\n";
  
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mousePress\"   arguments=\"1,1,33554432,25,5,/0:0\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mouseRelease\" arguments=\"1,0,33554432,25,5,/0:0\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mousePress\"   arguments=\"1,1,33554432,50,5,/0:0/"<<_nprocs<<":0\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mouseRelease\" arguments=\"1,0,33554432,50,5,/0:0/"<<_nprocs<<":0\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"keyEvent\"     arguments=\"7,16777248,0,,0,1\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mousePress\"   arguments=\"1,1,0,10,5,/0:0/"<<_nprocs<<":1\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mouseRelease\" arguments=\"1,0,0,10,5,/0:0/"<<_nprocs<<":1\" />\n";
  
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mousePress\"   arguments=\"1,1,0,50,5,/0:0/"<<_nprocs<<":1\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mouseRelease\" arguments=\"1,0,0,50,5,/0:0/"<<_nprocs<<":1\" />\n";
  
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mousePress\"   arguments=\"1,1,0,10,5,/0:0/"<<_nprocs<<":1\" />\n";
  Pfout<< "  <pqevent object=\"pqClientMainWindow/pipelineBrowserDock/pipelineBrowser\" command=\"mouseRelease\" arguments=\"1,0,0,10,5,/0:0/"<<_nprocs<<":1\" />\n";
  
  Pfout<< "</pqevents>\n";
  Pfout.close();
    
  
  
  return;   
}


} //end namespace femus




