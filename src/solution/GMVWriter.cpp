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


namespace femus {



GMVWriter::GMVWriter(MultiLevelSolution & ml_probl): Writer(ml_probl)
{
  _debugOutput = false;
}

GMVWriter::~GMVWriter()
{
  
}

void GMVWriter::write_system_solutions(const std::string output_path, const char order[], std::vector<std::string>& vars, const unsigned time_step) 
{ 
  unsigned igridn = _gridn; // aggiunta da me
      
  if (igridn==0) igridn=_gridn;
  
  unsigned igridr=(_gridr <= igridn)?_gridr:igridn;

  // ********** linear -> index==0 *** quadratic -> index==1 **********
  unsigned index=(strcmp(order,"linear"))?1:0;

  std::ostringstream filename;
  filename << output_path << "/sol.level" << _gridn << "." << time_step << "." << order << ".gmv"; 
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
    unsigned nvt_ig=_ml_sol._ml_msh->GetLevel(ig)->MetisOffset[index][_nprocs];
    nvt_max=(nvt_max>nvt_ig)?nvt_max:nvt_ig;
    nvt+=nvt_ig;
  }
  
  double *var_nd=new double [nvt_max+1]; //TO FIX Valgrind complaints! In reality it should be only nvt
  vector <NumericVector*> Mysol(igridn);
  for(unsigned ig=igridr-1u; ig<_gridn; ig++) {
    Mysol[ig] = NumericVector::build().release();
    Mysol[ig]->init(_ml_sol._ml_msh->GetLevel(ig)->MetisOffset[index][_nprocs],_ml_sol._ml_msh->GetLevel(ig)->own_size[index][_iproc],true,AUTOMATIC);
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
      Mysol[ig]->matrix_mult(*_ml_sol._ml_msh->GetLevel(ig)->_coordinate->_Sol[i],*_ProlQitoQj[index][2][ig]);
      vector <double> v_local;
      Mysol[ig]->localize_to_one(v_local,0);
      unsigned nvt_ig=_ml_sol._ml_msh->GetLevel(ig)->MetisOffset[index][_nprocs];      
      if(_iproc==0){ 
	for (unsigned ii=0; ii<nvt_ig; ii++) 
	  var_nd[ii]= v_local[ii];
      }
      if (_moving_mesh) {
	unsigned indDXDYDZ=_ml_sol.GetIndex(_moving_vars[i].c_str());
	Mysol[ig]->matrix_mult(*_ml_sol.GetSolutionLevel(ig)->_Sol[indDXDYDZ],*_ProlQitoQj[index][_ml_sol.GetSolutionType(indDXDYDZ)][ig]);
	Mysol[ig]->localize_to_one(v_local,0);
	unsigned nvt_ig=_ml_sol._ml_msh->GetLevel(ig)->MetisOffset[index][_nprocs];      
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
    nel+=( _ml_sol._ml_msh->GetLevel(ig)->GetNumberOfElements() - _ml_sol._ml_msh->GetLevel(ig)->el->GetRefinedElementNumber());
  nel+=_ml_sol._ml_msh->GetLevel(igridn-1u)->GetNumberOfElements();
  fout.write((char *)&nel,sizeof(unsigned));

  unsigned topology[27];
  unsigned offset=1;
  
  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
    for (unsigned ii=0; ii<_ml_sol._ml_msh->GetLevel(ig)->GetNumberOfElements(); ii++) {
      if ( ig==igridn-1u || 0==_ml_sol._ml_msh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
        short unsigned ielt=_ml_sol._ml_msh->GetLevel(ig)->el->GetElementType(ii);
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
	  
	  unsigned jnode=_ml_sol._ml_msh->GetLevel(ig)->el->GetElementVertexIndex(ii,j)-1u;
	  unsigned jnode_Metis = _ml_sol._ml_msh->GetLevel(ig)->GetMetisDof(jnode,index);
	  	  
	  topology[j]=jnode_Metis+offset;
	}
	fout.write((char *)topology,sizeof(unsigned)*NVE[ielt][index]);
      }
    }
    offset+=_ml_sol._ml_msh->GetLevel(ig)->MetisOffset[index][_nprocs];
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
    for (unsigned ii=0; ii<_ml_sol._ml_msh->GetLevel(ig)->GetNumberOfElements(); ii++) {
      if ( ig==igridn-1u || 0==_ml_sol._ml_msh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
	var_el[icount]=_ml_sol._ml_msh->GetLevel(ig)->el->GetElementGroup(ii);
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
      for (unsigned ii=0; ii<_ml_sol._ml_msh->GetLevel(ig)->GetNumberOfElements(); ii++) {
	if ( ig==igridn-1u || 0==_ml_sol._ml_msh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
	  var_el[icount]=_ml_sol._ml_msh->GetLevel(ig)->epart[ii];
	  icount++;
	}
      }
    }
    fout.write((char *)&var_el[0],nel*sizeof(double));
  }
  
  // ********** End printing Regions **********
  
  // ********** Start printing Solution **********
  bool printAll = 0;
  for (unsigned ivar=0; ivar < vars.size(); ivar++){
    printAll += !(vars[ivar].compare("All")) + !(vars[ivar].compare("all")) + !(vars[ivar].compare("ALL"));
  }
   
  for (unsigned ivar=0; ivar< !printAll*vars.size() + printAll*_ml_sol.GetSolutionSize(); ivar++) {
    unsigned i = ( printAll == 0 ) ? _ml_sol.GetIndex( vars[ivar].c_str()) : ivar;
  
    for(int name=0;name<4;name++){
      if (name==0){
	sprintf(det,"%s", _ml_sol.GetSolutionName(i));
      }
      else if (name==1){
	sprintf(det,"%s %s","Bdc",_ml_sol.GetSolutionName(i));
      }
      else if (name==2){
	sprintf(det,"%s %s","Res",_ml_sol.GetSolutionName(i));
      }
      else{
	sprintf(det,"%s %s","Eps",_ml_sol.GetSolutionName(i));
      }
      if(name==0 || ( _debugOutput  && _ml_sol.GetSolutionLevel(igridn-1u)->_ResEpsBdcFlag[i])){
	if (_ml_sol.GetSolutionType(i)<3) {  // **********  on the nodes **********
	  fout.write((char *)det,sizeof(char)*8);
	  fout.write((char *)&one,sizeof(unsigned));
	  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
	    if (name==0){
	      Mysol[ig]->matrix_mult(*_ml_sol.GetSolutionLevel(ig)->_Sol[i],*_ProlQitoQj[index][_ml_sol.GetSolutionType(i)][ig]);  
	    }
	    else if (name==1){
	      Mysol[ig]->matrix_mult(*_ml_sol.GetSolutionLevel(ig)->_Bdc[i],*_ProlQitoQj[index][_ml_sol.GetSolutionType(i)][ig]);
	    }
	    else if (name==2){
	      Mysol[ig]->matrix_mult(*_ml_sol.GetSolutionLevel(ig)->_Res[i],*_ProlQitoQj[index][_ml_sol.GetSolutionType(i)][ig]);
	    }
	    else{
	      Mysol[ig]->matrix_mult(*_ml_sol.GetSolutionLevel(ig)->_Eps[i],*_ProlQitoQj[index][_ml_sol.GetSolutionType(i)][ig]);
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
	      _ml_sol.GetSolutionLevel(ig)->_Sol[i]->localize_to_one(v_local,0); 
	    }
	    else if (name==1){
	      _ml_sol.GetSolutionLevel(ig)->_Bdc[i]->localize_to_one(v_local,0); 
	    }
	    else if (name==2){
	      _ml_sol.GetSolutionLevel(ig)->_Res[i]->localize_to_one(v_local,0);
	    }
	    else{
	      _ml_sol.GetSolutionLevel(ig)->_Eps[i]->localize_to_one(v_local,0);
	    }
	    for (unsigned ii=0; ii<_ml_sol._ml_msh->GetLevel(ig)->GetNumberOfElements(); ii++) {
	      if ( ig==igridn-1u || 0==_ml_sol._ml_msh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
		unsigned iel_Metis = _ml_sol._ml_msh->GetLevel(ig)->GetMetisDof(ii,_ml_sol.GetSolutionType(i));
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


} //end namespace femus


