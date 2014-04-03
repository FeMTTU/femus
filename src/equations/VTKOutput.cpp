/*=========================================================================

 Program: FEMUS
 Module: VTKOutput
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
#include "VTKOutput.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include <b64/b64.h>
#include "iostream"
#include <fstream>
#include "string.h"
#include "stdio.h"
#include <iomanip>
#include <algorithm>   


VTKOutput::VTKOutput(MultiLevelProblem& ml_probl): Output(ml_probl)
{
  
}

VTKOutput::~VTKOutput()
{
  
}


void VTKOutput::write_system_solutions(const char order[], std::vector<std::string>& vars, const unsigned time_step) 
{ 
  bool test_all=!(vars[0].compare("All"));
  
  int icount;
  unsigned index=0;
  unsigned index_nd=0;
  if (!strcmp(order,"linear")) {   //linear
    index=0;
    index_nd=0;
  } else if (!strcmp(order,"quadratic")) { //quadratic
    index=1;
    index_nd=1;
  } else if (!strcmp(order,"biquadratic")) { //biquadratic
    index=3;
    index_nd=2;
  }

  const int eltp[4][6]= {{12,10,13,9,5,3},{25,24,26,23,22,21},{},{29,1,1,28,22,1}};
  
  char *filename= new char[60];
  sprintf(filename,"./output/mesh.level%d.%d.%s.vtu",_gridn,time_step,order);
  std::ofstream fout;
  
  if(_iproc!=0) {
    fout.rdbuf();   //redirect to dev_null
  }
  else {
    fout.open(filename);
    if (!fout) {
      std::cout << "Output mesh file "<<filename<<" cannot be opened.\n";
      exit(0);
    }
  }
  // haed ************************************************
  fout<<"<?xml version=\"1.0\"?>" << std::endl;
  fout<<"<VTKFile type = \"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  fout << " <UnstructuredGrid>" << std::endl;
  //-----------------------------------------------------------------------------------------------------------

  //----------------------------------------------------------------------------------------------------------
  
  unsigned nvt=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    unsigned nvt_ig=_ml_probl._ml_msh->_level[ig]->MetisOffset[index_nd][_nprocs];
    nvt+=nvt_ig;
  }  

  unsigned nel=0;
  unsigned counter=0;
  for (unsigned ig=_gridr-1u; ig<_gridn-1u; ig++) {
    nel+=( _ml_probl._ml_msh->_level[ig]->GetElementNumber() - _ml_probl._ml_msh->_level[ig]->el->GetRefinedElementNumber());
    counter+=(_ml_probl._ml_msh->_level[ig]->el->GetElementNumber("Hex")-_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementNumber("Hex"))*NVE[0][index];
    counter+=(_ml_probl._ml_msh->_level[ig]->el->GetElementNumber("Tet")-_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementNumber("Tet"))*NVE[1][index];
    counter+=(_ml_probl._ml_msh->_level[ig]->el->GetElementNumber("Wedge")-_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementNumber("Wedge"))*NVE[2][index];
    counter+=(_ml_probl._ml_msh->_level[ig]->el->GetElementNumber("Quad")-_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementNumber("Quad"))*NVE[3][index];
    counter+=(_ml_probl._ml_msh->_level[ig]->el->GetElementNumber("Triangle")-_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementNumber("Triangle"))*NVE[4][index];
    counter+=(_ml_probl._ml_msh->_level[ig]->el->GetElementNumber("Line")-_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementNumber("Line"))*NVE[5][index];
  }
  nel+=_ml_probl._ml_msh->_level[_gridn-1u]->GetElementNumber();
  counter+=_ml_probl._ml_msh->_level[_gridn-1u]->el->GetElementNumber("Hex")*NVE[0][index];
  counter+=_ml_probl._ml_msh->_level[_gridn-1u]->el->GetElementNumber("Tet")*NVE[1][index];
  counter+=_ml_probl._ml_msh->_level[_gridn-1u]->el->GetElementNumber("Wedge")*NVE[2][index];
  counter+=_ml_probl._ml_msh->_level[_gridn-1u]->el->GetElementNumber("Quad")*NVE[3][index];
  counter+=_ml_probl._ml_msh->_level[_gridn-1u]->el->GetElementNumber("Triangle")*NVE[4][index];
  counter+=_ml_probl._ml_msh->_level[_gridn-1u]->el->GetElementNumber("Line")*NVE[5][index]; 
  
  const unsigned dim_array_coord [] = { nvt*3*sizeof(float) };  
  const unsigned dim_array_conn[]   = { counter*sizeof(int) };
  const unsigned dim_array_off []   = { nel*sizeof(int) };
  const unsigned dim_array_type []  = { nel*sizeof(short unsigned) };
  const unsigned dim_array_reg []   = { nel*sizeof(short unsigned) };
  const unsigned dim_array_elvar [] = { nel*sizeof(float) };
  const unsigned dim_array_ndvar [] = { nvt*sizeof(float) };

  // initialize common buffer_void memory 
  unsigned buffer_size=(dim_array_coord[0]>dim_array_conn[0])?dim_array_coord[0]:dim_array_conn[0];
  void *buffer_void=new char [buffer_size];
  char *buffer_char=static_cast <char *>(buffer_void);
  
  size_t cch;
  cch = b64::b64_encode(&buffer_char[0], buffer_size , NULL, 0);  
  vector <char> enc;
  enc.resize(cch);
  char *pt_char;
  
  fout << "  <Piece NumberOfPoints= \"" << nvt << "\" NumberOfCells= \"" << nel << "\" >" << std::endl;
  
  //-----------------------------------------------------------------------------------------------
  // print coordinates *********************************************Solu*******************************************
  fout << "   <Points>" << std::endl;
  fout << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;
  
//   unsigned indXYZ[3];
//   indXYZ[0]=_ml_probl.GetIndex("X");
//   indXYZ[1]=_ml_probl.GetIndex("Y");
//   indXYZ[2]=_ml_probl.GetIndex("Z");
  
  vector <NumericVector*> mysol(_gridn);
  for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    mysol[ig] = NumericVector::build().release();
    mysol[ig]->init(_ml_probl._ml_msh->_level[ig]->MetisOffset[index_nd][_nprocs],_ml_probl._ml_msh->_level[ig]->own_size[index_nd][_iproc],true,AUTOMATIC);
  }
  
  // point pointer to common mamory area buffer of void type;
  float *var_coord= static_cast<float*>(buffer_void);

  unsigned offset_nvt3=0;
  for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    std::vector<double> v_local;
    unsigned nvt_ig=_ml_probl._ml_msh->_level[ig]->MetisOffset[index_nd][_nprocs];
    for(int kk=0;kk<3;kk++) {
      mysol[ig]->matrix_mult(*_ml_probl._ml_msh->_level[ig]->_coordinate->_Sol[kk],*_ProlQitoQj[index_nd][2][ig]);
      mysol[ig]->localize_to_one(v_local,0);
      if(_iproc==0) { 
	for (unsigned i=0; i<nvt_ig; i++) {
	  var_coord[offset_nvt3+i*3+kk] = v_local[i];
	}
      } //if iproc
    } //loop over dimension
    offset_nvt3+=3*nvt_ig;
  }
  
  if (_moving_mesh) {
    
    unsigned offset_nvt3=0;
    unsigned indDXDYDZ[3];
    indDXDYDZ[0]=_ml_probl.GetIndex(_moving_vars[0].c_str());
    indDXDYDZ[1]=_ml_probl.GetIndex(_moving_vars[1].c_str());
    if(_ml_probl._ml_msh->_level[0]->GetDimension() == 3) {
      indDXDYDZ[2]=_ml_probl.GetIndex(_moving_vars[2].c_str());
    }
      
    for(unsigned ig=_gridr-1u; ig<_gridn; ig++){
      std::vector<double> v_local;
      unsigned nvt_ig=_ml_probl._ml_msh->_level[ig]->MetisOffset[index_nd][_nprocs];
      for(int kk=0;kk<_ml_probl._ml_msh->_level[0]->GetDimension();kk++) {
	mysol[ig]->matrix_mult(*_ml_probl._solution[ig]->_Sol[indDXDYDZ[kk]],*_ProlQitoQj[index_nd][_ml_probl.SolType[indDXDYDZ[kk]]][ig]);
        mysol[ig]->localize_to_one(v_local,0);
	if(_iproc==0) { 
	  for (unsigned i=0; i<nvt_ig; i++) {
	    var_coord[offset_nvt3+i*3+kk] += v_local[i];
	  }
	} //if iproc
      } //loop over dimension
      offset_nvt3+=3*nvt_ig;
    }
  }
  
  if(_iproc==0) {
    //print coordinates dimension
    cch = b64::b64_encode(&dim_array_coord[0], sizeof(dim_array_coord), NULL, 0);  
    b64::b64_encode(&dim_array_coord[0], sizeof(dim_array_coord), &enc[0], cch);
    pt_char=&enc[0];
    for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
    //print coordinates array
    cch = b64::b64_encode(&var_coord[0], dim_array_coord[0] , NULL, 0);  
    b64::b64_encode(&var_coord[0], dim_array_coord[0], &enc[0], cch);
    pt_char=&enc[0];
    for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    fout << std::endl;
  }
  fout << "    </DataArray>" << std::endl;
  fout << "   </Points>" << std::endl;
  //-----------------------------------------------------------------------------------------------
  
  //-----------------------------------------------------------------------------------------------
  // Printing of element connectivity - offset - format type  *
  fout << "   <Cells>" << std::endl;
  
  //-----------------------------------------------------------------------------------------------
  //print connectivity
  fout << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">" << std::endl;
  
  // point pointer to common mamory area buffer of void type;
  int *var_conn = static_cast <int*> (buffer_void);
  icount = 0;
  unsigned offset_nvt=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned iel=0; iel<_ml_probl._ml_msh->_level[ig]->GetElementNumber(); iel++) {
      if (_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementIndex(iel)==0 || ig==_gridn-1u) {
        for (unsigned j=0; j<_ml_probl._ml_msh->_level[ig]->el->GetElementDofNumber(iel,index); j++) {
	  unsigned jnode=_ml_probl._ml_msh->_level[ig]->el->GetElementVertexIndex(iel,j)-1u;
	  unsigned jnode_Metis = _ml_probl._ml_msh->_level[ig]->GetMetisDof(jnode,index_nd);
	  var_conn[icount] = offset_nvt+jnode_Metis;
	  icount++;
	}
      }
    }
    offset_nvt+=_ml_probl._ml_msh->_level[ig]->MetisOffset[index_nd][_nprocs];
  }
  
  //print connectivity dimension
  cch = b64::b64_encode(&dim_array_conn[0], sizeof(dim_array_conn), NULL, 0);  
  b64::b64_encode(&dim_array_conn[0], sizeof(dim_array_conn), &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  //print connectivity array
  cch = b64::b64_encode(&var_conn[0], dim_array_conn[0] , NULL, 0);  
  b64::b64_encode(&var_conn[0], dim_array_conn[0], &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  fout << std::endl;
  
  fout << "    </DataArray>" << std::endl;
  //------------------------------------------------------------------------------------------------
  
  //-------------------------------------------------------------------------------------------------
  //printing offset
  fout << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">" << std::endl;
 
  // point pointer to common mamory area buffer of void type;
  int *var_off=static_cast <int*>(buffer_void);
  icount = 0;
  int offset_el=0;
  //print offset array
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned iel=0; iel<_ml_probl._ml_msh->_level[ig]->GetElementNumber(); iel++) {
      if (_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementIndex(iel)==0 || ig==_gridn-1u) {
	unsigned iel_Metis = _ml_probl._ml_msh->_level[ig]->GetMetisDof(iel,3);
        offset_el += _ml_probl._ml_msh->_level[ig]->el->GetElementDofNumber(iel_Metis,index);
        var_off[icount] = offset_el;
	icount++;
      }
    }
  }
  
  //print offset dimension
  cch = b64::b64_encode(&dim_array_off[0], sizeof(dim_array_off), NULL, 0);
  b64::b64_encode(&dim_array_off[0], sizeof(dim_array_off), &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  //print offset array
  cch = b64::b64_encode(&var_off[0], dim_array_off[0] , NULL, 0);  
  b64::b64_encode(&var_off[0], dim_array_off[0], &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  fout << std::endl;
  
  fout << "    </DataArray>" << std::endl;
  //--------------------------------------------------------------------------------------------------
  
  //--------------------------------------------------------------------------------------------------
  //Element format type : 23:Serendipity(8-nodes)  28:Quad9-Biquadratic
  fout << "    <DataArray type=\"UInt16\" Name=\"types\" format=\"binary\">" << std::endl;
   
  // point pointer to common mamory area buffer of void type;
  unsigned short *var_type = static_cast <unsigned short*> (buffer_void);
  icount=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned ii=0; ii<_ml_probl._ml_msh->_level[ig]->GetElementNumber(); ii++) {
      if (_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementIndex(ii)==0 || ig==_gridn-1u) {
	unsigned iel_Metis = _ml_probl._ml_msh->_level[ig]->GetMetisDof(ii,3);
        short unsigned ielt= _ml_probl._ml_msh->_level[ig]->el->GetElementType(iel_Metis);
	var_type[icount] = (short unsigned)(eltp[index][ielt]);
	icount++;
      }
    }
  }
  
  //print element format dimension
  cch = b64::b64_encode(&dim_array_type[0], sizeof(dim_array_type), NULL, 0);  
  b64::b64_encode(&dim_array_type[0], sizeof(dim_array_type), &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  //print element format array
  cch = b64::b64_encode(&var_type[0], dim_array_type[0] , NULL, 0);  
  b64::b64_encode(&var_type[0], dim_array_type[0], &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  fout << std::endl;
  fout << "    </DataArray>" << std::endl;
  //----------------------------------------------------------------------------------------------------
  
  fout << "   </Cells>" << std::endl;
  //--------------------------------------------------------------------------------------------------

  // /Print Cell Data ****************************************************************************
  fout << "   <CellData Scalars=\"scalars\">" << std::endl;

  //--------------------------------------------------------------------------------------------
  // Print Regions
  fout << "    <DataArray type=\"UInt16\" Name=\"Regions\" format=\"binary\">" << std::endl;
  
  // point pointer to common mamory area buffer of void type;
  unsigned short* var_reg=static_cast <unsigned short*> (buffer_void);
  icount=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned ii=0; ii<_ml_probl._ml_msh->_level[ig]->GetElementNumber(); ii++) {
      if ( ig==_gridn-1u || 0==_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementIndex(ii)) {
	unsigned iel_Metis = _ml_probl._ml_msh->_level[ig]->GetMetisDof(ii,3);
	var_reg[icount]= _ml_probl._ml_msh->_level[ig]->el->GetElementGroup(ii);
	icount++;
      }
    }
  }
  
  //print regions dimension
  cch = b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), NULL, 0);  
  b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  
  //print regions array
  cch = b64::b64_encode(&var_reg[0], dim_array_reg[0] , NULL, 0);  
  b64::b64_encode(&var_reg[0], dim_array_reg[0], &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  fout << std::endl;
  fout << "    </DataArray>" << std::endl;
  //-----------------------------------------------------------------------------------------------------   
  // Print Metis Partitioning
  fout << "    <DataArray type=\"UInt16\" Name=\"Domain_partition\" format=\"binary\">" << std::endl;
  
  // point pointer to common mamory area buffer of void type;
  unsigned short* var_proc=static_cast <unsigned short*> (buffer_void);
  icount=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned ii=0; ii<_ml_probl._ml_msh->_level[ig]->GetElementNumber(); ii++) {
      if ( ig==_gridn-1u || 0==_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementIndex(ii)) {
	unsigned iel_Metis = _ml_probl._ml_msh->_level[ig]->GetMetisDof(ii,3);
	var_proc[icount]=(unsigned short)(_ml_probl._ml_msh->_level[ig]->epart[ii]);
	icount++;
      }
    }
  }
  
  //print regions dimension
  cch = b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), NULL, 0);  
  b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  
  //print regions array
  cch = b64::b64_encode(&var_proc[0], dim_array_reg[0] , NULL, 0);  
  b64::b64_encode(&var_proc[0], dim_array_reg[0], &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  fout << std::endl;
  fout << "    </DataArray>" << std::endl;
  

  //Print Solution (on element) ***************************************************************
  for (unsigned i=0; i<(1-test_all)*vars.size()+test_all*_ml_probl.SolName.size(); i++) {
    unsigned indx=(test_all==0)?_ml_probl.GetIndex(vars[i].c_str()):i;
    if (3 <= _ml_probl.SolType[indx]) {
      fout << "    <DataArray type=\"Float32\" Name=\"" << _ml_probl.SolName[indx] <<"\" format=\"binary\">" << std::endl;
      // point pointer to common memory area buffer of void type;
      float *var_el = static_cast< float*> (buffer_void);
      icount=0;
      for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
	vector<double> sol_local;
	_ml_probl._solution[ig]->_Sol[indx]->localize_to_one(sol_local,0);
	for (unsigned ii=0; ii<_ml_probl._ml_msh->_level[ig]->GetElementNumber(); ii++) {
	  if (ig==_gridn-1u || 0==_ml_probl._ml_msh->_level[ig]->el->GetRefinedElementIndex(ii)) {
	    unsigned iel_Metis = _ml_probl._ml_msh->_level[ig]->GetMetisDof(ii,_ml_probl.SolType[indx]);
	    var_el[icount]=sol_local[iel_Metis];
	    icount++;
	  }
	}
      }
      
      if(_iproc==0) {
        //print solution on element dimension
        cch = b64::b64_encode(&dim_array_elvar[0], sizeof(dim_array_elvar), NULL, 0);  
        b64::b64_encode(&dim_array_elvar[0], sizeof(dim_array_elvar), &enc[0], cch);
        pt_char=&enc[0];
        for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
        //print solution on element array
        cch = b64::b64_encode(&var_el[0], dim_array_elvar[0] , NULL, 0);  
        b64::b64_encode(&var_el[0], dim_array_elvar[0], &enc[0], cch);
        pt_char=&enc[0];
        for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
        fout << std::endl;
        fout << "    </DataArray>" << std::endl;
      }
      //----------------------------------------------------------------------------------------------------
    }
  }
  fout << "   </CellData>" << std::endl;
  //   //------------------------------------------------------------------------------------------------
  // 
  //   //------------------------------------------------------------------------------------------------
  // / Print Solution (on nodes) ********************************************************************
  fout<< " <PointData Scalars=\"scalars\"> " << std::endl;
  //Loop on variables
   
  // point pointer to common mamory area buffer of void type;
  float* var_nd = static_cast<float*>(buffer_void);
  for (unsigned i=0; i<(1-test_all)*vars.size()+test_all*_ml_probl.SolName.size(); i++) {
    unsigned indx=(test_all==0)?_ml_probl.GetIndex(vars[i].c_str()):i;
    if (_ml_probl.SolType[indx]<3) {
      fout << " <DataArray type=\"Float32\" Name=\"" << _ml_probl.SolName[indx] <<"\" format=\"binary\">" << std::endl;
      //print solutions on nodes dimension
      cch = b64::b64_encode(&dim_array_ndvar[0], sizeof(dim_array_ndvar), NULL, 0);  
      b64::b64_encode(&dim_array_ndvar[0], sizeof(dim_array_ndvar), &enc[0], cch);
      pt_char=&enc[0];
      for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
      
      unsigned offset_nvt=0;
      for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
	mysol[ig]->matrix_mult(*_ml_probl._solution[ig]->_Sol[indx],*_ProlQitoQj[index_nd][_ml_probl.SolType[indx]][ig]);
	vector<double> sol_local;
	mysol[ig]->localize_to_one(sol_local,0);
	unsigned nvt_ig=_ml_probl._ml_msh->_level[ig]->MetisOffset[index_nd][_nprocs];
	for (unsigned ii=0; ii<nvt_ig; ii++) {
	  var_nd[ii+offset_nvt] = sol_local[ii];
	}
	offset_nvt+=nvt_ig;
      }
      
      if(_iproc==0) {
        cch = b64::b64_encode(&var_nd[0], dim_array_ndvar [0], NULL, 0);  
        b64::b64_encode(&var_nd[0], dim_array_ndvar [0], &enc[0], cch);
        pt_char=&enc[0];
        for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
        fout << std::endl;

        fout << "    </DataArray>" << std::endl;
      }
    } //endif
  } // end for sol
  fout << "   </PointData>" << std::endl;
  
  //------------------------------------------------------------------------------------------------
  
  fout << "  </Piece>" << std::endl;
  fout << " </UnstructuredGrid>" << std::endl;
  fout << "</VTKFile>" << std::endl;
  fout.close();
  
  //-----------------------------------------------------------------------------------------------------
  //free memory
  for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    delete mysol[ig];
  }
  delete [] filename;
  delete [] var_nd;  
  
  //--------------------------------------------------------------------------------------------------------
  return;   
}
  
  
  