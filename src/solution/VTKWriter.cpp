/*=========================================================================

 Program: FEMUS
 Module: VTKWriter
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
#include "VTKWriter.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include <b64/b64.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <iomanip>
#include <algorithm>
#include "Files.hpp"

namespace femus {


 short unsigned int VTKWriter::femusToVtkCellType[3][6]= {{12,10,13,9,5,3},{25,24,26,23,22,21},{29,24,32,28,22,21}};


VTKWriter::VTKWriter(MultiLevelSolution * ml_sol): Writer(ml_sol) {}

VTKWriter::VTKWriter(MultiLevelMesh * ml_mesh): Writer(ml_mesh) {}

VTKWriter::~VTKWriter() {}


void VTKWriter::Pwrite(const std::string output_path, const char order[], const std::vector < std::string > & vars, const unsigned time_step) const {

  // *********** open vtu files *************
  std::ofstream fout;

  std::string dirnamePVTK = "VTKParallelFiles/";
  Files files;
  files.CheckDir(output_path,dirnamePVTK);

  int icount;
  unsigned index=0;
  if 	  (!strcmp(order,"linear")) 	 index=0; //linear
  else if (!strcmp(order,"quadratic")) 	 index=1; //quadratic
  else if (!strcmp(order,"biquadratic")) index=2; //biquadratic

  std::string filename_prefix;
  if( _ml_sol != NULL ) filename_prefix = "sol";
  else filename_prefix = "mesh";

  std::ostringstream filename;
  filename << output_path << "/"<<dirnamePVTK<< filename_prefix << ".level" << _gridn << "." <<_iproc<<"."<< time_step << "." << order << ".vtu";

  fout.open(filename.str().c_str());
  if (!fout.is_open()) {
    std::cout << std::endl << " The output file "<< filename.str() <<" cannot be opened.\n";
    abort();
  }

  // *********** write vtu header ************
  fout<<"<?xml version=\"1.0\"?>" << std::endl;
  fout<<"<VTKFile type = \"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  fout << "  <UnstructuredGrid>" << std::endl;

  // *********** open pvtu file *************
  std::ofstream Pfout;
  if(_iproc!=0) {
    Pfout.rdbuf();   //redirect to dev_null
  }
  else {
    std::ostringstream Pfilename;
    Pfilename << output_path << "/" << filename_prefix << ".level" << _gridn <<"."<< time_step << "." << order << ".pvtu";
    Pfout.open(Pfilename.str().c_str());
    if (Pfout.is_open()) {
      std::cout << std::endl << " The output is printed to file " << Pfilename.str() << " in parallel VTK-XML (64-based) format" << std::endl;
    }
    else {
      std::cout << std::endl << " The output file "<< Pfilename.str() <<" cannot be opened.\n";
      abort();
    }
  }

  // *********** write pvtu header ***********
  Pfout<<"<?xml version=\"1.0\"?>" << std::endl;
  Pfout<<"<VTKFile type = \"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  Pfout<< "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
  for(int jproc=0;jproc<_nprocs;jproc++){
    Pfout<<"    <Piece Source=\""<<dirnamePVTK
         << filename_prefix << ".level" << _gridn << "." <<jproc<<"."<< time_step << "." << order << ".vtu"
	 <<"\"/>" << std::endl;
  }
  // ****************************************


  //count the own node dofs on all levels
  unsigned nvt = 0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    unsigned nvt_ig = _ml_mesh->GetLevel(ig)->own_size[index][_iproc];
    nvt += nvt_ig;
  }

  map < unsigned, unsigned > ghostMap;

  // count the ghost node dofs and the own element dofs element on all levels
  unsigned ghostMapCounter = 0;
  unsigned gridOffset = 0;

  //unsigned ghost_counter = 0;
  unsigned counter=0;
  unsigned nel=0;
  for (unsigned ig = _gridr-1u; ig<_gridn; ig++) {
    unsigned offset_iprc = _ml_mesh->GetLevel(ig)->MetisOffset[index][_iproc];
    for (int iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
      unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
      if ( ig == _gridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
	nel++;
	short unsigned ielt = _ml_mesh->GetLevel(ig)->el->GetElementType(kel);
	for (unsigned j=0; j<_ml_mesh->GetLevel(ig)->el->GetElementDofNumber(kel,index); j++) {
          counter++;
	  unsigned loc_vtk_conn = FemusToVTKorToXDMFConn[j];
	  unsigned jnode=_ml_mesh->GetLevel(ig)->el->GetMeshDof(kel, loc_vtk_conn, index);
	  unsigned jnodeMetis = _ml_mesh->GetLevel(ig)->GetMetisDof(jnode, index);
	  if( jnodeMetis < offset_iprc ){ // check if jnodeMetis is a ghost node
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

  unsigned nvtOwned = nvt;
  nvt += ghostMap.size(); // total node dofs (own + ghost)

  const unsigned dim_array_coord [] = { nvt*3*sizeof(float) };
  const unsigned dim_array_conn[]   = { counter*sizeof(int) };
  const unsigned dim_array_off []   = { nel*sizeof(int) };
  const unsigned dim_array_type []  = { nel*sizeof(short unsigned) };
  const unsigned dim_array_reg []   = { nel*sizeof(short unsigned) };
  const unsigned dim_array_elvar [] = { nel*sizeof(float) };
  const unsigned dim_array_ndvar [] = { nvt*sizeof(float) };

  // initialize common buffer_void memory
  unsigned buffer_size=(dim_array_coord[0]>dim_array_conn[0])? dim_array_coord[0] : dim_array_conn[0];
  void *buffer_void=new char [buffer_size];
  char *buffer_char=static_cast <char *>(buffer_void);

  size_t cch;
  cch = b64::b64_encode(&buffer_char[0], buffer_size , NULL, 0);
  vector <char> enc;
  enc.resize(cch);
  char *pt_char;

  fout  << "    <Piece NumberOfPoints= \"" << nvt << "\" NumberOfCells= \"" << nel << "\" >" << std::endl;

  //-----------------------------------------------------------------------------------------------
  // print coordinates *********************************************Solu*******************************************
  fout  << "      <Points>" << std::endl;
  fout  << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;

  Pfout << "    <PPoints>" << std::endl;
  Pfout << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\"/>" << std::endl;

  vector < NumericVector* > mysol(_gridn);
  for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    mysol[ig] = NumericVector::build().release();

    if(n_processors()==1) { // IF SERIAL
      mysol[ig]->init(_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs],
		      _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs],false,SERIAL);
    }
    else{ // IF PARALLEL
      mysol[ig]->init(_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs],
		      _ml_mesh->GetLevel(ig)->own_size[index][_iproc],
		      _ml_mesh->GetLevel(ig)->ghost_nd_mts[index][_iproc],false,GHOSTED );
    }
  }

  // point pointer to common mamory area buffer of void type;
  float *var_coord= static_cast<float*>(buffer_void);
  unsigned offset_ig = 0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    unsigned offset_iprc = _ml_mesh->GetLevel(ig)->MetisOffset[index][_iproc];
    unsigned nvt_ig = _ml_mesh->GetLevel(ig)->own_size[index][_iproc];
    for (int i = 0; i < 3; i++) {
      mysol[ig]->matrix_mult(*_ml_mesh->GetLevel(ig)->_coordinate->_Sol[i],
			     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,2) );
      if( _surface && i == 2 ){
        unsigned indSurfVar=_ml_sol->GetIndex(_surfaceVariable.c_str());
        mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Sol[indSurfVar],
                               *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,_ml_sol->GetSolutionType(indSurfVar)) );
        mysol[ig]->close();
      }
      for (unsigned ii = 0; ii < nvt_ig; ii++) {
	var_coord[ offset_ig + ii*3 + i] = (*mysol[ig])(ii + offset_iprc);
      }
      if (_ml_sol != NULL && _moving_mesh  && _ml_mesh->GetLevel(0)->GetDimension() > i)  { // if moving mesh
	unsigned indDXDYDZ=_ml_sol->GetIndex(_moving_vars[i].c_str());
	mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Sol[indDXDYDZ],
			       *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,_ml_sol->GetSolutionType(indDXDYDZ)) );
	for (unsigned ii=0; ii<nvt_ig; ii++)
 	  var_coord[ offset_ig + ii*3 + i] += (*mysol[ig])(ii + offset_iprc);
      }
    }
    offset_ig += 3 * nvt_ig;
  }
  //print ghost nodes
  for (int i=0; i<3; i++) {
    for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
      mysol[ig]->matrix_mult(*_ml_mesh->GetLevel(ig)-> _coordinate->_Sol[i],
			     *_ml_mesh->GetLevel(ig)-> GetQitoQjProjection(index,2) );
      if(i==2 && _surface){
        unsigned indSurfVar=_ml_sol->GetIndex(_surfaceVariable.c_str());
        mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Sol[indSurfVar],
                               *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,_ml_sol->GetSolutionType(indSurfVar)) );
        mysol[ig]->close();
      }
    }
    gridOffset = 0;
    unsigned ig = _gridr-1u;
    for (std::map <unsigned, unsigned>::iterator it=ghostMap.begin(); it!=ghostMap.end(); ++it){
      while( it->first >= gridOffset + _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs] ) {
	gridOffset += _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
	ig++;
      }
      var_coord[ offset_ig + 3*it->second + i ] = (*mysol[ig])( it->first - gridOffset);
    }
  }
  for (int i=0; i<3; i++) {  // if moving mesh
    if (_ml_sol != NULL && _moving_mesh  && _ml_mesh->GetLevel(0)->GetDimension() > i)  {
      unsigned indDXDYDZ=_ml_sol->GetIndex(_moving_vars[i].c_str());
      for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
	mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Sol[indDXDYDZ],
			       *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,_ml_sol->GetSolutionType(indDXDYDZ)) );
      }
      gridOffset = 0;
      unsigned ig = _gridr-1u;
      for (std::map <unsigned, unsigned>::iterator it=ghostMap.begin(); it!=ghostMap.end(); ++it){
	while( it->first >= gridOffset + _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs] ) {
	  gridOffset += _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
	  ig++;
	}
	var_coord[ offset_ig + 3*it->second + i ] += (*mysol[ig])( it->first - gridOffset);
      }
    }
  }

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

  fout  << "        </DataArray>" << std::endl;
  fout  << "      </Points>" << std::endl;
  Pfout << "    </PPoints>" << std::endl;
  //-----------------------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------------------
  // Printing of element connectivity - offset - format type  *
  fout  << "      <Cells>" << std::endl;
  Pfout << "    <PCells>" << std::endl;
  //-----------------------------------------------------------------------------------------------
  //print connectivity
  fout  << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">" << std::endl;
  Pfout << "      <PDataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\"/>" << std::endl;

  // point pointer to common mamory area buffer of void type;
  int *var_conn = static_cast <int*> (buffer_void);
  icount = 0;
  //ghost_counter = 0;
  gridOffset = 0;
  unsigned offset_nvt=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    unsigned offset_iprc = _ml_mesh->GetLevel(ig)->MetisOffset[index][_iproc];
    unsigned nvt_ig= _ml_mesh->GetLevel(ig)->own_size[index][_iproc];
    for (int iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
      unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
      if ( ig == _gridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
        for (unsigned j=0; j<_ml_mesh->GetLevel(ig)->el->GetElementDofNumber(kel,index); j++) {
	  unsigned loc_vtk_conn = FemusToVTKorToXDMFConn[j];
	  unsigned jnode=_ml_mesh->GetLevel(ig)->el->GetMeshDof(kel, loc_vtk_conn, index);
	  unsigned jnodeMetis = _ml_mesh->GetLevel(ig)->GetMetisDof(jnode, index);
	  var_conn[icount] = (jnodeMetis >= offset_iprc )? offset_nvt + jnodeMetis - offset_iprc :
							   nvtOwned + ghostMap[gridOffset+jnodeMetis];
	  icount++;
	}
      }
    }
    offset_nvt+= nvt_ig;
    gridOffset += _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
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
  fout << "        </DataArray>" << std::endl;
  //------------------------------------------------------------------------------------------------

  //-------------------------------------------------------------------------------------------------
  //printing offset
  fout  << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">" << std::endl;
  Pfout << "      <PDataArray type=\"Int32\" Name=\"offsets\" format=\"binary\"/>" << std::endl;



  // point pointer to common mamory area buffer of void type;
  int *var_off=static_cast <int*>(buffer_void);
  icount = 0;
  int offset_el=0;
  //print offset array

  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (int iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
      unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
      if ( ig == _gridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
  	offset_el += _ml_mesh->GetLevel(ig)->el->GetElementDofNumber(kel,index);
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

  fout  << std::endl;

  fout  << "        </DataArray>" << std::endl;
  //--------------------------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------------------------
  //Element format type : 23:Serendipity(8-nodes)  28:Quad9-Biquadratic
  fout  << "        <DataArray type=\"UInt16\" Name=\"types\" format=\"binary\">" << std::endl;
  Pfout << "      <PDataArray type=\"UInt16\" Name=\"types\" format=\"binary\"/>" << std::endl;

  // point pointer to common mamory area buffer of void type;
  unsigned short *var_type = static_cast <unsigned short*> (buffer_void);
  icount=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (int iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
      unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
      if ( ig == _gridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
  	short unsigned ielt= _ml_mesh->GetLevel(ig)->el->GetElementType(kel);
	var_type[icount] = femusToVtkCellType[index][ielt];
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

  fout  << std::endl;
  fout  << "        </DataArray>" << std::endl;
  //----------------------------------------------------------------------------------------------------
//
  fout  << "      </Cells>" << std::endl;
  Pfout << "    </PCells>" << std::endl;
  //--------------------------------------------------------------------------------------------------

  // /Print Cell Data ****************************************************************************
  fout  << "      <CellData Scalars=\"scalars\">" << std::endl;
  Pfout << "    <PCellData Scalars=\"scalars\">" << std::endl;
  //--------------------------------------------------------------------------------------------
  // Print Regions
  fout  << "        <DataArray type=\"UInt16\" Name=\"Regions\" format=\"binary\">" << std::endl;
  Pfout << "      <PDataArray type=\"UInt16\" Name=\"Regions\" format=\"binary\"/>" << std::endl;


  // point pointer to common mamory area buffer of void type;
  unsigned short* var_reg=static_cast <unsigned short*> (buffer_void);

  icount=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (int iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
      unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
      if ( ig == _gridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
  	var_reg[icount]= _ml_mesh->GetLevel(ig)->el->GetElementGroup(kel);
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

  fout  << std::endl;
  fout  << "        </DataArray>" << std::endl;

  //-----------------------------------------------------------------------------------------------------
  // Print Metis Partitioning
  fout  << "        <DataArray type=\"UInt16\" Name=\"Domain_partition\" format=\"binary\">" << std::endl;
  Pfout << "      <PDataArray type=\"UInt16\" Name=\"Domain_partition\" format=\"binary\"/>" << std::endl;

  // point pointer to common mamory area buffer of void type;
  unsigned short* var_proc=static_cast <unsigned short*> (buffer_void);

  icount=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (int iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
      unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
      if ( ig == _gridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
  	var_proc[icount]=(unsigned short)(_ml_mesh->GetLevel(ig)->epart[kel]);
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

  fout  << std::endl;
  fout  << "        </DataArray>" << std::endl;
  if (_ml_sol == NULL) {
    delete [] var_proc;
  }

  bool print_all = 0;
    for (unsigned ivar=0; ivar < vars.size(); ivar++){
      print_all += !(vars[ivar].compare("All")) + !(vars[ivar].compare("all")) + !(vars[ivar].compare("ALL"));
  }
  if (_ml_sol != NULL) {
    //Print Solution (on element) ***************************************************************
    for (unsigned i=0; i<(!print_all)*vars.size() + print_all*_ml_sol->GetSolutionSize(); i++) {
      unsigned solIndex=( print_all == 0 ) ? _ml_sol->GetIndex(vars[i].c_str()):i;
      if (3 <= _ml_sol->GetSolutionType(solIndex)) {
	fout  << "        <DataArray type=\"Float32\" Name=\"" << _ml_sol->GetSolutionName(solIndex) <<"\" format=\"binary\">" << std::endl;
	Pfout << "      <PDataArray type=\"Float32\" Name=\"" << _ml_sol->GetSolutionName(solIndex) <<"\" format=\"binary\"/>" << std::endl;
	// point pointer to common memory area buffer of void type;
	float *var_el = static_cast< float*> (buffer_void);
	icount=0;
	for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
	  for (unsigned iel=_ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc]; iel < _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
	    unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
	    if ( ig == _gridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
	      unsigned iel_Metis = _ml_mesh->GetLevel(ig)->GetMetisDof(kel,_ml_sol->GetSolutionType(i));
	      var_el[icount] = (*_ml_sol->GetSolutionLevel(ig)->_Sol[i])(iel_Metis);
	      icount++;
	    }
	  }
	}

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
	fout << "        </DataArray>" << std::endl;
	//----------------------------------------------------------------------------------------------------
      }
    } //end _ml_sol != NULL
  }
  fout  << "      </CellData>" << std::endl;
  Pfout << "    </PCellData>" << std::endl;
  //   //------------------------------------------------------------------------------------------------

  if (_ml_sol != NULL) {
  //   //------------------------------------------------------------------------------------------------
  // / Print Solution (on nodes) ********************************************************************
  fout  << "      <PointData Scalars=\"scalars\"> " << std::endl;
  Pfout << "    <PPointData Scalars=\"scalars\"> " << std::endl;
  //Loop on variables

  // point pointer to common memory area buffer of void type;
  float* var_nd = static_cast<float*>(buffer_void);
  for (unsigned i=0; i< (!print_all)*vars.size() + print_all*_ml_sol->GetSolutionSize(); i++) {
    unsigned solIndex=( print_all == 0 )?_ml_sol->GetIndex(vars[i].c_str()):i;
    if (_ml_sol->GetSolutionType(solIndex)<3) {
      fout  << "        <DataArray type=\"Float32\" Name=\"" << _ml_sol->GetSolutionName(solIndex) <<"\" format=\"binary\">" << std::endl;
      Pfout << "      <PDataArray type=\"Float32\" Name=\"" << _ml_sol->GetSolutionName(solIndex) <<"\" format=\"binary\"/>" << std::endl;

      //print solutions on nodes dimension
      cch = b64::b64_encode(&dim_array_ndvar[0], sizeof(dim_array_ndvar), NULL, 0);
      b64::b64_encode(&dim_array_ndvar[0], sizeof(dim_array_ndvar), &enc[0], cch);
      pt_char=&enc[0];
      for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;

      unsigned offset_ig = 0;
      for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
	unsigned offset_iprc = _ml_mesh->GetLevel(ig)->MetisOffset[index][_iproc];
	unsigned nvt_ig = _ml_mesh->GetLevel(ig)->own_size[index][_iproc];

	mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Sol[solIndex],
			       *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,_ml_sol->GetSolutionType(solIndex)) );

	for (unsigned ii = 0; ii < nvt_ig; ii++) {
	  var_nd[ offset_ig + ii ] = (*mysol[ig])(ii + offset_iprc);
	}
	offset_ig += nvt_ig;
      }


      gridOffset = 0;
      unsigned ig = _gridr-1u;
      for (std::map <unsigned, unsigned>::iterator it=ghostMap.begin(); it!=ghostMap.end(); ++it){
	while( it->first >= gridOffset + _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs] ) {
	  gridOffset += _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
	  ig++;
	}
	var_nd[ offset_ig + it->second ] = (*mysol[ig])( it->first - gridOffset);
      }




//       gridOffset = 0;
//       for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
// 	unsigned offset_iprc = _ml_mesh->GetLevel(ig)->MetisOffset[index][_iproc];
// 	for (int iel = _ml_mesh->GetLevel(ig)-> IS_Mts2Gmt_elem_offset[_iproc];
// 		 iel < _ml_mesh->GetLevel(ig)-> IS_Mts2Gmt_elem_offset[_iproc+1]; iel++) {
// 	  unsigned kel = _ml_mesh->GetLevel(ig)->IS_Mts2Gmt_elem[iel];
// 	  if ( ig == _gridn-1u || 0 == _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(kel)) {
// 	    short unsigned ielt=_ml_mesh->GetLevel(ig)->el->GetElementType(kel);
// 	    for (unsigned j=0; j<_ml_mesh->GetLevel(ig)->el->GetElementDofNumber(kel,index); j++) {
// 	      unsigned loc_vtk_conn = FemusToVTKorToXDMFConn[j];
// 	      unsigned jnode=_ml_mesh->GetLevel(ig)->el->GetMeshDof(kel, loc_vtk_conn, index);
// 	      unsigned jnodeMetis = _ml_mesh->GetLevel(ig)->GetMetisDof(jnode, index);
// 	      if( jnodeMetis < offset_iprc ){
// 		//var_nd[ offset_ig + icounter] = (*mysol[ig])(jnodeMetis);
// 		var_nd[ offset_ig +  ghostMap[gridOffset+jnodeMetis] ] = (*mysol[ig])(jnodeMetis);
// 		//icounter++;
// 	      }
// 	    }
// 	  }
// 	}
// 	gridOffset += _ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
//       }

      cch = b64::b64_encode(&var_nd[0], dim_array_ndvar [0], NULL, 0);
      b64::b64_encode(&var_nd[0], dim_array_ndvar [0], &enc[0], cch);
      pt_char=&enc[0];
      for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
      fout << std::endl;

      fout  << "        </DataArray>" << std::endl;
      } //endif
    } // end for sol
    fout  << "      </PointData>" << std::endl;
    Pfout << "    </PPointData>" << std::endl;
    delete [] var_nd;
  }  //end _ml_sol != NULL

  //------------------------------------------------------------------------------------------------

  fout << "    </Piece>" << std::endl;
  fout << "  </UnstructuredGrid>" << std::endl;
  fout << "</VTKFile>" << std::endl;
  fout.close();

  Pfout << "  </PUnstructuredGrid>" << std::endl;
  Pfout << "</VTKFile>" << std::endl;
  Pfout.close();


  //-----------------------------------------------------------------------------------------------------
  //free memory
  for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    delete mysol[ig];
  }

  //--------------------------------------------------------------------------------------------------------
  return;
}


// void VTKWriter::write(const std::string output_path, const char order[], const std::vector < std::string > & vars, const unsigned time_step) const {
//
//   ParallelWrite(output_path, order, vars, time_step);
//   bool print_all = 0;
//   for (unsigned ivar=0; ivar < vars.size(); ivar++){
//     print_all += !(vars[ivar].compare("All")) + !(vars[ivar].compare("all")) + !(vars[ivar].compare("ALL"));
//   }
//
//    int icount;
//   unsigned index=0;
//   if 	  (!strcmp(order,"linear")) 	 index=0; //linear
//   else if (!strcmp(order,"quadratic")) 	 index=1; //quadratic
//   else if (!strcmp(order,"biquadratic")) index=2; //biquadratic
//
//
//   std::string filename_prefix;
//   if( _ml_sol != NULL ) filename_prefix = "sol";
//   else filename_prefix = "mesh";
//
//   std::ostringstream filename;
//   filename << output_path << "/" << filename_prefix << ".level" << _gridn << "." << time_step << "." << order << ".vtu";
//
//   std::ofstream fout;
//
//   if(_iproc!=0) {
//     fout.rdbuf();   //redirect to dev_null
//   }
//   else {
//     fout.open(filename.str().c_str());
//     if (fout.is_open()) {
//       std::cout << std::endl << " The output is printed to file " << filename.str() << " in VTK-XML (64-based) format" << std::endl;
//     }
//     else {
//       std::cout << std::endl << " The output file "<< filename.str() <<" cannot be opened.\n";
//       abort();
//     }
//   }
//
//   // head ************************************************
//   fout<<"<?xml version=\"1.0\"?>" << std::endl;
//   fout<<"<VTKFile type = \"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
//   fout << " <UnstructuredGrid>" << std::endl;
//   //-----------------------------------------------------------------------------------------------------------
//
//   //----------------------------------------------------------------------------------------------------------
//
//   unsigned nvt=0;
//   for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     unsigned nvt_ig=_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
//     nvt+=nvt_ig;
//   }
//
//   unsigned nel=0;
//   unsigned counter=0;
//   for (unsigned ig=_gridr-1u; ig<_gridn-1u; ig++) {
//     nel+=( _ml_mesh->GetLevel(ig)->GetNumberOfElements() - _ml_mesh->GetLevel(ig)->el->GetRefinedElementNumber());
//     counter+=(_ml_mesh->GetLevel(ig)->el->GetElementNumber("Hex")-_ml_mesh->GetLevel(ig)->el->GetRefinedElementNumber("Hex"))*NVE[0][index];
//     counter+=(_ml_mesh->GetLevel(ig)->el->GetElementNumber("Tet")-_ml_mesh->GetLevel(ig)->el->GetRefinedElementNumber("Tet"))*NVE[1][index];
//     counter+=(_ml_mesh->GetLevel(ig)->el->GetElementNumber("Wedge")-_ml_mesh->GetLevel(ig)->el->GetRefinedElementNumber("Wedge"))*NVE[2][index];
//     counter+=(_ml_mesh->GetLevel(ig)->el->GetElementNumber("Quad")-_ml_mesh->GetLevel(ig)->el->GetRefinedElementNumber("Quad"))*NVE[3][index];
//     counter+=(_ml_mesh->GetLevel(ig)->el->GetElementNumber("Triangle")-_ml_mesh->GetLevel(ig)->el->GetRefinedElementNumber("Triangle"))*NVE[4][index];
//     counter+=(_ml_mesh->GetLevel(ig)->el->GetElementNumber("Line")-_ml_mesh->GetLevel(ig)->el->GetRefinedElementNumber("Line"))*NVE[5][index];
//   }
//   nel+=_ml_mesh->GetLevel(_gridn-1u)->GetNumberOfElements();
//   counter+=_ml_mesh->GetLevel(_gridn-1u)->el->GetElementNumber("Hex")*NVE[0][index];
//   counter+=_ml_mesh->GetLevel(_gridn-1u)->el->GetElementNumber("Tet")*NVE[1][index];
//   counter+=_ml_mesh->GetLevel(_gridn-1u)->el->GetElementNumber("Wedge")*NVE[2][index];
//   counter+=_ml_mesh->GetLevel(_gridn-1u)->el->GetElementNumber("Quad")*NVE[3][index];
//   counter+=_ml_mesh->GetLevel(_gridn-1u)->el->GetElementNumber("Triangle")*NVE[4][index];
//   counter+=_ml_mesh->GetLevel(_gridn-1u)->el->GetElementNumber("Line")*NVE[5][index];
//
//   const unsigned dim_array_coord [] = { nvt*3*sizeof(float) };
//   const unsigned dim_array_conn[]   = { counter*sizeof(int) };
//   const unsigned dim_array_off []   = { nel*sizeof(int) };
//   const unsigned dim_array_type []  = { nel*sizeof(short unsigned) };
//   const unsigned dim_array_reg []   = { nel*sizeof(short unsigned) };
//   const unsigned dim_array_elvar [] = { nel*sizeof(float) };
//   const unsigned dim_array_ndvar [] = { nvt*sizeof(float) };
//
//   // initialize common buffer_void memory
//   unsigned buffer_size=(dim_array_coord[0]>dim_array_conn[0])?dim_array_coord[0]:dim_array_conn[0];
//   void *buffer_void=new char [buffer_size];
//   char *buffer_char=static_cast <char *>(buffer_void);
//
//   size_t cch;
//   cch = b64::b64_encode(&buffer_char[0], buffer_size , NULL, 0);
//   vector <char> enc;
//   enc.resize(cch);
//   char *pt_char;
//
//   fout << "  <Piece NumberOfPoints= \"" << nvt << "\" NumberOfCells= \"" << nel << "\" >" << std::endl;
//
//   //-----------------------------------------------------------------------------------------------
//   // print coordinates *********************************************Solu*******************************************
//   fout << "   <Points>" << std::endl;
//   fout << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;
//
//
//   vector <NumericVector*> mysol(_gridn);
//   for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     mysol[ig] = NumericVector::build().release();
//     mysol[ig]->init(_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs],_ml_mesh->GetLevel(ig)->own_size[index][_iproc],true,AUTOMATIC);
//   }
//
//   // point pointer to common mamory area buffer of void type;
//   float *var_coord= static_cast<float*>(buffer_void);
//
//   unsigned offset_nvt3=0;
//   for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     std::vector<double> v_local;
//     unsigned nvt_ig=_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
//     for(int kk=0;kk<3;kk++) {
//       mysol[ig]->matrix_mult(*_ml_mesh->GetLevel(ig)->_coordinate->_Sol[kk],
// 			     *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,2));
//       mysol[ig]->localize_to_one(v_local,0);
//       if(_iproc==0) {
// 	for (unsigned i=0; i<nvt_ig; i++) {
// 	  var_coord[offset_nvt3+i*3+kk] = v_local[i];
// 	}
//       } //if iproc
//     } //loop over dimension
//     offset_nvt3+=3*nvt_ig;
//   }
//
//
//   if (_ml_sol != NULL && _moving_mesh == true) {
//
//     unsigned offset_nvt3=0;
//     unsigned indDXDYDZ[3];
//     indDXDYDZ[0]=_ml_sol->GetIndex(_moving_vars[0].c_str());
//     indDXDYDZ[1]=_ml_sol->GetIndex(_moving_vars[1].c_str());
//     if(_ml_mesh->GetLevel(0)->GetDimension() == 3) {
//       indDXDYDZ[2]=_ml_sol->GetIndex(_moving_vars[2].c_str());
//     }
//
//     for(unsigned ig=_gridr-1u; ig<_gridn; ig++){
//       std::vector<double> v_local;
//       unsigned nvt_ig=_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
//       for(int kk=0;kk<_ml_mesh->GetLevel(0)->GetDimension();kk++) {
// 	mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Sol[indDXDYDZ[kk]],
// 			       *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,_ml_sol->GetSolutionType(indDXDYDZ[kk])));
//         mysol[ig]->localize_to_one(v_local,0);
// 	if(_iproc==0) {
// 	  for (unsigned i=0; i<nvt_ig; i++) {
// 	    var_coord[offset_nvt3+i*3+kk] += v_local[i];
// 	  }
// 	} //if iproc
//       } //loop over dimension
//       offset_nvt3+=3*nvt_ig;
//     }
//
//   } //end moving mesh and ml_sol
//
//   if(_iproc==0) {
//     //print coordinates dimension
//     cch = b64::b64_encode(&dim_array_coord[0], sizeof(dim_array_coord), NULL, 0);
//     b64::b64_encode(&dim_array_coord[0], sizeof(dim_array_coord), &enc[0], cch);
//     pt_char=&enc[0];
//     for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//     //print coordinates array
//     cch = b64::b64_encode(&var_coord[0], dim_array_coord[0] , NULL, 0);
//     b64::b64_encode(&var_coord[0], dim_array_coord[0], &enc[0], cch);
//     pt_char=&enc[0];
//     for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//     fout << std::endl;
//   }
//   fout << "    </DataArray>" << std::endl;
//   fout << "   </Points>" << std::endl;
//   //-----------------------------------------------------------------------------------------------
//
//   //-----------------------------------------------------------------------------------------------
//   // Printing of element connectivity - offset - format type  *
//   fout << "   <Cells>" << std::endl;
//
//   //-----------------------------------------------------------------------------------------------
//   //print connectivity
//   fout << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">" << std::endl;
//
//   // point pointer to common mamory area buffer of void type;
//   int *var_conn = static_cast <int*> (buffer_void);
//   icount = 0;
//   unsigned offset_nvt=0;
//   for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     for (unsigned iel=0; iel<_ml_mesh->GetLevel(ig)->GetNumberOfElements(); iel++) {
//       if (ig==_gridn-1u || _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(iel)==0) {
//         for (unsigned j=0; j<_ml_mesh->GetLevel(ig)->el->GetElementDofNumber(iel,index); j++) {
// 	  unsigned loc_vtk_conn = FemusToVTKorToXDMFConn[j];
// 	  unsigned jnode=_ml_mesh->GetLevel(ig)->el->GetElementVertexIndex(iel,loc_vtk_conn)-1u;
// 	  unsigned jnodeMetis = _ml_mesh->GetLevel(ig)->GetMetisDof(jnode,index);
// 	  var_conn[icount] = offset_nvt + jnodeMetis;
// 	  icount++;
// 	}
//       }
//     }
//     offset_nvt+=_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
//   }
//
//   //print connectivity dimension
//   cch = b64::b64_encode(&dim_array_conn[0], sizeof(dim_array_conn), NULL, 0);
//   b64::b64_encode(&dim_array_conn[0], sizeof(dim_array_conn), &enc[0], cch);
//   pt_char=&enc[0];
//   for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//   //print connectivity array
//   cch = b64::b64_encode(&var_conn[0], dim_array_conn[0] , NULL, 0);
//   b64::b64_encode(&var_conn[0], dim_array_conn[0], &enc[0], cch);
//   pt_char=&enc[0];
//   for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//   fout << std::endl;
//
//   fout << "    </DataArray>" << std::endl;
//   //------------------------------------------------------------------------------------------------
//
//   //-------------------------------------------------------------------------------------------------
//   //printing offset
//   fout << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">" << std::endl;
//
//   // point pointer to common mamory area buffer of void type;
//   int *var_off=static_cast <int*>(buffer_void);
//   icount = 0;
//   int offset_el=0;
//   //print offset array
//   for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     for (unsigned iel=0; iel<_ml_mesh->GetLevel(ig)->GetNumberOfElements(); iel++) {
//       if ( ig==_gridn-1u || 0==_ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(iel)) {
//         offset_el += _ml_mesh->GetLevel(ig)->el->GetElementDofNumber(iel,index);
//         var_off[icount] = offset_el;
// 	icount++;
//       }
//     }
//   }
//
//   //print offset dimension
//   cch = b64::b64_encode(&dim_array_off[0], sizeof(dim_array_off), NULL, 0);
//   b64::b64_encode(&dim_array_off[0], sizeof(dim_array_off), &enc[0], cch);
//   pt_char=&enc[0];
//   for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//   //print offset array
//   cch = b64::b64_encode(&var_off[0], dim_array_off[0] , NULL, 0);
//   b64::b64_encode(&var_off[0], dim_array_off[0], &enc[0], cch);
//   pt_char=&enc[0];
//   for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//   fout << std::endl;
//
//   fout << "    </DataArray>" << std::endl;
//   //--------------------------------------------------------------------------------------------------
//
//   //--------------------------------------------------------------------------------------------------
//   //Element format type : 23:Serendipity(8-nodes)  28:Quad9-Biquadratic
//   fout << "    <DataArray type=\"UInt16\" Name=\"types\" format=\"binary\">" << std::endl;
//
//   // point pointer to common mamory area buffer of void type;
//   unsigned short *var_type = static_cast <unsigned short*> (buffer_void);
//   icount=0;
//   for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     for (unsigned ii=0; ii<_ml_mesh->GetLevel(ig)->GetNumberOfElements(); ii++) {
//       if (ig==_gridn-1u || _ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(ii)==0) {
// 	short unsigned ielt= _ml_mesh->GetLevel(ig)->el->GetElementType(ii);
// 	var_type[icount] = femusToVtkCellType[index][ielt];
// 	icount++;
//       }
//     }
//   }
//
//   //print element format dimension
//   cch = b64::b64_encode(&dim_array_type[0], sizeof(dim_array_type), NULL, 0);
//   b64::b64_encode(&dim_array_type[0], sizeof(dim_array_type), &enc[0], cch);
//   pt_char=&enc[0];
//   for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//   //print element format array
//   cch = b64::b64_encode(&var_type[0], dim_array_type[0] , NULL, 0);
//   b64::b64_encode(&var_type[0], dim_array_type[0], &enc[0], cch);
//   pt_char=&enc[0];
//   for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//   fout << std::endl;
//   fout << "    </DataArray>" << std::endl;
//   //----------------------------------------------------------------------------------------------------
//
//   fout << "   </Cells>" << std::endl;
//   //--------------------------------------------------------------------------------------------------
//
//   // /Print Cell Data ****************************************************************************
//   fout << "   <CellData Scalars=\"scalars\">" << std::endl;
//
//   //--------------------------------------------------------------------------------------------
//   // Print Regions
//   fout << "    <DataArray type=\"UInt16\" Name=\"Regions\" format=\"binary\">" << std::endl;
//
//   // point pointer to common mamory area buffer of void type;
//   unsigned short* var_reg=static_cast <unsigned short*> (buffer_void);
//   icount=0;
//   for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     for (unsigned ii=0; ii<_ml_mesh->GetLevel(ig)->GetNumberOfElements(); ii++) {
//       if ( ig==_gridn-1u || 0==_ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
// 	var_reg[icount]= _ml_mesh->GetLevel(ig)->el->GetElementGroup(ii);
// 	icount++;
//       }
//     }
//   }
//
//   //print regions dimension
//   cch = b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), NULL, 0);
//   b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), &enc[0], cch);
//   pt_char=&enc[0];
//   for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//
//   //print regions array
//   cch = b64::b64_encode(&var_reg[0], dim_array_reg[0] , NULL, 0);
//   b64::b64_encode(&var_reg[0], dim_array_reg[0], &enc[0], cch);
//   pt_char=&enc[0];
//   for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//   fout << std::endl;
//   fout << "    </DataArray>" << std::endl;
//
//   //-----------------------------------------------------------------------------------------------------
//   // Print Metis Partitioning
//   fout << "    <DataArray type=\"UInt16\" Name=\"Domain_partition\" format=\"binary\">" << std::endl;
//
//   // point pointer to common mamory area buffer of void type;
//   unsigned short* var_proc=static_cast <unsigned short*> (buffer_void);
//   icount=0;
//   for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     for (unsigned ii=0; ii<_ml_mesh->GetLevel(ig)->GetNumberOfElements(); ii++) {
//       if ( ig==_gridn-1u || 0==_ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
// 	var_proc[icount]=(unsigned short)(_ml_mesh->GetLevel(ig)->epart[ii]);
// 	icount++;
//       }
//     }
//   }
//
//   //print regions dimension
//   cch = b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), NULL, 0);
//   b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), &enc[0], cch);
//   pt_char=&enc[0];
//   for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//
//   //print regions array
//   cch = b64::b64_encode(&var_proc[0], dim_array_reg[0] , NULL, 0);
//   b64::b64_encode(&var_proc[0], dim_array_reg[0], &enc[0], cch);
//   pt_char=&enc[0];
//   for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//   fout << std::endl;
//   fout << "    </DataArray>" << std::endl;
//
//   if (_ml_sol != NULL) {
//
//   //Print Solution (on element) ***************************************************************
//   for (unsigned i=0; i<(!print_all)*vars.size() + print_all*_ml_sol->GetSolutionSize(); i++) {
//     unsigned solIndex=( print_all == 0 ) ? _ml_sol->GetIndex(vars[i].c_str()):i;
//     if (3 <= _ml_sol->GetSolutionType(solIndex)) {
//       fout << "    <DataArray type=\"Float32\" Name=\"" << _ml_sol->GetSolutionName(solIndex) <<"\" format=\"binary\">" << std::endl;
//       // point pointer to common memory area buffer of void type;
//       float *var_el = static_cast< float*> (buffer_void);
//       icount=0;
//       for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
// 	vector<double> sol_local;
// 	_ml_sol->GetSolutionLevel(ig)->_Sol[solIndex]->localize_to_one(sol_local,0);
// 	for (unsigned ii=0; ii<_ml_mesh->GetLevel(ig)->GetNumberOfElements(); ii++) {
// 	  if (ig==_gridn-1u || 0==_ml_mesh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
// 	    unsigned iel_Metis = _ml_mesh->GetLevel(ig)->GetMetisDof(ii,_ml_sol->GetSolutionType(solIndex));
// 	    var_el[icount]=sol_local[iel_Metis];
// 	    icount++;
// 	  }
// 	}
//       }
//
//       if(_iproc==0) {
//         //print solution on element dimension
//         cch = b64::b64_encode(&dim_array_elvar[0], sizeof(dim_array_elvar), NULL, 0);
//         b64::b64_encode(&dim_array_elvar[0], sizeof(dim_array_elvar), &enc[0], cch);
//         pt_char=&enc[0];
//         for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//         //print solution on element array
//         cch = b64::b64_encode(&var_el[0], dim_array_elvar[0] , NULL, 0);
//         b64::b64_encode(&var_el[0], dim_array_elvar[0], &enc[0], cch);
//         pt_char=&enc[0];
//         for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//         fout << std::endl;
//         fout << "    </DataArray>" << std::endl;
//       }
//       //----------------------------------------------------------------------------------------------------
//     }
//
//   } //end mlsol
//
//   }
//   fout << "   </CellData>" << std::endl;
//   //   //------------------------------------------------------------------------------------------------
//
//   if (_ml_sol != NULL) {
//   //   //------------------------------------------------------------------------------------------------
//   // / Print Solution (on nodes) ********************************************************************
//   fout<< " <PointData Scalars=\"scalars\"> " << std::endl;
//   //Loop on variables
//
//   // point pointer to common memory area buffer of void type;
//   float* var_nd = static_cast<float*>(buffer_void);
//   for (unsigned i=0; i<(!print_all)*vars.size()+ print_all*_ml_sol->GetSolutionSize(); i++) {
//     unsigned solIndex=( print_all == 0 )?_ml_sol->GetIndex(vars[i].c_str()):i;
//     if (_ml_sol->GetSolutionType(solIndex)<3) {
//       fout << " <DataArray type=\"Float32\" Name=\"" << _ml_sol->GetSolutionName(solIndex) <<"\" format=\"binary\">" << std::endl;
//       //print solutions on nodes dimension
//       cch = b64::b64_encode(&dim_array_ndvar[0], sizeof(dim_array_ndvar), NULL, 0);
//       b64::b64_encode(&dim_array_ndvar[0], sizeof(dim_array_ndvar), &enc[0], cch);
//       pt_char=&enc[0];
//       for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//
//       unsigned offset_nvt=0;
//       for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
// 	mysol[ig]->matrix_mult(*_ml_sol->GetSolutionLevel(ig)->_Sol[solIndex],
// 			       *_ml_mesh->GetLevel(ig)->GetQitoQjProjection(index,_ml_sol->GetSolutionType(solIndex)) );
// 	vector<double> sol_local;
// 	mysol[ig]->localize_to_one(sol_local,0);
// 	unsigned nvt_ig=_ml_mesh->GetLevel(ig)->MetisOffset[index][_nprocs];
// 	for (unsigned ii=0; ii<nvt_ig; ii++) {
// 	  var_nd[ii+offset_nvt] = sol_local[ii];
// 	}
// 	offset_nvt+=nvt_ig;
//       }
//
//       if(_iproc==0) {
//         cch = b64::b64_encode(&var_nd[0], dim_array_ndvar [0], NULL, 0);
//         b64::b64_encode(&var_nd[0], dim_array_ndvar [0], &enc[0], cch);
//         pt_char=&enc[0];
//         for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char;
//         fout << std::endl;
//
//         fout << "    </DataArray>" << std::endl;
//       }
//     } //endif
//   } // end for sol
//   fout << "   </PointData>" << std::endl;
//
//     delete [] var_nd;
//
//   }  //end _ml_sol != NULL
//
//   //------------------------------------------------------------------------------------------------
//
//   fout << "  </Piece>" << std::endl;
//   fout << " </UnstructuredGrid>" << std::endl;
//   fout << "</VTKFile>" << std::endl;
//   fout.close();
//
//   //-----------------------------------------------------------------------------------------------------
//   //free memory
//   for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     delete mysol[ig];
//   }
//
//
//
//   //--------------------------------------------------------------------------------------------------------
//   return;
//}






} //end namespace femus


