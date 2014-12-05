/*=========================================================================

 Program: FEMUS
 Module: XDMFWriter
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
#include "FEMTTUConfig.h"
#include "XDMFWriter.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "stdio.h"
#include "fstream"
#include "iostream"
#include <algorithm>  
#include <cstring>
#include <sstream>

#ifdef HAVE_HDF5
  #include "hdf5.h"


namespace femus {


#endif


XDMFWriter::XDMFWriter(MultiLevelSolution& ml_probl): Writer(ml_probl)
{
  
}

XDMFWriter::~XDMFWriter()
{
  
}

void XDMFWriter::write_system_solutions(const char order[], std::vector<std::string>& vars, const unsigned time_step) 
{ 
#ifdef HAVE_HDF5
  
  bool test_all=!(vars[0].compare("All"));
    
  unsigned index=0;
  unsigned index_nd=0;
  if(!strcmp(order,"linear")) {    //linear
    index=0;
    index_nd=0;
  }
  else if(!strcmp(order,"quadratic")) {  //quadratic
    index=1;
    index_nd=1;
  }
  else if(!strcmp(order,"biquadratic")) { //biquadratic
    index=3;
    index_nd=2;
  }

  const std::string type_el[4][6] = {{"Hexahedron","Tetrahedron","Wedge","Quadrilateral","Triangle","Edge"},
                                {"Hexahedron_20","Tetrahedron_10","Not_implemented","Quadrilateral_8","Triangle_6","Edge_3"},
			        {"Not_implemented","Not_implemented","Not_implemented","Not_implemented","Not_implemented","Not_implemented"},
                                {"Hexahedron_27","Not_implemented","Not_implemented","Quadrilateral_9","Triangle_6","Edge_3"}};
			 
  
  //I assume that the mesh is not mixed
  std::string type_elem;
  unsigned elemtype = _ml_sol._ml_msh->GetLevel(_gridn-1u)->el->GetElementType(0);
  type_elem = type_el[index][elemtype];
  
  if (type_elem.compare("Not_implemented") == 0) 
  {
    std::cerr << "XDMF-Writer error: element type not supported!" << std::endl;
    exit(1);
  }
  
  unsigned nvt=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    unsigned nvt_ig=_ml_sol._ml_msh->GetLevel(ig)->GetDofNumber(index_nd);
    nvt+=nvt_ig;
  } 
  
  // Printing connectivity
  unsigned nel=0;
  for(unsigned ig=0;ig<_gridn-1u;ig++) {
    nel+=( _ml_sol._ml_msh->GetLevel(ig)->GetNumberOfElements() - _ml_sol._ml_msh->GetLevel(ig)->el->GetRefinedElementNumber());
  }
  nel+=_ml_sol._ml_msh->GetLevel(_gridn-1u)->GetNumberOfElements();
  
  unsigned icount;
  unsigned el_dof_number  = _ml_sol._ml_msh->GetLevel(_gridn-1u)->el->GetElementDofNumber(0,index);
  int *var_int            = new int [nel*el_dof_number];
  float *var_el_f         = new float [nel];
  float *var_nd_f         = new float [nvt];

  char *filename= new char[60];
  std::ofstream fout;
  
  //--------------------------------------------------------------------------------------------------
  // Print The Xdmf wrapper
  sprintf(filename,"./output/mesh.level%d.%d.%s.xmf",_gridn,time_step,order);
  
  if(_iproc!=0) {
    fout.rdbuf();   //redirect to dev_null
  }
  else {
    fout.open(filename);
    if (!fout) {
      std::cout << std::endl << " The output mesh file "<<filename<<" cannot be opened.\n";
      exit(0);
    }
    else {
      std::cout << std::endl << " The output is printed to file " << filename << " in XDMF-HDF5 format" << std::endl;   
    }
  }

  // Print The HDF5 file
  sprintf(filename,"./mesh.level%d.%d.%s.h5",_gridn,time_step,order);
  // haed ************************************************
  fout<<"<?xml version=\"1.0\" ?>" << std::endl;
  fout<<"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd []\">"<< std::endl;
  fout<<"<Xdmf>"<<std::endl;
  fout<<"<Domain>"<<std::endl;
  fout<<"<Grid Name=\"Mesh\">"<<std::endl;
  fout<<"<Time Value =\""<< time_step<< "\" />"<<std::endl;
  fout<<"<Topology TopologyType=\""<< type_elem <<"\" NumberOfElements=\""<< nel <<"\">"<<std::endl;
  //Connectivity
  fout<<"<DataStructure DataType=\"Int\" Dimensions=\""<< nel*el_dof_number <<"\"" << "  Format=\"HDF\">" << std::endl;
  fout << filename << ":CONNECTIVITY" << std::endl;
  fout <<"</DataStructure>" << std::endl;
  fout << "</Topology>" << std::endl;
  fout << "<Geometry Type=\"X_Y_Z\">" << std::endl;
  //Node_X
  fout<<"<DataStructure DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << std::endl;
  fout << filename << ":NODES_X1" << std::endl;
  fout <<"</DataStructure>" << std::endl;
  //Node_Y
  fout<<"<DataStructure DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << std::endl;
  fout << filename << ":NODES_X2" << std::endl;
  fout <<"</DataStructure>" << std::endl;
  //Node_Z
  fout<<"<DataStructure DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << std::endl;
  fout << filename << ":NODES_X3" << std::endl;
  fout <<"</DataStructure>" << std::endl;
  fout <<"</Geometry>" << std::endl;
  //Regions
  fout << "<Attribute Name=\""<< "Regions"<<"\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
  fout << "<DataItem DataType=\"Int\" Dimensions=\""<< nel << "\"" << "  Format=\"HDF\">" << std::endl;
  fout << filename << ":REGIONS" << std::endl;
  fout << "</DataItem>" << std::endl;
  fout << "</Attribute>" << std::endl;
  // Solution Variables
  for (unsigned i=0; i<vars.size(); i++) {
    unsigned indx=_ml_sol.GetIndex(vars[i].c_str());  
    //Printing biquadratic solution on the nodes
    if(_ml_sol.GetSolutionType(indx)<3) {  
      fout << "<Attribute Name=\""<< _ml_sol.GetSolutionName(indx)<<"\" AttributeType=\"Scalar\" Center=\"Node\">" << std::endl;
      fout << "<DataItem DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << std::endl;
      fout << filename << ":" << _ml_sol.GetSolutionName(indx) << std::endl;
      fout << "</DataItem>" << std::endl;
      fout << "</Attribute>" << std::endl;
    }
    else if (_ml_sol.GetSolutionType(indx)>=3) {  //Printing picewise constant solution on the element
      fout << "<Attribute Name=\""<< _ml_sol.GetSolutionName(indx)<<"\" AttributeType=\"Scalar\" Center=\"Cell\">" << std::endl;
      fout << "<DataItem DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nel << "\"  Format=\"HDF\">" << std::endl;
      fout << filename << ":" << _ml_sol.GetSolutionName(indx) << std::endl;
      fout << "</DataItem>" << std::endl;
      fout << "</Attribute>" << std::endl;
    }
  }

  fout <<"</Grid>" << std::endl;
  fout <<"</Domain>" << std::endl;
  fout <<"</Xdmf>" << std::endl;
  fout.close();
  //----------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------------------------------------
  hid_t file_id;
  sprintf(filename,"./output/mesh.level%d.%d.%s.h5",_gridn,time_step,order);
  file_id = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  hsize_t dimsf[2];
  herr_t status;
  hid_t dataspace;
  hid_t dataset;
  
  //-----------------------------------------------------------------------------------------------------------
  // Printing nodes coordinates 
  
  PetscScalar *MYSOL[1]; //TODO

  for (int i=0; i<3; i++) {
    unsigned offset_nvt=0;
    for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
      NumericVector* mysol;
      mysol = NumericVector::build().release();
      //mysol->init(_ml_sol._ml_msh->GetLevel(ig)->GetDofNumber(index_nd),_ml_sol._ml_msh->GetLevel(ig)->GetDofNumber(index_nd),true,AUTOMATIC);
      mysol->init(_ml_sol._ml_msh->GetLevel(ig)->MetisOffset[index_nd][_nprocs],_ml_sol._ml_msh->GetLevel(ig)->own_size[index_nd][_iproc],true,AUTOMATIC);
      mysol->matrix_mult(*_ml_sol._ml_msh->GetLevel(ig)->_coordinate->_Sol[i],*Writer::_ProlQitoQj[index_nd][2][ig]);
      unsigned nvt_ig=_ml_sol._ml_msh->GetLevel(ig)->GetDofNumber(index_nd);
      for (unsigned ii=0; ii<nvt_ig; ii++) var_nd_f[ii+offset_nvt] = (*mysol)(ii);
      if (_moving_mesh) {
	unsigned varind_DXDYDZ=_ml_sol.GetIndex(_moving_vars[i].c_str());
	mysol->matrix_mult(*_ml_sol.GetSolutionLevel(ig)->_Sol[varind_DXDYDZ],*Writer::_ProlQitoQj[index_nd][_ml_sol.GetSolutionType(varind_DXDYDZ)][ig]);
	for (unsigned ii=0; ii<nvt_ig; ii++) var_nd_f[ii+offset_nvt] += (*mysol)(ii);
      }
      offset_nvt+=nvt_ig;
      delete mysol;
    }
    
    dimsf[0] = nvt ;  dimsf[1] = 1;
    std::ostringstream Name; Name << "/NODES_X" << i+1;
    dataspace = H5Screate_simple(2,dimsf, NULL);
    dataset   = H5Dcreate(file_id,Name.str().c_str(),H5T_NATIVE_FLOAT,
			  dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_nd_f[0]);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    
  }

  //-------------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------------------------------
  //connectivity
  icount = 0;
  unsigned offset_conn=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned iel=0; iel<_ml_sol._ml_msh->GetLevel(ig)->GetNumberOfElements(); iel++) {
      if (_ml_sol._ml_msh->GetLevel(ig)->el->GetRefinedElementIndex(iel)==0 || ig==_gridn-1u) {
        for (unsigned j=0; j<_ml_sol._ml_msh->GetLevel(ig)->el->GetElementDofNumber(iel,index); j++) {
	  unsigned vtk_loc_conn = map_pr[j];
	  unsigned jnode=_ml_sol._ml_msh->GetLevel(ig)->el->GetElementVertexIndex(iel,vtk_loc_conn)-1u;
	  unsigned jnode_Metis = _ml_sol._ml_msh->GetLevel(ig)->GetMetisDof(jnode,index_nd);
	  var_int[icount] = offset_conn + jnode_Metis;
	  icount++;
	}
      }
    }
    offset_conn += _ml_sol._ml_msh->GetLevel(ig)->GetDofNumber(index_nd);
  }
  
  dimsf[0] = nel*el_dof_number ;  dimsf[1] = 1;
  dataspace = H5Screate_simple(2,dimsf, NULL);
  dataset   = H5Dcreate(file_id,"/CONNECTIVITY",H5T_NATIVE_INT,
			dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status   = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_int[0]);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  //------------------------------------------------------------------------------------------------------
  
  
  //-------------------------------------------------------------------------------------------------------
  // print regions
  icount=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned ii=0; ii<_ml_sol._ml_msh->GetLevel(ig)->GetNumberOfElements(); ii++) {
      if (ig==_gridn-1u || 0==_ml_sol._ml_msh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
	unsigned iel_Metis = _ml_sol._ml_msh->GetLevel(ig)->GetMetisDof(ii,3);
	var_int[icount] = _ml_sol._ml_msh->GetLevel(ig)->el->GetElementGroup(iel_Metis);
	icount++;
      }
    }
  } 
   
  dimsf[0] = nel;  dimsf[1] = 1;
  dataspace = H5Screate_simple(2,dimsf, NULL);
  dataset   = H5Dcreate(file_id,"/REGIONS",H5T_NATIVE_INT,
			dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status   = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_int[0]);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  
  //-------------------------------------------------------------------------------------------------------
  // printing element variables
  for (unsigned i=0; i<(1-test_all)*vars.size()+test_all*_ml_sol.GetSolutionSize(); i++) {
    unsigned indx=(test_all==0)?_ml_sol.GetIndex(vars[i].c_str()):i;
    if (_ml_sol.GetSolutionType(indx)>=3) {
      icount=0;
      for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
	for (unsigned ii=0; ii<_ml_sol._ml_msh->GetLevel(ig)->GetNumberOfElements(); ii++) {
	  if (ig==_gridn-1u || 0==_ml_sol._ml_msh->GetLevel(ig)->el->GetRefinedElementIndex(ii)) {
	    unsigned iel_Metis = _ml_sol._ml_msh->GetLevel(ig)->GetMetisDof(ii,_ml_sol.GetSolutionType(indx));
	    var_el_f[icount]=(*_ml_sol.GetSolutionLevel(ig)->_Sol[indx])(iel_Metis);
	    icount++;
	  }
	}
      } 
     
      dimsf[0] = nel;  dimsf[1] = 1;
      dataspace = H5Screate_simple(2,dimsf, NULL);
      dataset   = H5Dcreate(file_id,_ml_sol.GetSolutionName(indx),H5T_NATIVE_FLOAT,
			    dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status   = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_el_f[0]);
      H5Sclose(dataspace);
      H5Dclose(dataset);
     
    }
  }
  
  //-------------------------------------------------------------------------------------------------------
  // printing nodes variables
  for (unsigned i=0; i<(1-test_all)*vars.size()+test_all*_ml_sol.GetSolutionSize(); i++) {
    unsigned indx=(test_all==0)?_ml_sol.GetIndex(vars[i].c_str()):i;
    if (_ml_sol.GetSolutionType(indx) < 3) {
      unsigned offset_nvt=0;
      for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
        NumericVector* mysol;
	mysol = NumericVector::build().release();
        //mysol->init(_ml_sol._ml_msh->GetLevel(ig)->GetDofNumber(index_nd),_ml_sol._ml_msh->GetLevel(ig)->GetDofNumber(index_nd),true,AUTOMATIC);
	mysol->init(_ml_sol._ml_msh->GetLevel(ig)->MetisOffset[index_nd][_nprocs],_ml_sol._ml_msh->GetLevel(ig)->own_size[index_nd][_iproc],true,AUTOMATIC);
	mysol->matrix_mult(*_ml_sol.GetSolutionLevel(ig)->_Sol[indx],*_ProlQitoQj[index_nd][_ml_sol.GetSolutionType(indx)][ig]);
	unsigned nvt_ig=_ml_sol._ml_msh->GetLevel(ig)->GetDofNumber(index_nd);
	for (unsigned ii=0; ii<nvt_ig; ii++) var_nd_f[ii+offset_nvt] = (*mysol)(ii);
	offset_nvt+=nvt_ig;
	delete mysol;
      }
     
      dimsf[0] = nvt;  dimsf[1] = 1;
      dataspace = H5Screate_simple(2,dimsf, NULL);
      dataset   = H5Dcreate(file_id,_ml_sol.GetSolutionName(indx),H5T_NATIVE_FLOAT,
			    dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status   = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_nd_f[0]);
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }
  //-------------------------------------------------------------------------------------------------------
    
  // Close the file -------------
  H5Fclose(file_id);
 
  //free memory
  delete [] filename;
  delete [] var_int;
  delete [] var_el_f;
  delete [] var_nd_f;
  
#endif
  
  return;   
}

void XDMFWriter::write_solution_wrapper(const char type[]) const {
  
#ifdef HAVE_HDF5 
  
  // to add--> _time_step0, _ntime_steps
  int time_step0 = 0;
  int ntime_steps = 1;
  int print_step = 1;
  
  char *filename= new char[60];
  // Print The Xdmf transient wrapper
  sprintf(filename,"./output/mesh.level%d.%s.xmf",_gridn,type);
  std::ofstream ftr_out;
  ftr_out.open(filename);
  if (!ftr_out) {
    std::cout << "Transient Output mesh file "<<filename<<" cannot be opened.\n";
    exit(0);
  }
  
  ftr_out << "<?xml version=\"1.0\" ?> \n";
  ftr_out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd []\">"<<std::endl;
  ftr_out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> " << std::endl;
  ftr_out << "<Domain> " << std::endl;
  ftr_out << "<Grid Name=\"Mesh\" GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
  // time loop for grid sequence
  for ( unsigned time_step = time_step0; time_step < time_step0 + ntime_steps; time_step++) {
    if ( !(time_step%print_step) ) {
      sprintf(filename,"./mesh.level%d.%d.%s.xmf",_gridn,time_step,type);
      ftr_out << "<xi:include href=\"" << filename << "\" xpointer=\"xpointer(//Xdmf/Domain/Grid["<< 1 <<"])\">\n";
      ftr_out << "<xi:fallback/>\n";
      ftr_out << "</xi:include>\n";
    }
  }
  ftr_out << "</Grid> \n";
  ftr_out << "</Domain> \n";
  ftr_out << "</Xdmf> \n";
  ftr_out.close();
  ftr_out.close();  
  //----------------------------------------------------------------------------------------------------------
  delete [] filename;


} //end namespace femus


  
#endif
 
}


