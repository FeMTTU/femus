/*=========================================================================

 Program: FEMUS
 Module: XDMFOutput
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
#include "XDMFOutput.hpp"

XDMFOutput::XDMFOutput(MultiLevelProblem& ml_probl): Output(ml_probl)
{
  
}

XDMFOutput::~XDMFOutput()
{
  
}

void XDMFOutput::write_system_solutions(const char order[], std::vector<std::string>& vars, const unsigned time_step) 
{ 
//     // to add time_step
//   int time_step = 0;
//   
//   if(_iproc!=0) return;
//   
//   bool test_all=!(vars[0].compare("All"));
//     
//   unsigned index=0;
//   unsigned index_nd=0;
//   if(!strcmp(type,"linear")) {    //linear
//     index=0;
//     index_nd=0;
//   }
//   else if(!strcmp(type,"quadratic")) {  //quadratic
//     index=1;
//     index_nd=1;
//   }
//   else if(!strcmp(type,"biquadratic")) { //biquadratic
//     index=3;
//     index_nd=2;
//   }
// 
//   const string type_el[4][6] = {{"Hexahedron","Tetrahedron","Wedge","Quadrilateral","Triangle","Edge"},
//                                 {"Hexahedron_20","Tetrahedron_10","Not_implemented","Quadrilateral_8","Triangle_6","Edge_3"},
// 			        {"Not_implemented","Not_implemented","Not_implemented","Not_implemented",
// 				 "Not_implemented","Not_implemented"},
//                                 {"Not_implemented","Not_implemented","Not_implemented","Quadrilateral_9",
// 				 "Not_implemented","Not_implemented"}};
// 			 
//   
//   //I assume that the mesh is not mixed
//   string type_elem;
//   type_elem = type_el[index][_msh[_gridn-1u]->el->GetElementType(0)];
//   
//   if (type_elem.compare("Not_implemented") == 0) exit(1);
//   
//   unsigned nvt=0;
//   for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     unsigned nvt_ig=_msh[ig]->GetDofNumber(index_nd);
//     nvt+=nvt_ig;
//   } 
//   
//   // Printing connectivity
//   unsigned nel=0;
//   for(unsigned ig=0;ig<_gridn-1u;ig++) {
//     nel+=( _msh[ig]->GetElementNumber() - _msh[ig]->el->GetRefinedElementNumber());
//   }
//   nel+=_msh[_gridn-1u]->GetElementNumber();
//   
//   unsigned icount;
//   unsigned el_dof_number  = _msh[_gridn-1u]->el->GetElementDofNumber(0,index);
//   int *var_int             = new int [nel*el_dof_number];
//   float *var_el_f         = new float [nel];
//   float *var_nd_f         = new float [nvt];
// 
//   char *filename= new char[60];
//   std::ofstream fout;
//   
//   //--------------------------------------------------------------------------------------------------
//   // Print The Xdmf wrapper
//   sprintf(filename,"./output/mesh.level%d.%d.%s.xmf",_gridn,time_step,type);
//   //std::ofstream fout;
//   fout.open(filename);
//   if (!fout) {
//     cout << "Output mesh file "<<filename<<" cannot be opened.\n";
//     exit(0);
//   }
//   
//   // Print The HDF5 file
//   sprintf(filename,"./mesh.level%d.%d.%s.h5",_gridn,time_step,type);
//   // haed ************************************************
//   fout<<"<?xml version=\"1.0\" ?>" << endl;
//   fout<<"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd []\">"<<endl;
//   fout<<"<Xdmf>"<<endl;
//   fout<<"<Domain>"<<endl;
//   fout<<"<Grid Name=\"Mesh\">"<<endl;
//   fout<<"<Time Value =\""<< time_step<< "\" />"<<endl;
//   fout<<"<Topology TopologyType=\""<< type_elem <<"\" NumberOfElements=\""<< nel <<"\">"<<endl;
//   //Connectivity
//   fout<<"<DataStructure DataType=\"Int\" Dimensions=\""<< nel*el_dof_number <<"\"" << "  Format=\"HDF\">" << endl;
//   fout << filename << ":CONNECTIVITY" << endl;
//   fout <<"</DataStructure>" << endl;
//   fout << "</Topology>" << endl;
//   fout << "<Geometry Type=\"X_Y_Z\">" << endl;
//   //Node_X
//   fout<<"<DataStructure DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << endl;
//   fout << filename << ":NODES_X1" << endl;
//   fout <<"</DataStructure>" << endl;
//   //Node_Y
//   fout<<"<DataStructure DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << endl;
//   fout << filename << ":NODES_X2" << endl;
//   fout <<"</DataStructure>" << endl;
//   //Node_Z
//   fout<<"<DataStructure DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << endl;
//   fout << filename << ":NODES_X3" << endl;
//   fout <<"</DataStructure>" << endl;
//   fout <<"</Geometry>" << endl;
//   //Regions
//   fout << "<Attribute Name=\""<< "Regions"<<"\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl;
//   fout << "<DataItem DataType=\"Int\" Dimensions=\""<< nel << "\"" << "  Format=\"HDF\">" << endl;
//   fout << filename << ":REGIONS" << endl;
//   fout << "</DataItem>" << endl;
//   fout << "</Attribute>" << endl;
//   // Solution Variables
//   for (unsigned i=0; i<vars.size(); i++) {
//     unsigned indx=GetIndex(vars[i].c_str());  
//     //Printing biquadratic solution on the nodes
//     if(SolType[indx]<3) {  
//       fout << "<Attribute Name=\""<< SolName[indx]<<"\" AttributeType=\"Scalar\" Center=\"Node\">" << endl;
//       fout << "<DataItem DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << endl;
//       fout << filename << ":" << SolName[indx] << endl;
//       fout << "</DataItem>" << endl;
//       fout << "</Attribute>" << endl;
//     }
//     else if (SolType[indx]>=3) {  //Printing picewise constant solution on the element
//       fout << "<Attribute Name=\""<< SolName[indx]<<"\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl;
//       fout << "<DataItem DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nel << "\"  Format=\"HDF\">" << endl;
//       fout << filename << ":" << SolName[indx] << endl;
//       fout << "</DataItem>" << endl;
//       fout << "</Attribute>" << endl;
//     }
//   }
// 
//   fout <<"</Grid>" << endl;
//   fout <<"</Domain>" << endl;
//   fout <<"</Xdmf>" << endl;
//   fout.close();
//   //----------------------------------------------------------------------------------------------------------
//   
//   //----------------------------------------------------------------------------------------------------------
//   hid_t file_id;
//   sprintf(filename,"./output/mesh.level%d.%d.%s.h5",_gridn,time_step,type);
//   file_id = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
//   hsize_t dimsf[2];
//   herr_t status;
//   hid_t dataspace;
//   hid_t dataset;
//   
//   //-----------------------------------------------------------------------------------------------------------
//   // Printing nodes coordinates 
//   
//   PetscScalar *MYSOL[1]; //TODO
//   unsigned varind[3];
//   varind[0]=GetIndex("X");
//   varind[1]=GetIndex("Y");
//   varind[2]=GetIndex("Z");
// 
//   
//   for (int i=0; i<3; i++) {
//     unsigned offset_nvt=0;
//     for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//       NumericVector* mysol;
//       mysol = NumericVector::build().release();
//       mysol->init(_msh[ig]->GetDofNumber(index_nd),_msh[ig]->GetDofNumber(index_nd),true,SERIAL);
//       mysol->matrix_mult(*_solution[ig]->_Sol[varind[i]],*ProlQitoQj_[index_nd][SolType[varind[i]]][ig]);
//       unsigned nvt_ig=_msh[ig]->GetDofNumber(index_nd);
//       for (unsigned ii=0; ii<nvt_ig; ii++) var_nd_f[ii+offset_nvt] = (*mysol)(ii);
//       if (_moving_mesh) {
// 	unsigned varind_DXDYDZ=GetIndex(_moving_vars[i].c_str());
// 	mysol->matrix_mult(*_solution[ig]->_Sol[varind_DXDYDZ],*ProlQitoQj_[index_nd][SolType[varind_DXDYDZ]][ig]);
// 	for (unsigned ii=0; ii<nvt_ig; ii++) var_nd_f[ii+offset_nvt] += (*mysol)(ii);
//       }
//       offset_nvt+=nvt_ig;
//       delete mysol;
//     }
//     
//     dimsf[0] = nvt ;  dimsf[1] = 1;
//     std::ostringstream Name; Name << "/NODES_X" << i+1;
//     dataspace = H5Screate_simple(2,dimsf, NULL);
//     dataset   = H5Dcreate(file_id,Name.str().c_str(),H5T_NATIVE_FLOAT,
// 			  dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//     status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_nd_f[0]);
//     H5Sclose(dataspace);
//     H5Dclose(dataset);
//     
//   }
//   
//   //-------------------------------------------------------------------------------------------------------------
// 
//   //------------------------------------------------------------------------------------------------------
//   //connectivity
//   icount = 0;
//   unsigned offset_conn=0;
//   for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     for (unsigned iel=0; iel<_msh[ig]->GetElementNumber(); iel++) {
//       if (_msh[ig]->el->GetRefinedElementIndex(iel)==0 || ig==_gridn-1u) {
//         for (unsigned j=0; j<_msh[ig]->el->GetElementDofNumber(iel,index); j++) {
// 	  unsigned jnode=_msh[ig]->el->GetElementVertexIndex(iel,j)-1u;
// 	  unsigned jnode_Metis = _msh[ig]->GetMetisDof(jnode,index_nd);
// 	  var_int[icount] = offset_conn + jnode_Metis;
// 	  icount++;
// 	}
//       }
//     }
//     offset_conn += _msh[ig]->GetDofNumber(index_nd);
//   }
//   
//   dimsf[0] = nel*el_dof_number ;  dimsf[1] = 1;
//   dataspace = H5Screate_simple(2,dimsf, NULL);
//   dataset   = H5Dcreate(file_id,"/CONNECTIVITY",H5T_NATIVE_INT,
// 			dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//   status   = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_int[0]);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//   //------------------------------------------------------------------------------------------------------
//   
//   
//   //-------------------------------------------------------------------------------------------------------
//   // print regions
//   icount=0;
//   for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//     for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
//       if (ig==_gridn-1u || 0==_msh[ig]->el->GetRefinedElementIndex(ii)) {
// 	unsigned iel_Metis = _msh[ig]->GetMetisDof(ii,3);
// 	var_int[icount] = _msh[ig]->el->GetElementGroup(iel_Metis);
// 	icount++;
//       }
//     }
//   } 
//    
//   dimsf[0] = nel;  dimsf[1] = 1;
//   dataspace = H5Screate_simple(2,dimsf, NULL);
//   dataset   = H5Dcreate(file_id,"/REGIONS",H5T_NATIVE_INT,
// 			dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//   status   = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_int[0]);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//   
//   //-------------------------------------------------------------------------------------------------------
//   // printing element variables
//   for (unsigned i=0; i<(1-test_all)*vars.size()+test_all*SolName.size(); i++) {
//     unsigned indx=(test_all==0)?GetIndex(vars[i].c_str()):i;
//     if (SolType[indx]>=3) {
//       icount=0;
//       for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
// 	for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
// 	  if (ig==_gridn-1u || 0==_msh[ig]->el->GetRefinedElementIndex(ii)) {
// 	    unsigned iel_Metis = _msh[ig]->GetMetisDof(ii,SolType[indx]);
// 	    var_el_f[icount]=(*_solution[ig]->_Sol[indx])(iel_Metis);
// 	    icount++;
// 	  }
// 	}
//       } 
//      
//       dimsf[0] = nel;  dimsf[1] = 1;
//       dataspace = H5Screate_simple(2,dimsf, NULL);
//       dataset   = H5Dcreate(file_id,SolName[indx],H5T_NATIVE_FLOAT,
// 			    dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//       status   = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_el_f[0]);
//       H5Sclose(dataspace);
//       H5Dclose(dataset);
//      
//     }
//   }
//   
//   //-------------------------------------------------------------------------------------------------------
//   // printing nodes variables
//   for (unsigned i=0; i<(1-test_all)*vars.size()+test_all*SolName.size(); i++) {
//     unsigned indx=(test_all==0)?GetIndex(vars[i].c_str()):i;
//     if (SolType[indx] < 3) {
//       unsigned offset_nvt=0;
//       for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
//         NumericVector* mysol;
// 	mysol = NumericVector::build().release();
//         mysol->init(_msh[ig]->GetDofNumber(index_nd),_msh[ig]->GetDofNumber(index_nd),true,SERIAL);
// 	mysol->matrix_mult(*_solution[ig]->_Sol[indx],*ProlQitoQj_[index_nd][SolType[indx]][ig]);
// 	unsigned nvt_ig=_msh[ig]->GetDofNumber(index_nd);
// 	for (unsigned ii=0; ii<nvt_ig; ii++) var_nd_f[ii+offset_nvt] = (*mysol)(ii);
// 	offset_nvt+=nvt_ig;
// 	delete mysol;
//       }
//      
//       dimsf[0] = nvt;  dimsf[1] = 1;
//       dataspace = H5Screate_simple(2,dimsf, NULL);
//       dataset   = H5Dcreate(file_id,SolName[indx],H5T_NATIVE_FLOAT,
// 			    dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//       status   = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_nd_f[0]);
//       H5Sclose(dataspace);
//       H5Dclose(dataset);
//     }
//   }
//   //-------------------------------------------------------------------------------------------------------
//     
//   // Close the file -------------
//   H5Fclose(file_id);
//  
//   //free memory
//   delete [] filename;
//   delete [] var_int;
//   delete [] var_el_f;
//   delete [] var_nd_f;
  
  return;   
}

void XDMFOutput::write_solution_wrapper(const char type[]) const {
  
//   // to add--> _time_step0, _ntime_steps
//   int time_step0 = 0;
//   int ntime_steps = 1;
//   int print_step = 1;
//   
//   char *filename= new char[60];
//   // Print The Xdmf transient wrapper
//   sprintf(filename,"./output/mesh.level%d.%s.xmf",_gridn,type);
//   std::ofstream ftr_out;
//   ftr_out.open(filename);
//   if (!ftr_out) {
//     cout << "Transient Output mesh file "<<filename<<" cannot be opened.\n";
//     exit(0);
//   }
//   
//   ftr_out << "<?xml version=\"1.0\" ?> \n";
//   ftr_out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd []\">"<<endl;
//   ftr_out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> " << endl;
//   ftr_out << "<Domain> " << endl;
//   ftr_out << "<Grid Name=\"Mesh\" GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
//   // time loop for grid sequence
//   for ( unsigned time_step = time_step0; time_step < time_step0 + ntime_steps; time_step++) {
//     if ( !(time_step%print_step) ) {
//       sprintf(filename,"./mesh.level%d.%d.%s.xmf",_gridn,time_step,type);
//       ftr_out << "<xi:include href=\"" << filename << "\" xpointer=\"xpointer(//Xdmf/Domain/Grid["<< 1 <<"])\">\n";
//       ftr_out << "<xi:fallback/>\n";
//       ftr_out << "</xi:include>\n";
//     }
//   }
//   ftr_out << "</Grid> \n";
//   ftr_out << "</Domain> \n";
//   ftr_out << "</Xdmf> \n";
//   ftr_out.close();
//   ftr_out.close();  
//   //----------------------------------------------------------------------------------------------------------
//   delete [] filename;
 
}


