/*=========================================================================

 Program: FEMUS
 Module: SalomeIO
 Authors: Sureka Pathmanathan, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//local include
#include "SalomeIO.hpp"
#include "Mesh.hpp"

//C++ include
#include <cassert>
#include <cstdio>
#include <fstream>


namespace femus {
  
   const std::string SalomeIO::group_name_begin = "FAS";
   const std::string SalomeIO::group_name_end   = "ELEME";
   const std::string SalomeIO::mesh_ensemble = "ENS_MAA";
   const std::string SalomeIO::aux_zeroone = "-0000000000000000001-0000000000000000001";
   const std::string SalomeIO::connectivity = "MAI";
   const std::string SalomeIO::node_coord = "NOE/COO";
   const uint SalomeIO::max_length = 100;  ///@todo this length of the menu string is conservative enough...

  
  
 const unsigned SalomeIO::GambitToFemusVertexIndex[6][27]= 
   {
    {
      4,16,0,15,23,11,7,19,3,
      12,20,8,25,26,24,14,22,10,
      5,17,1,13,21,9,6,18,2
    },
    {
      0,4,1,6,5,
      2,7,8,9,3
    },
    {
      3, 11,5, 9, 10,4,
      12,17,14,15,16,13,
      0, 8, 2, 6, 7, 1
    },
    {0,4,1,5,2,6,3,7,8},
    {0,3,1,4,2,5},
    {0,2,1}
  };

  
const unsigned SalomeIO::GambitToFemusFaceIndex[6][6]= 
  {
    {0,4,2,5,3,1},
    {0,1,2,3},
    {2,1,0,4,3},
    {0,1,2,3},
    {0,1,2},
    {0,1}
  };

  
void SalomeIO::read(const std::string& name, vector < vector < double> > &coords, const double Lref, std::vector<bool> &type_elem_flag) {
   
    Mesh& mesh = GetMesh();
    mesh.SetGridNumber(0);

    hsize_t dims[2];
    std::string el_fem_type_vol(""); 
    std::string el_fem_type_bd("");
    
    hid_t  file_id = H5Fopen(name.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
    
    hid_t  gid = H5Gopen(file_id,mesh_ensemble.c_str(),H5P_DEFAULT);
      
    hsize_t     n_menus;
    hid_t status= H5Gget_num_objs(gid, &n_menus);  // number of menus
    if(status !=0) { std::cout << "Number of mesh menus not found"; abort(); }
    
    std::vector<int>  itype_vol;
     itype_vol.resize(n_menus);
    menu_names.resize(n_menus);
    for (unsigned j=0; j<menu_names.size(); j++) {
      menu_names[j] = new char[max_length];
     H5Gget_objname_by_idx(gid,j,menu_names[j],max_length); ///@deprecated see the HDF doc to replace this
     std::string tempj(menu_names[j]);
     itype_vol[j] =  ReadFE(file_id,el_fem_type_vol,el_fem_type_bd,tempj);
    }

    //   // read control data ******************** A
//   mesh.SetDimension(dim);  // this is determined in the other routine already
//   mesh.SetNumberOfElements(nel);
//   mesh.SetNumberOfNodes(nvt);

    
    
    //   // end read control data **************** A

    
    
    
    
    
    
    //   // read ELEMENT/cell ******************** B

    
    //   // end read  ELEMENT/CELL **************** B

    
    
    
    
    //   // read NODAL COORDINATES **************** C

    //   // end read NODAL COORDINATES ************* C

    
    
    
    
//   // read GROUP **************** E
    
//   // end read GROUP **************** E






//   // read boundary **************** D

//   // end read boundary **************** D


    for (unsigned j=0; j<menu_names.size(); j++) { delete [] menu_names[j]; }
    
    status = H5Fclose(file_id);
 
    exit(0);   // just to make the test exit successfully

// 
//   std::ifstream inf;
//   std::string str2;
//   unsigned ngroup;
//   unsigned nbcd;
//   unsigned dim;
//   double x,y,z;
//   unsigned nvt;
//   unsigned nel;
// 
//   mesh.SetGridNumber(0);
//   
//   
//   
//   // read control data ******************** A
//   inf.open(name.c_str());
//   if (!inf) {
//     std::cout<<"Generic-mesh file "<< name << " can not read parameters\n";
//     exit(0);
//   }
//   str2="0";
//   while (str2.compare("NDFVL") != 0) inf >> str2;
//   inf >> nvt >> nel >>  ngroup >> nbcd >> dim >> str2 ;
//   mesh.SetDimension(dim);
//   mesh.SetNumberOfElements(nel);
//   mesh.SetNumberOfNodes(nvt);
//   inf >> str2;
//   if (str2.compare("ENDOFSECTION") != 0) {
//     std::cout<<"error control data mesh"<<std::endl;
//     exit(0);
//   }
// //   std::cout << "***************" << _dimension << std::endl;
//   inf.close();
//   // end read control data **************** A
//   
//   
//   
//   
//   // read ELEMENT/cell ******************** B
//   inf.open(name.c_str());
//   if (!inf) {
//     std::cout<<"Generic-mesh file "<< name << " cannot read elements\n";
//     exit(0);
//   }
//   mesh.el= new elem(nel);
//   while (str2.compare("ELEMENTS/CELLS") != 0) inf >> str2;
//   inf >> str2;
//   for (unsigned iel=0; iel<nel; iel++) {
//     mesh.el->SetElementGroup(iel,1);
//     unsigned nve;
//     inf >> str2 >> str2 >> nve;
//     if (nve==27) {
//       type_elem_flag[0]=type_elem_flag[3]=true;
//       mesh.el->AddToElementNumber(1,"Hex");
//       mesh.el->SetElementType(iel,0);
//     } else if (nve==10) {
//       type_elem_flag[1]=type_elem_flag[4]=true;
//       mesh.el->AddToElementNumber(1,"Tet");
//       mesh.el->SetElementType(iel,1);
//     } else if (nve==18) {
//       type_elem_flag[2]=type_elem_flag[3]=type_elem_flag[4]=true;
//       mesh.el->AddToElementNumber(1,"Wedge");
//       mesh.el->SetElementType(iel,2);
//     } else if (nve==9) {
//       type_elem_flag[3]=true;
//       mesh.el->AddToElementNumber(1,"Quad");
//       mesh.el->SetElementType(iel,3);
//     }
//     else if (nve==6 && mesh.GetDimension()==2) {
//       type_elem_flag[4]=true;
//       mesh.el->AddToElementNumber(1,"Triangle");
//       mesh.el->SetElementType(iel,4);
//     }
//     else if (nve==3 && mesh.GetDimension()==1) {
//       mesh.el->AddToElementNumber(1,"Line");
//       mesh.el->SetElementType(iel,5);
//     } else {
//       std::cout<<"Error! Invalid element type in reading Gambit File!"<<std::endl;
//       std::cout<<"Error! Use a second order discretization"<<std::endl;
//       exit(0);
//     }
//     for (unsigned i=0; i<nve; i++) {
//       unsigned inode=SalomeIO::GambitToFemusVertexIndex[mesh.el->GetElementType(iel)][i];
//       double value;
//       inf>>value;
//       mesh.el->SetElementVertexIndex(iel,inode,value);
//     }
//   }
//   inf >> str2;
//   if (str2.compare("ENDOFSECTION") != 0) {
//     std::cout<<"error element data mesh"<<std::endl;
//     exit(0);
//   }
//   inf.close();
// 
//   // end read  ELEMENT/CELL **************** B
//   
//   
//   
// 
//   // read NODAL COORDINATES **************** C
//   inf.open(name.c_str());
//   if (!inf) {
//     std::cout<<"Generic-mesh file "<< name << " cannot read nodes\n";
//     exit(0);
//   }
//   while (str2.compare("COORDINATES") != 0) inf >> str2;
//   inf >> str2;  // 2.0.4
//   coords[0].resize(nvt);
//   coords[1].resize(nvt);
//   coords[2].resize(nvt);
// 
//   if (mesh.GetDimension()==3) {
//     for (unsigned j=0; j<nvt; j++) {
//       inf >> str2 >> x >> y >> z;
//       coords[0][j] = x/Lref;
//       coords[1][j] = y/Lref;
//       coords[2][j] = z/Lref;
//     }
//   }
// 
//   else if (mesh.GetDimension()==2) {
//     for (unsigned j=0; j<nvt; j++) {
//       inf >> str2 >> x >> y;
//       coords[0][j] = x/Lref;
//       coords[1][j] = y/Lref;
//       coords[2][j] = 0.;
//     }
//   }
// 
//   else if (mesh.GetDimension()==1) {
//     for (unsigned j=0; j<nvt; j++) {
//       inf >> str2 >> x;
//       coords[0][j] = x/Lref;
//       coords[1][j]=0.;
//       coords[2][j]=0.;
//     }
//   }
//   inf >> str2; // "ENDOFSECTION"
//   if (str2.compare("ENDOFSECTION") != 0) {
//     std::cout<<"error node data mesh 1"<<std::endl;
//     exit(0);
//   }
//   inf.close();
//   // end read NODAL COORDINATES ************* C
//   
//   
//   
// 
//   // read GROUP **************** E
//   inf.open(name.c_str());
//   if (!inf) {
//     std::cout<<"Generic-mesh file "<< name << " cannot read group\n";
//     exit(0);
//   }
//   mesh.el->SetElementGroupNumber(ngroup);
//   for (unsigned k=0; k<ngroup; k++) {
//     int ngel;
//     int name;
//     int mat;
//     while (str2.compare("GROUP:") != 0) inf >> str2;
//     inf >> str2 >> str2 >> ngel >> str2 >> mat >>str2 >> str2 >>name>> str2;
//     for (int i=0; i<ngel; i++) {
//       int iel;
//       inf >> iel;
//       mesh.el->SetElementGroup(iel-1,name);
//       mesh.el->SetElementMaterial(iel-1,mat);
//     }
//     inf >> str2;
//     if (str2.compare("ENDOFSECTION") != 0) {
//       std::cout<<"error group data mesh"<<std::endl;
//       exit(0);
//     }
//   }
//   inf.close();
//   // end read boundary **************** E
//   
//   
// 
//   // read boundary **************** D
//   inf.open(name.c_str());
//   if (!inf) {
//     std::cout<<"Generic-mesh file "<< name << " cannot read boudary\n";
//     exit(0);
//   }
//   for (unsigned k=0; k<nbcd; k++) {
//     while (str2.compare("CONDITIONS") != 0) inf >> str2;
//     inf >> str2;
//     int value;
//     unsigned nface;
//     inf >>value>>str2>>nface>>str2>>str2;
//     value=-value-1;
//     for (unsigned i=0; i<nface; i++) {
//       unsigned iel,iface;
//       inf>>iel>>str2>>iface;
//       iel--;
//       iface=SalomeIO::GambitToFemusFaceIndex[mesh.el->GetElementType(iel)][iface-1u];
//       mesh.el->SetFaceElementIndex(iel,iface,value);
//     }
//     inf >> str2;
//     if (str2.compare("ENDOFSECTION") != 0) {
//       std::cout<<"error boundary data mesh"<<std::endl;
//       exit(0);
//     }
//   }
//   inf.close();
//   // end read boundary **************** D
 
}




int  SalomeIO::ReadFE(
  hid_t file_id,
  std::string & el_fem_type_vol,
  std::string & el_fem_type_bd,
  const  std::string menu_name
) {

    Mesh& mesh = GetMesh();
    uint mydim = 1;  //this is the initial value, then it will be updated below
    mesh.SetDimension(mydim);

  // ************  determine mesh dimension  ************************************
    /// @todo this determination of the dimension from the mesh file would not work with a 2D mesh embedded in 3D
  std::string my_mesh_name_dir = mesh_ensemble +  "/" + menu_name + "/" +  aux_zeroone + "/" + connectivity + "/";  ///@todo here we have to loop
  
  hid_t       gid=H5Gopen(file_id,my_mesh_name_dir.c_str(),H5P_DEFAULT); // group identity
  hsize_t     n_fem_type;
  hid_t status= H5Gget_num_objs(gid, &n_fem_type);  // number of links
  if(status !=0) {std::cout << "SalomeIO::read_fem_type:   H5Gget_num_objs not found"; abort();}

  std::vector<char*> elem_list;  elem_list.resize(n_fem_type);
    for (unsigned j=0; j < elem_list.size(); j++) {
      elem_list[j] = new char[max_length];
      H5Gget_objname_by_idx(gid,j,elem_list[j],max_length); ///@deprecated see the HDF doc to replace this
      std::string tempj(elem_list[j]);
      
      if ( tempj.compare("HE8") == 0 || 
	   tempj.compare("H20") == 0 ||  
	   tempj.compare("H27") == 0 || 
	   tempj.compare("TE4") == 0 || 
	   tempj.compare("T10") == 0    )      { mydim = 3; } 
      else if ( tempj.compare("QU4") == 0 || 
	        tempj.compare("QU8") == 0 ||  
	        tempj.compare("QU9") == 0 || 
	        tempj.compare("TR3") == 0 || 
	        tempj.compare("TR6") == 0    ) { mydim = 2; } 
      
      if ( mydim > mesh.GetDimension() ) mesh.SetDimension(mydim); 
	
    }  //end for
    
  // ---------  Reading Element type (From salome to LIBMESH (itype_vol))------
  std::map< std::string, int > fem_type_vol;// Salome fem table name (vol)
  fem_type_vol["HE8"] = 5;  fem_type_vol["H20"] = 20; fem_type_vol["H27"] = 12;
  fem_type_vol["TE4"] = 4;  fem_type_vol["T10"] = 11;
  fem_type_vol["QU4"] = 3;  fem_type_vol["QU8"] = 19; fem_type_vol["QU9"] = 10;
  fem_type_vol["TR3"] = 2;  fem_type_vol["TR6"] = 9;
  fem_type_vol["SE2"] = 0;  fem_type_vol["SE3"] = 0;  // not valid in 3D
  if(mesh.GetDimension()==3) {    // not valid in 3D as boundary
    fem_type_vol["QU4"] = 0; fem_type_vol["QU8"] = 0; fem_type_vol["QU9"] = 0;
    fem_type_vol["TR3"] = 0; fem_type_vol["TR6"] = 0;
  }
  std::map< std::string, int > fem_type_bd; // Salome fem table name (surface)
  fem_type_bd["HE8"] = 0;  fem_type_bd["H20"] = 0;  fem_type_bd["H27"] = 0;
  fem_type_bd["TE4"] = 0;  fem_type_bd["T10"] = 0;  // not valid in 2D
  fem_type_bd["QU4"] = 3;  fem_type_bd["QU8"] = 19;  fem_type_bd["QU9"] = 10;
  fem_type_bd["TR3"] = 2;  fem_type_bd["TR6"] = 9;
  fem_type_bd["SE2"] = 1;  fem_type_bd["SE3"] = 8;
  if(mesh.GetDimension()==3) {    // not valid in 3D as boundary
    fem_type_bd["SE2"] = 0;  fem_type_bd["SE3"] = 0;
  }


  // Get the element name from MESH_NAME_DIR in the file med (el_name)
  char **el_fem_type=new char*[n_fem_type];
  int index_vol=0;  int index_bd=0;
  for(int i=0; i<(int)n_fem_type; i++) {
    el_fem_type[i]=new char[4];
    H5Lget_name_by_idx(file_id,my_mesh_name_dir.c_str(), H5_INDEX_NAME, H5_ITER_INC,i,el_fem_type[i],4, H5P_DEFAULT);
    if(fem_type_vol[el_fem_type[i]]!=0) index_vol=i;
    if(fem_type_bd[el_fem_type[i]]!=0) index_bd=i;
  }

  // LIBMESH volume fem type (itype_vol) and MEDfem (el_fem_type_vol-el_fem_type_bd)
  int itype_vol= fem_type_vol[el_fem_type[index_vol]]; assert(itype_vol!=0);
  el_fem_type_vol=el_fem_type[index_vol];
  el_fem_type_bd=el_fem_type[index_bd];
  // clean
  for(int i=0; i<(int)n_fem_type; i++) delete[] el_fem_type[i];
  delete[] el_fem_type; fem_type_vol.clear(); fem_type_bd.clear();
  return itype_vol;
}



} //end namespace femus

