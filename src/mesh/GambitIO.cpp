/*=========================================================================

 Program: FEMUS
 Module: GambitIO
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//local include
#include "GambitIO.hpp"
#include "Mesh.hpp"

//C++ include
#include "cstdio"
#include "fstream"


namespace femus {

  
void GambitIO::read(const std::string& name, vector < vector < double> > &vt, const double Lref, std::vector<bool> &type_elem_flag) {
  
 // set the file
  const unsigned GambitVertexIndex[6][27]= {{
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

  const unsigned GambitFaceIndex[6][6]= {{0,4,2,5,3,1},
    {0,1,2,3},
    {2,1,0,4,3},
    {0,1,2,3},
    {0,1,2},
    {0,1}
  };
  
  mesh& mesh = GetMesh();

  std::ifstream inf;
  std::string str2;
  unsigned ngroup;
  unsigned nbcd;
  unsigned dim;
  double x,y,z;
  unsigned nvt;
  unsigned nel;

  mesh.SetGridNumber(0);
  // read control data ******************** A
  inf.open(name.c_str());
  if (!inf) {
    std::cout<<"Generic-mesh file "<< name << " can not read parameters\n";
    exit(0);
  }
  str2="0";
  while (str2.compare("NDFVL") != 0) inf >> str2;
  inf >> nvt >> nel >>  ngroup >> nbcd >> dim >> str2 ;
  mesh.SetDimension(dim);
  mesh.SetElementNumber(nel);
  mesh.SetNumberOfNodes(nvt);
  inf >> str2;
  if (str2.compare("ENDOFSECTION") != 0) {
    std::cout<<"error control data mesh"<<std::endl;
    exit(0);
  }
//   std::cout << "***************" << _dimension << std::endl;
  inf.close();
  // end read control data **************** A
  // read ELEMENT/cell ******************** B
  inf.open(name.c_str());
  if (!inf) {
    std::cout<<"Generic-mesh file "<< name << " cannot read elements\n";
    exit(0);
  }
  mesh.el= new elem(nel);
  while (str2.compare("ELEMENTS/CELLS") != 0) inf >> str2;
  inf >> str2;
  for (unsigned iel=0; iel<nel; iel++) {
    mesh.el->SetElementGroup(iel,1);
    unsigned nve;
    inf >> str2 >> str2 >> nve;
    if (nve==27) {
      type_elem_flag[0]=type_elem_flag[3]=true;
      mesh.el->AddToElementNumber(1,"Hex");
      mesh.el->SetElementType(iel,0);
    } else if (nve==10) {
      type_elem_flag[1]=type_elem_flag[4]=true;
      mesh.el->AddToElementNumber(1,"Tet");
      mesh.el->SetElementType(iel,1);
    } else if (nve==18) {
      type_elem_flag[2]=type_elem_flag[3]=type_elem_flag[4]=true;
      mesh.el->AddToElementNumber(1,"Wedge");
      mesh.el->SetElementType(iel,2);
    } else if (nve==9) {
      type_elem_flag[3]=true;
      mesh.el->AddToElementNumber(1,"Quad");
      mesh.el->SetElementType(iel,3);
    }
    else if (nve==6 && mesh.GetDimension()==2) {
      type_elem_flag[4]=true;
      mesh.el->AddToElementNumber(1,"Triangle");
      mesh.el->SetElementType(iel,4);
    }
    else if (nve==3 && mesh.GetDimension()==1) {
      mesh.el->AddToElementNumber(1,"Line");
      mesh.el->SetElementType(iel,5);
    } else {
      std::cout<<"Error! Invalid element type in reading Gambit File!"<<std::endl;
      std::cout<<"Error! Use a second order discretization"<<std::endl;
      exit(0);
    }
    for (unsigned i=0; i<nve; i++) {
      unsigned inode=GambitVertexIndex[mesh.el->GetElementType(iel)][i];
      double value;
      inf>>value;
      mesh.el->SetElementVertexIndex(iel,inode,value);
    }
  }
  inf >> str2;
  if (str2.compare("ENDOFSECTION") != 0) {
    std::cout<<"error element data mesh"<<std::endl;
    exit(0);
  }
  inf.close();

  // end read  ELEMENT/CELL **************** B

  // read NODAL COORDINATES **************** C
  inf.open(name.c_str());
  if (!inf) {
    std::cout<<"Generic-mesh file "<< name << " cannot read nodes\n";
    exit(0);
  }
  while (str2.compare("COORDINATES") != 0) inf >> str2;
  inf >> str2;  // 2.0.4
  vt[0].resize(nvt);
  vt[1].resize(nvt);
  vt[2].resize(nvt);

  if (mesh.GetDimension()==3) {
    for (unsigned j=0; j<nvt; j++) {
      inf >> str2 >> x >> y >> z;
      vt[0][j] = x/Lref;
      vt[1][j] = y/Lref;
      vt[2][j] = z/Lref;
    }
  }

  else if (mesh.GetDimension()==2) {
    for (unsigned j=0; j<nvt; j++) {
      inf >> str2 >> x >> y;
      vt[0][j] = x/Lref;
      vt[1][j] = y/Lref;
      vt[2][j] = 0.;
    }
  }

  else if (mesh.GetDimension()==1) {
    for (unsigned j=0; j<nvt; j++) {
      inf >> str2 >> x;
      vt[0][j] = x/Lref;
      vt[1][j]=0.;
      vt[2][j]=0.;
    }
  }
  inf >> str2; // "ENDOFSECTION"
  if (str2.compare("ENDOFSECTION") != 0) {
    std::cout<<"error node data mesh 1"<<std::endl;
    exit(0);
  }
  inf.close();
  // end read NODAL COORDINATES ************* C

  // read GROUP **************** E
  inf.open(name.c_str());
  if (!inf) {
    std::cout<<"Generic-mesh file "<< name << " cannot read group\n";
    exit(0);
  }
  mesh.el->SetElementGroupNumber(ngroup);
  for (unsigned k=0; k<ngroup; k++) {
    int ngel;
    int name;
    int mat;
    while (str2.compare("GROUP:") != 0) inf >> str2;
    inf >> str2 >> str2 >> ngel >> str2 >> mat >>str2 >> str2 >>name>> str2;
//     std::cout<<ngel<<std::endl;
//     std::cout<<name<<std::endl;
//     std::cout << mat << std::endl;
    for (int i=0; i<ngel; i++) {
      int iel;
      inf >> iel;
      mesh.el->SetElementGroup(iel-1,name);
      mesh.el->SetElementMaterial(iel-1,mat);
    }
    inf >> str2;
    if (str2.compare("ENDOFSECTION") != 0) {
      std::cout<<"error group data mesh"<<std::endl;
      exit(0);
    }
  }
  inf.close();
  // end read boundary **************** E

  // read boundary **************** D
  inf.open(name.c_str());
  if (!inf) {
    std::cout<<"Generic-mesh file "<< name << " cannot read boudary\n";
    exit(0);
  }
  for (unsigned k=0; k<nbcd; k++) {
    while (str2.compare("CONDITIONS") != 0) inf >> str2;
    inf >> str2;
    int value;
    unsigned nface;
    inf >>value>>str2>>nface>>str2>>str2;
    value=-value-1;
    for (unsigned i=0; i<nface; i++) {
      unsigned iel,iface;
      inf>>iel>>str2>>iface;
      iel--;
      iface=GambitFaceIndex[mesh.el->GetElementType(iel)][iface-1u];
      mesh.el->SetFaceElementIndex(iel,iface,value);
    }
    inf >> str2;
    if (str2.compare("ENDOFSECTION") != 0) {
      std::cout<<"error boundary data mesh"<<std::endl;
      exit(0);
    }
  }
  inf.close();
  // end read boundary **************** D

  
  
  //*************** start reorder mesh dofs **************
  //(1)linear (2)quadratic (3)biquaratic
  
  vector <unsigned> dof_index;
  dof_index.resize(nvt);
  for(unsigned i=0;i<nvt;i++){
    dof_index[i]=i+1;
  }
  //reorder vertices and mid-points vs central points
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode<mesh.el->GetElementDofNumber(iel,1); inode++) {
      for (unsigned jel=0; jel<nel; jel++) {
	for (unsigned jnode=mesh.el->GetElementDofNumber(jel,1); jnode<mesh.el->GetElementDofNumber(jel,3); jnode++) { 
	  unsigned ii=mesh.el->GetElementVertexIndex(iel,inode)-1;
	  unsigned jj=mesh.el->GetElementVertexIndex(jel,jnode)-1;
	  unsigned i0=dof_index[ii];
          unsigned i1=dof_index[jj];
	  if(i0>i1){
	    dof_index[ii]=i1;
	    dof_index[jj]=i0; 
	  }
	}
      }
    }
  }
  //reorder vertices vs mid-points
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode<mesh.el->GetElementDofNumber(iel,0); inode++) {
      for (unsigned jel=0; jel<nel; jel++) {
        for (unsigned jnode=mesh.el->GetElementDofNumber(jel,0); jnode<mesh.el->GetElementDofNumber(jel,1); jnode++) {
          unsigned ii=mesh.el->GetElementVertexIndex(iel,inode)-1;
	  unsigned jj=mesh.el->GetElementVertexIndex(jel,jnode)-1;
	  unsigned i0=dof_index[ii];
          unsigned i1=dof_index[jj];
	  if(i0>i1){
	    dof_index[ii]=i1;
	    dof_index[jj]=i0; 
	  }
	}
      }
    }
  }
  
  // update all
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode<mesh.el->GetElementDofNumber(iel,3); inode++) {
      unsigned ii=mesh.el->GetElementVertexIndex(iel,inode)-1;
      mesh.el->SetElementVertexIndex(iel,inode,dof_index[ii]);
    }
  }
  vector <double> vt_temp;
  for(int i=0;i<3;i++){
    vt_temp=vt[i];
    for(unsigned j=0;j<nvt;j++){
      vt[i][dof_index[j]-1]=vt_temp[j];
    }
  }
  // **************  end reoreder mesh dofs **************
 
  mesh.el->SetNodeNumber(nvt);

  unsigned nv0=0;
  for (unsigned iel=0; iel<nel; iel++)
    for (unsigned inode=0; inode<mesh.el->GetElementDofNumber(iel,0); inode++) {
      unsigned i0=mesh.el->GetElementVertexIndex(iel,inode);
      if (nv0<i0) nv0=i0;
  }
  mesh.el->SetVertexNodeNumber(nv0);

  unsigned nv1=0;
  for (unsigned iel=0; iel<nel; iel++)
    for (unsigned inode=mesh.el->GetElementDofNumber(iel,0); inode<mesh.el->GetElementDofNumber(iel,1); inode++) {
      unsigned i1=mesh.el->GetElementVertexIndex(iel,inode);
      if (nv1<i1) nv1=i1;
  }
  mesh.el->SetMidpointNodeNumber(nv1-nv0);

  mesh.el->SetCentralNodeNumber(nvt-nv1);
  
};

}