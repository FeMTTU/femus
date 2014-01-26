#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include "mpi.h"
#include <algorithm>
using std::cout;
using std::endl;
using std::min;

#include "Mesh.hpp"
#include "metis.h"
using std::sort;

const unsigned mesh::_END_IND[5]= {0,1,3,4,5};

unsigned mesh::_dimension=2;
unsigned mesh::_ref_index=4;  // 8*DIM[2]+4*DIM[1]+2*DIM[0];
unsigned mesh::_face_index=2; // 4*DIM[2]+2*DIM[1]+1*DIM[0];


/**
 *  This constructur generates the coarse mesh level, $l_0$, from the gambit data file
 **/
mesh::mesh(const char infile[], vector < vector < double> > &vt, const double Lref) {
  
  MPI_Comm_rank(MPI_COMM_WORLD, &_iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &_nprocs);
  
  grid=0;
  // mesh
  //now we are reading only 2D or 3D
  //   if(DIM[0]==1)
  //     Read1D(infile,vt);
  //   else
  ReadGambit(infile,vt,Lref);

  // connectivity: find all the element near the vertices
  BuildAdjVtx();
  Buildkel();
  // some output
  cout <<"grid\tnel\tnvt"<< endl;
  cout <<grid<<"\t"<<nel<<"\t"<<nvt<<endl;
  cout <<"nv0\tnv1\tnv2"<< endl;
  cout <<el->GetVertexNodeNumber()<<"\t"
       <<el->GetMidpointNodeNumber()<<"\t"
       <<el->GetCentralNodeNumber()<<endl
       <<"-------------------------"<<endl;
       
 if (_nprocs>=1) generate_metis_mesh_partition();
  vector <double> vt_temp;
  for(int i=0;i<3;i++){
    vt_temp=vt[i];
    for(unsigned j=0;j<nvt;j++) {
      vt[i][GetMetisDof(j,2)]=vt_temp[j];
    }
  }
};

//------------------------------------------------------------------------------------------------------
void mesh::generate_metis_mesh_partition(){
  
  int *aux_vec = new int[nvt];
  for(unsigned i=0; i<nvt; i++) aux_vec[i]=-1;
  
  unsigned eind_size = el->GetElementNumber("Hex")*NVE[0][3] + el->GetElementNumber("Tet")*NVE[1][3] 
                     + el->GetElementNumber("Wedge")*NVE[2][3] + el->GetElementNumber("Quad")*NVE[3][3] 
                     + el->GetElementNumber("Triangle")*NVE[4][3] + el->GetElementNumber("Line")*NVE[5][3];

  idx_t *eptr=new idx_t [nel+1];
  idx_t *eind=new idx_t [eind_size];
  
  epart=new idx_t [nel];
  npart=new idx_t [nvt];
  
  idx_t objval;
  
  idx_t options[METIS_NOPTIONS]; 
  
  METIS_SetDefaultOptions(options);
  //for(int k=0;k<METIS_NOPTIONS;k++) cout<<options[k]<<endl;
  
  options[METIS_OPTION_NUMBERING]=0;
  options[METIS_OPTION_DBGLVL]   =0;
  
 
//   options[METIS_OPTION_CCORDER]  =0;   
//   options[METIS_OPTION_CTYPE]    = METIS_CTYPE_RM; 
//   options[METIS_OPTION_IPTYPE]   = METIS_IPTYPE_RANDOM;
//   METIS_SetDefaultOptions(options);
//   options[METIS_OPTION_PTYPE]   = params->ptype;
//   options[METIS_OPTION_OBJTYPE] = params->objtype;
//   options[METIS_OPTION_CTYPE]   = params->ctype;
//   options[METIS_OPTION_IPTYPE]  = params->iptype;
//   options[METIS_OPTION_RTYPE]   = params->rtype;
//   options[METIS_OPTION_DBGLVL]  = params->dbglvl;
//   options[METIS_OPTION_UFACTOR] = params->ufactor;
//   options[METIS_OPTION_MINCONN] = params->minconn;
    
//   options[METIS_OPTION_SEED]    = params->seed;
//   options[METIS_OPTION_NITER]   = params->niter;
//   options[METIS_OPTION_NCUTS]   = params->ncuts;
  
     
     options[METIS_OPTION_CTYPE]    = METIS_CTYPE_SHEM; 
     //cout<<options[METIS_OPTION_CTYPE]; 
     options[METIS_OPTION_PTYPE]= METIS_PTYPE_KWAY;
     //cout<<options[METIS_OPTION_PTYPE]<<endl; //exit(0);
  
     //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
  
     options[METIS_OPTION_IPTYPE]   = METIS_IPTYPE_RANDOM;
     //cout<<options[METIS_OPTION_IPTYPE]<<endl; 
     options[METIS_OPTION_CONTIG]  = 1;//params->contig;
     //cout<< options[METIS_OPTION_CONTIG]<<endl; //exit(0);
     options[METIS_OPTION_MINCONN] = 1;
     options[METIS_OPTION_NITER]   = 10;
     options[METIS_OPTION_UFACTOR] = 100;
  
  eptr[0]=0;
  unsigned counter=0;
  for (unsigned iel=0; iel<nel; iel++) {
    unsigned ielt=el->GetElementType(iel);
    eptr[iel+1]=eptr[iel]+NVE[ielt][3];
//     cout << eptr[iel+1] << " ";
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,3); inode++){
      eind[counter]=el->GetElementVertexIndex(iel,inode)-1;
//       cout << eind[counter] << " ";
      counter ++;
   }
  }
  
  printf("NEL : %d \n",nel);
  printf("NVT : %d \n",nvt);
  cout << "NVT : " << nvt << endl;
  
  idx_t mnel = nel; 
  idx_t mnvt = nvt;
  idx_t ncommon = _dimension+1;
  nsubdom = _nprocs;
  
  if(nsubdom!=1) {
  //I call the mesh partioning function of Metis library (output is epart(own elem) and npart (own nodes))
  int err = METIS_PartMeshDual(&mnel, &mnvt, eptr, eind, NULL, NULL, &ncommon, &nsubdom, NULL, options, &objval, epart, npart);
  
  if(err==METIS_OK) {
    cout << " METIS PARTITIONING IS OK " << endl;
  }
  else if(err==METIS_ERROR_INPUT) {
    cout << " METIS_ERROR_INPUT " << endl;
    exit(1);
  }
  else if (err==METIS_ERROR_MEMORY) {
    cout << " METIS_ERROR_MEMORY " << endl;
    exit(2);
  }
  else {
    cout << " METIS_GENERIC_ERROR " << endl;
    exit(3);
   }
  }
  else {
    //serial computation
    for(unsigned i=0;i<nel;i++) {
      epart[i]=0;
    } 
  }
  
  delete [] eptr;
  delete [] eind;
  delete [] aux_vec;

  //dof map: piecewise liner 0, quadratic 1, biquadratic 2, piecewise constant 3, picewise discontinous linear 4 
  
  //resize the vector IS_Gmt2Mts_dof and dof
  for(int k=0;k<5;k++) {
    IS_Gmt2Mts_dof[k].resize(GetDofNumber(k)); 
    IS_Gmt2Mts_dof_offset[k].resize(nsubdom+1);
  }
  IS_Mts2Gmt_elem.resize(nel);
  IS_Mts2Gmt_elem_offset.resize(nsubdom+1);
  
  // I 
  for(unsigned i=0;i<nvt;i++) {
    npart[i]=nsubdom;
  }
  
  IS_Mts2Gmt_elem_offset[0]=0;
  vector <unsigned> IS_Gmt2Mts_dof_counter(5,0);
  
   for(int k=0;k<5;k++) {
     IS_Gmt2Mts_dof[k].assign(GetDofNumber(k),GetDofNumber(k)-1); 
     //TODO for domain decomposition pourposes! the non existing dofs point to the last dof!!!!!!
   }
  
  
  IS_Gmt2Mts_dof_counter[3]=0;
  IS_Gmt2Mts_dof_counter[4]=0;
  
  for(int isdom=0;isdom<nsubdom;isdom++) {
    for(unsigned iel=0;iel<nel;iel++){
      if(epart[iel]==isdom){
	//filling the piecewise IS_Mts2Gmt_elem metis->gambit
	IS_Mts2Gmt_elem[ IS_Gmt2Mts_dof_counter[3] ]=iel;
	IS_Gmt2Mts_dof[3][iel]=IS_Gmt2Mts_dof_counter[3];
	IS_Gmt2Mts_dof_counter[3]++;
	IS_Mts2Gmt_elem_offset[isdom+1]=IS_Gmt2Mts_dof_counter[3];
	// linear+quadratic+biquadratic
	for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
	  unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	  if(npart[ii]>isdom) {
	    npart[ii]=isdom;
	    IS_Gmt2Mts_dof[0][ii]=IS_Gmt2Mts_dof_counter[0];
	    IS_Gmt2Mts_dof_counter[0]++;
	    IS_Gmt2Mts_dof[1][ii]=IS_Gmt2Mts_dof_counter[1];
	    IS_Gmt2Mts_dof_counter[1]++;
	    IS_Gmt2Mts_dof[2][ii]=IS_Gmt2Mts_dof_counter[2];
	    IS_Gmt2Mts_dof_counter[2]++;
	  }
	}
	// quadratic+biquadratic
	for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++) {
	  unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	  if(npart[ii]>isdom){
	    npart[ii]=isdom;
	    IS_Gmt2Mts_dof[1][ii]=IS_Gmt2Mts_dof_counter[1];
	    IS_Gmt2Mts_dof_counter[1]++;
	    IS_Gmt2Mts_dof[2][ii]=IS_Gmt2Mts_dof_counter[2];
	    IS_Gmt2Mts_dof_counter[2]++;
	  }
	}
	// biquadratic
	for (unsigned inode=el->GetElementDofNumber(iel,1); inode<el->GetElementDofNumber(iel,3); inode++) {
	  unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	  if(npart[ii]>isdom){
	    npart[ii]=isdom;
	    IS_Gmt2Mts_dof[2][ii]=IS_Gmt2Mts_dof_counter[2];
	    IS_Gmt2Mts_dof_counter[2]++;
	  }
	}	
      }
    }
    for(unsigned k_dim=0;k_dim<_dimension+1;k_dim++){
      for(unsigned iel=0;iel<nel;iel++){
     	if(epart[iel]==isdom){
	  IS_Gmt2Mts_dof[4][iel+k_dim*nel]=IS_Gmt2Mts_dof_counter[4];
	  IS_Gmt2Mts_dof_counter[4]++;
	}
      }
    }
  }  
   
 cout << " 1 " << endl; 
   
//  // a lot of time!!!!!!!!!!!!!!!!!
//  
//  
//   counter=0;
//  
//   // reorder the vertices and midpoints vs the central nodes    
//   for(int isdom=0;isdom<nsubdom;isdom++){
//     for(unsigned i=IS_Mts2Gmt_elem_offset[isdom];i<IS_Mts2Gmt_elem_offset[isdom+1];i++){
//       unsigned iel=IS_Mts2Gmt_elem[i]; 
//       // cout << "i: " << i << "  dof_map: " << iel << endl;
//       for (unsigned inode=0; inode<el->GetElementDofNumber(iel,3); inode++) {
// // 	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
// // 	 if(npart[ii]==isdom){
// //  	   unsigned i0=IS_Gmt2Mts_dof[2][ii];
// // 	   if(i0>counter){
// // 	     IS_Gmt2Mts_dof[2][ii]=counter;
// // 	     IS_Gmt2Mts_dof[2][counter]=i0; 
// // 	     counter++;
// // 	   }
// // 	 }
// 	
// 	for(unsigned j=IS_Mts2Gmt_elem_offset[isdom];j<IS_Mts2Gmt_elem_offset[isdom+1];j++){
// 	  unsigned jel=IS_Mts2Gmt_elem[j];
// 	  for (unsigned jnode=el->GetElementDofNumber(jel,1); jnode<el->GetElementDofNumber(jel,3); jnode++) { 
// 	    unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
// 	    unsigned jj=el->GetElementVertexIndex(jel,jnode)-1;
// 	    if(npart[ii]==isdom && npart[jj]==isdom){
// 	      unsigned i0=IS_Gmt2Mts_dof[2][ii];
// 	      unsigned i1=IS_Gmt2Mts_dof[2][jj];
// 	      if(i0>i1){
// 		IS_Gmt2Mts_dof[2][ii]=i1;
// 		IS_Gmt2Mts_dof[2][jj]=i0; 
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
// 
//   cout << " 2 " << endl; 
//   
//   
//   // a lot of time!!!!!!!!!!!!!!!!!
//   
//   // reoreder the vertices vs the midpoints
//   for(int isdom=0;isdom<nsubdom;isdom++){
//     for(unsigned i=IS_Mts2Gmt_elem_offset[isdom];i<IS_Mts2Gmt_elem_offset[isdom+1];i++){
//       unsigned iel=IS_Mts2Gmt_elem[i]; 
//       for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
// 	for(unsigned j=IS_Mts2Gmt_elem_offset[isdom];j<IS_Mts2Gmt_elem_offset[isdom+1];j++){
// 	  unsigned jel=IS_Mts2Gmt_elem[j];
// 	  for (unsigned jnode=el->GetElementDofNumber(jel,0); jnode<el->GetElementDofNumber(jel,1); jnode++) { 
// 	    unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
// 	    unsigned jj=el->GetElementVertexIndex(jel,jnode)-1;
// 	    if(npart[ii]==isdom && npart[jj]==isdom){
// 	      { // for quadratic
// 		unsigned i0=IS_Gmt2Mts_dof[1][ii];
// 		unsigned i1=IS_Gmt2Mts_dof[1][jj];
// 		if(i0>i1){
// 		  IS_Gmt2Mts_dof[1][ii]=i1;
// 		  IS_Gmt2Mts_dof[1][jj]=i0; 
// 		}
// 	      }
// 	      { // for biquadratic
// 		unsigned i0=IS_Gmt2Mts_dof[2][ii];
// 		unsigned i1=IS_Gmt2Mts_dof[2][jj];
// 		if(i0>i1){
// 		  IS_Gmt2Mts_dof[2][ii]=i1;
// 		  IS_Gmt2Mts_dof[2][jj]=i0; 
// 		}
// 	      }  
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
//    
//   cout << " 3 " << endl;  
   
  // ghost vs own nodes
  vector <unsigned short> node_count(nvt,0);
  
  for(unsigned k=0;k<5;k++){
    ghost_size[k].assign(nsubdom,0);
    own_size[k].assign(nsubdom,0);
  }
  //unsigned ghost_counter[5]={0,0,0,0,0};
  
  for(int isdom=0;isdom<nsubdom;isdom++){
    
    own_size[3][isdom] = IS_Mts2Gmt_elem_offset[isdom+1]-IS_Mts2Gmt_elem_offset[isdom];
    own_size[4][isdom] = (IS_Mts2Gmt_elem_offset[isdom+1]-IS_Mts2Gmt_elem_offset[isdom])*(_dimension+1);
    
    for(unsigned i=IS_Mts2Gmt_elem_offset[isdom];i<IS_Mts2Gmt_elem_offset[isdom+1];i++){
      unsigned iel=IS_Mts2Gmt_elem[i]; 
      
       for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if( npart[ii] != isdom){
	    ghost_size[0][isdom]++;
	    ghost_size[1][isdom]++;
	    ghost_size[2][isdom]++;
	  }else {
	    own_size[0][isdom]++;
	    own_size[1][isdom]++;
	    own_size[2][isdom]++;
	  }
	}
      }
      
      
       for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if( npart[ii] != isdom){
	    ghost_size[1][isdom]++;
	    ghost_size[2][isdom]++;
	  }else {
	    own_size[1][isdom]++;
	    own_size[2][isdom]++;
	  }
	}
      }
      
      
      for (unsigned inode=el->GetElementDofNumber(iel,1); inode<el->GetElementDofNumber(iel,3); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if( npart[ii] != isdom){
	    ghost_size[2][isdom]++;
	  }else {
	    own_size[2][isdom]++;
	  }
	}
      }
      
      
    }
  }
  

  for(int k=0; k<5; k++) {
    ghost_nd[k].resize(nsubdom);
    ghost_nd_mts[k].resize(nsubdom);
    for(int isdom=0; isdom<nsubdom; isdom++) {
      ghost_nd[k][isdom].resize(ghost_size[k][isdom]);
      ghost_nd_mts[k][isdom].resize(ghost_size[k][isdom]);    
    }
  } 
  
  node_count.assign (nvt,0);  
  for(int k=0; k<5; k++) {
    ghost_size[k].assign(nsubdom,0);
  }
  
  for(int isdom=0;isdom<nsubdom;isdom++) {
    for(unsigned i=IS_Mts2Gmt_elem_offset[isdom];i<IS_Mts2Gmt_elem_offset[isdom+1];i++){
      unsigned iel=IS_Mts2Gmt_elem[i]; 
      
      
      for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if(npart[ii] != isdom){
 	    ghost_nd_mts[0][isdom][ghost_size[0][isdom]]=IS_Gmt2Mts_dof[0][ii];
            ghost_nd[0][isdom][ghost_size[0][isdom]]=ii;
	    ghost_size[0][isdom]++;
 	    ghost_nd_mts[1][isdom][ghost_size[1][isdom]]=IS_Gmt2Mts_dof[1][ii];
	    ghost_nd[1][isdom][ghost_size[1][isdom]]=ii;
	    ghost_size[1][isdom]++;
 	    ghost_nd_mts[2][isdom][ghost_size[2][isdom]]=IS_Gmt2Mts_dof[2][ii];
	    ghost_nd[2][isdom][ghost_size[2][isdom]]=ii;
	    ghost_size[2][isdom]++;
	  }
	}
      }
      
      for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if( npart[ii] != isdom){
 	    ghost_nd_mts[1][isdom][ghost_size[1][isdom]]=IS_Gmt2Mts_dof[1][ii];
	    ghost_nd[1][isdom][ghost_size[1][isdom]]=ii;
	    ghost_size[1][isdom]++;
	    ghost_nd_mts[2][isdom][ghost_size[2][isdom]]=IS_Gmt2Mts_dof[2][ii];
	    ghost_nd[2][isdom][ghost_size[2][isdom]]=ii;
	    ghost_size[2][isdom]++;
	  }
	}
      }
      
      for(unsigned inode=el->GetElementDofNumber(iel,1); inode<el->GetElementDofNumber(iel,3); inode++) {
	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	if(node_count[ii]<isdom+1){
	  node_count[ii]=isdom+1;
	  if( npart[ii] != isdom){
 	    ghost_nd_mts[2][isdom][ghost_size[2][isdom]]=IS_Gmt2Mts_dof[2][ii];
            ghost_nd[2][isdom][ghost_size[2][isdom]]=ii;
	    ghost_size[2][isdom]++;
	  }
	}
      } 
    }
  }
  
  
  MetisOffset.resize(5);
  for(int i=0;i<5;i++) 
    MetisOffset[i].resize(nsubdom+1);
  
  MetisOffset[0][0]=0;
  MetisOffset[1][0]=0;
  MetisOffset[2][0]=0;
  MetisOffset[3][0]=0;
  MetisOffset[4][0]=0;
  
  for(int i=1;i<=nsubdom;i++){
    MetisOffset[0][i]= MetisOffset[0][i-1]+own_size[0][i-1];
    MetisOffset[1][i]= MetisOffset[1][i-1]+own_size[1][i-1];
    MetisOffset[2][i]= MetisOffset[2][i-1]+own_size[2][i-1];
    MetisOffset[3][i]= IS_Mts2Gmt_elem_offset[i];
    MetisOffset[4][i]= IS_Mts2Gmt_elem_offset[i]*(_dimension+1);
    
  }
  
//   for(int i=0;i<5;i++){
//     for(int j=0;j<=nsubdom;j++){
//       cout<<MetisOffset[i][j]<<" ";
//     }
//     cout<<endl;
//   }
  
  /*
   for(unsigned i=0;i<3*nel;i++){
     IS_Gmt2Mts_dof[4][i]=i;
   }*/
  
//   for(int isdom=0;isdom<nsubdom;isdom++){
//     cout <<"nel: "<< IS_Mts2Gmt_elem_offset[isdom+1]-IS_Mts2Gmt_elem_offset[isdom] << endl;
//     cout<<"linear\n";
//     cout<<"ghost size: "<<ghost_size[0][isdom]<<endl;
//     cout<<"own size: "<<own_size[0][isdom]<<endl;
//     sort(&ghost_nd[0][isdom][0],&ghost_nd[0][isdom][ghost_size[0][isdom]]);
//     cout << "ghost nodes: ";
//     for(unsigned i=0;i<ghost_size[0][isdom];i++){
//       cout<<  ghost_nd_mts[0][isdom][i] <<" ";
//     }
//     cout<<endl;
//     for(unsigned i=IS_Mts2Gmt_elem_offset[isdom];i<IS_Mts2Gmt_elem_offset[isdom+1];i++){
//       unsigned iel=IS_Mts2Gmt_elem[i];
//       for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
// 	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
// 	cout<<IS_Gmt2Mts_dof[0][ii]<<" ";
//       }
//       cout<<endl;
//     }
//     cout<<"quadratic\n";
//     cout<<"ghost size: "<<ghost_size[1][isdom]<<endl;
//     cout<<"own size: "<<own_size[1][isdom]<<endl;
//     sort(&ghost_nd[1][isdom][0],&ghost_nd[1][isdom][ghost_size[1][isdom]]);
//     cout << "ghost nodes: ";
//     for(unsigned i=0;i<ghost_size[1][isdom];i++){
//       cout<<  ghost_nd_mts[1][isdom][i] <<" ";
//     }
//     cout<<endl;
//     for(unsigned i=IS_Mts2Gmt_elem_offset[isdom];i<IS_Mts2Gmt_elem_offset[isdom+1];i++){
//       unsigned iel=IS_Mts2Gmt_elem[i];
//       for (unsigned inode=0; inode<el->GetElementDofNumber(iel,1); inode++) {
// 	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
// 	cout<<IS_Gmt2Mts_dof[1][ii]<<" ";
//       }
//       cout<<endl;
//     }
//    
//     
//     cout<<"biquadratic\n";
// 
//     cout<<"ghost size: "<<ghost_size[2][isdom]<<endl;
//     cout<<"own size: "<<own_size[2][isdom]<<endl;
//     sort(&ghost_nd[2][isdom][0],&ghost_nd[2][isdom][ghost_size[2][isdom]]);
//     cout << "ghost nodes: ";
//     for(unsigned i=0;i<ghost_size[2][isdom];i++){
//       //cout<<isdom<<" "<<i<<endl;
//       cout<<  ghost_nd_mts[2][isdom][i] <<" ";
//     }
//     cout<<endl;
//     
//     for(unsigned i=IS_Mts2Gmt_elem_offset[isdom];i<IS_Mts2Gmt_elem_offset[isdom+1];i++){
//       unsigned iel=IS_Mts2Gmt_elem[i];
//       for (unsigned inode=0; inode<el->GetElementDofNumber(iel,3); inode++) {
// 	unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
// 	cout<<IS_Gmt2Mts_dof[2][ii]<<" ";
//       }
//       cout<<endl;
//     }
//     cout<<endl;
//   }
    
//   cout<<nel<<" "<<counter_el<<endl;
//   for(unsigned i=0;i<3*nel;i++){
//     cout<<IS_Gmt2Mts_dof[4][i]<<" ";
//   }
//   cout<<endl;
  
// /*  
//   for(unsigned i=0;i<nel;i++){
//     cout<<IS_Gmt2Mts_dof[3][IS_Mts2Gmt_elem[i]]<<" ";
//   }
//   cout<<endl;*/
//   
//    for(unsigned i=0;i<nel;i++){
//     cout<<IS_Mts2Gmt_elem[IS_Gmt2Mts_dof[3][i]]<<" ";
//   }
//   cout<<endl;
//   
//   for(int isdom=0;isdom<nsubdom;isdom++){
//     for(unsigned k_dim=0;k_dim<_dimension+1;k_dim++){
//       for(unsigned i=IS_Mts2Gmt_elem_offset[isdom];i<IS_Mts2Gmt_elem_offset[isdom+1];i++){
// 	cout<<IS_Gmt2Mts_dof[4][nel*k_dim+IS_Mts2Gmt_elem[i]]<<" ";
//       }
//     }
//   }
//   cout<<endl;
//   
//   
//   for(int isdom=0;isdom<nsubdom;isdom++){
//     for(unsigned k_dim=0;k_dim<_dimension+1;k_dim++){
//       for(unsigned i=IS_Mts2Gmt_elem_offset[isdom];i<IS_Mts2Gmt_elem_offset[isdom+1];i++){
// 	cout<<IS_Gmt2Mts_dof[4][nel*k_dim+IS_Mts2Gmt_elem[i]]<<" ";
//       }
//     }
//   }
//   cout<<endl;
//   
  
  //exit(0);
  return; 
  
}

/**
 *  This constructur generates a fine mesh level, $l_i$, from a coarse mesh level $l_{i-1}$, $i>0$
 **/
//------------------------------------------------------------------------------------------------------
mesh::mesh(const unsigned & igrid,elem *elc) {
  
  MPI_Comm_rank(MPI_COMM_WORLD, &_iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &_nprocs);
  
  const unsigned vertex_index[6][8][8]= {{{1,9,25,12,17,21,27,24},{9,2,10,25,21,18,22,27},
      {25,10,3,11,27,22,19,23},{12,25,11,4,24,27,23,20},
      {17,21,27,24,5,13,26,16},{21,18,22,27,13,6,14,26},
      {27,22,19,23,26,14,7,15},{24,27,23,20,16,26,15,8}
    },
    { {1,5,7,8},{2,6,5,9},{3,7,6,10},{8,9,10,4},
      {5,6,8,9},{6,10,8,9},{5,8,6,7},{6,8,10,7}
    },
    { {1,7,9,13,16,18},{2,8,7,14,17,16},
      {3,9,8,15,18,17},{7,8,9,16,17,18},
      {13,16,18,4,10,12},{14,17,16,5,11,10},
      {15,18,17,6,12,11},{16,17,18,10,11,12}
    },
    {{1,5,9,8},{5,2,6,9},{9,6,3,7},{8,9,7,4}},
    {{1,4,6},{2,5,4},{3,6,5},{4,5,6}},
    {{1,3},{3,2}}
  };

  const unsigned face_index[6][6][4][2]= {{{{0,0},{1,0},{4,0},{5,0}},
      {{1,1},{2,1},{5,1},{6,1}},
      {{2,2},{3,2},{6,2},{7,2}},
      {{3,3},{0,3},{7,3},{4,3}},
      {{0,4},{1,4},{2,4},{3,4}},
      {{4,5},{5,5},{6,5},{7,5}}
    },
    { {{0,0},{1,0},{2,0},{6,3}},
      {{0,1},{1,3},{3,1},{4,3}},
      {{1,1},{2,3},{3,2},{5,1}},
      {{2,1},{0,3},{3,3},{7,2}}
    },
    { {{0,0},{1,2},{4,0},{5,2}},
      {{1,0},{2,2},{5,0},{6,2}},
      {{2,0},{0,2},{6,0},{4,2}},
      {{0,3},{1,3},{2,3},{3,3}},
      {{4,4},{5,4},{6,4},{7,4}}
    },
    { {{0,0},{1,0}},
      {{1,1},{2,1}},
      {{2,2},{3,2}},
      {{3,3},{0,3}}
    },
    { {{0,0},{1,2}},
      {{1,0},{2,2}},
      {{2,0},{0,2}}
    },
    { {{0,0}},
      {{1,1}}
    }
  };

  const unsigned midpoint_index[6][12][2]= {{{0,1},{1,2},{2,3},{3,0},
      {4,5},{5,6},{6,7},{7,4},
      {0,4},{1,5},{2,6},{3,7}
    },
    { {0,1},{1,2},{2,0},
      {0,3},{1,3},{2,3}
    },
    { {0,1},{1,2},{2,0},
      {3,4},{4,5},{5,3},
      {0,3},{1,4},{2,5}
    },
    {{0,1},{1,2},{2,3},{3,0}},
    {{0,1},{1,2},{2,0}},
    {{0,1}}
  };

  const unsigned midpoint_index2[6][7][8]= {{{0,8,0,11,16},
      {0,0,9,0,0,17},
      {0,0,0,10,0,0,18},
      {0,0,0,0,0,0,0,19},
      {0,0,0,0,0,12,0,15},
      {0,0,0,0,0,0,13},
      {0,0,0,0,0,0,0,14}
    },
    { {0,4,6,7},
      {0,0,5,8},
      {0,0,0,9}
    },
    { {0,6,8,12},
      {0,0,7,0,13},
      {0,0,0,0,0,14},
      {0,0,0,0,9,11},
      {0,0,0,0,0,10}
    },
    { {0,4,0,7},
      {0,0,5},
      {0,0,0,6}
    },
    { {0,3,5},
      {0,0,4}
    },
    {{0,2}}
  };
  grid=igrid;

//   nel=elc->GetRefinedElementNumber()*REF_INDEX;
  nel=elc->GetRefinedElementNumber()*_ref_index;
  el=new elem(elc,mesh::_ref_index);



  unsigned jel=0;
  //divide each coarse element in 8(3D), 4(2D) or 2(1D) fine elemets and find all the vertices

  el->SetElementGroupNumber(elc->GetElementGroupNumber());
  el->SetNumberElementFather(elc->GetElementNumber());

  for (unsigned iel=0; iel<elc->GetElementNumber(); iel++) {
    if ( elc->GetRefinedElementIndex(iel) ) {
      elc->SetRefinedElementIndex(iel,jel+1u);
      unsigned elt=elc->GetElementType(iel);
      // project element type
//       for(unsigned j=0;j<REF_INDEX;j++){
      for (unsigned j=0; j<_ref_index; j++) {
        el->SetElementType(jel+j,elt);
        el->SetElementFather(jel+j,iel);
	elc->SetChildElement(iel,j,jel+j);
      }

      unsigned elg=elc->GetElementGroup(iel);
      unsigned elmat=elc->GetElementMaterial(iel);
      // project element group
//       for(unsigned j=0;j<REF_INDEX;j++){
      for (unsigned j=0; j<_ref_index; j++) {
        el->SetElementGroup(jel+j,elg);
	el->SetElementMaterial(jel+j,elmat);
      }

      // project vertex indeces
//       for(unsigned j=0;j<REF_INDEX;j++)
      for (unsigned j=0; j<_ref_index; j++)
        for (unsigned inode=0; inode<elc->GetElementDofNumber(iel,0); inode++)
          el->SetElementVertexIndex(jel+j,inode,elc->GetElementVertexIndex(iel,vertex_index[elt][j][inode]-1u));
      // project face indeces
      for (unsigned iface=0; iface<elc->GetElementFaceNumber(iel); iface++) {
        int value=elc->GetFaceElementIndex(iel,iface);
        if (0>value)
// 	  for(unsigned jface=0;jface<FACE_INDEX;jface++)
          for (unsigned jface=0; jface<_face_index; jface++)
            el->SetFaceElementIndex(jel+face_index[elt][iface][jface][0],face_index[elt][iface][jface][1], value);
      }
      // update element numbers
//       jel+=REF_INDEX;
      jel+=_ref_index;
//       el->AddToElementNumber(REF_INDEX,elt);
      el->AddToElementNumber(_ref_index,elt);
    }
  }

  nvt=elc->GetNodeNumber();
  el->SetVertexNodeNumber(nvt);

  
  //find all the elements near each vertex
  BuildAdjVtx();
  //initialize to zero all the middle edge points
  for (unsigned iel=0; iel<nel; iel++)
    for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++)
      el->SetElementVertexIndex(iel,inode,0);
  //find all the middle edge points
  for (unsigned iel=0; iel<nel; iel++) {
    unsigned ielt=el->GetElementType(iel);
    unsigned istart=el->GetElementDofNumber(iel,0);
    unsigned iend=el->GetElementDofNumber(iel,1);
    for (unsigned inode=istart; inode<iend; inode++)
      if (0==el->GetElementVertexIndex(iel,inode)) {
        nvt++;
        el->SetElementVertexIndex(iel,inode,nvt);
        unsigned im=el->GetElementVertexIndex(iel,midpoint_index[ielt][inode-istart][0]);
        unsigned ip=el->GetElementVertexIndex(iel,midpoint_index[ielt][inode-istart][1]);
        //find all the near elements which share the same middle edge point
        for (unsigned j=0; j<el->GetVertexElementNumber(im-1u); j++) {
          unsigned jel=el->GetVertexElementIndex(im-1u,j)-1u;
          if (jel>iel) {
            unsigned jm=0,jp=0;
            unsigned jelt=el->GetElementType(jel);
            for (unsigned jnode=0; jnode<el->GetElementDofNumber(jel,0); jnode++) {
              if (el->GetElementVertexIndex(jel,jnode)==im) {
                jm=jnode+1u;
                break;
              }
            }
            if (jm!=0) {
              for (unsigned jnode=0; jnode<el->GetElementDofNumber(jel,0); jnode++) {
                if (el->GetElementVertexIndex(jel,jnode)==ip) {
                  jp=jnode+1u;
                  break;
                }
              }
              if (jp!=0) {
                if (jp<jm) {
                  unsigned tp=jp;
                  jp=jm;
                  jm=tp;
                }
                el->SetElementVertexIndex(jel,midpoint_index2[jelt][--jm][--jp],nvt);
              }
            }
          }
        }
      }
  }
  el->SetMidpointNodeNumber(nvt - el->GetVertexNodeNumber());
  // Now build kmid
  Buildkmid();
  Buildkel();
  
  //for parallel computations
  cout << "  start partioning " << endl;
  if (_nprocs>=1) generate_metis_mesh_partition();
  cout << " end partitioning " << endl;
  // some output
  cout <<"grid\tnel\tnvt"<< endl;
  cout <<grid<<"\t"<<nel<<"\t"<<nvt<<endl;
  cout <<"nv0\tnv1\tnv2"<< endl;
  cout <<el->GetVertexNodeNumber()<<"\t"
       <<el->GetMidpointNodeNumber()<<"\t"
       <<el->GetCentralNodeNumber()<<endl
       <<"-------------------------"<<endl;
}

/**
 * This function read the data form the custom 1D file.
 * It is used in the constructor of the coarse mesh.
 **/
void mesh::Read1D(const char infile [], vector < vector < double> > &vt) {
  // set the file
  std::ifstream inf;
  std::string str2;

  grid=0;
  // read NUMBER OF ELEMENTS ******************** A
  inf.open(infile);
  if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " can not read parameters\n";
    exit(0);
  }
  str2="0";
  while (str2.compare("NEL") != 0) inf >> str2;
  inf >> nel;
  inf.close();

  nvt=2*nel+1;
  el= new elem(nel);
  el->AddToElementNumber(nel,"Line");

  el->SetElementGroupNumber(1);
  for (unsigned iel=0; iel<nel; iel++) {
    el->SetElementGroup(iel,0);
    el->SetElementType(iel,5);
    el->SetElementVertexIndex(iel,0,iel+1);
    el->SetElementVertexIndex(iel,1,iel+2);
    el->SetElementVertexIndex(iel,2,(nel+1)+(iel+1));
  }
  // end read NUMBER OF ELEMENTS **************** A

  // read NODAL COORDINATES **************** B
  inf.open(infile);
  if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " cannot read nodes\n";
    exit(0);
  }
  double x0,x1,h;
  while (str2.compare("COORDINATES") != 0) inf >> str2;
  inf >> x0>>x1;
  inf.close();

  if (x1<x0) {
    double temp=x1;
    x1=x0;
    x0=temp;
  }
  vt[0].resize(nvt);
  vt[1].resize(nvt);
  vt[2].resize(nvt);

    
  h=(x1-x0)/nel;
  for (unsigned iel=0; iel<nel; iel++) {
    vt[0][iel]=x0+iel*h;
    vt[1][iel]=0.;
    vt[2][iel]=0.;
    vt[0][nel+(iel+1)]=vt[0][iel]+h/2.;
    vt[1][nel+(iel+1)]=0.;
    vt[2][nel+(iel+1)]=0.;
  }
  vt[0][nel]=x1;
  vt[1][nel]=0.;
  vt[2][nel]=0.;
  // end read NODAL COORDINATES ************* B

  // set boundary **************** C
  el->SetFaceElementIndex(0,0,-2);
  el->SetFaceElementIndex(nel-1,1,-3);
  // end set boundary ****************
  el->SetNodeNumber(nvt);
  el->SetVertexNodeNumber(nel+1);
  el->SetMidpointNodeNumber(nel);
  el->SetCentralNodeNumber(0);
  
 
  
}



/**
 * This function read the data form the gambit file.
 * It is used in the constructor of the coarse mesh.
 **/
void mesh::ReadGambit(const char infile [], vector < vector < double> > &vt, const double Lref) {
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

  std::ifstream inf;
  std::string str2;
  unsigned ngroup;
  unsigned nbcd;
  unsigned dim;
  double x,y,z;
//   double _lref = 1.;

  grid=0;
  // read control data ******************** A
  inf.open(infile);
  if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " can not read parameters\n";
    exit(0);
  }
  str2="0";
  while (str2.compare("NDFVL") != 0) inf >> str2;
  inf >> nvt >> nel >>  ngroup >> nbcd >> dim >> str2 ;
  mesh::_dimension = dim;
  mesh::_ref_index = pow(2,mesh::_dimension);  // 8*DIM[2]+4*DIM[1]+2*DIM[0];
  mesh::_face_index = pow(2,mesh::_dimension-1u); // 4*DIM[2]+2*DIM[1]+1*DIM[0];
//   cout << " ref index " << _ref_index << "    " << REF_INDEX << endl;
//   cout << " face index " << _face_index << "    " << FACE_INDEX << endl;
//   cout << " mesh dim " << mesh::_dimension <<  endl;
  inf >> str2;
  if (str2.compare("ENDOFSECTION") != 0) {
    cout<<"error control data mesh"<<endl;
    exit(0);
  }
//   cout << "***************" << _dimension << endl;
  inf.close();
  // end read control data **************** A
  // read ELEMENT/cell ******************** B
  inf.open(infile);
  if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " cannot read elements\n";
    exit(0);
  }
  el= new elem(nel);
  while (str2.compare("ELEMENTS/CELLS") != 0) inf >> str2;
  inf >> str2;
  for (unsigned iel=0; iel<nel; iel++) {
    el->SetElementGroup(iel,1);
    unsigned nve;
    inf >> str2 >> str2 >> nve;
    if (nve==27) {
      el->AddToElementNumber(1,"Hex");
      el->SetElementType(iel,0);
    } else if (nve==10) {
      el->AddToElementNumber(1,"Tet");
      el->SetElementType(iel,1);
    } else if (nve==18) {
      el->AddToElementNumber(1,"Wedge");
      el->SetElementType(iel,2);
    } else if (nve==9) {
      el->AddToElementNumber(1,"Quad");
      el->SetElementType(iel,3);
    }
//     else if(nve==6 && DIM[1]==1){
    else if (nve==6 && mesh::_dimension==2) {
      el->AddToElementNumber(1,"Triangle");
      el->SetElementType(iel,4);
    }
//     else if(nve==3 && DIM[0]==1){
    else if (nve==3 && mesh::_dimension==1) {
      el->AddToElementNumber(1,"Line");
      el->SetElementType(iel,5);
    } else {
      cout<<"Error! Invalid element type in reading Gambit File!"<<endl;
      cout<<"Error! Use a second order discretization"<<endl;
      exit(0);
    }
    for (unsigned i=0; i<nve; i++) {
      unsigned inode=GambitVertexIndex[el->GetElementType(iel)][i];
      double value;
      inf>>value;
      el->SetElementVertexIndex(iel,inode,value);
    }
  }
  inf >> str2;
  if (str2.compare("ENDOFSECTION") != 0) {
    cout<<"error element data mesh"<<endl;
    exit(0);
  }
  inf.close();

  // end read  ELEMENT/CELL **************** B

  // read NODAL COORDINATES **************** C
  inf.open(infile);
  if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " cannot read nodes\n";
    exit(0);
  }
  while (str2.compare("COORDINATES") != 0) inf >> str2;
  inf >> str2;  // 2.0.4
  vt[0].resize(nvt);
  vt[1].resize(nvt);
  vt[2].resize(nvt);

  if (mesh::_dimension==3) {
    for (unsigned j=0; j<nvt; j++) {
      inf >> str2 >> x >> y >> z;
      vt[0][j] = x/Lref;
      vt[1][j] = y/Lref;
      vt[2][j] = z/Lref;
    }
  }

  else if (mesh::_dimension==2) {
    for (unsigned j=0; j<nvt; j++) {
      inf >> str2 >> x >> y;
      vt[0][j] = x/Lref;
      vt[1][j] = y/Lref;
      vt[2][j] = 0.;
    }
  }

  else if (mesh::_dimension==1) {
    for (unsigned j=0; j<nvt; j++) {
      inf >> str2 >> x;
      vt[0][j] = x/Lref;
      vt[1][j]=0.;
      vt[2][j]=0.;
    }
  }
  inf >> str2; // "ENDOFSECTION"
  if (str2.compare("ENDOFSECTION") != 0) {
    cout<<"error node data mesh 1"<<endl;
    exit(0);
  }
  inf.close();
  // end read NODAL COORDINATES ************* C

  // read GROUP **************** E
  inf.open(infile);
  if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " cannot read group\n";
    exit(0);
  }
  el->SetElementGroupNumber(ngroup);
  for (unsigned k=0; k<ngroup; k++) {
    int ngel;
    int name;
    int mat;
    while (str2.compare("GROUP:") != 0) inf >> str2;
    inf >> str2 >> str2 >> ngel >> str2 >> mat >>str2 >> str2 >>name>> str2;
//     cout<<ngel<<endl;
//     cout<<name<<endl;
//     cout << mat << endl;
    for (int i=0; i<ngel; i++) {
      int iel;
      inf >> iel;
      el->SetElementGroup(iel-1,name);
      el->SetElementMaterial(iel-1,mat);
    }
    inf >> str2;
    if (str2.compare("ENDOFSECTION") != 0) {
      cout<<"error group data mesh"<<endl;
      exit(0);
    }
  }
  inf.close();
  // end read boundary **************** E

  // read boundary **************** D
  inf.open(infile);
  if (!inf) {
    cout<<"Generic-mesh file "<<infile<< " cannot read boudary\n";
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
      iface=GambitFaceIndex[el->GetElementType(iel)][iface-1u];
      el->SetFaceElementIndex(iel,iface,value);
    }
    inf >> str2;
    if (str2.compare("ENDOFSECTION") != 0) {
      cout<<"error boundary data mesh"<<endl;
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
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,1); inode++) {
      for (unsigned jel=0; jel<nel; jel++) {
	for (unsigned jnode=el->GetElementDofNumber(jel,1); jnode<el->GetElementDofNumber(jel,3); jnode++) { 
	  unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	  unsigned jj=el->GetElementVertexIndex(jel,jnode)-1;
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
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
      for (unsigned jel=0; jel<nel; jel++) {
        for (unsigned jnode=el->GetElementDofNumber(jel,0); jnode<el->GetElementDofNumber(jel,1); jnode++) {
          unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
	  unsigned jj=el->GetElementVertexIndex(jel,jnode)-1;
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
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,3); inode++) {
      unsigned ii=el->GetElementVertexIndex(iel,inode)-1;
      el->SetElementVertexIndex(iel,inode,dof_index[ii]);
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
 
  el->SetNodeNumber(nvt);

  unsigned nv0=0;
  for (unsigned iel=0; iel<nel; iel++)
    for (unsigned inode=0; inode<el->GetElementDofNumber(iel,0); inode++) {
      unsigned i0=el->GetElementVertexIndex(iel,inode);
      if (nv0<i0) nv0=i0;
  }
  el->SetVertexNodeNumber(nv0);

  unsigned nv1=0;
  for (unsigned iel=0; iel<nel; iel++)
    for (unsigned inode=el->GetElementDofNumber(iel,0); inode<el->GetElementDofNumber(iel,1); inode++) {
      unsigned i1=el->GetElementVertexIndex(iel,inode);
      if (nv1<i1) nv1=i1;
  }
  el->SetMidpointNodeNumber(nv1-nv0);

  el->SetCentralNodeNumber(nvt-nv1);
  
};




/**
 * This function searches all the elements around all the vertices
 **/
void mesh::BuildAdjVtx() {
  el->AllocateVertexElementMemory();
  for (unsigned iel=0; iel<nel; iel++) {
    for (unsigned inode=0; inode < el->GetElementDofNumber(iel,0); inode++) {
      unsigned ii=el->GetElementVertexIndex(iel,inode)-1u;
      unsigned jj=0;
      while ( 0 != el->GetVertexElementIndex(ii,jj) ) jj++;
      el->SetVertexElementIndex(ii,jj,iel+1u);
    }
  }
}


/**
 * This function generates kmid for hex and wedge elements
 **/
void mesh::Buildkmid() {
  for (unsigned iel=0; iel<el->GetElementNumber(); iel++)
    for (unsigned inode=el->GetElementDofNumber(iel,1); inode<el->GetElementDofNumber(iel,2); inode++)
      el->SetElementVertexIndex(iel,inode,0);

  for (unsigned iel=0; iel<el->GetElementNumber(); iel++) {
    for (unsigned iface=0; iface<el->GetElementFaceNumber(iel,0); iface++) {
      unsigned inode=el->GetElementDofNumber(iel,1)+iface;
      if ( 0==el->GetElementVertexIndex(iel,inode) ) {
        el->SetElementVertexIndex(iel,inode,++nvt);
        unsigned i1=el->GetFaceVertexIndex(iel,iface,0);
        unsigned i2=el->GetFaceVertexIndex(iel,iface,1);
        unsigned i3=el->GetFaceVertexIndex(iel,iface,2);
        for (unsigned j=0; j< el->GetVertexElementNumber(i1-1u); j++) {
          unsigned jel= el->GetVertexElementIndex(i1-1u,j)-1u;
          if (jel>iel) {
            for (unsigned jface=0; jface<el->GetElementFaceNumber(jel,0); jface++) {
              unsigned jnode=el->GetElementDofNumber(jel,1)+jface;
              if ( 0==el->GetElementVertexIndex(jel,jnode) ) {
                unsigned j1=el->GetFaceVertexIndex(jel,jface,0);
                unsigned j2=el->GetFaceVertexIndex(jel,jface,1);
                unsigned j3=el->GetFaceVertexIndex(jel,jface,2);
                unsigned j4=el->GetFaceVertexIndex(jel,jface,3);
                if ((i1==j1 || i1==j2 || i1==j3 ||  i1==j4 )&&
                    (i2==j1 || i2==j2 || i2==j3 ||  i2==j4 )&&
                    (i3==j1 || i3==j2 || i3==j3 ||  i3==j4 )) {
                  el->SetElementVertexIndex(jel,jnode,nvt);
                }
              }
            }
          }
        }
      }
    }
  }

  for (unsigned iel=0; iel<el->GetElementNumber(); iel++) {
    if (0==el->GetElementType(iel)) {
      el->SetElementVertexIndex(iel,26,++nvt);
    }
    if (3==el->GetElementType(iel)) {
      el->SetElementVertexIndex(iel,8,++nvt);
    }
  }
  el->SetNodeNumber(nvt);

  unsigned nv0= el->GetVertexNodeNumber();
  unsigned nv1= el->GetMidpointNodeNumber();
  el->SetCentralNodeNumber(nvt-nv0-nv1);

}

/**
 * This function stores the element adiacent to the element face (iel,iface)
 * and stores it in kel[iel][iface]
 **/
void mesh::Buildkel() {
  for (unsigned iel=0; iel<el->GetElementNumber(); iel++) {
    for (unsigned iface=0; iface<el->GetElementFaceNumber(iel); iface++) {
      if ( el->GetFaceElementIndex(iel,iface) <= 0) {
        unsigned i1=el->GetFaceVertexIndex(iel,iface,0);
        unsigned i2=el->GetFaceVertexIndex(iel,iface,1);
        unsigned i3=el->GetFaceVertexIndex(iel,iface,2);
        for (unsigned j=0; j< el->GetVertexElementNumber(i1-1u); j++) {
          unsigned jel= el->GetVertexElementIndex(i1-1u,j)-1u;
          if (jel>iel) {
            for (unsigned jface=0; jface<el->GetElementFaceNumber(jel); jface++) {
              if ( el->GetFaceElementIndex(jel,jface) <= 0) {
                unsigned j1=el->GetFaceVertexIndex(jel,jface,0);
                unsigned j2=el->GetFaceVertexIndex(jel,jface,1);
                unsigned j3=el->GetFaceVertexIndex(jel,jface,2);
                unsigned j4=el->GetFaceVertexIndex(jel,jface,3);
// 		if((DIM[2]==1 &&
                if ((mesh::_dimension==3 &&
                     (i1==j1 || i1==j2 || i1==j3 ||  i1==j4 )&&
                     (i2==j1 || i2==j2 || i2==j3 ||  i2==j4 )&&
                     (i3==j1 || i3==j2 || i3==j3 ||  i3==j4 ))||
// 		   (DIM[1]==1 &&
                    (mesh::_dimension==2 &&
                     (i1==j1 || i1==j2 )&&
                     (i2==j1 || i2==j2 ))||
// 		   (DIM[0]==1 &&
                    (mesh::_dimension==1 &&
                     (i1==j1))
                   ) {
                  el->SetFaceElementIndex(iel,iface,jel+1u);
                  el->SetFaceElementIndex(jel,jface,iel+1u);
                }
              }
            }
          }
        }
      }
    }
  }
}


unsigned mesh::GetDimension() {
  return mesh::_dimension;
}

unsigned mesh::GetRefIndex() {
  return mesh::_ref_index;
}


/**
 * This function returns the number of mesh nodes for different type of elemets
 **/
unsigned mesh::GetDofNumber(const unsigned type) const {
// cout << "Dim:" << DIMENSION << "   " << mesh::_dimension-1 << endl;
// if(mesh::_dimension != 2) exit(1);
  switch (type) {
  case 0:
    return el->GetVertexNodeNumber();
    break;
  case 1:
    return el->GetVertexNodeNumber()+el->GetMidpointNodeNumber();
    break;
  case 2:
    return nvt;
    break;
  case 3:
    return nel;
    break;
  case 4:
//     return nel*(2+DIMENSION);
    return nel*(2+mesh::_dimension-1);
    break;
  }
  return 0;
}

unsigned mesh::GetElementNumber() const {
  return nel;
}
unsigned mesh::GetGridNumber() const {
  return grid;
}


/**
 * This function copies the refined element index vector in other_vector
 **/
void mesh::copy_elr(vector <unsigned> &other_vec) const {
  for (unsigned i=0; i<nel; i++)
    other_vec[i]=el->GetRefinedElementIndex(i);
}


void mesh::AllocateAndMarkStructureNode() {
  el->AllocateNodeRegion();
  for (unsigned iel=0; iel<nel; iel++) {
//     bool structure = MarkStructureNode( el->GetElementGroup(iel) );
    int flag_mat = el->GetElementMaterial(iel);

//     cout << "structure: " << structure << "   flag_mat: " << flag_mat << endl;  
    
//     if (structure==1) {
    if (flag_mat==4) {
      unsigned nve=el->GetElementDofNumber(iel);
      for ( unsigned i=0; i<nve; i++) {
        unsigned inode=el->GetElementVertexIndex(iel,i)-1u;
        el->SetNodeRegion(inode, 1);
      }
    }
  }
  return;
}

unsigned mesh::GetEndIndex(const unsigned i) const{
  return _END_IND[i];
};
