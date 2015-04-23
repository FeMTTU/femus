/*=========================================================================

 Program: FEMUS
 Module: MeshMetisPartitioning
 Authors: Simone Bnà, Eugenio Aulisa
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MeshMetisPartitioning.hpp"
#include "Mesh.hpp"
#include "FemusConfig.hpp"

#ifdef HAVE_METIS
  #include "metis.h"
#endif

//C++ include
#include <iostream>


namespace femus {

using std::cout;
using std::endl;  
  
  
MeshMetisPartitioning::MeshMetisPartitioning(Mesh& mesh) : MeshPartitioning(mesh) {  
  
}


//------------------------------------------------------------------------------------------------------
void MeshMetisPartitioning::DoPartition() {
  
#ifndef HAVE_METIS
  
  std::cerr << "Fatal error: Metis library was not found. Metis partioning algorithm cannot be called!" << std::endl;
  exit(1);

#else  
   
  unsigned eind_size = _mesh.el->GetElementNumber("Hex")*NVE[0][2]      + _mesh.el->GetElementNumber("Tet")*NVE[1][2] 
                     + _mesh.el->GetElementNumber("Wedge")*NVE[2][2]    + _mesh.el->GetElementNumber("Quad")*NVE[3][2] 
                     + _mesh.el->GetElementNumber("Triangle")*NVE[4][2] + _mesh.el->GetElementNumber("Line")*NVE[5][2];

  int nelem = _mesh.GetNumberOfElements();
  
  vector < idx_t > eptr(nelem+1);
  vector < idx_t > eind(eind_size);
  
  idx_t objval;
  idx_t options[METIS_NOPTIONS]; 
  
  METIS_SetDefaultOptions(options);
    
  options[METIS_OPTION_NUMBERING]= 0;
  options[METIS_OPTION_DBGLVL]   = 0;
  options[METIS_OPTION_CTYPE]    = METIS_CTYPE_SHEM; 
  options[METIS_OPTION_PTYPE]	 = METIS_PTYPE_KWAY;
  options[METIS_OPTION_IPTYPE]   = METIS_IPTYPE_RANDOM;
  options[METIS_OPTION_CONTIG]   = 0;
  options[METIS_OPTION_MINCONN]  = 1;
  options[METIS_OPTION_NITER]    = 10;
  options[METIS_OPTION_UFACTOR]  = 100;
  
  eptr[0]=0;
  unsigned counter=0;
  for (unsigned iel=0; iel<nelem; iel++) {
    unsigned ielt=_mesh.el->GetElementType(iel);
    eptr[iel+1]=eptr[iel]+NVE[ielt][2];
    
    for (unsigned inode=0; inode<_mesh.el->GetElementDofNumber(iel,2); inode++){
      eind[counter]=_mesh.el->GetElementVertexIndex(iel,inode)-1;
      counter++;
    }
    
  }
  
  int nnodes = _mesh.GetNumberOfNodes();
  int ncommon = ( _mesh.GetDimension() == 1 ) ? 1 : _mesh.GetDimension()+1;
  _mesh.nsubdom = _nprocs;
    
  _mesh.epart.resize(nelem);
  _mesh.npart.resize(nnodes);
  
  if(_mesh.nsubdom!=1) {
    
  //I call the Mesh partioning function of Metis library (output is epart(own elem) and npart (own nodes))
  int err = METIS_PartMeshDual(&nelem, &nnodes, &eptr[0], &eind[0], NULL, NULL, &ncommon, &_mesh.nsubdom, NULL, options, &objval, &_mesh.epart[0], &_mesh.npart[0]);
  
  if(err == METIS_OK) {
    std::cout << " METIS PARTITIONING IS OK " << std::endl;
  }
  else if(err == METIS_ERROR_INPUT) {
    cout << " METIS_ERROR_INPUT " << endl;
    exit(1);
  }
  else if (err == METIS_ERROR_MEMORY) {
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
    for(unsigned i=0;i<nelem;i++) {
      _mesh.epart[i]=0;
    } 
  }
    
  return;
  
#endif
  
}

void MeshMetisPartitioning::DoPartition(const Mesh& meshc){
  _mesh.nsubdom = _nprocs;
  _mesh.epart.resize( _mesh.GetNumberOfElements() );
  _mesh.npart.resize( _mesh.GetNumberOfNodes() );
  unsigned refIndex = _mesh.GetRefIndex();
  for (unsigned iel=0; iel < meshc.GetNumberOfElements(); iel++) {
    for (unsigned j=0; j < refIndex; j++) {
      _mesh.epart[ iel * refIndex + j ] = meshc.epart[iel];
    }
  }
}

}