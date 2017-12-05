/*=========================================================================

 Program: FEMUS
 Module: MeshASMPartitioning
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
#include "MeshASMPartitioning.hpp"
#include "Mesh.hpp"

//C++ include



namespace femus {

MeshASMPartitioning::MeshASMPartitioning(Mesh& mesh) : MeshPartitioning(mesh) {


}

//----------------------------------------------------------------------------------------------------------------
void MeshASMPartitioning::DoPartitionOld( const unsigned *block_size, vector < vector< unsigned > > &block_elements,
					 vector <unsigned> &block_type_range){

  unsigned iproc=processor_id();
  unsigned ElemOffset    = _mesh._elementOffset[iproc];
  unsigned ElemOffsetp1  = _mesh._elementOffset[iproc+1];
  unsigned OwnedElements = ElemOffsetp1 - ElemOffset;

  unsigned counter[2]={0,0};
  for (unsigned iel = ElemOffset; iel < ElemOffsetp1; iel++) {
    unsigned flag_mat   = _mesh.GetElementMaterial(iel);

    if(flag_mat < 4){
      counter[1]++;
    }
  }
  counter[0]=OwnedElements - counter[1];

  block_type_range.resize(2);

  unsigned flag_block[2]={4,2};

  unsigned block_start = 0;
  unsigned iblock = 0;
  while(iblock < 2){
    if(counter[iblock] !=0 ){ //material of this type is there
      unsigned reminder = counter[iblock] % block_size[iblock];
      unsigned blocks = (0 == reminder)? counter[iblock]/block_size[iblock] : counter[iblock]/block_size[iblock] + 1 ;
      block_elements.resize(block_start+blocks);

      for(int i = 0; i < blocks; i++)
	block_elements[block_start + i].resize(block_size[iblock]);
      if  (0 != reminder ){
	block_elements[block_start + (blocks-1u)].resize(reminder);
      }

      unsigned counter=0;
      for (unsigned iel = ElemOffset; iel < ElemOffsetp1; iel++) {
	unsigned flag_mat   = _mesh.GetElementMaterial(iel);
	if( flag_block[iblock] == flag_mat || flag_block[iblock] + 1 == flag_mat){
	  block_elements[ block_start + (counter / block_size[iblock]) ][ counter % block_size[iblock] ]=iel;
	  counter++;
	}
      }
      block_type_range[iblock]=block_start+blocks;
      block_start += blocks;
    }
    else{
      block_type_range[iblock]=block_start;
    }
    iblock++;
  }

}

void MeshASMPartitioning::DoPartition( const unsigned *block_size, vector < vector< unsigned > > &block_elements,
                                         vector <unsigned> &block_type_range){

    unsigned iproc = processor_id();
    unsigned ElemOffset    = _mesh._elementOffset[iproc];
    unsigned ElemOffsetp1  = _mesh._elementOffset[iproc + 1];
    unsigned OwnedElements = ElemOffsetp1 - ElemOffset;

    unsigned counter[3] = {0, 0, 0};
    unsigned flag_block[3] = {4, 3, 2};
    
    for (unsigned iel = ElemOffset; iel < ElemOffsetp1; iel++) {
      unsigned flag_mat   = _mesh.GetElementMaterial(iel);

      if (flag_mat == flag_block[0]) {
        counter[0]++;
      }
      else if (flag_mat == flag_block[1]) {
        counter[1]++;
      }
    }
    counter[2] = OwnedElements - counter[0] - counter[1];

    block_type_range.resize(3);

  

    unsigned block_start = 0;
    unsigned iMaterial = 0;
    while (iMaterial < 3) {
      if (counter[iMaterial] != 0 ) { //material of this type is there
        unsigned reminder = counter[iMaterial] % block_size[iMaterial];
        unsigned blocks = (0 == reminder) ? counter[iMaterial] / block_size[iMaterial] : counter[iMaterial] / block_size[iMaterial] + 1 ;
        block_elements.resize(block_start + blocks);

        for (int i = 0; i < blocks; i++) {
          block_elements[block_start + i].resize(block_size[iMaterial]);
        }
        if  (0 != reminder ) {
          block_elements[block_start + (blocks - 1u)].resize(reminder);
        }

        unsigned counter = 0;
        for (unsigned iel = ElemOffset; iel < ElemOffsetp1; iel++) {
          unsigned flag_mat   = _mesh.GetElementMaterial(iel);
          if ( flag_block[iMaterial] == flag_mat ) {
            block_elements[ block_start + (counter / block_size[iMaterial]) ][ counter % block_size[iMaterial] ] = iel;
            counter++;
          }
        }
        block_type_range[iMaterial] = block_start + blocks;
        block_start += blocks;
      } 
      else {
        block_type_range[iMaterial] = block_start;
      }
      iMaterial++;
    }
  }


}