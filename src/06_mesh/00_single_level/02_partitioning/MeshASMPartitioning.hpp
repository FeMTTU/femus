/*=========================================================================

 Program: FEMuS
 Module: MeshASMPartitioning
 Authors: Simone Bnà, Eugenio Aulisa

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_mesh_MeshASMPartitioning_hpp__
#define __femus_mesh_MeshASMPartitioning_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MeshPartitioning.hpp"

#include <vector>


namespace femus {


class Mesh;

/**
 * This is the \p MeshASMPartitioning class.  This class calls a built-in mesh partitioner algorithm 
 * for ASM domain decomposition algorithms.
*/

class MeshASMPartitioning : public MeshPartitioning {

public:

    /** Constructor */
    MeshASMPartitioning(const Mesh& mesh);

    /** destructor */
    ~MeshASMPartitioning() {} 
    
    /** To be added */
    void DoPartition(const unsigned *block_size, std::vector < std::vector< unsigned > > &block_elements, std::vector <unsigned> &block_type_range);
    
    void DoPartitionOld(const unsigned *block_size, std::vector < std::vector< unsigned > > &block_elements, std::vector <unsigned> &block_type_range);
    
    
private:


};


}

#endif
