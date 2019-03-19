/*=========================================================================

 Program: FEMuS
 Module: MeshMetisPartitioning
 Authors: Simone Bnà, Eugenio Aulisa

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_mesh_MeshMetisPartitioning_hpp__
#define __femus_mesh_MeshMetisPartitioning_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <vector>
#include "MeshPartitioning.hpp"

namespace femus {


class Mesh;


/**
 * This is the \p MeshMetisPartitioning class.  This class call the metis algorithm for
 * mesh partitioning.
*/

class MeshMetisPartitioning : public MeshPartitioning {

public:

    /** Constructor */
    MeshMetisPartitioning(Mesh& mesh);

    /** destructor */
    ~MeshMetisPartitioning() {};

    /** New Metis parallel partitioning:
     *  for coarse and AMR mesh */
    void DoPartition( std::vector < unsigned > &partition, const bool &AMR );

    /** Parallel partitioning imported from coarser mesh partition:
     *  for uniformed refined meshes */
    void DoPartition( std::vector < unsigned > &partition, const Mesh &meshc );

private:


};


}

#endif
