/*=========================================================================

 Program: FEMuS
 Module: MeshPartitioning
 Authors: Simone Bnà, Eugenio Aulisa

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_mesh_MeshPartitioning_hpp__
#define __femus_mesh_MeshPartitioning_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "ParallelObject.hpp"

namespace femus {


class Mesh;


/**
 * This is the \p MeshPartitioning class.  This class is the abstract interface to
 * mesh partitioning algorithms.
*/

class MeshPartitioning : public ParallelObject {

public:

    /** Constructor */
    MeshPartitioning(Mesh& mesh);

    /** destructor */
    ~MeshPartitioning() {}; 
    
    /** Refinement functions */
    
protected:

   Mesh& _mesh;                 //< reference to the mesh which is built by refinement

};


}

#endif