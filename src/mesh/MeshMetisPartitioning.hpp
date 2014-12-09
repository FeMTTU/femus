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

#ifndef __meshmetispartitioning_hpp__
#define __meshmetispartitioning_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
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
    
    /** Refinement functions */
    void DoPartition();
    
private:



};


}

#endif