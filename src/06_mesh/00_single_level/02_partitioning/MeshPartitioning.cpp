/*=========================================================================

 Program: FEMUS
 Module: MeshPartitioning
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
#include "MeshPartitioning.hpp"
#include "Mesh.hpp"

namespace femus {

  
MeshPartitioning::MeshPartitioning(const Mesh &mesh) : _mesh(mesh) {

}

}
