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
#include "MultiLevelMesh.hpp"


//C++ include
#include <iostream>


namespace femus {

  
MeshPartitioning::MeshPartitioning(Mesh &mesh) : _mesh(mesh) {

}

}