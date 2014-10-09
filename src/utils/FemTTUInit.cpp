/*=========================================================================

 Program: FEMUS
 Module: FemTTUInit
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <iostream>
#include "FemTTUInit.hpp"


namespace femus {



// =======================================================
/// This function initializes the libraries if it is parallel
FemTTUInit::FemTTUInit(
    int & argc,            // integer program input
    char** & argv,         // char program input
    MPI_Comm comm_world_in // communicator for MPI direct
) {// ======================================================

#ifdef HAVE_PETSC
    int ierr = PetscInitialize (&argc, &argv, NULL, NULL);
    CHKERRABORT(PETSC_COMM_WORLD,ierr);

    std::cout << " FemTTUInit(): PETSC_COMM_WORLD initialized" << std::endl << std::endl;

#endif

#ifdef HAVE_MPI
    // redirect libMesh::out to nothing on all
    // other processors unless explicitly told
    // not to via the --keep-cout command-line argument.
    int i;
    MPI_Comm_rank(MPI_COMM_WORLD, &i);

    if ( i != 0 ) {
        std::cout.rdbuf(NULL);
    }
#endif

    return;
}


FemTTUInit::~FemTTUInit() { 

#ifdef HAVE_PETSC
    PetscFinalize();
    std::cout << std::endl << " ~FemTTUInit(): PETSC_COMM_WORLD ends" << std::endl;
#endif

    return;
}


} //end namespace femus


