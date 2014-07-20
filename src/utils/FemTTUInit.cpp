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
 
#ifdef HAVE_MPI 
//   	  MPI_Init (&argc, &argv);
//          int i;
//           MPI_Comm_size(comm_world_in, &i);
//           assert (i >= 0);
//           _size = static_cast<unsigned int>(i);

//           MPI_Comm_rank(comm_world_in, &i);
//           assert (i >= 0);
//           _rank = static_cast<unsigned int>(i);
#endif	  

#ifdef HAVE_PETSC  
  int ierr = PetscInitialize (&argc, &argv, NULL, NULL); CHKERRABORT(PETSC_COMM_WORLD,ierr);
      
std::cout << " FemTTUInit(): PETSC_COMM_WORLD initialized" << std::endl << std::endl;      

#endif 

  // redirect libMesh::out to nothing on all
  // other processors unless explicitly told
  // not to via the --keep-cout command-line argument.
  int i;  MPI_Comm_rank(MPI_COMM_WORLD, &i);

 // std::cout << " FemTTUInit: PETSC_COMM_WORLD initialized from proc " << i << "\n";

  if ( i != 0 ){std::cout.rdbuf(NULL);}
  
  return;
}


// =======================================================
/// This function initializes the libraries if it is parallel
FemTTUInit::~FemTTUInit() { // ========================

// 1 proc= nothing to do
#ifdef HAVE_PETSC 
  PetscFinalize();
  std::cout << std::endl << " ~FemTTUInit(): PETSC_COMM_WORLD ends" << std::endl;    
#endif
 
#ifdef HAVE_MPI
//       MPI_Comm_free (MPI_COMM_WORLD);
//       MPI_Finalize();
//   std::cout << " ~FemTTUInit(): MPI_COMM_WORLD ends \n";  
#endif
  
 return;
} 


} //end namespace femus


