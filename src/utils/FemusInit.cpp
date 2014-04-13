#include <iostream>

#include "FemusInit.hpp"

#include "FEMTTUConfig.h"


#ifdef HAVE_MPI
# include "mpi.h"
#endif


#ifdef HAVE_PETSC
# include "PetscMacro.hpp"
EXTERN_C_FOR_PETSC_BEGIN
# include "petsc.h"
# include "petscerror.h"
EXTERN_C_FOR_PETSC_END
#endif


#include "paral.hpp"

FemusInit::FemusInit(int & argc, char** & argv/*, MPI_Comm comm_world_in*/) {
 
#ifdef HAVE_MPI
  std::cout <<  "  MPI_COMM_WORLD first " << MPI_COMM_WORLD << std::endl;
//clearly, at this point MPI_Init was not called yet, so the different processors 
// have different MPI_COMM_WORLD values.
//but, after, they all must have the SAME COMMUNICATOR in common


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
std::cout << "PETSC_COMM_WORLD first " << PETSC_COMM_WORLD << std::endl;

//      PETSC_COMM_WORLD = comm_world_in;
  int ierr = PetscInitialize (&argc, &argv, NULL, NULL);  //TODO VALGRIND
              CHKERRABORT(PETSC_COMM_WORLD,ierr);
      
 std::cout <<   "PETSC_COMM_WORLD after" << PETSC_COMM_WORLD << std::endl;
   std::cout << "  MPI_COMM_WORLD after" <<  MPI_COMM_WORLD << std::endl;
#endif
    // redirect libMesh::out to nothing on all
  // other processors unless explicitly told
  // not to via the --keep-cout command-line argument.


//after initializing MPI, we can fill proc_rank and proc_size in namespace Paral
// paral::set_rank();
// paral::set_size();


//  std::cout << "My processor is " << paral::proc_rank <<std::endl;

//  std::cout << "The number of processors is " << paral::proc_size <<std::endl;
  
  //****** STD::COUT
  //1)either copy std::cout to a filestream
  //so that you have BOTH cout and YOUR FILE
  //2)or just REDIRECT std::cout to YOUR FILE
  //2) would be easy just with 2> in the run command  
//BUT the fact is that i want to do that INSIDE MY OUTTIME_FOLDER,
  //but i dont know what its name will be until i compose it!
  //so after composing i can redirect the std::cout
  //before that it will go 
  //clearly,this is a case where the interaction between SHELL and PROGRAM
//   cannot be avoided and is important.
//And this can be considered only in the context of the program,
//which means that you have to deal with c++ streams and you cant try
// to just OUTSOURCE everything to SYSTEM CALLS in BASH
//so i have to find a solution for the redirecting within c++
//and not by using shell symbols ("2>")
  
  
}




FemusInit::~FemusInit() {


#ifdef HAVE_PETSC
      PetscFinalize();
#endif

#ifdef HAVE_MPI
//       MPI_Comm_free (MPI_COMM_WORLD);
//       MPI_Finalize();
#endif


}
  
  
  