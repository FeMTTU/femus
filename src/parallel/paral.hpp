#ifndef __femus_parallel_paral_hpp__
#define __femus_parallel_paral_hpp__


#include "FemusConfig.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#endif


namespace femus {



///This namespace will hold global data and functions for parallel purposes
///In this way you dont have to pass these values through a particular instantiation of a class,
///but everyone can access these


namespace paral {

//a namespace is a list of data and functions, NOT OF CALLS. A list of calls is A FUNCTION
//it is not a class (type), in fact it cannot have different instantiations
//it may hold also data. Since namespaces do not provide private data, you have to enclose
//them in a namespace and call them like that


//MULTIPLE DEFINITIONS -> error at LINK TIME
// RE-DEFINITION        -> error at COMPILE TIME
// multiple definitions of the same object are not allowed, but of the same typedef or class are allowed

//  namespace  {  //if you put an ANONYMOUS (not general) namespace it doesnt bother for multiple definitions
//but an anonymous namespace has only FILE SCOPE so its used only inside a single file
// int proc_rank=0;
// int proc_size=1;
//  }
//const variables can be redefined multiple times as global objects!
//without const it doesnt work
//with extern it works

//if the function is inline it doesnt bother, otherwise multiple definitions occur!
//inline doesnt give problems because the code is substituted everywhere you put the call,
//so there is no linking between objects
//so the linker does not have to link anything
//so if you want to have global functions, with INLINE you can include them in more files
//what can you do to have global DATA? extern, or sthg

inline int get_rank() {
    int proc_rank=0;
#ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
#endif
    return proc_rank;
}

inline int get_size() {
    int proc_size=1;
#ifdef HAVE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
#endif
    return  proc_size;
}



}




} //end namespace femus



#endif
