/*=========================================================================

 Program: FEMUS
 Module: ParallelObject
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_parallel_ParallelObject_hpp__
#define __femus_parallel_ParallelObject_hpp__

#include "mpi.h"


namespace femus {



/**
   * This class forms the base class for all other classes
   * that are expected to be implemented in paralel.
*/


class ParallelObject
{

public:

    /** Constructor. Requires a reference to the communicator
     * that defines the object's parallel decomposition. */
    ParallelObject () {
        MPI_Comm_rank(MPI_COMM_WORLD, &_iproc);
        MPI_Comm_size(MPI_COMM_WORLD, &_nprocs);
    }

    /** Destructor. Virtual because we are a base class. */
    virtual ~ParallelObject () {}

    /** @returns the number of processors in the group. */
    int n_processors() const {
        return _nprocs;
    }

    /** @returns the rank of this processor in the group. */
    int processor_id() const {
        return _iproc;
    }


protected:

    int _nprocs;
    int _iproc;

};


} //end namespace femus



#endif
