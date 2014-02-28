#include "FemusExtLib_conf.hpp"

#ifdef FEMUS_HAVE_PETSC

// C++ includes

// Local Includes
// #include "auto_ptr.h"
#include "petsc_preconditionerM.hpp"
#include "petsc_macroM.hpp"
#include "petsc_matrixM.hpp"
#include "petsc_vectorM.hpp"

#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "mpi.h"  //for MPI_COMM_WORLD
#pragma GCC diagnostic warning "-Wunused-parameter"


// =============================================================================
void PetscPreconditionerM::apply(const NumericVectorM & x, NumericVectorM & y){
  PetscVectorM & x_pvec = static_cast<PetscVectorM&>(const_cast<NumericVectorM&>(x));
  PetscVectorM & y_pvec = static_cast<PetscVectorM&>(const_cast<NumericVectorM&>(y));
  Vec x_vec = x_pvec.vec();  Vec y_vec = y_pvec.vec();
  PCApply(_pc,x_vec,y_vec);
}
// =================================================
void PetscPreconditionerM::init (){
  if(!this->_matrix) {
    std::cerr << "ERROR: No matrix set for PetscPreconditioner, but init() called" << std::endl;
    abort();
  }
  //Clear the preconditioner in case it has been created in the past
  if(!this->_is_initialized)  {
    //Create the preconditioning object
    PCCreate(MPI_COMM_WORLD,&_pc);
    //Set the PCType
    set_petsc_preconditioner_type(this->_preconditioner_type, _pc);
// #ifdef LIBMESH_HAVE_PETSC_HYPRE
//     if(this->_preconditioner_type == AMG_PRECOND)
//       PCHYPRESetType(this->_pc, "boomerang");
// #endif
    PetscMatrixM * pmatrix = libmeshM_cast_ptr<PetscMatrixM*, SparseMatrixM >(this->_matrix);
    _mat = pmatrix->mat();
  }
  PCSetOperators(_pc,_mat,_mat,SAME_NONZERO_PATTERN);
  this->_is_initialized = true;
}

// =====================================================
void PetscPreconditionerM::set_petsc_preconditioner_type
                (const PreconditionerTypeM & preconditioner_type, PC & pc){
  int ierr = 0;
  switch (preconditioner_type)  {
  case IDENTITY_PRECONDM:
    ierr = PCSetType (pc, (char*) PCNONE);      CHKERRABORT(MPI_COMM_WORLD,ierr); break;
	
  case CHOLESKY_PRECONDM:
    ierr = PCSetType (pc, (char*) PCCHOLESKY);  CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case ICC_PRECONDM:
    ierr = PCSetType (pc, (char*) PCICC);       CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case ILU_PRECONDM:
    ierr = PCSetType (pc, (char*) PCILU);       CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case LU_PRECONDM:
    ierr = PCSetType (pc, (char*) PCLU);        CHKERRABORT(MPI_COMM_WORLD,ierr); break;
      
  case ASM_PRECONDM:
    ierr = PCSetType (pc, (char*) PCASM);       CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case JACOBI_PRECONDM:
    ierr = PCSetType (pc, (char*) PCJACOBI);    CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case BLOCK_JACOBI_PRECONDM:
    ierr = PCSetType (pc, (char*) PCBJACOBI);   CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case SOR_PRECONDM:
    ierr = PCSetType (pc, (char*) PCSOR);       CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case EISENSTAT_PRECONDM:
    ierr = PCSetType (pc, (char*) PCEISENSTAT); CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  case AMG_PRECONDM:
    ierr = PCSetType (pc, (char*) PCHYPRE);     CHKERRABORT(MPI_COMM_WORLD,ierr); break;

#if !(PETSC_VERSION_LESS_THAN(2,1,2))
    // Only available for PETSC >= 2.1.2      
  case USER_PRECONDM:
    ierr = PCSetType (pc, (char*) PCMAT);       CHKERRABORT(MPI_COMM_WORLD,ierr); break;
#endif

  case SHELL_PRECONDM:
    ierr = PCSetType (pc, (char*) PCSHELL);     CHKERRABORT(MPI_COMM_WORLD,ierr); break;

  default:
    std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
              << preconditioner_type       << std::endl
              << "Continuing with PETSC defaults" << std::endl;
  }

  //Let the commandline override stuff
  if( preconditioner_type != AMG_PRECONDM )   PCSetFromOptions(pc);
}

//------------------------------------------------------------------

#endif // #ifdef LIBMESH_HAVE_PETSC
