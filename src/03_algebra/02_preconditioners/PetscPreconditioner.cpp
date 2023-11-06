/*=========================================================================

 Program: FEMUS
 Module: PetscPreconditioner
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "FemusConfig.hpp"

#ifdef HAVE_PETSC

// Local Includes
#include "PetscPreconditioner.hpp"
#include "PetscMacro.hpp"
#include "PetscMatrix.hpp"
#include "PetscVector.hpp"

#include <mpi.h>
#include "petscpc.h" 
#include <petscksp.h>

namespace femus {



// =============================================================================
  void PetscPreconditioner::apply(const NumericVector & x, NumericVector & y) {
    PetscVector & x_pvec = static_cast<PetscVector&>(const_cast<NumericVector&>(x));
    PetscVector & y_pvec = static_cast<PetscVector&>(const_cast<NumericVector&>(y));
    Vec x_vec = x_pvec.vec();
    Vec y_vec = y_pvec.vec();
    PCApply(_pc, x_vec, y_vec);
  }
// =================================================
  void PetscPreconditioner::init() {
    if(!this->_matrix) {
      std::cerr << "ERROR: No matrix set for PetscPreconditioner, but init() called" << std::endl;
      abort();
    }
    //Clear the preconditioner in case it has been created in the past
    if(!this->_is_initialized)  {
      //Create the preconditioning object
      PCCreate(MPI_COMM_WORLD, &_pc);
      //Set the PCType
      set_petsc_preconditioner_type(this->_preconditioner_type, _pc);
// #ifdef LIBMESH_HAVE_PETSC_HYPRE
//     if(this->_preconditioner_type == AMG_PRECOND)
//       PCHYPRESetType(this->_pc, "boomerang");
// #endif
      PetscMatrix * pmatrix = libmeshM_cast_ptr<PetscMatrix*, SparseMatrix >(this->_matrix);
      _mat = pmatrix->mat();
    }
    //PCSetOperators(_pc,_mat,_mat,SAME_NONZERO_PATTERN);
    PCSetOperators(_pc, _mat, _mat); //PETSC3p5
    this->_is_initialized = true;
  }

// =====================================================
  void PetscPreconditioner::set_petsc_preconditioner_type
  (const PreconditionerType & preconditioner_type, PC & pc,  const int &parallelOverlapping) {
    int ierr = 0;
    switch(preconditioner_type)  {

      case IDENTITY_PRECOND:
        ierr = PCSetType(pc, (char*) PCNONE);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

      case CHOLESKY_PRECOND:
        ierr = PCSetType(pc, (char*) PCCHOLESKY);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

      case ICC_PRECOND:
        ierr = PCSetType(pc, (char*) PCICC);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;


      case ILU_PRECOND:
      {
        int nprocs;
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs); //TODO
        // In serial, just set the ILU preconditioner type
        if(nprocs == 1)
        {
          ierr = PCSetType(pc, (char*) PCILU);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
        }
        else
        {
//        But PETSc has no truly parallel ILU, instead you have to set
//        an actual parallel preconditioner (e.g. block Jacobi (parlleloverlapping 0) or ASM (parlleloverlapping >0))
//	  and then assign ILU sub-preconditioners.

          set_petsc_preconditioner_type(ASM_PRECOND, pc);
          PCASMSetOverlap(pc, parallelOverlapping);
          PCSetUp(pc);

          // Set ILU as the sub preconditioner type
          set_petsc_subpreconditioner_type(PCILU, pc);
        }
        break;
      }
      case LU_PRECOND: {
        int nprocs;
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        if(nprocs == 1) {
          ierr = PCSetType(pc, (char*) PCLU);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
        }
        else {
          ierr = PCSetType(pc, (char*) PCLU);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
	  
	  ierr = PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          ierr = PCFactorSetUpMatSolverType(pc);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          Mat       F;
          ierr = PCFactorGetMatrix(pc, &F);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
          ierr = MatMumpsSetIcntl(F, 14, 30);
          CHKERRABORT(MPI_COMM_WORLD, ierr);
        }
        break;
      }

      case SLU_PRECOND:
        ierr = PCSetType(pc, (char*) PCLU);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        ierr = PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;    //here we set the SuperLU_dist solver package

      case MLU_PRECOND: //here we set the MUMPS parallel direct solver package
        ierr = PCSetType(pc, (char*) PCLU);
        CHKERRABORT(MPI_COMM_WORLD, ierr);

        ierr = PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        ierr = PCFactorSetUpMatSolverType(pc);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        Mat       F;
        ierr = PCFactorGetMatrix(pc, &F);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        ierr = MatMumpsSetIcntl(F, 14, 30);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

      case ULU_PRECOND: //here we set the Umfpack serial direct solver package
        ierr = PCSetType(pc, (char*) PCLU);
        CHKERRABORT(MPI_COMM_WORLD, ierr);

        ierr = PCFactorSetMatSolverType(pc, MATSOLVERUMFPACK);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        ierr = PCFactorSetUpMatSolverType(pc);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

      case MCC_PRECOND:
        ierr = PCSetType(pc, (char*) PCCHOLESKY);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        ierr = PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
        CHKERRABORT(MPI_COMM_WORLD, ierr);                   //here we set the MUMPS parallel direct solver package
        break;

      case ASM_PRECOND:
      case ASM_MULTIPLICATIVE_PRECOND:
        PCSetType(pc, (char*) PCASM);
        PCASMSetType (pc,  PC_ASM_BASIC);
        PCASMSetLocalType (pc, PC_COMPOSITE_MULTIPLICATIVE);
        break;
      case ASM_ADDITIVE_PRECOND:
        PCSetType(pc, (char*) PCASM);
        PCASMSetType (pc,  PC_ASM_BASIC);
        PCASMSetLocalType (pc, PC_COMPOSITE_ADDITIVE);
        break;  

      case FIELDSPLIT_PRECOND:
      case FIELDSPLIT_ADDITIVE_PRECOND:
        PCSetType(pc, (char*) PCFIELDSPLIT);
        PCFieldSplitSetType (pc, PC_COMPOSITE_ADDITIVE);
        break;
      case FIELDSPLIT_MULTIPLICATIVE_PRECOND:
        PCSetType(pc, (char*) PCFIELDSPLIT);
        PCFieldSplitSetType (pc, PC_COMPOSITE_MULTIPLICATIVE);
        break;
      case FIELDSPLIT_SYMMETRIC_MULTIPLICATIVE_PRECOND:
        PCSetType(pc, (char*) PCFIELDSPLIT);
        PCFieldSplitSetType (pc, PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE);
        break;
      case FIELDSPLIT_SCHUR_PRECOND:
        PCSetType(pc, (char*) PCFIELDSPLIT);
        PCFieldSplitSetType (pc, PC_COMPOSITE_SCHUR);
        break;  
        
      case JACOBI_PRECOND:
        ierr = PCSetType(pc, (char*) PCJACOBI);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

      case BLOCK_JACOBI_PRECOND:
        ierr = PCSetType(pc, (char*) PCBJACOBI);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

      case SOR_PRECOND:
        ierr = PCSetType(pc, (char*) PCSOR);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

      case EISENSTAT_PRECOND:
        ierr = PCSetType(pc, (char*) PCEISENSTAT);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

      case AMG_PRECOND:
        ierr = PCSetType(pc, (char*) PCHYPRE);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

      case MG_PRECOND:
        ierr = PCSetType(pc, (char*) PCMG);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

      case LSC_PRECOND:
        ierr = PCSetType(pc, (char*) PCLSC);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

#if !(PETSC_VERSION_LESS_THAN(2,1,2))
        // Only available for PETSC >= 2.1.2
      case USER_PRECOND:
        ierr = PCSetType(pc, (char*) PCMAT);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;
#endif

      case SHELL_PRECOND:
        ierr = PCSetType(pc, (char*) PCSHELL);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        break;

      default:
        std::cerr << "ERROR:  Unsupported PETSC Preconditioner: "
                  << preconditioner_type << std::endl
                  << "Continuing with PETSC defaults" << std::endl;
    }

    //Let the commandline override stuff
    if(preconditioner_type != AMG_PRECOND && preconditioner_type != MG_PRECOND)   PCSetFromOptions(pc);   //!!!!!!
  }



  void PetscPreconditioner::set_petsc_subpreconditioner_type(const PCType type, PC& pc)  {

    int ierr;
    KSP* subksps;
    int nlocal;

    ierr = PCASMGetSubKSP(pc, &nlocal, PETSC_NULL, &subksps);
    CHKERRABORT(MPI_COMM_WORLD, ierr);

    PetscReal epsilon = 1.e-16;

    for(int i = 0; i < nlocal; i++) {
      PC subpc;

      ierr = KSPGetPC(subksps[i], &subpc);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      ierr = KSPSetTolerances(subksps[i], PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, 1);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      ierr = KSPSetFromOptions(subksps[i]);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      ierr = PCSetType(subpc, type);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      ierr = PCFactorSetZeroPivot(subpc, epsilon);
      CHKERRABORT(MPI_COMM_WORLD, ierr);

      ierr = PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
      CHKERRABORT(MPI_COMM_WORLD, ierr);
    }
  }













//   // For catching PETSc error return codes
//   int ierr = 0;
//
//   // All docs say must call KSPSetUp or PCSetUp before calling PCBJacobiGetSubKSP.
//   // You must call PCSetUp after the preconditioner operators have been set, otherwise you get the:
//   //
//   // "Object is in wrong state!"
//   // "Matrix must be set first."
//   //
//   // error messages...
//   ierr = PCSetUp(pc);
//   CHKERRABORT(MPI_COMM_WORLD, ierr);
//
//   // To store array of local KSP contexts on this processor
//   KSP* subksps;
// //
// //   // the number of blocks on this processor
//   int n_local;
// //
// //   // The global number of the first block on this processor.
// //   // This is not used, so we just pass PETSC_NULL instead.
// //   // int first_local;
// //
// //   // Fill array of local KSP contexts
//   ierr = PCBJacobiGetSubKSP(pc, &n_local, PETSC_NULL, &subksps);
//   CHKERRABORT(MPI_COMM_WORLD, ierr);
// //
// //   // Loop over sub-ksp objects, set ILU preconditioner
//   for(int i = 0; i < n_local; ++i) {
//     // Get pointer to sub KSP object's PC
//     PC subpc;
//     ierr = KSPGetPC(subksps[i], &subpc);
//     CHKERRABORT(MPI_COMM_WORLD, ierr);
//
//     // Set requested type on the sub PC
//     ierr = PCSetType(subpc, type);
//     CHKERRABORT(MPI_COMM_WORLD, ierr);
//   }
//}


//------------------------------------------------------------------


} //end namespace femus



#endif // #ifdef LIBMESH_HAVE_PETSC
