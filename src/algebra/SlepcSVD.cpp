/*=========================================================================

 Program: FEMUS
 Module: SlepcSVD
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "SlepcSVD.hpp"
#include "iostream"
#include "PetscMatrix.hpp"

#if defined(HAVE_PETSC) && defined(HAVE_SLEPC)

namespace femus {
  
  
SlepcSVD::SlepcSVD() {
 
  SlepcInitialize(NULL,NULL,(char*)0, NULL);
  
  SVDCreate(PETSC_COMM_WORLD,&m_svd);

  std::cout << "debug: slepc has been initialized" << std::endl;
}


SlepcSVD::~SlepcSVD() {
 
  SVDDestroy(&m_svd);
  
  SlepcFinalize();

  std::cout << "debug: slepc has been finalized" << std::endl;
}


void SlepcSVD::set_operator(Mat A) {
  SVDSetOperator(m_svd, A);
}

void SlepcSVD::set_operator(SparseMatrix* Amat) {
  Mat A = (static_cast< PetscMatrix* >(Amat))->mat();
  SVDSetOperator(m_svd, A);
}


void SlepcSVD::init() {
  SVDSetFromOptions(m_svd);
  SVDSetDimensions(m_svd, 1, PETSC_DEFAULT, PETSC_DEFAULT);
}


double SlepcSVD::compute_2norm_condition_number() {
  
  PetscInt       nconv1,nconv2;
  PetscReal      sigma_1,sigma_n;

  // Solve the singular value problem
  // First request a singular value from one end of the spectrum
  SVDSetWhichSingularTriplets(m_svd, SVD_LARGEST);
  SVDSolve(m_svd);
  
  // Get number of converged singular values
  SVDGetConverged(m_svd, &nconv1);

  // Get converged singular values: largest singular value is stored in sigma_1.
  // In this example, we are not interested in the singular vectors
  if (nconv1 > 0) {
    SVDGetSingularTriplet(m_svd,0,&sigma_1,NULL,NULL);
  } else {
    PetscPrintf(PETSC_COMM_WORLD," Unable to compute large singular value!\n\n");
  }

  // Request a singular value from the other end of the spectrum
  SVDSetWhichSingularTriplets(m_svd, SVD_SMALLEST);
  SVDSolve(m_svd);

  // Get number of converged eigenpairs
  SVDGetConverged(m_svd, &nconv2);

  // Get converged singular values: smallest singular value is stored in sigma_n.
  // As before, we are not interested in the singular vectors
  if (nconv2 > 0) {
    SVDGetSingularTriplet(m_svd, 0, &sigma_n, NULL, NULL);
  } else {
    PetscPrintf(PETSC_COMM_WORLD," Unable to compute small singular value!\n\n");
  }

  // return the estimated condition number
  if (nconv1 > 0 && nconv2 > 0) {
    double cond_numb = (double)(sigma_1/sigma_n);
    PetscPrintf(PETSC_COMM_WORLD, " Computed singular values: sigma_1=%6f, sigma_n=%6f\n", (double)sigma_1, (double)sigma_n);
    PetscPrintf(PETSC_COMM_WORLD, " Estimated condition number: sigma_1/sigma_n=%6f\n\n", cond_numb);
    return cond_numb;
  }
  else {
    return -1; 
  }
  
}

}

#endif