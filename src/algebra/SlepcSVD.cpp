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

double SlepcSVD::get_2norm_condition_number() {

return 1.;

}

}

#endif