/*=========================================================================

 Program: FEMUS
 Module: VTKOutput
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __vtkoutput_h_
#define __vtkoutput_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Output.hpp"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class MultiLevelProblem;


class VTKOutput : public Output {

public:

  /** Constructor. */
  VTKOutput(MultiLevelProblem& ml_probl);

  /** Destructor */
  virtual ~VTKOutput();
  
  /** write output function */
  virtual void write_system_solutions();
  
};

#endif