/*=========================================================================

 Program: FEMUS
 Module: XDMFOutput
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __xdmfoutput_h_
#define __xdmfoutput_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Output.hpp"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class MultiLevelProblem;


class XDMFOutput : public Output {

public:

  /** Constructor. */
  XDMFOutput(MultiLevelProblem& ml_probl);

  /** Destructor */
  virtual ~XDMFOutput();
  
  /** write output function */
  virtual void write_system_solutions(const char order[], std::vector<std::string>& vars, const unsigned time_step = 0);
  
  /** write a wrapper file for paraview to open all the files of an history toghether */
  void write_solution_wrapper(const char type[]) const;
  
};

#endif