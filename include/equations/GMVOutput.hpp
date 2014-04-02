/*=========================================================================

 Program: FEMUS
 Module: GMVOutput
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __gmvoutput_h_
#define __gmvoutput_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Output.hpp"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class MultiLevelProblem;


class GMVOutput : public Output {

public:

  /** Constructor. */
  GMVOutput(MultiLevelProblem& ml_probl);

  /** Destructor */
  virtual ~GMVOutput();
  
  /** write output function */
  virtual void write_system_solutions(const char order[], std::vector<std::string>& vars, const unsigned time_step=0);
  
};

#endif