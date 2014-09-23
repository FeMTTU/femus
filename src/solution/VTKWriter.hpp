/*=========================================================================

 Program: FEMUS
 Module: VTKWriter
 Authors: Eugenio Aulisa, Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __vtkwriter_h_
#define __vtkwriter_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Writer.hpp"


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class MultiLevelProblem;


class VTKWriter : public Writer {

public:

    /** Constructor. */
    VTKWriter(MultiLevelSolution& ml_sol);

    /** Destructor */
    virtual ~VTKWriter();

    /** write output function */
    virtual void write_system_solutions(const char order[], std::vector<std::string>& vars, const unsigned time_step=0);

};


} //end namespace femus



#endif
