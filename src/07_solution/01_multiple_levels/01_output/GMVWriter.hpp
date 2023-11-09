/*=========================================================================

 Program: FEMUS
 Module: GMVWriter
 Authors: Eugenio Aulisa, Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_solution_GMVWriter_hpp__
#define __femus_solution_GMVWriter_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Writer.hpp"


namespace femus {



  //------------------------------------------------------------------------------
  // Forward declarations
  //------------------------------------------------------------------------------
  class MultiLevelProblem;


  class GMVWriter : public Writer {

  public:

// === Constructors / Destructor  - BEGIN =================
      /** Constructor. */
      GMVWriter(MultiLevelSolution * ml_sol);

      /** Constructor. */
      GMVWriter(MultiLevelMesh * ml_mesh);

      /** Destructor */
      virtual ~GMVWriter();
// === Constructors / Destructor  - END =================

      /** write output function */
      void Write(const std::string output_path, const char order[], const std::vector < std::string > & vars = std::vector < std::string > (), const unsigned time_step=0) ;
    

    protected:
      
      void Mysol();
};


} //end namespace femus



#endif
