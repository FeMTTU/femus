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

#ifndef __femus_solution_GMVWriter_one_level_hpp__
#define __femus_solution_GMVWriter_one_level_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Writer_one_level.hpp"


namespace femus {



  //------------------------------------------------------------------------------
  // Forward declarations
  //------------------------------------------------------------------------------


  class GMVWriter_one_level : public Writer_one_level {

  public:

// === Constructors / Destructor  - BEGIN =================
      /** Constructor. */
      GMVWriter_one_level(const MultiLevelSolution * ml_sol);

      /** Constructor. */
      GMVWriter_one_level(const MultiLevelMesh * ml_mesh);

      /** Destructor */
      virtual ~GMVWriter_one_level();
// === Constructors / Destructor  - END =================

// === Write - BEGIN =================
  public:

      /** write output function */
      void Write(const std::string output_path, 
                 const std::string order,
                 const std::vector < std::string > & vars = std::vector < std::string > (), 
                 const unsigned time_step = _time_step_index_default);

  private:
        
      /** write output function */
      void Write(const unsigned level_in,
                 const std::string output_path, 
                 const std::string order,
                 const std::vector < std::string > & vars = std::vector < std::string > (), 
                 const unsigned time_step = _time_step_index_default);
// === Write - END =================
    

};


} //end namespace femus



#endif
