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


  class GMVWriter : public Writer {

  public:

// === Constructors / Destructor  - BEGIN =================
      /** Constructor. */
      GMVWriter(const MultiLevelSolution * ml_sol);

      /** Constructor. */
      GMVWriter(const MultiLevelMesh * ml_mesh);

// === Constructors / Destructor  - END =================

// === Write, at finest level - BEGIN =================
  public:

      /** write output function */
      void Write(const std::string output_path, 
                 const std::string order,
                 const std::vector < std::string > & vars = std::vector < std::string > (), 
                 const unsigned time_step = Writer_one_level::_time_step_index_default);

      /** write output function */
      void Write(const std::string filename_prefix, 
                 const std::string output_path, 
                 const std::string order,
                 const std::vector < std::string > & vars = std::vector < std::string > (), 
                 const unsigned time_step = Writer_one_level::_time_step_index_default);

private:
  
      /** at finest level:    */
     void Write(const std::string filename_prefix,
                       const std::string output_path, 
                       const std::string suffix_pre_extension, 
                       const std::string order,
                       const std::vector < std::string >& vars = std::vector < std::string > (), 
                       const unsigned time_step = Writer_one_level::_time_step_index_default)
     { abort(); }
// === Write, at finest level - END =================


// === Write, at arbitrary level - BEGIN =================
  private:        
        
      /** write output function */
      void Write(const unsigned level_in,
                 const std::string filename_prefix, 
                 const std::string output_path, 
                 const std::string suffix_pre_extension, 
                 const std::string order,
                 const std::vector < std::string > & vars = std::vector < std::string > (), 
                 const unsigned time_step = Writer_one_level::_time_step_index_default);
// === Write, at arbitrary level - END =================
    

};


} //end namespace femus



#endif
