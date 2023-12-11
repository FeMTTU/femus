/*=========================================================================

 Program: FEMUS
 Module: Writer
 Authors: Eugenio Aulisa, Simone Bnà, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_solution_Writer_hpp__
#define __femus_solution_Writer_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Writer_one_level.hpp"
#include "ParallelObject.hpp"
#include "WriterEnum.hpp"

#include "fe_projection_matrices_Lagrange_continuous.hpp"

#include <vector>
#include <string>
#include <memory>
#include <iostream>


namespace femus {

  //------------------------------------------------------------------------------
  // Forward declarations
  //------------------------------------------------------------------------------
  class MultiLevelMesh;
  class MultiLevelSolution;
  class Solution;
  

  class Writer : public Writer_one_level {


// === Constructors / Destructor  - BEGIN =================
  public:
    
    /** Constructor. */
    Writer(const MultiLevelSolution * ml_sol);

    /** Constructor. */
    Writer(const MultiLevelMesh * ml_mesh);


    /** runtime selection of writer for MLsol */
    static std::unique_ptr<Writer> build(const WriterEnum format, const MultiLevelSolution * ml_sol);

    /** runtime selection of writer for MLmesh */
    static std::unique_ptr<Writer> build(const WriterEnum format, const MultiLevelMesh * ml_mesh);

// === Constructors / Destructor  - END =================



// === Write - BEGIN =================
    
  public:
    /** at finest level: write output function */
    virtual void Write(const std::string output_path,
                       const std::string order,
                       const std::vector < std::string > & vars = std::vector < std::string > (), 
                       const unsigned time_step = _time_step_index_default)  = 0;
    
    /** at finest level: write output function with fixed level and arbitrary initial string */
    virtual void Write(const std::string init_string,
                       const std::string output_path, 
                       const std::string order,
                       const std::vector < std::string >& vars = std::vector < std::string > (), 
                       const unsigned time_step = _time_step_index_default)  = 0;
                       
    /** at finest level:    */
    virtual void Write(const std::string init_string,
                       const std::string output_path, 
                       const std::string suffix_pre_extension, 
                       const std::string order,
                       const std::vector < std::string >& vars = std::vector < std::string > (), 
                       const unsigned time_step = _time_step_index_default)  = 0;

    /**  at arbitrary level: write output function with arbitrary level */
    virtual void Write(const unsigned my_level, 
                       const std::string output_path, 
                       const std::string order,
                       const std::vector < std::string >& vars = std::vector < std::string > (), 
                       const unsigned time_step = _time_step_index_default)  = 0;
  
    /**  at arbitrary level: write output function with arbitrary level and arbitrary initial string and arbitrary suffix before the extension */
    virtual void Write(const unsigned my_level, 
                       const std::string init_string, 
                       const std::string output_path,
                       const std::string suffix_pre_extension, 
                       const std::string order,
                       const std::vector < std::string >& vars = std::vector < std::string > (), 
                       const unsigned time_step = _time_step_index_default)  = 0;

// === Write - END =================




// === Solution - BEGIN =================
  protected:

    const Solution * get_solution(const unsigned level_in) const;
    
    /** the multilevelsolution pointer: it is const, so it does not modify the object that is printed */
    const  MultiLevelSolution* _ml_sol;
// === Solution - END =================

// === Mesh - BEGIN =================
  protected:

    const Mesh * get_mesh(const unsigned level_in) const;
    
    /** the multilevel mesh: it is const, so it does not modify the object that is printed */
    const MultiLevelMesh* _ml_mesh;
// === Mesh - END =================


// === Mesh, Level - BEGIN =================
  protected:

    /** Number of mesh levels */
    const unsigned int _gridn;
// === Mesh, Level - END =================

  };

  
  
  
  
} //end namespace femus



#endif
