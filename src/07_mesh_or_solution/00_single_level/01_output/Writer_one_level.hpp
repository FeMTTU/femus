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

#ifndef __femus_solution_Writer_one_level_hpp__
#define __femus_solution_Writer_one_level_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
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
  class Mesh;
  class Solution;
  

  class Writer_one_level : public ParallelObject {


// === Constructors / Destructor  - BEGIN =================
  public:
    
    Writer_one_level();

// // //     /** Constructor. */
// // //     Writer_one_level(const Solution * ml_sol);
// // // 
// // //     /** Constructor. */
// // //     Writer_one_level(const Mesh * ml_mesh);
// // // 
// // //     
// // //     /** runtime selection of writer for MLsol */
// // //     static std::unique_ptr< Writer_one_level > build(const WriterEnum format, const Solution * ml_sol);
// // // 
// // //     /** runtime selection of writer for MLmesh */
// // //     static std::unique_ptr< Writer_one_level > build(const WriterEnum format, const Mesh * ml_mesh);

  private:
    
    void initialize_flags();
// === Constructors / Destructor  - END =================

    
// === Write - BEGIN =================
  public:
    
    // /** write output function with arbitrary level */
    // virtual void Write(const unsigned my_level, 
    //                    const std::string output_path, 
    //                    const std::string order,
    //                    const std::vector < std::string >& vars = std::vector < std::string > (), 
    //                    const unsigned time_step = _time_step_index_default)  = 0;
  
    /** write output function with arbitrary level and arbitrary initial string and arbitrary suffix before the extension */
    virtual void Write(const unsigned my_level, 
                       const std::string filename_prefix, 
                       const std::string output_path,
                       const std::string suffix_pre_extension, 
                       const std::string order,
                       const std::vector < std::string >& vars = std::vector < std::string > (), 
                       const unsigned time_step = _time_step_index_default)  = 0;

// === Write - END =================



// === Mesh or Solution - BEGIN =================
  protected:
    std::string get_filename_prefix(const Solution * solution) const; 
// === Mesh or Solution - END =================



// === Mesh - BEGIN =================
  protected:
    /** the mesh: it is const, so it does not modify the object that is printed */
    const Mesh* _mesh;
// === Mesh - END =================


// === Mesh, Geometric Element, Connectivities - BEGIN =================
  protected:
    /** map from femus connectivity to vtk-connectivity for paraview visualization */
    static const unsigned FemusToVTKorToXDMFConn[27];
// === Mesh, Geometric Element, Connectivities - END =================



// === Mesh, displacement - BEGIN =================


// === Mesh displacement: Moving Mesh - BEGIN =================
  public:
    
    /** set moving mesh */
    void SetMovingMesh(std::vector<std::string>& movvars_in);

  protected:
    
    /** a flag to move the output mesh */
    bool _moving_mesh;
    
    /** the displacement variables for moving mesh */
    std::vector<std::string> _moving_vars;
// === Mesh displacement: Moving Mesh - END =================
    

// === Mesh displacement: Graph Variables - BEGIN =================
  public:
    
    void SetGraphVariable(const std::string &GraphVaraible);
    void UnsetGraphVariable(){ _graph = false;};

  protected:
    
    bool _graph;
    std::string _graphVariable;

// === Mesh displacement: Graph Variables - END =================

    
// === Mesh displacement: Surface Variable - BEGIN =================
  public:
    
    ///@todo seems to be unused
    void SetSurfaceVariables( std::vector < std::string > &surfaceVariable );
    ///@todo seems to be unused
    void UnsetSurfaceVariables(){ _surface = false;};
    
  protected:
    
    bool _surface;
    std::vector < std::string > _surfaceVariables;
// === Mesh displacement: Surface Variable - END =================


// === Mesh, displacement - END =================


// === Solution - BEGIN =================
  protected:
    
// === Solution, FE index for printing - BEGIN =================
    unsigned fe_index(const std::string & order_str) const;
// === Solution, FE index for printing - END =================

    /** thesolution pointer: it is const, so it does not modify the object that is printed */
    const  Solution* _sol;
// === Solution - END =================

// === Solution, FE DOFMAP & PROJECTION at SAME LEVEL (needed for node-based printing, Only Lagrange) - BEGIN =================
  protected:

    FE_Proj_Matrices   _fe_proj_matrices;
   
// === Solution, FE DOFMAP & PROJECTION at SAME LEVEL (needed for node-based printing, Only Lagrange) - END =================

    
// === Solution, Debug, Solutions that are Unknowns - BEGIN =================
  public:
    
      /** Set if to print or not to print the debugging variables */
      void SetDebugOutput( bool value ) {
        _debugOutput = value;
      }
    
  protected:
    
      bool _debugOutput;
      
      
    unsigned compute_sol_bdc_res_eps_size(const Solution * solution, const unsigned i) const;
    
    std::string print_sol_bdc_res_eps_name(const std::string solName, const unsigned name) const;

    static const std::string _name_bdc;
    static const std::string _name_res;
    static const std::string _name_eps;

    static constexpr unsigned _index_sol = 0;
    static constexpr unsigned _index_bdc = 1;
    static constexpr unsigned _index_res = 2;
    static constexpr unsigned _index_eps = 3;

// === Solution, Debug, Solutions that are Unknowns - END =================
    


// === Time step, default - BEGIN =================
 protected:
    
    static constexpr unsigned _time_step_index_default = 0;
    
// === Time step, default - END =================
  


  };

  
  
  
  
} //end namespace femus



#endif
