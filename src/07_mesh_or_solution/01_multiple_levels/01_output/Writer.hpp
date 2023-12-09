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
  

  class Writer : public ParallelObject {


// === Constructors / Destructor  - BEGIN =================
  public:
    
    /** Constructor. */
    Writer(const MultiLevelSolution * ml_sol);

    /** Constructor. */
    Writer(const MultiLevelMesh * ml_mesh);

    /** Destructor */
    virtual ~Writer();

    
    /** runtime selection of writer for MLsol */
    static std::unique_ptr<Writer> build(const WriterEnum format, const MultiLevelSolution * ml_sol);

    /** runtime selection of writer for MLmesh */
    static std::unique_ptr<Writer> build(const WriterEnum format, const MultiLevelMesh * ml_mesh);

  private:
    
    void initialize_flags();
// === Constructors / Destructor  - END =================

    
// === Write - BEGIN =================
  public:
    
    /** write output function */
    virtual void Write(const std::string output_path,
                       const std::string order,
                       const std::vector < std::string > & vars = std::vector < std::string > (), 
                       const unsigned time_step = 0)  = 0;
    
    /** write output function with arbitrary level */
    virtual void Write(const unsigned my_level, 
                       const std::string output_path, 
                       const std::string order,
                       const std::vector < std::string >& vars = std::vector < std::string > (), 
                       const unsigned time_step = 0)
        { abort(); };
  
    /** write output function with fixed level and arbitrary initial string */
    virtual void Write(const std::string init_string,
                       const std::string output_path, 
                       const std::string order,
                       const std::vector < std::string >& vars = std::vector < std::string > (), 
                       const unsigned time_step = 0) 
        { abort(); };
  
    /** write output function with arbitrary level and arbitrary initial string and arbitrary suffix before the extension */
    virtual void Write(const unsigned my_level, 
                       const std::string init_string, 
                       const std::string output_path,
                       const std::string suffix_pre_extension, 
                       const std::string order,
                       const std::vector < std::string >& vars = std::vector < std::string > (), 
                       const unsigned time_step = 0)
        { abort(); };
// === Write - END =================

    
    
// === Debug - BEGIN =================
  public:
    
      /** Set if to print or not to print the debugging variables */
      void SetDebugOutput( bool value ) {
        _debugOutput = value;
      }
    
  protected:
    
      bool _debugOutput;
      
// === Debug, Solutions that are Unknowns - BEGIN =================
    unsigned compute_sol_bdc_res_eps_size(const Solution * solution, const unsigned i) const;
    
    std::string print_sol_bdc_res_eps_name(const std::string solName, const unsigned name) const;

    static const std::string _name_bdc;
    static const std::string _name_res;
    static const std::string _name_eps;

    static constexpr unsigned _index_sol = 0;
    static constexpr unsigned _index_bdc = 1;
    static constexpr unsigned _index_res = 2;
    static constexpr unsigned _index_eps = 3;
    
// === Debug, Solutions that are Unknowns - END =================

// === Debug - END =================
    
    
// === Moving Mesh - BEGIN =================
  public:
    
    /** set moving mesh */
    void SetMovingMesh(std::vector<std::string>& movvars_in);

  protected:
    
    /** a flag to move the output mesh */
    bool _moving_mesh;
    
    /** the displacement variables for moving mesh */
    std::vector<std::string> _moving_vars;
// === Moving Mesh - END =================
    

// === Graph Variables - BEGIN =================
  public:
    
    void SetGraphVariable(const std::string &GraphVaraible);
    void UnsetGraphVariable(){ _graph = false;};

  protected:
    
    bool _graph;
    std::string _graphVariable;

// === Graph Variables - END =================

    
// === Surface Variable - BEGIN =================
  public:
    
    void SetSurfaceVariables( std::vector < std::string > &surfaceVariable );
    void UnsetSurfaceVariables(){ _surface = false;};
    
  protected:
    
    bool _surface;
    std::vector < std::string > _surfaceVariables;
// === Surface Variable - END =================
    


  protected:

// === Mesh or Solution - BEGIN =================
    std::string get_filename_prefix() const; 
// === Mesh or Solution - END =================

    
// === Solution, FE index for printing - BEGIN =================
    unsigned fe_index(const std::string & order_str) const;
// === Solution, FE index for printing - END =================

// === Solution - BEGIN =================
    /** the multilevelsolution pointer: it is const, so it does not modify the object that is printed */
    const  MultiLevelSolution* _ml_sol;
// === Solution - END =================

// === Mesh - BEGIN =================
    /** the multilevel mesh: it is const, so it does not modify the object that is printed */
    const MultiLevelMesh* _ml_mesh;
// === Mesh - END =================


// === Mesh, Level - BEGIN =================
    int _gridn;
// === Mesh, Level - END =================


// === Geometric Element, Connectivities - BEGIN =================
    /** map from femus connectivity to vtk-connectivity for paraview visualization */
    static const unsigned FemusToVTKorToXDMFConn[27];
// === Geometric Element, Connectivities - END =================



// === FE DOFMAP & PROJECTION at SAME LEVEL (needed for node-based printing, Only Lagrange) - BEGIN =================

    FE_Proj_Matrices   _fe_proj_matrices;
   
// === FE DOFMAP & PROJECTION at SAME LEVEL (needed for node-based printing, Only Lagrange) - END =================

    
  };

  
  
  
  
} //end namespace femus



#endif
