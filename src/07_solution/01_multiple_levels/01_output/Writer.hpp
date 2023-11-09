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
  

  class Writer : public ParallelObject {


// === Constructors / Destructor  - BEGIN =================
  public:
    
    /** Constructor. */
    Writer(MultiLevelSolution * ml_sol);

    /** Constructor. */
    Writer(MultiLevelMesh * ml_mesh);

    /** Destructor */
    virtual ~Writer();

    
    /** runtime selection of writer for MLsol */
    static std::unique_ptr<Writer> build(const WriterEnum format, MultiLevelSolution * ml_sol);

    /** runtime selection of writer for MLmesh */
    static std::unique_ptr<Writer> build(const WriterEnum format, MultiLevelMesh * ml_mesh);

// === Constructors / Destructor  - END =================

    
// === Write - BEGIN =================
  public:
    
    /** write output function */
    virtual void Write(const std::string output_path, const char order[], const std::vector < std::string > & vars = std::vector < std::string > (), const unsigned time_step = 0)  = 0;
    
    /** write output function with arbitrary level */
    virtual void Write(const unsigned my_level, const std::string output_path, const char order[], const std::vector < std::string >& vars = std::vector < std::string > (), const unsigned time_step = 0) { abort(); };
  
    /** write output function with fixed level and arbitrary initial string */
    virtual void Write(const std::string init_string, const std::string output_path, const char order[], const std::vector < std::string >& vars = std::vector < std::string > (), const unsigned time_step = 0) { abort(); };
  
    /** write output function with arbitrary level and arbitrary initial string and arbitrary suffix before the extension */
    virtual void Write(const unsigned my_level, const std::string init_string, const std::string output_path, const std::string suffix_pre_extension, const char order[], const std::vector < std::string >& vars = std::vector < std::string > (), const unsigned time_step = 0) { abort(); };
// === Write - END =================

    
    
// === Debug - BEGIN =================
  public:
    
      /** Set if to print or not to print the debugging variables */
      void SetDebugOutput( bool value ) {
        _debugOutput = value;
      }
    
  protected:
    
      bool _debugOutput;
// === Debug - END =================
    
    
// === Moving Mesh - BEGIN =================
  public:
    
    /** set moving mesh */
    void SetMovingMesh(std::vector<std::string>& movvars_in);

  protected:
    
    /** a flag to move the output mesh */
    int _moving_mesh;
    
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

// === Solution - BEGIN =================
    /** the multilevelsolution pointer */
    MultiLevelSolution* _ml_sol;
// === Solution - END =================

// === Mesh - BEGIN =================
    /** the multilevel mesh */
    MultiLevelMesh* _ml_mesh;
// === Mesh - END =================

    int _gridn;

    /** map from femus connectivity to vtk-connectivity for paraview visualization */
    static const unsigned FemusToVTKorToXDMFConn[27];



  private:

  };

} //end namespace femus



#endif
