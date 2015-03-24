/*=========================================================================

 Program: FEMUS
 Module: Writer
 Authors: Eugenio Aulisa, Simone Bnà

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
#include <vector>
#include <string>
#include <memory>
#include "ParallelObject.hpp"
#include "WriterEnum.hpp"

namespace femus {

  // map from our connectivity to vtk-connectivity for paraview visualization  //TODO move this to the appropriate place
  const unsigned map_pr[27] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,23,21,20,22,24,25,26};

  //------------------------------------------------------------------------------
  // Forward declarations
  //------------------------------------------------------------------------------
  class MultiLevelMesh;
  class MultiLevelSolution;
  class SparseMatrix;
  class Vector;


  class Writer : public ParallelObject {

  public:

    /** Constructor. */
    Writer(MultiLevelSolution * ml_sol);

    /** Constructor. */
    Writer(MultiLevelMesh * ml_mesh);

    /** Destructor */
    virtual ~Writer();

    /** write output function */
    virtual void write(const std::string output_path, const char order[], const std::vector < std::string > & vars = std::vector < std::string > (), const unsigned time_step = 0)  const = 0;
    virtual void Pwrite(const std::string output_path, const char order[], const std::vector < std::string > & vars = std::vector < std::string > (), const unsigned time_step = 0) const = 0;
    /** set moving mesh */
    void SetMovingMesh(std::vector<std::string>& movvars_in);
    
    /** runtime selection of writer for MLsol */
    static std::auto_ptr<Writer> build(const WriterEnum format, MultiLevelSolution * ml_sol);

    /** runtime selection of writer for MLmesh */
    static std::auto_ptr<Writer> build(const WriterEnum format, MultiLevelMesh * ml_mesh);
    
  protected:

    /** a flag to move the output mesh */
    int _moving_mesh;

    /** the displacement variables for mesh moving */
    std::vector<std::string> _moving_vars;

    /** the multilevelsolution pointer */
    MultiLevelSolution* _ml_sol;

    /** the multilevel mesh */
    MultiLevelMesh* _ml_mesh;

    int _gridn;

    int _gridr;


  private:

  };
  
} //end namespace femus



#endif
