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

#ifndef __writer_h_
#define __writer_h_

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
class MultiLevelSolution;
class SparseMatrix;
class Vector;


class Writer : public ParallelObject {

public:

    /** Constructor. */
    Writer(MultiLevelSolution & ml_probl);

    /** Destructor */
    virtual ~Writer();

    /** write output function */
    virtual void write_system_solutions(const std::string output_path, const char order[], std::vector<std::string>& vars, const unsigned time_step = 0) = 0;

    /** set moving mesh */
    void SetMovingMesh(std::vector<std::string>& movvars_in);
    
    static std::auto_ptr<Writer> build(const WriterEnum format, MultiLevelSolution * ml_sol);

protected:

    /** a set of matrix for the projection of the solution */
    static std::vector<SparseMatrix*> _ProlQitoQj[3][3];

    /** a flag to move the output mesh */
    int _moving_mesh;

    /** the displacement variables for mesh moving */
    std::vector<std::string> _moving_vars;

    /** the multilevelsolution reference */
    MultiLevelSolution& _ml_sol;

    int _gridn;

    int _gridr;


private:

    /** This routine generates the matrices for the projection of the solutions from different FE spaces */
    void BuildProlongatorMatrices();

};



} //end namespace femus



#endif
