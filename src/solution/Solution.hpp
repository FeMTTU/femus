/*=========================================================================

 Program: FEMuS
 Module: Solution
 Authors: Simone Bnà

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * This class is useful to handle multilevel variables. It's going to be splitted
 * into meshvariables and solutionvariables
*/


#ifndef __solution_hpp__
#define __solution_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Mesh.hpp"
#include "petscmat.h"
#include "FElemTypeEnum.hpp"
#include "ParallelObject.hpp"

namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class elem_type;
class NumericVector;
class SparseMatrix;
class mesh;

using std::vector;

class Solution : public ParallelObject {

public:

    /** Constructor */
    Solution(mesh *other_msh);

    /** Destructor */
    ~Solution();

    /** Add a new variable called 'name' */
    void AddSolution( const char name[], const FEFamily fefamily, const FEOrder order, const unsigned& tmorder=0, const bool &Pde_type=1);

    /** Resize the solution vector */
    void ResizeSolutionVector(const char name[]);

    /** Free the solution vectors */
    void FreeSolutionVectors();

    /** Set the coarse coordinates */
    void SetCoarseCoordinates( vector < vector < double> > &vt);

    /** Sum to Solution vector the Epsilon vector. It is used inside the multigrid cycle */
    void SumEpsToSol(const vector <unsigned> &_SolPdeIndex,  NumericVector* EPS, NumericVector* RES, const vector <vector <unsigned> > &KKoffset);

    /** Update the solution */
    void UpdateSolution();

    /** Get the solution (Numeric Vector) by name */
    const NumericVector* GetSolutionName(const char* var) {
        return _Sol[GetIndex(var)];
    };

    /** member data - one for each variable - */
    vector <NumericVector*> _Sol;
    vector <NumericVector*> _SolOld;
    vector <NumericVector*> _Res;
    vector <NumericVector*> _Eps;
    vector <NumericVector*> _Bdc;
    vector <bool> _ResEpsBdcFlag;

    /** one for every type of variable */
    SparseMatrix* _ProjMat[5];
    bool _ProjMatFlag[5];


private:

    /** Get the index of the variable -name- */
    unsigned GetIndex(const char name[]) const;

    //member data
    vector <int> _SolType;
    vector <char*> _SolName;
    vector <unsigned> _SolTmOrder;
    vector <FEFamily> _family;
    vector <FEOrder> _order;
    mesh *_msh;

};


} //end namespace femus



#endif
