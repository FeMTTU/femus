/*=========================================================================

Program: FEMuS
Module: MultiLevelProblem
Authors: Eugenio Aulisa, Simone Bnà, Giorgio Bornia

Copyright (c) FEMuS
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_solution_MultiLevelSolution_hpp__
#define __femus_solution_MultiLevelSolution_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <vector>
#include <memory>
#include "MultiLevelMesh.hpp"
#include "Solution.hpp"
#include "ParallelObject.hpp"
#include "FElemTypeEnum.hpp"
#include "BDCTypeEnum.hpp"
#include "FunctionBase.hpp"
#include "Writer.hpp"

namespace femus {


class MultiLevelProblem;

/**
 * This class is a black box container to handle multilevel solutions.
 */

class MultiLevelSolution : public ParallelObject {

private:

    /** Initial condition function pointer typedef */
    typedef double (*InitFunc) (const std::vector < double >& x);

    /** duplicate */
    typedef double (*InitFuncMLProb) (const MultiLevelProblem * ml_prob, const double &x, const double &y, const double &z, const char * name);

    /** Boundary condition function pointer typedef */
    typedef bool (*BoundaryFunc) (const std::vector < double >& x, const char name[], double &value, const int FaceName, const double time);

    /** duplicate */
    typedef bool (*BoundaryFuncMLProb) (const MultiLevelProblem * ml_prob, const double &x, const double &y, const double &z,const char name[], double &value, const int FaceName, const double time);

public:

    /** Constructor */
    MultiLevelSolution(MultiLevelMesh *ml_msh);

    /** Destructor */
    ~MultiLevelSolution();

    /** To be Added */
    void AddSolution(const char name[], const FEFamily fefamily, const FEOrder order, unsigned tmorder=0, const bool &Pde_type=1);

    /** If you want to add a vector whose components are treated the same way */
    void AddSolutionVector(const unsigned n_components, const std::string name, const FEFamily fefamily, const FEOrder order, unsigned tmorder=0, const bool &Pde_type=1);

    /** To be Added */
    void AddSolutionLevel();

    /** To be Added */
    void AssociatePropertyToSolution(const char solution_name[], const char solution_property[]);

    /** To be Added */
    void PairSolution(const char solution_name[], const char solution_pair[]);

    /** To be Added */
    void ResizeSolutionVector( const char name[]);

    /** To be Added */
    void Initialize(const char name[], InitFunc func = NULL);

    /** duplicate */
    void InitializeMLProb(const MultiLevelProblem * ml_prob, const char * name, InitFuncMLProb func = NULL);

    /** To be Added */
    unsigned GetIndex(const char name[]) const;

    /** To be Added */
    unsigned GetSolType(const char name[]);

    /** To be Added */
    unsigned GetSolutionSize() {
        return _SolType.size();
    };

    /** To be Added */
    vector <char*>  GetSolName() {
        return _SolName;
    };

    /** To be Added */
    vector <int>  GetSolType() {
        return _SolType;
    };

    /** To be Added */
    void AttachSetBoundaryConditionFunction (BoundaryFunc SetBoundaryConditionFunction );
    void FixSolutionAtOnePoint( const char sol[] ){ _FixSolutionAtOnePoint[GetIndex(sol)] = true ;};


    /** To be Added */
    void GenerateBdc(const char name[], const char bdc_type[]="Steady", const MultiLevelProblem * ml_prob = NULL);

    /** To be Added */
    void InitializeBdc();

    /** To be Added */
    void UpdateBdc(const double time);

    /** duplicate */
    void AttachSetBoundaryConditionFunctionMLProb ( BoundaryFuncMLProb SetBoundaryConditionFunction_in );

    /** duplicate */
    void GenerateBdcMLProb(const MultiLevelProblem * ml_prob, const unsigned int k, const unsigned int grid0, const double time);

    /** duplicate */
    BoundaryFuncMLProb _SetBoundaryConditionFunctionMLProb;

    /** To be Added */
    void GenerateBdc(const unsigned int k, const unsigned grid0, const double time);

    /** To be Added */
    void GenerateBdc_new(const unsigned int k, const unsigned grid0, const double time);

    /** To be Added */
    BDCType GetBoundaryCondition(const std::string varname, const unsigned int facename) const;

    /** To be Added */
    bool Ishomogeneous(const std::string varname, const unsigned int facename) const;

    /** To be Added */
    void SetBoundaryCondition_new(const std::string name, const std::string facename, const BDCType bdctype = DIRICHLET,
                              const bool istimedependent = false, FunctionBase* func = NULL);

    /** To be Added */
    FunctionBase* GetBdcFunction(const std::string varname, const unsigned int facename) const;

    /** To be Added */
    Solution* GetSolutionLevel(const unsigned i) {
        return _solution[i];
    };

    /** To be Added */
    char* GetSolutionName(unsigned i) {
        return _SolName[i];
    };

    /** To be Added */
    int   GetSolutionType(unsigned i) {
        return _SolType[i];
    };

    /** To be Added */
    unsigned GetSolutionType(const char name[]);

    /** To be Added */
    char* GetBdcType(unsigned i) {
        return _BdcType[i];
    };

    /** To be Added */
    int   GetSolutionTimeOrder(unsigned i) {
        return _SolTmorder[i];
    };

    /** To be Added */
    bool  TestIfSolutionIsPressure(unsigned i) {
        return _TestIfPressure[i];
    };

    /** To be Added */
    unsigned GetSolutionPairIndex(const unsigned& i) const{
      return _SolPairIndex[i];
    }

    // member data
    MultiLevelMesh* _ml_msh; //< Multilevel mesh

    /** boundary condition function pointer */
    BoundaryFunc _SetBoundaryConditionFunction;

    void build();

    bool _Use_GenerateBdc_new;

    /** To be Added */
    Writer* GetWriter() {return _writer; }

    /** To be Added */
    const Writer* GetWriter() const {return _writer; }

    /** To be Added */
    void SetWriter(const WriterEnum format) { _writer = Writer::build(format,this).release(); }

private:

    /** To be Added */
    BDCType GetBoundaryCondition(const unsigned int var, const unsigned int facename) const;

    /** To be Added */
    bool Ishomogeneous(const unsigned int var, const unsigned int facename) const;

    /** To be Added */
    FunctionBase* GetBdcFunction(const unsigned int var, const unsigned int facename) const;

    /** Array of solution, dimension number of levels */
    vector <Solution*>  _solution;

    /** Flag to tell whether the BC function has been set */
    bool _bdc_func_set;

    /** This group of vectors has the size of the number of added solutions */
    vector< vector <BDCType> > _boundaryconditions;
    vector< vector <bool> > _ishomogeneous;
    vector< vector <FunctionBase *> > _nonhomogeneousbcfunction;

    unsigned short  _gridn;
    vector < int >    _SolType;    /* Tells the FE index */
    vector < FEFamily > _family;
    vector < FEOrder > _order;
    vector < char* >  _SolName;
    vector < char* >  _BdcType;
    vector < int >    _SolTmorder;
    vector < bool >   _PdeType;    /*Tells whether the Solution is an unknown of a PDE or not*/
    vector < bool >   _TestIfPressure;
    vector < bool > _FixSolutionAtOnePoint;
    vector <unsigned> _SolPairIndex;

    /** Multilevel solution writer */
    Writer* _writer;

};


inline
BDCType MultiLevelSolution::GetBoundaryCondition(const unsigned int var, const unsigned int facename) const {
    return _boundaryconditions[var][facename];
}

inline
bool MultiLevelSolution::Ishomogeneous(const unsigned int var, const unsigned int facename) const {
    return _ishomogeneous[var][facename];
}

inline
FunctionBase* MultiLevelSolution::GetBdcFunction(const unsigned int var, const unsigned int facename) const {
    return _nonhomogeneousbcfunction[var][facename];
}

inline
BDCType MultiLevelSolution::GetBoundaryCondition(const std::string varname, const unsigned int facename) const {
    unsigned int var = GetIndex(varname.c_str());
    return _boundaryconditions[var][facename];
}

inline
bool MultiLevelSolution::Ishomogeneous(const std::string varname, const unsigned int facename) const {
    unsigned int var = GetIndex(varname.c_str());
    return _ishomogeneous[var][facename];
}

inline
FunctionBase* MultiLevelSolution::GetBdcFunction(const std::string varname, const unsigned int facename) const {
    unsigned int var = GetIndex(varname.c_str());
    return _nonhomogeneousbcfunction[var][facename];
}


} //end namespace femus



#endif

