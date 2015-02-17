/*=========================================================================

Program: FEMuS
Module: MultiLevelProblem
Authors: Eugenio Aulisa, Simone Bnà

Copyright (c) FEMuS
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __MultiLevelSolution_hpp__
#define __MultiLevelSolution_hpp__

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




/**
 * This class is a black box container to handle multilevel solutions.
 */

class MultiLevelSolution : public ParallelObject {

private:
  
    typedef double (*initfunc) (const double &x, const double &y, const double &z);

public:

    /** Constructor */
    MultiLevelSolution(MultiLevelMesh *ml_msh);

    /** Destructor */
    ~MultiLevelSolution();

    /** To be Added */
    void AddSolution(const char name[], const FEFamily fefamily, const FEOrder order, unsigned tmorder=0, const bool &Pde_type=1);

    /** To be Added */
    void AddSolutionLevel();
    
    /** To be Added */
    void AssociatePropertyToSolution(const char solution_name[], const char solution_property[]);

    /** To be Added */
    void ResizeSolutionVector( const char name[]);

    /** To be Added */
    void Initialize(const char name[], initfunc func = NULL);

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
    void BuildProlongatorMatrix(unsigned gridf, unsigned SolIndex);

    /** To be Added */
    void AttachSetBoundaryConditionFunction ( bool (* SetBoundaryConditionFunction) (const double &x, const double &y, const double &z,const char name[],
            double &value, const int FaceName, const double time) );

    /** To be Added */
    void GenerateBdc(const char name[], const char bdc_type[]="Steady");
    
    /** To be Added */
    void InitializeBdc();

    /** To be Added */
    void UpdateBdc(const double time);

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
    bool  TestIfSolutionIsDisplacemenet(unsigned i) {
        return _TestIfDisplacement[i];
    };

    // member data
    MultiLevelMesh* _ml_msh; //< Multilevel mesh

    bool (*_SetBoundaryConditionFunction) (const double &x, const double &y, const double &z,const char name[],
                                           double &value, const int FaceName, const double time); //< boundary condition function pointer
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
    
    /** This group of vectors has the size of the number of added solutions */
    vector< vector <BDCType> > _boundaryconditions;
    vector< vector <bool> > _ishomogeneous;
    vector< vector <FunctionBase *> > _nonhomogeneousbcfunction; 
    
    bool _bdc_func_set;
    
    unsigned short  _gridn;
    vector <int>    _SolType;
    vector <FEFamily> _family;
    vector <FEOrder> _order;
    vector <char*>  _SolName;
    vector <char*>  _BdcType;
    vector <int>    _SolTmorder;
    vector <bool>   _PdeType;
    vector <bool>   _TestIfPressure;
    vector <bool>   _TestIfDisplacement;
    
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

