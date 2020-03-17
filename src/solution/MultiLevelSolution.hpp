/*=========================================================================

Program: FEMuS
Module: MultiLevelSolution
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


public:
    
    /** Initial condition function pointer typedef */
    typedef double (*InitFunc) (const std::vector < double >& x);

    /** duplicate */
    typedef double (*InitFuncMLProb) (const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]);

    /** Boundary condition function pointer typedef */
    typedef bool (*BoundaryFunc) (const std::vector < double >& x, const char name[], double &value, const int FaceName, const double time);
    
    /** duplicate */
    typedef bool (*BoundaryFuncMLProb) (const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double &value, const int FaceName, const double time);

    /** Constructor */
    MultiLevelSolution(MultiLevelMesh *ml_msh);

    /** Destructor */
    ~MultiLevelSolution();

    /** this is the destructor that can be called explicitly, instead of the automatic destructor */
    void clear();
 
    /** To be Added */
    void AddSolution(const char name[], const FEFamily fefamily, const FEOrder order, unsigned tmorder = 0, const bool &Pde_type = 1);

    /** If you want to add a vector whose components are treated the same way */
    void AddSolutionVector(const unsigned n_components, const std::string name, const FEFamily fefamily, const FEOrder order, unsigned tmorder=0, const bool &Pde_type=1);

    /** To be Added */
    void AddSolutionLevel();

    /** To be Added */
    void ResizeSolution_par(const unsigned new_size);
    
    /** To be Added */
    void AssociatePropertyToSolution(const char solution_name[], const char solution_property[], const bool &bool_property = true);

    /** To be Added */
    void PairSolution(const char solution_name[], const char solution_pair[]);

    /** To be Added */
    void ResizeSolutionVector( const char name[]);

    /** To be Added */
    void Initialize(const char name[], InitFunc func = NULL);

    void Initialize(const char * name, InitFuncMLProb func, const MultiLevelProblem * ml_prob);
    
    void Initialize(const char name[], InitFunc func, InitFuncMLProb funcMLProb, const MultiLevelProblem *ml_prob);

    inline void Set(const char name[], InitFuncMLProb funcMLProb, const MultiLevelProblem *ml_prob);
    
    /** To be Added */
    unsigned GetIndex(const char name[]) const;

    /** To be Added */
    unsigned GetSolType(const char name[]);

    /** To be Added */
    unsigned GetSolutionSize() {
        return _solType.size();
    };

    /** To be Added */
    const unsigned GetSolutionSize() const {
        return _solType.size();
    };
    
    /** To be Added */
    vector <char*>  GetSolName() {
        return _solName;
    };

    /** To be Added */
    vector <int>  GetSolType() {
        return _solType;
    };

    /** To be Added */
    void AttachSetBoundaryConditionFunction( BoundaryFunc SetBoundaryConditionFunction );
    
    void AttachSetBoundaryConditionFunction( BoundaryFuncMLProb SetBoundaryConditionFunction );
    
    void FixSolutionAtOnePoint( const char sol[] ){
      _fixSolutionAtOnePoint[GetIndex(sol)] = true ;
      for(unsigned ig = 1; ig < _gridn; ig++){
        _solution[ig]->RemoveNullSpace(GetIndex(sol));
      }
    };

    
    

    /** To be Added */
    void GenerateBdc(const char name[], const char bdc_type[]="Steady", const MultiLevelProblem * ml_prob = NULL);

    /** To be Added */
    void InitializeBdc();

    /** To be Added */
    void UpdateBdc(const double time);

    /** To be Added */
    void GenerateBdc( const unsigned int k, const unsigned grid0, const double time );
    void GenerateRKBdc(const unsigned int &solIndex, const std::vector<unsigned> &solKiIndex, 
                       const unsigned int &grid0, const std::vector < double> & time,  const double &time0, 
                       const double &dt, const double* AI);
    
    //for NONLOCAL problems, _Bdc must be 0 on all the volume constraint
    void GenerateBdcOnVolumeConstraint(const std::vector<unsigned> &volumeConstraintFlags, const unsigned &solIndex, const unsigned &grid0);

    /** To be Added */
    BDCType GetBoundaryCondition(const std::string varname, const unsigned int facename) const;

    /** To be Added */
    bool Ishomogeneous(const std::string varname, const unsigned int facename) const;

    /** To be Added */
    void SetBoundaryCondition_new(const std::string name, const std::string facename, const BDCType bdctype = DIRICHLET,
                              const bool istimedependent = false, FunctionBase* func = NULL);

    /** To be Added */
    FunctionBase* GetBdcFunction(const std::string varname, const unsigned int facename) const;

       /** duplicate of GetSolutionLevel, to be removed @todo */
    Solution* GetLevel(const unsigned i) {
        return _solution[i];
    };
    
       /** To be Added */
    Solution* GetSolutionLevel(const unsigned i) {
        return _solution[i];
    };

    /** To be Added */
    const Solution* GetSolutionLevel(const unsigned i) const {
        return _solution[i];
    };
    
    /** To be Added */
    char* GetSolutionName(unsigned i) {
        return _solName[i];
    };

    /** To be Added */
    int GetSolutionType(unsigned i) {
        return _solType[i];
    };

    /** To be Added */
    const int GetSolutionType(unsigned i) const {
        return _solType[i];
    };
    
    /** To be Added */
    unsigned GetSolutionType(const char name[]);

    /** To be Added */
    char* GetBdcType(unsigned i) {
        return _bdcType[i];
    };

    /** To be Added */
    int   GetSolutionTimeOrder(unsigned i) {
        return _solTimeOrder[i];
    };

    const int   GetSolutionTimeOrder(const std::string & sol_name) const;

    /** To be Added */
    bool  TestIfSolutionIsPressure(unsigned i) {
        return _testIfPressure[i];
    };

    /** To be Added */
    unsigned GetSolutionPairIndex(const unsigned& i) const{
      return _solPairIndex[i];
    }
    
    unsigned GetSolutionPairInverseIndex(const unsigned& i) const{
      return _solPairInverseIndex[i];
    }
    
    
    FEFamily GetSolutionFamily(const unsigned& i){
      return _family[i];  
    };
    
    const FEFamily GetSolutionFamily(const std::string & sol_name) const;
    
    FEOrder GetSolutionOrder(const unsigned& i){
      return _order[i];    
    }
    
    const FEOrder GetSolutionOrder(const std::string & sol_name) const;

    void build();

    /** To be Added */
    Writer* GetWriter() {return _writer; }

    /** To be Added */
    const Writer* GetWriter() const {return _writer; }

    /** To be Added */
    void SetWriter(const WriterEnum format) { _writer = Writer::build(format,this).release(); }

    // member data
    MultiLevelMesh* _mlMesh; //< Multilevel mesh

    BoundaryFunc GetBdcFunction(){
      return _SetBoundaryConditionFunction;
    }

    BoundaryFuncMLProb GetBdcFunctionMLProb(){
      return _SetBoundaryConditionFunctionMLProb;
    }

    bool _useParsedBCFunction;

    void SaveSolution(const char* filename, const double &time=0.);
    void SaveSolution(const char* filename, const unsigned &iteration);
    void LoadSolution(const char* filename);
    void LoadSolution(const unsigned &level, const char* filename);
    
     // *******************************************************

    void RefineSolution( const unsigned &gridf );
    void CoarsenSolutionByOneLevel_wrong( const unsigned &gridf );
    void CoarsenSolutionByOneLevel( const unsigned &gridf );

  // ********************************************
    
    void UpdateSolution(const char name[], InitFunc func, const double& time);
    
    void fill_at_level_from_level(const unsigned lev_out, const unsigned lev_in, const MultiLevelSolution & ml_sol_in);
        
    void CopySolutionToOldSolution();
    
    void SetIfFSI(const bool &FSI = true){
	_FSI = FSI; 
	for(unsigned i=0;i<_gridn;i++){
	  _solution[i]->SetIfFSI(FSI);
	}
    }
      
    bool GetIfFSI(){
      return _FSI; 
    }
    
  
private:
    
    /** boundary condition function pointer */

    BoundaryFunc _SetBoundaryConditionFunction;
    BoundaryFuncMLProb _SetBoundaryConditionFunctionMLProb;
    /** Flag to tell whether the BC function has been set */
    bool _bdcFuncSet;
    bool _bdcFuncSetMLProb;

    /** To be Added */
    BDCType GetBoundaryCondition(const unsigned int var, const unsigned int facename) const;

    /** To be Added */
    bool Ishomogeneous(const unsigned int var, const unsigned int facename) const;

    /** To be Added */
    FunctionBase* GetBdcFunction(const unsigned int var, const unsigned int facename) const;

    /** Array of solution, dimension number of levels */
    vector < Solution* >  _solution;
    unsigned short  _gridn;
    

    /** This group of vectors has the size of the number of added solutions */
    vector < vector <BDCType> >         _boundaryConditions;
    vector < vector <bool> >            _isHomogeneous;
    vector < vector <FunctionBase *> >  _nonHomogeneousBCFunction;
    vector < int >                      _solType;    /* Tells the FE index */
    vector < FEFamily >                 _family;
    vector < FEOrder >                  _order;
    vector < char* >                    _solName;
    vector < char* >                    _bdcType;
    vector < int >                      _solTimeOrder;  //0 = steady, 2 = time-dependent
    vector < bool >                     _pdeType;    /*Tells whether the Solution is an unknown of a PDE or not*/
    vector < bool >                     _testIfPressure;
    vector < bool >                     _addAMRPressureStability;
    vector < bool >                     _fixSolutionAtOnePoint;

    vector <unsigned>                   _solPairIndex;
    vector <unsigned>                   _solPairInverseIndex;

    /** Multilevel solution writer */
    Writer* _writer;

    const MultiLevelProblem* _mlBCProblem;
    bool _FSI;

};


inline
BDCType MultiLevelSolution::GetBoundaryCondition(const unsigned int var, const unsigned int facename) const {
    return _boundaryConditions[var][facename];
}

inline
bool MultiLevelSolution::Ishomogeneous(const unsigned int var, const unsigned int facename) const {
    return _isHomogeneous[var][facename];
}

inline
FunctionBase* MultiLevelSolution::GetBdcFunction(const unsigned int var, const unsigned int facename) const {
    return _nonHomogeneousBCFunction[var][facename];
}

inline
BDCType MultiLevelSolution::GetBoundaryCondition(const std::string varname, const unsigned int facename) const {
    unsigned int var = GetIndex(varname.c_str());
    return _boundaryConditions[var][facename];
}

inline
bool MultiLevelSolution::Ishomogeneous(const std::string varname, const unsigned int facename) const {
    unsigned int var = GetIndex(varname.c_str());
    return _isHomogeneous[var][facename];
}

inline
FunctionBase* MultiLevelSolution::GetBdcFunction(const std::string varname, const unsigned int facename) const {
    unsigned int var = GetIndex(varname.c_str());
    return _nonHomogeneousBCFunction[var][facename];
}

inline 
void MultiLevelSolution::Set(const char * name, InitFuncMLProb funcMLProb, const MultiLevelProblem * ml_prob) {
    Initialize(name, funcMLProb, ml_prob);
}


} //end namespace femus



#endif

