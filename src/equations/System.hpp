/*=========================================================================

 Program: FEMUS
 Module: System
 Authors: Simone Bnà, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_System_hpp__
#define __femus_equations_System_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelProblem.hpp"


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class System;
class MultiLevelProblem;
class String;

/**
 * The system abstract class
 */

class System {

protected:

    /** Function pointer type, easiest way to declare function pointer instantiations */
    typedef void (* AssembleFunctionType) (MultiLevelProblem &ml_prob);//, unsigned level, const unsigned &gridn, const bool &assemble_matrix);

public:

    /** Constructor.  Optionally initializes required data structures. */
    System (MultiLevelProblem& ml_prob, const std::string& name, const unsigned int number, const MgSmoother & smoother_type);

    /** destructor */
    virtual ~System();

    /** To be Added */
    unsigned int number() const;

    /** To be Added */
    const std::string & name() const;

    /**
     * @returns the type of system, helpful in identifying
     * which system type to use when reading equation system
     * data from file.  Should be overloaded in derived classes.
    */
    virtual std::string system_type () const {
        return "Basic";
    }

    /** Add a vector of solution variables to the system PDE */
    void AddSolutionToSystemPDEVector(const unsigned n_components,  const std::string name);

    /** Associate the solution variables to the system PDE */
    void AddSolutionToSystemPDE(const char solname[]);

    /** Register a user function to use in assembling the system matrix and RHS. */
    void SetAssembleFunction (AssembleFunctionType );

    AssembleFunctionType  GetAssembleFunction();

    /** Solves the system.  Should be overloaded in derived systems. */
    virtual void solve () {};

    /** Clear all the data structures associated with the system. */
    virtual void clear();

    /** Init the system PDE structures */
    virtual void init();

    /** @deprecated Init the system PDE structures */
    virtual void init_two(){};

    /** Get the index of the Solution "solname" for this system */
    unsigned GetSolPdeIndex(const char solname[]);

    /** Get MultiLevelProblem */
    const MultiLevelProblem &  GetMLProb() const { return _equation_systems; }

    /** Get MultiLevelProblem */
    MultiLevelProblem &  GetMLProb() { return _equation_systems; }

    /** Get Number of Levels */
    inline const unsigned GetGridn() const { return _gridn; }

    inline unsigned GetLevelToAssemble(){ return _levelToAssemble; }
    inline unsigned GetLevelMax(){ return _levelMax; }
    inline bool GetAssembleMatrix(){ return _assembleMatrix; }

    inline unsigned SetLevelToAssemble(const unsigned &level){ _levelToAssemble = level; }
    inline bool SetAssembleMatrix(const bool& assembleMatrix){ _assembleMatrix = assembleMatrix; }
    inline unsigned SetLevelMax(const unsigned &levelMax){ _levelMax = levelMax; }

protected:

    /** Constant reference to the \p EquationSystems object used for the simulation. */
    MultiLevelProblem& _equation_systems;

    /** Mesh vector, dimension _gridn */
    vector<Mesh*> _msh;

    /** Solution vector, dimension _gridn */
    vector<Solution*> _solution;

    /** pointer */
    MultiLevelSolution* _ml_sol;

    /** pointer */
    MultiLevelMesh* _ml_msh;

    /** indices of the solutions, dynamical dimension */
    vector <unsigned> _SolSystemPdeIndex;

    /** Number of Levels */
    unsigned _gridn;

    /** Number of Totally Refined Levels */
    unsigned _gridr;

    bool _assembleMatrix;
    unsigned _levelToAssemble;
    unsigned _levelMax;





    /** Function that assembles the system. */
    AssembleFunctionType _assemble_system_function;

    /** The number associated with this system */
    const unsigned int _sys_number;

    /** A name associated with this system. */
    const std::string _sys_name;

};

// System inline methods
inline
const std::string & System::name() const
{
    return _sys_name;
}

inline
unsigned int System::number() const
{
    return _sys_number;
}


} //end namespace femus



#endif
