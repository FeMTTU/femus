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
#include "SolvertypeEnum.hpp"
#include "MgTypeEnum.hpp"
#include "LinearEquationSolverEnum.hpp"
#include "Mesh.hpp"


namespace femus {

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class MultiLevelProblem;
class MultiLevelMesh;
class Unknown;

/**
 * The system abstract class
 */

class System {



//==== Constructors / Destructor - BEGIN ========
public:
  
    /** Constructor.  Optionally initializes required data structures. */
    System (MultiLevelProblem& ml_prob, const std::string& name, const unsigned int number, const LinearEquationSolverType & smoother_type);

    /** destructor */
    virtual ~System();

  
//==== Constructors / Destructor - END ========
  
  
//==== Basic - BEGIN ========
public:

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
    
    
    /** Init the system PDE structures */
    virtual void init();

protected:
  
    /** The number associated with this system */
    const unsigned int _sys_number;

    /** A name associated with this system. */
    const std::string _sys_name;

    
    
//==== Basic - END ========

    
//==== 
//==== LINEAR ALGEBRA - BEGIN ========
//==== 


//==== Assemble function (Residual and Jacobian) - BEGIN ========

public:
  
    /** Function pointer type, easiest way to declare function pointer instantiations */
    typedef void (* AssembleFunctionType) (MultiLevelProblem &ml_prob);//, unsigned level, const unsigned &gridn, const bool &assemble_matrix);
    
    /** Register a user function to use in assembling the system matrix and RHS. */
    void SetAssembleFunction (AssembleFunctionType );

    AssembleFunctionType  GetAssembleFunction();
    
    
    inline unsigned GetLevelToAssemble() const { return _levelToAssemble; }

    inline void SetLevelToAssemble(const unsigned &level){ _levelToAssemble = level; }

    /** Only call assemble function */
    virtual void assemble_call_before_boundary_conditions(const unsigned int n_times);

protected:
  
    /** Function that assembles the system. */
    AssembleFunctionType _assemble_system_function;

    /** Flag to assemble the matrix or not */
    bool _buildSolver;
    
    unsigned _levelToAssemble;

//==== Assemble function (Residual and Jacobian) - END ========

//==== Unknowns - BEGIN ========
public:
  
    /** Add a vector of solution variables to the system PDE */
    void AddSolutionToSystemPDEVector(const unsigned n_components,  const std::string name);

    /** Associate the solution variables to the system PDE */
    virtual void AddSolutionToSystemPDE(const char solname[]);

    /** Get the index of the Solution "solname" for this system */
    unsigned GetSolPdeIndex(const char solname[]);
    
    /** Get the index of the Solution "solname" for this system */
    const unsigned GetSolPdeIndex(const char solname[]) const;

    std::vector <unsigned> & GetSolPdeIndex() {
      return _SolSystemPdeIndex;
    }

    const std::vector <unsigned> & GetSolPdeIndex() const {
      return _SolSystemPdeIndex;
    }
    
    /** Set Unknown list for the current System */
    void set_unknown_list_for_assembly(const std::vector< Unknown > unknown_in );
    
    /** Get Unknown list for the current System */
    const std::vector< Unknown > get_unknown_list_for_assembly() const;
    
    
protected:

    /** indices of the solutions, dynamical dimension */
    std::vector <unsigned> _SolSystemPdeIndex;

    /** List of unknowns for the assembly routine */
    std::vector< Unknown > _unknown_list_for_assembly;

//==== Unknowns - END ========




//==== 
//==== LINEAR ALGEBRA - END ========
//==== 

    
//==== Solver - BEGIN ========
public:
  
    virtual void SetOuterSolver (const SolverType & mgOuterSolver) {};
    
    virtual void MGsolve (const MgSmootherType& mgSmootherType = MULTIPLICATIVE){
      //solve(mgSmootherType);
    };

// protected:
// virtual void solve( const MgSmootherType& mgSmootherType = MULTIPLICATIVE ){};
    
//==== Solver - END ========

    
//==== Problem - BEGIN ======== 
public:

    /** Get MultiLevelProblem */
    const MultiLevelProblem &  GetMLProb() const { return _equation_systems; }

    /** Get MultiLevelProblem */
    MultiLevelProblem &  GetMLProb() { return _equation_systems; }

protected:

    /** Constant reference to the \p EquationSystems object used for the simulation. */
    MultiLevelProblem & _equation_systems;

//==== Problem - END ======== 


//==== Mesh - BEGIN ========
public:
  
    /** Get Number of Levels */
    inline const unsigned GetGridn() const { return _gridn; }

protected:

    /** Mesh vector, dimension _gridn */
    std::vector<Mesh*> _msh;

    /** pointer */
    MultiLevelMesh* _ml_msh;

    /** Number of Levels */
    unsigned _gridn;

//==== Mesh - END ========
    
    
//==== Solution (may or may not be Unknowns) - BEGIN ========

protected:


    /** Solution vector, dimension _gridn */
    std::vector<Solution*> _solution;

    /** pointer */
    MultiLevelSolution* _ml_sol;

//==== Solution (may or may not be Unknowns) - END ========


};



//==== Basic - BEGIN ========
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
//==== Basic - END ========


} //end namespace femus



#endif
