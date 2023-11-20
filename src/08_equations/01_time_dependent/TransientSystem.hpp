/*=========================================================================

 Program: FEMUS
 Module: TransientSystem
 Authors: Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_TransientSystem_hpp__
#define __femus_equations_TransientSystem_hpp__

#include "LinearEquationSolverEnum.hpp"
#include "MgTypeEnum.hpp"

#include <string>
#include <vector>

namespace femus {


//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class System;
class LinearImplicitSystem;
class NonLinearImplicitSystem;


class MultiLevelProblem;


/**
 * This class provides a specific system class.  It aims
 * at transient systems, offering nothing more than just
 * the essentials needed to solve a system.
 */

// ------------------------------------------------------------
// TransientSystem class definition
template <class Base>
class TransientSystem : public Base {


//==== Constructors / Destructor - BEGIN ========
public:
    /** Constructor.  Initializes required data structures.  */
    TransientSystem (MultiLevelProblem& ml_probl,
                     const std::string& name,
                     const unsigned int number, 
                     const LinearEquationSolverType & smoother_type);

    /** Destructor. */
    virtual ~TransientSystem ();
//==== Constructors / Destructor - END ========


//==== Basic - BEGIN ========
public:
    
    /** The type of system. */
    typedef TransientSystem<Base> sys_type;

    /** @returns a clever pointer to the system. */
    sys_type & system () {
        return *this;
    }

    /**
     * @returns \p "Transient" prepended to T::system_type().
     * Helps in identifying the system type in an equation
     * system file.
    */
    virtual std::string system_type () const;
//==== Basic - END ========


//==== Time - BEGIN ========
public:
    
    /** Set up before calling the parent solve */
    void SetUpForSolve();
    
    /** attach the GetTimeInterval Function for selective interval time */
    void AttachGetTimeIntervalFunction (double (* get_time_interval_function)(const double time));


    /** Set the interval time */
    void SetIntervalTime(const double dt) {
        _dt = dt;
    };


    /** Get the interval time */
    double GetIntervalTime() const {
        return _dt;
    };

    /** Get the time */
    double GetTime() const {
        return _time;
    };

    /** Get the time */
    void SetTime(const double time) {
        _time = time;
    };

protected:

    double _dt;
    
    double _time;

private:

    bool _is_selective_timestep;

    unsigned int _time_step;

    /** pointer function to the set time step function */
    double (* _get_time_interval_function)(const double time);

    /** Count how many times the assemble is performed */
    unsigned _assembleCounter;
    
//==== Time - END ========

   
//==== Time, Solver of 1 iteration - BEGIN ========
public:
    
    /** calling the parent solve */
    virtual void MGsolve( const MgSmootherType& mgSmootherType = MULTIPLICATIVE );
    
//==== Time, Solver of 1 iteration - END ========


//==== Time, Solution Update - BEGIN ========
public:
    
    /** Update the old solution with new ones. It calls the update solution function of the Solution class */
    void CopySolutionToOldSolution();
//==== Time, Solution Update - END ========

    
//==== Time, Newmark - BEGIN ========
public:
    
    /** update the Newmark variables */
    void NewmarkAccUpdate(const std::vector< std::string > acceleration_name,
                          const std::vector< std::string > velocity_name);

//==== Time, Newmark - END ========


};


// -----------------------------------------------------------
// Useful typedefs
typedef TransientSystem<System> TransientBaseSystem;
typedef TransientSystem<LinearImplicitSystem> TransientLinearImplicitSystem;
typedef TransientSystem<NonLinearImplicitSystem> TransientNonlinearImplicitSystem;




// ------------------------------------------------------------
// TransientSystem inline methods
template <class Base>
inline
std::string TransientSystem<Base>::system_type () const
{
    std::string type = "Transient";
    type += Base::system_type ();

    return type;
}


} //end namespace femus



#endif
