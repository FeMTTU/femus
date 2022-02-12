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

#include <string>

#include "LinearEquationSolverEnum.hpp"
#include "MgTypeEnum.hpp"

namespace femus {


//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class MonolithicFSINonLinearImplicitSystem;
class LinearImplicitSystem;
class NonLinearImplicitSystem;
class ExplicitSystem;
class MultiLevelProblem;
class System;


/**
 * This class provides a specific system class.  It aims
 * at transient systems, offering nothing more than just
 * the essentials needed to solve a system.
 */

// ------------------------------------------------------------
// TransientSystem class definition
template <class Base>
class TransientSystem : public Base {

public:

    /** Constructor.  Initializes required data structures.  */
    TransientSystem (MultiLevelProblem& ml_probl,
                     const std::string& name,
                     const unsigned int number, 
                     const LinearEquationSolverType & smoother_type);

    /** Destructor. */
    virtual ~TransientSystem ();

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


    /** Update the old solution with new ones. It calls the update solution function of the Solution class */
    virtual void CopySolutionToOldSolution();

    /** Set up before calling the parent solve */
    void SetUpForSolve();
   
    /** calling the parent solve */
    virtual void MGsolve( const MgSmootherType& mgSmootherType = MULTIPLICATIVE );

    /** update the Newmark variables */
    void NewmarkAccUpdate();


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

    unsigned _assembleCounter;

};


// -----------------------------------------------------------
// Useful typedefs
typedef TransientSystem<LinearImplicitSystem> TransientImplicitSystem;
typedef TransientSystem<LinearImplicitSystem> TransientLinearImplicitSystem;
typedef TransientSystem<NonLinearImplicitSystem> TransientNonlinearImplicitSystem;
typedef TransientSystem<MonolithicFSINonLinearImplicitSystem> TransientMonolithicFSINonlinearImplicitSystem;
typedef TransientSystem<ExplicitSystem> TransientExplicitSystem;
typedef TransientSystem<System> TransientBaseSystem;


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
