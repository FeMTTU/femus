#ifndef __transient_system_h__
#define __transient_system_h__

// Local Includes

// C++ includes

// Forward declarations
class LinearImplicitSystem;
class NonLinearImplicitSystem;
class ExplicitSystem;
class NonLinearMultiLevelProblem;

/**
 * This class provides a specific system class.  It aims
 * at transient systems, offering nothing more than just
 * the essentials needed to solve a system.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the parent classes.
 */

// ------------------------------------------------------------
// TransientSystem class definition
template <class Base>
class TransientSystem : public Base
{
public:

  /**
   * Constructor.  Initializes required
   * data structures.
   */
  TransientSystem (MultiLevelProblem& ml_probl,
		   const std::string& name,
		   const unsigned int number);

  /**
   * Destructor.
   */
  virtual ~TransientSystem ();

  /**
   * The type of system.
   */
  typedef TransientSystem<Base> sys_type;

  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear ();

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit ();

  /**
   * @returns \p "Transient" prepended to T::system_type().
   * Helps in identifying the system type in an equation
   * system file.
   */
  virtual std::string system_type () const;


  //-----------------------------------------------------------------
  // access to the solution data fields

  /**
   * @returns the old solution (at the previous timestep)
   * for the specified global DOF.
   */
  double old_solution (const int global_dof_number) const;

  /**
   * @returns the older solution (two timesteps ago)
   * for the specified global DOF.
   */
  double older_solution (const int global_dof_number) const;

  /**
   * All the values I need to compute my contribution
   * to the simulation at hand.  Think of this as the
   * current solution with any ghost values needed from
   * other processors.
   */
  NumericVector* old_local_solution;

  /**
   * All the values I need to compute my contribution
   * to the simulation at hand.  Think of this as the
   * current solution with any ghost values needed from
   * other processors.
   */
  NumericVector* older_local_solution;


protected:


  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void init_data ();

  /**
   * Re-update the local values when the mesh has changed.
   * This method takes the data updated by \p update() and
   * makes it up-to-date on the current mesh.
   */
  virtual void re_update ();
};



// -----------------------------------------------------------
// Useful typedefs
 typedef TransientSystem<LinearImplicitSystem> TransientImplicitSystem;
 typedef TransientSystem<LinearImplicitSystem> TransientLinearImplicitSystem;
 typedef TransientSystem<NonLinearImplicitSystem> TransientNonlinearImplicitSystem;
 typedef TransientSystem<ExplicitSystem> TransientExplicitSystem;
 typedef TransientSystem<System> TransientBaseSystem;



// ------------------------------------------------------------
// TransientSystem inline methods
template <class Base>
inline
std::string TransientSystem<Base>::system_type () const
{
  std::string type = "Transient";
  //type += Base::system_type ();  // ancora non esiste

  return type;
}



#endif 
