// C++ includes

// Local includes
#include "MultiLevelProblem.hpp"
#include "TransientSystem.hpp"
#include "ExplicitSystem.hpp"
#include "LinearImplicitSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

// ------------------------------------------------------------
// TransientSystem implementation
template <class Base>
TransientSystem<Base>::TransientSystem (MultiLevelProblem& ml_probl,
					const std::string& name_in,
					const unsigned int number_in) :

  Base                 (ml_probl, name_in, number_in)
{
// #ifdef LIBMESH_ENABLE_GHOSTED
//   old_local_solution =
//     AutoPtr<NumericVector<Number> >
//       (&(this->add_vector("_transient_old_local_solution", true, GHOSTED)));
//   older_local_solution =
//     AutoPtr<NumericVector<Number> >
//       (&(this->add_vector("_transient_older_local_solution", true, GHOSTED)));
// #else
//   old_local_solution =
//     AutoPtr<NumericVector<Number> >
//       (&(this->add_vector("_transient_old_local_solution", true, SERIAL)));
//   older_local_solution =
//     AutoPtr<NumericVector<Number> >
//       (&(this->add_vector("_transient_older_local_solution", true, SERIAL)));
// #endif
}



template <class Base>
TransientSystem<Base>::~TransientSystem ()
{
  this->clear();

  // We still have AutoPtrs for API compatibility, but
  // now that we're System::add_vector()ing these, we can trust
  // the base class to handle memory management
  old_local_solution;
  older_local_solution;
}



template <class Base>
void TransientSystem<Base>::clear ()
{
  // clear the parent data
  //Base::clear(); // ancora non esiste

  // the old & older local solutions
  // are now deleted by System!
  // old_local_solution->clear();
  // older_local_solution->clear();

//   // FIXME: This preserves maximum backwards compatibility,
//   // but is probably grossly unnecessary:
//   old_local_solution.release();
//   older_local_solution.release();
// 
//   old_local_solution =
//     AutoPtr<NumericVector<Number> >
//       (&(this->add_vector("_transient_old_local_solution")));
//   older_local_solution =
//     AutoPtr<NumericVector<Number> >
//       (&(this->add_vector("_transient_older_local_solution")));
}




template <class Base>
void TransientSystem<Base>::init_data ()
{
  // initialize parent data
  //Base::init_data(); ancora non esiste

//   // Initialize the old & older solutions
//   // Using new ghosted vectors if enabled
// #ifdef LIBMESH_ENABLE_GHOSTED
//   old_local_solution->init   (this->n_dofs(), this->n_local_dofs(),
//                               this->get_dof_map().get_send_list(), false,
//                               GHOSTED);
//   older_local_solution->init (this->n_dofs(), this->n_local_dofs(),
//                               this->get_dof_map().get_send_list(), false,
//                               GHOSTED);
// #else
//   old_local_solution->init   (this->n_dofs(), false, SERIAL);
//   older_local_solution->init (this->n_dofs(), false, SERIAL);
// #endif
}



template <class Base>
void TransientSystem<Base>::reinit ()
{
  // initialize parent data
  //Base::reinit();  // ancora non esiste

  // Project the old & older vectors to the new mesh
  // The System::reinit handles this now
  // this->project_vector (*old_local_solution);
  // this->project_vector (*older_local_solution);
}



template <class Base>
void TransientSystem<Base>::re_update ()
{
  // re_update the parent system
//   Base::re_update ();
// 
//   const std::vector<dof_id_type>& send_list = this->get_dof_map().get_send_list ();
// 
//   const dof_id_type first_local_dof = Base::get_dof_map().first_dof();
//   const dof_id_type end_local_dof  = Base::get_dof_map().end_dof();
// 
//   // Check sizes
//   libmesh_assert_greater_equal (end_local_dof, first_local_dof);
//   libmesh_assert_greater_equal (older_local_solution->size(), send_list.size());
//   libmesh_assert_greater_equal (old_local_solution->size(), send_list.size());
// 
//   // Even if we don't have to do anything ourselves, localize() may
//   // use parallel_only tools
//   // if (first_local_dof == end_local_dof)
//   //   return;
// 
//   // Update the old & older solutions with the send_list,
//   // which may have changed since their last update.
//   older_local_solution->localize (first_local_dof,
// 				  end_local_dof-1,
// 				  send_list);
// 
//   old_local_solution->localize (first_local_dof,
// 				end_local_dof-1,
// 				send_list);
}




template <class Base>
double TransientSystem<Base>::old_solution (const int global_dof_number) const
{
  // Check the sizes
//   libmesh_assert_less (global_dof_number, this->get_dof_map().n_dofs());
//   libmesh_assert_less (global_dof_number, old_local_solution->size());
// 
//   return (*old_local_solution)(global_dof_number);
}



template <class Base>
double TransientSystem<Base>::older_solution (const int global_dof_number) const
{
  // Check the sizes
//   libmesh_assert_less (global_dof_number, this->get_dof_map().n_dofs());
//   libmesh_assert_less (global_dof_number, older_local_solution->size());
// 
//   return (*older_local_solution)(global_dof_number);
}




// ------------------------------------------------------------
// TransientSystem instantiations
template class TransientSystem<LinearImplicitSystem>;
template class TransientSystem<NonLinearImplicitSystem>;
template class TransientSystem<ExplicitSystem>;
template class TransientSystem<System>;

