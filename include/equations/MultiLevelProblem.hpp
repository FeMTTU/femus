/*=========================================================================

 Program: FEMUS
 Module: MultiLevelProblem
 Authors: Eugenio Aulisa, Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __MultiLevelProblem_hpp__
#define __MultiLevelProblem_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelMesh.hpp"
#include "MultiLevelSolution.hpp"
#include "Solution.hpp"
#include "Parameters.hpp"
#include "ParallelObject.hpp"
#include "MgSmootherEnum.hpp"
#include <vector>
#include <map>


namespace femus {

using std::map;

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class System;

typedef double (*initfunc) (const double &x, const double &y, const double &z);

/**
* This class is a black box container to handle multilevel problems.
*/

class MultiLevelProblem {

public:  

  /** Constructor */
  MultiLevelProblem(MultiLevelMesh *ml_msh, MultiLevelSolution *ml_sol);

  /** Destructor */
  ~MultiLevelProblem() {};
  
  /** Multilevel solution pointer */
  MultiLevelSolution *_ml_sol;
  
  /** Multilevel mesh pointer */
  MultiLevelMesh *_ml_msh;
  
  /** Data structure holding arbitrary parameters. */
  Parameters parameters;
  
   /** Data structure holding the systems. */
  std::map<std::string, System*> _systems;

  /** Typedef for system iterators */
  typedef std::map<std::string, System*>::iterator       system_iterator;

  /** Typedef for constatnt system iterators */
  typedef std::map<std::string, System*>::const_iterator const_system_iterator; 
     
  /** Add the system of type \p system_type named \p name to the systems array. */
  virtual System & add_system (const std::string& system_type, const std::string& name);

  /** Add the system named \p name to the systems array. */
  template <typename T_sys> T_sys & add_system (const std::string& name, const MgSmoother & smoother_type = GMRES_SMOOTHER);
  
  /**
   * @returns a constant reference to the system named \p name.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> ("sys");
   * is an example of how the method might be used
   */
  template <typename T_sys>
  const T_sys & get_system (const std::string& name) const;

  /**
   * @returns a writeable referene to the system named \p name.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> ("sys");
   * is an example of how the method might be used
   */
  template <typename T_sys>
  T_sys & get_system (const std::string& name);

  /**
   * @returns a constant reference to system number \p num.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> (0);
   * is an example of how the method might be used
   */
  template <typename T_sys>
  const T_sys & get_system (const unsigned int num) const;

  /**
   * @returns a writeable referene to the system number \p num.
   * The template argument defines the return type.  For example,
   * const SteadySystem& sys = eq.get_system<SteadySystem> (0);
   * is an example of how the method might be used
   */
  template <typename T_sys>
  T_sys & get_system (const unsigned int num);

  /**
   * @returns a constant reference to the system named \p name.
   */
  const System & get_system (const std::string& name) const;

  /**
   * @returns a writeable referene to the system named \p name.
   */
  System & get_system (const std::string& name);

  /**
   * @returns a constant reference to system number \p num.
   */
  const System & get_system (const unsigned int num) const;

  /**
   * @returns a writeable referene to the system number \p num.
   */
  System & get_system (const unsigned int num);
    
  /** @returns the number of equation systems. */
  unsigned int n_systems() const;
  
  /** Clear all the Sytems PDE structures */
  void clear();
  
  /** init the system pde structures */
  void init();
  
  /** Get the total number of grid, both totally refined and partial refined for DomainDecomposition */
  const unsigned GetNumberOfGrid() const {return _gridn; };
  
  /** Get the number of grid totally refined; it is a subset of _gridn */
  const unsigned GetNumberOfGridTotallyRefined() const {return _gridr; };

  /** This function is not supposed to be placed here!!! */
  int ComputeBdIntegral(const char pdename[],const char var_name[], const unsigned & kel, 
                         const unsigned & jface, unsigned level, unsigned dir);
  

private:
    
  vector < map <unsigned,bool> > index;
  
  unsigned short _gridn;
  
  unsigned short _gridr;
  
};

template <typename T_sys>
inline
T_sys & MultiLevelProblem::add_system (const std::string& name,const MgSmoother & smoother_type )
{
  T_sys* ptr = NULL;

  if (!_systems.count(name))
    {
      ptr = new T_sys(*this, name, this->n_systems(),smoother_type);

      _systems.insert (std::make_pair(name, ptr));

    }
  else
    {
      // We now allow redundant add_system calls, to make it
      // easier to load data from files for user-derived system
      // subclasses
      std::cerr << "ERROR: There was already a system"
              << " named " << name
              << std::endl;

      ptr = &(this->get_system<T_sys>(name));
    }

  // Return a dynamically casted reference to the newly added System.
  return *ptr;
}


template <typename T_sys>
inline
const T_sys & MultiLevelProblem::get_system (const std::string& name) const
{
  const_system_iterator pos = _systems.find(name);

  // Check for errors
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: no system named \"" << name << "\" found!"
                    << std::endl;
    }

  // Attempt dynamic cast
  return *static_cast<T_sys*>(pos->second);
}

template <typename T_sys>
inline
T_sys & MultiLevelProblem::get_system (const std::string& name)
{
  system_iterator pos = _systems.find(name);

  // Check for errors
  if (pos == _systems.end())
    {
      std::cerr << "ERROR: no system named " << name << " found!"
                    << std::endl;
    }

  // Attempt dynamic cast
  return *static_cast<T_sys*>(pos->second);
}

inline
const System & MultiLevelProblem::get_system (const std::string& name) const
{
  return this->get_system<System>(name);
}


inline
System & MultiLevelProblem::get_system (const std::string& name)
{
  return this->get_system<System>(name);
}


inline
const System & MultiLevelProblem::get_system (const unsigned int num) const
{
  return this->get_system<System>(num);
}


inline
System & MultiLevelProblem::get_system (const unsigned int num)
{
  return this->get_system<System>(num);
}

inline
unsigned int MultiLevelProblem::n_systems () const
{
  return static_cast<unsigned int>(_systems.size());  //libmesh static_cast_int
}


} //end namespace femus



#endif

