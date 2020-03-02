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

#ifndef __femus_equations_MultiLevelProblem_hpp__
#define __femus_equations_MultiLevelProblem_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <vector>
#include <map>
#include "MultiLevelMesh.hpp"
#include "Parameters.hpp"
#include "ParallelObject.hpp"
#include "LinearEquationSolverEnum.hpp"
#include "GaussPoints.hpp"
#include "FemusInputParser.hpp"
#include "Files.hpp"
#include "System.hpp"
#include "ElemType_template.hpp"

namespace femus {

using std::map;

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class MultiLevelSolution;
class MultiLevelMeshTwo;
class elem_type;
class QuantityMap;



/**
* This class is a black box container to handle multilevel problems.
*/

class MultiLevelProblem {

public:

    /** Constructor */
    MultiLevelProblem();

    MultiLevelProblem(MultiLevelSolution *ml_sol);

    /** Destructor */

    ~MultiLevelProblem();

    /** Multilevel solution pointer */
    MultiLevelSolution *_ml_sol;

    /** Multilevel mesh pointer */
    MultiLevelMesh *_ml_msh;

    /** Data structure holding arbitrary parameters. */
    Parameters parameters;

    /** Typedef for system iterators */
    typedef std::map<std::string, System*>::iterator       system_iterator;

    /** Typedef for constatnt system iterators */
    typedef std::map<std::string, System*>::const_iterator const_system_iterator;

    /** Add the system of type \p system_type named \p name to the systems array. */
    virtual System & add_system (const std::string& system_type, const std::string& name);

    /** Add the system named \p name to the systems array. */
    template <typename T_sys> T_sys & add_system (const std::string& name, const LinearEquationSolverType & smoother_type = FEMuS_DEFAULT);

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

//   /** init the system pde structures */
//   void init();

    /** Get the total number of grid, both totally refined and partial refined for DomainDecomposition */
    const unsigned GetNumberOfLevels() const {
        return _gridn;
    };

    /** Increase of one the number of levels */
    void AddLevel(){
        _gridn++;
    };

    /** returns the beginning of the systems map */
  inline system_iterator       begin()       { return _systems.begin(); }

    /**  */
  inline system_iterator         end()       { return _systems.end(); }

    /**  */
  inline const_system_iterator begin() const { return _systems.begin(); }

    /**  */
  inline const_system_iterator   end() const { return _systems.end(); }

    /**  New get/set for new data */
  inline void SetMeshTwo(const MultiLevelMeshTwo * mesh_in)  {   _mesh = mesh_in; return; }

  inline const  MultiLevelMeshTwo & GetMeshTwo() const { return  *_mesh; }

    /** Quantity Map */
  inline void SetQtyMap(const QuantityMap * qtymap_in) { _qtymap = qtymap_in; return; }

  inline const QuantityMap & GetQtyMap() const { return  *_qtymap; }

    /** ElemType and Quadrature rule */
  inline const std::vector<const elem_type*>  & GetElemType(const unsigned dim) const { return  _elem_type[dim - 1]; }

  inline const std::vector< std::vector<const elem_type*> >  & GetElemType() const { return  _elem_type; }

  inline const std::vector<Gauss> & GetQuadratureRuleAllGeomElems() const { return _qrule; }
  
  inline const Gauss & GetQuadratureRule(const unsigned geom_elem_type) const { return _qrule[geom_elem_type]; }

  void SetQuadratureRuleAllGeomElems(const std::string quadr_order_in);
  
    /** Input Parser */
  inline const FemusInputParser<double> &  GetInputParser() const { return *_phys; }

  void SetInputParser(const FemusInputParser<double> * parser_in) { _phys = parser_in; return; }

    /** Files Handler */
  void SetFilesHandler(const Files * files_in) { _files = files_in; return; }
  
  inline const Files * GetFilesHandler() const { return  _files; }

    /** Input Parser */
  const int get_current_system_number() const { return _current_system_number; }

  void set_current_system_number(const unsigned current_system_number_in) { _current_system_number = current_system_number_in; }

  void SetMultiLevelMeshAndSolution(MultiLevelMesh * ml_mesh, MultiLevelSolution * ml_sol);
  
  //   Returns a non-const reference to the map of Systems
  std::map<std::string, System*> & get_systems_map() { return _systems; }
  
  
  void get_all_abstract_fe(std::vector < std::vector < /*const*/ elem_type_templ_base< double, double > *  > > & elem_all_in)                 /*const*/ { elem_all_in = _elem_all_dd; }
  
  void get_all_abstract_fe(std::vector < std::vector < /*const*/ elem_type_templ_base< adept::adouble, double > *  > > & elem_all_in)         /*const*/ { elem_all_in = _elem_all_ad; }
  
  void get_all_abstract_fe(std::vector < std::vector < /*const*/ elem_type_templ_base< adept::adouble, adept::adouble > *  > > & elem_all_in) /*const*/ { elem_all_in = _elem_all_aa; }

  void get_all_abstract_fe(std::vector < std::vector < /*const*/ elem_type_templ_base< double, double > *  > > & elem_all_in)                 const { elem_all_in = _elem_all_dd; }
  
  void get_all_abstract_fe(std::vector < std::vector < /*const*/ elem_type_templ_base< adept::adouble, double > *  > > & elem_all_in)         const { elem_all_in = _elem_all_ad; }
  
  void get_all_abstract_fe(std::vector < std::vector < /*const*/ elem_type_templ_base< adept::adouble, adept::adouble > *  > > & elem_all_in) const { elem_all_in = _elem_all_aa; }

  
  
  void set_all_abstract_fe() {
      
       set_all_abstract_fe<double, double>(_elem_all_dd);
       set_all_abstract_fe<adept::adouble, double>(_elem_all_ad);
       set_all_abstract_fe<adept::adouble, adept::adouble>(_elem_all_aa);
    
}  
  
 template <class type, class type_mov>
  void set_all_abstract_fe(std::vector < std::vector < /*const*/ elem_type_templ_base<type, type_mov> *  > > & elem_all_in) const {

//this function performs the initialization of all abstract FE families on all abstract Geometric Elements      
      
//     clock_t start_evals = clock();
  
  //prepare Abstract quantities for all fe fams for all geom elems: perform all quadrature evaluations beforehand
        elem_all_in.resize( femus::geom_elems.size() );
  
         for (unsigned int g = 0; g < femus::geom_elems.size(); g++) {
             elem_all_in[g].resize(femus::fe_fams.size());
             const std::string quad_order = this->GetQuadratureRule(g).GetGaussOrderString();  ///@todo what if you choose different quadrature orders on different geom elems?

         for (unsigned int fe = 0; fe < femus::fe_fams.size(); fe++) {
            elem_all_in[g][fe] = elem_type_templ_base<type, type_mov>::build(femus::geom_elems[g], femus::fe_fams[fe], quad_order.c_str(), 3);          
           }
       }
       
//   clock_t end_evals = clock();
//    std::cout << " FE Evals time " << static_cast<double>(end_evals - start_evals) / CLOCKS_PER_SEC << std::endl;

   
   }
   
   
private:


    /** Data structure holding the systems. */
    std::map<std::string, System*> _systems;

    // member data
    vector < map <unsigned,bool> > index;
    unsigned short _gridn;

    std::vector< std::vector<const elem_type*> >  _elem_type;  ///@deprecated 
    std::vector<Gauss>                      _qrule;            //over all Geom Elems
    const FemusInputParser<double>        * _phys;
    const QuantityMap                     * _qtymap;
    const MultiLevelMeshTwo               * _mesh;

    const Files                           * _files;
    unsigned int _current_system_number;

    // attempt to handle templated classes from non-templated class
    std::vector< std::vector< /*const*/ elem_type_templ_base< double, double > * > >  _elem_all_dd;
    std::vector< std::vector< /*const*/ elem_type_templ_base< adept::adouble, double > * > >  _elem_all_ad;
    std::vector< std::vector< /*const*/ elem_type_templ_base< adept::adouble, adept::adouble > * > >  _elem_all_aa;
    
    
    
};


template <typename T_sys>
inline
const T_sys & MultiLevelProblem::get_system (const unsigned int num) const
{
  assert(num < this->n_systems());

  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    if (pos->second->number() == num)
      break;

  // Check for errors
  if (pos == end)
  {
    std::cerr << "ERROR: no system number " << num << " found!" << std::endl;
  }

  // Attempt dynamic cast
  return *static_cast<T_sys*>(pos->second);
}

template <typename T_sys>
inline
T_sys & MultiLevelProblem::get_system (const unsigned int num)
{
  assert(num < this->n_systems());

  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    if (pos->second->number() == num)
      break;

  // Check for errors
  if (pos == end)
  {
    std::cerr << "ERROR: no system number " << num << " found!" << std::endl;
  }

  // Attempt dynamic cast
  return *static_cast<T_sys*>(pos->second);
}



template <typename T_sys>
inline
T_sys & MultiLevelProblem::add_system (const std::string& name,const LinearEquationSolverType & smoother_type )
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

