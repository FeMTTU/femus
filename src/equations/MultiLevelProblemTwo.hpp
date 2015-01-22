/*=========================================================================

 Program: FEMUS
 Module: MultiLevelProblemTwo
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __mgequationsmap_h__
#define __mgequationsmap_h__

#include <map>
#include <string>
using namespace std;

#include "Typedefs.hpp"


#include "SystemTwo.hpp"
#include "GaussPoints.hpp"


namespace femus {



class Physics;
class MeshTwo;
class FEElemBase;
class elem_type;
class TimeLoop;

class QuantityMap;


class MultiLevelProblemTwo  {

public:
  
    Files&       _files;
    FemusInputParser<double> &  _phys;
    QuantityMap& _qtymap;
    MeshTwo&     _mesh;
    std::vector< std::vector<elem_type*> >  &  _elem_type;
    std::vector<Gauss>       _qrule;

  /// Constructor
    MultiLevelProblemTwo( Files& files_in,
		  FemusInputParser<double> & phys_in,
		  QuantityMap& qtymap_in,
		  MeshTwo& mesh_in,
                  std::vector< std::vector<elem_type*> > & elem_type_in,
		  std::vector<Gauss> qrule_in
		);

    
    /// Destructor
  ~MultiLevelProblemTwo(){};
  
  void clean(); ///< Clean all substructures

  // equation get/set
  inline          void  set_eqs(SystemTwo* value)            {_equations.insert(make_pair(value->_eqname,value));}
  inline       SystemTwo* get_eqs(const string & name)       {return _equations.find(name)->second;}
  inline const SystemTwo* get_eqs(const string & name) const {return _equations.find(name)->second;}

  typedef std::map<string, SystemTwo*>::iterator iterator;
  typedef std::map<string, SystemTwo*>::const_iterator const_iterator;

  inline iterator       begin()       { return _equations.begin();}
  inline iterator         end()       { return _equations.end();}
  inline const_iterator begin() const { return _equations.begin();}
  inline const_iterator   end() const { return _equations.end();}


  void setDofBcOpIc() ;


  // read-print functions -------------------------------------------
  void PrintSol(const uint t_step,const double curr_time) const; ///< Print solution 
  void ReadSol(const uint t_step,double& time_out) const;  ///< Read solution //TODO must be updated
  void PrintCase(const uint t_init) const; ///< Print ic and bc
  
private:
  
 map<string,SystemTwo*> _equations;   // system map
    
 void PrintSolXDMF(const uint t_step,const double curr_time) const;
 void PrintSolHDF5(const uint t_flag) const;
 
 void PrintCaseXDMF(const uint t_init) const;
 void PrintCaseHDF5(const uint t_init) const;

 void PrintXDMFTopologyGeometry(std::ofstream& out,const uint Level, const uint vb) const;

 
};


} //end namespace femus



#endif
