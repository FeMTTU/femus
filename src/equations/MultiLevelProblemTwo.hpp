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

#include "MultiLevelProblem.hpp"

namespace femus {



class MultiLevelMeshTwo;
class elem_type;
class QuantityMap;


class MultiLevelProblemTwo : public MultiLevelProblem {

public:
  

  /// Constructor
    MultiLevelProblemTwo(const MultiLevelMeshTwo& mesh_in,
		  const std::string quadr_order_in
		);

  inline const  MultiLevelMeshTwo & GetMeshTwo() const { return  _mesh; }
  
  inline void SetQtyMap(const QuantityMap * qtymap_in) { _qtymap = qtymap_in; return; }
   
  inline const QuantityMap & GetQtyMap() const { return  *_qtymap; }

  inline const std::vector<elem_type*>  & GetElemType(const unsigned dim) const { return  _elem_type[dim - 1]; }
    
  inline const std::vector< std::vector<elem_type*> >  & GetElemType() const { return  _elem_type; }

  inline const Gauss & GetQrule(const unsigned dim) const { return _qrule[dim - 1]; }
  
  void SetQruleOrder(const std::string order_in) { _quadrature_order = order_in; return; }
  
  inline const FemusInputParser<double> &  GetInputParser() const { return *_phys; }

  void SetInputParser(const FemusInputParser<double> * parser_in) { _phys = parser_in; return; }
  
  ~MultiLevelProblemTwo(){};
  
  void clean();

  // equation get/set
  inline          void  add_system(SystemTwo* value)            {_equations.insert(make_pair(value->_eqname,value));}
  inline       SystemTwo* get_system(const string & name)       {return _equations.find(name)->second;}
  inline const SystemTwo* get_system(const string & name) const {return _equations.find(name)->second;}

  typedef std::map<string, SystemTwo*>::iterator iterator;
  typedef std::map<string, SystemTwo*>::const_iterator const_iterator;

  inline iterator       begin()       { return _equations.begin();}
  inline iterator         end()       { return _equations.end();}
  inline const_iterator begin() const { return _equations.begin();}
  inline const_iterator   end() const { return _equations.end();}

private:
  
    map<string,SystemTwo*> _equations;   // system map
    
    std::vector< std::vector<elem_type*> >  _elem_type;
    
    std::vector<Gauss>       _qrule;
    
    std::string _quadrature_order;
    
    const FemusInputParser<double>  * _phys;

    const QuantityMap * _qtymap;
    
    const MultiLevelMeshTwo &     _mesh;

};


} //end namespace femus



#endif