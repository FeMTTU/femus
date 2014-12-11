#ifndef __quantity_h
#define __quantity_h

#include <string>
#include <map>
#include <iostream>
#include <cstdlib>

#include "Typedefs.hpp"
#include "QuantityLocal.hpp"


namespace femus {


class EqnBase;
class Physics;
class QuantityMap;




class Quantity { 
  
public:
  
   Quantity(std::string name_in,QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Quantity();  

  /*virtual*/ void FunctionDof(const uint vb, QuantityLocal& myvect, const double t,const double* xx) const/* =0*/;
      virtual void Function_txyz(const double t, const double* xp, double* temp) const   = 0;  
      
  void set_eqn(EqnBase*);
  
 inline       void SetPosInAssocEqn(uint pos_in) { _pos = pos_in; return;}
 
  std::string  _name;      //quantity name, to retrieve it
  uint         _dim;       //number of scalar components
  uint         _FEord;     //FEorder
  double *     _refvalue;  //ref values for the scalar components (_dim)
  QuantityMap& _qtymap;
  EqnBase *      _eqn;
  uint           _pos;     //block position in the associated equation

  
};



///This class is the GLOBAL QuantityMap of my simulation,
///i.e. it holds the list of all the involved physical quantities



class QuantityMap { 
  
public:

   QuantityMap(Physics& phys_in);
  ~QuantityMap(){};
  
  inline           void  set_qty(Quantity* value)          {_QuantMap.insert(make_pair(value->_name,value));}
  
  inline       Quantity* get_qty(const std::string & name)      {
 
    std::map<std::string,Quantity*>::iterator myit = _QuantMap.find(name);
      if ( myit == _QuantMap.end() ) { 
       std::cout << "QuantityMap::get_qty: Sorry but there is no ---> "
                 << name  << " <--- element in the Global Quantity Map" << std::endl; abort();}

    return myit->second;
    
  }

  
 std::map<std::string,Quantity*> _QuantMap;
 Physics& _phys; 
  
};







} //end namespace femus



#endif