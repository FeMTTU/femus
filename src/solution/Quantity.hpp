/*=========================================================================

 Program: FEMUS
 Module: Quantity
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __quantity_h
#define __quantity_h

#include <string>
#include <map>
#include <iostream>
#include <cstdlib>

#include "Typedefs.hpp"
#include "CurrentQuantity.hpp"


namespace femus {


class SystemTwo;
class Physics;
class QuantityMap;




class Quantity { 
  
public:
  
   Quantity(std::string name_in,QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Quantity();  

  /*virtual*/ void FunctionDof(CurrentQuantity& myvect, const double t,const double* xx) const/* =0*/;
      virtual void Function_txyz(const double t, const double* xp, double* temp) const   = 0;  
      
  void set_eqn(SystemTwo*);
  
 inline       void SetPosInAssocEqn(uint pos_in) { _pos = pos_in; return;}
 
  std::string  _name;      //quantity name, to retrieve it
  uint         _dim;       //number of scalar components
  uint         _FEord;     //FEorder
  double *     _refvalue;  //ref values for the scalar components (_dim)
  QuantityMap& _qtymap;
  SystemTwo *      _eqn;
  uint           _pos;     //block position in the associated equation

  
};



///This class is the GLOBAL QuantityMap of my simulation,
///i.e. it holds the list of all the involved physical quantities



class QuantityMap { 
  
public:

   QuantityMap(const MeshTwo & mesh, const FemusInputParser<double> * map_in);
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
 const MeshTwo & _mesh;
 const FemusInputParser<double> * _physmap;
 
  
};







} //end namespace femus



#endif