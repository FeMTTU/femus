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

#ifndef __femus_solution_Quantity_hpp__
#define __femus_solution_Quantity_hpp__

#include <string>
#include <map>
#include <iostream>
#include <cstdlib>

#include "Typedefs.hpp"
#include "CurrentQuantity.hpp"


namespace femus {


class SystemTwo;
class QuantityMap;


//a quantity may or may not have an equation
//if it has an equation, it must have some boundary condition flag

//The Solution is like a vector of Quantities

class Quantity { 
  
public:
  
   Quantity(std::string name_in,QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Quantity();  

      virtual void Function_txyz(const double t, const double* xp, double* temp) const  { std::cout <<  "A quantity that calls this function must have an implementation of it" << std::endl; abort();  }  
      virtual void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const  { std::cout <<  "A quantity that calls this function must have an implementation of it" << std::endl; abort(); } 
      virtual void initialize_xyz(const double* xp, std::vector<double> & value) const  { std::cout <<  "A quantity that calls this function must have an implementation of it" << std::endl; abort(); } 
     
  void set_eqn(SystemTwo*);
  
 inline       void SetPosInAssocEqn(uint pos_in) { _pos = pos_in; return;}
 
  std::string  _name;      //quantity name, to retrieve it
  uint         _dim;       //number of scalar components
  uint         _FEord;     //FEorder
  std::vector<double>     _refvalue;  //ref values for the scalar components (_dim)
  QuantityMap& _qtymap;
  SystemTwo *      _eqn;
  uint           _pos;     //block position in the associated equation

  
};



///This class is the GLOBAL QuantityMap of my simulation,
///i.e. it holds the list of all the involved physical quantities



class QuantityMap { 
  
public:

   QuantityMap() {};
  ~QuantityMap() {};
  
  inline           void  AddQuantity(Quantity* value)          { _QuantMap.insert(make_pair(value->_name,value)); }
  
  inline       Quantity* GetQuantity(const std::string & name)  const    {
 
    std::map<std::string,Quantity*>::const_iterator myit = _QuantMap.find(name);
    
      if ( myit == _QuantMap.end() ) { 
       std::cout << "QuantityMap::GetQuantity: Sorry but there is no ---> "
                 << name  << " <--- element in the Global Quantity Map" << std::endl; abort();
      }

    return myit->second;
    
  }

    /** get/set */
  inline void SetMeshTwo(const MultiLevelMeshTwo * mesh_in)  {   _mesh = mesh_in; return; }
  
  inline const  MultiLevelMeshTwo * GetMeshTwo() const { return  _mesh; }
  
  inline const FemusInputParser<double> *  GetInputParser() const { return _physmap; }

  void SetInputParser(const FemusInputParser<double> * parser_in) { _physmap = parser_in; return; }

private:
  
 std::map<std::string,Quantity*> _QuantMap;
 
 const MultiLevelMeshTwo * _mesh;
 const FemusInputParser<double> * _physmap;
  
};



} //end namespace femus



#endif