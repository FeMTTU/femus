#ifndef __quantity_h
#define __quantity_h

#include <string>

#include "Typedefs.hpp"
#include "QuantityLocal.hpp"

class EqnBase;
class Physics;
class QuantityMap;

class Quantity { 
  
public:
  
  std::string  _name;      //quantity name, to retrieve it
  uint         _dim;       //number of scalar components
  uint         _FEord;     //FEorder
  double *     _refvalue;  //ref values for the scalar components (_dim)
  QuantityMap& _qtymap;
  EqnBase *      _eqn;
  uint           _pos;     //block position in the associated equation

   Quantity(std::string name_in,QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Quantity();  

  /*virtual*/ void FunctionDof(const uint vb, QuantityLocal& myvect, const double t,const double* xx) const/* =0*/;
      virtual void Function_txyz(const double t, const double* xp, double* temp) const   = 0;  
//its good to have them as pure virtual, so that the user has to implement them somehow
//in case no equation is inserted
//on the other hand these functions are not used if the quantities are not taken for the dofs  
//but for instance for temperature dependence  
//i decide to put them as virtual so the child that doesnt need them may not worry about having them... No,
// you better do PURE virtual so that the child doesnt risk picking a bad default from the father
//i wanna do a COMMON FunctionDof, not even virtual. It is just the same for all, and it cant be 
  //overwritten
  
  void set_eqn(EqnBase*);
 inline       void SetPosInAssocEqn(uint pos_in) { _pos = pos_in; return;}
  
};


////////////////////////
//C++ includes
#include <map>
#include <string>
#include <iostream>
#include <cstdlib>

///This class is the GLOBAL QuantityMap of my simulation,
///i.e. it holds the list of all the involved physical quantities

class QuantityMap { 
  
public:

  std::map<std::string,Quantity*> _QuantMap;   // system map
 Physics& _phys; 
  
   QuantityMap(Physics& phys_in);
  ~QuantityMap(){};

  
  inline           void  set_qty(Quantity* value)          {_QuantMap.insert(make_pair(value->_name,value));}
  inline       Quantity* get_qty(const std::string & name)      {

//let me write this so that I have a more friendly message even in optimized mode
//I will no longer have the message "attempting to dereference a past-the-end iterator..."
//because i'll put a check before dereferencing it!
//what if i just do std::iterator? what is the difference with std::map<>::iterator?
//it says "missing template arguments"
//it means that if i take values from a templated class, like QuantMap,
//then the return of the function find is templated also, and so must be the variable where this return is stored.
// iterator is a type defined within the templated class std::map
//that is why you define this type like that
//here, the first :: is for a "NAMESPACE SCOPE", the second is for a "CLASS SCOPE"

std::map<std::string,Quantity*>::iterator myit = _QuantMap.find(name);
      if ( myit == _QuantMap.end() ) { 
       std::cout << "QuantityMap::get_qty: Sorry but there is no ---> "
                 << name  << " <--- element in the Global Quantity Map" << std::endl; abort();}

return /*_QuantMap.find(name)*/myit->second;
    
  }

//   inline const EqnBase* get_eqs(const string & name) const {return _equations.find(name)->second;}

//why does anyone need the const version of that?
//if you return a const thing, then you cannot modify it
//if you return a non-const thing, then you may modify it
//the problem about modifying is that you may modify the DATA of the class,
//not clearly the functions... (what if you could overload functions at runtime... :) )
//but what is the difference between making a class non-const and with protected members
//or making a class const with public members?
//const is a cv-qualifier
//public,protected,private are...
  //they are used to RESTRICT the CONTEXT where some data/functions of the class can be called

  //can CONST or private/public/protected help me in avoiding WRITING ON TOP OF a VARIABLE at RUN-TIME?
  //let us try that with a small struct!
};





#endif