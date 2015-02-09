#ifndef __physicsuser__
#define __physicsuser__

#include <iostream>
#include <cstdlib>

#include "Typedefs.hpp"
#include "Quantity.hpp"

#include "OptLoop.hpp"
#include "TimeLoop.hpp"


using namespace femus;

class Temperature : public Quantity {

  public:
   Temperature(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Temperature(){};
  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 
//specific function
  //this is the function of the IMPOSED DERIVATIVE of TEMPERATURE, aka heat flux
  void heatflux_txyz(const double t,const double* xyz, double* qflux) const;
  

};

//===============
class TempLift : public Quantity {

  public:
   TempLift(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in, const TimeLoop & time_loop_in);
  ~TempLift(){};
  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
  
    const TimeLoop & _my_timeloop;
  

};


//===============
class TempAdj : public Quantity {

  public:
   TempAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~TempAdj(){};
  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
  

};

//===============
class TempDes : public Quantity {

  public:
   TempDes(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~TempDes(){};
  void Function_txyz(const double t, const double* xp,double* temp) const;  
  

};



class Pressure : public Quantity {

  public:
   Pressure(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Pressure(){};
  
  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;

 

};


class Velocity : public Quantity {

  public:
  Velocity(std::string name_in, QuantityMap& qtymap_in);
    Velocity(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  ~Velocity(){};

  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;

  //specific function
  //this is the STRAIN DERIVATIVE of VELOCITY, so it must stay here
  //from the physical and also mathematical point of view
  //the shape funcs of the same order as v will then be used and so on
//multi-dimensional arrays must have bounds for all dimensions except the first
    void strain_txyz(const double t, const double* xp,double strain[][DIMENSION]) const;  //TODO convert this double array

};


class Pressure2 : public Quantity {

  public:
  Pressure2(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Pressure2(){};
  
  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;

};


#endif