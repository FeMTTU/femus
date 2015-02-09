#ifndef __physicsuser__
#define __physicsuser__

#include <iostream>
#include <cstdlib>

#include "Typedefs.hpp"
#include "Quantity.hpp"

#include "OptLoop.hpp"

using namespace femus;



class MagnFieldHom : public Quantity {

  public:
    
   MagnFieldHom(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~MagnFieldHom(){};
  
  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldHomAdj : public Quantity {

  public:
    
   MagnFieldHomAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~MagnFieldHomAdj(){};
  
  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldExt : public Quantity {

  public:
    
   MagnFieldExt(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~MagnFieldExt(){};
  
  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldHomLagMult : public Quantity {

  public:
    
   MagnFieldHomLagMult(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~MagnFieldHomLagMult(){};
  
  void Function_txyz(const double t, const double* xp,double* temp) const;
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};


class MagnFieldHomLagMultAdj : public Quantity {

  public:
    
   MagnFieldHomLagMultAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~MagnFieldHomLagMultAdj(){};
  
  void Function_txyz(const double t, const double* xp,double* temp) const;
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};


class MagnFieldExtLagMult : public Quantity {

  public:
    
   MagnFieldExtLagMult(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~MagnFieldExtLagMult(){};
   
  void Function_txyz(const double t, const double* xp,double* temp) const;
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};


class Pressure : public Quantity {

  public:
    
   Pressure(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Pressure(){};

  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};


class PressureAdj : public Quantity {

  public:
    
   PressureAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~PressureAdj(){};

  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};


class Velocity : public Quantity {

  public:
    
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
    void strain_txyz(const double t, const double* xp,double strain[][DIMENSION]) const;  //TODO remove this DIMENSION and use double pointers or std vector

  

};

class VelocityAdj : public Quantity {

  public:
    
   VelocityAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~VelocityAdj(){};

  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;


};


class DesVelocity : public Quantity {

  public:
    
   DesVelocity(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~DesVelocity(){};

  void Function_txyz(const double t, const double* xp,double* temp) const;  
 

};

//===============================
//temp-dep =========
//===============================

class Temperature : public Quantity {

  public:
    
   Temperature(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Temperature(){};
  
  void Function_txyz(const double t, const double* xp,double* temp) const;  
 
//specific function
  //this is the function of the IMPOSED DERIVATIVE of TEMPERATURE, aka heat flux
  void heatflux_txyz(const double t,const double* xyz, double* qflux) const;
  

};


class Density : public Quantity {

  public:
    
   Density(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Density(){};

  void Function_txyz(const double t, const double* xp,double* temp) const{};  
 
void Temp_dep(const double temp_in, double& rho_out) const {rho_out = 1.;return;} 

};

class Viscosity : public Quantity {

  public:
    
   Viscosity(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~Viscosity(){};

  void Function_txyz(const double t, const double* xp,double* temp) const{};  
 
void Temp_dep(const double temp_in, double& mu_out) const {mu_out = 1.;return;} 

};


class HeatConductivity : public Quantity {

  public:
    
   HeatConductivity(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~HeatConductivity(){};

  void Function_txyz(const double t, const double* xp,double* temp) const{};  
 
void Temp_dep(const double temp_in, double& mu_out) const {mu_out = 1.;return;} 

};

class SpecificHeatP : public Quantity {

  public:
    
   SpecificHeatP(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  ~SpecificHeatP(){};

  void Function_txyz(const double t, const double* xp,double* temp) const{};  
 
void Temp_dep(const double temp_in, double& mu_out) const {mu_out = 1.;return;} 

};



#endif

