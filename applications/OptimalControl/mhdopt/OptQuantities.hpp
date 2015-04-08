#ifndef __femus_mhdopt_OptQuantities_hpp__
#define __femus_mhdopt_OptQuantities_hpp__

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


class VelocityX : public Quantity {

  public:
    
   VelocityX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;

};


class VelocityY : public Quantity {

  public:
    
   VelocityY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;

};


class VelocityZ : public Quantity {

  public:
    
   VelocityZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;

};




class VelocityAdjX : public Quantity {

  public:
    
  VelocityAdjX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;


};


class VelocityAdjY : public Quantity {

  public:
    
  VelocityAdjY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;


};


class VelocityAdjZ : public Quantity {

  public:
    
  VelocityAdjZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;


};


class DesVelocityX : public Quantity {

  public:
    
   DesVelocityX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void Function_txyz(const double t, const double* xp,double* temp) const;  
 

};

class DesVelocityY : public Quantity {

  public:
    
   DesVelocityY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void Function_txyz(const double t, const double* xp,double* temp) const;  
 

};

class DesVelocityZ : public Quantity {

  public:
    
   DesVelocityZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void Function_txyz(const double t, const double* xp,double* temp) const;  
 

};

#endif

