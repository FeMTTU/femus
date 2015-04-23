#ifndef __femus_mhdopt_OptQuantities_hpp__
#define __femus_mhdopt_OptQuantities_hpp__

#include <iostream>
#include <cstdlib>

#include "Typedefs.hpp"
#include "Quantity.hpp"

#include "OptLoop.hpp"

using namespace femus;



class MagnFieldHomX : public Quantity {

  public:
    
   MagnFieldHomX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldHomY : public Quantity {

  public:
    
   MagnFieldHomY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldHomZ : public Quantity {

  public:
    
   MagnFieldHomZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldHomAdjX : public Quantity {

  public:
    
   MagnFieldHomAdjX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldHomAdjY : public Quantity {

  public:
    
   MagnFieldHomAdjY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldHomAdjZ : public Quantity {

  public:
    
   MagnFieldHomAdjZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldExtX : public Quantity {

  public:
    
   MagnFieldExtX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldExtY : public Quantity {

  public:
    
   MagnFieldExtY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldExtZ : public Quantity {

  public:
    
   MagnFieldExtZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};

class MagnFieldHomLagMult : public Quantity {

  public:
    
   MagnFieldHomLagMult(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};


class MagnFieldHomLagMultAdj : public Quantity {

  public:
    
   MagnFieldHomLagMultAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};


class MagnFieldExtLagMult : public Quantity {

  public:
    
   MagnFieldExtLagMult(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);
   
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};


class Pressure : public Quantity {

  public:
    
   Pressure(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void Function_txyz(const double t, const double* xp,double* temp) const;  
  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};


class PressureAdj : public Quantity {

  public:
    
   PressureAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;
 

};


class VelocityX : public Quantity {

  public:
    
   VelocityX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;

};


class VelocityY : public Quantity {

  public:
    
   VelocityY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;

};


class VelocityZ : public Quantity {

  public:
    
   VelocityZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;

};




class VelocityAdjX : public Quantity {

  public:
    
  VelocityAdjX(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;


};


class VelocityAdjY : public Quantity {

  public:
    
  VelocityAdjY(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;


};


class VelocityAdjZ : public Quantity {

  public:
    
  VelocityAdjZ(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in);

  void bc_flag_txyz(const double t, const double* xp, std::vector<int> & flag) const;
  void initialize_xyz(const double* xp, std::vector<double> & value) const;


};

#endif

