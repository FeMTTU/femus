#ifndef NSOPT_PARAMETERS
#define NSOPT_PARAMETERS


#include   "../boundary_control_inequality/param.hpp"

namespace femus {

//*********************** Physics *******************************************************
#define FLUID_DENSITY  1

 double force[3] = {0., 0., 0.};

 //******************************************* Desired Target*******************************************************
 std::vector<double> DesiredTargetVel() {
     
    std::vector<double>  Vel_desired(3, 0.);
    
   const unsigned int axis_dir = axis_direction_Gamma_control(FACE_FOR_CONTROL);
   
    Vel_desired[axis_dir] = 1.;
    
   return Vel_desired;
}

 
//******************************** switch between stokes and navier stokes *********************************************
 
 const int advection_flag = 0;
 const int advection_Picard = 0;
 
//  Stokes: advection_flag = 0; advection_Picard = 0;
//  Newton: advection_flag = 1; advection_Picard = 0;
//  Picard: advection_flag = 1; advection_Picard = 1;
 
 
} //end namespace

 
#endif
