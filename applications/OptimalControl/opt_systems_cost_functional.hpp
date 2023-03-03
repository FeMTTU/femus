#ifndef OPT_SYSTEMS_COST_FUNCTIONAL_HPP
#define OPT_SYSTEMS_COST_FUNCTIONAL_HPP




namespace femus {

    

namespace ctrl {

    
    
namespace cost_functional_Square_or_Cube {

const unsigned int axis_direction_target_reg(const unsigned int face_index) {

    unsigned int axis_dir;

        if (face_index == 1 || face_index == 2) { axis_dir = 0; }
   else if (face_index == 3 || face_index == 4) { axis_dir = 1; }
   else if (face_index == 5 || face_index == 6) { axis_dir = 2; }

    return axis_dir;

}




//*********************** Find volume elements that contain a  Target domain element **************************************

int ElementTargetFlag(const std::vector<double> & elem_center) {

    const double target_line_position_along_coordinate = TARGET_LINE_ORTHOGONAL_DISTANCE_FROM_FACE_ATTACHED_TO_TARGET_REG;
 //***** set target domain flag ******
  int target_flag = 0; //set 0 to 1 to get the entire domain

  const double offset_to_include_line = OFFSET_TO_INCLUDE_LINE;

  const unsigned int axis_dir = axis_direction_target_reg(FACE_FOR_TARGET);

  const int  target_line_sign =  ctrl:: Square_or_Cube ::sign_function_for_delimiting_region(FACE_FOR_TARGET);

   const double target_line = target_line_position_along_coordinate + target_line_sign * offset_to_include_line;



      if ((  target_line_sign * elem_center[axis_dir] < target_line_sign * target_line ) &&
          (  target_line_sign * elem_center[axis_dir] > - target_line_position_along_coordinate + target_line_sign * (target_line_position_along_coordinate - target_line_sign * offset_to_include_line)))
          {  target_flag = 1;  }

     return target_flag;

}



//******************************************* Desired Target *******************************************************

double DesiredTarget() {
   return 0.9;
}


 std::vector<double> DesiredTargetVec() {

    std::vector<double>  Vel_desired(3, 0.);

   const unsigned int axis_dir = 0;

    Vel_desired[axis_dir] = 1.;

   return Vel_desired;
    }




  
 }  //end namespace 
 

}  //end namespace ctrl

}



#include "00_cost_functional.hpp"


#endif
