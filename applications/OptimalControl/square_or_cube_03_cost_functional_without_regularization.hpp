#ifndef _SQUARE_OR_CUBE_COST_FUNCTIONAL_HPP_
#define _SQUARE_OR_CUBE_COST_FUNCTIONAL_HPP_



#include "00_cost_functional_without_regularization.hpp"


//*********************** Control, cost functional, target region - BEGIN *******************************************************
  /* Rectangular/Hexahedral domain:  1-2 x coords, 3-4 y coords, 5-6 z coords */
  /* L-shaped domain (2d):  1-2 x coords, 3-4 y coords, 5 indent between 1 and 2, 6 indent between 3 and 4 */
#define FACE_FOR_TARGET         2

#define  TARGET_LINE_ORTHOGONAL_DISTANCE_FROM_FACE_ATTACHED_TO_TARGET_REG  0.5
//*********************** Control, cost functional, target region - END *******************************************************





namespace femus {

    

namespace ctrl {

  namespace square_or_cube {

class cost_functional_without_regularization : public  femus::ctrl::cost_functional_without_regularization {

public:

static const unsigned int axis_direction_target_reg(const unsigned int face_index) {

    unsigned int axis_dir;

        if (face_index == 1 || face_index == 2) { axis_dir = 0; }
   else if (face_index == 3 || face_index == 4) { axis_dir = 1; }
   else if (face_index == 5 || face_index == 6) { axis_dir = 2; }

    return axis_dir;

}




//*********************** Find volume elements that contain a  Target domain element **************************************

static int ElementTargetFlag(const std::vector<double> & elem_center) {

    const double target_line_position_along_coordinate = TARGET_LINE_ORTHOGONAL_DISTANCE_FROM_FACE_ATTACHED_TO_TARGET_REG;
 //***** set target domain flag ******
  int target_flag = 0; //set 0 to 1 to get the entire domain

  const double offset_to_include_line = OFFSET_TO_INCLUDE_LINE;

  const unsigned int axis_dir = axis_direction_target_reg(FACE_FOR_TARGET);

  const int  target_line_sign =  ctrl:: square_or_cube:: List_of_faces ::sign_function_for_delimiting_region(FACE_FOR_TARGET);

   const double target_line = target_line_position_along_coordinate + target_line_sign * offset_to_include_line;



      if ((  target_line_sign * elem_center[axis_dir] < target_line_sign * target_line ) &&
          (  target_line_sign * elem_center[axis_dir] > - target_line_position_along_coordinate + target_line_sign * (target_line_position_along_coordinate - target_line_sign * offset_to_include_line)))
          {  target_flag = 1;  }

     return target_flag;

}



//******************************************* Desired Target *******************************************************

 static std::vector<double> DesiredTargetVec() {

    std::vector<double>  Vel_desired(3, 0.);

   const unsigned int axis_dir = 0;

    Vel_desired[axis_dir] = 1.;

   return Vel_desired;
    }




  
 };


}

}  //end namespace ctrl

}


#endif
