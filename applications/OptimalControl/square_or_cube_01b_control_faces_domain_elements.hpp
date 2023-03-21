#ifndef CONTROL_FACES_SQUARE_OR_CUBE_DOMAIN_ELEMENTS_HPP
#define CONTROL_FACES_SQUARE_OR_CUBE_DOMAIN_ELEMENTS_HPP

#include "square_or_cube_01_control_faces.hpp"


//*********************** Control, boundary - BEGIN *******************************************************
#define BOUNDARY_ORTHOGONAL_DISTANCE_FROM_GAMMA_C  1.   //how far it goes orthogonally to the Control piece of the Boundary 
//*********************** Control, boundary - END *******************************************************


//*********************** Control, lifting internal - BEGIN *******************************************************
#define LIFTING_INTERNAL_ORTHOGONAL_DISTANCE_FROM_GAMMA_C  1.   //how far it goes orthogonally to the Control piece of the Boundary 
#define LIFTING_INTERNAL_WIDTH_LOWER  0.
#define LIFTING_INTERNAL_WIDTH_UPPER  1.
//*********************** Control, lifting internal - END *******************************************************



namespace femus {

namespace ctrl {


namespace  square_or_cube {
 
template < class T >
class Domain_elements_containing_Gamma_control  {

  public:


static const double face_coordinate_extreme_position_normal_to_Gamma_control(const unsigned int face_index) {

  double extreme_pos;

        if (face_index == 1 || face_index == 3 || face_index == 5) {  extreme_pos = 0.; }
   else if (face_index == 2 || face_index == 4 || face_index == 6) {  extreme_pos = 1.; }

   return extreme_pos;

}





//*********************** Find volume elements that contain a Control domain element *********************************

static int ControlDomainFlag_external_restriction(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
 // flag = 1: we are in the lifting nonzero domain
  int exterior_el_flag = 0.;
   if ( elem_center[0] >  1. -   OFFSET_TO_INCLUDE_LINE ) { exterior_el_flag = 1; }

     return exterior_el_flag;

}



//*********************** Find volume elements that contain a Control domain element *********************************

// template < class T >
static int ControlDomainFlag_internal_restriction(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
 // flag = 1: we are in the lifting nonzero domain
  int control_el_flag = 0.;

  const double offset_to_include_line =  OFFSET_TO_INCLUDE_LINE;

  const double control_domain_depth = LIFTING_INTERNAL_ORTHOGONAL_DISTANCE_FROM_GAMMA_C;

  const double control_domain_width_lower = LIFTING_INTERNAL_WIDTH_LOWER;
  const double control_domain_width_upper = LIFTING_INTERNAL_WIDTH_UPPER;


	  for(unsigned f = 0; f < /*ctrl::*/ T ::_face_with_extremes_index_size; f++) {

   const int  line_sign =  /*ctrl::*/ T ::  sign_function_for_delimiting_region(/*ctrl::*/ T ::_face_with_extremes_index[f]);  //1,3,5 = 1 || 2,4,6 = -1

   const double extreme_pos = face_coordinate_extreme_position_normal_to_Gamma_control(/*ctrl::*/ T ::_face_with_extremes_index[f]);        //1,3,5 = 0 || 2,4,6 = +1

   const unsigned int axis_dir = T ::tangential_direction_to_Gamma_control(/*ctrl::*/ T ::_face_with_extremes_index[f]);                  // 1-2 , 5-6 = 1 || 3-4 = 0


   if ( ( line_sign * elem_center[1 - axis_dir] <   line_sign * ( extreme_pos + line_sign * control_domain_depth ) )
       && ( elem_center[axis_dir] > control_domain_width_lower - offset_to_include_line )
       && ( elem_center[axis_dir] < control_domain_width_upper + offset_to_include_line ) )
      { control_el_flag = 1; }

      }

     return control_el_flag;

}



 static bool volume_elem_contains_a_Gamma_control_face( const std::vector<double> & elem_center ) {

      int control_flag_jel = 0;
        control_flag_jel =  ControlDomainFlag_bdry(elem_center);

        return (control_flag_jel == 1);


  }


//*********************** Find volume elements that contain a Control Face element *********************************

// template < class T >
static int ControlDomainFlag_bdry(const std::vector<double> & elem_center) {


  int control_el_flag = 0;

  const double offset_to_include_line = OFFSET_TO_INCLUDE_LINE;

  const double control_domain_depth = BOUNDARY_ORTHOGONAL_DISTANCE_FROM_GAMMA_C; //this picks a lot more elements, but then the if on the faces only gets the control boundary


	  for(unsigned f = 0; f < /*ctrl::*/ T ::_face_with_extremes_index_size; f++) {

   const int  line_sign =  /*ctrl::*/ T ::sign_function_for_delimiting_region(/*ctrl::*/ T ::_face_with_extremes_index[f]);

   const double extreme_pos = face_coordinate_extreme_position_normal_to_Gamma_control(/*ctrl::*/ T ::_face_with_extremes_index[f]);

   const unsigned int Gamma_c_dir_tangential = T ::tangential_direction_to_Gamma_control(/*ctrl::*/ T ::_face_with_extremes_index[f]);


   if ( ( line_sign * elem_center[1 - Gamma_c_dir_tangential] <   line_sign * (  extreme_pos  + line_sign * control_domain_depth) )
       && ( elem_center[Gamma_c_dir_tangential] >/* ctrl::*/ T ::_face_with_extremes_extremes[f][0] - offset_to_include_line )
       && ( elem_center[Gamma_c_dir_tangential] < /*ctrl::*/ T ::_face_with_extremes_extremes[f][1] + offset_to_include_line ) )
      { control_el_flag = 1; }

   }



     return control_el_flag;
}







};



}


}

}


#endif
