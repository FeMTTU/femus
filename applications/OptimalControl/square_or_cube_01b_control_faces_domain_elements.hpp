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
 
//*********************** Find volume elements that contain a Control domain element *********************************
template < class LIST_OF_CTRL_FACES >
class Domain_elements_containing_Gamma_control  {

};

//*********************** Pure boundary - BEGIN *********************************
template < class LIST_OF_CTRL_FACES >
class pure_boundary : public Domain_elements_containing_Gamma_control<LIST_OF_CTRL_FACES>  {
    
  public:
      


  static const std::string _cont_reg_name;    

 static bool volume_elem_contains_a_Gamma_control_face(const Solution * sol, const Mesh * msh, const unsigned int iel ) {
    
     
     const std::string sol_for_cont_reg = pure_boundary<LIST_OF_CTRL_FACES>::_cont_reg_name;
     const unsigned sol_for_cont_reg_index = sol->GetIndex(sol_for_cont_reg.c_str());
     const unsigned sol_for_cont_reg_fe_type = sol->GetSolutionType( sol_for_cont_reg_index );
     
     const unsigned placeholder_index = 0;
     
     const unsigned sol_for_cont_reg_dof = msh->GetSolutionDof(placeholder_index, iel, sol_for_cont_reg_fe_type);
     
     const bool does_iel_contain_Gamma_c = (bool) (*sol->_Sol[sol_for_cont_reg_index])(sol_for_cont_reg_dof);

//      return true;
        return (does_iel_contain_Gamma_c == true);
        
  }

};

template < class LIST_OF_CTRL_FACES >
  const std::string pure_boundary<LIST_OF_CTRL_FACES>::_cont_reg_name = "ContReg";    
//*********************** Pure boundary - END *********************************


//*********************** Lifting internal - BEGIN *********************************

template < class LIST_OF_CTRL_FACES >
class lifting_internal : public Domain_elements_containing_Gamma_control<LIST_OF_CTRL_FACES>  {
    
  public:
      

static const double face_coordinate_extreme_position_normal_to_Gamma_control(const unsigned int face_index) {

  double extreme_pos;

        if (face_index == 1 || face_index == 3 || face_index == 5) {  extreme_pos = 0.; }
   else if (face_index == 2 || face_index == 4 || face_index == 6) {  extreme_pos = 1.; }

   return extreme_pos;

}



static int ControlDomainFlag_internal_restriction(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
 // flag = 1: we are in the lifting nonzero domain
  int control_el_flag = 0.;

  const double offset_to_include_line =  OFFSET_TO_INCLUDE_LINE;

  const double control_domain_depth = LIFTING_INTERNAL_ORTHOGONAL_DISTANCE_FROM_GAMMA_C;

  const double control_domain_width_lower = LIFTING_INTERNAL_WIDTH_LOWER;
  const double control_domain_width_upper = LIFTING_INTERNAL_WIDTH_UPPER;


	  for(unsigned f = 0; f < /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

   const int  line_sign =  /*ctrl::*/ LIST_OF_CTRL_FACES ::  sign_function_for_delimiting_region(/*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]);  //1,3,5 = 1 || 2,4,6 = -1

   const double extreme_pos = face_coordinate_extreme_position_normal_to_Gamma_control(/*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]);        //1,3,5 = 0 || 2,4,6 = +1

   const unsigned int normal_dir =  LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(/*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]);

   const unsigned int tang_dir = 1 - normal_dir;

   if ( ( line_sign * elem_center[normal_dir] <   line_sign * ( extreme_pos + line_sign * control_domain_depth ) )
       /*&& ( elem_center[tang_dir] > control_domain_width_lower - offset_to_include_line )
       && ( elem_center[tang_dir] < control_domain_width_upper + offset_to_include_line )*/ )
      { control_el_flag = 1; }

      }

     return control_el_flag;

}

    
};
//*********************** Lifting internal - END *********************************


//*********************** Lifting external - BEGIN *********************************

template < class LIST_OF_CTRL_FACES >
class lifting_external : public Domain_elements_containing_Gamma_control<LIST_OF_CTRL_FACES>  {
    
  public:
      


static int ControlDomainFlag_external_restriction(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
 // flag = 1: we are in the lifting nonzero domain
  int exterior_el_flag = 0.;
   if ( elem_center[0] >  1. -   OFFSET_TO_INCLUDE_LINE ) { exterior_el_flag = 1; }

     return exterior_el_flag;

}



};
//*********************** Lifting external - END *********************************


}


}

}


#endif
