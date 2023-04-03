#ifndef SQUARE_OR_CUBE_BOUNDARY_CONDITIONS_HPP
#define SQUARE_OR_CUBE_BOUNDARY_CONDITIONS_HPP

#include "square_or_cube_01_control_faces.hpp"

namespace femus {

namespace ctrl {

  namespace square_or_cube {

class Boundary_Condition {

    public:


 static double opposite_face_ctrl_or_state_value(const unsigned int face_index, const double domain_length) {

            double opposite_value = 0.;

                 if (face_index == 1 || face_index == 3 || face_index == 5) { opposite_value = 0.; }
            else if (face_index == 2 || face_index == 4 || face_index == 6) { opposite_value = domain_length; }

             return opposite_value;

         }
         
    static constexpr double _scale_factor_max_value_on_line = 5.;     


};
  }
  
  
  namespace square {
  
  
namespace boundary_control_between_extreme {

template < class LIST_OF_CTRL_FACES >
 class Boundary_Condition_control_between_extreme: public square_or_cube :: Boundary_Condition {

    public:

      // template < class T >
 static bool ctrl_or_state_set_dirichlet_flags(const MultiLevelProblem * ml_prob,
                                        const int faceName,
                                        const std::vector < double > & x,
                                        bool &  dirichlet)  {


//      assert( _face_with_extremes_index_size == 1 );
//     const unsigned face_for_control = _face_with_extremes_index[0];

      /*static*/ bool  is_facename_a_control_face = false;
      for(unsigned f = 0; f < /*ctrl::*/LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

           if (faceName == /*ctrl::*/LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]) { is_facename_a_control_face = true; break; }

           }

           if  (is_facename_a_control_face == false) { dirichlet = true; }


           for(unsigned f = 0; f < /*ctrl::*/LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

           if (faceName == /*ctrl::*/LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]) {

              if ( !(x[ /*ctrl::*/LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName) ] >
                  /*ctrl::*/LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][LIST_OF_CTRL_FACES :: _num_of_tang_components_per_face_2d-1][0] +  OFFSET_TO_INCLUDE_LINE  &&
                     x[ /*ctrl::*/LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName) ] <
                  /*ctrl::*/LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][LIST_OF_CTRL_FACES :: _num_of_tang_components_per_face_2d-1][1] -  OFFSET_TO_INCLUDE_LINE
                 ) ) {
                      dirichlet = true;
                 }
              }
            }


          return dirichlet;
        }

 };

template < class LIST_OF_CTRL_FACES >
 class Multiple_controls_and_homogeneous_boundary_conditions : public Boundary_Condition_control_between_extreme < LIST_OF_CTRL_FACES > {

 public:



 static double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {


   value = 0 + PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

   return value;
            }
};




template < class LIST_OF_CTRL_FACES >
 class Multiple_controls_in_front_constant : public Boundary_Condition_control_between_extreme < LIST_OF_CTRL_FACES > {

 public:

static double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

//      if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

//      assert( _face_with_extremes_index_size == 2 );

     const unsigned face_for_control_principal = /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_index[0];

    const double domain_length = 1.;

      const double gamma = square_or_cube :: Boundary_Condition::_scale_factor_max_value_on_line;


	  for(unsigned f = 0; f < /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

        if (faceName == /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]) {
              if(faceName == face_for_control_principal) {  value = 0.; }
              else if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){
                     if ( !(x[ /*ctrl::*/ LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName) ] > /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][LIST_OF_CTRL_FACES :: _num_of_tang_components_per_face_2d-1][0] +  OFFSET_TO_INCLUDE_LINE  &&
                            x[ /*ctrl::*/ LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName) ] <
                     /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][LIST_OF_CTRL_FACES :: _num_of_tang_components_per_face_2d-1][1] -  OFFSET_TO_INCLUDE_LINE ) )   { value =  gamma * domain_length; }
                     else {value = 0.; }
              }
              else {
                     if ( !(x[ /*ctrl::*/ LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName) ] > /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][LIST_OF_CTRL_FACES :: _num_of_tang_components_per_face_2d-1][0] +  OFFSET_TO_INCLUDE_LINE  &&
                            x[ /*ctrl::*/ LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName) ] <
                     /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][LIST_OF_CTRL_FACES :: _num_of_tang_components_per_face_2d-1][1] -  OFFSET_TO_INCLUDE_LINE ) ) { value =
                            gamma * (
                                      /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(face_for_control_principal) *  x[      /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(face_for_control_principal) ] ); }
                     else {value = 0.; }
              }
        }
        else if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)) { value =  gamma * domain_length; }
        else      { value = gamma * (
                                    /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(face_for_control_principal) *  x[      /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(face_for_control_principal) ] ); }

      }

   value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

   return value;
             }
};




template < class LIST_OF_CTRL_FACES >
 class Single_control_in_front_linear : public Boundary_Condition_control_between_extreme < LIST_OF_CTRL_FACES > {

 public:


 static const unsigned int adjacent_face(const unsigned int face_index) {

   unsigned int adjacent_face = 0;

   if (face_index < 3 )    { adjacent_face = face_index + 2; }
   else                    { adjacent_face = face_index - 2; }

    return adjacent_face;

}

 static double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

    assert( /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size == 1 );

    const unsigned face_for_control = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[0];
    const unsigned adjacent_face_for_control = adjacent_face(face_for_control);


    const double domain_length = 1.;

    const double gamma = square_or_cube :: Boundary_Condition::_scale_factor_max_value_on_line;


    if(faceName == face_for_control) {  value = 0.; }
    else if(faceName == adjacent_face_for_control) {  value = 0.; }
    else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control)){ value = gamma * (
             /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(adjacent_face_for_control, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(adjacent_face_for_control) *  x[
             /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(adjacent_face_for_control) ] );  }
    else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){ value = gamma * (
             /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(face_for_control, domain_length) +
              /*ctrl::*/  LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(face_for_control) *  x[
             /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(face_for_control) ] );  }

   value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

   return value;
                }
};




template < class LIST_OF_CTRL_FACES >
 class Double_controls_adjacent_in_front_linear : public Boundary_Condition_control_between_extreme < LIST_OF_CTRL_FACES > {

 public:



 static double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

     if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

     assert( /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size == 2 );

     const unsigned face_for_control_principal = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[0];
     const unsigned adjacent_face_for_control = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[1];

     const double domain_length = 1.;

     const double gamma = square_or_cube :: Boundary_Condition::_scale_factor_max_value_on_line;


     for(unsigned f = 0; f < /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

        if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f])     {  value = 0.; }
        else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){ value = gamma * (
             /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(adjacent_face_for_control, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(adjacent_face_for_control) *  x[
             /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(adjacent_face_for_control) ] );  }
        else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){ value = gamma * (
             /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(face_for_control_principal) *  x[
             /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(face_for_control_principal) ] );  }
     }

     value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

     return value;
}
 };



template < class LIST_OF_CTRL_FACES >
 class Triple_controls_adjacent_in_front_linear : public Boundary_Condition_control_between_extreme < LIST_OF_CTRL_FACES > {

 public:



 static double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

     if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

     assert( /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size == 3 );

     const unsigned face_for_control_principal = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[0];
     const unsigned adjacent_face_for_control = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[1];

     const double domain_length = 1.;

     const double gamma = square_or_cube :: Boundary_Condition::_scale_factor_max_value_on_line;


     for(unsigned f = 0; f < /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

         if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]) {

                if ( (faceName == face_for_control_principal) || (faceName == adjacent_face_for_control) )     {  value = 0.; }
                else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){
                        if ( !(x[ /*ctrl::*/ LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName) ] > /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][LIST_OF_CTRL_FACES :: _num_of_tang_components_per_face_2d-1][0] +  OFFSET_TO_INCLUDE_LINE  &&
                               x[ /*ctrl::*/ LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName) ] < /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][LIST_OF_CTRL_FACES :: _num_of_tang_components_per_face_2d-1][1] -  OFFSET_TO_INCLUDE_LINE ) ) {
                                  value = gamma * (
                                  /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(adjacent_face_for_control, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(adjacent_face_for_control) *  x[
                                  /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(adjacent_face_for_control) ] );  }
                        else {value = 0.; }
                }
                else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){
                        if ( !(x[ /*ctrl::*/ LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName) ] > /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][LIST_OF_CTRL_FACES :: _num_of_tang_components_per_face_2d-1][0] +  OFFSET_TO_INCLUDE_LINE  &&
                               x[ /*ctrl::*/ LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName) ] < /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][LIST_OF_CTRL_FACES :: _num_of_tang_components_per_face_2d-1][1] -  OFFSET_TO_INCLUDE_LINE ) ) {
                               value = gamma * (
                               /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(face_for_control_principal) *  x[
                               /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(face_for_control_principal) ] );  }
                               else {value = 0.; }
                }
         }
         else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){ value = gamma * (
                                  /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(adjacent_face_for_control, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(adjacent_face_for_control) *  x[
                                  /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(adjacent_face_for_control) ] );  }

         else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){ value = gamma * (
                            /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(face_for_control_principal) *  x[
                            /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(face_for_control_principal) ] );  }

     }

     value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

     return value;
}






 };





}

namespace boundary_control_full_face {

template < class LIST_OF_CTRL_FACES >
 class Boundary_Condition_control_full_face: public square_or_cube :: Boundary_Condition {

    public:

      // template < class T >
 static bool ctrl_or_state_set_dirichlet_flags(const MultiLevelProblem * ml_prob,
                                        const int faceName,
                                        const std::vector < double > & x,
                                        bool &  dirichlet)  {


//      assert( _face_with_extremes_index_size == 1 );
//     const unsigned face_for_control = _face_with_extremes_index[0];

      /*static*/ bool  is_facename_a_control_face = false;
      for(unsigned f = 0; f < /*ctrl::*/LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

           if (faceName == /*ctrl::*/LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]) { is_facename_a_control_face = true; break; }

           }

           if  (is_facename_a_control_face == false) { dirichlet = true; }



          return dirichlet;
        }

 };

template < class LIST_OF_CTRL_FACES >
 class Multiple_controls_and_homogeneous_boundary_conditions : public Boundary_Condition_control_full_face < LIST_OF_CTRL_FACES > {

 public:



 static double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {


   value = 0 + PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

   return value;
            }
};




template < class LIST_OF_CTRL_FACES >
 class Multiple_controls_in_front_constant : public Boundary_Condition_control_full_face < LIST_OF_CTRL_FACES > {

 public:

static double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

//      if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

//      assert( _face_with_extremes_index_size == 2 );

     const unsigned face_for_control_principal = /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_index[0];

    const double domain_length = 1.;

      const double gamma = square_or_cube :: Boundary_Condition::_scale_factor_max_value_on_line;


	  for(unsigned f = 0; f < /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

        if (faceName == /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]) {value = 0.;}

        else if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)) { value =  gamma * domain_length; }
        else      { value = gamma * (
                                    /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(face_for_control_principal) *  x[      /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(face_for_control_principal) ] ); }

      }

   value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

   return value;
             }
};




template < class LIST_OF_CTRL_FACES >
 class Single_control_in_front_linear : public Boundary_Condition_control_full_face < LIST_OF_CTRL_FACES > {

 public:


 static const unsigned int adjacent_face(const unsigned int face_index) {

   unsigned int adjacent_face = 0;

   if (face_index < 3 )    { adjacent_face = face_index + 2; }
   else                    { adjacent_face = face_index - 2; }

    return adjacent_face;

}

 static double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

    assert( /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size == 1 );

    const unsigned face_for_control = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[0];
    const unsigned adjacent_face_for_control = adjacent_face(face_for_control);


    const double domain_length = 1.;

    const double gamma = square_or_cube :: Boundary_Condition::_scale_factor_max_value_on_line;


    if(faceName == face_for_control) {  value = 0.; }
    else if(faceName == adjacent_face_for_control) {  value = 0.; }
    else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control)){ value = gamma * (
             /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(adjacent_face_for_control, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(adjacent_face_for_control) *  x[
             /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(adjacent_face_for_control) ] );  }
    else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){ value = gamma * (
             /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(face_for_control, domain_length) +
              /*ctrl::*/  LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(face_for_control) *  x[
             /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(face_for_control) ] );  }

   value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

   return value;
                }
};




template < class LIST_OF_CTRL_FACES >
 class Double_controls_adjacent_in_front_linear : public Boundary_Condition_control_full_face < LIST_OF_CTRL_FACES > {

 public:



 static double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

     if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

     assert( /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size == 2 );

     const unsigned face_for_control_principal = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[0];
     const unsigned adjacent_face_for_control = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[1];

     const double domain_length = 1.;

     const double gamma = square_or_cube :: Boundary_Condition::_scale_factor_max_value_on_line;


     for(unsigned f = 0; f < /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

        if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f])     {  value = 0.; }
        else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){ value = gamma * (
             /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(adjacent_face_for_control, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(adjacent_face_for_control) *  x[
             /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(adjacent_face_for_control) ] );  }
        else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){ value = gamma * (
             /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(face_for_control_principal) *  x[
             /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(face_for_control_principal) ] );  }
     }

     value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

     return value;
}
 };



template < class LIST_OF_CTRL_FACES >
 class Triple_controls_adjacent_in_front_linear : public Boundary_Condition_control_full_face < LIST_OF_CTRL_FACES > {

 public:



 static double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

     if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

     assert( /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size == 3 );

     const unsigned face_for_control_principal = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[0];
     const unsigned adjacent_face_for_control = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[1];

     const double domain_length = 1.;

     const double gamma = square_or_cube :: Boundary_Condition::_scale_factor_max_value_on_line;


     for(unsigned f = 0; f < /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

         if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]) { value = 0.; }

         else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){ value = gamma * (
                                  /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(adjacent_face_for_control, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(adjacent_face_for_control) *  x[
                                  /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(adjacent_face_for_control) ] );  }

         else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){ value = gamma * (
                            /*ctrl::boundary_conditions::*/square_or_cube :: Boundary_Condition::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) +  /*ctrl::*/ LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(face_for_control_principal) *  x[
                            /*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(face_for_control_principal) ] );  }

     }

     value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

     return value;
}






 };





}

  }
  
  
  namespace cube {
      
    namespace boundary_control_between_extreme {
        
template < class LIST_OF_CTRL_FACES >
 class Multiple_controls_and_homogeneous_boundary_conditions : public femus:: ctrl :: square :: boundary_control_between_extreme :: Multiple_controls_and_homogeneous_boundary_conditions< LIST_OF_CTRL_FACES > {

 };
 
        
        
      }
    namespace boundary_control_full_face {
    }  
      
  }  
  
  
  

   }

}


#endif
