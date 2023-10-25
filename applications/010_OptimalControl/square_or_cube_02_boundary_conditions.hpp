#ifndef SQUARE_OR_CUBE_BOUNDARY_CONDITIONS_HPP
#define SQUARE_OR_CUBE_BOUNDARY_CONDITIONS_HPP

#include "square_or_cube_01_control_faces.hpp"

namespace femus {

namespace ctrl {

  namespace square_or_cube {

template < class LIST_OF_CTRL_FACES >
class Boundary_Condition {

    public:


 static double opposite_face_ctrl_or_state_value(const unsigned int face_index, const double domain_length) {

            double opposite_value = 0.;

                 if (face_index == 1 || face_index == 3 || face_index == 5) { opposite_value = 0.; }
            else if (face_index == 2 || face_index == 4 || face_index == 6) { opposite_value = domain_length; }

             return opposite_value;

         }
         

 static constexpr double _scale_factor_max_value_on_line = 5.;



 static  std::pair< bool, bool >  is_boundary_of_boundary_region (/*const MultiLevelProblem * ml_prob,*/
                                        const int faceName,
                                        const std::vector < double > & x,
                                        unsigned number_of_tangential_direction_components
                                        ){

    bool is_facename_a_control_face = false;
    unsigned int is_bdry_bdry_region = 0;


    for(unsigned f = 0; f < /*ctrl::*/LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

        if (faceName == /*ctrl::*/LIST_OF_CTRL_FACES :: _face_with_extremes_index[f]) {

            is_facename_a_control_face = true;
            is_bdry_bdry_region = 1;

            for(unsigned t = 0; t < number_of_tangential_direction_components; t++) {

                          if (
                                  (x[ /*ctrl::*/LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName, number_of_tangential_direction_components)[t] ] >
                                            /*ctrl::*/LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][t][0] +  OFFSET_TO_INCLUDE_LINE &&
                                   x[ /*ctrl::*/LIST_OF_CTRL_FACES ::tangential_direction_to_Gamma_control(faceName, number_of_tangential_direction_components)[t] ] <
                                            /*ctrl::*/LIST_OF_CTRL_FACES ::_face_with_extremes_extremes_on_tang_surface[f][t][1] -  OFFSET_TO_INCLUDE_LINE
                                  )
                             ) { is_bdry_bdry_region *= 1; }

                          else { is_bdry_bdry_region *= 0; }

            }

            break;
        }

    }

    return std::pair< bool, unsigned int >/*std::make_pair*/( is_facename_a_control_face, (is_bdry_bdry_region == 1) );

}


static double linear_function_normal_direction(const double domain_length, const double gamma, const unsigned int face_index, const std::vector < double > & x){

    double value = 0;

    value = gamma * ( /*ctrl::boundary_conditions::*//*square_or_cube ::*/
    Boundary_Condition<LIST_OF_CTRL_FACES >::opposite_face_ctrl_or_state_value(face_index, domain_length) +  /*ctrl::*/
    LIST_OF_CTRL_FACES ::sign_function_for_delimiting_region(face_index) *
    x[/*ctrl::*/ LIST_OF_CTRL_FACES ::normal_direction_to_Gamma_control(face_index) ]
                    );

    return value;
}

}; // end class Boundary_Condition


}// end namespace square_or_cube
  
  
namespace boundary_control_between_extreme {
    
  namespace square {
  
  

template < class LIST_OF_CTRL_FACES >
 class Boundary_Condition_control_between_extreme: public square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES > {

    public:

      // template < class T >
 static bool ctrl_or_state_set_dirichlet_flags(const MultiLevelProblem * ml_prob,
                                        const int faceName,
                                        const std::vector < double > & x,
                                        bool &  dirichlet)  {

             const unsigned int number_of_tangential_direction_components = LIST_OF_CTRL_FACES ::  _num_of_tang_components_per_face ;
             std::pair< bool, bool > boundary_and_control_flags_region = square_or_cube:: Boundary_Condition<LIST_OF_CTRL_FACES >::is_boundary_of_boundary_region(faceName, x, number_of_tangential_direction_components);

             /*static*/ bool is_facename_a_control_face = boundary_and_control_flags_region.first;
             /*static*/ bool is_bndry_cntrl_region      = boundary_and_control_flags_region.second;

             if  (is_facename_a_control_face == false) { dirichlet = true; }
        else if  ( (is_facename_a_control_face == true) && !(is_bndry_cntrl_region) ) { dirichlet = true; }

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

    const unsigned int number_of_tangential_direction_components = LIST_OF_CTRL_FACES ::  _num_of_tang_components_per_face ;
    std::pair< bool, bool > boundary_and_control_flags_region = square_or_cube:: Boundary_Condition<LIST_OF_CTRL_FACES >:: is_boundary_of_boundary_region(faceName, x, number_of_tangential_direction_components);

    /*static*/ bool  is_facename_a_control_face = boundary_and_control_flags_region.first;
    /*static*/ bool  is_bndry_cntrl_region      = boundary_and_control_flags_region.second;



     const unsigned face_for_control_principal = /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_index[0];

    const double domain_length = 1.;

      const double gamma = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::_scale_factor_max_value_on_line;

        if (is_facename_a_control_face) {

              if(faceName == face_for_control_principal) {  value = 0.; }

              else if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){
                       if ( !(is_bndry_cntrl_region) )   { value =  gamma * domain_length; }
                       else {value = 0.; }
              }

              else {
                     if ( !(is_bndry_cntrl_region) ) {
                         value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,face_for_control_principal, x);
                         }
                     else {value = 0.; }
              }
        }
        else if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)) { value =  gamma * domain_length; }
        else      { value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,face_for_control_principal, x);}

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

    const double gamma = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::_scale_factor_max_value_on_line;


    if(faceName == face_for_control) {  value = 0.; }

    else if(faceName == adjacent_face_for_control) {  value = 0.; }

    else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control)){
        value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,adjacent_face_for_control, x);
           }
    else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){
        value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,face_for_control, x);
           }

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

     const unsigned int number_of_tangential_direction_components = LIST_OF_CTRL_FACES ::  _num_of_tang_components_per_face ;
     std::pair< bool, bool > boundary_and_control_flags_region =  square_or_cube:: Boundary_Condition<LIST_OF_CTRL_FACES >::is_boundary_of_boundary_region(faceName, x, number_of_tangential_direction_components);

     /*static*/ bool is_facename_a_control_face = boundary_and_control_flags_region.first;
     /*static*/ bool is_bndry_cntrl_region      = boundary_and_control_flags_region.second;

     const unsigned face_for_control_principal = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[0];
     const unsigned adjacent_face_for_control = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[1];

     const double domain_length = 1.;

     const double gamma = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::_scale_factor_max_value_on_line;


     if (is_facename_a_control_face)     {  value = 0.; }

     else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){
         value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,adjacent_face_for_control, x);
             }

     else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){
         value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,face_for_control_principal, x);
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

     const unsigned int number_of_tangential_direction_components = LIST_OF_CTRL_FACES ::  _num_of_tang_components_per_face ;
     std::pair< bool, bool > boundary_and_control_flags_region =  square_or_cube:: Boundary_Condition<LIST_OF_CTRL_FACES >::is_boundary_of_boundary_region(faceName, x, number_of_tangential_direction_components);

     /*static*/ bool is_facename_a_control_face = boundary_and_control_flags_region.first;
     /*static*/ bool is_bndry_cntrl_region      = boundary_and_control_flags_region.second;

     assert( /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size == 3 );

     const unsigned face_for_control_principal = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[0];
     const unsigned adjacent_face_for_control = /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[1];

     const double domain_length = 1.;

     const double gamma = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::_scale_factor_max_value_on_line;



         if (is_facename_a_control_face) {

                if ( (faceName == face_for_control_principal) || (faceName == adjacent_face_for_control) )     {  value = 0.; }

                else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){

                        if ( !(is_bndry_cntrl_region) ) {
                                  value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,adjacent_face_for_control, x);
                                  }
                        else {value = 0.; }
                }
                else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){
                        if ( !(is_bndry_cntrl_region ) ) {
                               value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,face_for_control_principal, x);
                                }
                               else {value = 0.; }
                }

         }
         else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){
             value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,adjacent_face_for_control, x);
              }

         else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){
             value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,face_for_control_principal, x);
               }

     value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

     return value;
}

 };



}



        
  namespace cube {
      
        
template < class LIST_OF_CTRL_FACES >
 class Multiple_controls_and_homogeneous_boundary_conditions : public femus:: ctrl :: boundary_control_between_extreme :: square :: Multiple_controls_and_homogeneous_boundary_conditions< LIST_OF_CTRL_FACES > {

 };
 
        
        
      }
      
  }
  
  
  
 namespace boundary_control_full_face {
     
 namespace square {
  


template < class LIST_OF_CTRL_FACES >
 class Boundary_Condition_control_full_face: public square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES > {

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

      const double gamma = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::_scale_factor_max_value_on_line;


	  for(unsigned f = 0; f < /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

        if (faceName == /*ctrl:: */ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]) {value = 0.;}

        else if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)) { value =  gamma * domain_length; }
        else      {
            value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,face_for_control_principal, x);
           }
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

    const double gamma = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::_scale_factor_max_value_on_line;


    if(faceName == face_for_control) {  value = 0.; }
    else if(faceName == adjacent_face_for_control) {  value = 0.; }
    else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control)){
        value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,adjacent_face_for_control, x);
         }
    else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){
        value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,face_for_control, x);
        }

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

     const double gamma = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::_scale_factor_max_value_on_line;


     for(unsigned f = 0; f < /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

        if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f])     {  value = 0.; }
        else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){
            value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,adjacent_face_for_control, x);
            }
        else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){
            value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,face_for_control_principal, x);
            }
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

     const double gamma = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::_scale_factor_max_value_on_line;


     for(unsigned f = 0; f < /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index_size; f++) {

         if (faceName == /*ctrl::*/ LIST_OF_CTRL_FACES ::_face_with_extremes_index[f]) { value = 0.; }

         else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(face_for_control_principal)){
             value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,adjacent_face_for_control, x);
              }

         else if(faceName == /*ctrl::*/ LIST_OF_CTRL_FACES :: opposite_face(adjacent_face_for_control)){
             value = square_or_cube :: Boundary_Condition<LIST_OF_CTRL_FACES >::linear_function_normal_direction(domain_length, gamma,face_for_control_principal, x);
              }

     }

     value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

     return value;
}






 };





}

        
   namespace cube {
     
    }  
      
  }  
  
  
  

   }

}


#endif
