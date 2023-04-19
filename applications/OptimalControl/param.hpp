#ifndef OPT_CONTROL_PARAMETERS
#define OPT_CONTROL_PARAMETERS





//*******************************************************************************************
//*********************** Domain and Mesh Dependent, Parameters - BEGIN *****************************************
//*******************************************************************************************

//*********************** Control, Gamma_c specification - BEGIN  *******************************************************
  /* Rectangular/Hexahedral domain:  1-2 x coords, 3-4 y coords, 5-6 z coords */
  /* L-shaped domain (2d):  1-2 x coords, 3-4 y coords, 5 indent between 1 and 2, 6 indent between 3 and 4 */
#define FACE_FOR_CONTROL        2



#define SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN 0.25  /* 0. */
#define SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX 0.75  /* 1. */

//*********************** Control, Gamma_c specification - END *******************************************************

//*******************************************************************************************
//*********************** Domain and Mesh Dependent, Parameters - END *****************************************
//*******************************************************************************************




#include "opt_systems_boundary_control_eqn_sobolev_fractional.hpp"









//*******************************************************************************************
//*********************** Domain and Mesh Independent, Parameters - BEGIN *****************************************
//*******************************************************************************************

namespace femus {




   
namespace ctrl {
    

//*********************** 
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = false;
//*********************** 

  
//*********************** Mesh - BEGIN *****************************************

//*********************** Mesh, Number of refinements - BEGIN *****************************************
// #define N_UNIFORM_LEVELS 7// for 2D applications
#define N_UNIFORM_LEVELS 3 // for 3D bdry application

#define N_ERASED_LEVELS   N_UNIFORM_LEVELS - 1

  //with 0 it only works in serial, you must put 2 to make it work in parallel...: that's because when you fetch the dofs from _topology you get the wrong indices
//*********************** Mesh, Number of refinements - END *****************************************


//*********************** Mesh, offset for point inclusion - BEGIN *******************************************************
#define OFFSET_TO_INCLUDE_LINE  1.e-5  
//*********************** Mesh, offset for point inclusion - END *******************************************************

//*********************** Mesh - END *****************************************

  
//*********************** Control, cost functional without regularization - BEGIN *******************************************************
#define COST_FUNCTIONAL_TYPE  0 /*1  *//*------[0: target ; 1: gradient]---------*/

#define COST_FUNCTIONAL_COEFF 1 
//*********************** Control, cost functional without regularization - END *******************************************************

//*********************** Control, cost functional, regularization - BEGIN *******************************************************
// for pure boundary approaches
#define ALPHA_CTRL_BDRY 1.e-7
#define BETA_CTRL_BDRY   ALPHA_CTRL_BDRY

// for lifting approaches (both internal and external)
#define ALPHA_CTRL_VOL ALPHA_CTRL_BDRY 
#define BETA_CTRL_VOL ALPHA_CTRL_VOL
//*********************** Control, cost functional, regularization - END *******************************************************

  

//*********************** Control, cost functional, regularization (almost), Boundary, Fractional or Integer - BEGIN  *******************************************************
#define CONTROL_FORMULATION   0  /* 0 = pure boundary; 1 = lifting internal; 2 = lifting external */  

     #define OP_H1       1 /*0*/

//***** Operator-related - BEGIN ****************** 
#if CONTROL_FORMULATION == 0

  
/* 0: integer norm, 1: fractional norm */
//   #define IS_CTRL_FRACTIONAL_SOBOLEV  0 // INTEGER NORM
  #define IS_CTRL_FRACTIONAL_SOBOLEV  1  // FRACTIONAL NORM


  //for L2 norm test
#if IS_CTRL_FRACTIONAL_SOBOLEV == 0
     #define OP_L2       1

     #define OP_Hhalf    0

     #define UNBOUNDED    0

      #define USE_Cns     0
#endif

#if IS_CTRL_FRACTIONAL_SOBOLEV == 1

// #define NORM_GIR_RAV 0  /* Leave it at 0 --> 0 = Gagliardo seminorm   */
#define NORM_GIR_RAV  1    /* Leave it at 0 --> 1 = Girault Raviart norm */

#if NORM_GIR_RAV == 0

  #define OP_L2       0 /*1*/
  
  #define OP_Hhalf    1

  #define UNBOUNDED   1

  #define USE_Cns     1

#elif NORM_GIR_RAV == 1 

  #define OP_L2       1
  
  #define OP_Hhalf    1 /*0*/

  #define UNBOUNDED   0

  #define USE_Cns     0
#endif
  
#endif
  
  #endif

  
//***** Operator-related - END ****************** 

  
//***** Quadrature-related - BEGIN ****************** 
// for integrations in the same element
#define Nsplit 0
#define Quadrature_split_index  0

//for semi-analytical integration in the unbounded domain
#define N_DIV_FACE_OF_FACE_FOR_UNBOUNDED_INTEGRAL  10

#define QRULE_I   0
#define QRULE_J   1
#define QRULE_K   QRULE_I
//***** Quadrature-related - END ****************** 

  


//******** Penalties for equations - BEGIN ******************************
#define PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY           1.e50      
#define PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY           0.      
#define PENALTY_DIRICHLET_BC_U_EQUAL_Q_BOUNDARY           1.e10         // penalty for u = q
//******** Penalties for equations - END ******************************


//*********************** Control, cost functional, regularization (almost), Boundary, Fractional or Integer - END  *******************************************************

///********domain dependent -BEGIN
//*********************** Control, cost functional and equation, Lifting External - BEGIN  *******************************************************
#define GROUP_INTERNAL  12
#define GROUP_EXTERNAL  13
//*********************** Control, cost functional and equation, Lifting External - END *******************************************************
///********domain dependent -END


} //end namespace 







} //end namespace femus

  
  
//*******************************************************************************************
//*********************** Domain and Mesh Independent, Parameters - END *****************************************
//*******************************************************************************************






//*******************************************************************************************
//*********************** Domain and Mesh Dependent: Square or Cube - BEGIN *****************************************
//*******************************************************************************************


#include "square_or_cube_00_mesh_files_for_control.hpp"
#include "square_or_cube_01_control_faces.hpp"
#include "square_or_cube_01b_control_faces_domain_elements.hpp"
#include "square_or_cube_02_boundary_conditions.hpp"
#include "square_or_cube_03_cost_functional_without_regularization.hpp"
#include "square_or_cube_04_cost_functional_regularization.hpp"
#include "square_or_cube_05_opt_system_inequalities.hpp"



namespace femus {



namespace ctrl {


//*********************** Gamma_c, Choice of Domain,
//*********************** Gamma_c, Choice of Mesh File,
//******************************** Choice of List of control faces, ***********************
//******************************** Choice of corresponding boundary conditions and *****************************************
//******************************** Choice of Cost Functional-related stuff - BEGIN *********************** 









#define  TYPE_OF_BOUNDARY_CONTROL   boundary_control_between_extreme
// #define  TYPE_OF_BOUNDARY_CONTROL   boundary_control_full_face




//------------------------------------ Square: BEGIN ------------------------------------
// // // #define  DOMAIN_NAMESPACE_CHILD   square

// // //    const std::string mesh_input = ctrl:: DOMAIN_NAMESPACE_CHILD ::mesh::_2d_square_1x1;

//------------------------------------ single: BEGIN ------------------------------------

// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_One
// // // // // // // // // // // #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Two
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Three
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Four

// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Single_control_in_front_linear
// // // // // // // // // // // // // #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Multiple_controls_and_homogeneous_boundary_conditions
// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Multiple_controls_in_front_constant

//------------------------------------ single: END ------------------------------------


//------------------------------------ double: BEGIN ------------------------------------

// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_One_Three
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_One_Four
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Two_One
// // // #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Two_Three
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Two_Four
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Three_One
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Three_Two
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Three_Four
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Four_One
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Four_Two
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Four_Three

//  #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS   TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Double_controls_adjacent_in_front_linear
// // // #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Multiple_controls_and_homogeneous_boundary_conditions
// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Multiple_controls_in_front_constant

//------------------------------------ double: END ------------------------------------


//------------------------------------ Triple: BEGIN ------------------------------------

// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_One_Three_Two
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_One_Four_Two
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_One_Three_Four
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_One_Four_Three
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Two_Three_One
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Two_Four_One
// // // // // // // #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Two_Three_Four
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Two_Four_Three
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Three_One_Four
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Three_Two_Four
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Three_One_Two
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Three_Two_One
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Four_One_Three
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Four_Two_Three
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Four_One_Two

// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Triple_controls_adjacent_in_front_linear
// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Multiple_controls_and_homogeneous_boundary_conditions
// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Multiple_controls_in_front_constant

//------------------------------------ Triple: END ------------------------------------


//------------------------------------ Quadruple: BEGIN ------------------------------------

// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES      TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_One_Three_Two_Four

// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Multiple_controls_and_homogeneous_boundary_conditions
// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Multiple_controls_in_front_constant

//------------------------------------ Quadruple: END ------------------------------------


//------------------------------------ Square: END ------------------------------------






//------------------------------------ Cube: BEGIN ------------------------------------

#define  DOMAIN_NAMESPACE_CHILD   cube

//       const std::string mesh_input = ctrl::DOMAIN_NAMESPACE_CHILD:: mesh::_3d_cube_single_face_control_2_old_coarser;
//       const std::string mesh_input = ctrl::DOMAIN_NAMESPACE_CHILD:: mesh::_3d_cube_single_face_control_2_old_coarser_between;
      const std::string mesh_input = ctrl::DOMAIN_NAMESPACE_CHILD:: mesh::_3d_cube_all_faces_between_between;

// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_One
#define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Two
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Three
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Four
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Five
// #define  GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: List_of_Gamma_control_faces_Six

#define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    TYPE_OF_BOUNDARY_CONTROL :: DOMAIN_NAMESPACE_CHILD :: Multiple_controls_and_homogeneous_boundary_conditions

//------------------------------------ Cube: END ------------------------------------









//------------------------------------ Cost functional without regularization: BEGIN ------------------------------------
#define  COST_FUNCTIONAL_WITHOUT_REG    DOMAIN_NAMESPACE_CHILD :: cost_functional_without_regularization
//------------------------------------ Cost functional without regularization: END ------------------------------------

      
      
//*********************** Gamma_c, Choice of List of control faces and *********************** 
//******************************** Choice of possible corresponding boundary conditions - END *****************************************








}


}


//*******************************************************************************************
//*********************** Domain and Mesh Dependent: Square or Cube - END *****************************************
//*******************************************************************************************





  
#endif
