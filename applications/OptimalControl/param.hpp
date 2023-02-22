#ifndef OPT_CONTROL_PARAMETERS
#define OPT_CONTROL_PARAMETERS



#include "Mesh.hpp"
#include "CurrentElem.hpp"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "MED_IO.hpp"
#include "paral.hpp"

#include "fractional_functions.hpp"


#include <functional>

#include <boost/mpi.hpp>


// using namespace femus;



namespace femus {




    
//*******************************************************************************************
//*********************** Domain and Mesh Independent - BEGIN *****************************************
//*******************************************************************************************






 
void el_dofs_quantities_vol(const Solution*                sol,
                        const Mesh * msh,
                        const unsigned int iel,
                        const    vector < unsigned > & SolFEType,
                        vector < unsigned > & Sol_n_el_dofs 
     ) {
    
        for (unsigned  k = 0; k < Sol_n_el_dofs.size(); k++) {
            unsigned  ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
            Sol_n_el_dofs[k] = ndofs_unk;
        }   
    
}


 
void el_dofs_unknowns_vol(const Solution*                sol,
                      const Mesh * msh,
                      const  LinearEquationSolver* pdeSys,
                      const unsigned int iel,
                      const    vector < unsigned > & SolFEType,
                      const vector < unsigned > & SolIndex,
                      const vector < unsigned > & SolPdeIndex,
                      vector < unsigned > & Sol_n_el_dofs, 
                      vector < vector < double > > & sol_eldofs,  
                      vector < vector < int > > & L2G_dofmap ) {
    
    assert(Sol_n_el_dofs.size() == sol_eldofs.size());
    
        //all vars###################################################################
        for (unsigned  k = 0; k < Sol_n_el_dofs.size(); k++) {
            unsigned  ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
            Sol_n_el_dofs[k] = ndofs_unk;
            sol_eldofs[k].resize(ndofs_unk);
            L2G_dofmap[k].resize(ndofs_unk);
            for (unsigned i = 0; i < ndofs_unk; i++) {
                unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);
                sol_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);
                L2G_dofmap[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);
            }
        }
        //all vars###################################################################

}




void  print_global_residual_jacobian(const bool print_algebra_global,
                                     const MultiLevelProblem& ml_prob,
                                     const NonLinearImplicitSystem * mlPdeSys,
                                       LinearEquationSolver* pdeSys,
                                     NumericVector*           RES,
                                     SparseMatrix*             JAC,
                                     const unsigned iproc,
                                     const bool assembleMatrix)  {


  if (print_algebra_global) {
    const unsigned nonlin_iter = mlPdeSys->GetNonlinearIt();
    assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, JAC, nonlin_iter);
//     assemble_jacobian< double, double >::print_global_residual(ml_prob, RES,  mlPdeSys->GetNonlinearIt());

    RES->close();
    std::ostringstream res_out; res_out << ml_prob.GetFilesHandler()->GetOutputPath() << "./" << "res_" << nonlin_iter  << ".txt";
    pdeSys->print_with_structure_matlab_friendly(iproc, res_out.str().c_str(), RES);

     } 
  
  }


  
    
  
  bool check_if_same_elem(const unsigned iel, const unsigned jel) {
      
   return (iel == jel);
      
  }
  
  
  
  bool check_if_same_elem_bdry(const unsigned iel, const unsigned jel, const unsigned iface, const unsigned jface) {

   return (iel == jel && iface == jface);
      
  }


  
  void  set_dense_pattern_for_unknowns(NonLinearImplicitSystem/*WithPrimalDualActiveSetMethod*/  & system, const std::vector < Unknown > unknowns)  {
  
///     system.init();   ///@todo Understand why I cannot put this here but it has to be in the main(), there must be some objects that get destroyed, passing the reference is not enough

  const MultiLevelProblem &  ml_prob = system.GetMLProb();
  const MultiLevelMesh *  ml_mesh = ml_prob.GetMLMesh();
  
  unsigned n_levels = ml_mesh->GetNumberOfLevels();
  std::ostringstream sp_out_base; sp_out_base << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "sp_";
  system._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base.str(), "on");
  system._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base.str(), "off");

  
  const Mesh* msh = ml_mesh->GetLevel(n_levels - 1);
  unsigned nprocs = msh->n_processors();
  unsigned iproc = msh->processor_id();

  
    for(int ivar = 0; ivar < unknowns.size(); ivar++) {
  
        if (unknowns[ivar]._is_sparse == false ) { 
        
  unsigned variable_index = system.GetSolPdeIndex(unknowns[ivar]._name.c_str());
  
  
  
  unsigned n_dofs_var_all_procs = 0;
  for (int ip = 0; ip < nprocs; ip++) {
     n_dofs_var_all_procs += system._LinSolver[n_levels - 1]->KKoffset[variable_index + 1][ip] - system._LinSolver[n_levels - 1]->KKoffset[variable_index][ip];
  // how does this depend on the number of levels and the number of processors? 
  // For the processors I summed over them and it seems to work fine
  // For the levels... should I pick the coarsest level instead of the finest one, or is it the same?
   } 

  unsigned n_vars = system.GetSolPdeIndex().size();   //assume all variables are dense: we should sum 3 sparse + 1 dense... or better n_components_ctrl * dense and the rest sparse
  //check that the dofs are picked correctly, it doesn't seem so 
  
  system.SetSparsityPatternMinimumSize (n_dofs_var_all_procs * n_vars, unknowns[ivar]._name);  ///@todo this is like AddSolution: it increases a vector
        }
    }
    
    
    
  return; 
  
  }


//*******************************************************************************************
//*********************** Domain and Mesh Independent - END *****************************************
//*******************************************************************************************




   
namespace ctrl {
    

//*******************************************************************************************
//*********************** Domain and Mesh Independent, Parameters - BEGIN *****************************************
//*******************************************************************************************

//*********************** 
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = false;
//*********************** 

  
//*********************** Mesh - BEGIN *****************************************

//*********************** Mesh, Number of refinements - BEGIN *****************************************
#define N_UNIFORM_LEVELS 5
#define N_ERASED_LEVELS   N_UNIFORM_LEVELS - 1

#define FE_DOMAIN  2 //with 0 it only works in serial, you must put 2 to make it work in parallel...: that's because when you fetch the dofs from _topology you get the wrong indices
//*********************** Mesh, Number of refinements - END *****************************************


//*********************** Mesh, offset for point inclusion - BEGIN *******************************************************
#define OFFSET_TO_INCLUDE_LINE  1.e-5  
//*********************** Mesh, offset for point inclusion - END *******************************************************

//*********************** Mesh - END *****************************************

  
//*********************** Control, cost functional - BEGIN *******************************************************
#define COST_FUNCTIONAL_TYPE  0 /*1  *//*------[0: target ; 1: gradient]---------*/

#define COST_FUNCTIONAL_COEFF 1 

// for pure boundary approaches
#define ALPHA_CTRL_BDRY 1.e-6
#define BETA_CTRL_BDRY   ALPHA_CTRL_BDRY

// for lifting approaches (both internal and external)
#define ALPHA_CTRL_VOL ALPHA_CTRL_BDRY 
#define BETA_CTRL_VOL ALPHA_CTRL_VOL
//*********************** Control, cost functional - END *******************************************************

  

//*********************** Control, Boundary, Fractional or Integer - BEGIN  *******************************************************


//***** Operator-related - BEGIN ****************** 
#define IS_CTRL_FRACTIONAL_SOBOLEV  1      /* 0: integer norm, 1: fractional norm */


#define RHS_ONE             0.

#define KEEP_ADJOINT_PUSH   1

#define S_FRAC 0.5       /* Order of fractional derivative */

#define NORM_GIR_RAV  0  /* Leave it at 0 */

#if NORM_GIR_RAV == 0

  #define OP_L2       0
  #define OP_H1       0
  #define OP_Hhalf    1

  #define UNBOUNDED   1

  #define USE_Cns     1

#elif NORM_GIR_RAV == 1 

  #define OP_L2       1
  #define OP_H1       0
  #define OP_Hhalf    1

  #define UNBOUNDED   0

  #define USE_Cns     0
#endif
//***** Operator-related - END ****************** 

  
//***** Quadrature-related - BEGIN ****************** 
// for integrations in the same element
#define Nsplit 0
#define Quadrature_split_index  0

//for semi-analytical integration in the unbounded domain
#define N_DIV_UNBOUNDED  10

#define QRULE_I   0
#define QRULE_J   1
#define QRULE_K   QRULE_I
//***** Quadrature-related - END ****************** 

  
//***** Implementation-related: where are L2 and H1 norms implemented - BEGIN ****************** 
#define IS_BLOCK_DCTRL_CTRL_INSIDE_MAIN_BIG_ASSEMBLY   0  // 1 internal routine; 0 external routine
//***** Implementation-related: where are L2 and H1 norms implemented - END ****************** 


//******** Penalties for equations - BEGIN ******************************
#define PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY           1.e50      
#define PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY           0.      
#define PENALTY_DIRICHLET_BC_U_EQUAL_Q_BOUNDARY           1.e10         // penalty for u = q
//******** Penalties for equations - END ******************************


//*********************** Domain Dependent - BEGIN *****************************************

//***** How to identify boundary of boundary in the 3D case - BEGIN ****************** 
#define  NODE_BASED_BDRY_BDRY   "node_based_bdry_bdry_flag"
//***** How to identify boundary of boundary in the 3D case - END ****************** 

//*********************** Domain Dependent - END *****************************************



//*********************** Control, Boundary, Fractional or Integer - END  *******************************************************


//*********************** Inequality - BEGIN *******************************************************
#define  INEQ_FLAG 1
#define  C_COMPL 1.
//*********************** Inequality - END *******************************************************

  
//*******************************************************************************************
//*********************** Domain and Mesh Independent, Parameters - END *****************************************
//*******************************************************************************************




  
//*******************************************************************************************
//*********************** Domain and Mesh Dependent, Parameters - BEGIN *****************************************
//*******************************************************************************************

//*********************** Control, Gamma_c specification - BEGIN  *******************************************************
  /* Rectangular/Hexahedral domain:  1-2 x coords, 3-4 y coords, 5-6 z coords */
  /* L-shaped domain (2d):  1-2 x coords, 3-4 y coords, 5 indent between 1 and 2, 6 indent between 3 and 4 */
#define FACE_FOR_CONTROL        3



#define GAMMA_CONTROL_LOWER 0.25  /*0.*/
#define GAMMA_CONTROL_UPPER 0.75  /*1.*/
//***** Domain-related ****************** 
#define EX_1        GAMMA_CONTROL_LOWER
#define EX_2        GAMMA_CONTROL_UPPER
#define EY_1        0.
#define EY_2        1.   ///@todo  see here

#define DOMAIN_EX_1 0
#define DOMAIN_EX_2 1
//**************************************
//*********************** Control, Gamma_c specification - END *******************************************************


//*********************** Control, cost functional, target region - BEGIN *******************************************************
  /* Rectangular/Hexahedral domain:  1-2 x coords, 3-4 y coords, 5-6 z coords */
  /* L-shaped domain (2d):  1-2 x coords, 3-4 y coords, 5 indent between 1 and 2, 6 indent between 3 and 4 */
#define FACE_FOR_TARGET         4

#define  TARGET_LINE_ORTHOGONAL_DISTANCE_FROM_FACE_ATTACHED_TO_TARGET_REG  0.5
//*********************** Control, cost functional, target region - END *******************************************************



//*********************** Control, boundary - BEGIN *******************************************************
#define BOUNDARY_ORTHOGONAL_DISTANCE_FROM_GAMMA_C  1.   //how far it goes orthogonally to the Control piece of the Boundary 

//*********************** Control, boundary - END *******************************************************



//*********************** Control, Lifting internal - BEGIN *******************************************************
#define LIFTING_INTERNAL_ORTHOGONAL_DISTANCE_FROM_GAMMA_C  1.   //how far it goes orthogonally to the Control piece of the Boundary 
#define LIFTING_INTERNAL_WIDTH_LOWER  0. /*GAMMA_CONTROL_LOWER*/
#define LIFTING_INTERNAL_WIDTH_UPPER  1. /*GAMMA_CONTROL_UPPER*/

//******** Penalties for equations - BEGIN ******************************
#define PENALTY_OUTSIDE_CONTROL_DOMAIN_LIFTING_INTERNAL   1.e20         // penalty for zero control outside
//******** Penalties for equations - END ******************************

//*********************** Control, Lifting internal - END *******************************************************



//*******************************************************************************************
//*********************** Domain and Mesh Dependent, Parameters - END *****************************************
//*******************************************************************************************






//*******************************************************************************************
//*********************** Domain and Mesh Dependent - BEGIN *****************************************
//*******************************************************************************************
    
 namespace Gamma_control_list_of_faces_with_extremes {
 
//*********************** Mesh dependent, Gamma_c, List of Faces - BEGIN *****************************************
    
//     class face_with_extremes {
//         
//         
//     public:
        
     static const unsigned face_with_extremes_index_size = /*1*/ /*2*/ 3 ;
    
     static const unsigned face_with_extremes_index[ ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size ] = {
         FACE_FOR_CONTROL
       , 1
       , 2
  //     , FACE_FOR_CONTROL + 2
    };
     
     static const bool     face_with_extremes_extract_subface[ ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size ] = {
         true
       , true
       , true
  //     , true
    };
        
     static const double   face_with_extremes_extremes[ ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size ][2] = {
       { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }
     , { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }
     , { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }
  //   , { GAMMA_CONTROL_LOWER, GAMMA_CONTROL_UPPER }
    };
        
     
//     };


#define  FACE_CONTROL_ADJACENT   2/*3*/// ***for "Gamma_c_single_control_in_front_linear"***


// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    Gamma_c_double_adjacent
//  #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    Gamma_c_single


// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    Multiple_controls_and_homogenous_boundary_conditions

// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    Multiple_controls_in_front_constant

// #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    Gamma_c_single_control_in_front_linear

//  #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    Gamma_c_double_adjacent_control_in_front_linear

 #define NAMESPACE_FOR_GAMMA_C_BOUNDARY_CONDITIONS    Gamma_c_triple_adjacent_control_in_front_linear

//*********************** Mesh dependent, Gamma_c, List of Faces - END *****************************************

 }


namespace mesh {

//   const std::string input = "Study7_GroupofNode_Face2.med";

  const std::string input = "parametric_square_1x1.med";

//      const std::string input = "Mesh_3_groups_with_bdry_nodes_coarser.med";

//   std::string input_file = "parametric_square_1x1.med";
//   std::string input_file = "parametric_square_1x2.med";
//   std::string input_file = "parametric_square_2x2.med";
//   std::string input_file = "parametric_square_4x5.med";
//   std::string input_file = "Mesh_3_groups_with_bdry_nodes.med";

     
  std::string file_with_prefix()  {
     

     const std::string input_file = ctrl::mesh::input;
     
     const std::string mesh_location = "../../../";
     
      std::ostringstream mystream; mystream << mesh_location /*"./"*/ << DEFAULT_INPUTDIR << "/" << input_file;
      const std::string infile = mystream.str();
      
      return infile;

   }
     
     
}
      

namespace boundary_conditions_or_cost_functional {

    
const  int sign_function_for_delimiting_region(const unsigned int face_index) {
    
   int  target_line_sign;
  
        if (face_index == 1 || face_index == 3 || face_index == 5) { target_line_sign = 1;  }
   else if (face_index == 2 || face_index == 4 || face_index == 6) { target_line_sign = -1; }
   
   return target_line_sign;
   
}


}



namespace cost_functional {
 
    
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

  const int  target_line_sign = boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(FACE_FOR_TARGET);
  
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

     
    
    
    
}  




namespace boundary_conditions {


    
 //direction of the line that contains \Gamma_c    
const unsigned int normal_direction_to_Gamma_control(const unsigned int face_index) {
    
    unsigned int axis_dir;
    
        if (face_index == 1 || face_index == 2) { axis_dir = 0; }
   else if (face_index == 3 || face_index == 4) { axis_dir = 1; }
   else if (face_index == 5 || face_index == 6) { axis_dir = 2; }

    return axis_dir;
    
}




const unsigned int tangential_direction_to_Gamma_control(const unsigned int face_index) {
    
    unsigned int axis_dir;
    
        if (face_index == 1 || face_index == 2) { axis_dir = 1; }
   else if (face_index == 3 || face_index == 4) { axis_dir = 0; }
   else if (face_index == 5 || face_index == 6) { /*abort();*/ axis_dir = 1; }  ///@todo this depends on the mesh file 

    return axis_dir;
    
}




double opposite_face_ctrl_or_state_value(const unsigned int face_index, const double domain_length) {

   double opposite_value = 0.;

        if (face_index == 1 || face_index == 3 || face_index == 5) { opposite_value = 0.; }
   else if (face_index == 2 || face_index == 4 || face_index == 6) { opposite_value = domain_length; }

    return opposite_value;

}



const unsigned int opposite_face(const unsigned int face_index) {

   unsigned int face_opposite = 0;

        if (face_index == 1 || face_index == 3 || face_index == 5) { face_opposite = face_index + 1; }
   else if (face_index == 2 || face_index == 4 || face_index == 6) { face_opposite = face_index - 1; }

    return face_opposite;

}


namespace Gamma_c_single {
    
    

 double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob, 
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {
     
     assert( ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size == 1 );

    const unsigned face_for_control = ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[0];

    const double domain_length = 1.;

      const double gamma = 5.;

        if (faceName == face_for_control)     {  value = 0.; }
   else if (faceName == ctrl::boundary_conditions::opposite_face(face_for_control)) { value =  gamma * domain_length; }
   else                                       { value = gamma * ( boundary_conditions::opposite_face_ctrl_or_state_value(face_for_control, domain_length) + boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(face_for_control) *  x[ boundary_conditions::normal_direction_to_Gamma_control(face_for_control) ] ); }

   value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;
   
   return value;
}


 bool ctrl_or_state_set_dirichlet_flags(const MultiLevelProblem * ml_prob, 
                                        const int faceName,
                                        const std::vector < double > & x,
                                        bool &  dirichlet)  {
     

     assert( ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size == 1 );
    const unsigned face_for_control = ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[0];

     
     if (faceName == face_for_control) {
        if ( !(x[ boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] > ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[0][0] + 1.e-5 &&
               x[ boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] < ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[0][1] - 1.e-5) ) {
                dirichlet = true;
           }
    }
    else {
          dirichlet = true;
    }

    return dirichlet;
}



}

namespace Gamma_c_double_adjacent {
    
    

    
 double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob, 
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

     if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

     assert( ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size == 2 );

    const double domain_length = 1.;

      const double gamma = 5.;
      
      
	  for(unsigned f = 0; f < ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

        if (faceName == ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f])     {  value = 0.; }
   else if (faceName == ctrl::boundary_conditions::opposite_face(ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f])) { value = gamma * 
       ( boundary_conditions::opposite_face_ctrl_or_state_value(ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f], domain_length) + boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) *  x[ boundary_conditions::tangential_direction_to_Gamma_control(ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) ] );  }

      }
   
   value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;
   
   return value;
}


 bool ctrl_or_state_set_dirichlet_flags(const MultiLevelProblem * ml_prob, 
                                        const int faceName,
                                        const std::vector < double > & x,
                                        bool &  dirichlet)  {

     if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();
     
     assert( ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size == 2 );
     
     
     bool  is_facename_a_control_face = false;
     
     for(unsigned f = 0; f < ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

     if (faceName == ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) { is_facename_a_control_face = true; break; }
          
     }

     if  (is_facename_a_control_face == false) { dirichlet = true; }

      
     for(unsigned f = 0; f < ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

     if (faceName == ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) {
         
        if ( !(x[ boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] > ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][0] + 1.e-5 &&
               x[ boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] < ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][1] - 1.e-5) ) {
                dirichlet = true;
           }
        }
      }    
    

    return dirichlet;
}



}



namespace Multiple_controls_and_homogenous_boundary_conditions {



 double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

//      assert( face_with_extremes_index_size == 1 );

//     const unsigned face_for_control = face_with_extremes_index[0];

//     const double domain_length = 1.;

//       const double gamma = 5.;

//         if (faceName == face_for_control)     {  value = 0.; }
//    else if (faceName == ctrl::opposite_face(face_for_control)) { value =  gamma * domain_length; }
//    else                                       { value = gamma * (
//                                                 ctrl::opposite_face_ctrl_or_state_value(face_for_control, domain_length) + ctrl::sign_function_for_delimiting_region(face_for_control) *  x[      ctrl::normal_direction_to_Gamma_control(face_for_control) ] ); }

   value = 0 + PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

   return value;
}


 bool ctrl_or_state_set_dirichlet_flags(const MultiLevelProblem * ml_prob,
                                        const int faceName,
                                        const std::vector < double > & x,
                                        bool &  dirichlet)  {


//      assert( face_with_extremes_index_size == 1 );
//     const unsigned face_for_control = face_with_extremes_index[0];

      bool  is_facename_a_control_face = false;
      for(unsigned f = 0; f < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

           if (faceName == femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) { is_facename_a_control_face = true; break; }

           }

           if  (is_facename_a_control_face == false) { dirichlet = true; }


           for(unsigned f = 0; f < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

           if (faceName == femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) {

              if ( !(x[ femus::ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] > femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][0] + 1.e-5 &&
                     x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][1] - 1.e-5) ) {
                      dirichlet = true;
                 }
              }
            }


          return dirichlet;
}



}


namespace Multiple_controls_in_front_constant {


 double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

//      if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

//      assert( face_with_extremes_index_size == 2 );

     const unsigned face_for_control_principal = femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[0];

    const double domain_length = 1.;

      const double gamma = 5.;


	  for(unsigned f = 0; f < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

        if (faceName == femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) {
              if(faceName == face_for_control_principal) {  value = 0.; }
              else if (faceName == ctrl::boundary_conditions::opposite_face(face_for_control_principal)){
                     if ( !(x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] > femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][0] + 1.e-5 &&
                            x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][1] - 1.e-5) )   { value =  gamma * domain_length; }
                     else {value = 0.; }
              }
              else {
                     if ( !(x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] > femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][0] + 1.e-5 &&
                            x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][1] - 1.e-5) ) { value = gamma * (
                                      ctrl::boundary_conditions::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) + ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(face_for_control_principal) *  x[      ctrl::boundary_conditions::normal_direction_to_Gamma_control(face_for_control_principal) ] ); }
                     else {value = 0.; }
              }
        }
        else if (faceName == ctrl::boundary_conditions::opposite_face(face_for_control_principal)) { value =  gamma * domain_length; }
        else      { value = gamma * (
                                    ctrl::boundary_conditions::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) + ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(face_for_control_principal) *  x[      ctrl::boundary_conditions::normal_direction_to_Gamma_control(face_for_control_principal) ] ); }

      }

   value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

   return value;
}


 bool ctrl_or_state_set_dirichlet_flags(const MultiLevelProblem * ml_prob,
                                        const int faceName,
                                        const std::vector < double > & x,
                                        bool &  dirichlet)  {

 //    if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

  //   assert( face_with_extremes_index_size == 2 );


     bool  is_facename_a_control_face = false;

     for(unsigned f = 0; f < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

     if (faceName == femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) { is_facename_a_control_face = true; break; }

     }

     if  (is_facename_a_control_face == false) { dirichlet = true; }


     for(unsigned f = 0; f < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

     if (faceName == femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) {

        if ( !(x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] > femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][0] + 1.e-5 &&
               x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][1] - 1.e-5) ) {
                dirichlet = true;
           }
        }
      }


    return dirichlet;
}



}


namespace Gamma_c_single_control_in_front_linear {



 double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

    assert( femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size == 1 );

    const unsigned face_for_control = femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[0];
    const unsigned adjacent_face_for_control = FACE_CONTROL_ADJACENT; /*face_with_extremes_index[1]*/


    const double domain_length = 1.;

    const double gamma = 5.;


    if(faceName == face_for_control) {  value = 0.; }
    else if(faceName == adjacent_face_for_control) {  value = 0.; }
    else if(faceName == ctrl::boundary_conditions::opposite_face(face_for_control)){ value = gamma * (
             ctrl::boundary_conditions::opposite_face_ctrl_or_state_value(adjacent_face_for_control, domain_length) + ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(adjacent_face_for_control) *  x[
             ctrl::boundary_conditions::normal_direction_to_Gamma_control(adjacent_face_for_control) ] );  }
    else if(faceName == ctrl::boundary_conditions::opposite_face(adjacent_face_for_control)){ value = gamma * (
             ctrl::boundary_conditions::opposite_face_ctrl_or_state_value(face_for_control, domain_length) +
             ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(face_for_control) *  x[
             ctrl::boundary_conditions::normal_direction_to_Gamma_control(face_for_control) ] );  }

   value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

   return value;
}


 bool ctrl_or_state_set_dirichlet_flags(const MultiLevelProblem * ml_prob,
                                        const int faceName,
                                        const std::vector < double > & x,
                                        bool &  dirichlet)  {


     assert( femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size == 1 );
    const unsigned face_for_control = femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[0];


     if (faceName == face_for_control) {
        if ( !(x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] > femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[0][0] + 1.e-5 &&
               x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[0][1] - 1.e-5) ) {
                dirichlet = true;
           }
    }
    else {
          dirichlet = true;
    }

    return dirichlet;
}



}


namespace Gamma_c_double_adjacent_control_in_front_linear {




 double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

     if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

     assert( femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size == 2 );

     const unsigned face_for_control_principal = femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[0];
     const unsigned adjacent_face_for_control = femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[1];

     const double domain_length = 1.;

     const double gamma = 5.;


     for(unsigned f = 0; f < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

        if (faceName == femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f])     {  value = 0.; }
        else if(faceName == ctrl::boundary_conditions::opposite_face(face_for_control_principal)){ value = gamma * (
             ctrl::boundary_conditions::opposite_face_ctrl_or_state_value(adjacent_face_for_control, domain_length) + ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(adjacent_face_for_control) *  x[
             ctrl::boundary_conditions::normal_direction_to_Gamma_control(adjacent_face_for_control) ] );  }
        else if(faceName == ctrl::boundary_conditions::opposite_face(adjacent_face_for_control)){ value = gamma * (
             ctrl::boundary_conditions::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) + ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(face_for_control_principal) *  x[
             ctrl::boundary_conditions::normal_direction_to_Gamma_control(face_for_control_principal) ] );  }
     }

     value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

     return value;
}


 bool ctrl_or_state_set_dirichlet_flags(const MultiLevelProblem * ml_prob,
                                        const int faceName,
                                        const std::vector < double > & x,
                                        bool &  dirichlet)  {

     if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

     assert( femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size == 2 );


     bool  is_facename_a_control_face = false;

     for(unsigned f = 0; f < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

     if (faceName == femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) { is_facename_a_control_face = true; break; }

     }

     if  (is_facename_a_control_face == false) { dirichlet = true; }


     for(unsigned f = 0; f < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

     if (faceName == femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) {

        if ( !(x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] > femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][0] + 1.e-5 &&
               x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][1] - 1.e-5) ) {
                dirichlet = true;
           }
        }
      }


    return dirichlet;
}



}


namespace Gamma_c_triple_adjacent_control_in_front_linear {




 double ctrl_or_state_set_dirichlet_fixed_values(const MultiLevelProblem * ml_prob,
                                                 const int faceName,
                                                 const std::vector < double > & x,
                                                 double &  value)  {

     if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

     assert( femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size == 3 );

     const unsigned face_for_control_principal = femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[0];
     const unsigned adjacent_face_for_control = femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[1];

     const double domain_length = 1.;

     const double gamma = 5.;


     for(unsigned f = 0; f < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

         if (faceName == femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) {

                if ( (faceName == face_for_control_principal) || (faceName == adjacent_face_for_control) )     {  value = 0.; }
                else if(faceName == ctrl::boundary_conditions::opposite_face(face_for_control_principal)){
                        if ( !(x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] > femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][0] + 1.e-5 &&
                               x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][1] - 1.e-5) ) { value = gamma * (
                                  ctrl::boundary_conditions::opposite_face_ctrl_or_state_value(adjacent_face_for_control, domain_length) + ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(adjacent_face_for_control) *  x[
                                  ctrl::boundary_conditions::normal_direction_to_Gamma_control(adjacent_face_for_control) ] );  }
                        else {value = 0.; }
                }
                else if(faceName == ctrl::boundary_conditions::opposite_face(adjacent_face_for_control)){
                        if ( !(x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] > femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][0] + 1.e-5 &&
                               x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][1] - 1.e-5) ) { value = gamma * (
                               ctrl::boundary_conditions::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) + ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(face_for_control_principal) *  x[
                               ctrl::boundary_conditions::normal_direction_to_Gamma_control(face_for_control_principal) ] );  }
                               else {value = 0.; }
                }
         }
         else if(faceName == ctrl::boundary_conditions::opposite_face(face_for_control_principal)){ value = gamma * (
                                  ctrl::boundary_conditions::opposite_face_ctrl_or_state_value(adjacent_face_for_control, domain_length) + ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(adjacent_face_for_control) *  x[
                                  ctrl::boundary_conditions::normal_direction_to_Gamma_control(adjacent_face_for_control) ] );  }

         else if(faceName == ctrl::boundary_conditions::opposite_face(adjacent_face_for_control)){ value = gamma * (
                            ctrl::boundary_conditions::opposite_face_ctrl_or_state_value(face_for_control_principal, domain_length) + ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(face_for_control_principal) *  x[
                            ctrl::boundary_conditions::normal_direction_to_Gamma_control(face_for_control_principal) ] );  }

     }

     value += PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY;

     return value;
}


 bool ctrl_or_state_set_dirichlet_flags(const MultiLevelProblem * ml_prob,
                                        const int faceName,
                                        const std::vector < double > & x,
                                        bool &  dirichlet)  {

     if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

     assert( femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size == 3 );


     bool  is_facename_a_control_face = false;

     for(unsigned f = 0; f < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

     if (faceName == femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) { is_facename_a_control_face = true; break; }

     }

     if  (is_facename_a_control_face == false) { dirichlet = true; }


     for(unsigned f = 0; f < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {

     if (faceName == femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) {

        if ( !(x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] > femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][0] + 1.e-5 &&
               x[ ctrl::boundary_conditions::tangential_direction_to_Gamma_control(faceName) ] < femus::ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][1] - 1.e-5) ) {
                dirichlet = true;
           }
        }
      }


    return dirichlet;
}



}


}




namespace Gamma_control {
    

    
const double face_coordinate_extreme_position_normal_to_Gamma_control(const unsigned int face_index) {
    
  double extreme_pos;
  
        if (face_index == 1 || face_index == 3 || face_index == 5) {  extreme_pos = 0.; }
   else if (face_index == 2 || face_index == 4 || face_index == 6) {  extreme_pos = 1.; }
   
   return extreme_pos;
   
}





//*********************** Find volume elements that contain a Control domain element *********************************

int ControlDomainFlag_external_restriction(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
 // flag = 1: we are in the lifting nonzero domain
  int exterior_el_flag = 0.;
   if ( elem_center[0] >  1. -  1.e-5) { exterior_el_flag = 1; }

     return exterior_el_flag;

}



//*********************** Find volume elements that contain a Control domain element *********************************

int ControlDomainFlag_internal_restriction(const std::vector<double> & elem_center) {

 //***** set target domain flag ******
 // flag = 1: we are in the lifting nonzero domain
  int control_el_flag = 0.;
  
  const double offset_to_include_line =  OFFSET_TO_INCLUDE_LINE;
  
  const double control_domain_depth = LIFTING_INTERNAL_ORTHOGONAL_DISTANCE_FROM_GAMMA_C;
  
  const double control_domain_width_lower = LIFTING_INTERNAL_WIDTH_LOWER;
  const double control_domain_width_upper = LIFTING_INTERNAL_WIDTH_UPPER;
   
  
	  for(unsigned f = 0; f < ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {
          
   const int  line_sign = ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]);

   const double extreme_pos = Gamma_control::face_coordinate_extreme_position_normal_to_Gamma_control(ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]);

   const unsigned int axis_dir = boundary_conditions::tangential_direction_to_Gamma_control(ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]);

   
   if ( ( line_sign * elem_center[1 - axis_dir] <   line_sign * ( extreme_pos + line_sign * control_domain_depth ) )
       && ( elem_center[axis_dir] > control_domain_width_lower - offset_to_include_line ) 
       && ( elem_center[axis_dir] < control_domain_width_upper + offset_to_include_line ) )
      { control_el_flag = 1; }
   
      }
   
     return control_el_flag;

}





//*********************** Find volume elements that contain a Control Face element *********************************

int ControlDomainFlag_bdry(const std::vector<double> & elem_center) {

   
  int control_el_flag = 0;
  
  const double offset_to_include_line = OFFSET_TO_INCLUDE_LINE;

  const double control_domain_depth = BOUNDARY_ORTHOGONAL_DISTANCE_FROM_GAMMA_C; //this picks a lot more elements, but then the if on the faces only gets the control boundary

  
	  for(unsigned f = 0; f < ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {
          
   const int  line_sign = ctrl::boundary_conditions_or_cost_functional::sign_function_for_delimiting_region(ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]);

   const double extreme_pos = Gamma_control::face_coordinate_extreme_position_normal_to_Gamma_control(ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]);
   
   const unsigned int Gamma_c_dir_tangential = boundary_conditions::tangential_direction_to_Gamma_control(ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]);

  
   if ( ( line_sign * elem_center[1 - Gamma_c_dir_tangential] <   line_sign * (  extreme_pos  + line_sign * control_domain_depth) )
       && ( elem_center[Gamma_c_dir_tangential] > ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][0] - offset_to_include_line ) 
       && ( elem_center[Gamma_c_dir_tangential] < ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_extremes[f][1] + offset_to_include_line ) )
      { control_el_flag = 1; }
      
   }
                  
                  
                  
     return control_el_flag;
}

  
    
    
    

  
}



//*******************************************************************************************
//*********************** Domain and Mesh Dependent - END *****************************************
//*******************************************************************************************






//*******************************************************************************************
//*********************** Mesh independent - BEGIN *****************************************
//*******************************************************************************************

 namespace Gamma_control {
 

    
    
    
      std::vector< std::vector< int > > is_dof_associated_to_Gamma_control_equation(
      const Mesh * msh,
      /*const*/ MultiLevelSolution * ml_sol,
         const MultiLevelProblem *    ml_prob,
      const unsigned iel,
      CurrentElem < double > & geom_element_iel,
      const unsigned solType_coords,
      std::vector < std::string > Solname_Mat,      
      std::vector < unsigned > SolFEType_Mat,    
      std::vector < unsigned > Sol_n_el_dofs_Mat,    
      const unsigned pos_mat_ctrl,
      const unsigned n_components_ctrl
) {
//=============== construct control node flag field  =========================    
	      /* For every component:
           * (control_node_flag[c][i])       picks nodes on \Gamma_c
           * (1 - control_node_flag[c][i])   picks nodes on \Omega \setminus \Gamma_c
	       */
          
       std::vector< std::vector< int > > control_node_flag_iel_all_faces(n_components_ctrl);
       
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
              control_node_flag_iel_all_faces[c].resize(Sol_n_el_dofs_Mat[pos_mat_ctrl + c]);
              std::fill(control_node_flag_iel_all_faces[c].begin(), control_node_flag_iel_all_faces[c].end(), 0);   
         }
       
          
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

          
	    // look for boundary faces
            const int bdry_index = msh->el->GetFaceElementIndex(iel, iface);
            
	    if( bdry_index < 0) {
	      const unsigned int face_in_rectangle_domain = - ( msh->el->GetFaceElementIndex(iel, iface) + 1);
         
         
         
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
          
	     double tau = 0.;
         
	      const bool  dir_bool_c = ml_sol->GetBdcFunctionMLProb()(ml_prob, geom_element_iel.get_elem_center_bdry_3d(), Solname_Mat[pos_mat_ctrl + c].c_str(), tau, face_in_rectangle_domain, 0.);

	      if (dir_bool_c == false) {
              
          const unsigned ndofs_ctrl_bdry = msh->GetElementFaceDofNumber(iel, iface, SolFEType_Mat[pos_mat_ctrl + c]);
		  for(unsigned i_bdry = 0; i_bdry < ndofs_ctrl_bdry; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
		
				  control_node_flag_iel_all_faces[c][i_vol] = 1;
			    
              }
         }
      } 
        
        
      }   
    }
    
    return   control_node_flag_iel_all_faces;
    
  }

  
    
  bool face_is_a_Gamma_control_face(/*const*/ elem * el, const unsigned jel, const unsigned jface) {
      
  	    // look for boundary faces
            const int bdry_index_j = el->GetFaceElementIndex(jel, jface);
	    // look for face number equal to control face
	      const unsigned int face_in_rectangle_domain_j = - ( el->GetFaceElementIndex(jel,jface) + 1);

	    // look for boundary faces && look for control faces
		
          bool  is_face_for_control = false;
          
          		  for(unsigned f = 0; f < ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index_size; f++) {
                      if (face_in_rectangle_domain_j == ctrl::Gamma_control_list_of_faces_with_extremes::face_with_extremes_index[f]) { is_face_for_control = true; }
                  }

          
   return ( bdry_index_j < 0 && is_face_for_control );
      
  }

  
   

  
  
  

  
  bool volume_elem_contains_a_Gamma_control_face( const std::vector<double> & elem_center ) {

      int control_flag_jel = 0;
        control_flag_jel = ControlDomainFlag_bdry(elem_center);
        
        return (control_flag_jel == 1);

      
  }
  
    



  
}
  

 
  
 namespace Gamma_control_equation_fractional {
//*********************** Mesh independent, ALMOST - BEGIN *****************************************
 
  
  void mixed_integral(const unsigned unbounded,
                      const unsigned dim,
                      const unsigned dim_bdry,
                      const std::vector < std::vector < double > > & ex_control,
                      const unsigned div, 
                      const double weight_iqp_bdry,
                      const std::vector < double > & x_iqp_bdry,
                      const std::vector < double > & phi_ctrl_iel_bdry_iqp_bdry,
                      const std::vector<double> sol_ctrl_iqp_bdry,   //////////
                      const double s_frac,
                      const double check_limits,
                      const double C_ns,
                      const unsigned int operator_Hhalf,
                      const double beta,
                      const std::vector<unsigned int> nDof_vol_iel,  //////////
                      std::vector < double > & KK_local_iel,
                      std::vector < double > & Res_local_iel,
                      std::vector < double > & KK_local_iel_mixed_num,
                      std::vector < double > & Res_local_iel_mixed_num,
                      const Mesh * msh,
                      const Solution *    sol,
                      const MultiLevelSolution *    ml_sol,
                      const int iel,
                      const int jel,
                      const unsigned iface,
                      std::vector <int> bdry_bdry,
                      CurrentElem < double > & geom_element_jel,
                      const unsigned jelGeom_bdry,
                      unsigned solType_coords,
                      double & integral
                     ) {
      
     
  const unsigned  sol_node_flag_index =  ml_sol->GetIndex(NODE_BASED_BDRY_BDRY);
  const unsigned  group_salome = 2;   ///@todo fix here, maybe pass it in the args
  
  const unsigned int n_components_ctrl = nDof_vol_iel.size();
  
  unsigned nDof_iel_vec = 0;
        for (unsigned c = 0; c < n_components_ctrl; c++) {nDof_iel_vec  +=  nDof_vol_iel[c];    }

  
      
  if(unbounded == 1) {
      
      //============ Mixed Integral 1D - Analytical ==================      
      if (dim_bdry == 1) {
      
//               double ex_1 = EX_1;
//               double ex_2 = EX_2;
//               std::vector < double > ex_1_vec(dim);
//               std::vector < double > ex_2_vec(dim);
//               ex_1_vec[0] = EX_1;
//               ex_1_vec[1] = 0.;
//               ex_2_vec[0] = EX_2;
//               ex_2_vec[1] = 0.;
// //               ex_1_vec[0] = 1.;
// //               ex_1_vec[1] = EX_1;
// //               ex_2_vec[0] = 1.;
// //               ex_2_vec[1] = EX_2;
              
              double dist2_1 = 0.;
              double dist2_2 = 0.;
              double dist_1 = 0.;  //distance from node to extreme 1
              double dist_2 = 0.;  //distance from node to extreme 2
              
//               for(int d = 0; d < dim; d++) {
//                 dist2_1 += (x_iqp_bdry[d] - ex_1_vec[d]) * (x_iqp_bdry[d] - ex_1_vec[d]);
//                 dist2_2 += (x_iqp_bdry[d] - ex_2_vec[d]) * (x_iqp_bdry[d] - ex_2_vec[d]);
//               }
//               // ex_control Built as  ( x1 y1 ) ( x2 y2 ) 
//               
              for(int d = 0; d < dim; d++) {
                dist2_1 += (x_iqp_bdry[d] - ex_control[d][0]) * (x_iqp_bdry[d] - ex_control[d][0]);
                dist2_2 += (x_iqp_bdry[d] - ex_control[d][1]) * (x_iqp_bdry[d] - ex_control[d][1]);
//                 std::cout<< ex_control[d][0] << "  " << ex_control[d][1] << "\n";
              }
              
              
              dist_1 = sqrt( dist2_1 );
              dist_2 = sqrt( dist2_2 );
              
              double mixed_denominator = pow(dist_1, -2. * s_frac) + pow(dist_2, - 2. * s_frac);

              
              for (unsigned c = 0; c < n_components_ctrl; c++) {
                  integral +=  0.5 * C_ns * check_limits * operator_Hhalf  * beta * sol_ctrl_iqp_bdry[c] * sol_ctrl_iqp_bdry[c] * weight_iqp_bdry * mixed_denominator * (1. / s_frac);
              }   
              
              
              
              for (unsigned c = 0; c < n_components_ctrl; c++) {

              for(unsigned i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); i_bdry++) {
                unsigned int i_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
                
                const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_vol_iel, c, i_vol_iel);

                Res_local_iel[ res_pos/*i_vol_iel*/ ] += - 0.5 * C_ns * check_limits * operator_Hhalf  * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry[c] * weight_iqp_bdry * mixed_denominator * (1. / s_frac);
                
                for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) { 
                     for(unsigned j_bdry = 0; j_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); j_bdry++) {
                  unsigned int j_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);
                  const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_vol_iel, nDof_iel_vec, c, e, i_vol_iel, j_vol_iel);         

                  KK_local_iel[ jac_pos /*i_vol_iel * nDof_vol_iel + j_vol_iel*/ ] += 0.5 * C_ns * check_limits * operator_Hhalf * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry] * weight_iqp_bdry * mixed_denominator * (1. / s_frac);
                     }
                
                  }
                }
                
              }
           }
      }
      
    //============ Mixed Integral 2D - Numerical ==================      
      else if (dim_bdry == 2) {           
          
            unsigned iel_geom_type = msh->GetElementType(iel);
            unsigned iel_geom_type_face = msh->GetElementFaceType(iel, iface);

            unsigned f_n_faces_faces  =  msh->el->GetNFC(iel_geom_type, iel_geom_type_face); /* ElementFaceFaceNumber */

         double mixed_denominator_numerical = 0.;
         
            // *** Face Gauss point loop (Integral over 2d object) ***
            for(unsigned e_bdry_bdry = 0; e_bdry_bdry < f_n_faces_faces/*bdry_bdry.size()*/; e_bdry_bdry++) {  ///@todo I think this has to be fixed

              // look for boundary faces
            

//               unsigned n_dofs_bdry_bdry = msh->el->GetNFC(LINE, solType_coords);
              
              unsigned n_dofs_bdry_bdry =  msh->el->GetNFACENODES(iel_geom_type_face/*jelGeom_bdry*/, e_bdry_bdry, solType_coords); //TODO ///@todo this is only taking linear nodes, do we want to do that for the coordinates too?

              vector  < vector  <  double> > delta_coordinates_bdry_bdry(dim);    // A matrix holding the face coordinates rowwise.
              for(int k = 0; k < dim; k++) {
                delta_coordinates_bdry_bdry[k].resize(n_dofs_bdry_bdry);
              }
              
                  std::vector < int > nodes_face_face_flags(n_dofs_bdry_bdry, 0); 
              
              for(unsigned i_bdry_bdry = 0; i_bdry_bdry < n_dofs_bdry_bdry; i_bdry_bdry++) {
                  
                unsigned inode_bdry = msh->el->GetIG(jelGeom_bdry, e_bdry_bdry/*bdry_bdry[e_bdry_bdry]*/, i_bdry_bdry); // face-to-element local node mapping. TODO: verify jelGeom_bdry
                unsigned inode_vol    = msh->el->GetIG(iel_geom_type, iface, inode_bdry); //from n-1 to n

                unsigned node_global = msh->el->GetElementDofIndex(jel, inode_vol);
                
                nodes_face_face_flags[i_bdry_bdry] = (*sol->_Sol[sol_node_flag_index])(node_global);
                
              // delta coords  -----
                for(unsigned k = 0; k < dim; k++) {
                  delta_coordinates_bdry_bdry[k][i_bdry_bdry] = geom_element_jel.get_coords_at_dofs_3d()[k][inode_vol] - x_iqp_bdry[k];  ///@todo// TODO We extract the local coordinates on the face from local coordinates on the element.
                }
              }
              
              
                bool is_face_bdry_bdry  =  MED_IO::boundary_of_boundary_3d_check_face_of_face_via_nodes( nodes_face_face_flags, group_salome);

              if(is_face_bdry_bdry) {
                
              // delta coords - refinement -----
              
// // //               std::vector  <  double > delta_coordinates_bdry_bdry_refined( (div + 1) * dim);
              
              vector  < vector  <  double > > delta_coordinates_bdry_bdry_refined(dim);
              for(int k = 0; k < dim; k++) {
                delta_coordinates_bdry_bdry_refined[k].resize(div + 1); // set "4" as a parameter
              }
              
              for(unsigned n = 0; n <= div; n++) {
                for(int k = 0; k < dim; k++) {
// // //                   delta_coordinates_bdry_bdry_refined[n + k * div] = delta_coordinates_bdry_bdry[k][0] + n * (delta_coordinates_bdry_bdry[k][1] - delta_coordinates_bdry_bdry[k][0]) /  div ;
                  delta_coordinates_bdry_bdry_refined[k][n] = delta_coordinates_bdry_bdry[k][0] + n * (delta_coordinates_bdry_bdry[k][1] - delta_coordinates_bdry_bdry[k][0]) /  div ;
                }
              }
              for(unsigned n = 0; n < div; n++) {
                  
                const unsigned dir_x_for_atan = ( ( (FACE_FOR_CONTROL - 1) / 2 ) + 1 ) % 3;  ///@todo I think needs to be changed
                const unsigned dir_y_for_atan = ( dir_x_for_atan + 1 ) % 3 ;  ///@todo I think needs to be changed
// // //                 double teta2 = atan2(delta_coordinates_bdry_bdry_refined[(n+1) + dir_y_for_atan * div], delta_coordinates_bdry_bdry_refined[(n+1) + dir_x_for_atan * div]);
// // //                 double teta1 = atan2(delta_coordinates_bdry_bdry_refined[n + dir_y_for_atan * div], delta_coordinates_bdry_bdry_refined[n + dir_x_for_atan * div]);
                double teta2 = atan2(delta_coordinates_bdry_bdry_refined[dir_y_for_atan][n + 1], delta_coordinates_bdry_bdry_refined[dir_x_for_atan][n + 1]);
                double teta1 = atan2(delta_coordinates_bdry_bdry_refined[dir_y_for_atan][n], delta_coordinates_bdry_bdry_refined[dir_x_for_atan][n]);

// // //                 double delta_teta = 0.;
// // //                 if(teta2 < teta1) delta_teta = std::min(teta1 - teta2, 2. * M_PI + teta2 - teta1);
// // //                 else delta_teta = std::min(teta2 - teta1, 2. * M_PI + teta1 - teta2);
                if(teta2 < teta1) {
                    teta2 += 2. * M_PI;
                }
                
                double delta_teta = teta2 - teta1;

                vector <double> mid_point;
                mid_point.resize(dim);
                for(unsigned k = 0; k < dim; k++) {
// // //                   mid_point[k] = (delta_coordinates_bdry_bdry_refined[(n+1) + k * div] + delta_coordinates_bdry_bdry_refined[n + k * div]) * 0.5;
                  mid_point[k] = (delta_coordinates_bdry_bdry_refined[k][n + 1] + delta_coordinates_bdry_bdry_refined[k][n]) * 0.5;
                }
                double dist2 = 0;
                for(int k = 0; k < dim; k++) {
                  dist2 += mid_point[k] * mid_point[k];
                }
                double dist = sqrt(dist2);
                
                mixed_denominator_numerical += pow(dist, -  2. * s_frac) * delta_teta;
              }
//               delta coords - refinement -----

          
              }
              
            }
            

            
            
            
//              for(unsigned i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); i_bdry++) {
//                 unsigned int i_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
//                 
//                 Res_local_iel_mixed_num[ i_vol_iel ] += - 0.5 * C_ns * check_limits * operator_Hhalf * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry * weight_iqp_bdry * mixed_denominator_numerical  * (1. / s_frac);
//                 
//                 for(unsigned j_bdry = 0; j_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); j_bdry++) {
//                   unsigned int j_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);
//                   KK_local_iel_mixed_num[ i_vol_iel * nDof_vol_iel + j_vol_iel ] += 0.5 * C_ns * check_limits * operator_Hhalf * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry]  * weight_iqp_bdry * mixed_denominator_numerical * (1. / s_frac);
//                 }
//                 
//               }
            
            /// @todo ONLY DIFFERENCES: the mixed_denominator is numerical, and so also the corresponding Res and Jac. It could be done with a single function
            
              for (unsigned c = 0; c < n_components_ctrl; c++) {
                 integral +=  0.5 * C_ns * check_limits * operator_Hhalf  * beta  * sol_ctrl_iqp_bdry[c] * sol_ctrl_iqp_bdry[c] * weight_iqp_bdry * mixed_denominator_numerical * (1. / s_frac);
              }
              
              
              
              for (unsigned c = 0; c < n_components_ctrl; c++) {

              for(unsigned i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); i_bdry++) {
                unsigned int i_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
                
                const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_vol_iel, c, i_vol_iel);

                Res_local_iel_mixed_num[ res_pos/*i_vol_iel*/ ] += - 0.5 * C_ns * check_limits * operator_Hhalf  * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry[c] * weight_iqp_bdry * mixed_denominator_numerical * (1. / s_frac);
                
                for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) { 
                     for(unsigned j_bdry = 0; j_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); j_bdry++) {
                  unsigned int j_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);
                  const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_vol_iel, nDof_iel_vec, c, e, i_vol_iel, j_vol_iel);         

                  KK_local_iel_mixed_num[ jac_pos /*i_vol_iel * nDof_vol_iel + j_vol_iel*/ ] += 0.5 * C_ns * check_limits * operator_Hhalf * beta * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry] * weight_iqp_bdry * mixed_denominator_numerical * (1. / s_frac);
                     }
                
                  }
                }
                
              }
           }
            
          
      }          
      
  }
  
  
 }
  
  
    void   bdry_bdry_flag_copy_and_delete(   MultiLevelProblem & ml_prob,
                                             MultiLevelSolution & ml_sol,
                                           const MultiLevelMesh & ml_mesh, 
                                              const unsigned erased_levels,
                                             MultiLevelSolution *  ml_sol_bdry_bdry_flag,
                                             const std::string node_based_bdry_bdry_flag_name,
                                             const unsigned  steady_flag,
                                             const bool      is_an_unknown_of_a_pde,
                                             const FEFamily node_bdry_bdry_flag_fe_fam,
                                             const FEOrder node_bdry_bdry_flag_fe_ord) {
        
 ml_sol.AddSolution(node_based_bdry_bdry_flag_name.c_str(), node_bdry_bdry_flag_fe_fam, node_bdry_bdry_flag_fe_ord, steady_flag, is_an_unknown_of_a_pde);
  ml_sol.Initialize(node_based_bdry_bdry_flag_name.c_str());
  // copy ml_sol_bdry_bdry_flag at the non-removed levels into ml_sol
  for(unsigned l = 0; l < ml_mesh.GetNumberOfLevels(); l++) {
      *(ml_sol.GetSolutionLevel(l)->_Sol[ ml_sol.GetIndex(node_based_bdry_bdry_flag_name.c_str()) ]) =
      *(ml_sol_bdry_bdry_flag->GetSolutionLevel(l + erased_levels)->_Sol[ ml_sol_bdry_bdry_flag->GetIndex(node_based_bdry_bdry_flag_name.c_str()) ]);
  }
  
  delete ml_sol_bdry_bdry_flag;
     
     }  

  
     MultiLevelSolution *  bdry_bdry_flag(  Files & files,
                                            MultiLevelMesh & ml_mesh, 
                                               const std::string infile,
                                               std::vector < unsigned > & node_mapping_from_mesh_file_to_new,
                                             const std::string node_based_bdry_bdry_flag_name,
                                             const unsigned  steady_flag,
                                             const bool      is_an_unknown_of_a_pde,
                                             const FEFamily node_bdry_bdry_flag_fe_fam,
                                             const FEOrder node_bdry_bdry_flag_fe_ord) {
     
     
  MultiLevelSolution * ml_sol_bdry_bdry_flag = new MultiLevelSolution(&ml_mesh);
  ml_sol_bdry_bdry_flag->SetWriter(VTK);
  ml_sol_bdry_bdry_flag->GetWriter()->SetDebugOutput(true);
  

  ml_sol_bdry_bdry_flag->AddSolution(node_based_bdry_bdry_flag_name.c_str(), node_bdry_bdry_flag_fe_fam, node_bdry_bdry_flag_fe_ord, steady_flag, is_an_unknown_of_a_pde);
  ml_sol_bdry_bdry_flag->Initialize(node_based_bdry_bdry_flag_name.c_str());

  // ======= COARSE READING and REFINEMENT ========================
  ml_sol_bdry_bdry_flag->GetSolutionLevel(0)->GetSolutionName(node_based_bdry_bdry_flag_name.c_str()) = MED_IO(*ml_mesh.GetLevel(0)).node_based_flag_read_from_file(infile, node_mapping_from_mesh_file_to_new);

  ml_mesh.GetLevelZero(0)->deallocate_node_mapping(node_mapping_from_mesh_file_to_new);

  for(unsigned l = 1; l < ml_mesh.GetNumberOfLevels(); l++) {
     ml_sol_bdry_bdry_flag->RefineSolution(l);
  }
  
  std::vector < std::string > variablesToBePrinted_aux;
  variablesToBePrinted_aux.push_back("all");
  for(unsigned l = 0; l < ml_mesh.GetNumberOfLevels(); l++) {
  ml_sol_bdry_bdry_flag->GetWriter()->Write(l+1, "aux", files.GetOutputPath(), "", "biquadratic", variablesToBePrinted_aux);
   }
 
 
 return ml_sol_bdry_bdry_flag;
 
}


  
  
  double  compute_C_ns(const unsigned dim_bdry,
                       const double s_frac,
                       const unsigned use_Cns) {
      
  const double C_ns = 2 * (1 - use_Cns) + use_Cns * s_frac * pow(2, (2. * s_frac)) * tgamma((dim_bdry + 2. * s_frac) / 2.) / (pow(M_PI, dim_bdry / 2.) * tgamma(1 -  s_frac)) ;
  
  return C_ns;
}
  
    //********** FRAC CONTROL - BEGIN *****************************************

  void control_eqn_bdry_fractional_sobolev_differentiability_index(const unsigned iproc,
                                   const unsigned nprocs,
                        MultiLevelProblem &    ml_prob,
                        MultiLevelSolution*    ml_sol,
                        const Solution*        sol,
                        const Mesh * msh,
                        const  LinearEquationSolver* pdeSys,
                        //-----------
                        CurrentElem < double > & geom_element_iel,
                        const unsigned int solType_coords,
                        const unsigned int dim,
                        const unsigned int space_dim,
                        const unsigned int dim_bdry,
                        //-----------
                        const unsigned int n_unknowns,
                        const    vector < std::string > & Solname_Mat,
                        const    vector < unsigned > & SolFEType_Mat,
                        const vector < unsigned > & SolIndex_Mat,
                        const vector < unsigned > & SolPdeIndex,
                        vector < unsigned > & Sol_n_el_dofs_Mat, 
                        vector < vector < double > > & sol_eldofs_Mat,  
                        vector < vector < int > > & L2G_dofmap_Mat,
                        const unsigned int max_size,
                        //-----------
                         const unsigned int n_quantities,
                        vector < unsigned > SolFEType_quantities,
                        //-----------
                        std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all,
                        //-----------
                        std::vector < std::vector < double > >  Jac_iel_bdry_iqp_bdry,
                        std::vector < std::vector < double > >  JacI_iel_bdry_iqp_bdry,
                        double detJac_iel_bdry_iqp_bdry,
                        double weight_iqp_bdry,
                        vector <double> phi_ctrl_iel_bdry_iqp_bdry,
                        vector <double> phi_ctrl_x_iel_bdry_iqp_bdry, 
                        //---- Control unknown -------
                        const unsigned int n_components_ctrl,
                        const unsigned int pos_mat_ctrl,
                        const unsigned int pos_sol_ctrl,
                        const unsigned int is_block_dctrl_ctrl_inside_main_big_assembly,
                        //-----------
                        SparseMatrix*  KK,
                        NumericVector* RES,
                        const bool assembleMatrix,
                        //-----------
                        const double alpha,
                        const double beta,
                        //-----------
                        const double s_frac,
                        const double check_limits,
                        const unsigned use_Cns,
                        const unsigned int operator_Hhalf,
                        const unsigned int operator_L2,
                        const double rhs_one,
                        const unsigned int unbounded,
                        //--- Quadrature --------
                        const unsigned qrule_i,
                        const unsigned qrule_j,
                        const unsigned qrule_k,
                        const unsigned int integration_num_split,
                        const unsigned int integration_split_index,
                        const unsigned int N_div_unbounded,
                        //-----------
                        const bool print_algebra_local
                       ) {

// --- Fractional - BEGIN
const double C_ns =    compute_C_ns(dim_bdry, s_frac, use_Cns);  
// --- Fractional - END
      
      
// --- Geometry - BEGIN
  CurrentElem < double > geom_element_jel(dim, msh);            // must be adept if the domain is moving, otherwise double
// --- Geometry - END
      
// --- Quadrature - BEGIN
     std::vector < std::vector < double > >  JacI_jel_bdry_jqp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_jel_bdry_jqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_jel_bdry_jqp_bdry.size(); d++) {   Jac_jel_bdry_jqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_jel_bdry_jqp_bdry.size(); d++) { JacI_jel_bdry_jqp_bdry[d].resize(dim-1); }
    
    double detJac_jel_bdry_jqp_bdry;
// --- Quadrature - END


//--- Control domain - BEGIN -------
//*************************************************** 
  unsigned n_max = pow(2,dim_bdry);
  std::vector < double > extremes(n_max);
  extremes[0] = EX_1;
  extremes[1] = EX_2;
  if(dim_bdry == 2){
    extremes[2] = EY_1;
    extremes[3] = EY_2;
  }
  
  std::vector < std::vector < double > > ex_control(dim); // control zone extremes. Built as: 2D) ( x1 y1 ) ( x2 y2 )    3D) ( x1 y1 z1 ) ( x2 y2 z2 ) ( x3 y3 z3 ) 
  for(unsigned d = 0; d < dim; d++) {
      ex_control[d].reserve(n_max);
  }
  int control_xyz = (FACE_FOR_CONTROL - 1) / 2;
  bool ctrl_min_max = (FACE_FOR_CONTROL - 1) % 2;
  for(unsigned d = 0; d < dim; d++) {
    for(unsigned n_e = 0; n_e < n_max; n_e++){
      if(control_xyz == d) ex_control[d][n_e] =  (ctrl_min_max)? DOMAIN_EX_2:DOMAIN_EX_1;
      else ex_control[d][n_e] = extremes[n_e];
    }
  }
//*************************************************** 


//--- Control domain - END -------
  


 //********************* bdry cont *******************
 //*************************************************** 
  std::vector <double> phi_coords_iel_bdry_iqp_bdry;  
  std::vector <double> phi_coords_x_iel_bdry_iqp_bdry; 

  phi_coords_iel_bdry_iqp_bdry.reserve(max_size);
  phi_coords_x_iel_bdry_iqp_bdry.reserve(max_size * space_dim);

 //*************************************************** 

  
  std::vector <unsigned> solType_ctrl(n_components_ctrl);
  std::vector <unsigned> solIndex_ctrl(n_components_ctrl);
  std::vector <unsigned> solPdeIndex_ctrl(n_components_ctrl);
  
  for (unsigned c = 0; c < n_components_ctrl; c++) {
          solType_ctrl[c]     = SolFEType_Mat[pos_mat_ctrl + c];
          solIndex_ctrl[c]    = SolIndex_Mat[pos_mat_ctrl + c];
          solPdeIndex_ctrl[c] = SolPdeIndex[pos_mat_ctrl + c];
  }
  
  std::vector < std::vector < double > > sol_ctrl_iel(n_components_ctrl);  ///@todo 1
  std::vector < std::vector < double > > sol_ctrl_jel(n_components_ctrl);
    
    for (unsigned c = 0; c < n_components_ctrl; c++) {
        sol_ctrl_iel[c].reserve(max_size);
        sol_ctrl_jel[c].reserve(max_size);
    }
    
    
   const unsigned  elem_dof_size_max = n_components_ctrl * max_size;
 
 
  //-------- local to global mappings - BEGIN --------------
  vector< vector< int > > l2gMap_iel(n_components_ctrl); 
  vector< vector< int > > l2gMap_jel(n_components_ctrl);

    for (unsigned c = 0; c < n_components_ctrl; c++) {
        l2gMap_iel[c].reserve(max_size);
        l2gMap_jel[c].reserve(max_size);
    }
  
  vector< int > l2gMap_iel_vec;  l2gMap_iel_vec.reserve(elem_dof_size_max);
  vector< int > l2gMap_jel_vec;  l2gMap_jel_vec.reserve(elem_dof_size_max);
  //-------- local to global mappings - END --------------


  //-------- Local matrices and rhs - BEGIN --------------
  vector < double > Res_local_iel; Res_local_iel.reserve(elem_dof_size_max);
  vector < double > KK_local_iel;  KK_local_iel.reserve(elem_dof_size_max * elem_dof_size_max);

//   Local matrices and rhs for adaptive quadrature (iel == jel)
  vector < double > Res_local_iel_refined; Res_local_iel_refined.reserve(elem_dof_size_max);
  vector < double > KK_local_iel_refined;   KK_local_iel_refined.reserve(elem_dof_size_max * elem_dof_size_max);

//   Local matrices and rhs for the mixed internal-external integral term (both adaptive and non-adaptive)
  vector < double > Res_local_iel_mixed_num;  Res_local_iel_mixed_num.reserve(elem_dof_size_max);
  vector < double > KK_local_iel_mixed_num;   KK_local_iel_mixed_num.reserve(elem_dof_size_max * elem_dof_size_max);

//   Non local matrices and vectors for H^s laplacian operator
  vector< double >         Res_nonlocal_iel;  Res_nonlocal_iel.reserve(elem_dof_size_max);
  vector< double >         Res_nonlocal_jel;  Res_nonlocal_jel.reserve(elem_dof_size_max);

  vector < double > KK_nonlocal_iel_iel;  KK_nonlocal_iel_iel.reserve(elem_dof_size_max * elem_dof_size_max);
  vector < double > KK_nonlocal_iel_jel;  KK_nonlocal_iel_jel.reserve(elem_dof_size_max * elem_dof_size_max);
  vector < double > KK_nonlocal_jel_iel;  KK_nonlocal_jel_iel.reserve(elem_dof_size_max * elem_dof_size_max);
  vector < double > KK_nonlocal_jel_jel;  KK_nonlocal_jel_jel.reserve(elem_dof_size_max * elem_dof_size_max); 
  //-------- Local matrices and rhs - END --------------

  
  //-------- Global matrices and rhs - BEGIN --------------
  KK->zero();
  RES->zero(); 
  //-------- Global matrices and rhs - END --------------

  
//   phi_x_placeholder_unused as unused input of certain functions
  vector < double > phi_x_placeholder_unused;
  phi_x_placeholder_unused.reserve(max_size * dim);
  
 
///   boost::mpi::communicator world(MPI_COMM_WORLD, boost::mpi::comm_attach);  /// @todo future solution: broadcast whole class instances


  unsigned count_visits_of_boundary_faces = 0;


 // integral - BEGIN ************
  double integral  = 0.;
// integral - END ************
      

  
   for(int kproc = 0; kproc < nprocs; kproc++) {
       
  const int proc_to_bcast_from = kproc;
  
    for(int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++) {

// --- geometry
   // all these little vectors are filled in one proc and broadcast to all       

        // --- 
        unsigned nDof_jel_coords;
         
      if (iproc == kproc) {
        nDof_jel_coords = msh->GetElementDofNumber(jel, solType_coords);
      }
     MPI_Bcast(& nDof_jel_coords, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
    // ---       
        
    // --- 
         unsigned short  jel_geommm;
      if (iproc == kproc) {
       jel_geommm = msh->el->GetElementType(jel);
//         std::cout  << " current_proc " << iproc << " from external_proc " << kproc << " elem " << jel << " (before bcast) "  << jel_geommm << std::endl;
//         geom_element_jel.set_geom_type(jel);
//         jel_geom = geom_element_jel.geom_type();
     }
      MPI_Bcast(& jel_geommm, 1, MPI_UNSIGNED_SHORT, proc_to_bcast_from, MPI_COMM_WORLD);
//         std::cout << " current_proc " << iproc << " from external_proc " << kproc << " elem " << jel << " (after  bcast) "  << jel_geommm << std::endl;
        // --- 

        // --- coords - other way
        geom_element_jel.allocate_coords_at_dofs_3d(jel, nDof_jel_coords, solType_coords);
        
      if(kproc == iproc) {
        geom_element_jel.fill_coords_at_dofs_3d(jel, solType_coords);
      }
      for(unsigned k = 0; k < dim; k++) {
        MPI_Bcast(& geom_element_jel.get_coords_at_dofs()[k][0], nDof_jel_coords, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
      }
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element_jel.get_coords_at_dofs_3d()[k][0], nDof_jel_coords, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
      }
// --- coords - other way

      if(kproc == iproc) {
        geom_element_jel.set_elem_center_3d(jel, solType_coords);
      }
        MPI_Bcast(& geom_element_jel.get_elem_center_3d()[0], space_dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);

// --- geometry        

        
// // // all of this is not used right now in this routine        
// // //       if(kproc == iproc) {
// // //  //***************************************************
// // //    el_dofs_unknowns_vol(sol, msh, pdeSys, jel,
// // //                         SolFEType_Mat,
// // //                         SolIndex_Mat,
// // //                         SolPdeIndex,
// // //                         Sol_n_el_dofs_Mat, 
// // //                         sol_eldofs_Mat,  
// // //                         L2G_dofmap_Mat);  //all unknowns here, perhaps we could restrict it to the ctrl components only
// // //       }
// // //         MPI_Bcast(& SolFEType_Mat[0], , MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// // //         MPI_Bcast(& SolIndex_Mat[0], , MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// // //         MPI_Bcast(& SolPdeIndex[0], , MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// // //         MPI_Bcast(& Sol_n_el_dofs_Mat[0], , MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// // //   //***************************************************

   
	// Perform face loop over elements that contain some control face
        
        
	if ( ctrl::Gamma_control::volume_elem_contains_a_Gamma_control_face(geom_element_jel.get_elem_center_3d()) ) {

      
// ***************************************
// ******* jel-related stuff - BEGIN *************
// ***************************************

// --- 2 - solution -----------------
      vector< unsigned > nDof_jel(n_components_ctrl);

      if(kproc == iproc) {
        for (unsigned c = 0; c < n_components_ctrl; c++) {
           nDof_jel[c]  = msh->GetElementDofNumber(jel, solType_ctrl[c]);    // number of solution element dofs
           }
    }

      MPI_Bcast(& nDof_jel[0], n_components_ctrl, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
      
      
      unsigned nDof_jel_vec = 0;
        for (unsigned c = 0; c < n_components_ctrl; c++) { nDof_jel_vec  +=  nDof_jel[c];    }

      
      
        for (unsigned c = 0; c < n_components_ctrl; c++) {
           sol_ctrl_jel[c].resize(nDof_jel[c]);
        }
        
      if(kproc == iproc) {
          
      for (unsigned c = 0; c < n_components_ctrl; c++) {
       for(unsigned j = 0; j < nDof_jel[c]; j++) {
          unsigned jDof  = msh->GetSolutionDof(j, jel, solType_ctrl[c]);
          sol_ctrl_jel[c][j] = (*sol->_Sol[solIndex_ctrl[c]])(jDof);
           }
         }
      }
      
      for (unsigned c = 0; c < n_components_ctrl; c++) {
         MPI_Bcast(& sol_ctrl_jel[c][0], nDof_jel[c], MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
      }
    // --- 2 - solution -----------------

      
// --- 3 - l2GMap -----------------
      for (unsigned c = 0; c < n_components_ctrl; c++) {
          l2gMap_jel[c].resize(nDof_jel[c]);
      }
           
           
      // local storage of global mapping and solution ********************
      if(kproc == iproc) {
       for (unsigned c = 0; c < n_components_ctrl; c++) {
         for(unsigned j = 0; j < nDof_jel[c]; j++) {
          l2gMap_jel[c][j] = pdeSys->GetSystemDof(solIndex_ctrl[c], solPdeIndex_ctrl[c], j, jel);  // global to global mapping between solution node and pdeSys dof
         }
       }
    }
    
      for (unsigned c = 0; c < n_components_ctrl; c++) {
         MPI_Bcast(& l2gMap_jel[c][0], nDof_jel[c], MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
      }
      // ******************************************************************
// --- 3 - l2GMap -----------------



// --- l2GMapVec -----------------
    l2gMap_jel_vec.resize(0);
      for (unsigned c = 0; c < n_components_ctrl; c++) {
          l2gMap_jel_vec.insert(l2gMap_jel_vec.end(), l2gMap_jel[c].begin(), l2gMap_jel[c].end());
      }
// --- l2GMapVec -----------------



// ***************************************
// ******* jel-related stuff - END *************
// ***************************************


      
//------------------------------------        
//------------ jface opening ---------        
//------------------------------------        
      
// --- 
      unsigned n_faces_jel;
      if(kproc == iproc) {
          n_faces_jel = msh->GetElementFaceNumber(jel); 
    }
      MPI_Bcast(& n_faces_jel, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
// ---
      
	  // loop on faces of the current element
	  for(unsigned jface = 0; jface < n_faces_jel; jface++) {
          
          
// --- geometry        
       unsigned jelGeom_bdry;
       if(kproc == iproc) {
           jelGeom_bdry = msh->GetElementFaceType(jel, jface);
        }  
      MPI_Bcast(& jelGeom_bdry, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);

        unsigned nve_bdry;
       if(kproc == iproc) {
          nve_bdry  = msh->GetElementFaceDofNumber(jel, jface, solType_coords);
        }  
      MPI_Bcast(& nve_bdry, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);

      
       geom_element_jel.allocate_coords_at_dofs_bdry_3d(jel, jface, nve_bdry);
      unsigned coords_at_dofs_bdry_3d_size = 0;
      
       if(kproc == iproc) {
         geom_element_jel.fill_coords_at_dofs_bdry_3d(jel, jface, solType_coords);
         coords_at_dofs_bdry_3d_size = geom_element_jel.get_coords_at_dofs_bdry_3d()[0].size();  ///@todo all coords have the same ndofs
        }  
      MPI_Bcast(& coords_at_dofs_bdry_3d_size, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
      
      for(unsigned k = 0; k < space_dim; k++) {
         MPI_Bcast(& geom_element_jel.get_coords_at_dofs_bdry_3d()[k][0], coords_at_dofs_bdry_3d_size, MPI_DOUBLE, kproc, MPI_COMM_WORLD);
      }
            
      
           std::fill(geom_element_jel.get_elem_center_bdry_3d().begin(), geom_element_jel.get_elem_center_bdry_3d().end(), 0.);

       if(kproc == iproc) {
       geom_element_jel.set_elem_center_bdry_3d();
        }  
      for(unsigned k = 0; k < space_dim; k++) {
        MPI_Bcast(& geom_element_jel.get_elem_center_bdry_3d()[0], space_dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
      }
// --- geometry        

     /*bool*/int jface_is_a_boundary_control;
       if(kproc == iproc) {
           jface_is_a_boundary_control = ctrl::Gamma_control::face_is_a_Gamma_control_face(msh->el, jel, jface);
       }
      MPI_Bcast(& jface_is_a_boundary_control, 1, MPI_INTEGER, proc_to_bcast_from, MPI_COMM_WORLD);

      
	    if( jface_is_a_boundary_control ) {
//------------------------------------        
//------------ jface opening ---------    
//------------------------------------        

            

// // // //---- Quadrature in jqp_bdry, preparation right before iel - BEGIN ------- 

// Wait a second... If I want to prepare the jqp_bdry loop, I must be inside jface as well...
// So now it seems to me that I have to do jel jface iel iface instead...
// Previously it was jel iel iqp jqp
// Now, it has to be jel jface  - iel iface - iqp_bdry jqp_bdry
// The two quadrature loops must be the innermost. In this way you exclude all non-needed volume elements and all non-needed faces, so that you minimize the number of inner ifs. You keep them outside as much as possible
// There will be a storage of jqp_bdry

            
      const unsigned n_jqp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_j, jelGeom_bdry).GetGaussPointsNumber();

      std::vector < std::vector < double > > x_jqp_bdry(n_jqp_bdry);
      std::vector < double > weight_jqp_bdry(n_jqp_bdry);
      std::vector < std::vector < double > > phi_coords_jel_bdry_jqp_bdry(n_jqp_bdry);
      std::vector < std::vector < double > > phi_coords_x_jel_bdry_jqp_bdry(n_jqp_bdry);
      
      std::vector < std::vector < double > > phi_ctrl_jel_bdry_jqp_bdry(n_jqp_bdry);   /// @todo assume all ctrl components have the same FE family
      std::vector < std::vector < double > > phi_ctrl_x_jel_bdry_jqp_bdry(n_jqp_bdry);   /// @todo assume all ctrl components have the same FE family

      std::vector< std::vector< double > > sol_ctrl_jqp_bdry(n_components_ctrl);
      
        for (unsigned c = 0; c < n_components_ctrl; c++) {
             sol_ctrl_jqp_bdry[c].resize(n_jqp_bdry);
                 std::fill(sol_ctrl_jqp_bdry[c].begin(), sol_ctrl_jqp_bdry[c].end(), 0.);
            }

      
         for(unsigned jqp_bdry = 0; jqp_bdry < n_jqp_bdry; jqp_bdry++) {

    elem_all[qrule_j][jelGeom_bdry][solType_coords]->JacJacInv(geom_element_jel.get_coords_at_dofs_bdry_3d(), jqp_bdry, Jac_jel_bdry_jqp_bdry, JacI_jel_bdry_jqp_bdry, detJac_jel_bdry_jqp_bdry, space_dim);
    
    weight_jqp_bdry[jqp_bdry] = detJac_jel_bdry_jqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_j, jelGeom_bdry).GetGaussWeightsPointer()[jqp_bdry];

    elem_all[qrule_j][jelGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(jqp_bdry, JacI_jel_bdry_jqp_bdry, phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry], phi_ctrl_x_jel_bdry_jqp_bdry[jqp_bdry], boost::none, space_dim);
            
    elem_all[qrule_j][jelGeom_bdry][solType_coords] ->shape_funcs_current_elem(jqp_bdry, JacI_jel_bdry_jqp_bdry, phi_coords_jel_bdry_jqp_bdry[jqp_bdry], phi_coords_x_jel_bdry_jqp_bdry[jqp_bdry], boost::none, space_dim);

//========== compute gauss quantities on the boundary ===============================================
//--- geom
           x_jqp_bdry[jqp_bdry].assign(dim, 0.);
          
//        if(kproc == iproc) {
         for(unsigned d = 0; d < dim; d++) {
	      for (int j_bdry = 0; j_bdry < geom_element_jel.get_coords_at_dofs_bdry_3d()[d].size(); j_bdry++)  {
			
              x_jqp_bdry[jqp_bdry][d] += geom_element_jel.get_coords_at_dofs_bdry_3d()[d][j_bdry] * phi_coords_jel_bdry_jqp_bdry[jqp_bdry][j_bdry];

		      }
            }
//        }
//       MPI_Bcast(& x_jqp_bdry[jqp_bdry][0], dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
//--- geom
    
//--- solution
      for (unsigned c = 0; c < n_components_ctrl; c++) {
         sol_ctrl_jqp_bdry[c][jqp_bdry] = 0.;
//        if(kproc == iproc) {
	      for (int j_bdry = 0; j_bdry < phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry].size()/*Sol_n_el_dofs_quantities[pos_sol_ctrl]*/; j_bdry++)  {
		    unsigned int j_vol = msh->GetLocalFaceVertexIndex_PassElemType(jel_geommm, jface, j_bdry);
			
			sol_ctrl_jqp_bdry[c][jqp_bdry] +=  /*sol_eldofs_Mat[pos_mat_ctrl]*/sol_ctrl_jel[c][j_vol] * phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][j_bdry];

	      }
       }
//        }
//       MPI_Bcast(& sol_ctrl_jqp_bdry[dim], 1, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
//--- solution
//========== compute gauss quantities on the boundary ================================================


        }  //jqp_bdry

// // //         //we can do the broadcast after the loop, faster
// // //         for(unsigned jqp_bdry = 0; jqp_bdry < n_jqp_bdry; jqp_bdry++) {
// // //             MPI_Bcast(& x_jqp_bdry[jqp_bdry][0], dim, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
// // //         }
// // //             MPI_Bcast(& sol_ctrl_jqp_bdry[0], n_jqp_bdry, MPI_DOUBLE, proc_to_bcast_from, MPI_COMM_WORLD);
  
        
// // // //---- Quadrature in jqp_bdry, preparation right before iel - END ------- 

    
// ---- boundary faces in jface: compute and broadcast - BEGIN ----
// This is only needed for when the boundary is a 2D face. We'll look at it later on
// look for what face of jface are on the boundary of the domain
// I believe we have to see deeply how we can extend this to the boundary case
// TODO: instead of faceIndex we want to have the new group condition that we impose 
// on the boundary of the boundary.


       std::vector< int > bdry_bdry(0);
///       unsigned n_faces;
/// 
///       if(iproc == kproc) {
///         for(unsigned j_bd_face = 0; j_bd_face < msh->GetElementFaceNumber_PassElemType(jelGeom_bdry); j_bd_face++) {
///           int faceIndex = msh->el->GetBoundaryIndex(jface, j_bd_face); /// TODO find a new condition and correct msh->GetElementFaceNumber ///@todo this is wrong
/// 
///           // look for boundary faces of the boundary
///           if(faceIndex >= 1) {
///             unsigned i = bdry_bdry.size();
///             bdry_bdry.resize(i + 1);
///             bdry_bdry[i] = jface;
///           }
///         }
///         n_faces = bdry_bdry.size();
///       }
/// 
///       MPI_Bcast(& n_faces, 1, MPI_UNSIGNED, proc_to_bcast_from, MPI_COMM_WORLD);
/// 
///       bdry_bdry.resize(n_faces);
///       MPI_Bcast(& bdry_bdry[0], n_faces, MPI_INT, proc_to_bcast_from, MPI_COMM_WORLD);  

// ---- boundary faces in jface: compute and broadcast - END ----    


              
       for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
           
           
// --- geometry        
        geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);

        short unsigned ielGeom = geom_element_iel.geom_type();
              
        geom_element_iel.set_elem_center_3d(iel, solType_coords);
// --- geometry        

      
	// Perform face loop over elements that contain some control face
	if ( ctrl::Gamma_control::volume_elem_contains_a_Gamma_control_face( geom_element_iel.get_elem_center_3d() ) ) {
        
// ***************************************
// ******* iel-related stuff - BEGIN *************
// ***************************************
        
      
// --- l2GMap
        vector< unsigned > nDof_iel(n_components_ctrl);
        
       for (unsigned c = 0; c < n_components_ctrl; c++) {
         nDof_iel[c]  = msh->GetElementDofNumber(iel, solType_ctrl[c]);
       }
              
      for (unsigned c = 0; c < n_components_ctrl; c++) {
        l2gMap_iel[c].resize(nDof_iel[c]);
        for(unsigned i = 0; i < l2gMap_iel[c].size(); i++) {
          l2gMap_iel[c][i] = pdeSys->GetSystemDof(solIndex_ctrl[c], solPdeIndex_ctrl[c], i, iel);
           }
         }
// --- l2GMap

// --- l2GMapVec -----------------
    l2gMap_iel_vec.resize(0);
      for (unsigned c = 0; c < n_components_ctrl; c++) {
          l2gMap_iel_vec.insert(l2gMap_iel_vec.end(), l2gMap_iel[c].begin(), l2gMap_iel[c].end());
      }
// --- l2GMapVec -----------------


// --- element matrix and vector resizing
unsigned nDof_iel_vec = 0;
        for (unsigned c = 0; c < n_components_ctrl; c++) {  nDof_iel_vec  +=  nDof_iel[c];    }


        KK_nonlocal_iel_iel.assign(nDof_iel_vec * nDof_iel_vec, 0.);   //resize
        KK_nonlocal_iel_jel.assign(nDof_iel_vec * nDof_jel_vec, 0.);   //resize
        KK_nonlocal_jel_iel.assign(nDof_jel_vec * nDof_iel_vec, 0.);   //resize
        KK_nonlocal_jel_jel.assign(nDof_jel_vec * nDof_jel_vec, 0.);   //resize
        Res_nonlocal_iel.assign(nDof_iel_vec, 0.);    //resize
        Res_nonlocal_jel.assign(nDof_jel_vec, 0.);    //resize

        Res_local_iel_mixed_num.assign(nDof_iel_vec, 0.);    //resize
        KK_local_iel_mixed_num.assign(nDof_iel_vec * nDof_iel_vec, 0.);

        if( check_if_same_elem(iel, jel) ) {
          Res_local_iel.assign(nDof_iel_vec, 0.);    //resize
          KK_local_iel.assign(nDof_iel_vec * nDof_iel_vec, 0.);
          if(integration_num_split != 0) {
//             Vectors and matrices for adaptive quadrature
            Res_local_iel_refined.assign(nDof_iel_vec, 0.);    //resize
            KK_local_iel_refined.assign(nDof_iel_vec * nDof_iel_vec, 0.);
          }
        }
// --- element matrix and vector resizing 


// --- solution        

        for (unsigned c = 0; c < n_components_ctrl; c++) {
           sol_ctrl_iel[c].resize(nDof_iel[c]);

           
        for(unsigned i = 0; i < sol_ctrl_iel[c].size(); i++) {
          unsigned iDof  = msh->GetSolutionDof(i, iel, solType_ctrl[c]);  // global to global mapping between coordinates node and coordinate dof
          sol_ctrl_iel[c][i] = (*sol->_Sol[solIndex_ctrl[c]])(iDof);  // global extraction and local storage for the element coordinates
          }
        }
// --- solution        

// ***************************************
// ******* iel-related stuff - END *************
// ***************************************
        


        
//------------ iface opening ---------        
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
// --- geom          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();
// --- geom          

   
	    if( ctrl::Gamma_control::face_is_a_Gamma_control_face(msh->el, iel, iface) ) {
//------------ iface opening ---------        
		
//                 count_visits_of_boundary_faces++;


             
     //**** Adaptive preparation - BEGIN ******** 
     std::vector < std::vector < double > >  Jac_kel_bdry_kqp_bdry;
     std::vector < std::vector < double > >  JacI_kel_bdry_kqp_bdry;
     double detJac_kel_bdry_kqp_bdry;
      vector < double >  x_kqp_bdry;
      double  weight_kqp_bdry;
      vector < double >  phi_coords_kel_bdry_kqp_bdry;
      vector < double >  phi_coords_x_kel_bdry_kqp_bdry;
      
      vector < double >  phi_ctrl_kel_bdry_kqp_bdry; /// @todo assume all ctrl components have the same FE family
      vector < double >  phi_ctrl_x_kel_bdry_kqp_bdry;
      
       std::vector< double > sol_ctrl_kqp_bdry(n_components_ctrl);
            
        for (unsigned c = 0; c < n_components_ctrl; c++) {   sol_ctrl_kqp_bdry[c] = 0.;   }


      
//         double weight3;
//         vector < double > phi3;

//        Evaluating coarse FE functions on Quadrature Points of the "sub-elements"
       std::vector < std::vector < std::vector <double > > > aP(3);  //[NFE_FAMS][DIM==3][N_DOFS]
        if(integration_num_split != 0) {
          for(unsigned fe_type = 0; fe_type < /*solType*/solType_coords + 1; fe_type++) { //loop up to the FE type + 1 of the unknown
            ProjectNodalToPolynomialCoefficients(aP[fe_type], geom_element_iel.get_coords_at_dofs_bdry_3d(), ielGeom_bdry, fe_type) ;         ///@todo check this!!!!!!!!!!!!
          }
        }                      
     //**** Adaptive preparation - END ********  
                      
              //Quadrature loop initialization - BEGIN
        const unsigned n_iqp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber();
        
        
       std::vector< double > sol_ctrl_iqp_bdry(n_components_ctrl);
        for (unsigned c = 0; c < n_components_ctrl; c++) {   sol_ctrl_iqp_bdry[c] = 0.;   }
              //Quadrature loop initialization - END
    
         
		for(unsigned iqp_bdry = 0; iqp_bdry < n_iqp_bdry; iqp_bdry++) {
            
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), iqp_bdry, Jac_iel_bdry_iqp_bdry, JacI_iel_bdry_iqp_bdry, detJac_iel_bdry_iqp_bdry, space_dim);
    
    weight_iqp_bdry = detJac_iel_bdry_iqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[iqp_bdry];

    elem_all[qrule_i][ielGeom_bdry][solType_coords] ->shape_funcs_current_elem(iqp_bdry, JacI_iel_bdry_iqp_bdry, phi_coords_iel_bdry_iqp_bdry, phi_coords_x_iel_bdry_iqp_bdry, boost::none, space_dim);

    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(iqp_bdry, JacI_iel_bdry_iqp_bdry, phi_ctrl_iel_bdry_iqp_bdry, phi_ctrl_x_iel_bdry_iqp_bdry, boost::none, space_dim);
            
//========== compute gauss quantities on the boundary - BEGIN ===============================================
//--- geom
          std::vector < double > x_iqp_bdry(dim, 0.);  ///@todo is this dim or dim_bdry?

            for(unsigned d = 0; d < x_iqp_bdry.size(); d++) {
	      for (int i_bdry = 0; i_bdry < geom_element_iel.get_coords_at_dofs_bdry_3d()[d].size(); i_bdry++)  {
			
              x_iqp_bdry[d] += geom_element_iel.get_coords_at_dofs_bdry_3d()[d][i_bdry] * phi_coords_iel_bdry_iqp_bdry[i_bdry];

		      }
            }
//--- geom
    
//--- solution
        for (unsigned c = 0; c < n_components_ctrl; c++) {
          sol_ctrl_iqp_bdry[c] = 0.;
	        for (int i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_ctrl_iqp_bdry[c] +=  /*sol_eldofs_Mat[pos_mat_ctrl]*/sol_ctrl_iel[c][i_vol] * phi_ctrl_iel_bdry_iqp_bdry[i_bdry];

		      }
       }
//--- solution
//========== compute gauss quantities on the boundary - END ================================================


  
      //============  Non-fractional assembly - BEGIN ==================
          
            if( check_if_same_elem_bdry(iel, jel, iface, jface) ) {
              
                
       //============  Mass assembly - BEGIN ==================
    for (unsigned c = 0; c < n_components_ctrl; c++) {
        
        integral += operator_L2 * alpha * weight_iqp_bdry * sol_ctrl_iqp_bdry[c] * sol_ctrl_iqp_bdry[c];
    }
  
    
    
    for (unsigned c = 0; c < n_components_ctrl; c++) {
          for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); l_bdry++) {
               		    unsigned int l_vol = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);

              double mass_res_i = phi_ctrl_iel_bdry_iqp_bdry[l_bdry] * sol_ctrl_iqp_bdry[c];
              const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_iel, c, l_vol);
              Res_local_iel[ res_pos/*l_vol*/ ] += operator_L2 * alpha * weight_iqp_bdry * mass_res_i ;
              Res_local_iel[ res_pos/*l_vol*/ ] += - rhs_one * weight_iqp_bdry * (phi_ctrl_iel_bdry_iqp_bdry[l_bdry] * (-1.) /** ( sin(2 * acos(0.0) * x1[0][l_bdry])) * ( sin(2 * acos(0.0) * x1[1][l_bdry]))*/);
              
              for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) { 
            for(unsigned m_bdry = 0; m_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); m_bdry++) {
               		    unsigned int m_vol = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
                const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol, m_vol);         
                KK_local_iel[ jac_pos/*l_vol * nDof_iel + m_vol*/ ] += operator_L2 * alpha * phi_ctrl_iel_bdry_iqp_bdry[l_bdry] * phi_ctrl_iel_bdry_iqp_bdry[m_bdry] * weight_iqp_bdry;
                 }
               }
             }
             
            } 
       }
        //============  Mass assembly - END ==================
         
      //============  Laplacian assembly - BEGIN ==================
      //============  Laplacian assembly - END ==================
            }
    
      //============  Non-fractional assembly - END ==================
             
             
      //============  Fractional assembly - BEGIN ==================
        if(operator_Hhalf != 0) {
                 
      //============ Same elem, && Adaptive quadrature - BEGIN ==================
        if( check_if_same_elem_bdry(iel, jel, iface, jface) && integration_num_split != 0) {
                
          /*const*/ short unsigned kelGeom_bdry = ielGeom_bdry;
           
          const unsigned n_kqp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_k, kelGeom_bdry).GetGaussPointsNumber();
                
                
          std::cout.precision(14);
          std::vector< std::vector< std::vector<double> > > x3;

          for(unsigned split = 0; split <= integration_num_split; split++) {

            
            if (dim_bdry/*dim*/ == 1) GetElementPartition1D(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, integration_num_split, x3, space_dim); //TODO space_dim or dim?
            else if (dim_bdry/*dim*/ == 2) {
              //TODO need to be checked !!!
              //GetElementPartition2D(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, integration_num_split, x3);
              GetElementPartitionQuad(x_iqp_bdry, geom_element_iel.get_coords_at_dofs_bdry_3d(), split, integration_num_split, x3);
            }

            //for(unsigned r = 0; r < size_part; r++) {
            for(unsigned r = 0; r < x3.size(); r++) {


              for(unsigned k_qp_bdry = 0; k_qp_bdry < n_kqp_bdry; k_qp_bdry++) {

// ********* PREPARATION PART - BEGIN ***************
                elem_all[qrule_k][kelGeom_bdry][solType_coords]->JacJacInv(x3[r]/*geom_element_iel.get_coords_at_dofs_bdry_3d()*/, k_qp_bdry, Jac_kel_bdry_kqp_bdry, JacI_kel_bdry_kqp_bdry, detJac_kel_bdry_kqp_bdry, space_dim);
    
//                 weight_kqp_bdry = detJac_kel_bdry_kqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_k, kelGeom_bdry).GetGaussWeightsPointer()[k_qp_bdry];
                
                msh->_finiteElement[kelGeom_bdry][solType_coords]->Jacobian(x3[r], k_qp_bdry, weight_kqp_bdry, phi_coords_kel_bdry_kqp_bdry, phi_x_placeholder_unused);
                
//                 elem_all[qrule_k][kelGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(k_qp_bdry, JacI_kel_bdry_kqp_bdry, phi_ctrl_kel_bdry_kqp_bdry, phi_ctrl_x_kel_bdry_kqp_bdry, boost::none, space_dim);
            
//                 elem_all[qrule_k][kelGeom_bdry][solType_coords] ->shape_funcs_current_elem(k_qp_bdry, JacI_kel_bdry_kqp_bdry, phi_coords_kel_bdry_kqp_bdry, phi_coords_x_kel_bdry_kqp_bdry, boost::none, space_dim);

//--- geom
                vector < double > x_kqp_bdry(dim, 0.);  ///@todo is this dim or dim_bdry?

                for(unsigned d = 0; d < x_kqp_bdry.size(); d++) {
                  for (int k_bdry = 0; k_bdry < /*geom_element_iel.get_coords_at_dofs_bdry_3d()[d].size()*/phi_coords_kel_bdry_kqp_bdry.size(); k_bdry++)  {
			
                    x_kqp_bdry[d] += x3[r][d][k_bdry]/*geom_element_iel.get_coords_at_dofs_bdry_3d()[d][i_vol]*/ * phi_coords_kel_bdry_kqp_bdry[k_bdry];

                  }
                }
//--- geom
    
                std::vector<double> xi3(dim, 0.);

                GetClosestPointInReferenceElement(geom_element_iel.get_coords_at_dofs_bdry_3d(), x_kqp_bdry, kelGeom_bdry, xi3);
//                 GetInverseMapping(solType_coords, kelGeom_bdry, aP, x_kqp_bdry, xi3, 1000);  ///@todo generalize to rectangular Jacobian TODO is needed?

                msh->_finiteElement[kelGeom_bdry][/*solType*/ solType_coords]->GetPhi(phi_ctrl_kel_bdry_kqp_bdry, xi3); //TODO solType or solType_coords?

                 std::vector<double> solY3(n_components_ctrl, 0.);

                for (unsigned c = 0; c < n_components_ctrl; c++) {
     for(unsigned i_bdry = 0; i_bdry < phi_ctrl_kel_bdry_kqp_bdry.size(); i_bdry++) {
                  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
                  solY3[c] += sol_ctrl_iel[c][i_vol] * phi_ctrl_kel_bdry_kqp_bdry[i_bdry];
                }
               }
// // // 
// // // 
// // //                     msh->_finiteElement[ielGeom][solType]->Jacobian(x3[r], k_qp_bdry, weight_kqp_bdry, phi3, phi_x_placeholder_unused);
// // // 
// // //                     vector < double > xg3(dim, 0.);
// // // 
// // //                     for(unsigned d = 0; d < dim; d++) {
// // //                       for(unsigned i = 0; i < nDof_iel; i++) {
// // //                         xg3[d] += x3[r][d][i] * phi3[i];
// // //                       }
// // //                     }
// // // 
// // //                     std::vector<double> xi3(dim, 0.);
// // // 
// // //                     GetClosestPointInReferenceElement(geom_element_iel.get_coords_at_dofs_bdry_3d(), xg3, ielGeom, xi3);
// // //                     GetInverseMapping(solType, ielGeom, aP, xg3, xi3, 1000);
// // // 
// // //                     msh->_finiteElement[ielGeom][solType]->GetPhi(phi3, xi3);
// // // 
// // //                     double solY3 = 0.;
// // //                     for(unsigned i = 0; i < nDof_iel; i++) {
// // //                       solY3 += sol_ctrl_iel[i] * phi3[i];
// // //                     }
// ********* PREPARATION PART - END ***************

// ********* BOUNDED PART - BEGIN ***************
// // // 
                double dist_xyz3 = 0;
                for(unsigned k = 0; k < x_iqp_bdry.size(); k++) {
                  dist_xyz3 += (x_iqp_bdry[k] - x_kqp_bdry[k]) * (x_iqp_bdry[k] - x_kqp_bdry[k]);
                }

                const double denom_ik = pow(dist_xyz3, (double)( 0.5 * dim_bdry + s_frac));
                
                const double common_weight =  0.5 * C_ns * operator_Hhalf * beta * check_limits * weight_iqp_bdry * weight_kqp_bdry  / denom_ik;


     for (unsigned c = 0; c < n_components_ctrl; c++) {
         integral    +=       common_weight * (sol_ctrl_iqp_bdry[c] - solY3[c]) * (sol_ctrl_iqp_bdry[c] - solY3[c]);
     }
     
                for (unsigned c = 0; c < n_components_ctrl; c++) {
               for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); l_bdry++) { //dofs of test function
                  unsigned int l_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);

              const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(nDof_iel, c, l_vol_iel);
                  Res_local_iel_refined[res_pos/* l_vol_iel*/ ]    +=      - common_weight * (sol_ctrl_iqp_bdry[c] - solY3[c]) * (phi_ctrl_iel_bdry_iqp_bdry[l_bdry] - phi_ctrl_kel_bdry_kqp_bdry[l_bdry]);

              for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) {
                for(unsigned m_bdry = 0; m_bdry < phi_ctrl_kel_bdry_kqp_bdry.size(); m_bdry++) { //dofs of unknown function
                    unsigned int m_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
                    
                    const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_iel, m_vol_iel);         
                KK_local_iel_refined[ jac_pos/*l_vol_iel * nDof_jel + m_vol_iel*/ ] += common_weight * (phi_ctrl_iel_bdry_iqp_bdry[m_bdry] - phi_ctrl_kel_bdry_kqp_bdry[m_bdry]) * 
                                                        (phi_ctrl_iel_bdry_iqp_bdry[l_bdry] - phi_ctrl_kel_bdry_kqp_bdry[l_bdry]);

                  }
                }
                  
              }
            }
     }
// ********* BOUNDED PART - END ***************
// ********* UNBOUNDED PART - BEGIN ***************
                if( iqp_bdry == integration_split_index ) { ///@todo is there a way to put this outside of the quadrature loop?
                    
              mixed_integral(unbounded,
                              dim,
                              dim_bdry,
                              ex_control,
                              N_div_unbounded,
                              weight_kqp_bdry,
                              x_kqp_bdry,
                              phi_ctrl_kel_bdry_kqp_bdry,
                              sol_ctrl_kqp_bdry,
                              s_frac,
                              check_limits,
                              C_ns,
                              operator_Hhalf,
                              beta,
                              nDof_jel, ///@todo
//                               KK_local_iel_refined,
//                               Res_local_iel_refined,
                              KK_local_iel,
                              Res_local_iel,
                              KK_local_iel_mixed_num,
                              Res_local_iel_mixed_num,
                              msh,
                              sol,
                              ml_sol,
                              iel,
                             jel,
                              iface,
                              bdry_bdry,
                              geom_element_jel,
                              jelGeom_bdry,
                              solType_coords,
                              integral
                             ); 
                }


// ********* UNBOUNDED PART - END ***************
                                         
                  } //end k_qp_bdry
                } //end r
              }  //end split
              
              
                
                
                
            }  //end iel == jel && integration_num_split != 0
      //============ Same elem, && Adaptive quadrature - END ==================
              
             
      //============ Either different elements, or lack of adaptivity (so all elements) - BEGIN ==================
        else {  //  if(iel != jel || integration_num_split == 0) 
            
// ********* UNBOUNDED PART - BEGIN ***************
//           if(check_if_same_elem_bdry(iel, jel, iface, jface)) { //TODO I removed this since we don't want iel==jel here
              
               mixed_integral(unbounded,
                              dim,
                              dim_bdry,
                              ex_control,
                              N_div_unbounded,
                              weight_iqp_bdry,
                              x_iqp_bdry,
                              phi_ctrl_iel_bdry_iqp_bdry,
                              sol_ctrl_iqp_bdry,
                              s_frac,
                              check_limits,
                              C_ns,
                              operator_Hhalf,
                              beta,
                              nDof_iel,
                              KK_local_iel,
                              Res_local_iel,
                              KK_local_iel_mixed_num,
                              Res_local_iel_mixed_num,
                              msh,
                              sol,
                              ml_sol,
                              iel,
                              jel,
                              iface,
                              bdry_bdry,
                              geom_element_jel,
                              jelGeom_bdry,
                              solType_coords,
                              integral
                             );
               
//           }
              
// ********* UNBOUNDED PART - END ***************
            
// ********* BOUNDED PART - BEGIN ***************
            for(unsigned jqp_bdry = 0; jqp_bdry < n_jqp_bdry; jqp_bdry++) {

              double dist_xyz = 0.;
              for(unsigned d = 0; d < x_iqp_bdry.size(); d++) {
                dist_xyz += (x_iqp_bdry[d] - x_jqp_bdry[jqp_bdry][d]) * (x_iqp_bdry[d] - x_jqp_bdry[jqp_bdry][d]);
              }

              const double denom = pow(dist_xyz, (double)(  0.5 * /*dim*/dim_bdry + s_frac));
              
              const double common_weight = (0.5 * C_ns) * operator_Hhalf * beta * check_limits * weight_iqp_bdry * weight_jqp_bdry[jqp_bdry]  / denom;

              for (unsigned c = 0; c < n_components_ctrl; c++) {
                  
                  integral +=  common_weight * (sol_ctrl_iqp_bdry[c] - sol_ctrl_jqp_bdry[c][jqp_bdry]) * sol_ctrl_iqp_bdry[c];
                  integral +=  common_weight * (sol_ctrl_iqp_bdry[c] - sol_ctrl_jqp_bdry[c][jqp_bdry]) * (- sol_ctrl_jqp_bdry[c][jqp_bdry]);
                  
              }              
              
              
              for (unsigned c = 0; c < n_components_ctrl; c++) {
                 for(unsigned l_bdry = 0; l_bdry < phi_ctrl_iel_bdry_iqp_bdry.size(); l_bdry++) { //dofs of test function
               		    unsigned int l_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, l_bdry);
               		    unsigned int l_vol_jel = msh->el->GetIG(jel_geommm, jface, l_bdry)/*msh->GetLocalFaceVertexIndex(jel, jface, l_bdry)*/;

              const unsigned res_pos_iel = assemble_jacobian<double,double>::res_row_index(nDof_iel, c, l_vol_iel);
              const unsigned res_pos_jel = assemble_jacobian<double,double>::res_row_index(nDof_jel, c, l_vol_jel);
                Res_nonlocal_iel[ res_pos_iel /*l_vol_iel*/ ]      +=      - common_weight * (sol_ctrl_iqp_bdry[c] - sol_ctrl_jqp_bdry[c][jqp_bdry]) * (phi_ctrl_iel_bdry_iqp_bdry[l_bdry]);

                Res_nonlocal_jel[ res_pos_jel /*l_vol_jel*/ ]      +=      - common_weight * (sol_ctrl_iqp_bdry[c] - sol_ctrl_jqp_bdry[c][jqp_bdry]) * (- phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry]);

               
//                 for(unsigned j = 0; j < nDof_jel; j++) {
                for (unsigned e = 0; e < n_components_ctrl; e++) {
                  if (e == c) {
                      for(unsigned m_bdry = 0; m_bdry < phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry].size(); m_bdry++) { //dofs of unknown function
               		    unsigned int m_vol_iel = msh->GetLocalFaceVertexIndex(iel, iface, m_bdry);
               		    unsigned int m_vol_jel = msh->el->GetIG(jel_geommm, jface, m_bdry)/*msh->GetLocalFaceVertexIndex(jel, jface, m_bdry)*/;

                    const unsigned jac_pos_iel_iel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_iel, m_vol_iel);         
                    const unsigned jac_pos_iel_jel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_iel, m_vol_jel);
                    const unsigned jac_pos_jel_iel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_jel, m_vol_iel);
                    const unsigned jac_pos_jel_jel = assemble_jacobian< double, double >::jac_row_col_index(nDof_iel, nDof_iel_vec, c, e, l_vol_jel, m_vol_jel);
                    
             /*  u(x) v(x)*/     KK_nonlocal_iel_iel[ jac_pos_iel_iel /*l_vol_iel * nDof_jel + m_vol_iel*/ ] += common_weight *          phi_ctrl_iel_bdry_iqp_bdry[m_bdry]            *    phi_ctrl_iel_bdry_iqp_bdry[l_bdry];

             /*- u(y) v(x)*/     KK_nonlocal_iel_jel[ jac_pos_iel_jel /*l_vol_iel * nDof_jel + m_vol_jel*/ ] += common_weight * (- 1.) * phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][m_bdry]  *    phi_ctrl_iel_bdry_iqp_bdry[l_bdry];

             /*- u(x) v(y)*/     KK_nonlocal_jel_iel[ jac_pos_jel_iel /*l_vol_jel * nDof_jel + m_vol_iel*/ ] += common_weight * (- 1.) * phi_ctrl_iel_bdry_iqp_bdry[m_bdry]            *   phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry];

             /*  u(y) v(y)*/     KK_nonlocal_jel_jel[ jac_pos_jel_jel /*l_vol_jel * nDof_jel + m_vol_jel*/ ] += common_weight *          phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][m_bdry]  *    phi_ctrl_jel_bdry_jqp_bdry[jqp_bdry][l_bdry];


                     }
                   }
                 }
                }
                
               }
               
              } //endl jqp_bdry loop
// ********* BOUNDED PART - END ***************
            
            
         } //end if(iel != jel || integration_num_split == 0)
      //============ Either different elements, or lack of adaptivity (so all elements) - END ==================

        
         } //operator_Hhalf != 0              
      //============  Fractional assembly - END ==================
            
      }   //end iqp_bdry
      
//                std::vector<unsigned> Sol_n_el_dofs_Mat_vol2(1, nDof_jel);
//          assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
//          assemble_jacobian<double,double>::print_element_jacobian(iel, KK_local_iel_mixed_num, Sol_n_el_dofs_Mat_vol2, 10, 5);
      
      
              
//----- iface ---        
        } //end face for control
        
      } //end iface
//----- iface ---        




                
//============ add to global - BEGIN ==================
// // multiply everything by -1.? Don't think so
// std::transform(KK_local_iel.begin(), KK_local_iel.end(), KK_local_iel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_local_iel.begin(), Res_local_iel.end(), Res_local_iel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(KK_local_iel_refined.begin(), KK_local_iel_refined.end(), KK_local_iel_refined.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_local_iel_refined.begin(), Res_local_iel_refined.end(), Res_local_iel_refined.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(KK_local_iel_mixed_num.begin(), KK_local_iel_mixed_num.end(), KK_local_iel_mixed_num.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_local_iel_mixed_num.begin(), Res_local_iel_mixed_num.end(), Res_local_iel_mixed_num.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(KK_nonlocal_iel_iel.begin(), KK_nonlocal_iel_iel.end(), KK_nonlocal_iel_iel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(KK_nonlocal_iel_jel.begin(), KK_nonlocal_iel_jel.end(), KK_nonlocal_iel_jel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(KK_nonlocal_jel_iel.begin(), KK_nonlocal_jel_iel.end(), KK_nonlocal_jel_iel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(KK_nonlocal_jel_jel.begin(), KK_nonlocal_jel_jel.end(), KK_nonlocal_jel_jel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// 
// std::transform(Res_nonlocal_iel.begin(), Res_nonlocal_iel.end(), Res_nonlocal_iel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// std::transform(Res_nonlocal_jel.begin(), Res_nonlocal_jel.end(), Res_nonlocal_jel.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, -1.));
// // multiply everything by -1.


 if( check_if_same_elem(iel, jel)/*check_if_same_elem_bdry(iel, jel, iface, jface) if you put it inside the face loops */ ) {
     

         RES->add_vector_blocked(Res_local_iel, l2gMap_iel_vec);
          KK->add_matrix_blocked(KK_local_iel, l2gMap_iel_vec, l2gMap_iel_vec);
        
        
          if(integration_num_split != 0) {
            RES->add_vector_blocked(Res_local_iel_refined, l2gMap_iel_vec);
            KK->add_matrix_blocked(KK_local_iel_refined, l2gMap_iel_vec, l2gMap_iel_vec);
          }
          
       }

        RES->add_vector_blocked(Res_local_iel_mixed_num, l2gMap_iel_vec);
        KK->add_matrix_blocked(KK_local_iel_mixed_num, l2gMap_iel_vec, l2gMap_iel_vec);
        

        RES->add_vector_blocked(Res_nonlocal_iel, l2gMap_iel_vec);
        RES->add_vector_blocked(Res_nonlocal_jel, l2gMap_jel_vec);
        
        KK->add_matrix_blocked(KK_nonlocal_iel_iel, l2gMap_iel_vec, l2gMap_iel_vec);
        KK->add_matrix_blocked(KK_nonlocal_iel_jel, l2gMap_iel_vec, l2gMap_jel_vec);
        KK->add_matrix_blocked(KK_nonlocal_jel_iel, l2gMap_jel_vec, l2gMap_iel_vec);
        KK->add_matrix_blocked(KK_nonlocal_jel_jel, l2gMap_jel_vec, l2gMap_jel_vec);
// Since 1 is dense and 3 are sparse, and the dense dofs are 30, we should have at most 3x9 + 30 = 57, but in the sparsity print it shows 30. That's the problem

//============ add to global - END ==================

        
          
// //              if (print_algebra_local) {
//          std::vector<unsigned> Sol_n_el_dofs_Mat_vol2(1, nDof_jel);
// //          assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
//          assemble_jacobian<double,double>::print_element_jacobian(iel, KK_local_iel_mixed_num, Sol_n_el_dofs_Mat_vol2, 10, 5);
// //      }
         

        
//----- iel ---        
    } //end control elem flag i (control_flag_iel == 1)
  } //end iel
//----- iel ---        


//----- jface ---        
     
} //end face for control


      } //end jface

 //----- jface ---   
 
     
    }  //end control elem flag jel
    



   } //end jel
//----- jel ---        
   
 } //end kproc
    
    
    // integral - BEGIN ************
  std::cout << std::endl;
  std::cout <<  "&&&&&& Integrals from previous iteration &&&&&&&&&&" << std::endl;
  double integral_parallel = 0.; MPI_Allreduce( &integral, &integral_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  std::cout << "The value of the integral  is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_parallel << std::endl;
  std::cout <<  "&&&&&& Integrals from previous iteration &&&&&&&&&&" << std::endl;
  std::cout << std::endl;
                  
       std::cout << "SQUARE ROOOOOOOOTTTTTTTTTTTTTTTTTTTT AFTER" << std::endl;
       
// integral - END ************

    
    
    
    
    
  
  }
    
  //**********FRAC CONTROL - END *****************************************
  

//*********************** Mesh independent, ALMOST - END *****************************************
    
}

 

 

  
namespace Gamma_control_equation_integer {
  
   
 
  void control_eqn_bdry_integer_sobolev_differentiability_index(const unsigned iproc,
                        MultiLevelProblem &    ml_prob,
                        MultiLevelSolution*    ml_sol,
                        const Solution*        sol,
                        const Mesh * msh,
                        //----- Eqn ------
                        const  LinearEquationSolver* pdeSys,
                        //----- Geom Element ------
                        CurrentElem < double > & geom_element_iel,
                        const unsigned int solType_coords,
                        const unsigned int space_dim,
                        //----- Mat ------
                        const unsigned int n_unknowns,
                        const    vector < std::string > & Solname_Mat,
                        const    vector < unsigned > & SolFEType_Mat,
                        const vector < unsigned > & SolIndex_Mat,
                        const vector < unsigned > & SolPdeIndex,
                        vector < unsigned > & Sol_n_el_dofs_Mat, 
                        vector < vector < double > > & sol_eldofs_Mat,  
                        vector < vector < int > > & L2G_dofmap_Mat,
                        //--- Equation, local --------
                        const unsigned max_size,
                        //----- Sol ------
                        const unsigned int n_quantities,
                        vector < unsigned > SolFEType_quantities,
                        //---- Quadrature - FE Evaluations -------
                        std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all,
                        //---- Quadrature, Geometry ------
                        std::vector < std::vector < double > >  Jac_iel_bdry_iqp_bdry,
                        std::vector < std::vector < double > >  JacI_iel_bdry_iqp_bdry,
                        double detJac_iel_bdry_iqp_bdry,
                        double weight_iqp_bdry,
                        //---- Quadrature, Control ------
                        vector <double> phi_ctrl_iel_bdry_iqp_bdry,
                        vector <double> phi_ctrl_x_iel_bdry_iqp_bdry, 
                        //---- Control -------
                        const unsigned int n_components_ctrl,
                        const unsigned int pos_mat_ctrl,
                        const unsigned int pos_sol_ctrl,
                        //-----------
                        const unsigned int is_block_dctrl_ctrl_inside_main_big_assembly,
                        //-----------
                        SparseMatrix *  KK,
                        NumericVector * RES,
                        const bool assembleMatrix,
                        //-----------
                        const double alpha,
                        const double beta,
                        const double rhs_one,
                        //-----------
                        const unsigned qrule_i,
                        //-----------
                        const bool print_algebra_local
                       ) {
      
   
   const unsigned  elem_dof_size_max = n_components_ctrl * max_size;
   
   std::vector < double >  Res_ctrl_only;      Res_ctrl_only.reserve( elem_dof_size_max );                         //should have Mat order
   std::vector < double >  Jac_ctrl_only;      Jac_ctrl_only.reserve( elem_dof_size_max * elem_dof_size_max);   //should have Mat order

   
// integral - BEGIN ************
  double integral_alpha  = 0.;
  double integral_beta   = 0.;
// integral - END ************

 
    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
        
// --- geometry        
            geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);

            const short unsigned ielGeom = geom_element_iel.geom_type();

           geom_element_iel.set_elem_center_3d(iel, solType_coords);
// --- geometry        

           
 //***************************************************
   el_dofs_unknowns_vol(sol, msh, pdeSys, iel,
                        SolFEType_Mat,
                        SolIndex_Mat,
                        SolPdeIndex,
                        Sol_n_el_dofs_Mat, 
                        sol_eldofs_Mat,  
                        L2G_dofmap_Mat);  //all unknowns here, perhaps we could restrict it to the ctrl components only
  //***************************************************
      
 //***************************************************
 //***************************************************
   //extract a subvector containing only the control components, starting from zero      
      
   std::vector< unsigned >::const_iterator first  = Sol_n_el_dofs_Mat.begin() + pos_mat_ctrl;
   std::vector< unsigned >::const_iterator last   = Sol_n_el_dofs_Mat.begin() + pos_mat_ctrl + n_components_ctrl;
   std::vector< unsigned > Sol_n_el_dofs_Mat_ctrl_only(first, last);
   
   unsigned int sum_Sol_n_el_dofs_ctrl_only = ElementJacRes< double >::compute_sum_n_dofs(Sol_n_el_dofs_Mat_ctrl_only);
 //***************************************************
 //***************************************************
   //create a subvector containing only the control components, starting from zero  
   std::vector< unsigned > L2G_dofmap_Mat_ctrl_only; 
   L2G_dofmap_Mat_ctrl_only.resize(0);
      for (unsigned  k = 0; k < n_components_ctrl; k++)     L2G_dofmap_Mat_ctrl_only.insert(L2G_dofmap_Mat_ctrl_only.end(), L2G_dofmap_Mat[pos_mat_ctrl + k].begin(), L2G_dofmap_Mat[pos_mat_ctrl + k].end());
 //***************************************************
 //***************************************************

   
 //***************************************************
    const unsigned int res_length =  n_components_ctrl * Sol_n_el_dofs_Mat[pos_mat_ctrl];

    Res_ctrl_only.resize(res_length);                 std::fill(Res_ctrl_only.begin(), Res_ctrl_only.end(), 0.);

    Jac_ctrl_only.resize(res_length * res_length);    std::fill(Jac_ctrl_only.begin(), Jac_ctrl_only.end(), 0.);
 //***************************************************


    if ( ctrl::Gamma_control::volume_elem_contains_a_Gamma_control_face(geom_element_iel.get_elem_center_3d()) ) {
        
  
 //************ set control flag *********************
  std::vector< std::vector< int > > control_node_flag_iel_all_faces = 
       ctrl::Gamma_control::is_dof_associated_to_Gamma_control_equation(msh, ml_sol, & ml_prob, iel, geom_element_iel, solType_coords, Solname_Mat, SolFEType_Mat, Sol_n_el_dofs_Mat, pos_mat_ctrl, n_components_ctrl);
       
       ///@todo here I have to do it "on the go", for each boundary dof!!!
  //*************************************************** 
      
	  
	  // loop on faces of the current element
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
// ---
       std::vector<unsigned int> Sol_el_n_dofs_current_face(n_quantities); ///@todo the active flag is not an unknown!

       for (unsigned  k = 0; k < Sol_el_n_dofs_current_face.size(); k++) {
                 if (SolFEType_quantities[k] < 3) Sol_el_n_dofs_current_face[k] = msh->GetElementFaceDofNumber(iel, iface, SolFEType_quantities[k]);  ///@todo fix this absence
       }
       
// ---     
       
// --- geometry        
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

// --- geometry        
         
	    if( ctrl::Gamma_control::face_is_a_Gamma_control_face(msh->el, iel, iface) ) {
              

//========= initialize gauss quantities on the boundary ============================================
                std::vector< double > sol_ctrl_iqp_bdry(n_components_ctrl);
                std::vector< std::vector< double > > sol_ctrl_x_iqp_bdry(n_components_ctrl);
                
           for (unsigned c = 0; c < n_components_ctrl; c++) {
               sol_ctrl_iqp_bdry[c] = 0.;                
               sol_ctrl_x_iqp_bdry[c].assign(space_dim, 0.); 
            }
//========= initialize gauss quantities on the boundary ============================================
		
        const unsigned n_qp_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber();
        
    
		for(unsigned iqp_bdry = 0; iqp_bdry < n_qp_bdry; iqp_bdry++) {
    
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), iqp_bdry, Jac_iel_bdry_iqp_bdry, JacI_iel_bdry_iqp_bdry, detJac_iel_bdry_iqp_bdry, space_dim);
// 	elem_all[qrule_i][ielGeom_bdry][solType_coords]->compute_normal(Jac_qp_bdry, normal);
    
    weight_iqp_bdry = detJac_iel_bdry_iqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[iqp_bdry];

    
   for (unsigned c = 0; c < n_components_ctrl; c++) {
    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl + c]] ->shape_funcs_current_elem(iqp_bdry, JacI_iel_bdry_iqp_bdry, phi_ctrl_iel_bdry_iqp_bdry, phi_ctrl_x_iel_bdry_iqp_bdry, boost::none, space_dim);
   }  
   
//========== compute gauss quantities on the boundary ===============================================
   for (unsigned c = 0; c < n_components_ctrl; c++) {
          sol_ctrl_iqp_bdry[c] = 0.;
                  std::fill(sol_ctrl_x_iqp_bdry[c].begin(), sol_ctrl_x_iqp_bdry[c].end(), 0.);
		      for (int i_bdry = 0; i_bdry < phi_ctrl_iel_bdry_iqp_bdry.size()/*Sol_n_el_dofs_quantities[pos_sol_ctrl + c]*/; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_ctrl_iqp_bdry[c] +=  sol_eldofs_Mat[pos_mat_ctrl + c][i_vol] * phi_ctrl_iel_bdry_iqp_bdry[i_bdry];
                            for (int d = 0; d < space_dim; d++) {
			      sol_ctrl_x_iqp_bdry[c][d] += sol_eldofs_Mat[pos_mat_ctrl + c][i_vol] * phi_ctrl_x_iel_bdry_iqp_bdry[i_bdry * space_dim + d];
			    }
		      }
   }
//========== compute gauss quantities on the boundary ================================================

// integral - BEGIN ************
   for (unsigned c = 0; c < n_components_ctrl; c++) {
                  integral_alpha +=  weight_iqp_bdry * sol_ctrl_iqp_bdry[c] * sol_ctrl_iqp_bdry[c]; 
                 for (unsigned d = 0; d < space_dim; d++) {
                  integral_beta  +=  weight_iqp_bdry * sol_ctrl_x_iqp_bdry[c][d] * sol_ctrl_x_iqp_bdry[c][d];
      }
   }
  // integral - END ************


//equation - BEGIN ************
   for (unsigned c = 0; c < n_components_ctrl; c++) {

		  // *** phi_i loop ***
		  for(unsigned i_bdry = 0; i_bdry < Sol_el_n_dofs_current_face[pos_sol_ctrl + c] ; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);

                 double lap_rhs_dctrl_ctrl_bdry_gss_i_c = 0.;
                 for (unsigned d = 0; d < space_dim; d++) {
                       /*if ( i_vol < Sol_n_el_dofs_Mat[pos_mat_ctrl] )*/  lap_rhs_dctrl_ctrl_bdry_gss_i_c +=  phi_ctrl_x_iel_bdry_iqp_bdry[i_bdry * space_dim + d] * sol_ctrl_x_iqp_bdry[c][d];
                 }
                 
		 
//============ Bdry Residuals - BEGIN ==================	
           const unsigned res_pos = assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_ctrl_only, /*pos_mat_ctrl +*/ c, i_vol);
           
                Res_ctrl_only[ res_pos ]  +=  - control_node_flag_iel_all_faces[c][i_vol] *  weight_iqp_bdry *
                                         (     ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * alpha * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * sol_ctrl_iqp_bdry[c]
							                +  ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * beta  * lap_rhs_dctrl_ctrl_bdry_gss_i_c
							                         );  //boundary optimality condition
                Res_ctrl_only[ res_pos ] += - rhs_one * weight_iqp_bdry * (phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * (-1.) /** ( sin(2 * acos(0.0) * x1[0][l_bdry])) * ( sin(2 * acos(0.0) * x1[1][l_bdry]))*/);
                          
//============ Bdry Residuals - END ==================    
		    
//============ Bdry Jacobians - BEGIN ==================	
   for (unsigned e = 0; e < n_components_ctrl; e++) {
        if (e == c) {
		    for(unsigned j_bdry = 0; j_bdry < Sol_el_n_dofs_current_face[pos_sol_ctrl + e]; j_bdry++) {
		         unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);

//========block delta_control / control ========
              double  lap_mat_dctrl_ctrl_bdry_gss = 0.;
		      for (unsigned d = 0; d < space_dim; d++) {  
                  lap_mat_dctrl_ctrl_bdry_gss += phi_ctrl_x_iel_bdry_iqp_bdry[i_bdry * space_dim + d] * phi_ctrl_x_iel_bdry_iqp_bdry[j_bdry * space_dim + d];    
               }

          
           const unsigned jac_pos = assemble_jacobian< double, double >::jac_row_col_index(Sol_n_el_dofs_Mat_ctrl_only, sum_Sol_n_el_dofs_ctrl_only, /*pos_mat_ctrl +*/ c, /*pos_mat_ctrl +*/ e, i_vol, j_vol);
              Jac_ctrl_only[ jac_pos ]   +=  control_node_flag_iel_all_faces[c][i_vol] *  weight_iqp_bdry * (
                                    ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * alpha * phi_ctrl_iel_bdry_iqp_bdry[i_bdry] * phi_ctrl_iel_bdry_iqp_bdry[j_bdry] 
			                      + ( 1 - is_block_dctrl_ctrl_inside_main_big_assembly ) * beta  * lap_mat_dctrl_ctrl_bdry_gss);   
				
	        }   //end j loop
	      }
   }
//============ Bdry Jacobians - END ==================	
	      
		  }  //end i loop
		  
   }  //end components ctrl	  
//equation - END ************
		  
		}  //end iqp_bdry loop
	  }    //end if control face
	      
	  }    //end loop over faces
	  
	} //end if control element flag

	
  /*is_res_control_only*/
                     RES->add_vector_blocked(Res_ctrl_only, L2G_dofmap_Mat_ctrl_only);
              
                  if (assembleMatrix) {
                     KK->add_matrix_blocked(Jac_ctrl_only, L2G_dofmap_Mat_ctrl_only, L2G_dofmap_Mat_ctrl_only);
                  }   
         
         
         
        if (print_algebra_local) {
  
         assemble_jacobian<double,double>::print_element_residual(iel, Res_ctrl_only, Sol_n_el_dofs_Mat_ctrl_only, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac_ctrl_only, Sol_n_el_dofs_Mat_ctrl_only, 10, 5);
        }
     
    }  //iel

    
// integral - BEGIN ************
  std::cout << std::endl;
  std::cout <<  "&&&&&& Integrals from previous iteration &&&&&&&&&&" << std::endl;
  double integral_alpha_parallel = 0.; MPI_Allreduce( &integral_alpha, &integral_alpha_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_alpha_parallel << std::endl;
  double integral_beta_parallel = 0.;  MPI_Allreduce( &integral_beta, &integral_beta_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_beta_parallel << std::endl;
  std::cout <<  "&&&&&& Integrals from previous iteration &&&&&&&&&&" << std::endl;
  std::cout << std::endl;
// integral - END ************
   
    
    
}




}
 
 

 
 

 namespace cost_functional {

 
//************** how to retrieve theta from proc0 ************************************* 
const double get_theta_value(const unsigned int nprocs, const Solution * sol, const unsigned int sol_theta_index) {
    
NumericVector* local_theta_vec;

          local_theta_vec = NumericVector::build().release();
        local_theta_vec->init(*sol->_Sol[sol_theta_index]);
        sol->_Sol[sol_theta_index]->localize(*local_theta_vec);
        
        PetscScalar* values;
        VecGetArray(static_cast<PetscVector*>(local_theta_vec)->vec(), &values);
        double theta_value = values[0];
        if (nprocs == 1) {
            if ( (*sol->_Sol[sol_theta_index])(0) != theta_value) abort();
        }
        
        return theta_value;
}
//*************************************************** 


  
  

void compute_cost_functional_regularization_bdry_vec(const MultiLevelProblem& ml_prob, 
                     const unsigned level,
                     const unsigned iteration,
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars  
) {

    
    
  const double cost_functional_coeff = COST_FUNCTIONAL_COEFF;
  const double alpha = ALPHA_CTRL_BDRY;
  const double beta  = BETA_CTRL_BDRY;

  
  
  const Mesh*          msh          	= ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         	= msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  ml_sol    = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        	= ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object
  
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  unsigned solType_coords = coordXType;  //FE_DOMAIN = 0; //we do linear FE this time // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = /*dim*/  3  /*2*/    ;
 
  std::vector<double> normal_iqp(dim_offset_grad /*space_dim*/, 0.);

  vector < vector < double > > coordX(dim);    // local coordinates
  vector< vector < double> > coordX_bd(dim);

  for (unsigned  k = 0; k < dim; k++) { 
        coordX[k].reserve(max_size);
        coordX_bd[k].reserve(max_size); 
  }
  
  double AbsDetJxWeight_iqp = 0.;
  double AbsDetJxWeight_iqp_bdry = 0.;
  
  
  //geometry *******************************

//STATE######################################################################
  //velocity *******************************
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = ml_sol->GetIndex("u_0");
  solVIndex[1] = ml_sol->GetIndex("u_1");

  if (dim == 3) solVIndex[2] = ml_sol->GetIndex("u_2");

  unsigned solVType = ml_sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
  
  vector < vector < double > >  solV(dim);    // local solution
  vector <double >  V_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(max_size);
  }

  
  vector <double> phiV_gss;  // local test function
  vector <double> phiV_x_gss; // local test function first order partial derivatives

  phiV_gss.reserve(max_size);
  phiV_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
   
//STATE######################################################################
  

//CONTROL_@bdry######################################################################
  vector < unsigned > solVctrlIndex(dim);
  solVctrlIndex[0] = ml_sol->GetIndex("ctrl_0");
  solVctrlIndex[1] = ml_sol->GetIndex("ctrl_1");
  if (dim == 3) solVctrlIndex[2] = ml_sol->GetIndex("ctrl_2");

  unsigned solVctrlType = ml_sol->GetSolutionType(solVctrlIndex[0]);    // get the finite element type for "u"
  
  vector < vector < double > >  solVctrl(dim);    // local solution
  vector < double >   Vctrl_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solVctrl[k].reserve(max_size);
  }

  
  vector <double> phiVctrl_gss_bd;  // local test function
  vector <double> phiVctrl_x_gss_bd; // local test function first order partial derivatives

  phiVctrl_gss_bd.reserve(max_size);
  phiVctrl_x_gss_bd.reserve(max_size * dim_offset_grad );
  
//CONTROL_@bdry######################################################################
  
//Theta value ######################################################################
   const unsigned solThetaIndex = ml_sol->GetIndex("theta");
   const unsigned solThetaType = ml_sol->GetSolutionType(solThetaIndex);
   
//    double solTheta = (*sol->_Sol[solThetaIndex])(0)/*0.*/;
   //************** how to retrieve theta from proc0 ************************************* 
 double solTheta = ctrl::cost_functional::get_theta_value(msh->n_processors(), sol, solThetaIndex);
//*************************************************** 
// 		     solTheta = (*sol->_Sol[solThetaIndex])(0);
//Theta value ######################################################################


// Vel_desired##################################################################
  vector <double> phiVdes_gss;  // local test function
  vector <double> phiVdes_x_gss; // local test function first order partial derivatives

  phiVdes_gss.reserve(max_size);
  phiVdes_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);

//   vector< vector < double > >  solVdes(dim);    // local solution
  vector <double>  solVdes(dim,0.);
  vector<double> Vdes_gss(dim, 0.);  
  
//  for (unsigned  k = 0; k < dim; k++) {
//     solVdes[k].reserve(max_size);
//   }
//   


// Vel_desired##################################################################

  
vector<double> integral(dim);

double  integral_target_alpha = 0.;

double	integral_beta   = 0.;
double	integral_gamma  = 0.;

double integral_g_dot_n = 0.;
  

// double	integral_div_ctrl  = 0.;

//*************************************************** 
  //--- quadrature rules -------------------
  constexpr unsigned qrule_i = QRULE_I;
  
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
    double detJac_iqp;

     std::vector < std::vector < double > >  JacI_iqp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_iqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_iqp_bdry.size(); d++) {   Jac_iqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp_bdry.size(); d++) { JacI_iqp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
    
      //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all;
  ml_prob.get_all_abstract_fe_multiple(elem_all);
//*************************************************** 

  
  

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

   // geometry *****************************
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element_iel.geom_type();
  // geometry end *****************************

// equation
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
//     unsigned nDofsVdes = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsVctrl = msh->GetElementDofNumber(iel, solVctrlType);    // number of solution element dofs
    
    unsigned nDofsThetactrl = msh->GetElementDofNumber(iel,solThetaType);
    
     
    for (unsigned  k = 0; k < dim; k++)  {
      solV[k].resize(nDofsV);
      solVctrl[k].resize(nDofsVctrl);
//       solVdes[k].resize(nDofsVdes);
    }
  //*************************************** 
  
  //***** set target domain flag ********************************** 
  geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = femus::ctrl::cost_functional::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
//***************************************       
    
    
 //STATE###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
      }
    }
//STATE###################################################################

//CONTROL###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVctrlDof = msh->GetSolutionDof(i, iel, solVctrlType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solVctrl[k][i] = (*sol->_Sol[solVctrlIndex[k]])(solVctrlDof);      // global extraction and local storage for the solution
      }
    }
//CONTROL###################################################################




  //DESIRED VEL###################################################################  
    // velocity ************
//     for (unsigned i = 0; i < nDofsV; i++) {
//       unsigned solVdesDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < solVdes.size() /*dim*/; k++) {
        solVdes[k]/*[i]*/ = femus::ctrl::cost_functional::DesiredTargetVec()[k] /*(*sol->_Sol[solVIndex[k]])(solVdesDof)*/;      // global extraction and local storage for the solution
     }
//     }
 //DESIRED VEL###################################################################

 

//========BoundaryLoop=====================================================================

	if ( ctrl::Gamma_control::volume_elem_contains_a_Gamma_control_face( geom_element_iel.get_elem_center_3d() ) ) {

        
    vector<double> normal_iqp(dim_offset_grad /*space_dim*/, 0.);
	  
    for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {

       const unsigned ielGeom_bd = msh->GetElementFaceType(iel, jface);    
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

	    // look for boundary faces
	    if( ctrl::Gamma_control::face_is_a_Gamma_control_face(msh->el, iel, jface) ) {
	  
//=================================================== 
		
//========= initialize gauss quantities on the boundary ============================================
    vector < double >   Vctrl_bd_qp(dim, 0.);    //  solution@bdry
    vector < vector < double > > gradVctrl_bd_qp(dim);
      for (unsigned  k = 0; k < dim; k++) {
          gradVctrl_bd_qp[k].resize(dim_offset_grad /*space_dim*/);
          std::fill(gradVctrl_bd_qp[k].begin(), gradVctrl_bd_qp[k].end(), 0);
        }

//========= gauss_loop boundary===============================================================
	    for(unsigned iqp_bdry=0; iqp_bdry < ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussPointsNumber(); iqp_bdry++) {
    elem_all[qrule_i][ielGeom_bd][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), iqp_bdry, Jac_iqp_bdry, JacI_iqp_bdry, detJac_iqp_bdry, space_dim);
	elem_all[qrule_i][ielGeom_bd][solType_coords]->compute_normal(Jac_iqp_bdry, normal_iqp);
    
    AbsDetJxWeight_iqp_bdry = detJac_iqp_bdry * ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussWeightsPointer()[iqp_bdry];

    elem_all[qrule_i][ielGeom_bd][solVctrlType] ->shape_funcs_current_elem(iqp_bdry, JacI_iqp_bdry,phiVctrl_gss_bd,phiVctrl_x_gss_bd , boost::none , space_dim);

     
//========== compute gauss quantities on the boundary ===============================================
    for (unsigned  k = 0; k < dim; k++) {
	  Vctrl_bd_qp[k] = 0.;
	  for(unsigned ivar2 = 0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++) { gradVctrl_bd_qp[k][ivar2] = 0.; }
	  
	  for (unsigned i = 0; i < nDofsVctrl; i++) {
                 const unsigned ndof_bdry = msh->GetElementFaceDofNumber(iel, jface, solVctrlType);

		   for(int i_bd = 0; i_bd < ndof_bdry; i_bd++) {
		       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
		       Vctrl_bd_qp[k] += phiVctrl_gss_bd[i_bd] * solVctrl[k][i_vol];
		       for(unsigned ivar2 = 0; ivar2 < dim_offset_grad; ivar2++) {
			   gradVctrl_bd_qp[k][ivar2] += phiVctrl_x_gss_bd[i_bd * dim_offset_grad + ivar2 ] * solVctrl[k][i_vol]; 
		         }
		   }
      }
    }
 //end unknowns eval at gauss points ********************************
		  
//========== compute gauss quantities on the boundary ================================================
      for (unsigned  k = 0; k < dim; k++) {
	 integral_beta	+= ((Vctrl_bd_qp[k])*(Vctrl_bd_qp[k]) * AbsDetJxWeight_iqp_bdry);
	 integral_g_dot_n += Vctrl_bd_qp[k]*normal_iqp[k] * AbsDetJxWeight_iqp_bdry;
      }
      for (unsigned  k = 0; k < dim; k++) {
	for (unsigned  j = 0; j < dim; j++) {	
		integral_gamma	  += (gradVctrl_bd_qp[k][j])*(gradVctrl_bd_qp[k][j]) * AbsDetJxWeight_iqp_bdry;
	}
      }


                }  //end iqp_bdryry loop
	  
             }    //end if control face

        
    }  // loop over element faces //jface  
      
  } //end if control element flag

  
  
      // *** Gauss point loop ***
      for (unsigned iqp = 0; iqp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); iqp++) {
//STATE######## VolumeLoop #####################################################################	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), iqp, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);
    AbsDetJxWeight_iqp = detJac_iqp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[iqp];
   
    elem_all[qrule_i][ielGeom][solVType]->shape_funcs_current_elem(iqp, JacI_iqp, phiV_gss, phiV_x_gss, boost::none , space_dim);
    elem_all[qrule_i][ielGeom][solVType /*solVdes*/]->shape_funcs_current_elem(iqp, JacI_iqp, phiVdes_gss, phiVdes_x_gss, boost::none , space_dim);

    
      for (unsigned  k = 0; k < dim; k++) {
	           V_gss[k] = 0.;
	           Vdes_gss[k] = 0.;
	    for (unsigned i = 0; i < nDofsV; i++) {
	   	   V_gss[k] += solV[k][i] * phiV_gss[i];
		   Vdes_gss[k] += solVdes[k]/*[i]*/ * phiVdes_gss[i];
	    }
	  }
	
//       for (unsigned  k = 0; k < dim; k++) {
//           integral_div_ctrl +=  AbsDetJxWeight_iqp * gradVctrl_gss[k][k] /** phiVctrl_gss[i]*/;
//       }

      for (unsigned  k = 0; k < dim; k++) {
	      integral_target_alpha += (( target_flag ) *((V_gss[k]  - Vdes_gss[k]) * (V_gss[k]  - Vdes_gss[k])) * AbsDetJxWeight_iqp);
      }
      
      }// end gauss point loop
    } //end element loop  

       std::ostringstream filename_out; filename_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "Integral_computation"  << ".txt";

       std::ofstream intgr_fstream;
  if (paral::get_rank() == 0 ) {
      
      intgr_fstream.open(filename_out.str().c_str(),std::ios_base::app);
      
      intgr_fstream << " ***************************** Iteration "<< iteration << " *********************************** " <<  std::endl << std::endl;
      intgr_fstream << "The value of the target functional for " << "alpha " <<   std::setprecision(0) << std::scientific << cost_functional_coeff << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_target_alpha << std::endl;
      intgr_fstream << "The value of the L2 control for        " << "beta  " <<   std::setprecision(0) << std::scientific << alpha  << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_beta         << std::endl;
      intgr_fstream << "The value of the H1 control for        " << "gamma " <<   std::setprecision(0) << std::scientific << beta << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_gamma        << std::endl;
      intgr_fstream << "The value of the integral of g.n "<<    integral_g_dot_n << std::endl;
      intgr_fstream << "The value of the theta is                             " <<    std::setw(11) << std::setprecision(10) <<  solTheta << std::endl;
      intgr_fstream << "The value of the total integral is " << std::setw(11) << std::setprecision(10) <<  integral_target_alpha * cost_functional_coeff * 0.5  + integral_beta * alpha * 0.5 + integral_gamma * beta * 0.5 << std::endl;
      intgr_fstream <<  std::endl;
      
      intgr_fstream.close();  //you have to close to disassociate the file from the stream
}  
     
    return; 
	  
  
}

  
  
  /** This function computes a functional with a volume part and a boundary part
     We pass a 2 Solution objects: the first for the cost functional, the second for the regularization
    */
void compute_cost_functional_regularization_bdry(const MultiLevelProblem & ml_prob, 
                     const unsigned level, 
                     const unsigned iteration,
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars  
                    )  {
    
  std::cout << "=== Compute cost functional parts with boundary regularization =============" << std::endl;  
  
  
  
  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //=============== Geometry ========================================
   unsigned solType_coords = FE_DOMAIN;
 
  CurrentElem < double > geom_element_iel(dim, msh);
    
  constexpr unsigned int space_dim = 3;
  
  std::vector<double> normal(space_dim, 0.);
 //***************************************************

  //=============== Integration ========================================

 //***************************************************
  const double alpha = ALPHA_CTRL_BDRY;
  const double beta  = BETA_CTRL_BDRY;
  
 //*************** state ***************************** 
 //***************************************************
  const unsigned n_components_state = 1;
  
  vector <double> phi_u;     phi_u.reserve(max_size);
  vector <double> phi_u_x;   phi_u_x.reserve(max_size * space_dim);
//   vector <double> phi_u_xx;  phi_u_xx.reserve(max_size * dim2);
 
  unsigned solIndex_u = ml_sol->GetIndex( state_vars[ n_components_state - 1].c_str() );
  unsigned solType_u  = ml_sol->GetSolutionType(solIndex_u);

  vector < double >  sol_u;
  sol_u.reserve(max_size);
  
  double u_gss = 0.;
  double u_x_gss = 0.;
 //*************************************************** 
 //***************************************************

  
 //************** cont *******************************
 //***************************************************
  const unsigned n_components_ctrl = 1;
  
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(max_size);
  phi_ctrl_x_bdry.reserve(max_size * space_dim);

  unsigned solIndex_ctrl = ml_sol->GetIndex(  ctrl_vars[ n_components_ctrl - 1].c_str() );
  unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

   vector < double >  sol_ctrl;   sol_ctrl.reserve(max_size);
 //***************************************************
 //*************************************************** 
  
  
 //************** desired ****************************
 //***************************************************
  vector <double> phi_udes;
  vector <double> phi_udes_x;

    phi_udes.reserve(max_size);
    phi_udes_x.reserve(max_size * space_dim);
 

  vector < double >  sol_udes;
  sol_udes.reserve(max_size);

  double udes_gss = 0.;
 //***************************************************
 //***************************************************

 //********** DATA *********************************** 
  double u_des = femus::ctrl::cost_functional::DesiredTarget();
 //*************************************************** 
  
  double integral_target = 0.;
  double integral_alpha  = 0.;
  double integral_beta   = 0.;

  

 //*************************************************** 
  constexpr unsigned qrule_i = QRULE_I;
  
     std::vector < std::vector < double > >  JacI_qp(space_dim);
     std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < Jac_qp.size(); d++) {   Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp.size(); d++) { JacI_qp[d].resize(dim); }
    
    double detJac_iqp;
  double weight_iqp = 0.;

     std::vector < std::vector < double > >  JacI_qp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_qp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_qp_bdry.size(); d++) {   Jac_qp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp_bdry.size(); d++) { JacI_qp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
  double weight_iqp_bdry = 0.;
    
      //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all;
  ml_prob.get_all_abstract_fe_multiple(elem_all);
 //*************************************************** 
  
  
  
  
  
  
  
  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element_iel.geom_type();

  //************* set target domain flag **************
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = femus::ctrl::cost_functional::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
 //***************************************************

   
 //*********** state ********************************* 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);
    sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
      unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
    }
 //*********** state ********************************* 


 //*********** cont ********************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);
    sol_ctrl    .resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);
    } 

 //*********** cont ********************************** 
 
 
 //*********** udes ********************************** 
    unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u);
    sol_udes    .resize(nDof_udes);
    for (unsigned i = 0; i < sol_udes.size(); i++) {
            sol_udes[i] = u_des;  //dof value
    } 
 //*********** udes ********************************** 

 
 //********** ALL VARS ******************************* 
    int nDof_max    =  nDof_u;   //  TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    if(nDof_udes > nDof_max) 
    {
      nDof_max = nDof_udes;
      }
    
    if(nDof_ctrl > nDof_max)
    {
      nDof_max = nDof_ctrl;
    }
    
 //***************************************************

//=================== BOUNDARY PART - BEGIN ==================================================================================================  
  
	if ( ctrl::Gamma_control::volume_elem_contains_a_Gamma_control_face( geom_element_iel.get_elem_center_3d() ) ) {
	  
	       
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       const unsigned nve_bdry_ctrl = msh->GetElementFaceDofNumber(iel,iface,solType_ctrl);
       
// ----------
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();
// ----------

		
	    if( ctrl::Gamma_control::face_is_a_Gamma_control_face(msh->el, iel, iface) ) {

	
		//============ initialize gauss quantities on the boundary ==========================================
                double sol_ctrl_bdry_gss = 0.;
                std::vector<double> sol_ctrl_x_bdry_gss(space_dim);
		//============ initialize gauss quantities on the boundary ==========================================
		
		for(unsigned ig_bdry = 0; ig_bdry < ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber(); ig_bdry++) {
		  
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_qp_bdry, JacI_qp_bdry, detJac_iqp_bdry, space_dim);
    weight_iqp_bdry = detJac_iqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];
    elem_all[qrule_i][ielGeom_bdry][solType_ctrl] ->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_ctrl_bdry, phi_ctrl_x_bdry, boost::none, space_dim);

		  
		 //========== compute gauss quantities on the boundary ===============================================
		  sol_ctrl_bdry_gss = 0.;
                  std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);
		      for (int i_bdry = 0; i_bdry < nve_bdry_ctrl; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_ctrl_bdry_gss +=  sol_ctrl[i_vol] * phi_ctrl_bdry[i_bdry];
                            for (int d = 0; d < space_dim; d++) {
			      sol_ctrl_x_bdry_gss[d] += sol_ctrl[i_vol] * phi_ctrl_x_bdry[i_bdry * space_dim + d];
			    }
		      }
		      
		      double laplace_ctrl_surface = 0.;  for (int d = 0; d < space_dim; d++) { laplace_ctrl_surface += sol_ctrl_x_bdry_gss[d] * sol_ctrl_x_bdry_gss[d]; }

                 //========= compute gauss quantities on the boundary ================================================
                  integral_alpha +=  weight_iqp_bdry * sol_ctrl_bdry_gss * sol_ctrl_bdry_gss; 
                  integral_beta  +=  weight_iqp_bdry * laplace_ctrl_surface;
                 
             }
	      } //end face
	      
	  }  // loop over element faces   
	  
	} //end if control element flag

//=================== BOUNDARY PART - END ==================================================================================================  

  
  
   
//=================== VOLUME PART - BEGIN ==================================================================================================  

   for (unsigned ig = 0; ig < ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom).GetGaussPointsNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), ig, Jac_qp, JacI_qp, detJac_iqp, space_dim);
    weight_iqp = detJac_iqp * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom).GetGaussWeightsPointer()[ig];

    elem_all[qrule_i][ielGeom][solType_u]                 ->shape_funcs_current_elem(ig, JacI_qp, phi_u, phi_u_x, boost::none, space_dim);
    elem_all[qrule_i][ielGeom][solType_u/*solTypeTdes*/]  ->shape_funcs_current_elem(ig, JacI_qp, phi_udes, phi_udes_x, boost::none, space_dim);
    
	u_gss     = 0.;  for (unsigned i = 0; i < nDof_u; i++)        u_gss += sol_u[i]     * phi_u[i];
	udes_gss  = 0.;  for (unsigned i = 0; i < nDof_udes; i++)  udes_gss += sol_udes[i]  * phi_udes[i];

    u_x_gss  = 0.;
        for (unsigned i = 0; i < nDof_u; i++)  {
          for (unsigned idim = 0; idim < dim; idim ++) u_x_gss  += sol_u[i] * phi_u_x[i * space_dim + idim];
        }

               integral_target +=  weight_iqp * target_flag*
#if COST_FUNCTIONAL_TYPE == 0
               (u_gss  - udes_gss) * (u_gss - udes_gss)
#elif COST_FUNCTIONAL_TYPE == 1
               u_x_gss * u_x_gss
#endif
               ;
	  
      } // end gauss point loop
//=================== VOLUME PART - END ==================================================================================================  
      
  } //end element loop

  ////////////////////////////////////////
  double total_integral = 0.5 * integral_target + 0.5 * alpha * integral_alpha + 0.5 * beta * integral_beta;
  
  std::cout << "total integral on processor " << iproc << ": " << total_integral << std::endl;

  double integral_target_parallel = 0.; MPI_Allreduce( &integral_target, &integral_target_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_target_parallel << std::endl;
  double integral_alpha_parallel = 0.; MPI_Allreduce( &integral_alpha, &integral_alpha_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_alpha_parallel << std::endl;
  double integral_beta_parallel = 0.;  MPI_Allreduce( &integral_beta, &integral_beta_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_beta_parallel << std::endl;

  double total_integral_parallel = 0.; MPI_Allreduce( &total_integral, &total_integral_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << total_integral_parallel << std::endl;


  
 
return;
  
}
  

  
  

void compute_cost_functional_regularization_lifting_internal(const MultiLevelProblem & ml_prob, 
                     const unsigned level, 
                     const unsigned iteration,
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars  
                    )   {
  
  std::cout << "=== Compute cost functional parts with internal lifting regularization =============" << std::endl;  
  
  Mesh*                          msh = ml_prob._ml_msh->GetLevel(level);

  MultiLevelSolution*          ml_sol = ml_prob._ml_sol;
  Solution*                      sol = ml_prob._ml_sol->GetSolutionLevel(level);

  const unsigned     dim = msh->GetDimension();                                 // get the domain dimension of the problem
  const unsigned    dim2 = (3 * (dim - 1) + !(dim - 1));                        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  const unsigned   iproc = msh->processor_id(); 

 //********** Geometry ***************************************** 
 unsigned solType_coords = FE_DOMAIN;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
 //*************************************************** 

  
 //***************************************************  
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;

 //******************** state ************************ 
 //*************************************************** 
  vector <double> phi_u;
  vector <double> phi_u_x;
  vector <double> phi_u_xx;

  phi_u.reserve(max_size);
  phi_u_x.reserve(max_size * dim);
  phi_u_xx.reserve(max_size * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = ml_sol->GetIndex("state");
  unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);

  vector < double >  sol_u; // local solution
  sol_u.reserve(max_size);
  
  double u_gss = 0.;
double u_x_gss = 0.;
 //*************************************************** 

 //******************** control ********************** 
 //*************************************************** 
  vector <double> phi_ctrl;
  vector <double> phi_ctrl_x;
  vector <double> phi_ctrl_xx;

  phi_ctrl.reserve(max_size);
  phi_ctrl_x.reserve(max_size * dim);
  phi_ctrl_xx.reserve(max_size * dim2);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = ml_sol->GetIndex("control");
  unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

  vector < double >  sol_ctrl; // local solution
  sol_ctrl.reserve(max_size);
  
  double ctrl_gss = 0.;
  double ctrl_x_gss = 0.;
 //*************************************************** 
 //***************************************************  

  
 //********************* desired ********************* 
 //*************************************************** 
  vector <double> phi_udes;
  vector <double> phi_udes_x;
  vector <double> phi_udes_xx;

  phi_udes.reserve(max_size);
  phi_udes_x.reserve(max_size * dim);
  phi_udes_xx.reserve(max_size * dim2);
 
  
//  unsigned solIndex_udes;
//   solIndex_udes = ml_sol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solType_udes = ml_sol->GetSolutionType(solIndex_udes);    // get the finite element type for "state"

  vector < double >  sol_udes; // local solution
  sol_udes.reserve(max_size);

  double udes_gss = 0.;
 //*************************************************** 
 //*************************************************** 

 //********************* DATA ************************ 
  double u_des = femus::ctrl::cost_functional::DesiredTarget();
 //*************************************************** 
  
  double integral_target = 0.;
  double integral_alpha  = 0.;
  double integral_beta   = 0.;


 //***************************************************  
       //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all;
  ml_prob.get_all_abstract_fe_multiple(elem_all);
  
   constexpr unsigned qrule_i = QRULE_I;
 //***************************************************  


   // ====== Geometry at Quadrature points - BEGIN ==============================================================================
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
    double detJac_iqp;
    double AbsDetJxWeight_iqp = 0.;
   // ====== Geometry at Quadrature points - END ==============================================================================

    
  
    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
 
    geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element_iel.geom_type();

  //************* set target domain flag **************
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = femus::ctrl::cost_functional::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
 //*************************************************** 

 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ctrl::Gamma_control::ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
 //*************************************************** 
   
 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);
        sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);
            sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
    }
 //*************************************************** 


 //************** control **************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);
    sol_ctrl    .resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);
    } 
 //***************************************************  
 
 
 //**************** u_des **************************** 
    unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u);
    sol_udes    .resize(nDof_udes);
    for (unsigned i = 0; i < sol_udes.size(); i++) {
      sol_udes[i] = u_des;  //dof value
    } 
 //*************************************************** 

 
 //******************* ALL VARS ********************** 
    int nDof_max    =  nDof_u;   // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    if(nDof_udes > nDof_max) 
    {
      nDof_max = nDof_udes;
      }
    
    if(nDof_ctrl > nDof_max)
    {
      nDof_max = nDof_ctrl;
    }
    
 //*************************************************** 
   
      // *** Gauss point loop ***
      for (unsigned iqp = 0; iqp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); iqp++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), iqp, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);

    AbsDetJxWeight_iqp = detJac_iqp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[iqp];

    
    elem_all[qrule_i][ielGeom][solType_u]  ->shape_funcs_current_elem(iqp, JacI_iqp, phi_u, phi_u_x, boost::none, space_dim);
    elem_all[qrule_i][ielGeom][solType_ctrl]  ->shape_funcs_current_elem(iqp, JacI_iqp, phi_ctrl, phi_ctrl_x, boost::none, space_dim);
    elem_all[qrule_i][ielGeom][solType_u]  ->shape_funcs_current_elem(iqp, JacI_iqp, phi_udes, phi_udes_x, boost::none, space_dim);
    
//         // *** get gauss point weight, test function and test function partial derivatives ***
// 	    msh->_finiteElement[ielGeom][solType_u]               ->Jacobian(geom_element.get_coords_at_dofs(), iqp, weight, phi_u, phi_u_x, phi_u_xx);
//         msh->_finiteElement[ielGeom][solType_ctrl]            ->Jacobian(geom_element.get_coords_at_dofs(), iqp, weight, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);
//         msh->_finiteElement[ielGeom][solType_u/*solTypeudes*/]->Jacobian(geom_element.get_coords_at_dofs(), iqp, weight, phi_udes, phi_udes_x, phi_udes_xx);

	u_gss = 0.;  for (unsigned i = 0; i < nDof_u; i++) u_gss += sol_u[i] * phi_u[i];		
	ctrl_gss = 0.; for (unsigned i = 0; i < nDof_ctrl; i++) ctrl_gss += sol_ctrl[i] * phi_ctrl[i];  
	udes_gss  = 0.; for (unsigned i = 0; i < nDof_udes; i++)  udes_gss  += sol_udes[i]  * phi_udes[i]; 
        ctrl_x_gss  = 0.; 
        for (unsigned i = 0; i < nDof_ctrl; i++)  {
          for (unsigned idim = 0; idim < dim; idim ++) ctrl_x_gss  += sol_ctrl[i] * phi_ctrl_x[i * space_dim + idim];
        }

        u_x_gss  = 0.; 
        for (unsigned i = 0; i < nDof_u; i++)  {
          for (unsigned idim = 0; idim < dim; idim ++) u_x_gss  += sol_u[i] * phi_u_x[i * space_dim + idim];
        }
               integral_target += AbsDetJxWeight_iqp * target_flag * 
#if COST_FUNCTIONAL_TYPE == 0               
               (u_gss +  ctrl_gss - udes_gss) * (u_gss +  ctrl_gss - udes_gss)
#elif COST_FUNCTIONAL_TYPE == 1
               (u_x_gss + ctrl_x_gss) * (u_x_gss + ctrl_x_gss) 
#endif
               ;
               integral_alpha  += AbsDetJxWeight_iqp * control_el_flag * ctrl_gss * ctrl_gss;
               integral_beta   += AbsDetJxWeight_iqp * control_el_flag * ctrl_x_gss * ctrl_x_gss;
	  
      } // end gauss point loop
  } //end element loop

  ////////////////////////////////////////
  double total_integral = 0.5 * integral_target + 0.5 * alpha * integral_alpha + 0.5 * beta * integral_beta;
  
  std::cout << "total integral on processor " << iproc << ": " << total_integral << std::endl;

  double integral_target_parallel = 0.; MPI_Allreduce( &integral_target, &integral_target_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  double integral_alpha_parallel = 0.; MPI_Allreduce( &integral_alpha, &integral_alpha_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  double integral_beta_parallel = 0.;  MPI_Allreduce( &integral_beta, &integral_beta_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  double total_integral_parallel = 0.; MPI_Allreduce( &total_integral, &total_integral_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  
  std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_target_parallel << std::endl;
  std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_alpha_parallel << std::endl;
  std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_beta_parallel << std::endl;
  std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << total_integral_parallel << std::endl;

return;
  
}
  
  

void compute_cost_functional_regularization_lifting_internal_vec(
                     const MultiLevelProblem & ml_prob, 
                     const unsigned level,
                     const unsigned iteration,
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars) {

    
    
  const double cost_functional_coeff = COST_FUNCTIONAL_COEFF;
  const double alpha = ALPHA_CTRL_VOL;
  const double beta  = BETA_CTRL_VOL;

  
  Mesh*          msh          	= ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         	= msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  ml_sol    = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        	= ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object
  
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  unsigned solType_coords = coordXType;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = /*dim*/  3  /*2*/    ;
 
  vector < vector < double > > coordX(dim);    // local coordinates

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(max_size);
  }
  
  double AbsDetJxWeight_iqp;
  
  
  //geometry *******************************

//STATE######################################################################
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = ml_sol->GetIndex("u_0");
  solVIndex[1] = ml_sol->GetIndex("u_1");

  if (dim == 3) solVIndex[2] = ml_sol->GetIndex("u_2");

  unsigned solVType = ml_sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
  
  vector < vector < double > >  solV(dim);    // local solution
  vector <double >  V_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(max_size);
  }

  
  vector <double> phiV_gss;  // local test function
  vector <double> phiV_x_gss; // local test function first order partial derivatives

  phiV_gss.reserve(max_size);
  phiV_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  
//STATE######################################################################
  

//CONTROL######################################################################
  vector < unsigned > solVctrlIndex(dim);
  solVctrlIndex[0] = ml_sol->GetIndex("ctrl_0");
  solVctrlIndex[1] = ml_sol->GetIndex("ctrl_1");

  if (dim == 3) solVctrlIndex[2] = ml_sol->GetIndex("ctrl_2");

  unsigned solVctrlType = ml_sol->GetSolutionType(solVctrlIndex[0]);
  
  vector < vector < double > >  solVctrl(dim);    // local solution
  vector < double >   Vctrl_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solVctrl[k].reserve(max_size);
  }

  
  vector <double> phiVctrl_gss;  // local test function
  vector <double> phiVctrl_x_gss; // local test function first order partial derivatives

  phiVctrl_gss.reserve(max_size);
  phiVctrl_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  
//CONTROL######################################################################

// Vel_desired##################################################################
  vector <double> phiVdes_gss;  // local test function
  vector <double> phiVdes_x_gss; // local test function first order partial derivatives

  phiVdes_gss.reserve(max_size);
  phiVdes_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);

  vector <double>  solVdes(dim,0.);
  vector<double> Vdes_gss(dim, 0.);  
  
// Vel_desired##################################################################




double  integral_target_alpha = 0.;
double	integral_beta   = 0.;
double	integral_gamma  = 0.;
double  integral_div_ctrl = 0.;



//*************************************************** 
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
    double detJac_iqp;

    //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
//*************************************************** 

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

   // geometry *****************************
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element_iel.geom_type();
  // geometry end *****************************
   
// equation
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);
//     unsigned nDofsVdes = msh->GetElementDofNumber(iel, solVType); 
    unsigned nDofsVctrl = msh->GetElementDofNumber(iel, solVctrlType);
    
     
    for (unsigned  k = 0; k < dim; k++)  {
      solV[k].resize(nDofsV);
      solVctrl[k].resize(nDofsVctrl);
//       solVdes[k].resize(nDofsVdes);
    }
  //*************************************** 
  
  //***** set target domain flag ********************************** 
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = femus::ctrl::cost_functional::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
//***************************************       
    
 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ctrl::Gamma_control::ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
 //*************************************************** 
   
    
 //STATE###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
      }
    }
//STATE###################################################################

//CONTROL###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVctrlDof = msh->GetSolutionDof(i, iel, solVctrlType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solVctrl[k][i] = (*sol->_Sol[solVctrlIndex[k]])(solVctrlDof);      // global extraction and local storage for the solution
      }
    }
//CONTROL###################################################################




  //DESIRED VEL###################################################################  
    // velocity ************
//     for (unsigned i = 0; i < nDofsV; i++) {
//       unsigned solVdesDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < solVdes.size() /*dim*/; k++) {
        solVdes[k]/*[i]*/ = femus::ctrl::cost_functional::DesiredTargetVec()[k] /*(*sol->_Sol[solVIndex[k]])(solVdesDof)*/;      // global extraction and local storage for the solution
      }
//     }
 //DESIRED VEL###################################################################

 
 
 
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {

//STATE#############################################################################	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), ig, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);
    AbsDetJxWeight_iqp = detJac_iqp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[ig];
   
    elem_all[ielGeom][solVType]->shape_funcs_current_elem(ig, JacI_iqp, phiV_gss, phiV_x_gss, boost::none , space_dim);
    elem_all[ielGeom][solVctrlType]->shape_funcs_current_elem(ig, JacI_iqp, phiVctrl_gss, phiVctrl_x_gss,  boost::none, space_dim);
    elem_all[ielGeom][solVType /*solVdes*/]->shape_funcs_current_elem(ig, JacI_iqp, phiVdes_gss, phiVdes_x_gss,  boost::none, space_dim);

	
	  vector < vector < double > > gradVctrl_gss(dim);
      for (unsigned  k = 0; k < dim; k++) {
          gradVctrl_gss[k].resize(dim_offset_grad /*space_dim*/);
          std::fill(gradVctrl_gss[k].begin(), gradVctrl_gss[k].end(), 0);
        }
	
    for (unsigned  k = 0; k < dim; k++) {
      V_gss[k]       = 0.;
      Vdes_gss[k]    = 0.;
       Vctrl_gss[k]  = 0.;
    }
    
      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned  k = 0; k < dim; k++) {
            V_gss[k] += solV[k][i] * phiV_gss[i];
            Vdes_gss[k] += solVdes[k]/*[i]*/ * phiVdes_gss[i];
		}
      }
	
      for (unsigned i = 0; i < nDofsVctrl; i++) {
        for (unsigned  k = 0; k < dim; k++) {
            Vctrl_gss[k] += solVctrl[k][i] * phiVctrl_gss[i];
	 }
     for (unsigned j = 0; j < dim_offset_grad /*space_dim*/; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              gradVctrl_gss[k][j] += phiVctrl_x_gss[i * dim_offset_grad /*space_dim*/ + j] * solVctrl[k][i];
            }
          }
      }
          
//                 for (unsigned i = 0; i < nDofsV; i++) {

      for (unsigned  k = 0; k < dim; k++) {
          integral_div_ctrl +=  AbsDetJxWeight_iqp * gradVctrl_gss[k][k] /** phiVctrl_gss[i]*/;
      }
//       }
	
      for (unsigned  k = 0; k < dim; k++) {
	 integral_target_alpha +=  target_flag * (V_gss[k] + Vctrl_gss[k] - Vdes_gss[k]) * (V_gss[k] + Vctrl_gss[k] - Vdes_gss[k]) * AbsDetJxWeight_iqp; 
	 integral_beta	+=  control_el_flag * ((Vctrl_gss[k]) * (Vctrl_gss[k]) * AbsDetJxWeight_iqp);
      }
      for (unsigned  k = 0; k < dim; k++) {
	for (unsigned  j = 0; j < dim; j++) {	
		integral_gamma	  +=  control_el_flag * gradVctrl_gss[k][j] * gradVctrl_gss[k][j] * AbsDetJxWeight_iqp;
	}
      }
   
  
      }// end gauss point loop
    } //end element loop  

       std::ostringstream filename_out; filename_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "Integral_computation"  << ".txt";

       std::ofstream intgr_fstream;
  if (paral::get_rank() == 0 ) {
      intgr_fstream.open(filename_out.str().c_str(),std::ios_base::app);
      intgr_fstream << " ***************************** Iteration " << iteration << " *********************************** " <<  std::endl << std::endl;
      intgr_fstream << "The value of the target functional for " << "alpha " <<   std::setprecision(0) << std::scientific << cost_functional_coeff << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_target_alpha << std::endl;
      intgr_fstream << "The value of the L2 control for        " << "beta  " <<   std::setprecision(0) << std::scientific << alpha  << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_beta         << std::endl;
      intgr_fstream << "The value of the H1 control for        " << "gamma " <<   std::setprecision(0) << std::scientific << beta << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_gamma        << std::endl;
      intgr_fstream << "The value of the total integral is " << std::setw(11) << std::setprecision(10) <<  integral_target_alpha * cost_functional_coeff*0.5  + integral_beta *alpha*0.5 + integral_gamma *beta*0.5 << std::endl;
      intgr_fstream << "The value of the divergence of the control is " << std::setw(11) << std::setprecision(10) <<  integral_div_ctrl << std::endl;
      intgr_fstream <<  std::endl;
      intgr_fstream.close();  //you have to close to disassociate the file from the stream
}  

    
//     std::cout << "The value of the integral of target for alpha "<< std::setprecision(0)<< std::scientific<<  cost_functional_coeff<< " is " << std::setw(11) << std::setprecision(10) << std::fixed<< integral_target_alpha << std::endl;
//     std::cout << "The value of the integral of beta for beta "<<  std::setprecision(0)<<std::scientific<<alpha << " is " << std::setw(11) << std::setprecision(10) <<  std::fixed<< integral_beta << std::endl;
//     std::cout << "The value of the integral of gamma for gamma "<< std::setprecision(0)<<std::scientific<<beta<< " is " << std::setw(11) << std::setprecision(10) <<  std::fixed<< integral_gamma << std::endl; 
//     std::cout << "The value of the total integral is " << std::setw(11) << std::setprecision(10) <<  integral_target_alpha *(cost_functional_coeff*0.5)+ integral_beta *(alpha*0.5) + integral_gamma*(beta*0.5) << std::endl; 
   
    
    return; 
	  
  
}

  
 }
 
 


namespace ctrl_inequality {


 ///@@@@@@@@@@@@@@@@@@@todo I believe you need to pass other rows here.........
    // Also, I think you the active flag should be frozen if it was frozen once...!!!
    //I don't think I have to put other mu's for the other constraint equations... there is a mu per constrained variable... 
    //just make sure where they go
    //only try to constrain each component, and see how it behaves.
    //So, the active set is correctly found, but mu does not change, and it must! Maybe I didn't put all the pieces from the scalar to the vector routine?!
    //ctrl is not modified correctly...
    //I think the constraint is only acting on the interior nodes but not on the boundary, just check the boundary conditions!!!
 std::vector<double>  InequalityConstraint(const unsigned n_components_ctrl, const std::vector<double> & dof_obj_coord, const bool upper) {

     const unsigned dim = dof_obj_coord.size();
     
     std::vector<double> constr_value(n_components_ctrl, 0.);
     
     
     double constr_value_upper_0 =  1000000.; // dof_obj_coord[1]*(1. - dof_obj_coord[1]);
     double constr_value_lower_0 = -1000000.; //-3.e-13;
     assert(constr_value_lower_0 < constr_value_upper_0); 
     
     double constr_value_upper_1 =   1000.;
     double constr_value_lower_1 =  -1000.;
     assert(constr_value_lower_1 < constr_value_upper_1); 
     
     double constr_value_upper_2 =  1000.;
     double constr_value_lower_2 = -1000.;
     assert(constr_value_lower_2 < constr_value_upper_2); 
     
    if (upper)   { 
                      constr_value[0] = constr_value_upper_0; 
                      constr_value[1] = constr_value_upper_1; 
       if (dim == 3)  constr_value[2] = constr_value_upper_2;
    }
    else         { 
                      constr_value[0] = constr_value_lower_0; 
                      constr_value[1] = constr_value_lower_1; 
       if (dim == 3)  constr_value[2] = constr_value_lower_2; 
    }
    
    
  return constr_value;
     
}
   
   
    
    
    
 void update_active_set_flag_for_current_nonlinear_iteration(const femus::Mesh* msh,
                                                             const femus::Solution* sol,
                                                             const unsigned int iel,
                                                             const std::vector < std::vector < double > > & coords_at_dofs,
                                                             const std::vector < std::vector < double > > sol_eldofs,  ///@todo why not a reference?
                                                             const std::vector < unsigned int > & Sol_n_el_dofs,
                                                             const double c_compl,
                                                             const std::vector<unsigned int>  pos_mu,
                                                             const std::vector<unsigned int>  pos_ctrl,
                                                             const std::vector<unsigned int>  solIndex_act_flag,
                                                                 std::vector <   std::vector < double > > & ctrl_lower_dofs,
                                                                 std::vector <   std::vector < double > > & ctrl_upper_dofs,
                                                                 std::vector <   std::vector < double > > & sol_actflag_dofs) {
     
        const unsigned int dim = coords_at_dofs.size();
        
    const unsigned int   n_components_ctrl = pos_ctrl.size();
    
    
              for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 

        
// 0: inactive; 1: active_a; 2: active_b
        assert(Sol_n_el_dofs[pos_mu[kdim] ] == Sol_n_el_dofs[ pos_ctrl[kdim] ]);
        sol_actflag_dofs[kdim].resize(Sol_n_el_dofs[pos_mu[kdim] ]);
        ctrl_lower_dofs[kdim].resize(Sol_n_el_dofs[pos_mu[kdim] ]);
        ctrl_upper_dofs[kdim].resize(Sol_n_el_dofs[pos_mu[kdim] ]);
        std::fill(sol_actflag_dofs[kdim].begin(), sol_actflag_dofs[kdim].end(), 0);
        std::fill(ctrl_lower_dofs[kdim].begin(),   ctrl_lower_dofs[kdim].end(), 0.);
        std::fill(ctrl_upper_dofs[kdim].begin(),   ctrl_upper_dofs[kdim].end(), 0.);

// // //         std::cout << " mu dofs " << std::endl;
// // //                 for (unsigned i = 0; i < sol_actflag_dofs.size(); i++) {
// // //                 std::cout << sol_eldofs[pos_mu][i] << " ";
// // //         }
// // //         
// // //         std::cout << std::endl;
        
// // //         std::cout << " ctrl dofs " << std::endl;
// // //         for (unsigned i = 0; i < sol_actflag_dofs.size(); i++) {
// // //                 std::cout << sol_eldofs[pos_ctrl][i] << " ";
// // //         }
// // //         std::cout << std::endl;
        
        for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) {
            std::vector<double> node_coords_i(dim, 0.);
            for (unsigned d = 0; d < dim; d++) node_coords_i[d] = coords_at_dofs[d][i];
            
            ctrl_lower_dofs[kdim][i] = ctrl::ctrl_inequality::InequalityConstraint(n_components_ctrl, node_coords_i, false)[kdim];
            ctrl_upper_dofs[kdim][i] = ctrl::ctrl_inequality::InequalityConstraint(n_components_ctrl, node_coords_i, true)[kdim];
            
            const double lower_test_value = sol_eldofs[ pos_mu[kdim] ][i] + c_compl * ( sol_eldofs[ pos_ctrl[kdim] ][i] - ctrl_lower_dofs[kdim][i] );
            const double upper_test_value = sol_eldofs[ pos_mu[kdim] ][i] + c_compl * ( sol_eldofs[ pos_ctrl[kdim] ][i] - ctrl_upper_dofs[kdim][i] );

            if      ( lower_test_value < 0. )  {
//                 std::cout << "Found active node below" << std::endl;
//                 std::cout << "The current value of mu is " <<  sol_eldofs[ pos_mu[kdim] ][i] << std::endl;
                   sol_actflag_dofs[kdim][i] = 1;
            }
            else if ( upper_test_value > 0. )  {
//                 std::cout << "Found active node above" << std::endl;
//                 std::cout << "The current value of mu is " <<  sol_eldofs[ pos_mu[kdim] ][i] << std::endl;
                sol_actflag_dofs[kdim][i] = 2;
            }
        }

//************** local to global act flag ***************************
     const unsigned solFEType_act_flag = sol->GetSolutionType(solIndex_act_flag[kdim]);

        for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) {
            unsigned solDof_mu = msh->GetSolutionDof(i, iel, solFEType_act_flag);
            (sol->_Sol[ solIndex_act_flag[kdim] ])->set(solDof_mu, sol_actflag_dofs[kdim][i]);
        }


      }
      
                                                                 
   }





 void update_active_set_flag_for_current_nonlinear_iteration_bdry(const femus::Mesh* msh,
                                                             const femus::Solution* sol,
                                                             const unsigned int iel,
                                                             const unsigned int iface,
                                                             const std::vector < std::vector < double > > & coords_at_dofs,
                                                             const std::vector < std::vector < double > > sol_eldofs,  ///@todo why not a reference?
                                                             const std::vector < unsigned int > & Sol_n_el_dofs,
                                                             const double c_compl,
                                                             const std::vector<unsigned int> pos_mu,
                                                             const std::vector<unsigned int> pos_ctrl,
                                                             const std::vector<unsigned int>  solIndex_act_flag,
                                                              std::vector < std::vector < double > > & ctrl_lower_dofs, 
                                                              std::vector < std::vector < double > > & ctrl_upper_dofs,
                                                              std::vector < std::vector < double > > & sol_actflag_dofs) {
     
         const unsigned dim = coords_at_dofs.size();
         
      const unsigned int   n_components_ctrl = pos_ctrl.size();
    
       for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
     

     
     const unsigned solFEType_act_flag = sol->GetSolutionType(solIndex_act_flag[kdim]);
     
     const unsigned ndofs_bdry = msh->GetElementFaceDofNumber(iel, iface, solFEType_act_flag);
        
        
         // 0: inactive; 1: active_a; 2: active_b
        assert(Sol_n_el_dofs[ pos_mu[kdim] ] == Sol_n_el_dofs[ pos_ctrl[kdim] ]);///@todo More appropriately, 
            sol_actflag_dofs[kdim].resize(ndofs_bdry/*nDof_mu*/);
            ctrl_lower_dofs[kdim].resize(ndofs_bdry/*nDof_mu*/);
            ctrl_upper_dofs[kdim].resize(ndofs_bdry/*nDof_mu*/);
           std::fill(sol_actflag_dofs[kdim].begin(), sol_actflag_dofs[kdim].end(), 0);
           std::fill(ctrl_lower_dofs[kdim].begin(), ctrl_lower_dofs[kdim].end(), 0.);
           std::fill(ctrl_upper_dofs[kdim].begin(), ctrl_upper_dofs[kdim].end(), 0.);

           
      for (unsigned int i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
        std::vector<double> node_coords_i(dim, 0.);
        for (unsigned d = 0; d < dim; d++) node_coords_i[d] = coords_at_dofs[d][i_bdry];
        
        ctrl_lower_dofs[kdim][i_bdry] = ctrl::ctrl_inequality::InequalityConstraint(n_components_ctrl, node_coords_i, false)[kdim];
        ctrl_upper_dofs[kdim][i_bdry] = ctrl::ctrl_inequality::InequalityConstraint(n_components_ctrl, node_coords_i, true)[kdim];

        const double lower_test_value = sol_eldofs[ pos_mu[kdim] ][i_vol] + c_compl * ( sol_eldofs[ pos_ctrl[kdim] ][i_vol] - ctrl_lower_dofs[kdim][i_bdry] );
        const double upper_test_value = sol_eldofs[ pos_mu[kdim] ][i_vol] + c_compl * ( sol_eldofs[ pos_ctrl[kdim] ][i_vol] - ctrl_upper_dofs[kdim][i_bdry] );
        
        if      ( lower_test_value < 0. )  sol_actflag_dofs[kdim][i_bdry] = 1;
        else if ( upper_test_value > 0. )  sol_actflag_dofs[kdim][i_bdry] = 2;
            }
            
//************** local to global act flag ***************************
      for (unsigned int i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
      unsigned solDof_actflag = msh->GetSolutionDof(i_vol, iel, solFEType_act_flag); 
      (sol->_Sol[solIndex_act_flag[kdim]])->set(solDof_actflag, sol_actflag_dofs[kdim][i_bdry]);     
           }
    
       }
       
       
    }
    
 
 
 void node_insertion_bdry(const unsigned int iel,
                         const unsigned int iface,
                         const   Mesh* msh,
                         const     vector < vector < int > > & L2G_dofmap,
                                                             const std::vector<unsigned int> pos_mu,
                                                             const std::vector<unsigned int> pos_ctrl,
                         const std::vector < std::vector < double > > & sol_eldofs,
                         const std::vector < unsigned int > & Sol_n_el_dofs,
                         std::vector < std::vector < double > > & sol_actflag_dofs,
                         const std::vector < std::vector < double > > & ctrl_lower_dofs,
                         const std::vector < std::vector < double > > & ctrl_upper_dofs,
                         const double ineq_flag,
                         const double c_compl,
                         SparseMatrix*             KK,
                         NumericVector* RES,
                          const bool assembleMatrix
                        ) {

     
   const unsigned int   n_components_ctrl = pos_ctrl.size();
    
       for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
     
     
     
// Create the L2G boundary maps from the volume ones
  std::vector < int > L2G_dofmap_mu_bdry(sol_actflag_dofs[kdim].size());
  std::vector < int > L2G_dofmap_ctrl_bdry(sol_actflag_dofs[kdim].size());

      for (int i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
  L2G_dofmap_mu_bdry[i_bdry]   = L2G_dofmap[ pos_mu[kdim]  ][i_vol];
  L2G_dofmap_ctrl_bdry[i_bdry] = L2G_dofmap[ pos_ctrl[kdim] ][i_vol];
      }
      
 //============= delta_mu row - BEGIN ===============================
      std::vector<double> Res_mu_bdry (sol_actflag_dofs[kdim].size());     std::fill(Res_mu_bdry.begin(),Res_mu_bdry.end(), 0.);
//       std::vector<double> Res_mu (Sol_n_el_dofs[pos_mu]);       std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
      for (int i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
        
      if (sol_actflag_dofs[kdim][i_bdry] == 0) {  //inactive
//          Res_mu [i_vol]      = - ineq_flag * ( 1. * sol_eldofs[pos_mu[kdim]][i_vol] - 0. ); 
         Res_mu_bdry[i_bdry] = - ineq_flag * ( 1. * sol_eldofs[pos_mu[kdim]][i_vol] - 0. ); 
      }
      else if (sol_actflag_dofs[kdim][i_bdry] == 1) {  //active_a 
// 	 Res_mu [i_vol]      = - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl[kdim]][i_vol] - c_compl * ctrl_lower_dofs[kdim][i_bdry]);
     Res_mu_bdry[i_bdry] = - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl[kdim]][i_vol] - c_compl * ctrl_lower_dofs[kdim][i_bdry]);
          
    }
      else if (sol_actflag_dofs[kdim][i_bdry] == 2) {  //active_b 
// 	Res_mu [i_vol]      =  - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl[kdim]][i_vol] - c_compl * ctrl_upper_dofs[kdim][i_bdry]);
    Res_mu_bdry[i_bdry] =  - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl[kdim]][i_vol] - c_compl * ctrl_upper_dofs[kdim][i_bdry]);
      }
    }

    
//     RES->insert(Res_mu,  L2G_dofmap[pos_mu]);    
    RES->insert(Res_mu_bdry,  L2G_dofmap_mu_bdry);    
 //============= delta_mu row - END ===============================
    
 //============= delta_mu-delta_ctrl row  - BEGIN ===============================
 //auxiliary volume vector for act flag
//  unsigned nDof_actflag_vol  = msh->GetElementDofNumber(iel, solFEType_act_flag);
//  std::vector<double> sol_actflag_dofs_vol(nDof_actflag_vol); 


 for (unsigned i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++) if (sol_actflag_dofs[kdim][i_bdry] != 0 ) sol_actflag_dofs[kdim][i_bdry] = ineq_flag * c_compl;    
 
//  std::fill(sol_actflag_dofs_vol.begin(), sol_actflag_dofs_vol.end(), 0.);
//     for (int i_bdry = 0; i_bdry < sol_actflag_dofs.size(); i_bdry++)  {
//        unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
//        sol_actflag_dofs_vol[i_vol] = sol_actflag_dofs[i_bdry];
//     }
 
//  KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_ctrl], sol_actflag_dofs_vol);
 if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap_mu_bdry, L2G_dofmap_ctrl_bdry, sol_actflag_dofs[kdim]); }
 //============= delta_mu-delta_ctrl row - END ===============================

 //============= delta_mu-delta_mu row - BEGIN ===============================
 // Attention: this equation goes in contrast with \mu = 0 on \Omega \setminus \Gamma_c
 // In fact, here we shouldn't insert all VOLUME values, but only the BOUNDARY ones
 // The best way is to then do a L2G_map ON THE BOUNDARY only
 
  for (unsigned i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++) sol_actflag_dofs[kdim][i_bdry] =  ineq_flag * (1 - sol_actflag_dofs[kdim][i_bdry]/c_compl)  + (1-ineq_flag) * 1.;  //can do better to avoid division, maybe use modulo operator 

//  std::fill(sol_actflag_dofs_vol.begin(), sol_actflag_dofs_vol.end(), 0.);
//     for (int i_bdry = 0; i_bdry < sol_actflag_dofs.size(); i_bdry++)  {
//        unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
//        sol_actflag_dofs_vol[i_vol] = sol_actflag_dofs[i_bdry];
//     }
  
//   KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_mu], sol_actflag_dofs_vol );
  if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap_mu_bdry, L2G_dofmap_mu_bdry, sol_actflag_dofs[kdim]);  }
 //============= delta_mu-delta_mu row - END ===============================
  
       }//kdim
       
       
}




 void node_insertion(const unsigned int iel,
                         const   Mesh* msh,
                         const     vector < vector < int > > & L2G_dofmap_Mat,
                                                             const std::vector<unsigned int> pos_mat_mu,
                                                             const std::vector<unsigned int> pos_mat_ctrl,
                         const std::vector < std::vector < double > > & sol_eldofs_Mat,
                         const std::vector < unsigned int > & Sol_n_el_dofs,
                         std::vector < std::vector < double > > & sol_actflag_dofs,
                         const std::vector < std::vector < double > > & ctrl_lower_dofs,
                         const std::vector < std::vector < double > > & ctrl_upper_dofs,
                         const double ineq_flag,
                         const double c_compl,
                         SparseMatrix*             KK,
                         NumericVector* RES,
                          const bool assembleMatrix
                        ) {
     
      
   const unsigned int   n_components_ctrl = pos_mat_ctrl.size();
    
       for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
     
    
     
      //*************************************************** 
    std::vector < int > l2GMap_mu(sol_actflag_dofs[kdim].size());
    std::vector < int > l2GMap_ctrl(sol_actflag_dofs[kdim].size());
    for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) {
      l2GMap_mu[i]   = L2G_dofmap_Mat[ pos_mat_mu[kdim]  ][i];   //pdeSys->GetSystemDof(solIndex_mu, solPdeIndex_mu, i, iel);
      l2GMap_ctrl[i] = L2G_dofmap_Mat[ pos_mat_ctrl[kdim] ][i]; //pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);
    } 
 //*************************************************** 

 //============= delta_mu row - BEGIN ===============================
      std::vector<double> Res_mu (sol_actflag_dofs[kdim].size()); std::fill(Res_mu.begin(),Res_mu.end(), 0.);
      
    for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) {
      if (sol_actflag_dofs[kdim][i] == 0){  //inactive
         Res_mu [i] = - ineq_flag * ( 1. * sol_eldofs_Mat[pos_mat_mu[kdim]][i] - 0. ); 
// 	 Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i]; 
      }
      else if (sol_actflag_dofs[kdim][i] == 1){  //active_a 
	 Res_mu [i] = - ineq_flag * ( c_compl *  sol_eldofs_Mat[pos_mat_ctrl[kdim]][i] - c_compl * ctrl_lower_dofs[kdim][i]);
      }
      else if (sol_actflag_dofs[kdim][i] == 2){  //active_b 
	Res_mu [i]  =  - ineq_flag * ( c_compl *  sol_eldofs_Mat[pos_mat_ctrl[kdim]][i] - c_compl * ctrl_upper_dofs[kdim][i]);
      }
    }
//          Res[nDof_u + nDof_ctrl + nDof_adj + i]  = c_compl * (  (2 - sol_actflag_dofs[i]) * (ctrl_lower_dofs[i] - sol_ctrl[i]) + ( sol_actflag_dofs[i] - 1 ) * (ctrl_upper_dofs[i] - sol_ctrl[i])  ) ;
//          Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i] ;

    
    RES->insert(Res_mu, l2GMap_mu);
//     RES->insert(Res_ctrl, l2GMap_ctrl);
//     RES->insert(Res_u, l2GMap_u);
//     RES->insert(Res_adj, l2GMap_adj);
 //============= delta_mu row - END ===============================
    
//  //============= delta_state-delta_state row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_u, l2GMap_u, 1.);

//  //============= delta_ctrl-delta_ctrl row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_ctrl, l2GMap_ctrl, 1.);
 
//  //============= delta_adj-delta_adj row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_adj, l2GMap_adj, 1.);
  
 //============= delta_mu-delta_ctrl row  - BEGIN ===============================
 for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) if (sol_actflag_dofs[kdim][i] != 0 ) sol_actflag_dofs[kdim][i] = ineq_flag * c_compl;    
  
  if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_ctrl, sol_actflag_dofs[kdim]); }
 //============= delta_mu-delta_ctrl row  - END ===============================

 //============= delta_mu-delta_mu row - BEGIN ===============================
  for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) sol_actflag_dofs[kdim][i] =   ineq_flag * (1 - sol_actflag_dofs[kdim][i]/c_compl)  + (1-ineq_flag) * 1.;  //can do better to avoid division, maybe use modulo operator 

  if (assembleMatrix) {    KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_mu, sol_actflag_dofs[kdim] );  }
 //============= delta_mu-delta_mu row - END ===============================
  
  
       }//kdim
       

 }



///@todo This is being added to a weak form?
///this is the same for volume or boundary
 void add_one_times_mu_res_ctrl(const unsigned iproc,
                         const double ineq_flag,
                         const std::vector<unsigned int> pos_ctrl_in_Mat,
                         const std::vector<unsigned int> pos_mu_in_Mat,
                         const vector < unsigned > & SolIndex,
                         const Solution*                sol,
                         const NonLinearImplicitSystemWithPrimalDualActiveSetMethod * mlPdeSys,
                         const  LinearEquationSolver* pdeSys,
                         NumericVector* RES) {
  
     assert(pos_ctrl_in_Mat.size() == pos_mu_in_Mat.size());
     
    const unsigned int   n_components_ctrl = pos_ctrl_in_Mat.size();
    
       for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
     
 const unsigned int ctrl_index_in_Mat = pos_ctrl_in_Mat[kdim];
 const unsigned int mu_index_in_Mat   = pos_mu_in_Mat[kdim];
 const unsigned int mu_index_in_Sol   = pos_mu_in_Mat[kdim];

  unsigned int ctrl_size_iproc_in_Mat = pdeSys->KKoffset[ctrl_index_in_Mat + 1][iproc] - pdeSys->KKoffset[ctrl_index_in_Mat][iproc];
  unsigned int mu_size_iproc_in_Sol = (*sol->_Sol[ SolIndex[mu_index_in_Sol] ]).last_local_index() - (*sol->_Sol[ SolIndex[mu_index_in_Sol] ]).first_local_index(); // pdeSys->KKoffset[mu_index_in_Mat + 1][iproc] - pdeSys->KKoffset[mu_index_in_Mat][iproc];

  assert(ctrl_size_iproc_in_Mat == mu_size_iproc_in_Sol);

  std::vector<double>  one_times_mu(ctrl_size_iproc_in_Mat, 0.);
  std::vector<int>    positions_ctrl_in_Res(ctrl_size_iproc_in_Mat);
  std::vector<int>    positions_mu_in_Sol(ctrl_size_iproc_in_Mat);      

  for (unsigned i = 0; i < positions_ctrl_in_Res.size(); i++) {
    positions_ctrl_in_Res[i] = pdeSys->KKoffset[ctrl_index_in_Mat][iproc] + i;
    positions_mu_in_Sol[i]   = (*sol->_Sol[ SolIndex[mu_index_in_Sol] ]).first_local_index()/*pdeSys->KKoffset[mu_index_in_Mat][iproc]*/ + i;
    //this should not come from pdeSys but from Sol//actually I can take it from the Numeric Vector! 
    ///@todo put the Dof range for Sol 
//          unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);  //this needs iel, in fact it is only in an ELEMENT loop, but here I am in a NODE loop
    
    one_times_mu[i] = ineq_flag * 1. * (*sol->_Sol[ SolIndex[mu_index_in_Sol] ])(positions_mu_in_Sol[i]) ;
    }
    
    RES->add_vector_blocked(one_times_mu, positions_ctrl_in_Res);
    
    
       }
       
    
 }
 

 
  void store_act_flag_in_old(  const NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys,
                               const MultiLevelSolution *    ml_sol,
                              Solution *                sol,
                             std::vector<unsigned int>  & solIndex_act_flag) {
      
     const unsigned int   n_components_ctrl = mlPdeSys->GetActiveSetFlagName().size();
     
          for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
   
  const std::string act_flag_name = mlPdeSys->GetActiveSetFlagName()[kdim];
  
  solIndex_act_flag[kdim]  = ml_sol->GetIndex(act_flag_name.c_str());
  
     if(sol->GetSolutionTimeOrder(solIndex_act_flag[kdim]) == 2) {
       *(sol->_SolOld[solIndex_act_flag[kdim]]) = *(sol->_Sol[solIndex_act_flag[kdim]]);
     }
     
  }
  
  }
  
}



//*******************************************************************************************
//*********************** Mesh independent - END *****************************************
//*******************************************************************************************

  
   
 
  } //end namespace ctrl




  
  
  
} //end namespace
 
#endif
