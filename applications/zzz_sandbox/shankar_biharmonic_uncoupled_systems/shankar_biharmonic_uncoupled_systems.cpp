/** tutorial/Ex3
 * This example shows how to set and solve the weak form of the nonlinear problem
 *                     -\Delta^2 u = f(x) \text{ on }\Omega,
 *            u=0 \text{ on } \Gamma,
 *      \Delta u=0 \text{ on } \Gamma,
 * on a box domain $\Omega$ with boundary $\Gamma$,
 * by using a system of second order partial differential equation.
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"
#include "VTKWriter.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"
#include "CurrentElem.hpp"
#include "ElemType_template.hpp"
#include "adept.h"

#include "FE_convergence.hpp"

//Additional includes
#include "00_system_specifics.hpp"

#include "Solution_functions_over_domains_or_mesh_files.hpp"


using namespace femus;

//Boundary_integral-BEGIN
using namespace femus;

namespace shankar_biharmonic{
    class biharmonic_equation {

    public :

//=======1D_Case-BEGIN==============
static void natural_loop_u_1d (const MultiLevelProblem* ml_prob,
                             const Mesh * msh,
                             const MultiLevelSolution * ml_sol,
                             const unsigned iel,
                             CurrentElem <double> & geom_element,
                             const unsigned xType,
                             const std::string solname_u,
                             const unsigned solFEType_u,
                             std::vector <double> & Res
            ){
    double grad_u_dot_n = 0;
    for (unsigned jface=0; jface <msh->GetElementFaceNumber(iel); jface++) {
        geom_element.set_coords_at_dofs_bdry_3d(iel, jface, xType);
        geom_element.set_elem_center_bdry_3d();
        std::vector <double> xx_face_elem_center(3,0.);
         xx_face_elem_center = geom_element.get_elem_center_bdry_3d();

        const int boundary_index = msh->GetMeshElements()-> GetFaceElementIndex(iel,jface);

        if (boundary_index <0) {

            unsigned int face = -(boundary_index +1);

            bool is_dirichlet = ml_sol-> GetBdcFunctionMLProb()(ml_prob, xx_face_elem_center, solname_u.c_str(), grad_u_dot_n,face,0.);

            if ( !(is_dirichlet)  &&  (grad_u_dot_n != 0.) ) {  //dirichlet == false and nonhomogeneous Neumann



                   unsigned n_dofs_face = msh->GetElementFaceDofNumber(iel, jface, solFEType_u);

                  for (unsigned i = 0; i < n_dofs_face; i++) {

                 unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);

                 Res[i_vol] +=  grad_u_dot_n /* * phi[node] = 1. */;

                         }




                    }

              }

    }
}
static void natural_loop_v_1d (const MultiLevelProblem* ml_prob,
                             const Mesh * msh,
                             const MultiLevelSolution * ml_sol,
                             const unsigned iel,
                             CurrentElem <double> & geom_element,
                             const unsigned xType,
                             const std::string solname_v,
                             const unsigned solFEType_v,
                             std::vector <double> & Res
            ){
    double grad_v_dot_n = 0;
    for (unsigned jface=0; jface <msh->GetElementFaceNumber(iel); jface++) {
        geom_element.set_coords_at_dofs_bdry_3d(iel, jface, xType);
        geom_element.set_elem_center_bdry_3d();
        std:: vector <double> xx_face_elem_center(3,0.);
        xx_face_elem_center = geom_element.get_elem_center_bdry_3d();

        const int boundary_index = msh ->GetMeshElements() -> GetFaceElementIndex(iel, jface);

        if (boundary_index <0) {
            unsigned int face = -(boundary_index + 1);

            bool is_dirichlet = ml_sol -> GetBdcFunctionMLProb()(ml_prob, xx_face_elem_center, solname_v.c_str(), grad_v_dot_n, face,0.);

            if (!(is_dirichlet) && (grad_v_dot_n !=0.)) { //dirichlet == false and nonhomogenous Neumann condition

                unsigned n_dofs_face = msh -> GetElementFaceDofNumber(iel, jface, solFEType_v);

                for (unsigned i = 0; i < n_dofs_face; i++) {
                    unsigned int i_vol = msh -> GetLocalFaceVertexIndex(iel, jface, i);

                    Res[i_vol] += grad_v_dot_n;
                }
            }

        }

    }
  }
  //===========1D_Case-END====================
//===========2D_case-BEGIN==================
   //===========2D_U-BEGIN
template < class real_num, class real_num_mov >
static void natural_loop_u_2d3d(const MultiLevelProblem *    ml_prob,
                       const Mesh *                    msh,
                       const MultiLevelSolution *    ml_sol,
                       const unsigned iel,
                       CurrentElem < double > & geom_element,
                       const unsigned solType_coords,
                       const std::string  solname_u, /*<std::string>*/
                       const unsigned solFEType_u,
                       std::vector< double > & Res,
                       //-----------
                       std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > >  elem_all,
                       const unsigned dim,
                       const unsigned space_dim,
                       const unsigned maxSize
                    ) {
                std::vector <std::vector <double>> JacI_iqp_bdry(space_dim);
                std::vector <std::vector <double>> Jac_iqp_bdry(dim-1);
                for (unsigned d=0; d < Jac_iqp_bdry.size(); d++){
                    Jac_iqp_bdry[d].resize(space_dim);
                }
                for (unsigned d=0; d < JacI_iqp_bdry.size(); d++){
                    JacI_iqp_bdry[d].resize(space_dim);
                }
                double detJac_iqp_bdry;
                double weight_iqp_bdry = 0.;
                std::vector <double> phi_u_bdry;
                std::vector <double> phi_u_x_bdry;

                phi_u_bdry.reserve(maxSize);
                phi_u_x_bdry.reserve(maxSize *space_dim);

                std::vector <double> phi_coords_bdry;
                std::vector <double> phi_coords_x_bdry;

                phi_coords_bdry.reserve(maxSize);

                phi_coords_x_bdry.reserve(maxSize * space_dim);

                double grad_u_dot_n = 0.;
                for (unsigned jface = 0; jface < msh -> GetElementFaceNumber(iel); jface++) {
                    geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);

                    geom_element.set_elem_center_bdry_3d();

                    const unsigned ielGeom_bdry = msh-> GetElementFaceType(iel, jface);


                    std::vector <double> xx_face_elem_center(3,0);
                    xx_face_elem_center = geom_element.get_elem_center_bdry_3d();

                    const int boundary_index = msh->GetMeshElements() -> GetFaceElementIndex(iel, jface);


                    if (boundary_index <0) {

                        unsigned int face = -(boundary_index +1);


                        bool is_dirichlet =  ml_sol->GetBdcFunctionMLProb()(ml_prob, xx_face_elem_center, solname_u.c_str(), grad_u_dot_n, face, 0.);

                        if (!(is_dirichlet)) {
                            unsigned n_dofs_face_u = msh -> GetElementFaceDofNumber(iel, jface, solFEType_u);

                            std::vector <double> grad_u_dot_n_at_dofs(n_dofs_face_u);

                            for (unsigned i_bdry = 0; i_bdry <grad_u_dot_n_at_dofs.size(); i_bdry++){

                                std::vector <double> x_at_node(dim,0.);
                                for (unsigned jdim=0; jdim < x_at_node.size(); jdim++){
                                    x_at_node[jdim] = geom_element.get_coords_at_dofs_bdry_3d()[jdim][i_bdry];
                                }

                                double grad_u_dot_n_at_dofs_temp = 0.;

                                ml_sol-> GetBdcFunctionMLProb()(ml_prob, x_at_node, solname_u.c_str(), grad_u_dot_n_at_dofs_temp, face, 0.);

                                grad_u_dot_n_at_dofs[i_bdry] = grad_u_dot_n_at_dofs_temp;
                            }

                            const unsigned n_guass_bdry = ml_prob-> GetQuadratureRule(ielGeom_bdry).GetGaussPointsNumber();

                    for (unsigned ig_bdry = 0; ig_bdry < n_guass_bdry; ig_bdry++) {
                        elem_all[ielGeom_bdry][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_iqp_bdry, JacI_iqp_bdry, detJac_iqp_bdry, space_dim);

                        weight_iqp_bdry = detJac_iqp_bdry * ml_prob->GetQuadratureRule(ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];
    elem_all[ielGeom_bdry][solFEType_u ]->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_u_bdry, phi_u_x_bdry,  boost::none, space_dim);

                        elem_all[ielGeom_bdry][solType_coords ]->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_coords_bdry, phi_coords_x_bdry,  boost::none, space_dim);

                std::vector<double> x_qp_bdry(dim, 0.);

         for (unsigned i = 0; i < phi_coords_bdry.size(); i++) {
           	for (unsigned d = 0; d < dim; d++) {
 	                                                x_qp_bdry[d]    += geom_element.get_coords_at_dofs_bdry_3d()[d][i] * phi_coords_bdry[i]; // fetch of coordinate points
             }
         }


           double grad_u_dot_n_qp = 0.;
           for (unsigned i_bdry = 0; i_bdry < phi_u_bdry.size(); i_bdry ++) {
           grad_u_dot_n_qp +=  grad_u_dot_n_at_dofs[i_bdry] * phi_u_bdry[i_bdry];
           }

           for (unsigned i_bdry = 0; i_bdry < n_dofs_face_u; i_bdry++) {

                 unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

                 Res[i_vol] +=  weight_iqp_bdry * grad_u_dot_n_qp /*grad_u_dot_n*/  * phi_u_bdry[i_bdry];

                           }

                        }


                }


            }
                }
}

   //============2D_U-END===================


   //============2D_V-BEGIN=================
template <class real_num, class real_num_mov >
static void natural_loop_V_2d3d(const MultiLevelProblem *    ml_prob,
                       const Mesh *                    msh,
                       const MultiLevelSolution *    mlSol,
                       const unsigned iel,
                       CurrentElem < double > & geom_element,
                       const unsigned solType_coords,
                       /*const std::vector <std::string> solname_v*/
                       std::string solname_v,
                       const unsigned solvType,
                       std::vector< double > & Res,
                       //-----------
                       std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > >  elem_all,
                       const unsigned dim,
                       const unsigned int space_dim,
                       const unsigned maxSize
                    ) {

       std::vector < std::vector < double > >  JacI_iqp_bdry(space_dim);
       std::vector < std::vector < double > >  Jac_iqp_bdry(dim-1);
       for (unsigned d = 0; d < Jac_iqp_bdry.size(); d++) {
           Jac_iqp_bdry[d].resize(space_dim); }
       for (unsigned d = 0; d < JacI_iqp_bdry.size(); d++) {
           JacI_iqp_bdry[d].resize(dim-1); }

       double detJac_iqp_bdry;
       double weight_iqp_bdry = 0.;

       //boundary state shape functions
      std::vector <double> phi_v_bdry;
      std::vector <double> phi_v_x_bdry;

      phi_v_bdry.reserve(maxSize);
      phi_v_x_bdry.reserve(maxSize * space_dim);

      std::vector <double> phi_coords_bdry;
      std::vector <double> phi_coords_x_bdry;

      phi_coords_bdry.reserve(maxSize);
      phi_coords_x_bdry.reserve(maxSize * space_dim);
      //boundary state shape functions

      double grad_v_dot_n = 0.;

      for (unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {

       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);

       geom_element.set_elem_center_bdry_3d();

       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, jface);


       std::vector <  double > xx_face_elem_center(3, 0.);
       xx_face_elem_center = geom_element.get_elem_center_bdry_3d();

       const int boundary_index = msh->GetMeshElements()->GetFaceElementIndex(iel, jface);

       if ( boundary_index < 0) { //I am on the boundary

         unsigned int face = - (boundary_index + 1);

         bool is_dirichlet =  mlSol->GetBdcFunctionMLProb()(ml_prob, xx_face_elem_center, solname_v.c_str(), grad_v_dot_n, face, 0.);
         //we have to be careful here, because in GenerateBdc those coordinates are passed as NODE coordinates,
         //while here we pass the FACE ELEMENT CENTER coordinates.
         // So, if we use this for enforcing space-dependent Dirichlet or Neumann values, we need to be careful!

             if ( !(is_dirichlet) /* &&  (grad_u_dot_n != 0.)*/ ) {  //dirichlet == false and nonhomogeneous Neumann

    unsigned n_dofs_face_v = msh->GetElementFaceDofNumber(iel, jface, solvType);

// dof-based - BEGIN
     std::vector< double > grad_v_dot_n_at_dofs(n_dofs_face_v);


    for (unsigned i_bdry = 0; i_bdry < grad_v_dot_n_at_dofs.size(); i_bdry++) {
        std::vector<double> x_at_node(dim, 0.);
        for (unsigned jdim = 0; jdim < x_at_node.size(); jdim++) x_at_node[jdim] = geom_element.get_coords_at_dofs_bdry_3d()[jdim][i_bdry];

      double grad_v_dot_n_at_dofs_temp = 0.;
      mlSol->GetBdcFunctionMLProb()(ml_prob, x_at_node, solname_v.c_str(), grad_v_dot_n_at_dofs_temp, face, 0.);
     grad_v_dot_n_at_dofs[i_bdry] = grad_v_dot_n_at_dofs_temp;

    }

// dof-based - END
    const unsigned n_gauss_bdry = ml_prob->GetQuadratureRule(ielGeom_bdry).GetGaussPointsNumber();


		for(unsigned ig_bdry = 0; ig_bdry < n_gauss_bdry; ig_bdry++) {

     elem_all[ielGeom_bdry][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_iqp_bdry, JacI_iqp_bdry, detJac_iqp_bdry, space_dim);
//      elem_all[ielGeom_bdry][solType_coords]->compute_normal(Jac_iqp_bdry, normal);

    weight_iqp_bdry = detJac_iqp_bdry * ml_prob->GetQuadratureRule(ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];

    elem_all[ielGeom_bdry][solvType]->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_v_bdry, phi_v_x_bdry,  boost::none, space_dim);



//---------------------------------------------------------------------------------------------------------

     elem_all[ielGeom_bdry][solType_coords ]->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_coords_bdry, phi_coords_x_bdry,  boost::none, space_dim);

  std::vector<double> x_qp_bdry(dim, 0.);

         for (unsigned i = 0; i < phi_coords_bdry.size(); i++) {
           	for (unsigned d = 0; d < dim; d++) {
 	                                                x_qp_bdry[d]    += geom_element.get_coords_at_dofs_bdry_3d()[d][i] * phi_coords_bdry[i]; // fetch of coordinate points
             }
         }

           double grad_v_dot_n_qp = 0.;  ///@todo here we should do a function that provides the gradient at the boundary, and then we do "dot n" with the normal at qp

// dof-based
         for (unsigned i_bdry = 0; i_bdry < phi_v_bdry.size(); i_bdry ++) {
           grad_v_dot_n_qp +=  grad_v_dot_n_at_dofs[i_bdry] * phi_v_bdry[i_bdry];
         }

// quadrature point based
 // // // ml_sol->GetBdcFunctionMLProb()(ml_prob, x_qp_bdry, solname_u.c_str(), grad_u_dot_n_qp, face, 0.);

//---------------------------------------------------------------------------------------------------------





                  for (unsigned i_bdry = 0; i_bdry < n_dofs_face_v; i_bdry++) {

                 unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

                 Res[i_vol] +=  weight_iqp_bdry * grad_v_dot_n_qp /*grad_u_dot_n*/  * phi_v_bdry[i_bdry];

                           }



                        }
   }

       }
      }
   }



   //============2D_V-END====================

//===========2D_case-END==================
    };


};

//Boundary_integral-END



//================MANUFACTURED_SOLUTION-BEGIN===================
double GetExactSolutionValue(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char * name) {

    Math::Function< double > *  exact_sol =  ml_prob->get_ml_solution()->get_analytical_function(name);

    double value = exact_sol->value(x);

   return value;
};


void GetExactSolutionGradient(const MultiLevelProblem * ml_prob, std::vector < double >& x, std::vector < double >& gradientValue, const char * name) {

    Math::Function <double > * gradient_value = ml_prob -> get_ml_solution() -> get_analytical_function(name);

    gradientValue[0] = gradient_value -> gradient(x)[0];

    gradientValue[1] = gradient_value -> gradient(x)[1];


};


double GetExactSolutionLaplace(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char * name) {

    Math::Function< double > *  laplace_value =  ml_prob->get_ml_solution()->get_analytical_function(name);


    double LaplacianValue = laplace_value -> laplacian(x);

    return LaplacianValue;

};


/*
// for v - BEGIN ----
double LaplaceGetExactSolutionValue(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char * name) {

    double pi = acos(-1.);

    Math::Function< double > *  laplace_value =  ml_prob->get_ml_solution()->get_analytical_function(name);

    double LaplaceExactSolution = laplace_value -> laplacian(x);

    return -2.* pi * pi *LaplaceExactSolution;

};


void LaplaceGetExactSolutionGradient(const MultiLevelProblem * ml_prob, std::vector < double >& x, std::vector < double >& gradientValue, const char * name) {

    double pi = acos(-1.);

    Math::Function <double > * gradient_value = ml_prob -> get_ml_solution() -> get_analytical_function(name);

    gradientValue[0] = gradient_value -> gradient(x)[0];

    gradientValue[1] = gradient_value -> gradient(x)[1];

    std::vector <double> LaplacegradientValue;

    LaplacegradientValue[0] = -2.* pi * pi *gradientValue[0];
    LaplacegradientValue[1] = -2.* pi * pi *gradientValue[1];

};
//for v - END ----*/

//================MANUFACTURED_SOLUTION-END===================


//===========Boundary_Conditions-BEGIN================

//===========NONHOMOGENOUS_Dirichlet_BC-BEGIN====================
namespace boundaryConditions {

    namespace NonhomogenousCondition{

bool SetBoundaryConditionNonhomogenous(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& Value, const int facename, const double time) {
  bool dirichlet = true;

  if (!strcmp(SolName, "u")) {
      Math::Function <double> * u = ml_prob -> get_ml_solution() -> get_analytical_function(SolName);
      // strcmp compares two string in lexiographic sense.
    Value = u -> value(x);
  }
  else if (!strcmp(SolName, "v")) {
      Math::Function <double> * v = ml_prob -> get_ml_solution() -> get_analytical_function(SolName);
    Value = v -> value(x);
  }

  return dirichlet;
}
    }

}

//===========NONHOMOGENOUS_Dirichlet__BC-END====================


//===========Boundary_Conditions-END================


void AssembleU_AD(MultiLevelProblem& ml_prob);
void AssembleV_AD(MultiLevelProblem& ml_prob);


// template <class system_type, class real_num, class real_num_mov >


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  //Multilevelproblem
  MultiLevelProblem ml_prob;





  //=======Files-BEGIN=======
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = false;
  Files files;
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);


  //======Problem, files-BEGIN========
  ml_prob.SetFilesHandler(& files);
   //======Problem, files-END========


  //=========Problem, Quad Rule-BEGIN=======
  std::string fe_quad_rule("seventh");
  //=========Problem, Quad Rule-END=========

  //Start defining applications for each manufactured solution
//   //==========APP_SPECIFICS-BEGIN===========
//   std::string system_common_name = "Poisson";
//
//   std::vector <system_specifics>  my_specifics;
//
//   system_specifics app_square_01by01;
//
//   system_specifics app_square_m05p05;
//
//   system_specifics app_square_01by01_2;
//
//   system_specifics  app_cylinder;
//
// //   system_specifics app_circle;
//
//   const std::string relative_path_to_build_directory =  "../../../";
//
//   //==========app_square_01by01-BEGIN=======
//   app_square_01by01._system_name = "Poisson1";
//
//   app_square_01by01._mesh_files.push_back("square_-0p5-0p5x-0p5-0p5_divisions_2x2.med");
//
//   app_square_01by01._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/square/minus0p5-plus0p5_minus0p5-plus0p5/");
//
//   Domains::square_01by01::Function_Zero_on_boundary_1<>  app_square_function_zero_on_boundary_1;
//
//   app_square_01by01._true_solution_function = & app_square_function_zero_on_boundary_1;
//
//   app_square_01by01._assemble_function_for_rhs = & app_square_function_zero_on_boundary_1;
//
//   app_square_01by01._boundary_conditions_types_and_values = boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;
//
//   //==========app_square_01by01-END=======
//
//
//
//   //==========app_square_m05p05-BEGIN=======
//   app_square_m05p05._system_name = "Poisson2";
//
//   app_square_m05p05._mesh_files.push_back("square_-0p5-0p5x-0p5-0p5_divisions_2x2.med");
//
//   app_square_m05p05._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/square/minus0p5-plus0p5_minus0p5-plus0p5/");
//
//   Domains::square_01by01::Function_Zero_on_boundary_4<>  app_square_function_zero_on_boundary_4;
//
//   app_square_m05p05._true_solution_function = & app_square_function_zero_on_boundary_4;
//
//   app_square_m05p05._assemble_function_for_rhs = & app_square_function_zero_on_boundary_4;
//
//   app_square_m05p05._boundary_conditions_types_and_values = boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;
//
//   //==========app_square_01by01-END=======
//
//
//
//
//   //===========app_square_01by01_2-BEGIN======
//   app_square_01by01_2._system_name = "Poisson3";
//
//   app_square_01by01_2._mesh_files.push_back("square_0-1x0-1_divisions_2x2.med");
//
//   app_square_01by01_2._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/square/0-1x0-1/");
//
//   Domains::square_01by01::Function_Zero_on_boundary_2<>  app_square_function_zero_on_boundary_2;
//
//   app_square_01by01_2._true_solution_function = & app_square_function_zero_on_boundary_2;
//
//   app_square_01by01_2._assemble_function_for_rhs = & app_square_function_zero_on_boundary_2;
//
//   app_square_01by01_2._boundary_conditions_types_and_values = boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;
//
//   //===========app_square_01by01_2-END======
//
//
// /*
//   //===========App_circle-BEGIN========
//
//   app_circle._mesh_files.push_back("disk_quad.med");
//
//   app_circle._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/circle/");
//
//   Domains::circle::Function_Zero_on_boundary_1<>  app_circle_function_zero_on_boundary_1;
//
//   app_circle._true_solution_function = & app_circle_function_zero_on_boundary_1;
//
//   app_circle._assemble_function_for_rhs = & app_circle_function_zero_on_boundary_1;
//
//   //===========App_circle-END==========*/
//
//
//
//
// /*
//  //============App_cylinder-BEGIN========
//   app_cylinder._mesh_files.push_back("assignment_cylinder_hexahedral.med");
//
//   app_cylinder._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/3d/cylinder/");
//
//   Domains::cylinder_along_z_with_base_centered_at_1_by_1::Function_Zero_on_boundary_1<>  app_cylinder_function_zero_on_boundary_1;
//
//   app_cylinder._true_solution_function = & app_cylinder_function_zero_on_boundary_1;
//
//   app_cylinder._assemble_function_for_rhs = & app_cylinder_function_zero_on_boundary_1;
//
//   //===========App_cylinder-END==========*/
//
//
//
//
//
//   my_specifics.push_back(app_square_01by01);
//
//   my_specifics.push_back(app_square_m05p05);
//
//    my_specifics.push_back(app_square_01by01_2);
//
// //   my_specifics.push_back(app_cylinder);
//
// //   my_specifics.push_back(app_circle);
//
// //===========APP_SPECIFICS-END============



    // ======= App Specifics for all mesh files - BEGIN  ==================
//   std::string system_common_name = "Laplace";
  std::vector< system_specifics >   my_specifics;

  system_specifics  app_segment;   //me

  system_specifics  app_square;   //me
//   system_specifics  app_circle;   //Gayani
  system_specifics  app_semicircle;   //Himali
  system_specifics  app_quarter_circle;      //Max
  system_specifics  app_annulus;       //Jon
  system_specifics  app_semiannulus;        //Fahad

  system_specifics  app_cube;   //me
  system_specifics  app_cylinder;            //Aman
//   system_specifics  app_semicylinder;   //Rifat
  system_specifics  app_quarter_cylinder; //Armando   //coarser mesh to be done
  system_specifics  app_prism_annular_base;  //Abu

  const std::string relative_path_to_build_directory =  "../../../";

  //segment - BEGIN
  app_segment._system_name = "Segment";
  app_segment._mesh_files.push_back("segment_16_dir_neu.med");
  app_segment._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/1d/segment/0-1/");

  app_segment._boundary_conditions_types_and_values             = boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;

  Domains::segment_0x1::Function_Zero_on_boundary_1<>   app_segment_function_zero_on_boundary_1;
  app_segment._assemble_function_for_rhs   = & app_segment_function_zero_on_boundary_1;
  app_segment._true_solution_function      = & app_segment_function_zero_on_boundary_1;
  //segment - END


 //assignment_square - BEGIN
  app_square._system_name = "Square";
  app_square._mesh_files.push_back("square_0-1x0-1_divisions_2x2.med");
  app_square._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/square/0-1x0-1/");

  app_square._boundary_conditions_types_and_values             =boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;

  Domains::square_01by01::Function_Zero_on_boundary_1<>   app_square_function_zero_on_boundary_1;
  app_square._assemble_function_for_rhs        = & app_square_function_zero_on_boundary_1;
  app_square._true_solution_function           = & app_square_function_zero_on_boundary_1;
 //assignment_square - END


  //assignment_semicircle - BEGIN
  app_semicircle._system_name = "Semi_circle";
  // app_semicircle._mesh_files.push_back("assignment_semicircle_triangular.med");
  app_semicircle._mesh_files.push_back("assignment_semicircle_quadrilateral.med");
  app_semicircle._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/circle_semi/");

  app_semicircle._boundary_conditions_types_and_values             =boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;

  Domains::semicircle_centered_at_0_by_0::Function_Zero_on_boundary_1<>     app_semicircle_function_zero_on_boundary_1;
  app_semicircle._assemble_function_for_rhs = & app_semicircle_function_zero_on_boundary_1;
  app_semicircle._true_solution_function    = & app_semicircle_function_zero_on_boundary_1;
  //assignment_semicircle - END


  //assignment_quarter_circle - BEGIN
  app_quarter_circle._system_name = "Quarter_circle";
  app_quarter_circle._mesh_files.push_back("assignment_quarter_circle_triangular.med");
  app_quarter_circle._mesh_files.push_back("assignment_quarter_circle_quadrilateral.med");
  app_quarter_circle._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/circle_quarter/");
  app_quarter_circle._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/circle_quarter/");

  app_quarter_circle._boundary_conditions_types_and_values             =boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;

  Domains::quarter_circle_centered_at_0_by_0::Function_Zero_on_boundary_1<>     app_quarter_circle_function_zero_on_boundary_1;
  app_quarter_circle._assemble_function_for_rhs = & app_quarter_circle_function_zero_on_boundary_1;
  app_quarter_circle._true_solution_function    = & app_quarter_circle_function_zero_on_boundary_1;
  //assignment_quarter_circle - END


  //assignment_annulus - BEGIN
  app_annulus._system_name = "Annulus";
  //   app_annulus._mesh_files.push_back("assignment_annulus_triangular.med");
  app_annulus._mesh_files.push_back("assignment_annulus_quadrilateral.med");
  app_annulus._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/annulus/");

  app_annulus._boundary_conditions_types_and_values             =boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;

  Domains::annulus_centered_at_0_by_0::Function_Zero_on_boundary_1<>     app_annulus_function_zero_on_boundary_1;
  app_annulus._assemble_function_for_rhs = & app_annulus_function_zero_on_boundary_1;
  app_annulus._true_solution_function    = & app_annulus_function_zero_on_boundary_1;
  //assignment_annulus - END


  //assignment_semiannulus - BEGIN
  app_semiannulus._system_name = "Semi_annulus";
  app_semiannulus._mesh_files.push_back("assignment_semiannulus_triangular.med");
  app_semiannulus._mesh_files.push_back("assignment_semiannulus_quadrilateral.med");
  app_semiannulus._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/annulus_semi/");
  app_semiannulus._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/annulus_semi/");

  app_semiannulus._boundary_conditions_types_and_values             =boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;

  Domains::semiannulus_centered_at_0_by_0_cut_along_y::Function_Zero_on_boundary_1<>     app_semiannulus_function_zero_on_boundary_1;
  app_semiannulus._assemble_function_for_rhs = & app_semiannulus_function_zero_on_boundary_1;
  app_semiannulus._true_solution_function    = & app_semiannulus_function_zero_on_boundary_1;
  //assignment_semiannulus - END


 //assignment_cube - BEGIN
  app_cube._system_name = "Cube";
  app_cube._mesh_files.push_back("cube_hexahedral_divisions_2x2x2.med");
  app_cube._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/3d/cube/0-1x0-1x0-1/");

  app_cube._boundary_conditions_types_and_values             =boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;

  Domains::cube_01_by_01_by_01::Function_Zero_on_boundary_1<>   app_cube_function_zero_on_boundary_1;
  app_cube._assemble_function_for_rhs              = & app_cube_function_zero_on_boundary_1;
  app_cube._true_solution_function                 = & app_cube_function_zero_on_boundary_1;
 //assignment_cube - END



  //assignment_cylinder - BEGIN
  app_cylinder._system_name = "Cylinder";
  app_cylinder._mesh_files.push_back("assignment_cylinder_tetrahedral.med");
  app_cylinder._mesh_files.push_back("assignment_cylinder_hexahedral.med");
  app_cylinder._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/3d/cylinder/");
  app_cylinder._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/3d/cylinder/");

  app_cylinder._boundary_conditions_types_and_values             =boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;

  Domains::cylinder_along_z_with_base_centered_at_1_by_1::Function_Zero_on_boundary_1<>   app_cylinder_function_zero_on_boundary_1;
  app_cylinder._assemble_function_for_rhs                                    = & app_cylinder_function_zero_on_boundary_1;
  app_cylinder._true_solution_function                                       = & app_cylinder_function_zero_on_boundary_1;
  //assignment_cylinder - END


  //assignment_quarter_cylinder - BEGIN
  app_quarter_cylinder._system_name = "Quarter_circle";

  app_quarter_cylinder._mesh_files.push_back("assignment_quarter_cylinder_tetrahedral.med");
  app_quarter_cylinder._mesh_files.push_back("assignment_quarter_cylinder_hexahedral.med");
//   app_quarter_cylinder._mesh_files.push_back("assignment_quarter_cylinder_hexahedral_0.med");
//   app_quarter_cylinder._mesh_files.push_back("assignment_quarter_cylinder_hexahedral_1.med");
  app_quarter_cylinder._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/3d/cylinder_quarter/");
  app_quarter_cylinder._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/3d/cylinder_quarter/");

  app_quarter_cylinder._boundary_conditions_types_and_values             =boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;

  Domains::quarter_cylinder_along_z_with_base_centered_at_0_by_0::Function_Zero_on_boundary_1<>   app_quarter_cylinder_function_zero_on_boundary_1;
  app_quarter_cylinder._assemble_function_for_rhs                                    = & app_quarter_cylinder_function_zero_on_boundary_1;
  app_quarter_cylinder._true_solution_function                                       = & app_quarter_cylinder_function_zero_on_boundary_1;
  //assignment_quarter_cylinder - END




  //assignment_tetra_prism_annular_base - BEGIN
  app_prism_annular_base._system_name = "Annular_Base";
  app_prism_annular_base._mesh_files.push_back("assignment_prism_annular_base_tetrahedral.med");
  app_prism_annular_base._mesh_files.push_back("assignment_prism_annular_base_hexahedral.med");
  app_prism_annular_base._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/3d/prism_annular_base/");
  app_prism_annular_base._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/3d/prism_annular_base/");

  app_prism_annular_base._boundary_conditions_types_and_values             =boundaryConditions::NonhomogenousCondition::SetBoundaryConditionNonhomogenous;

  Domains::prism_annular_base_along_z_with_base_centered_at_0_by_0::Function_Zero_on_boundary_1<>   app_prism_annular_base_function_zero_on_boundary_1;
  app_prism_annular_base._assemble_function_for_rhs                                    = & app_prism_annular_base_function_zero_on_boundary_1;
  app_prism_annular_base._true_solution_function                                       = & app_prism_annular_base_function_zero_on_boundary_1;
  //assignment_tetra_prism_annular_base - END







  my_specifics.push_back(app_segment);
  // my_specifics.push_back(app_square);
  // my_specifics.push_back(app_semicircle);
  // my_specifics.push_back(app_quarter_circle);
  // my_specifics.push_back(app_annulus);
  // my_specifics.push_back(app_semiannulus);
  // my_specifics.push_back(app_cube);
  // my_specifics.push_back(app_cylinder);
  // my_specifics.push_back(app_quarter_cylinder);
  // my_specifics.push_back(app_prism_annular_base);

    // ======= App Specifics for all mesh files - END  ==================


  //Application loop
  for (unsigned int app = 0; app < my_specifics.size(); app++) {

      //=======Problem, App===============
      ml_prob.set_app_specs_pointer(&my_specifics[app]);

       // define multilevel mesh
      MultiLevelMesh mlMsh;

      //Mesh Reading loop
      for(unsigned int m=0; m< my_specifics[app]._mesh_files.size(); m++){
          // ======= Mesh, Coarse reading - BEGIN ==================
            double Lref = 1.;

            const bool read_groups = true; //with this being false, we don't read any group at all. Therefore, we cannot even read the boundary groups that specify what are the boundary faces, for the boundary conditions
            const bool read_boundary_groups = true;


            const std::string mesh_file = my_specifics[app]._mesh_files_path_relative_to_executable[m] + my_specifics[app]._mesh_files[m];

            mlMsh.ReadCoarseMesh(mesh_file, Lref, read_groups, read_boundary_groups);
            // ======= Mesh, Coarse reading - END ==================

            // read coarse level mesh and generate finers level meshes
            double scalingFactor = 1.;

            mlMsh.ReadCoarseMesh(mesh_file.c_str(), fe_quad_rule.c_str(), scalingFactor);

            unsigned maxNumberOfMeshes = 5;

            std::vector < std::vector < double > > l2Norm;
            l2Norm.resize(maxNumberOfMeshes);

            std::vector < std::vector < double > > semiNorm;
            semiNorm.resize(maxNumberOfMeshes);

                std::vector<FEOrder> feOrder;
                feOrder.push_back(FIRST);
                feOrder.push_back(SERENDIPITY);
                feOrder.push_back(SECOND);


  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i + 1;
    unsigned numberOfSelectiveLevels = 0;
    mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

    // print mesh info
    mlMsh.PrintInfo();

    l2Norm[i].resize( feOrder.size() );
    semiNorm[i].resize( feOrder.size() );

    for (unsigned j = 0; j < feOrder.size(); j++) {   // loop on the FE Order

      // define the multilevel solution and attach the mlMsh object to it
      MultiLevelSolution mlSol(&mlMsh);

      mlSol.AddSolution("u", LAGRANGE, feOrder[j]);

      mlSol.set_analytical_function("u",  my_specifics[app]._true_solution_function);

      mlSol.AddSolution("v", LAGRANGE, feOrder[j]);

      mlSol.set_analytical_function("v",my_specifics[app]._assemble_function_for_rhs);

      mlSol.Initialize("All");

      // define the multilevel problem attach the mlSol object to it
      MultiLevelProblem ml_prob(&mlSol);

      // Attach the boundary conditions
      mlSol.AttachSetBoundaryConditionFunction(my_specifics[app]._boundary_conditions_types_and_values);
      mlSol.GenerateBdc("u", "Steady", & ml_prob);
      mlSol.GenerateBdc("v", "Steady", & ml_prob);



      // ======= Problem, Quad Rule - BEGIN ========================
      ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule);
      ml_prob.set_all_abstract_fe_AD_or_not();
      // ======= Problem, Quad Rule - END  ========================

      // add system Poisson in ml_prob as a Linear Implicit Sygeom_elementstem
      NonLinearImplicitSystem& systemU = ml_prob.add_system < NonLinearImplicitSystem > ("PoissonU");
      NonLinearImplicitSystem& systemV = ml_prob.add_system < NonLinearImplicitSystem > ("PoissonV");

      // add solution "u" to system
      systemU.AddSolutionToSystemPDE("u");
      systemV.AddSolutionToSystemPDE("v");



      // attach the assembling function to system
      systemU.SetAssembleFunction(AssembleU_AD);
      systemV.SetAssembleFunction(AssembleV_AD);

      // initilaize and solve the system
      systemU.init();
      systemV.init();

      systemV.MGsolve(); //first solve for v
      systemU.MGsolve(); //then solve for u using v

      // convergence for u
//        std::pair< double , double > norm = GetErrorNorm_L2_H1_with_analytical_sol(&mlSol, "u",my_specifics[app]._true_solution_function );

      // // convergence for v
      // std::pair< double , double > norm = GetErrorNorm_L2_H1_with_analytical_sol(&mlSol, "v", LaplaceGetExactSolutionValue, LaplaceGetExactSolutionGradient );

//       l2Norm[i][j]  = norm.first;
//       semiNorm[i][j] = norm.second;

      // print solutions
      const std::string print_order = fe_fams_for_files[ FILES_CONTINUOUS_BIQUADRATIC ];
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");

      VTKWriter vtkIO(&mlSol);
      vtkIO.Write(my_specifics[app]._system_name + "_" + my_specifics[app]._mesh_files[m], files.GetOutputPath(), print_order, variablesToBePrinted);
    }
  }

/*
  // ======= L2 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < feOrder.size(); j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < feOrder.size(); j++) {
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= L2 - END  ========================

  // ======= H1 - BEGIN  ========================

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < feOrder.size(); j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < feOrder.size(); j++) {
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }

  // ======= H1 - END  ========================*/
  }

}
  return 0;
}



// - Delta v = f
void AssembleV_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("PoissonV");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  CurrentElem < double > geom_element(dim, msh);
  constexpr unsigned int space_dim = 3;

  //solution variable
  const std::string solname_v = "v"; //mlSol-> GetSolName_string_vec()[0];


  unsigned solvIndex = mlSol->GetIndex( solname_v.c_str() );    // get the position of "u" in the ml_sol object
  unsigned solvType = mlSol->GetSolutionType(solvIndex);    // get the finite element type for "u"

  unsigned solvPdeIndex;
  solvPdeIndex = mlPdeSys->GetSolPdeIndex( solname_v.c_str() );    // get the position of "u" in the pdeSys object

  std::vector < adept::adouble >  solv; // local solution

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < int > sysDof; // local to global pdeSys dofs
  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  std::vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  std::vector < double > Res; // local redidual vector
  std::vector < adept::adouble > aRes; // local redidual vector


  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solv.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  sysDof.reserve(maxSize);
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve(maxSize * dim2);
  Res.reserve(maxSize);
  aRes.reserve(maxSize);

  std::vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  Jac.reserve(maxSize * maxSize);
  std::vector < double > Jact; // local Jacobian matrix (ordered by raw, PETSC)
  Jact.reserve(maxSize * maxSize);


  KK->zero(); // Set to zero all the entries of the Global Matrix

  std::vector < std::vector < /*const*/ elem_type_templ_base <double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    geom_element.set_coords_at_dofs_and_geom_type(iel, xType);

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, solvType);    // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    sysDof.resize(nDofs);
    solv.resize(nDofs);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofs2);
    }

    Res.resize(nDofs);    //resize
    std::fill(Res.begin(), Res.end(), 0);    //set Res to zero

    aRes.resize(nDofs);    //resize
    std::fill(aRes.begin(), aRes.end(), 0);    //set aRes to zero

    Jact.resize(nDofs * nDofs);
    Jac.resize(nDofs * nDofs);

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel,  solvType);    // global to global mapping between solution node and solution dof
      solv[i] = (*sol->_Sol[solvIndex])(solDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(solvIndex, solvPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }

    //=========== BOUNDARY - BEGIN===========
    if (dim==1) shankar_biharmonic::biharmonic_equation::natural_loop_v_1d(& ml_prob, msh, mlSol, iel, geom_element, xType, solname_v, solvType, Res);

    if (dim==2 || dim==3) shankar_biharmonic::biharmonic_equation::natural_loop_V_2d3d<double, double>(& ml_prob, msh, mlSol,
                      iel, geom_element, xType,
                      solname_v, solvType,
                      Res,
                      elem_all,
                      dim,
                      space_dim,
                      maxSize
                     );;
    //=========== BOUNDARY-END===============


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solvType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solvType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble solvGauss = 0;
      std::vector < adept::adouble > solvGauss_x(dim, 0.);
      std::vector < double > xGauss(dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {
        solvGauss += phi[i] * solv[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          solvGauss_x[jdim] += phi_x[i * dim + jdim] * solv[i];
          xGauss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        adept::adouble Laplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          Laplace   +=  - phi_x[i * dim + jdim] * solvGauss_x[jdim];
        }

        double exactSolValue =  ml_prob. get_ml_solution()-> get_analytical_function("v")->laplacian(xGauss);

        double f = exactSolValue * phi[i] ;
        aRes[i] += (f -  Laplace) * weight;

      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    for (int i = 0; i < nDofs; i++) {
      Res[i] = aRes[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);



    // define the dependent variables
    s.dependent(&aRes[0], nDofs);

    // define the independent variables
    s.independent(&solv[0], nDofs);

    // get the jacobian matrix (ordered by column)
    s.jacobian(&Jac[0]);

    // get the jacobian matrix (ordered by raw, i.e. Jact=Jac^t)
    for (int inode = 0; inode < nDofs; inode++) {
      for (int jnode = 0; jnode < nDofs; jnode++) {
        Jact[inode * nDofs + jnode] = -Jac[jnode * nDofs + inode];
      }
    }

    //store Jact in the global matrix KK
    KK->add_matrix_blocked(Jact, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}


// - Delta u = v
void AssembleU_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("PoissonU");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned solvIndex;
  solvIndex = mlSol->GetIndex("v");    // get the position of "u" in the ml_sol object
  //unsigned solvType = mlSol->GetSolutionType(solvIndex);

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  CurrentElem < double > geom_element(dim, msh);

  const std::string  solname_u = mlSol -> GetSolName_string_vec()[0];

  constexpr unsigned int space_dim = 3;

  std::vector < adept::adouble >  solu; // local solution
  std::vector < double >  solv; // local solution

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector < int > sysDof; // local to global pdeSys dofs
  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  std::vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  std::vector < double > Res; // local redidual vector
  std::vector < adept::adouble > aRes; // local redidual vector


  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);
  solv.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  sysDof.reserve(maxSize);
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve(maxSize * dim2);
  Res.reserve(maxSize);
  aRes.reserve(maxSize);

  std::vector < double > Jac; // local Jacobian matrix (ordered by column, adept)
  Jac.reserve(maxSize * maxSize);
  std::vector < double > Jact; // local Jacobian matrix (ordered by raw, PETSC)
  Jact.reserve(maxSize * maxSize);



  KK->zero(); // Set to zero all the entries of the Global Matrix

  std::vector < std::vector < /*const*/ elem_type_templ_base <double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);



  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    geom_element.set_coords_at_dofs_and_geom_type(iel, xType);

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofs  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofs2 = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    sysDof.resize(nDofs);
    solu.resize(nDofs);
    solv.resize(nDofs);



    //=========Boundary-BEGIN===============
    if (dim==1) shankar_biharmonic::biharmonic_equation::natural_loop_u_1d(&ml_prob, msh, mlSol,
                                                                           iel, geom_element, xType,
                                                                           solname_u, soluType, Res);

    if (dim==2 || dim==3) shankar_biharmonic::biharmonic_equation::natural_loop_u_2d3d< double, double >(& ml_prob, msh, mlSol, iel, geom_element, xType, solname_u, soluType, Res, elem_all, dim, space_dim, maxSize);

    //=========Boundary-END===============

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofs2);
    }

    Res.resize(nDofs);    //resize
    std::fill(Res.begin(), Res.end(), 0);    //set Res to zero

    aRes.resize(nDofs);    //resize
    std::fill(aRes.begin(), aRes.end(), 0);    //set aRes to zero


    Jact.resize(nDofs * nDofs);
    Jac.resize(nDofs * nDofs);


    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofs; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      solv[i] = (*sol->_Sol[solvIndex])(solDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofs2; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble soluGauss = 0;
      std::vector < adept::adouble > soluGauss_x(dim, 0.);

      double solvGauss = 0;

      std::vector < double > xGauss(dim, 0.);

      for (unsigned i = 0; i < nDofs; i++) {
        soluGauss += phi[i] * solu[i];
        solvGauss += phi[i] * solv[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          soluGauss_x[jdim] += phi_x[i * dim + jdim] * solu[i];
          xGauss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofs; i++) {

        adept::adouble Laplace = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          Laplace   +=  - phi_x[i * dim + jdim] * soluGauss_x[jdim];
        }

        double exactSolValue =  ml_prob. get_ml_solution()-> get_analytical_function("u")->laplacian(xGauss);

        const double f = solvGauss/*exactSolValue*/  * phi[i] ;
        aRes[i] += (f - Laplace) * weight;

      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    for (int i = 0; i < nDofs; i++) {
      Res[i] = aRes[i].value();
    }

    RES->add_vector_blocked(Res, sysDof);



    // define the dependent variables
    s.dependent(&aRes[0], nDofs);

    // define the independent variables
    s.independent(&solu[0], nDofs);

    // get the jacobian matrix (ordered by column)
    s.jacobian(&Jac[0]);

    // get the jacobian matrix (ordered by raw, i.e. Jact=Jac^t)
    for (int inode = 0; inode < nDofs; inode++) {
      for (int jnode = 0; jnode < nDofs; jnode++) {
        Jact[inode * nDofs + jnode] = -Jac[jnode * nDofs + inode];
      }
    }

    //store Jact in the global matrix KK
    KK->add_matrix_blocked(Jact, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

}




