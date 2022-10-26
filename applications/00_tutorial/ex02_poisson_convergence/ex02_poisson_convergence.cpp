//compute convergence using Cauchy convergence test
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv
//

/** tutorial/Ex2
 * This example shows how to set and solve the weak form of the Poisson problem
 *                    $$ \Delta u = 1 \text{ on }\Omega, $$
 *          $$ u=0 \text{ on } \Gamma, $$
 * on a square domain $\Omega$ with boundary $\Gamma$;
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "LinearImplicitSystem.hpp"
#include "Files.hpp"
#include "FE_convergence.hpp"

#include "../tutorial_common.hpp"


#include "adept.h"


// command to view matrices
// ./tutorial_ex2_c -mat_view ::ascii_info_detail
// ./tutorial_ex2_c -mat_view > matview_print_in_file.txt

using namespace femus;



double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char * name) {

Math::Function< double > *  exact_sol =  ml_prob->get_ml_solution()->get_analytical_function(name);

double value = exact_sol->value(x);

   return value;   

}


bool Solution_set_boundary_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char * name, double& value, const int faceName, const double time) {

    bool dirichlet = true; //dirichlet

    Math::Function< double > *  exact_sol =  ml_prob->get_ml_solution()->get_analytical_function(name);
    value = exact_sol->value(x);
//     value = 0.;

    return dirichlet;
}


template < class system_type, class real_num, class real_num_mov >
void System_assemble_interface(MultiLevelProblem & ml_prob);



// template < class system_type, class real_num, class real_num_mov >
// void System_assemble_flexible_Laplacian_With_Manufactured_Sol(const std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > > & elem_all,
//                               const std::vector<Gauss> & quad_rules,
//                               system_type * mlPdeSys,
//                               MultiLevelMesh * ml_mesh,
//                               MultiLevelSolution * ml_sol,
//                               const std::vector< Unknown > &  unknowns,
//                               const std::vector< Math::Function< double > * > & exact_sol);





// this is the class that I pass to the FE_convergence functions. Basically I pass a class, instead of passing a function pointer.
// notice that this is a class template, not a class with a function template
template < class real_num >
class Solution_generation_1 : public Solution_generation_single_level {

public:

    const MultiLevelSolution  run_on_single_level(MultiLevelProblem & ml_prob,
                                                  MultiLevelMesh & ml_mesh,
                                                  const unsigned i,
                                                  const std::vector< Unknown > & unknowns,
                                                  const std::vector< Math::Function< double > * > &  exact_sol,
                                                  const MultiLevelSolution::InitFuncMLProb      SetInitialCondition_in,
                                                  const MultiLevelSolution::BoundaryFuncMLProb  SetBoundaryCondition_in,
                                                  const bool equation_solve
                                                  ) const;

};






int main(int argc, char** args) {

    // ======= Init ==========================
    FemusInit mpinit(argc, args, MPI_COMM_WORLD);

    // ======= Problem ========================
    MultiLevelProblem ml_prob;

    // ======= Files - BEGIN  =========================
    Files files;
    const bool use_output_time_folder = false;
    const bool redirect_cout_to_file = true;
    files.CheckIODirectories(use_output_time_folder);
    files.RedirectCout(redirect_cout_to_file);

    ml_prob.SetFilesHandler(&files);
    // ======= Files - END  =========================

    // ======= Quad Rule - BEGIN ========================
    std::string fe_quad_rule("seventh");
    /* this is the order of accuracy that is used in the gauss integration scheme
       In the future it is not going to be an argument of the mesh function   */
    
    ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule);
    ml_prob.set_all_abstract_fe_multiple();
    // ======= Quad Rule - END ========================


    // ======= Mesh, Coarse - BEGIN ========================
    MultiLevelMesh ml_mesh;

//     std::string input_file = "square_0-1x0-1_divisions_2x2.med"; // @todo does not work with biquadratic exact solution
//     std::string input_file = "square_0-1x0-1_divisions_5x4.med"; // @todo does not work with biquadratic exact solution 
//     std::string input_file = "../../../../unittests/test_mesh_read_write/input/salome_parametric_with_notebook/square_0-1x0-1_divisions_5x3.med";  // @todo does not work with biquadratic exact solution 
    std::string input_file = "../../../../unittests/test_mesh_read_write/input/salome_parametric_with_notebook/square_0-1x0-1_divisions_2x2_unstructured.med";  // @todo WORKS with biquadratic exact solution
//      std::string input_file = "assignment_square_regular_triangular.med";    // @todo WORKS with biquadratic exact solution

//     std::string input_file = "L_shaped_domain.med";
//     std::string input_file = "interval.med";
//     std::string input_file = "cylinder_hexahedral.med";
    std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
    const std::string infile = mystream.str();

    const double Lref = 1.;
  const bool read_groups = true;
  const bool read_boundary_groups = true;
    
//   ml_mesh.ReadCoarseMesh(infile.c_str(), fe_quad_rule.c_str(), Lref, read_groups, read_boundary_groups);

// // //   const unsigned mesh_file_type = 0;
// // //    std::string mesh_name = Domain_square_m05p05::quad9_all_mesh_generation_methods_Structured(mesh_file_type, ml_mesh, fe_quad_rule);

  
  
  
  ml_mesh.ReadCoarseMeshFileReadingBeforePartitioning(infile.c_str(), Lref, read_groups, read_boundary_groups);
    
  ml_mesh.GetLevelZero(0)->build_dofmap_all_fe_families_and_elem_and_node_structures();
 

  ml_mesh.BuildFETypesBasedOnExistingCoarseMeshGeomElements();
  
  ml_mesh.PrepareNewLevelsForRefinement();

    // ======= Mesh, Coarse - END ========================
  

    // ======= Unknowns - BEGIN ========================
    std::vector< Unknown > unknowns = systems__provide_list_of_unknowns_lagrangian();
    // ======= Unknowns - END ========================



    // ======= Normal run (without convergence study) ========================
    Solution_generation_1< double  /*adept::adouble*/ > my_solution_generation;
//     const unsigned int n_levels = 3;
//     my_solution_generation.run_on_single_level( ... );

    // ======= Convergence study - BEGIN ========================
    
    // Auxiliary mesh, all levels - BEGIN  ================
    unsigned max_number_of_meshes = 6;
    if (ml_mesh.GetDimension() == 3) max_number_of_meshes = 5;

    ///set coarse storage mesh (
    ///@todo should write the copy constructor or "=" operator to copy the previous mesh) ==================
    // If you try to use the default copy constructor it doesn't work.
    // In fact, the copy constructor will copy ALL THE POINTERS, and if there are pointers that were dynamically allocated with new, and destroyed with delete,
    //     new will be called only once but delete will be called twice...
//     MultiLevelMesh ml_mesh_all_levels( ml_mesh);
    MultiLevelMesh ml_mesh_all_levels;

    ml_mesh_all_levels.ReadCoarseMesh(infile.c_str(), fe_quad_rule.c_str(), Lref, read_groups, read_boundary_groups);
// // //     Domain_square_m05p05::quad9_all_mesh_generation_methods_Structured(mesh_file_type, ml_mesh_all_levels, fe_quad_rule);
    // Auxiliary mesh, all levels - END  ================


    // 1) Which exact solution - BEGIN ================
    std::vector< Math::Function< double > * > analytical_function( unknowns.size() );         ///@todo you have to switch it below too, or maybe pass it to MultiLevelProblem  provide exact solution, if available =

    Domain_square_01by01::Function_Zero_on_boundary_1< double >  analytical_function_1;
//     Domain_square_01by01::Function_Zero_on_boundary_2< double >  analytical_function_1;
//     Domain_square_m05p05::Function_Zero_on_boundary_4< double >  analytical_function_1;

    for (unsigned int u = 0; u < analytical_function.size(); u++) {
    analytical_function[u] =  & analytical_function_1;
    }
    // 1) Which exact solution - END ================

    // 2) Choose how to compute the convergence order - BEGIN ============== //0: incremental 1: absolute (with analytical sol)  2: absolute (with projection of finest sol)...    
     const unsigned   conv_order_flag = 0;
    // 2) Choose how to compute the convergence order - END ============== 

     // 3) Choose what norms to compute - BEGIN  (//0 = only L2: //1 = L2 + H1) ==============
     const unsigned norm_flag = 1;
     // 3) Choose what norms to compute - END  ==============

     // 4) Solve Equation or only Approximation Theory - BEGIN   ==============
       const bool equation_solve = true; 
     // 4) Solve Equation or only Approximation Theory  - END   ==============

     // object ================
    FE_convergence<>  fe_convergence;

// we are going to do one Convergence Study for each System. This will give more flexibility when we export this to an arbitrary Application   
//     for (unsigned int u = 0; u < unknowns.size(); u++) {

    for (unsigned int vb = 0; vb < 1; vb++) { //0: volume, 1: boundary, ...
        
    fe_convergence.convergence_study(ml_prob, 
                                     ml_mesh, 
                                     ml_mesh_all_levels, 
                                     max_number_of_meshes, 
                                     norm_flag,
                                     conv_order_flag,
                                     vb,
                                     equation_solve,
                                     my_solution_generation, 
                                     unknowns/*[u]*/,
                                     analytical_function,
                                     Solution_set_initial_conditions,
                                     Solution_set_boundary_conditions
                                    );

    }
    // ======= Convergence study - END ========================
//     }
    
    
    return 0;

}









template < class real_num >
const MultiLevelSolution  Solution_generation_1< real_num >::run_on_single_level(MultiLevelProblem & ml_prob,
                                                                                MultiLevelMesh & ml_mesh_single_level,
                                                                                const unsigned lev,                                                                                const std::vector< Unknown > &  unknowns,
                                                                                const std::vector< Math::Function< double > * > &  exact_sol,
                                                                                const MultiLevelSolution::InitFuncMLProb SetInitialCondition_in,
                                                                                const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition_in,
                                                                                const bool equation_solve

) const {


    //Mesh - BEGIN   ==================
    unsigned numberOfUniformLevels = lev + 1;
    unsigned numberOfSelectiveLevels = 0;
    ml_mesh_single_level.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
    ml_mesh_single_level.EraseCoarseLevels(numberOfUniformLevels - 1);

    ml_mesh_single_level.PrintInfo();
    //Mesh - END   ==================


    //Solution - BEGIN  ==================
    MultiLevelSolution ml_sol_single_level(& ml_mesh_single_level);

    ml_sol_single_level.SetWriter(VTK);
    ml_sol_single_level.GetWriter()->SetDebugOutput(true);

    // ======= Problem ========================
    ml_prob.SetMultiLevelMeshAndSolution(& ml_sol_single_level);
    //Solution - END  ==================
    
    ml_prob.get_systems_map().clear();  //at every lev we'll have a different map of systems

    //only one Problem,
    //only one Mesh per level,
    //only one Solution per level,
    // a separate System for each unknown (I want to do like this to see how to handle multiple coupled systems later on)

    for (unsigned int u = 0; u < unknowns.size(); u++) {

        // ======= Solution, II ==================
        ml_sol_single_level.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);
        ml_sol_single_level.set_analytical_function(unknowns[u]._name.c_str(), exact_sol[u]);   
        ml_sol_single_level.Initialize(unknowns[u]._name.c_str(), SetInitialCondition_in, & ml_prob);

        
// // If you just want an interpolation study, without equation, just initialize every Solution with some function - BEGIN      
       if (equation_solve)  {
       

        ml_sol_single_level.AttachSetBoundaryConditionFunction(SetBoundaryCondition_in);
        ml_sol_single_level.GenerateBdc(unknowns[u]._name.c_str(),  (unknowns[u]._time_order == 0) ? "Steady" : "Time_dependent", & ml_prob);

        // ======= Problem, System - BEGIN ========================
        std::ostringstream sys_name;
        sys_name << unknowns[u]._name;  //give to each system the name of the unknown it solves for!

        LinearImplicitSystem & system = ml_prob.add_system < LinearImplicitSystem > (sys_name.str());

        // ======= System, Unknowns ========================
        system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
        std::vector< Unknown > unknowns_vec(1);
        unknowns_vec[0] = unknowns[u]; //need to turn this into a vector
        system.set_unknown_list_for_assembly(unknowns_vec); //way to communicate to the assemble function, which doesn't belong to any class

         // ======= System, Assemble Function ========================
        system.SetAssembleFunction( System_assemble_interface< LinearImplicitSystem, real_num, double > );

        // ======= System, Exact Solution ========================
//         system.set_exact_solution();
        
       // ======= System, Current number ========================
        ml_prob.set_current_system_number(u);               //way to communicate to the assemble function, which doesn't belong to any class

        // initialize and solve the system
        system.init();
        system.ClearVariablesToBeSolved();
        system.AddVariableToBeSolved("All");

//             system.SetDebugLinear(true);
//             system.SetMaxNumberOfLinearIterations(6);
//             system.SetAbsoluteLinearConvergenceTolerance(1.e-4);

        system.SetOuterSolver(PREONLY/*GMRES*/);
        system.MGsolve();  //everything is stored into the Solution after this
        // ======= Problem, System - END ========================

       }
       
// // If you just want an interpolation study, without equation, just initialize every Solution with some function  - END      

        // ======= Print - BEGIN  ========================
        std::vector < std::string > variablesToBePrinted;
        variablesToBePrinted.push_back(unknowns[u]._name);
        ml_sol_single_level.GetWriter()->Write(unknowns[u]._name, ml_prob.GetFilesHandler()->GetOutputPath(), "biquadratic", variablesToBePrinted, lev);
        // ======= Print - END  ========================

    }


    return ml_sol_single_level;
}



template < class system_type, class real_num, class real_num_mov >
void System_assemble_interface(MultiLevelProblem& ml_prob) {
    //  ml_prob is the global object from/to where get/set all the data

// this is meant to be like a tiny addition to the main function, because we cannot pass these arguments through the function pointer
//what is funny is that this function is attached to a system, but I cannot retrieve the system to which it is attached unless I know the name or the number of it!
//all I have at hand is the MultiLevelProblem, which contains a Vector of Systems
// all I can do is put in the MultiLevelProblem a number that tells me what is the current system being solved

    
    
    // all the arguments here are retrieved through the Multilevel Problem
    
    
 
    const unsigned current_system_number = ml_prob.get_current_system_number();

        // ======= Unknowns - BEGIN  ========================
std::vector< Unknown >  unknowns = ml_prob.get_system< system_type >(current_system_number).get_unknown_list_for_assembly();
        // ======= Unknowns - END  ========================


        // ======= Exact sol - BEGIN  ========================
std::vector< Math::Function< double > * > exact_sol( unknowns.size() );

    for(int u = 0; u < exact_sol.size(); u++) {
        exact_sol[u] = ml_prob.get_ml_solution()->get_analytical_function( unknowns[u]._name.c_str() );
    }
        // ======= Exact sol - END  ========================


        // ======= FE Quadrature - BEGIN  ========================
   //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
        // ======= FE Quadrature - END  ========================

    System_assemble_flexible_Laplacian_With_Manufactured_Sol< system_type, real_num, real_num_mov > (elem_all,
                                                                     ml_prob.GetQuadratureRuleAllGeomElems(),
                                                                     & ml_prob.get_system< system_type >(current_system_number),
                                                                     ml_prob.GetMLMesh(),
                                                                     ml_prob.get_ml_solution(),
                                                                     unknowns,
                                                                     exact_sol);

}






///@todo: check bc for discontinuous FE
///@todo: compute error in L-\infty norm
///@todo: compute nonlinear convergence rate
///@todo: compute time convergence rate, pointwise and then in norms
///@todo: uncouple Gauss from Mesh
///@todo: make non-isoparametric Jacobian routines (abstract Jacobian) for all dims
///@todo: check solver and prec choices
///@todo: create a copy constructor/operator=() for the Mesh
///@todo: check FE convergence in 3D for various operators
///@todo: check face names in 3D
///@todo: test with linear compressible solid first
///@todo: check FE convergence with EXACT solution
///@todo: problems in the P0 or P1 should happen because of the boundary conditions not being enforced
