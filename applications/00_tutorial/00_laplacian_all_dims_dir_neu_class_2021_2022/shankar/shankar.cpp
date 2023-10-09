/**
 * Solve 
 *     - \Delta u = 1
*/

// library includes
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "NumericVector.hpp"
#include "cmath"

#include "CurrentElem.hpp"
#include "ElemType_template.hpp"

#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"


// application includes
#include "00_poisson_eqn.hpp"


#define LIBRARY_OR_USER   1 //0: library; 1: user

#if LIBRARY_OR_USER == 0
   #include "00_poisson_eqn_with_dirichlet_or_neumann_bc.hpp"
   #define NAMESPACE_FOR_POISSON  femus
#elif LIBRARY_OR_USER == 1
   #include "00_poisson_eqn.hpp"
   #define NAMESPACE_FOR_POISSON  shankar
#endif

 


double InitialValueDS(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {
    
  return 0.;
  
}




// SEGMENT - BEGIN
double segment_dir_neu_fine__laplacian__rhs(const std::vector<double> & x_qp){
    
    // for a 1d segment
    
    return  -2.;
}

 
bool segment_dir_neu_fine__laplacian__bc(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

  bool dirichlet = false;
  value = 0.;
  
  const double tolerance = 1.e-8;
  
 if (ml_prob->GetMLMesh()->GetDimension() == 1 )  {
  
  if (face_name == 1) {
      dirichlet = true;
        value = 0.; //Dirichlet value
    }
  else if (face_name == 2) {
      dirichlet = false;
        value = 1.; //Neumann value
    }

    
 }
 
 
  return dirichlet;
  
 }
// SEGMENT - END


// SQUARE - BEGIN
bool square__laplacian__bc(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

  if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();
  
   bool dirichlet = false;
  value = 0.;
    
     
  if (face_name == 1) {
      dirichlet = true;
        value = 1.5;
  }
  else if (face_name == 2) {
      dirichlet = true;
        value = 1.5;
  }

 else  if (face_name == 3) {
      dirichlet = true;
        value = 1.5;
  }
  else if (face_name == 4) {
      dirichlet = true;
        value = 0.;//x[0] * x[0];
  }
  
  

   return dirichlet;
   
}

double square__laplacian__rhs(const std::vector < double >& x) {
    
  return 8*(M_1_PI*M_1_PI)*cos(2*M_1_PI*x[0])*cos(2*M_1_PI*x[1]);
//return -2. * ( x[0] * (1. - x[0])  + x[1] * (1. - x[1]) );
}

double square__laplacian__true_solution(const std::vector < double >& x) {

    return cos(2*M_1_PI*x[0])*cos(2*M_1_PI*x[1]);
  //return x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
    
}
// SQUARE - END



 


int main(int argc, char** args) {

  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Problem ========================
  MultiLevelProblem ml_prob;

  // ======= Files - BEGIN  ========================
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = false;
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);

  // ======= Problem, Files ========================
  ml_prob.SetFilesHandler(&files);
  // ======= Files - END  ========================

  // ======= Problem, Quad Rule - BEGIN ========================
  std::string fe_quad_rule("seventh");

  ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule);
  ml_prob.set_all_abstract_fe_multiple();
  // ======= Problem, Quad Rule - END  ========================

    // ======= App Specifics - BEGIN  ==================
  app_specifics  app_segment;   //me

  //segment_dir_neu_fine
// <<<<<<< HEAD
  app_segment._mesh_files.push_back("Mesh_2_xy_boundaries_groups_4x4.med");
  

  app_segment._system_name = "Equation";
  app_segment._assemble_function = NAMESPACE_FOR_POISSON::poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;
  
// // //   app_segment._mesh_files.push_back("Mesh_1_x_dir_neu.med");
// // //   app_segment._boundary_conditions_types_and_values             = segment_dir_neu_fine__laplacian__bc;
// // //   app_segment._assemble_function_rhs = segment_dir_neu_fine__laplacian__rhs;
// // // //   app_segment._true_solution    = segment_dir_neu_fine__laplacian__true_solution;  
  
  app_segment._mesh_files.push_back("Mesh_2_xy_boundaries_groups_4x4.med");
  app_segment._boundary_conditions_types_and_values             = square__laplacian__bc;
  app_segment._assemble_function_rhs = square__laplacian__rhs;
//   app_segment._true_solution    = square__laplacian__true_solution; 
  
  ///@todo if this is not set, nothing should happen here

  
    // ======= App Specifics - END  ==================


  // ======= Mesh, files - BEGIN  ==================   
//    mesh_files.push_back("Mesh_1_x_dir_neu.med");
//    mesh_files.push_back("Mesh_2_xy_all_dir.med");
//    mesh_files.push_back("Mesh_2_xy_boundaries_groups_4x4.med");
//    mesh_files.push_back("Mesh_1_x_all_dir.med");
//    mesh_files.push_back("Mesh_1_y_all_dir.med");
//    mesh_files.push_back("Mesh_1_z_all_dir.med");
//    mesh_files.push_back("Mesh_2_xz_all_dir.med");
//    mesh_files.push_back("Mesh_2_yz_all_dir.med");
//    mesh_files.push_back("Mesh_3_xyz_all_dir.med");
//    mesh_files.push_back("dome_tri.med");
//    mesh_files.push_back("dome_quad.med");
//    mesh_files.push_back("disk_quad.med");
//    mesh_files.push_back("disk_quad_45x.med");
//    mesh_files.push_back("disk_quad_90x.med");
//    mesh_files.push_back("disk_tri.med");
//    mesh_files.push_back("disk_tri_45x.med");
//    mesh_files.push_back("disk_tri_90x.med");
  // ======= Mesh, files - END ==================



 for (unsigned int m = 0; m < app_segment._mesh_files.size(); m++)  {
   
  // ======= Problem, App (every Problem has 1 App for now) ========================
  ml_prob.set_app_specs_pointer(& app_segment);
  
  
  // ======= Mesh - BEGIN  ==================
  MultiLevelMesh ml_mesh;
  // ======= Mesh - END  ==================

  // ======= Mesh, Coarse reading - BEGIN ==================
  double scalingFactor = 1.;
 
  const bool read_groups = true; //with this being false, we don't read any group at all. Therefore, we cannot even read the boundary groups that specify what are the boundary faces, for the boundary conditions
  const bool read_boundary_groups = true;
  
  std::string mesh_file_tot = "./input/" + app_segment._mesh_files[m];
  
  ml_mesh.ReadCoarseMesh(mesh_file_tot.c_str(), fe_quad_rule.c_str(), scalingFactor, read_groups, read_boundary_groups);
//     ml_mesh.GenerateCoarseBoxMesh(2,0,0,0.,1.,0.,0.,0.,0.,EDGE3,fe_quad_rule.c_str());
//     ml_mesh.GenerateCoarseBoxMesh(0,2,0,0.,0.,0.,1.,0.,0.,EDGE3,fe_quad_rule.c_str());
  // ======= Mesh, Coarse reading - END ==================

  // ======= Mesh: Refinement - BEGIN  ==================
  unsigned numberOfUniformLevels = /*1*/4;
  unsigned numberOfSelectiveLevels = 0;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  // ======= Mesh: Refinement - END  ==================

  // ======= Mesh: Coarse erasing - BEGIN  ========================
  ml_mesh.EraseCoarseLevels(numberOfUniformLevels + numberOfSelectiveLevels - 1);
  ml_mesh.PrintInfo();
  // ======= Mesh: Coarse erasing - END  ========================

  
  // ======= Solution - BEGIN ==================
  MultiLevelSolution ml_sol(&ml_mesh);

  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);

  // ======= Problem, Mesh and Solution  ==================
  ml_prob.SetMultiLevelMeshAndSolution(& ml_sol);
  // ======= Solution - END ==================

  // ======= Solutions that are Unknowns - BEGIN ==================
  // add variables to ml_sol
  ml_sol.AddSolution("d_s", LAGRANGE, FIRST/*DISCONTINUOUS_POLYNOMIAL, ZERO*/);

  // ======= Solution: Initial Conditions ==================
  ml_sol.Initialize("All");    // initialize all variables to zero
  ml_sol.Initialize("d_s", InitialValueDS, & ml_prob);

  // ======= Solution: Boundary Conditions ==================
  ml_sol.AttachSetBoundaryConditionFunction(app_segment._boundary_conditions_types_and_values);
  ml_sol.GenerateBdc("d_s", "Steady",  & ml_prob);

  // ======= Solutions that are Unknowns - END ==================

  
//   std::vector < std::vector < const elem_type_templ_base<double, double> *  > > elem_all = ml_prob.evaluate_all_fe<double, double>();
  
    // ======= Problem, System - BEGIN ========================
  ml_prob.clear_systems();

  NonLinearImplicitSystem& system = ml_prob.add_system < NonLinearImplicitSystem > (app_segment._system_name);
  
  system.SetDebugNonlinear(true);
 
  system.AddSolutionToSystemPDE("d_s");
 
  // attach the assembling function to system
  system.SetAssembleFunction( app_segment._assemble_function );

//   system.SetMaxNumberOfLinearIterations(2);
  // initialize and solve the system
  system.SetMgType(V_CYCLE/*F_CYCLE*//*M_CYCLE*/); //it doesn't matter if I use only 1 level


  system.SetOuterSolver(GMRES);
 
  system.init();
  
  system.MGsolve();
    // ======= Problem, System - END ========================

    // ======= Print - BEGIN ========================
  const std::string print_order = "biquadratic"; //"linear", "quadratic", "biquadratic"
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
 
  ml_sol.GetWriter()->Write(app_segment._mesh_files[m], files.GetOutputPath(), print_order.c_str(), variablesToBePrinted);
    // ======= Print - END ========================

  }  //end mesh file loop
 
 
  return 0;
}


