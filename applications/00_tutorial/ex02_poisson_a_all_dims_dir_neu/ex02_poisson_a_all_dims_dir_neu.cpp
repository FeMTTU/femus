/**
 * Solve 
 *     - \Delta u = 1
*/


#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "NumericVector.hpp"

#include "CurrentElem.hpp"
#include "ElemType_template.hpp"

#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"


#include "00_poisson_eqn_with_dirichlet_or_neumann_bc.hpp"


using namespace femus;
 

// SEGMENT - BEGIN

namespace segment {

  
  namespace function_0 {
    

double value(const std::vector<double> & x) {
    
    // for a 1d segment
    
    return  x[0] * (1. - x[0]);
}


// user-made equation - accepts only coordinates
double laplacian(const std::vector<double> & x){
    
    // for a 1d segment
    
    return  -2.;
     }

  }

  
  
 
// This depends on: 
// mesh file (for the face flags)
// underlying exact solution (if provided)
// equation (name of the unknowns; also, if it is Laplace, biharmonic, Stokes, etc)
// 
// 
bool bc_all_dirichlet(const MultiLevelProblem * ml_prob,
                      const std::vector < double >& x,
                      const char name[], 
                      double& value,
                      const int face_name,
                      const double time) {

  bool dirichlet = false;
  value = 0.;
  
  
 if (ml_prob->GetMLMesh()->GetDimension() != 1 )  abort();
  
  if (face_name == 1) {
      dirichlet = true;
        value = 0.;
    }
  else if (face_name == 2) {
      dirichlet = true;
        value = 0.;
    }

 
  return dirichlet;
  
 }
 
  }
  
  

// SEGMENT - END




double InitialValueDS(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {
    
  return 0.;
  
}


 


int main(int argc, char** args) {
 // Very impotant to be there to run in parallel
  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Problem ========================
  MultiLevelProblem ml_prob;

  // ======= Files - BEGIN  ========================
  const bool use_output_time_folder = false; //not everytime this will go to output folder. This allows you to run the code multiple times without overwritting. This guy will generate output folder each time you run.
  const bool redirect_cout_to_file = false; // puts the output to a file instead of the command window or the konsole. In a log file.
  Files files; 
        files.CheckIODirectories(use_output_time_folder); //dot operator is calling or accessing the function in that class of that name.
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
  const std::string relative_path_to_build_directory =  "../../../";
  app_segment._mesh_files.push_back("segment_16_dir_neu.med");
  app_segment._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/01_1d/segment_0-1/");
  
  app_segment._system_name = "Equation";
  app_segment._assemble_function = femus::poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;
  
  app_segment._boundary_conditions_types_and_values             = segment::bc_all_dirichlet;
  app_segment._assemble_function_rhs = segment::function_0::laplacian;
//   app_segment._true_solution    = segment::function_0::value;  
  ///@todo if this is not set, nothing happens here. It is used to compute absolute errors

  
    // ======= App Specifics - END  ==================


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
  
  const std::string mesh_file =  app_segment._mesh_files_path_relative_to_executable[m] + 
                                 app_segment._mesh_files[m];
    
  ml_mesh.ReadCoarseMesh(mesh_file.c_str(), fe_quad_rule.c_str(), scalingFactor, read_groups, read_boundary_groups);
//     ml_mesh.GenerateCoarseBoxMesh(2,0,0,0.,1.,0.,0.,0.,0.,EDGE3,fe_quad_rule.c_str());
//     ml_mesh.GenerateCoarseBoxMesh(0,2,0,0.,0.,0.,1.,0.,0.,EDGE3,fe_quad_rule.c_str());
  // ======= Mesh, Coarse reading - END ==================

  // ======= Mesh: Refinement - BEGIN  ==================
  unsigned numberOfUniformLevels = /*1*/3;
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
  ml_prob.SetMultiLevelMeshAndSolution(& ml_sol); // & looking for the address of the ml_sol
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


