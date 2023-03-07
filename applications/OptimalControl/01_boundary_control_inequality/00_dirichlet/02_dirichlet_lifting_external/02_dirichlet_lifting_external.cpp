#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "NumericVector.hpp"


#include "../../../param.hpp"


using namespace femus;

#include  "../../../opt_systems_cost_functional.hpp"

#include  "../../../opt_systems_dirichlet.hpp"







double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

    if(!strcmp(name,"state")) {
        value = 0.;
    }
    else if(!strcmp(name,"control")) {
        value = 0.;
    }
    else if(!strcmp(name,"adjoint")) {
        value = 0.;
    }
    else if(!strcmp(name,"adjoint_ext")) {
        value = 0.;
    }
    else if(!strcmp(name,"mu")) {
        value = 0.;
    }
    else if(!strcmp(name,"TargReg")) {
        value = ctrl::cost_functional_without_regularization_Square_or_Cube::ElementTargetFlag(x);
    }
    else if(!strcmp(name,"ContReg")) {
        value = ctrl:: square_or_cube:: Domain_elements_containing_Gamma_control< ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ControlDomainFlag_external_restriction(x);
    }
    else if(!strcmp(name,"act_flag")) {
        value = 0.;
    }


    return value;
}


bool Solution_set_boundary_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

    bool dirichlet = false; //dirichlet
    value = 0.;

  if(!strcmp(name, "state")) {
       dirichlet = true;
       value = 0.;
  }
  else if(!strcmp(name, "control")) {
       dirichlet = true;
       value = 0.;
  }
  else if(!strcmp(name, "adjoint")) {
       dirichlet = true;
       value = 0.;
  }
  else if(!strcmp(name, "adjoint_ext")) {
       dirichlet = true;
       value = 0.;
  }
  else if(!strcmp(name, "mu")) {
        dirichlet = false;
    }
    
    
  else if(!strcmp(name,"act_flag")) {
        dirichlet = true;
    }
  else if(!strcmp(name,"TargReg")) {
    }
  else if(!strcmp(name,"ContReg")) {
        dirichlet = true;
    }
  

    return dirichlet;
}





int main(int argc, char** args) {

  // ======= Init ========================
    FemusInit mpinit(argc, args, MPI_COMM_WORLD);

    // ======= Problem ========================
    MultiLevelProblem ml_prob;

    // ======= Files ========================
  Files files; 
        files.CheckIODirectories(ctrl::use_output_time_folder);
        files.RedirectCout(ctrl::redirect_cout_to_file);
 
  // ======= Problem, Files ========================
  ml_prob.SetFilesHandler(&files);

  // ======= Problem, Quad Rule ========================
  std::string fe_quad_rule("seventh");
  ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule);
  ml_prob.set_all_abstract_fe_multiple();

    // ======= Mesh  ==================
    MultiLevelMesh ml_mesh;

    double scalingFactor = 1.;

    // read coarse level mesh and generate finers level meshes
    std::string mesh_file = "./input/extended_box_coarse.med";
//     std::string mesh_file = "./input/extended_box_coarse.med";
//     std::string mesh_file = "./input/extended_box.med";

    const double Lref = 1.;
    ml_mesh.ReadCoarseMesh(mesh_file.c_str(), fe_quad_rule.c_str(), Lref);

    //ml_mesh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,"seventh");
    unsigned numberOfUniformLevels = N_UNIFORM_LEVELS;
    const unsigned erased_levels = N_ERASED_LEVELS;
    unsigned numberOfSelectiveLevels = 0;
    ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
    ml_mesh.EraseCoarseLevels(erased_levels);
    ml_mesh.PrintInfo();

    // ======= Solution  ==================
    MultiLevelSolution ml_sol(&ml_mesh);

    ml_sol.SetWriter(VTK);
    ml_sol.GetWriter()->SetDebugOutput(true);
  
 // ======= Problem, Mesh and Solution  ==================
 ml_prob.SetMultiLevelMeshAndSolution(& ml_sol);


  // ======= Solutions that are Unknowns - BEGIN ==================
    ml_sol.AddSolution("state", LAGRANGE, FIRST);
    ml_sol.AddSolution("control", LAGRANGE, FIRST);
    ml_sol.AddSolution("adjoint", LAGRANGE, FIRST);
    ml_sol.AddSolution("adjoint_ext", LAGRANGE, FIRST);
    ml_sol.AddSolution("mu", LAGRANGE, FIRST);

  // ======= Solution: Initial Conditions ==================
    ml_sol.Initialize("state",       Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize("control",     Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize("adjoint",     Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize("adjoint_ext", Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize("mu",          Solution_set_initial_conditions, & ml_prob);

  // ======= Solution: Boundary Conditions ==================
    ml_sol.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
    ml_sol.GenerateBdc("state", "Steady", & ml_prob);
    ml_sol.GenerateBdc("control", "Steady", & ml_prob);
    ml_sol.GenerateBdc("adjoint", "Steady", & ml_prob);
    ml_sol.GenerateBdc("adjoint_ext", "Steady", & ml_prob);
    ml_sol.GenerateBdc("mu", "Steady", & ml_prob);
  // ======= Solutions that are Unknowns - END ==================

  // ======= Solutions that are not Unknowns - BEGIN  ==================
    ml_sol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
    ml_sol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
 
  //MU
  const bool      act_flag_is_an_unknown_of_a_pde = false;
  std::vector<std::string> act_set_flag_name(1);  act_set_flag_name[0] = "act_flag";
  const unsigned int act_set_fake_time_dep_flag = 2;
  ml_sol.AddSolution(act_set_flag_name[0].c_str(), LAGRANGE, /*FIRST*/SECOND, act_set_fake_time_dep_flag, act_flag_is_an_unknown_of_a_pde);
  //MU

    ml_sol.Initialize("TargReg",     Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize("ContReg",     Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize(act_set_flag_name[0].c_str(), Solution_set_initial_conditions, & ml_prob);

    ml_sol.GenerateBdc("TargReg", "Steady", & ml_prob);
    ml_sol.GenerateBdc("ContReg", "Steady", & ml_prob);
    ml_sol.GenerateBdc(act_set_flag_name[0].c_str(), "Steady", & ml_prob);
  // ======= Solutions that are not Unknowns - END  ==================

    
  //==== Solution: CHECK SOLUTION FE TYPES ==
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType("state")) abort();
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType("mu")) abort();
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType(act_set_flag_name[0].c_str())) abort();
  //==== Solution: CHECK SOLUTION FE TYPES ==

  

    
  // ======= Problem, System - BEGIN ========================
    NonLinearImplicitSystemWithPrimalDualActiveSetMethod & system = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("LiftRestr");

    system.SetActiveSetFlagName(act_set_flag_name);

    system.AddSolutionToSystemPDE("state");
    system.AddSolutionToSystemPDE("control");
    system.AddSolutionToSystemPDE("adjoint");
    system.AddSolutionToSystemPDE("adjoint_ext");
    system.AddSolutionToSystemPDE("mu");

    // attach the assembling function to system
    system.SetAssembleFunction( lifting_external::assemble_elliptic_dirichlet_control );

    system.SetDebugNonlinear(true);
//     system.SetDebugFunction( ctrl::cost_functional::compute_cost_functional_regularization_lifting_external );
    
    // system.SetMaxNumberOfNonLinearIterations(4);

    // initialize and solve the system
    system.init();
//     system.assemble_call_before_boundary_conditions(2);
    system.MGsolve();
  // ======= Problem, System  - END ========================

    // ======= Print ========================
    std::vector < std::string > variablesToBePrinted;
    variablesToBePrinted.push_back("all");
    ml_sol.GetWriter()->Write(files.GetOutputPath(), "biquadratic", variablesToBePrinted);

    return 0;

}



