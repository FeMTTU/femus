#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "NumericVector.hpp"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"



#include "../../../param.hpp"



using namespace femus;


#include  "../../../opt_systems.hpp"



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
        value = ctrl::cost_functional::cost_functional_Square_or_Cube::ElementTargetFlag(x);
    }
    else if(!strcmp(name,"ContReg")) {
        value = ctrl::Domain_elements_containing_Gamma_control< ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ControlDomainFlag_external_restriction(x);
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





void compute_cost_functional_regularization_lifting_external(const MultiLevelProblem& ml_prob, 
                     const unsigned level, 
                     const unsigned iteration,
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars  
);




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
    system.SetAssembleFunction(AssembleLiftExternalProblem);

    system.SetDebugNonlinear(true);
    system.SetDebugFunction(compute_cost_functional_regularization_lifting_external);
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



void compute_cost_functional_regularization_lifting_external(const MultiLevelProblem& ml_prob, 
                     const unsigned level, 
                     const unsigned iteration,
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars  
)    {



    Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);            // pointer to the mesh (level) object

    MultiLevelSolution*    ml_sol = ml_prob._ml_sol;                             // pointer to the multilevel solution object
    Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

    const unsigned           dim = msh->GetDimension();                         // get the domain dimension of the problem
    const unsigned          dim2 = (3 * (dim - 1) + !(dim - 1));                // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
    const unsigned       max_size = static_cast< unsigned >(ceil(pow(3, dim)));  // conservative: based on line3, quad9, hex27

    const unsigned         iproc = msh->processor_id();                         // get the process_id (for parallel computation)

//***************************************************
    vector < vector < double > > coords_at_dofs(dim);   // local coordinates
    unsigned solType_coords  = BIQUADR_FE;  // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
    for (unsigned i = 0; i < dim; i++) {
        coords_at_dofs[i].reserve(max_size);
    }
//***************************************************

//***************************************************
    double weight_qp; // gauss point weight

//***************************************************
    double alpha = ALPHA_CTRL_VOL;
    double beta  = BETA_CTRL_VOL;

//******************** state ************************
//***************************************************
    vector <double> phi_u;    // local test function
    vector <double> phi_u_x;  // local test function first order partial derivatives
    vector <double> phi_u_xx; // local test function second order partial derivatives

    phi_u.reserve(max_size);
    phi_u_x.reserve(max_size * dim);
    phi_u_xx.reserve(max_size * dim2);


    unsigned solIndex_u;
    solIndex_u = ml_sol->GetIndex("state");                    // get the position of "state" in the ml_sol object
    unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);  // get the finite element type for "state"

    vector < double >  sol_u; // local solution
    sol_u.reserve(max_size);

    double u_gss = 0.;
//***************************************************
//***************************************************

//******************** control **********************
//***************************************************
    vector <double> phi_ctrl;    // local test function
    vector <double> phi_ctrl_x;  // local test function first order partial derivatives
    vector <double> phi_ctrl_xx; // local test function second order partial derivatives

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
    vector <double> phi_udes;    // local test function
    vector <double> phi_udes_x;  // local test function first order partial derivatives
    vector <double> phi_udes_xx; // local test function second order partial derivatives

    phi_udes.reserve(max_size);
    phi_udes_x.reserve(max_size * dim);
    phi_udes_xx.reserve(max_size * dim2);

    vector < double >  sol_udes; // local solution
    sol_udes.reserve(max_size);

    double udes_gss = 0.;
//***************************************************
//***************************************************

//********************* DATA ************************
    double u_des = ctrl::cost_functional::cost_functional_Square_or_Cube::DesiredTarget();
//***************************************************

    double integral_target = 0.;
    double integral_alpha  = 0.;
    double integral_beta   = 0.;


    // element loop: each process loops only on the elements that owns
    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        int group_flag         = msh->GetElementGroup(iel);      // element group flag (Exterior = GROUP_EXTERNAL, Interior = GROUP_INTERNAL)
        short unsigned ielGeom = msh->GetElementType(iel);       // element geometry type

//***************** GEOMETRY ************************
        unsigned nDofx = msh->GetElementDofNumber(iel, solType_coords); // number of coordinate element dofs
        for (int i = 0; i < dim; i++)  coords_at_dofs[i].resize(nDofx);
        // local storage of coordinates
        for (unsigned i = 0; i < nDofx; i++) {
            unsigned xDof  = msh->GetSolutionDof(i, iel, solType_coords); // global to global mapping between coordinates node and coordinate dof

            for (unsigned jdim = 0; jdim < dim; jdim++) {
                coords_at_dofs[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);  // global extraction and local storage for the element coordinates
            }
        }

        // elem average point
        vector < double > elem_center(dim);
        for (unsigned j = 0; j < dim; j++) {
            elem_center[j] = 0.;
        }
        for (unsigned j = 0; j < dim; j++) {
            for (unsigned i = 0; i < nDofx; i++) {
                elem_center[j] += coords_at_dofs[j][i];
            }
        }

        for (unsigned j = 0; j < dim; j++) {
            elem_center[j] = elem_center[j]/nDofx;
        }
//***************************************************

//************** set target domain flag *************
        int target_flag = 0;
        target_flag = ctrl::cost_functional::cost_functional_Square_or_Cube::ElementTargetFlag(elem_center);
//***************************************************


//**************** state ****************************
        unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);     // number of solution element dofs
        sol_u    .resize(nDof_u);
        // local storage of global mapping and solution
        for (unsigned i = 0; i < sol_u.size(); i++) {
            unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);        // global to global mapping between solution node and solution dof
            sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);              // global extraction and local storage for the solution
        }
//***************************************************


//************** control ****************************
        unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);  // number of solution element dofs
        sol_ctrl    .resize(nDof_ctrl);
        for (unsigned i = 0; i < sol_ctrl.size(); i++) {
            unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl); // global to global mapping between solution node and solution dof
            sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);           // global extraction and local storage for the solution
        }
//***************************************************


//**************** u_des ****************************
        unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
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
        for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {

            // *** get gauss point weight, test function and test function partial derivatives ***
            //  ==== State
            msh->_finiteElement[ielGeom][solType_u]   ->Jacobian(coords_at_dofs, ig, weight_qp, phi_u, phi_u_x, phi_u_xx);
            //  ==== Adjoint
            msh->_finiteElement[ielGeom][solType_u/*solTypeTdes*/]->Jacobian(coords_at_dofs, ig, weight_qp, phi_udes, phi_udes_x, phi_udes_xx);
            //  ==== Control
            msh->_finiteElement[ielGeom][solType_ctrl]  ->Jacobian(coords_at_dofs, ig, weight_qp, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);

            u_gss = 0.;
            for (unsigned i = 0; i < nDof_u; i++) u_gss += sol_u[i] * phi_u[i];
            ctrl_gss = 0.;
            for (unsigned i = 0; i < nDof_ctrl; i++) ctrl_gss += sol_ctrl[i] * phi_ctrl[i];
            udes_gss  = 0.;
            for (unsigned i = 0; i < nDof_udes; i++)  udes_gss  += sol_udes[i]  * phi_udes[i];
            ctrl_x_gss  = 0.;
            for (unsigned i = 0; i < nDof_ctrl; i++)
            {
                for (unsigned idim = 0; idim < dim; idim ++) ctrl_x_gss  += sol_ctrl[i] * phi_ctrl_x[i + idim * nDof_ctrl];
            }

            integral_target += target_flag * weight_qp * (u_gss - udes_gss) * (u_gss - udes_gss);
            integral_alpha  += (group_flag - GROUP_INTERNAL) * weight_qp * ctrl_gss * ctrl_gss;
            integral_beta   += (group_flag - GROUP_INTERNAL) * weight_qp * ctrl_x_gss * ctrl_x_gss;

        } // end gauss point loop
    } //end element loop
    
//     std::ios_base::fmtflags f( std::cout.flags() );

    std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << integral_target << std::endl;
    std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << integral_alpha << std::endl;
    std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << integral_beta << std::endl;
    std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_target + 0.5 * alpha * integral_alpha + 0.5 * beta * integral_beta << std::endl;

//     std::cout.flags( f );  ///@todo attempt at restoring default cout flags

    return;

}



