#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "NumericVector.hpp"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"

#include "ElemType.hpp"


#define FACE_FOR_CONTROL             1

#include "../../param.hpp"


#define FE_DOMAIN  2 //with 0 it only works in serial, you must put 2 to make it work in parallel...: that's because when you fetch the dofs from _topology you get the wrong indices

///@todo do a very weak impl of Laplacian
///@todo Review the ordering for phi_ctrl_x_bdry
///@todo check computation of 2nd derivatives in elem_type_template
///@todo Implement rather fast way to add inequality constraint to a problem
///@todo If I have a mesh that has 1 element only at the coarse level, can I run it in parallel? I need to factorize the ReadCoarseMesh function
///@todo merge elliptic_nonlin in here
///@todo What if I did a Point domain, could I solve ODEs in time like this? :)
///@todo Re-double check that things are fine in elem_type_template, probably remove _gauss_bdry!
///@todo See if with Petsc you can enforce Dirichlet conditions using NEGATIVE indices
///@todo Do Parallel ComputeIntegral
///@todo Remove the prints, possible cause of slowing down (maybe do assert)
///@todo The \mu/actflag pieces are now basically separated, except for setting to zero on Omega minus Gamma_c (such as is done for control)
///@todo put assembleMatrix everywhere there is a filling of the matrix!
///@todo I tried the assembly alone and it seems fine (well, actually there are problems in parallel...). The problem seems to be in the Solve part... am I doing something weird there?
///@todo the compare function is probably responsible for parallel slowing down!


using namespace femus;



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
    else if(!strcmp(name,"mu")) {
        value = 0.;
    }
    else if(!strcmp(name,"TargReg")) {
        value = ElementTargetFlag(x);
    }
    else if(!strcmp(name,"ContReg")) {
        value = ControlDomainFlag_bdry(x);
    }
    else if(!strcmp(name,"act_flag")) {
        value = 0.;
    }


    return value;
}



///@todo notice that even if you set Dirichlet from the mesh file, here you can override it
bool Solution_set_boundary_conditions(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet; // = true; //dirichlet
  value = 0.;

  if(!strcmp(name,"control")) {
      
  if (faceName == FACE_FOR_CONTROL) {
     if (x[ axis_direction_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 && x[ axis_direction_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5)
     { dirichlet = false; }
     else { dirichlet = true;  }
  }
  else { dirichlet = true;  }
  
  }

  else if(!strcmp(name,"state")) {  //"state" corresponds to the first block row (u = q)
      
  if (faceName == FACE_FOR_CONTROL) {
      
     if (x[ axis_direction_Gamma_control(faceName) ] > GAMMA_CONTROL_LOWER - 1.e-5 && x[ axis_direction_Gamma_control(faceName) ] < GAMMA_CONTROL_UPPER + 1.e-5) 
     { dirichlet = false; }
     else { dirichlet = true;  }
  }
  else { dirichlet = true;  }
      
  }

  else if(!strcmp(name,"mu")) {
      
    dirichlet = false;

  }
  
  else { dirichlet = true; }
  
//     if(!strcmp(name,"adjoint")) { 
//     dirichlet = false;
//   }

  
  return dirichlet;
}


void ComputeIntegral(const MultiLevelProblem& ml_prob);

void AssembleOptSys(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {

  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();

  // ======= Quad Rule ========================
  std::string fe_quad_rule("seventh");

  // ======= Mesh  ==================
  MultiLevelMesh ml_mesh;

  
  std::string input_file = "square_4x5.med";
//   std::string input_file = "square_parametric.med";
//   std::string input_file = "Mesh_3_groups.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
  const std::string infile = mystream.str();
  const double Lref = 1.;
  
  
  ml_mesh.ReadCoarseMesh(infile.c_str(), fe_quad_rule.c_str(), Lref);
  
//   ml_mesh.GenerateCoarseBoxMesh(NSUB_X, NSUB_Y, 0, 0., 1., 0., 1., 0., 0., QUAD9, fe_quad_rule.c_str());  
//   ml_mesh.GenerateCoarseBoxMesh(NSUB_X, NSUB_Y, NSUB_Z, 0., 1., 0., 1., 0., 1., HEX27, fe_quad_rule.c_str());  
     ///@todo seems like GenerateCoarseBoxMesh doesn't assign flags to faces correctly, 
     //so I created a .med file with the following flags:
     //1: x = x_min  //2: x = x_max  //3: y = y_min  //4: y = y_max //5: z = z_min //6: z = z_max
     //1: x = x_min  //2: x = x_max  //3: y = y_min  //4: y = y_max                                (in 2d) in the new .med file
  
   //1: bottom  //2: right  //3: top  //4: left (in 2d) GenerateCoarseBoxMesh 
  

  unsigned numberOfUniformLevels = 4;
  unsigned numberOfSelectiveLevels = 0;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);
  ml_mesh.PrintInfo();

  // ======= Solution  ==================
  MultiLevelSolution ml_sol(&ml_mesh);

  ml_sol.AddSolution("state",   LAGRANGE, SECOND/*FIRST*/);
  ml_sol.AddSolution("control", LAGRANGE, SECOND/*FIRST*/);
  ml_sol.AddSolution("adjoint", LAGRANGE, SECOND/*FIRST*/);
  ml_sol.AddSolution("mu",      LAGRANGE, SECOND/*FIRST*/);  //MU
  ml_sol.AddSolution("TargReg", DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  ml_sol.AddSolution("ContReg", DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  //MU
  const unsigned int fake_time_dep_flag = 2;
  const std::string act_set_flag_name = "act_flag";
  ml_sol.AddSolution(act_set_flag_name.c_str(), LAGRANGE, SECOND/*FIRST*/, fake_time_dep_flag);               //this variable is not solution of any eqn, it's just a given field
  //MU

  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType("state")) abort();
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType("mu")) abort();
  if ( ml_sol.GetSolutionType("control") != ml_sol.GetSolutionType(act_set_flag_name.c_str())) abort();
      
  // ======= Problem  ==================
  MultiLevelProblem ml_prob(&ml_sol);
  
  ml_prob.SetFilesHandler(&files);
  ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule);
  ml_prob.set_all_abstract_fe();

  // ======= Solution: Initial Conditions ==================
  ml_sol.Initialize("All");    // initialize all varaibles to zero

  ml_sol.Initialize("state",       Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("control",     Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("adjoint",     Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("mu",          Solution_set_initial_conditions, & ml_prob);   //MU
  ml_sol.Initialize("TargReg",     Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize("ContReg",     Solution_set_initial_conditions, & ml_prob);
  ml_sol.Initialize(act_set_flag_name.c_str(), Solution_set_initial_conditions, & ml_prob);   //MU


  // ======= Solution: Boundary Conditions ==================
  ml_sol.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
  ml_sol.GenerateBdc("state");
  ml_sol.GenerateBdc("control");
  ml_sol.GenerateBdc("adjoint");
  ml_sol.GenerateBdc("mu");   //MU

  // ======= System ========================
  NonLinearImplicitSystemWithPrimalDualActiveSetMethod& system = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("BoundaryControl");
  
  system.SetActiveSetFlagName(act_set_flag_name);    //MU
//   system.SetMaxNumberOfNonLinearIterations(50);

  system.AddSolutionToSystemPDE("state");  
  system.AddSolutionToSystemPDE("control");  
  system.AddSolutionToSystemPDE("adjoint");  
  system.AddSolutionToSystemPDE("mu");     //MU
  
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleOptSys);
  
  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);

  system.SetDebugNonlinear(true);
  system.SetDebugFunction(ComputeIntegral);  //weird error if I comment this line, I expect nothing to happen but something in the assembly gets screwed up in memory I guess
   
//   // initialize and solve the system
  system.init();
  system.MGsolve();
//   system.assemble_call(1);
  
  // ======= Print ========================
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  ml_sol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/, "biquadratic", variablesToBePrinted);

  return 0;
}




  
  

//This Opt system is characterized by the following ways of setting matrix values:
// Add_values (Mat or Vec) in the volume loop
// Add_values (Mat or Vec) in the boundary loop
// Insert_values (Mat or Vec) in the boundary loop
// Insert_values (Mat or Vec) in the volume loop
// Insert_values (Mat or Vec) outside all loops
// We're going to split the two parts and add a close() at the end of each


void AssembleOptSys(MultiLevelProblem& ml_prob) {
    
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("BoundaryControl");
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;
  
  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  constexpr bool print_algebra_global = false;
  constexpr bool print_algebra_local = false;
  
  
  //=============== Integration ========================================
  double weight = 0.;
  double weight_bdry = 0.;


  //=============== Geometry ========================================
   unsigned solType_coords = FE_DOMAIN;
 
  CurrentElem < double > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  
  std::vector< double > normal(space_dim, 0.);
 //***************************************************  
  
  vector < double > coord_at_qp_bdry(space_dim);
  
  vector <double> phi_coords;
  vector <double> phi_coords_x;
  vector <double> phi_coords_xx; 

  phi_coords.reserve(max_size);
  phi_coords_x.reserve(max_size * space_dim);
  phi_coords_xx.reserve(max_size * dim2);
  

 //*************************************************** 

 //********************* state *********************** 
 //*************************************************** 
  vector <double> phi_u;
  vector <double> phi_u_x;
  vector <double> phi_u_xx;

  phi_u.reserve(max_size);
  phi_u_x.reserve(max_size * space_dim);
  phi_u_xx.reserve(max_size * dim2);
  
  
  //boundary state shape functions
  vector <double> phi_u_bdry;  
  vector <double> phi_u_x_bdry; 

  phi_u_bdry.reserve(max_size);
  phi_u_x_bdry.reserve(max_size * space_dim);
  
 //***************************************************  
 //***************************************************  

  
 //********************** adjoint ********************
 //*************************************************** 
  vector <double> phi_adj;
  vector <double> phi_adj_x;
  vector <double> phi_adj_xx;

  phi_adj.reserve(max_size);
  phi_adj_x.reserve(max_size * space_dim);
  phi_adj_xx.reserve(max_size * dim2);
 

  //boundary adjoint shape functions  
  vector <double> phi_adj_bdry;  
  vector <double> phi_adj_x_bdry; 

  phi_adj_bdry.reserve(max_size);
  phi_adj_x_bdry.reserve(max_size * space_dim);

  
  //volume shape functions at boundary
  vector <double> phi_adj_vol_at_bdry;
  vector <double> phi_adj_x_vol_at_bdry;
  phi_adj_vol_at_bdry.reserve(max_size);
  phi_adj_x_vol_at_bdry.reserve(max_size * space_dim);
  vector <double> sol_adj_x_vol_at_bdry_gss(space_dim);
 //*************************************************** 
 //*************************************************** 

  
 //********************* bdry cont *******************
 //*************************************************** 
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 
  vector <double> phi_xx_bdry_placeholder;

  phi_ctrl_bdry.reserve(max_size);
  phi_ctrl_x_bdry.reserve(max_size * space_dim);
  phi_xx_bdry_placeholder.reserve(max_size * dim2);
 //*************************************************** 

  //MU
  //************** act flag ****************************   
  unsigned int solIndex_act_flag_sol; 
  unsigned int solFEType_act_flag_sol;
  store_act_flag_in_old(mlPdeSys, ml_sol, sol,
                        solIndex_act_flag_sol, //this becomes a vector
                        solFEType_act_flag_sol //remove this one, only Index
                       );
  
  
  //********* variables for ineq constraints *****************
  vector < double/*int*/ >  sol_actflag;   sol_actflag.reserve(max_size); //flag for active set
  vector < double >  ctrl_lower;   ctrl_lower.reserve(max_size);
  vector < double >  ctrl_upper;   ctrl_upper.reserve(max_size);
  const int ineq_flag = INEQ_FLAG;
  const double c_compl = C_COMPL;
  //***************************************************  
  //MU
  
  
//***************************************************
//********* WHOLE SET OF VARIABLES ******************
    const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();
    const unsigned int n_quantities = ml_sol->GetSolutionSize();

//***************************************************
    enum Pos_in_matrix {pos_mat_state = 0, pos_mat_ctrl, pos_mat_adj, pos_mat_mu}; //these are known at compile-time 
                    ///@todo these are the positions in the MlSol object or in the Matrix? I'd say the matrix, but we have to check where we use it...

                    
    assert(pos_mat_state   == mlPdeSys->GetSolPdeIndex("state"));
    assert(pos_mat_ctrl    == mlPdeSys->GetSolPdeIndex("control"));
    assert(pos_mat_adj     == mlPdeSys->GetSolPdeIndex("adjoint"));
    assert(pos_mat_mu      == mlPdeSys->GetSolPdeIndex("mu"));
//***************************************************

    vector < std::string > Solname_Mat(n_unknowns);  //this coincides with Pos_in_matrix
    Solname_Mat[0] = "state";
    Solname_Mat[1] = "control";
    Solname_Mat[2] = "adjoint";
    Solname_Mat[3] = "mu";

 //***************************************************
   enum Pos_in_Sol {pos_sol_state = 0, pos_sol_ctrl, pos_sol_adj, pos_sol_mu, pos_sol_targreg, pos_sol_contreg, pos_sol_actflag}; //these are known at compile-time 

        assert(pos_sol_actflag == solIndex_act_flag_sol);
//***************************************************
    
    vector < std::string > Solname_quantities(n_quantities);
    
        for(unsigned ivar=0; ivar < Solname_quantities.size(); ivar++) {
            Solname_quantities[ivar] = ml_sol->GetSolutionName(ivar);
        }
        
    vector < unsigned > SolIndex_Mat(n_unknowns);      //should have Sol order
    vector < unsigned > SolFEType_Mat(n_unknowns);       //should have Mat order
    vector < unsigned > SolPdeIndex(n_unknowns);     //should have Mat order, of course

    vector < unsigned > SolIndex_quantities(n_quantities);      //should have Sol order
    vector < unsigned > SolFEType_quantities(n_quantities);     //should have Sol order
    vector < unsigned > Sol_n_el_dofs_quantities(n_quantities); //should have Sol order
 
  

    for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
        SolIndex_Mat[ivar]    = ml_sol->GetIndex        (Solname_Mat[ivar].c_str());
        SolFEType_Mat[ivar]   = ml_sol->GetSolutionType(SolIndex_Mat[ivar]);
        SolPdeIndex[ivar] = mlPdeSys->GetSolPdeIndex(Solname_Mat[ivar].c_str());
    }
    
    for(unsigned ivar=0; ivar < n_quantities; ivar++) {
        SolIndex_quantities[ivar]    = ml_sol->GetIndex        (Solname_quantities[ivar].c_str());
        SolFEType_quantities[ivar]   = ml_sol->GetSolutionType(SolIndex_quantities[ivar]);
    }    

    vector < unsigned > Sol_n_el_dofs_Mat(n_unknowns);    //should have Mat order

//***************************************************
    //----------- quantities (at dof objects) ------------------------------
    std::vector< int >       L2G_dofmap_Mat_AllVars;
    L2G_dofmap_Mat_AllVars.reserve( n_unknowns * max_size );
    vector < vector < int > >     L2G_dofmap_Mat(n_unknowns);     //should have Mat order
    for(int i = 0; i < n_unknowns; i++) {
        L2G_dofmap_Mat[i].reserve(max_size);
    }
    
    vector < vector < double > >  sol_eldofs_Mat(n_unknowns);  //should have Mat order
    for(int k = 0; k < n_unknowns; k++) {        sol_eldofs_Mat[k].reserve(max_size);    }


 //*************************************************** 
  std::vector< double > Res;   Res.reserve( n_unknowns*max_size);                         //should have Mat order
  std::vector < double > Jac;  Jac.reserve( n_unknowns*max_size * n_unknowns*max_size);   //should have Mat order
 //*************************************************** 

 
 //********************* DATA ************************ 
  double u_des = DesiredTarget();
  double alpha = ALPHA_CTRL_BDRY;
  double beta  = BETA_CTRL_BDRY;
  double penalty_outside_control_boundary = 1.e50;       // penalty for zero control outside Gamma_c and zero mu outside Gamma_c
  double penalty_strong_bdry = 1.e20;  // penalty for boundary equation on Gamma_c
  double penalty_ctrl = 1.e10;         //penalty for u=q
 //*************************************************** 
  
  
  if (assembleMatrix)  KK->zero();

 //*************************************************** 
     std::vector < std::vector < double > >  JacI_qp(space_dim);
     std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < Jac_qp.size(); d++) {   Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp.size(); d++) { JacI_qp[d].resize(dim); }
    
    double detJac_qp;

     std::vector < std::vector < double > >  JacI_qp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_qp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_qp_bdry.size(); d++) {   Jac_qp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp_bdry.size(); d++) { JacI_qp_bdry[d].resize(dim-1); }
    
    double detJac_qp_bdry;
    
    //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
//*************************************************** 
 
    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    geom_element.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element.geom_type();


   el_dofs_unknowns(sol, msh, pdeSys, iel,
                        SolFEType_Mat,
                        SolIndex_Mat,
                        SolPdeIndex,
                        Sol_n_el_dofs_Mat, 
                        sol_eldofs_Mat,  
                        L2G_dofmap_Mat);
  
        
 
 //***************************************************
    unsigned int nDof_max          = ElementJacRes<double>::compute_max_n_dofs(Sol_n_el_dofs_Mat);
    
    unsigned int sum_Sol_n_el_dofs = ElementJacRes<double>::compute_sum_n_dofs(Sol_n_el_dofs_Mat);

    Res.resize(sum_Sol_n_el_dofs);
    std::fill(Res.begin(), Res.end(), 0.);

    Jac.resize(sum_Sol_n_el_dofs * sum_Sol_n_el_dofs);
    std::fill(Jac.begin(), Jac.end(), 0.);
    
    L2G_dofmap_Mat_AllVars.resize(0);
      for (unsigned  k = 0; k < n_unknowns; k++)     L2G_dofmap_Mat_AllVars.insert(L2G_dofmap_Mat_AllVars.end(), L2G_dofmap_Mat[k].begin(), L2G_dofmap_Mat[k].end());
 //***************************************************

      
      
  //************* set target domain flag **************
   geom_element.set_elem_center(iel, solType_coords);

   int target_flag = 0;
   target_flag = ElementTargetFlag(geom_element.get_elem_center());
 //*************************************************** 
   

 //************ set control flag *********************
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag_bdry(geom_element.get_elem_center());
  std::vector<int> control_node_flag(Sol_n_el_dofs_Mat[pos_mat_ctrl],0);
 //*************************************************** 
 
  
  
 //===================================================   

	// Perform face loop over elements that contain some control face
	if (control_el_flag == 1) {
	  
	  double tau=0.;
	  std::vector<double> normal(space_dim,0.);
	       
	  // loop on faces of the current element

	  for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, jface);    
       
       std::vector<unsigned int> Sol_el_n_dofs_current_face(n_quantities); ///@todo the active flag is not an unknown!

       for (unsigned  k = 0; k < Sol_el_n_dofs_current_face.size(); k++) {
                 if (SolFEType_quantities[k] < 3) Sol_el_n_dofs_current_face[k] = msh->GetElementFaceDofNumber(iel, jface, SolFEType_quantities[k]);  ///@todo fix this absence
       }
       
       const unsigned nDof_max_bdry = ElementJacRes<double>::compute_max_n_dofs(Sol_el_n_dofs_current_face);
       
       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element.set_elem_center_bdry_3d();

	    // look for boundary faces
            const int bdry_index = el->GetFaceElementIndex(iel, jface);
            
	    if( bdry_index < 0) {
	      const unsigned int face_in_rectangle_domain = -( msh->el->GetFaceElementIndex(iel,jface)+1);
		
// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
	      if(  face_in_rectangle_domain == FACE_FOR_CONTROL) { //control face
              
 //=================================================== 
		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
	      bool  dir_bool = ml_sol->GetBdcFunction()(geom_element.get_elem_center_bdry(), Solname_Mat[pos_mat_ctrl].c_str(), tau, face_in_rectangle_domain, 0.);

 //=================================================== 
        
 
//========= initialize gauss quantities on the boundary ============================================
                double sol_ctrl_bdry_gss = 0.;
                double sol_adj_bdry_gss = 0.;
                std::vector<double> sol_ctrl_x_bdry_gss(space_dim);   std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);

//========= initialize gauss quantities on the boundary ============================================
		
        const unsigned n_gauss_bdry = ml_prob.GetQuadratureRule(ielGeom_bdry).GetGaussPointsNumber();
        
    
		for(unsigned ig_bdry = 0; ig_bdry < n_gauss_bdry; ig_bdry++) {
    
    elem_all[ielGeom_bdry][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_qp_bdry, JacI_qp_bdry, detJac_qp_bdry, space_dim);
	elem_all[ielGeom_bdry][solType_coords]->compute_normal(Jac_qp_bdry, normal);
    
    weight_bdry = detJac_qp_bdry * ml_prob.GetQuadratureRule(ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];

    elem_all[ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_ctrl_bdry, phi_ctrl_x_bdry, phi_xx_bdry_placeholder, space_dim);
    elem_all[ielGeom_bdry][SolFEType_quantities[pos_sol_state]]->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_u_bdry, phi_u_x_bdry,  phi_xx_bdry_placeholder, space_dim);
    elem_all[ielGeom_bdry][SolFEType_quantities[pos_sol_adj]]  ->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_adj_bdry, phi_adj_x_bdry,  phi_xx_bdry_placeholder, space_dim);


    elem_all[ielGeom][solType_coords]->JacJacInv_vol_at_bdry_new(geom_element.get_coords_at_dofs_3d(), ig_bdry, jface, Jac_qp/*not_needed_here*/, JacI_qp, detJac_qp/*not_needed_here*/, space_dim);
    elem_all[ielGeom][SolFEType_quantities[pos_sol_adj]]->shape_funcs_vol_at_bdry_current_elem(ig_bdry, jface, JacI_qp, phi_adj_vol_at_bdry, phi_adj_x_vol_at_bdry, boost::none, space_dim);
     
//     msh->_finiteElement[ielGeom][SolFEType_quantities[pos_sol_adj]]->fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem(geom_element.get_coords_at_dofs(), geom_element.get_coords_at_dofs_bdry_3d(), jface, ig_bdry, phi_adj_vol_at_bdry, phi_adj_x_vol_at_bdry);

		  
//========== compute gauss quantities on the boundary ===============================================
		  sol_ctrl_bdry_gss = 0.;
                  std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);
		      for (int i_bdry = 0; i_bdry < Sol_n_el_dofs_quantities[pos_sol_ctrl]; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
			
			sol_ctrl_bdry_gss +=  sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_bdry[i_bdry];
                            for (int d = 0; d < space_dim; d++) {
			      sol_ctrl_x_bdry_gss[d] += sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_x_bdry[i_bdry * space_dim + d];
			    }
		      }
		      
		      
		  sol_adj_bdry_gss = 0.;
		      for (int i_bdry = 0; i_bdry < Sol_n_el_dofs_quantities[pos_sol_adj]; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
			
			sol_adj_bdry_gss  +=  sol_eldofs_Mat[pos_mat_adj][i_vol] * phi_adj_bdry[i_bdry];
              }		      
		      
//=============== grad dot n for residual ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
           std::fill(sol_adj_x_vol_at_bdry_gss.begin(), sol_adj_x_vol_at_bdry_gss.end(), 0.);
		      for (int iv = 0; iv < Sol_n_el_dofs_Mat[pos_mat_adj]; iv++)  {
			
         for (int d = 0; d < space_dim; d++) {
			      sol_adj_x_vol_at_bdry_gss[d] += sol_eldofs_Mat[pos_mat_adj][iv] * phi_adj_x_vol_at_bdry[iv * space_dim + d];
			    }
		      }  
		      
    double grad_adj_dot_n_res = 0.;
        for(unsigned d=0; d < space_dim; d++) {
	  grad_adj_dot_n_res += sol_adj_x_vol_at_bdry_gss[d] * normal[d];  
	}
//=============== grad dot n  for residual =========================================       

//========== compute gauss quantities on the boundary ================================================

		  // *** phi_i loop ***
		  for(unsigned i_bdry=0; i_bdry < nDof_max_bdry; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

                 double lap_rhs_dctrl_ctrl_bdry_gss_i = 0.;
                 for (unsigned d = 0; d < space_dim; d++) {
                       if ( i_vol < Sol_n_el_dofs_Mat[pos_mat_ctrl] )  lap_rhs_dctrl_ctrl_bdry_gss_i +=  phi_ctrl_x_bdry[i_bdry * space_dim + d] * sol_ctrl_x_bdry_gss[d];
                 }
                 
//=============== construct control node flag field on the go  =========================================    
	      /* (control_node_flag[i])       picks nodes on \Gamma_c
	         (1 - control_node_flag[i])   picks nodes on \Omega \setminus \Gamma_c
	       */
	      if (dir_bool == false) { 
// 		std::cout << " found boundary control nodes ==== " << std::endl;
			for(unsigned k=0; k<control_node_flag.size(); k++) {
				  control_node_flag[i_vol] = 1;
			}
              }
//=============== construct control node flag field on the go  =========================================    

		 
//============ Bdry Residuals ==================	
                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat,pos_mat_state,i_vol) ] +=  - control_node_flag[i_vol] * penalty_ctrl * (   sol_eldofs_Mat[pos_mat_state][i_vol] - sol_eldofs_Mat[pos_mat_ctrl][i_vol] )
                    - control_node_flag[i_vol] *  weight_bdry * (grad_adj_dot_n_res * phi_u_bdry[i_bdry]);   // u = q


                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat,pos_mat_ctrl,i_vol) ]  +=  - control_node_flag[i_vol] *  weight_bdry *
                                                                                (    alpha * phi_ctrl_bdry[i_bdry] * sol_ctrl_bdry_gss
							                           +  beta * lap_rhs_dctrl_ctrl_bdry_gss_i 
							                           - grad_adj_dot_n_res * phi_ctrl_bdry[i_bdry]
// 							                           -         phi_ctrl_bdry[i_bdry]*sol_adj_bdry_gss // for Neumann control
							                         );  //boundary optimality condition
                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat,pos_mat_adj,i_vol) ]  += 0.; 
//============ Bdry Residuals ==================    
		    
		    for(unsigned j_bdry=0; j_bdry < nDof_max_bdry; j_bdry ++) {
		         unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);

//============ Bdry Jacobians ==================	
//============ Bdry Jacobians ==================	


// FIRST BLOCK ROW
//============ u = q ===========================	    
                 
if ( i_vol == j_vol )  {
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_state, i_vol, j_vol) ] += penalty_ctrl * ( control_node_flag[i_vol]);
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_ctrl, i_vol, j_vol) ]  += penalty_ctrl * ( control_node_flag[i_vol]) * (-1.);
		}
//============ u = q ===========================

		    

// SECOND BLOCK ROW
//=========== boundary control eqn =============	    

//========block delta_control / control ========
              double  lap_mat_dctrl_ctrl_bdry_gss = 0.;
		      for (unsigned d = 0; d < space_dim; d++) {  lap_mat_dctrl_ctrl_bdry_gss += phi_ctrl_x_bdry[i_bdry * space_dim + d] * phi_ctrl_x_bdry[j_bdry * space_dim + d];    }

          
              Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_ctrl, i_vol, j_vol) ] 
			+=  control_node_flag[i_vol] *  weight_bdry * (alpha * phi_ctrl_bdry[i_bdry] * phi_ctrl_bdry[j_bdry] 
			                                              + beta *  lap_mat_dctrl_ctrl_bdry_gss);   
    
		   
//============ End Bdry Jacobians ==================	
//============ End Bdry Jacobians ==================	
				
	      }  //end j loop
	      
//===================loop over j in the VOLUME (while i is in the boundary)	      
	for(unsigned j=0; j < nDof_max; j ++) {
		      
  //=============== grad dot n  =========================================    
    double grad_adj_dot_n_mat = 0.;
        for(unsigned d=0; d< space_dim; d++) {
	  grad_adj_dot_n_mat += phi_adj_x_vol_at_bdry[j * space_dim + d] * normal[d];  //notice that the convention of the orders x y z is different from vol to bdry
	}
//=============== grad dot n  =========================================    

		      
//==========block delta_control/adjoint ========
		     Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_adj, i_vol, j) ]  += 
		     control_node_flag[i_vol] * (-1.) * weight_bdry * grad_adj_dot_n_mat * phi_ctrl_bdry[i_bdry];    		      

//==========block delta_state/adjoint ========
		     Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_adj, i_vol, j) ] += 
		     control_node_flag[i_vol] * (1.) * weight_bdry * grad_adj_dot_n_mat * phi_u_bdry[i_bdry];  
		      
		    }   //end loop i_bdry // j_vol
	      
	      

		  }  //end i loop
		}  //end ig_bdry loop
	      }    //end if control face
	      
	    }  //end if boundary faces
	  }    //end loop over faces
	  
	} //end if control element flag
	

//========= gauss value quantities on the volume ==============  
	double sol_u_gss = 0.;
	double sol_adj_gss = 0.;
	std::vector<double> sol_u_x_gss(space_dim);     std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::vector<double> sol_adj_x_gss(space_dim);   std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
//=============================================== 
 
 
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[ielGeom][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_3d(), ig, Jac_qp, JacI_qp, detJac_qp, space_dim);
    weight = detJac_qp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[ig];

    elem_all[ielGeom][SolFEType_quantities[pos_sol_state]]->shape_funcs_current_elem(ig, JacI_qp, phi_u, phi_u_x, phi_u_xx, space_dim);
    elem_all[ielGeom][SolFEType_quantities[pos_sol_adj]]  ->shape_funcs_current_elem(ig, JacI_qp, phi_adj, phi_adj_x, phi_adj_xx, space_dim);

          
	sol_u_gss = 0.;
	sol_adj_gss = 0.;
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < Sol_n_el_dofs_Mat[pos_mat_state]; i++) {
	                                                sol_u_gss      += sol_eldofs_Mat[pos_mat_state][i] * phi_u[i];
                   for (unsigned d = 0; d < space_dim; d++)   sol_u_x_gss[d] += sol_eldofs_Mat[pos_mat_state][i] * phi_u_x[i * space_dim + d];
          }
	
	for (unsigned i = 0; i < Sol_n_el_dofs_Mat[pos_mat_adj]; i++) {
	                                                sol_adj_gss      += sol_eldofs_Mat[pos_mat_adj][i] * phi_adj[i];
                   for (unsigned d = 0; d < space_dim; d++)   sol_adj_x_gss[d] += sol_eldofs_Mat[pos_mat_adj][i] * phi_adj_x[i * space_dim + d];
        }

//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
              double laplace_rhs_du_adj_i = 0.;
              double laplace_rhs_dadj_u_i = 0.;
              for (unsigned kdim = 0; kdim < space_dim; kdim++) {
              if ( i < Sol_n_el_dofs_Mat[pos_mat_state] )  laplace_rhs_du_adj_i +=  phi_u_x   [i * space_dim + kdim] * sol_adj_x_gss[kdim];
              if ( i < Sol_n_el_dofs_Mat[pos_mat_adj] )    laplace_rhs_dadj_u_i +=  phi_adj_x [i * space_dim + kdim] * sol_u_x_gss[kdim];
	      }
	      
//============ Volume residuals ==================	    
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat,pos_mat_state,i) ] += - weight * ( target_flag * phi_u[i] * ( sol_u_gss - u_des)  - laplace_rhs_du_adj_i ); 
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat,pos_mat_ctrl,i) ]  += - penalty_outside_control_boundary * ( (1 - control_node_flag[i]) * (  sol_eldofs_Mat[pos_mat_ctrl][i] - 0.)  );
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat,pos_mat_adj,i) ]   += - weight * (-1.) * (laplace_rhs_dadj_u_i);
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat,pos_mat_mu,i) ]    += - penalty_outside_control_boundary * ( (1 - control_node_flag[i]) * (  sol_eldofs_Mat[pos_mat_mu][i] - 0.)  );  //MU
//============  Volume Residuals ==================	    
	      
	      
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
                
              double laplace_mat_dadj_u = 0.;
              double laplace_mat_du_adj = 0.;

              for (unsigned kdim = 0; kdim < space_dim; kdim++) {
              if ( i < Sol_n_el_dofs_Mat[pos_mat_adj] && j < Sol_n_el_dofs_Mat[pos_mat_state] )     laplace_mat_dadj_u        +=  (phi_adj_x [i * space_dim + kdim] * phi_u_x   [j * space_dim + kdim]);
              if ( i < Sol_n_el_dofs_Mat[pos_mat_state]   && j < Sol_n_el_dofs_Mat[pos_mat_adj] )   laplace_mat_du_adj        +=  (phi_u_x   [i * space_dim + kdim] * phi_adj_x [j * space_dim + kdim]);
		
	      }

              //============ delta_state row ============================
              // BLOCK delta_state / state
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_state, i, j) ]  += weight  * target_flag *  phi_u[i] * phi_u[j];   
              //BLOCK delta_state / adjoint
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_adj, i, j) ]  += weight * (-1) * laplace_mat_du_adj;
	      
	      
              //=========== delta_control row ===========================
              //enforce control zero outside the control boundary
	      if ( i==j )
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_ctrl, i, j) ]  += penalty_outside_control_boundary * ( (1 - control_node_flag[i]));    /*weight * phi_adj[i]*phi_adj[j]*/
              
	      //=========== delta_adjoint row ===========================
	      // BLOCK delta_adjoint / state
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_adj, pos_mat_state, i, j) ]  += weight * (-1) * laplace_mat_dadj_u;

	      
	      //============= delta_mu row ===============================
	        if ( i==j )   
		  Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat, sum_Sol_n_el_dofs, pos_mat_mu, pos_mat_mu, i, j) ]  += penalty_outside_control_boundary * ( (1 - control_node_flag[i]));    //MU
          
	         } // end phi_j loop
           } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop

  

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked(Res, L2G_dofmap_Mat_AllVars);

    if (assembleMatrix) {
      KK->add_matrix_blocked(Jac, L2G_dofmap_Mat_AllVars, L2G_dofmap_Mat_AllVars);
    }
    
    
    //========== dof-based part, without summation
 
     if (print_algebra_local) {
         assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs_Mat, 10, 5);
     }
     
     
  } //end element loop for each process
  

  //MU
add_one_times_mu_res_ctrl_bdry(iproc,
                               ineq_flag,
                               pos_mat_ctrl,
                               pos_mat_mu,
                               SolIndex_Mat,
                               sol,
                               mlPdeSys,
                               pdeSys,
                               RES);

    
  // ***************** END ASSEMBLY - ADD PART *******************

RES->close();
if (assembleMatrix) KK->close();  ///@todo is it needed? I think so


//   ***************** INSERT PART - BEGIN (must go AFTER the sum, clearly) *******************
 // One very important thing to consider: we have some PENALTIES that were set before during the SUMMATION part.
 // Now, if we do INSERT, we may end up OVERWRITING certain values, SUCH AS THOSE PENALTIES!!!
    
     //MU

   for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
       
     geom_element.set_coords_at_dofs_and_geom_type(iel, solType_coords);
      
    el_dofs_unknowns(sol, msh, pdeSys, iel,
                        SolFEType_Mat, SolIndex_Mat, SolPdeIndex,
                        Sol_n_el_dofs_Mat, sol_eldofs_Mat, L2G_dofmap_Mat);

   geom_element.set_elem_center(iel, solType_coords);

 //************ set control flag *********************
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag_bdry(geom_element.get_elem_center());
//  *************************************************** 

// Perform face loop over elements that contain some control face
	if (control_el_flag == 1) {

    	  for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {

       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
// 	    look for boundary faces
            const int bdry_index = el->GetFaceElementIndex(iel, jface);
   
	    if( bdry_index < 0) {
            	      const unsigned int face_in_rectangle_domain = -( msh->el->GetFaceElementIndex(iel,jface)+1);

	      if(  face_in_rectangle_domain == FACE_FOR_CONTROL) { //control face

       update_active_set_flag_for_current_nonlinear_iteration_bdry
   (msh, sol,
    iel, jface,
    geom_element.get_coords_at_dofs_bdry_3d(), 
    sol_eldofs_Mat, 
    Sol_n_el_dofs_Mat, 
    pos_mat_mu,               //this becomes a vector
    pos_mat_ctrl,             //this becomes a vector
    c_compl, 
    ctrl_lower, ctrl_upper,   //this becomes a vector
    sol_actflag,              //this becomes a vector
    solFEType_act_flag_sol, //remove this one, only Index
    solIndex_act_flag_sol);   //this becomes a vector
 

  node_insertion_bdry(iel, jface, 
                      msh,
                      L2G_dofmap_Mat,
                      pos_mat_mu, 
                      pos_mat_ctrl,
                      sol_eldofs_Mat,
                      Sol_n_el_dofs_Mat,
                      sol_actflag, 
                      solFEType_act_flag_sol,  //remove this one, only Index
                      ineq_flag,
                      c_compl,
                      ctrl_lower, ctrl_upper,
                      KK, RES,
                      assembleMatrix
                      );
  
             }
          }
       }
     }


     //============= delta_ctrl-delta_mu row ===============================
 if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap_Mat[pos_mat_ctrl],  L2G_dofmap_Mat[pos_mat_mu], ineq_flag * 1.); }   //this becomes a vector

   }
   
   
  // ***************** INSERT PART - END *******************
RES->close();
if (assembleMatrix) KK->close();  ///@todo is it needed? I think so
    
    
  if (print_algebra_global) {
    if (assembleMatrix) KK->close();
    std::ostringstream mat_out; mat_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "matrix_" << mlPdeSys->GetNonlinearIt()  << ".txt";
    KK->print_matlab(mat_out.str(),"ascii"); //  KK->print();

    RES->close();
    std::ostringstream res_out; res_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "res_" << mlPdeSys->GetNonlinearIt()  << ".txt";
    std::filebuf res_fb;
    res_fb.open (res_out.str().c_str(),std::ios::out);
    std::ostream  res_file_stream(&res_fb);
    RES->print(res_file_stream);
  }
     
     

  return;
}


 
  
void ComputeIntegral(const MultiLevelProblem& ml_prob)    {
  
  
  const NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("BoundaryControl");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

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
 
  CurrentElem < double > geom_element(dim, msh);
    
  constexpr unsigned int space_dim = 3;
  
  std::vector<double> normal(space_dim, 0.);
 //***************************************************

  //=============== Integration ========================================
  double weight; // gauss point weight
  double weight_bdry = 0.; // gauss point weight on the boundary

 //***************************************************
  double alpha = ALPHA_CTRL_BDRY;
  double beta  = BETA_CTRL_BDRY;
  
 //*************** state ***************************** 
 //*************************************************** 
  vector <double> phi_u;     phi_u.reserve(max_size);
  vector <double> phi_u_x;   phi_u_x.reserve(max_size * space_dim);
  vector <double> phi_u_xx;  phi_u_xx.reserve(max_size * dim2);

 
  unsigned solIndex_u;
  solIndex_u = ml_sol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

  vector < double >  sol_u; // local solution
  sol_u.reserve(max_size);
  
  double u_gss = 0.;
 //*************************************************** 
 //***************************************************

  
 //************** desired ****************************
 //***************************************************
  vector <double> phi_udes;
  vector <double> phi_udes_x;
  vector <double> phi_udes_xx;

    phi_udes.reserve(max_size);
    phi_udes_x.reserve(max_size * space_dim);
    phi_udes_xx.reserve(max_size * dim2);
 
  
//  unsigned solIndexTdes;
//   solIndexTdes = ml_sol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solTypeTdes = ml_sol->GetSolutionType(solIndexTdes);    // get the finite element type for "state"

  vector < double >  sol_udes;
  sol_udes.reserve(max_size);

  double udes_gss = 0.;
 //***************************************************
 //***************************************************

 //************** cont *******************************
 //***************************************************
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 
  vector <double> phi_xx_bdry_placeholder;

  phi_ctrl_bdry.reserve(max_size);
  phi_ctrl_x_bdry.reserve(max_size * space_dim);
  phi_xx_bdry_placeholder.reserve(max_size * dim2);

  unsigned solIndex_ctrl = ml_sol->GetIndex("control");
  unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

   vector < double >  sol_ctrl;   sol_ctrl.reserve(max_size);
 //***************************************************
 //*************************************************** 
  
 //********** DATA *********************************** 
  double u_des = DesiredTarget();
 //*************************************************** 
  
  double integral_target = 0.;
  double integral_alpha  = 0.;
  double integral_beta   = 0.;

  

 //*************************************************** 
     std::vector < std::vector < double > >  JacI_qp(space_dim);
     std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < Jac_qp.size(); d++) {   Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp.size(); d++) { JacI_qp[d].resize(dim); }
    
    double detJac_qp;

     std::vector < std::vector < double > >  JacI_qp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_qp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_qp_bdry.size(); d++) {   Jac_qp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp_bdry.size(); d++) { JacI_qp_bdry[d].resize(dim-1); }
    
    double detJac_qp_bdry;
    
      //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
 //*************************************************** 
  
  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    geom_element.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element.geom_type();

  //************* set target domain flag **************
   geom_element.set_elem_center(iel, solType_coords);

   int target_flag = 0;
   target_flag = ElementTargetFlag(geom_element.get_elem_center());
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

 // ==================================================
 //****** set control flag ***************************
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag_bdry(geom_element.get_elem_center());
 //***************************************************

  
  	if (control_el_flag == 1) {
	  
	  double tau=0.;
	       
	  // loop on faces of the current element

	  for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, jface);    
       const unsigned nve_bdry_ctrl = msh->GetElementFaceDofNumber(iel,jface,solType_ctrl);
       

       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element.set_elem_center_bdry_3d();

       
	    // look for boundary faces
            const int bdry_index = el->GetFaceElementIndex(iel,jface);
            
	    if( bdry_index < 0) {
	      unsigned int face = -( msh->el->GetFaceElementIndex(iel,jface)+1);
	      
		
// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
	      if(  face == FACE_FOR_CONTROL) { //control face

	
		//============ initialize gauss quantities on the boundary ==========================================
                double sol_ctrl_bdry_gss = 0.;
                std::vector<double> sol_ctrl_x_bdry_gss(space_dim);
		//============ initialize gauss quantities on the boundary ==========================================
		
		for(unsigned ig_bdry = 0; ig_bdry < ml_prob.GetQuadratureRule(ielGeom_bdry).GetGaussPointsNumber(); ig_bdry++) {
		  
    elem_all[ielGeom_bdry][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_qp_bdry, JacI_qp_bdry, detJac_qp_bdry, space_dim);
    weight_bdry = detJac_qp_bdry * ml_prob.GetQuadratureRule(ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];
    elem_all[ielGeom_bdry][solType_ctrl] ->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_ctrl_bdry, phi_ctrl_x_bdry, phi_xx_bdry_placeholder, space_dim);

		  
		 //========== compute gauss quantities on the boundary ===============================================
		  sol_ctrl_bdry_gss = 0.;
                  std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);
		      for (int i_bdry = 0; i_bdry < nve_bdry_ctrl; i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
			
			sol_ctrl_bdry_gss +=  sol_ctrl[i_vol] * phi_ctrl_bdry[i_bdry];
                            for (int d = 0; d < space_dim; d++) {
			      sol_ctrl_x_bdry_gss[d] += sol_ctrl[i_vol] * phi_ctrl_x_bdry[i_bdry * space_dim + d];
			    }
		      }
		      
		      double laplace_ctrl_surface = 0.;  for (int d = 0; d < space_dim; d++) { laplace_ctrl_surface += sol_ctrl_x_bdry_gss[d] * sol_ctrl_x_bdry_gss[d]; }

                 //========= compute gauss quantities on the boundary ================================================
                  integral_alpha +=  weight_bdry * sol_ctrl_bdry_gss * sol_ctrl_bdry_gss; 
                  integral_beta  +=  weight_bdry * laplace_ctrl_surface;
                 
             }
	      } //end face == 3
	      
	    } //end if boundary faces
	  }  // loop over element faces   
	  
	} //end if control element flag

//=====================================================================================================================  
//=====================================================================================================================  
//=====================================================================================================================  
  
  
   
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[ielGeom][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_3d(), ig, Jac_qp, JacI_qp, detJac_qp, space_dim);
    weight = detJac_qp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[ig];

    elem_all[ielGeom][solType_u]                 ->shape_funcs_current_elem(ig, JacI_qp, phi_u, phi_u_x, phi_u_xx, space_dim);
    elem_all[ielGeom][solType_u/*solTypeTdes*/]  ->shape_funcs_current_elem(ig, JacI_qp, phi_udes, phi_udes_x, phi_udes_xx, space_dim);
    
	u_gss     = 0.;  for (unsigned i = 0; i < nDof_u; i++)        u_gss += sol_u[i]     * phi_u[i];
	udes_gss  = 0.;  for (unsigned i = 0; i < nDof_udes; i++)  udes_gss += sol_udes[i]  * phi_udes[i];  

               integral_target += target_flag * weight * (u_gss  - udes_gss) * (u_gss - udes_gss);
	  
      } // end gauss point loop
      
  } //end element loop

  double total_integral = 0.5 * integral_target + 0.5 * alpha * integral_alpha + 0.5 * beta * integral_beta;
  
  
  ////////////////////////////////////////
       std::cout << "integral on processor: " << total_integral << std::endl;

   double J = 0.;
      MPI_Allreduce( &total_integral, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  //THIS IS THE RIGHT ONE!!


    std::cout << "@@@@@@@@@@@@@@@@ functional value: " << J << std::endl;
  
//   std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << integral_target << std::endl;
//   std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << integral_alpha << std::endl;
//   std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << integral_beta << std::endl;
//   std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << total_integral << std::endl;
 
return;
  
}
  
