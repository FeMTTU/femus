#ifndef OPT_SYSTEMS_DIRICHLET_HPP
#define OPT_SYSTEMS_DIRICHLET_HPP


#include "MultiLevelProblem.hpp"
#include "SparseMatrix.hpp"

// #include "Assemble_jacobian.hpp"
// #include "Assemble_unknown_jacres.hpp"
// #include "Assemble_unknown.hpp"

#include "01_opt_system.hpp"

#include  "03_opt_system_inequalities.hpp"

#include  "opt_systems_boundary_control_eqn_sobolev_integer.hpp"
  

//This Opt system is characterized by the following ways of setting matrix values:
// Add_values (Mat or Vec) in the volume loop
// Add_values (Mat or Vec) in the boundary loop
// Insert_values (Mat or Vec) in the boundary loop
// Insert_values (Mat or Vec) in the volume loop
// Insert_values (Mat or Vec) outside all loops
// We're going to split the two parts and add a close() at the end of each

namespace femus  {

namespace elliptic  {


template < bool IS_BOUNDARY_CONTROL_INTEGER_OR_FRACTIONAL >    
class pure_boundary : public femus::pure_boundary  {

private:
    
    static constexpr bool _is_boundary_control_integer_or_fractional = IS_BOUNDARY_CONTROL_INTEGER_OR_FRACTIONAL;

    
public: 

 static const std::vector< Unknown >  provide_list_of_unknowns(const unsigned int dimension) {
     
     
  std::vector< FEFamily > feFamily;
  std::vector< FEOrder >   feOrder;

                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
 
                        feOrder.push_back(/*FIRST*/SECOND);  //same
                        feOrder.push_back(/*FIRST*/SECOND);  //same
                        feOrder.push_back(/*FIRST*/SECOND);
                        feOrder.push_back(/*FIRST*/SECOND);  //same
 

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Unknown >  unknowns(feFamily.size());

   unknowns[0]._name      = "state";
   unknowns[1]._name      = "control";
   unknowns[2]._name      = "adjoint";
   unknowns[3]._name      = "mu";

   unknowns[0]._is_sparse = true;
   unknowns[1]._is_sparse = _is_boundary_control_integer_or_fractional ? false: true;
   unknowns[2]._is_sparse = true;
   unknowns[3]._is_sparse = true;
   
     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._fe_family  = feFamily[u];
              unknowns[u]._fe_order   = feOrder[u];
              unknowns[u]._time_order = 0;
              unknowns[u]._is_pde_unknown = true;
              
     }
 
   return unknowns;
     
}

    
 static void assemble_elliptic_dirichlet_control(MultiLevelProblem & ml_prob) {
    
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = & ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("BoundaryControl");
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->GetMeshElements();
  
  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*            JAC = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors();

  constexpr bool print_algebra_global = false;
  constexpr bool print_algebra_local = false;
  
  

  //=============== Geometry ========================================
   unsigned solType_coords = CONTINUOUS_BIQUADRATIC;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
  CurrentElem < double > geom_element_jel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  
  std::vector< double > normal(space_dim, 0.);
  
 //************ geom ***************************************  
  std::vector < double > coord_at_qp_bdry(space_dim);
  
  std::vector < double > phi_coords;
  std::vector < double > phi_coords_x;

  phi_coords.reserve(max_size);
  phi_coords_x.reserve(max_size * space_dim);
  

 //********************* state *********************** 
 //*************************************************** 
  std::vector <double> phi_u;
  std::vector <double> phi_u_x;

  phi_u.reserve(max_size);
  phi_u_x.reserve(max_size * space_dim);
  
  
  //boundary state shape functions
  std::vector <double> phi_u_bdry;  
  std::vector <double> phi_u_x_bdry; 

  phi_u_bdry.reserve(max_size);
  phi_u_x_bdry.reserve(max_size * space_dim);
  
 //***************************************************  
 //***************************************************  

  
 //********************** adjoint ********************
 //*************************************************** 
  std::vector <double> phi_adj;
  std::vector <double> phi_adj_x;

  phi_adj.reserve(max_size);
  phi_adj_x.reserve(max_size * space_dim);
 

  //boundary adjoint shape functions  
  std::vector <double> phi_adj_bdry;  
  std::vector <double> phi_adj_x_bdry; 

  phi_adj_bdry.reserve(max_size);
  phi_adj_x_bdry.reserve(max_size * space_dim);

  
  //volume shape functions at boundary
  std::vector <double> phi_adj_vol_at_bdry;
  std::vector <double> phi_adj_x_vol_at_bdry;
  phi_adj_vol_at_bdry.reserve(max_size);
  phi_adj_x_vol_at_bdry.reserve(max_size * space_dim);
  
  std::vector <double> sol_adj_x_vol_at_bdry_gss(space_dim);
 //*************************************************** 
 //*************************************************** 

  
 //********************* bdry cont *******************
 //*************************************************** 
  std::vector <double> phi_ctrl_bdry;  
  std::vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(max_size);
  phi_ctrl_x_bdry.reserve(max_size * space_dim);
 //*************************************************** 

    const unsigned int n_components_ctrl = 1;
    const unsigned int first_loc_comp_ctrl = 0;


    //MU
  //************** act flag ****************************   
    std::vector <unsigned int> solIndex_act_flag_sol(n_components_ctrl);

    femus::ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::store_act_flag_in_old(mlPdeSys, ml_sol, sol, solIndex_act_flag_sol);
  
  
  //********* variables for ineq constraints *****************
     std::vector <std::vector < double/*int*/ > > sol_actflag(n_components_ctrl);    //flag for active set
     std::vector <std::vector < double > > ctrl_lower(n_components_ctrl);  
     std::vector <std::vector < double > > ctrl_upper(n_components_ctrl);  
  
     for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {
          sol_actflag[kdim].reserve(max_size);
           ctrl_lower[kdim].reserve(max_size);
           ctrl_upper[kdim].reserve(max_size);
     }
     
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

    std::vector < std::string > Solname_Mat(n_unknowns);  //this coincides with Pos_in_matrix
    Solname_Mat[0] = "state";
    Solname_Mat[1] = "control";
    Solname_Mat[2] = "adjoint";
    Solname_Mat[3] = "mu";

 //***************************************************
   enum Pos_in_Sol {pos_sol_state = 0, pos_sol_ctrl, pos_sol_adj, pos_sol_mu, pos_sol_targreg, pos_sol_contreg, pos_sol_actflag}; //these are known at compile-time 

        assert(pos_sol_actflag == solIndex_act_flag_sol[0]);
        
        
  std::vector<unsigned int> ctrl_index_in_mat(1);  ctrl_index_in_mat[0] = pos_mat_ctrl;
  std::vector<unsigned int>   mu_index_in_mat(1);    mu_index_in_mat[0] = pos_mat_mu;

//***************************************************
    
 //***************************************************
    std::vector < std::string > Solname_quantities(n_quantities);
    
        for(unsigned ivar=0; ivar < Solname_quantities.size(); ivar++) {
            Solname_quantities[ivar] = ml_sol->GetSolName_from_index(ivar);
        }
 //***************************************************
        
    std::vector < unsigned > SolIndex_Mat(n_unknowns);      //should have Mat order
    std::vector < unsigned > SolFEType_Mat(n_unknowns);       //should have Mat order
    std::vector < unsigned > SolPdeIndex(n_unknowns);     //should have Mat order, of course

    std::vector < unsigned > SolIndex_quantities(n_quantities);      //should have Sol order
    std::vector < unsigned > SolFEType_quantities(n_quantities);     //should have Sol order
    std::vector < unsigned > Sol_n_el_dofs_quantities_vol(n_quantities); //should have Sol order
 
  

    for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
        SolIndex_Mat[ivar]    = ml_sol->GetIndex        (Solname_Mat[ivar].c_str());
        SolFEType_Mat[ivar]   = ml_sol->GetSolutionType(SolIndex_Mat[ivar]);
        SolPdeIndex[ivar] = mlPdeSys->GetSolPdeIndex(Solname_Mat[ivar].c_str());
    }
    
    for(unsigned ivar=0; ivar < n_quantities; ivar++) {
        SolIndex_quantities[ivar]    = ml_sol->GetIndex        (Solname_quantities[ivar].c_str());
        SolFEType_quantities[ivar]   = ml_sol->GetSolutionType(SolIndex_quantities[ivar]);
    }    

    std::vector < unsigned > Sol_n_el_dofs_Mat_vol(n_unknowns);    //should have Mat order

//***************************************************
    std::vector < std::vector < double > >  Sol_eldofs_Mat(n_unknowns);  //should have Mat order
    for(int k = 0; k < n_unknowns; k++) {        Sol_eldofs_Mat[k].reserve(max_size);    }


    //----------- quantities (at dof objects) ------------------------------
    std::vector< int >       L2G_dofmap_Mat_AllVars;
    L2G_dofmap_Mat_AllVars.reserve( n_unknowns * max_size );
    std::vector < std::vector < int > >     L2G_dofmap_Mat(n_unknowns);     //should have Mat order
    for(int i = 0; i < n_unknowns; i++) {
        L2G_dofmap_Mat[i].reserve(max_size);
    }
    
 //*************************************************** 
  std::vector < double > Res;   Res.reserve( n_unknowns*max_size);                         //should have Mat order
  std::vector < double > Jac;   Jac.reserve( n_unknowns*max_size * n_unknowns*max_size);   //should have Mat order
 //*************************************************** 

 
 //********************* DATA ************************ 
  const double u_des = femus::ctrl::square_or_cube :: cost_functional_without_regularization::DesiredTargetVec()[0];
  const double alpha = ALPHA_CTRL_BDRY;
  const double beta  = BETA_CTRL_BDRY;
  const double penalty_outside_control_domain_boundary = PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY;       // penalty for zero control outside Gamma_c and zero mu outside Gamma_c
  const double penalty_dirichlet_bc_u_equal_q = PENALTY_DIRICHLET_BC_U_EQUAL_Q_BOUNDARY;         //penalty for u=q
 //*************************************************** 
  
  RES->zero();  
  if (assembleMatrix)  JAC->zero();

 //*************************************************** 
// ---
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
    double detJac_iqp;
  double weight_iqp = 0.;
// ---

// ---
    std::vector < std::vector < double > >  JacI_iqp_bdry(space_dim);
    std::vector < std::vector < double > >  Jac_iqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_iqp_bdry.size(); d++) {   Jac_iqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp_bdry.size(); d++) { JacI_iqp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
  double weight_iqp_bdry = 0.;

// ---
    
//*************************************************** 

    //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all;
  ml_prob.get_all_abstract_fe_multiple(elem_all);
//*************************************************** 
//***************************************************
  
  const unsigned dim_bdry = dim - 1;
  
  const double s_frac = _s_frac;

  const double check_limits = 1.;//1./(1. - s_frac); // - s_frac;

  


  //--- quadrature rules -------------------
  constexpr unsigned qrule_i = QRULE_I;
  constexpr unsigned qrule_j = QRULE_J;
  constexpr unsigned qrule_k = QRULE_K;
  //----------------------


    
  if ( _is_boundary_control_integer_or_fractional ) {
     // fractional
     femus::ctrl::Gamma_control_equation_fractional_sobolev_differentiability_index<
                femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES, 
                femus::ctrl:: square_or_cube:: pure_boundary< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >
                >::control_eqn_bdry(
                    iproc,
                    nprocs,
                    ml_prob,
                    ml_sol,
                    sol,
                    msh,
                    pdeSys,
                    //-----------
                    geom_element_iel,
                    solType_coords,
                    dim,
                    space_dim,
                    dim_bdry,
                    //-----------
                    n_unknowns,
                    Solname_Mat,
                    SolFEType_Mat,
                    SolIndex_Mat,
                    SolPdeIndex,
                    Sol_n_el_dofs_Mat_vol, 
                    Sol_eldofs_Mat,  
                    L2G_dofmap_Mat,
                    max_size,
                    //-----------
                    n_quantities,
                    SolFEType_quantities,
//                     Sol_n_el_dofs_quantities, //filled inside
                    //-----------
                    elem_all,
                     //-----------
                    Jac_iqp_bdry,
                    JacI_iqp_bdry,
                    detJac_iqp_bdry,
                    weight_iqp_bdry,
                    phi_ctrl_bdry,
                    phi_ctrl_x_bdry, 
                    //-----------
                    n_components_ctrl,
                    pos_mat_ctrl,
                    pos_sol_ctrl,
                    _is_block_dctrl_ctrl_inside_main_big_assembly,
                    //-----------
                    JAC,
                    RES,
                    assembleMatrix,
                    //-----------
                    alpha,
                    beta,     
                    s_frac,
                    check_limits,
                    USE_Cns,
                    OP_Hhalf,
                    OP_L2,
                    _rhs_one,
                    UNBOUNDED,
                    _node_based_bdry_bdry,
                    //-----------
                    qrule_i,
                    qrule_j,
                    qrule_k,
                    Nsplit,
                    Quadrature_split_index,
                    N_DIV_FACE_OF_FACE_FOR_UNBOUNDED_INTEGRAL,
                    //-----------
                    print_algebra_local
                    );
                    
  }
  
  else {
   // integer
   femus::ctrl::Gamma_control_equation_integer_sobolev_differentiability_index<
                femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES, 
                femus::ctrl:: square_or_cube:: pure_boundary< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >
                >
                ::control_eqn_bdry(iproc,
                    ml_prob,
                    ml_sol,
                    sol,
                    msh,
                    pdeSys,
                    //-----------
                    geom_element_iel,
                    solType_coords,
                    space_dim,
                    //-----------
                    n_unknowns,
                    Solname_Mat,
                    SolFEType_Mat,
                    SolIndex_Mat,
                    SolPdeIndex,
                    Sol_n_el_dofs_Mat_vol, 
                    Sol_eldofs_Mat,  
                    L2G_dofmap_Mat,
                    max_size,
                    //-----------
                     n_quantities,
                    SolFEType_quantities,
                    //-----------
                    elem_all,
                     //-----------
                    Jac_iqp_bdry,
                    JacI_iqp_bdry,
                    detJac_iqp_bdry,
                    weight_iqp_bdry,
                    phi_ctrl_bdry,
                    phi_ctrl_x_bdry, 
                    //-----------
                    n_components_ctrl,
                    pos_mat_ctrl,
                    pos_sol_ctrl,
                    _is_block_dctrl_ctrl_inside_main_big_assembly,
                    //-----------
                    JAC,
                    RES,
                    assembleMatrix,
                    //-----------
                    alpha,
                    beta,
                    _rhs_one,
                    OP_L2,
                    OP_H1,
                    qrule_i,
                    //-----------
                    print_algebra_local
                    ) ;
  
  }
  
  

//**************************                    
// AAA do not close this because later they will be filled with the rest of the system!!!      
//    JAC->close();
//   RES->close();

   //print JAC and RES to files
   //You print only what is being sent to the global matrix. If nothing is sent, nothing is printed                 
   // I will keep this print here for later because it highlights what positions were filled in the matrix
   // If I remove everything above here, it seems like the very last diagonal position is filled... why? from where?
//   JAC->close(); //JAC->zero();
//   const unsigned nonlin_iter = 0/*mlPdeSys->GetNonlinearIt()*/;
//     assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, JAC, nonlin_iter);
// //     assemble_jacobian< double, double >::print_global_residual(ml_prob, RES, nonlin_iter);
//   std::cout << "****************************" << std::endl;
//**************************                    

  
                    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

// -------
    geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element_iel.geom_type();

    geom_element_iel.set_elem_center_3d(iel, solType_coords);
// -------
    

 //***************************************************
   el_dofs_unknowns_vol(sol, msh, pdeSys, iel,
                        SolFEType_Mat,
                        SolIndex_Mat,
                        SolPdeIndex,
                        Sol_n_el_dofs_Mat_vol, 
                        Sol_eldofs_Mat,  
                        L2G_dofmap_Mat);
  
   el_dofs_quantities_vol(sol, msh, iel, SolFEType_quantities, Sol_n_el_dofs_quantities_vol); 
  //***************************************************
   
 
    unsigned int nDof_max          = ElementJacRes<double>::compute_max_n_dofs(Sol_n_el_dofs_Mat_vol);
    
    unsigned int sum_Sol_n_el_dofs = ElementJacRes<double>::compute_sum_n_dofs(Sol_n_el_dofs_Mat_vol);

    Res.resize(sum_Sol_n_el_dofs);
    std::fill(Res.begin(), Res.end(), 0.);

    Jac.resize(sum_Sol_n_el_dofs * sum_Sol_n_el_dofs);
    std::fill(Jac.begin(), Jac.end(), 0.);
    
    L2G_dofmap_Mat_AllVars.resize(0);
      for (unsigned  k = 0; k < n_unknowns; k++)     L2G_dofmap_Mat_AllVars.insert(L2G_dofmap_Mat_AllVars.end(), L2G_dofmap_Mat[k].begin(), L2G_dofmap_Mat[k].end());
 //***************************************************

      
  //************* set target domain flag **************
   int target_flag = 0;
   target_flag = femus::ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
 //*************************************************** 
   

 //************ set control flag *********************
   std::vector< std::vector< int > > control_node_flag = 
       femus::is_dof_associated_to_Gamma_control_equation(msh, ml_sol, & ml_prob, iel, ielGeom, geom_element_iel, solType_coords, Solname_Mat, SolFEType_Mat, Sol_n_el_dofs_Mat_vol, pos_mat_ctrl, n_components_ctrl);
  //*************************************************** 
 

	if ( femus::ctrl:: square_or_cube:: pure_boundary< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES > ::volume_elem_contains_a_Gamma_control_face(sol, msh, iel) ) {
	  
	  std::vector<double> normal(space_dim, 0.);
	       
	  // loop on faces of the current element

	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
// -------
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
// -------
       
// -------
       std::vector<unsigned int> Sol_el_n_dofs_current_face(n_unknowns); ///@todo the active flag is not an unknown! However, if I add to the quantities something that has higher order, then I will have error below when I take the max number of dofs on a face. Since the act_flag must have the same FE family as the control, then I can limit this array to n_unknowns instead of n_quantities

       for (unsigned  k = 0; k < Sol_el_n_dofs_current_face.size(); k++) {
                 if (SolFEType_quantities[k] < 3) Sol_el_n_dofs_current_face[k] = msh->GetElementFaceDofNumber(iel, iface, SolFEType_quantities[k]);  ///@todo fix this absence
       }
       
       const unsigned nDof_max_bdry = ElementJacRes<double>::compute_max_n_dofs(Sol_el_n_dofs_current_face);
// -------

        std::pair< int, unsigned int > pair_control_iface = femus::face_is_a_Gamma_control_face_of_some_index< femus::ctrl:: GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >(msh->GetMeshElements(), iel, iface);

        const int  iface_is_a_boundary_control  = pair_control_iface.first;

	    if( iface_is_a_boundary_control ) {
		
//========= initialize gauss quantities on the boundary ============================================
                double sol_ctrl_bdry_gss = 0.;
                double sol_adj_bdry_gss = 0.;
                std::vector<double> sol_ctrl_x_bdry_gss(space_dim);   

//========= initialize gauss quantities on the boundary ============================================
		
        const unsigned n_gauss_bdry = ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber();
        
    
		for(unsigned ig_bdry = 0; ig_bdry < n_gauss_bdry; ig_bdry++) {
    
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_iqp_bdry, JacI_iqp_bdry, detJac_iqp_bdry, space_dim);
	elem_all[qrule_i][ielGeom_bdry][solType_coords]->compute_normal(Jac_iqp_bdry, normal);
    
    weight_iqp_bdry = detJac_iqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];

    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_ctrl]] ->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_ctrl_bdry, phi_ctrl_x_bdry, boost::none, space_dim);
    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_state]]->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_u_bdry, phi_u_x_bdry,  boost::none, space_dim);
    elem_all[qrule_i][ielGeom_bdry][SolFEType_quantities[pos_sol_adj]]  ->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_adj_bdry, phi_adj_x_bdry,  boost::none, space_dim);


    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv_vol_at_bdry_new(geom_element_iel.get_coords_at_dofs_3d(), ig_bdry, iface, Jac_iqp/*not_needed_here*/, JacI_iqp, detJac_iqp/*not_needed_here*/, space_dim);
    elem_all[qrule_i][ielGeom][SolFEType_quantities[pos_sol_adj]]->shape_funcs_vol_at_bdry_current_elem(ig_bdry, iface, JacI_iqp, phi_adj_vol_at_bdry, phi_adj_x_vol_at_bdry, boost::none, space_dim);
     
//     msh->_finiteElement[ielGeom][SolFEType_quantities[pos_sol_adj]]->fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem(geom_element_iel.get_coords_at_dofs(), geom_element_iel.get_coords_at_dofs_bdry_3d(), iface, ig_bdry, phi_adj_vol_at_bdry, phi_adj_x_vol_at_bdry);

		  
//========== compute gauss quantities on the boundary ===============================================
		  sol_ctrl_bdry_gss = 0.;
                  std::fill(sol_ctrl_x_bdry_gss.begin(), sol_ctrl_x_bdry_gss.end(), 0.);
                  
		      for (unsigned int i_bdry = 0; i_bdry < phi_ctrl_bdry.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_ctrl_bdry_gss +=  Sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_bdry[i_bdry];
                            for (unsigned int d = 0; d < space_dim; d++) {
			      sol_ctrl_x_bdry_gss[d] += Sol_eldofs_Mat[pos_mat_ctrl][i_vol] * phi_ctrl_x_bdry[i_bdry * space_dim + d];
			    }
		      }
		      
		      
		  sol_adj_bdry_gss = 0.;
		      for (unsigned int i_bdry = 0; i_bdry <  phi_adj_bdry.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_adj_bdry_gss  +=  Sol_eldofs_Mat[pos_mat_adj][i_vol] * phi_adj_bdry[i_bdry];
              }		      
		      
//=============== grad dot n for residual ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
           std::fill(sol_adj_x_vol_at_bdry_gss.begin(), sol_adj_x_vol_at_bdry_gss.end(), 0.);
		      for (unsigned int iv = 0; iv < Sol_n_el_dofs_Mat_vol[pos_mat_adj]; iv++)  {
			
         for (unsigned int d = 0; d < space_dim; d++) {
			      sol_adj_x_vol_at_bdry_gss[d] += Sol_eldofs_Mat[pos_mat_adj][iv] * phi_adj_x_vol_at_bdry[iv * space_dim + d];
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
              
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);

                 double lap_rhs_dctrl_ctrl_bdry_gss_i = 0.;
                 for (unsigned d = 0; d < space_dim; d++) {
                       if ( i_vol < Sol_n_el_dofs_Mat_vol[pos_mat_ctrl] )  lap_rhs_dctrl_ctrl_bdry_gss_i +=  phi_ctrl_x_bdry[i_bdry * space_dim + d] * sol_ctrl_x_bdry_gss[d];
                 }
                 
		 
//============ Bdry Residuals - BEGIN  ==================	
                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_state, i_vol) ] +=
                    - control_node_flag[first_loc_comp_ctrl][i_vol] * penalty_dirichlet_bc_u_equal_q     * _keep_adjoint_push * (   Sol_eldofs_Mat[pos_mat_state][i_vol] - Sol_eldofs_Mat[pos_mat_ctrl][i_vol] )
                    - control_node_flag[first_loc_comp_ctrl][i_vol] *  weight_iqp_bdry * _keep_adjoint_push * grad_adj_dot_n_res * phi_u_bdry[i_bdry];   // u = q


                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_ctrl, i_vol) ]  += 
                    - control_node_flag[first_loc_comp_ctrl][i_vol] *  weight_iqp_bdry *
                          (    _is_block_dctrl_ctrl_inside_main_big_assembly * OP_L2 * alpha * phi_ctrl_bdry[i_bdry]  * sol_ctrl_bdry_gss
							+  _is_block_dctrl_ctrl_inside_main_big_assembly * OP_H1 * beta * lap_rhs_dctrl_ctrl_bdry_gss_i
							                            - _keep_adjoint_push * grad_adj_dot_n_res * phi_ctrl_bdry[i_bdry]
// 							                           -         phi_ctrl_bdry[i_bdry]*sol_adj_bdry_gss // for Neumann control
						  );  //boundary optimality condition
                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_adj, i_vol) ]  += 0.; 
//============ Bdry Residuals - END  ==================    
		    
		    for(unsigned j_bdry=0; j_bdry < nDof_max_bdry; j_bdry ++) {
		         unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, iface, j_bdry);

//============ Bdry Jacobians - BEGIN ==================	


// FIRST BLOCK ROW
//============ u = q - BEGIN ===========================	    
                 
if ( i_vol == j_vol )  {
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_state, i_vol, j_vol) ] += 
		     penalty_dirichlet_bc_u_equal_q * _keep_adjoint_push * ( control_node_flag[first_loc_comp_ctrl][i_vol]) * ( 1.);
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_ctrl, i_vol, j_vol) ]  += 
		     penalty_dirichlet_bc_u_equal_q * _keep_adjoint_push * ( control_node_flag[first_loc_comp_ctrl][i_vol]) * (-1.);
		}
//============ u = q - END ===========================

		    

// SECOND BLOCK ROW
//=========== boundary control eqn =============	    

//========block delta_control / control ========
              double  lap_mat_dctrl_ctrl_bdry_gss = 0.;
		      for (unsigned d = 0; d < space_dim; d++) {  lap_mat_dctrl_ctrl_bdry_gss += phi_ctrl_x_bdry[i_bdry * space_dim + d] * phi_ctrl_x_bdry[j_bdry * space_dim + d];    }

          
              Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_ctrl, i_vol, j_vol) ] 
			+=  control_node_flag[first_loc_comp_ctrl][i_vol] *  weight_iqp_bdry * ( _is_block_dctrl_ctrl_inside_main_big_assembly * OP_L2 * alpha * phi_ctrl_bdry[i_bdry] * phi_ctrl_bdry[j_bdry]
			                                              +  _is_block_dctrl_ctrl_inside_main_big_assembly * OP_H1 * beta *  lap_mat_dctrl_ctrl_bdry_gss);
    
		   
//============ Bdry Jacobians - END ==================	
				
	      }  //end j loop
	      
//===================loop over j in the VOLUME (while i is in the boundary)	      
	for(unsigned j=0; j < nDof_max; j ++) {
		      
  //=============== grad dot n  =========================================    
    double grad_adj_dot_n_mat = 0.;
        for(unsigned d = 0; d< space_dim; d++) {
	  grad_adj_dot_n_mat += phi_adj_x_vol_at_bdry[j * space_dim + d] * normal[d];  //notice that the convention of the orders x y z is different from vol to bdry
	}
//=============== grad dot n  =========================================    

		      
//==========block delta_control/adjoint ========
		     Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_adj, i_vol, j) ]  += 
		     control_node_flag[first_loc_comp_ctrl][i_vol] * (-1.) * weight_iqp_bdry * _keep_adjoint_push * grad_adj_dot_n_mat * phi_ctrl_bdry[i_bdry];    		      

//==========block delta_state/adjoint ========
		     Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_adj, i_vol, j) ] += 
		     control_node_flag[first_loc_comp_ctrl][i_vol] * (1.) * weight_iqp_bdry * _keep_adjoint_push * grad_adj_dot_n_mat * phi_u_bdry[i_bdry];  
		      
		    }   //end loop i_bdry // j_vol
	      
	      

		  }  //end i loop
		}  //end ig_bdry loop
		
	      }    //end if control face
	      
	  }    //end loop over faces
	  
	} //end if control element flag
	

//========= gauss value quantities on the volume ==============  
	double sol_u_gss = 0.;
	double sol_adj_gss = 0.;
	std::vector<double> sol_u_x_gss(space_dim);     std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::vector<double> sol_adj_x_gss(space_dim);   std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
//=============================================== 
 
 
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom).GetGaussPointsNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), ig, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);
    weight_iqp = detJac_iqp * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom).GetGaussWeightsPointer()[ig];

    elem_all[qrule_i][ielGeom][SolFEType_quantities[pos_sol_state]]->shape_funcs_current_elem(ig, JacI_iqp, phi_u, phi_u_x, boost::none, space_dim);
    elem_all[qrule_i][ielGeom][SolFEType_quantities[pos_sol_adj]]  ->shape_funcs_current_elem(ig, JacI_iqp, phi_adj, phi_adj_x, boost::none, space_dim);

          
	sol_u_gss = 0.;
	sol_adj_gss = 0.;
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[pos_mat_state]; i++) {
	                                                sol_u_gss      += Sol_eldofs_Mat[pos_mat_state][i] * phi_u[i];
                   for (unsigned d = 0; d < space_dim; d++)   sol_u_x_gss[d] += Sol_eldofs_Mat[pos_mat_state][i] * phi_u_x[i * space_dim + d];
          }
	
	for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[pos_mat_adj]; i++) {
	                                                sol_adj_gss      += Sol_eldofs_Mat[pos_mat_adj][i] * phi_adj[i];
                   for (unsigned d = 0; d < space_dim; d++)   sol_adj_x_gss[d] += Sol_eldofs_Mat[pos_mat_adj][i] * phi_adj_x[i * space_dim + d];
        }

//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {

              double laplace_rhs_du_adj_i = 0.;
              double laplace_rhs_dadj_u_i = 0.;
              double laplace_rhs_du_u_i = 0.;

              for (unsigned kdim = 0; kdim < space_dim; kdim++) {
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_state] )  laplace_rhs_du_adj_i  +=  phi_u_x   [i * space_dim + kdim] * sol_adj_x_gss[kdim];
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_adj] )    laplace_rhs_dadj_u_i  +=  phi_adj_x [i * space_dim + kdim] * sol_u_x_gss[kdim];
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_state] )  laplace_rhs_du_u_i    +=  phi_u_x   [i * space_dim + kdim] * sol_u_x_gss[kdim];
	      }

//============ Volume residuals - BEGIN  ==================	    
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_state, i) ] += - weight_iqp * (
          - laplace_rhs_du_adj_i
          + target_flag *
          #if COST_FUNCTIONAL_TYPE == 0
                phi_u[i] * ( sol_u_gss - u_des)
          #elif COST_FUNCTIONAL_TYPE == 1
                laplace_rhs_du_u_i
          #endif
          );
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_ctrl, i) ]  += - penalty_outside_control_domain_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]) * (  Sol_eldofs_Mat[pos_mat_ctrl][i] - PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY)  );
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_adj, i) ]   += - weight_iqp * (-1.) * (laplace_rhs_dadj_u_i);
          Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mat_mu, i) ]    += - penalty_outside_control_domain_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]) * (  Sol_eldofs_Mat[pos_mat_mu][i] - 0.)  );  //MU
//============  Volume Residuals - END ==================	    
	      
	      
//============ Volume Jacobians - BEGIN  ==================	    
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
                
              double laplace_mat_dadj_u = 0.;
              double laplace_mat_du_adj = 0.;
              double laplace_mat_du_u = 0.;

              for (unsigned kdim = 0; kdim < space_dim; kdim++) {
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_adj]     && j < Sol_n_el_dofs_Mat_vol[pos_mat_state] )     laplace_mat_dadj_u        +=  (phi_adj_x [i * space_dim + kdim] * phi_u_x   [j * space_dim + kdim]);
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_state]   && j < Sol_n_el_dofs_Mat_vol[pos_mat_adj] )   laplace_mat_du_adj        +=  (phi_u_x   [i * space_dim + kdim] * phi_adj_x [j * space_dim + kdim]);
              if ( i < Sol_n_el_dofs_Mat_vol[pos_mat_state]   && j < Sol_n_el_dofs_Mat_vol[pos_mat_state] )     laplace_mat_du_u          += (phi_u_x  [i * space_dim + kdim] * phi_u_x  [j * space_dim + kdim]);
		
	      }

              //============ delta_state row ============================
              // BLOCK delta_state / state
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_state, i, j) ]  += weight_iqp  * target_flag *
		#if COST_FUNCTIONAL_TYPE == 0
             phi_u[i] * phi_u[j]
        #elif COST_FUNCTIONAL_TYPE == 1
             laplace_mat_du_u
        #endif
        ;

              //BLOCK delta_state / adjoint
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_state, pos_mat_adj, i, j) ]  += weight_iqp * (-1.) * laplace_mat_du_adj;
	      
	      
              //=========== delta_control row ===========================
              //enforce control zero outside the control boundary
	      if ( i==j )
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_ctrl, pos_mat_ctrl, i, j) ]  += penalty_outside_control_domain_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]));    /*weight * phi_adj[i]*phi_adj[j]*/
              
	      //=========== delta_adjoint row ===========================
	      // BLOCK delta_adjoint / state
		Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_adj, pos_mat_state, i, j) ]  += weight_iqp * (-1.) * laplace_mat_dadj_u;

	      
	      //============= delta_mu row ===============================
	        if ( i==j )   
		  Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mat_mu, pos_mat_mu, i, j) ]  += penalty_outside_control_domain_boundary * ( (1 - control_node_flag[first_loc_comp_ctrl][i]));    //MU
          
	         } // end phi_j loop
           } // endif assemble_matrix
//============ Volume Jacobians - END ==================	    

        } // end phi_i loop
        
      } // end gauss point loop

  

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked(Res, L2G_dofmap_Mat_AllVars);

    if (assembleMatrix) {
      JAC->add_matrix_blocked(Jac, L2G_dofmap_Mat_AllVars, L2G_dofmap_Mat_AllVars);
    }
    
    
    //========== dof-based part, without summation
 
     if (print_algebra_local) {
         assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs_Mat_vol, 10, 5);
     }
     
     
  } //end element loop for each process
  


  //MU in res ctrl - BEGIN  ***********************************
femus::ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::add_one_times_mu_res_ctrl(iproc,
                               ineq_flag,
                               ctrl_index_in_mat,
                               mu_index_in_mat,
                               SolIndex_Mat,
                               sol,
                               mlPdeSys,
                               pdeSys,
                               RES);
  //MU in res ctrl - END ***********************************

    
  // ***************** END ASSEMBLY - ADD PART *******************

RES->close();
if (assembleMatrix) JAC->close();  /// This is needed for the parallel, when splitting the add part from the insert part!!!


//   ***************** INSERT PART - BEGIN (must go AFTER the sum, clearly) *******************
 /// @todo One very important thing to consider: we have some PENALTIES that were set before during the SUMMATION part.
 // Now, if we do INSERT, we may end up OVERWRITING certain values, SUCH AS THOSE PENALTIES!!!
 // So you have to be very careful here!
    
     //MU

   for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {
       
// -------
   geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
      
   geom_element_iel.set_elem_center_3d(iel, solType_coords);
// -------
   
// -------
    el_dofs_unknowns_vol(sol, msh, pdeSys, iel,
                         SolFEType_Mat,
                         SolIndex_Mat,
                         SolPdeIndex,
                         Sol_n_el_dofs_Mat_vol, 
                         Sol_eldofs_Mat, 
                         L2G_dofmap_Mat);
// -------

	if ( ctrl::square_or_cube :: pure_boundary< ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >  ::volume_elem_contains_a_Gamma_control_face( sol, msh, iel ) ) {


    	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);


       std::pair< int, unsigned int > pair_control_iface = femus::face_is_a_Gamma_control_face_of_some_index< femus::ctrl:: GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >(msh->GetMeshElements(), iel, iface);

        const int  iface_is_a_boundary_control  = pair_control_iface.first;

	    if( iface_is_a_boundary_control ) {

       femus::ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::update_active_set_flag_for_current_nonlinear_iteration_bdry
   (msh, sol,
    iel, iface,
    geom_element_iel.get_coords_at_dofs_bdry_3d(), 
    Sol_eldofs_Mat, 
    Sol_n_el_dofs_Mat_vol, 
    c_compl,
    mu_index_in_mat,
    ctrl_index_in_mat,
    solIndex_act_flag_sol,
    ctrl_lower,
    ctrl_upper,
    sol_actflag);
 

  femus::ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::node_insertion_bdry(iel, iface, 
                      msh,
                      L2G_dofmap_Mat,
                      mu_index_in_mat, 
                      ctrl_index_in_mat,
                      Sol_eldofs_Mat,
                      Sol_n_el_dofs_Mat_vol,
                      sol_actflag, 
                      ctrl_lower, ctrl_upper,
                      ineq_flag,
                      c_compl,
                      JAC, 
                      RES,
                      assembleMatrix
                      );
  
             }
             
       }
     }


     //============= delta_ctrl-delta_mu row ===============================
 if (assembleMatrix) {
       for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
         JAC->matrix_set_off_diagonal_values_blocked(L2G_dofmap_Mat[ctrl_index_in_mat[kdim]],  L2G_dofmap_Mat[mu_index_in_mat[kdim]], ineq_flag * 1.);
       }
}

   }
   
   
  // ***************** INSERT PART - END *******************

  RES->close();
  if (assembleMatrix) JAC->close();  ///@todo is it needed? I think so
   
  
  print_global_residual_jacobian(print_algebra_global,
                                 ml_prob,
                                 mlPdeSys,
                                 pdeSys,
                                 RES,
                                 JAC,
                                 iproc,
                                 assembleMatrix);
  

  return;
}

};



class lifting_internal : public femus::lifting_internal   {



public:
    
 static const std::vector< Unknown >  provide_list_of_unknowns(const unsigned int dimension) {
     
     
  std::vector< FEFamily > feFamily;
  std::vector< FEOrder >   feOrder;

                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
 
                        feOrder.push_back(/*FIRST*/SECOND);  //same
                        feOrder.push_back(/*FIRST*/SECOND);  //same
                        feOrder.push_back(/*FIRST*/SECOND);
                        feOrder.push_back(/*FIRST*/SECOND);  //same
 

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Unknown >  unknowns(feFamily.size());

   unknowns[0]._name      = "state";
   unknowns[1]._name      = "control";
   unknowns[2]._name      = "adjoint";
   unknowns[3]._name      = "mu";

   unknowns[0]._is_sparse = true;
   unknowns[1]._is_sparse = true;
   unknowns[2]._is_sparse = true;
   unknowns[3]._is_sparse = true;
   
     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._fe_family  = feFamily[u];
              unknowns[u]._fe_order   = feOrder[u];
              unknowns[u]._time_order = 0;
              unknowns[u]._is_pde_unknown = true;
              
     }
 
 
   return unknowns;
     
}


    
    

    
 static void assemble_elliptic_dirichlet_control(MultiLevelProblem& ml_prob) {

  NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             JAC = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id();

  constexpr bool print_algebra_global = false;
  constexpr bool print_algebra_local = false;
  
  //=============== Geometry ========================================
  unsigned solType_coords = CONTINUOUS_BIQUADRATIC;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
   //***************************************************  
    

 //********************* state *********************** 
 //***************************************************  
  std::vector <double> phi_u; 
  std::vector <double> phi_u_x;

  phi_u.reserve(max_size);
  phi_u_x.reserve(max_size * dim);
  
 
  unsigned solIndex_u =  solIndex_u = ml_sol->GetIndex("state");
  unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);

  unsigned solPdeIndex_u  =  solPdeIndex_u = mlPdeSys->GetSolPdeIndex("state");

  std::vector < double >  sol_u;
  sol_u.reserve(max_size);
  std::vector < int > l2GMap_u;
  l2GMap_u.reserve(max_size);
 //***************************************************  
 //***************************************************  

  
 //******************** control ********************** 
 //***************************************************   
  std::vector <double> phi_ctrl;
  std::vector <double> phi_ctrl_x;

  phi_ctrl.reserve(max_size);
  phi_ctrl_x.reserve(max_size * dim);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = ml_sol->GetIndex("control");
  unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

  unsigned solPdeIndex_ctrl;
  solPdeIndex_ctrl = mlPdeSys->GetSolPdeIndex("control");

  std::vector < double >  sol_ctrl;
  sol_ctrl.reserve(max_size);
  std::vector < int > l2GMap_ctrl;
  l2GMap_ctrl.reserve(max_size);
 //***************************************************  
 //***************************************************  
  
  
 //********************* adjoint ********************* 
 //***************************************************  
  std::vector <double> phi_adj;
  std::vector <double> phi_adj_x;

  phi_adj.reserve(max_size);
  phi_adj_x.reserve(max_size * dim);
 
  
  unsigned solIndex_adj;
  solIndex_adj = ml_sol->GetIndex("adjoint");
  unsigned solType_adj = ml_sol->GetSolutionType(solIndex_adj);

  unsigned solPdeIndex_adj;
  solPdeIndex_adj = mlPdeSys->GetSolPdeIndex("adjoint");

  std::vector < double >  sol_adj;
    sol_adj.reserve(max_size);
  std::vector < int > l2GMap_adj;
    l2GMap_adj.reserve(max_size);
 //***************************************************  
 //***************************************************  

 //****************** mu ******************************  
 //***************************************************  
  std::vector <double> phi_mu;
  std::vector <double> phi_mu_x;

  phi_mu.reserve(max_size);
  phi_mu_x.reserve(max_size * dim);

  
  unsigned solIndex_mu;
  solIndex_mu = ml_sol->GetIndex("mu");
   
  unsigned solPdeIndex_mu;
  solPdeIndex_mu = mlPdeSys->GetSolPdeIndex("mu");
  
  unsigned solType_mu = ml_sol->GetSolutionType(solIndex_mu);
  std::vector < double >  sol_mu;   sol_mu.reserve(max_size);
  std::vector < int > l2GMap_mu;   l2GMap_mu.reserve(max_size);

  
    const unsigned int n_components_ctrl = 1;
    const unsigned int first_loc_comp_ctrl = 0;

  //************** variables for ineq constraints: act flag ****************************   
    std::vector <unsigned int> solIndex_act_flag_sol(n_components_ctrl);
  
  femus::ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::store_act_flag_in_old(mlPdeSys, ml_sol, sol, solIndex_act_flag_sol);
    
  
  
  //********* variables for ineq constraints *****************
     std::vector <std::vector < double/*int*/ > > sol_actflag(n_components_ctrl);    //flag for active set
     std::vector <std::vector < double > > ctrl_lower(n_components_ctrl);  
     std::vector <std::vector < double > > ctrl_upper(n_components_ctrl);  
  
     for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {
          sol_actflag[kdim].reserve(max_size);
           ctrl_lower[kdim].reserve(max_size);
           ctrl_upper[kdim].reserve(max_size);
     }
     
  const int ineq_flag = INEQ_FLAG;
  const double c_compl = C_COMPL;
  //***************************************************  

 //***************************************************  
 //********* WHOLE SET OF VARIABLES ****************** 

  const int n_unknowns =  mlPdeSys->GetSolPdeIndex().size();
 
  std::vector < int > l2GMap_AllVars;
  l2GMap_AllVars.reserve(n_unknowns*max_size);
  
  std::vector < double > Res;
  Res.reserve(n_unknowns*max_size);

  std::vector < double > Jac;
  Jac.reserve( n_unknowns*max_size * n_unknowns*max_size);
  

//***************************************************
    enum Pos_in_matrix {pos_mat_state = 0, pos_mat_ctrl, pos_mat_adj, pos_mat_mu}; //these are known at compile-time 
                    ///@todo these are the positions in the MlSol object or in the Matrix? I'd say the matrix, but we have to check where we use it...

                    
    assert(pos_mat_state   == mlPdeSys->GetSolPdeIndex("state"));
    assert(pos_mat_ctrl    == mlPdeSys->GetSolPdeIndex("control"));
    assert(pos_mat_adj     == mlPdeSys->GetSolPdeIndex("adjoint"));
    assert(pos_mat_mu      == mlPdeSys->GetSolPdeIndex("mu"));

 std::vector < unsigned >   pos_ctrl_in_mat(n_components_ctrl);  pos_ctrl_in_mat[0] = pos_mat_ctrl;
 std::vector < unsigned >   pos_mu_in_mat(n_components_ctrl);    pos_mu_in_mat[0] = pos_mat_mu;
    
    //***************************************************

  std::vector < std::string > Solname(n_unknowns);
  Solname[0] = "state";
  Solname[1] = "control";
  Solname[2] = "adjoint";
  Solname[3] = "mu";
  

 //***************************************************
   enum Pos_in_Sol {pos_sol_state = 0, pos_sol_ctrl, pos_sol_adj, pos_sol_mu, pos_sol_targreg, pos_sol_contreg, pos_sol_actflag}; //these are known at compile-time 

        assert(pos_sol_actflag == solIndex_act_flag_sol[0]);
//***************************************************
    
  


  std::vector < unsigned > SolPdeIndex(n_unknowns);
  std::vector < unsigned > SolIndex(n_unknowns);  
  std::vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(Solname[ivar].c_str());
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

    std::vector < unsigned > Sol_n_el_dofs_Mat_vol(n_unknowns);
 //***************************************************  

//***************************************************
    std::vector < std::vector < double > >  sol_eldofs_Mat(n_unknowns);  //should have Mat order
    for(int k = 0; k < n_unknowns; k++) {        sol_eldofs_Mat[k].reserve(max_size);    }


    //----------- quantities (at dof objects) ------------------------------
    std::vector< int >       L2G_dofmap_Mat_AllVars;
    L2G_dofmap_Mat_AllVars.reserve( n_unknowns * max_size );
    std::vector < std::vector < int > >     L2G_dofmap_Mat(n_unknowns);     //should have Mat order
    for(int i = 0; i < n_unknowns; i++) {
        L2G_dofmap_Mat[i].reserve(max_size);
    }

    
    
 //********************* DATA ************************ 
  double u_des = femus::ctrl::square_or_cube :: cost_functional_without_regularization::DesiredTargetVec()[0];
  double alpha = ALPHA_CTRL_VOL;
  double beta  = BETA_CTRL_VOL;
  double penalty_outside_control_domain = _lifting_internal_penalty_outside_control_domain;         // penalty for zero control outside
 //***************************************************  

 
 
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

  
 
      // ***************** ADD PART - BEGIN *******************
  RES->zero();
  if (assembleMatrix)  JAC->zero();

    

  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element_iel.geom_type();
    
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

  //************* set target domain flag **************
   int target_flag = 0;
   target_flag = femus::ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
 //*************************************************** 
   
    
 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);
    sol_u    .resize(nDof_u);
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
      l2GMap_u[i] = pdeSys->GetSystemDof(solIndex_u, solPdeIndex_u, i, iel);
    }
 //***************************************************  
 
 //***************** control ************************* 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);
    sol_ctrl    .resize(nDof_ctrl);
    l2GMap_ctrl.resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);
      l2GMap_ctrl[i] = pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);
    } 
 //*************************************************** 
 

 //************** adjoint **************************** 
    unsigned nDof_adj  = msh->GetElementDofNumber(iel, solType_adj);
        sol_adj    .resize(nDof_adj);
        l2GMap_adj.resize(nDof_adj);
    for (unsigned i = 0; i < sol_adj.size(); i++) {
      unsigned solDof_adj = msh->GetSolutionDof(i, iel, solType_adj);
      sol_adj[i] = (*sol->_Sol[solIndex_adj])(solDof_adj);
      l2GMap_adj[i] = pdeSys->GetSystemDof(solIndex_adj, solPdeIndex_adj, i, iel);
    } 
 //***************************************************  
 
 //************** mu **************************** 
    unsigned nDof_mu  = msh->GetElementDofNumber(iel, solType_mu);
    sol_mu   .resize(nDof_mu);
    l2GMap_mu.resize(nDof_mu);
    for (unsigned i = 0; i < sol_mu.size(); i++) {
      unsigned solDof_mu = msh->GetSolutionDof(i, iel, solType_mu);
      sol_mu[i] = (*sol->_Sol[solIndex_mu])(solDof_mu);
      l2GMap_mu[i] = pdeSys->GetSystemDof(solIndex_mu, solPdeIndex_mu, i, iel);
    }
 //***************************************************  
    
 
    
 //******************** ALL VARS ********************* 
    unsigned nDof_AllVars = nDof_u + nDof_ctrl + nDof_adj + nDof_mu; 
    int nDof_max    =  nDof_u;   // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    if(nDof_adj > nDof_max) 
    {
      nDof_max = nDof_adj;
      }
    
    if(nDof_ctrl > nDof_max)
    {
      nDof_max = nDof_ctrl;
    }
    
    
    Res.resize(nDof_AllVars);
    std::fill(Res.begin(), Res.end(), 0.);

    Jac.resize(nDof_AllVars * nDof_AllVars);
    std::fill(Jac.begin(), Jac.end(), 0.);
    
    l2GMap_AllVars.resize(0);
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_u.begin(),l2GMap_u.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_ctrl.begin(),l2GMap_ctrl.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_adj.begin(),l2GMap_adj.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_mu.begin(),l2GMap_mu.end());
    
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	Sol_n_el_dofs_Mat_vol[k] = ndofs_unk;
  }
    
 //*************************************************** 
    
 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = femus::ctrl:: square_or_cube:: lifting_internal< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES>::ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
  std::vector<int> control_node_flag(nDof_ctrl, 0);
  if (control_el_flag == 1) std::fill(control_node_flag.begin(), control_node_flag.end(), 1);
 //*************************************************** 
  
 //========= gauss value quantities ==================   
	double sol_u_gss = 0.;
	double sol_adj_gss = 0.;
	double sol_ctrl_gss = 0.;
	double sol_mu_gss = 0.;
	std::vector<double> sol_u_x_gss(dim);       std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	std::vector<double> sol_adj_x_gss(dim);     std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	std::vector<double> sol_ctrl_x_gss(dim);    std::fill(sol_ctrl_x_gss.begin(), sol_ctrl_x_gss.end(), 0.);
	std::vector<double> sol_mu_x_gss(dim);      std::fill(sol_mu_x_gss.begin(), sol_mu_x_gss.end(), 0.);
 //===================================================   

      // *** Gauss point loop ***
      for (unsigned iqp = 0; iqp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); iqp++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), iqp, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);

    AbsDetJxWeight_iqp = detJac_iqp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[iqp];

    
    elem_all[qrule_i][ielGeom][solType_u]  ->shape_funcs_current_elem(iqp, JacI_iqp, phi_u, phi_u_x, boost::none, space_dim);
    elem_all[qrule_i][ielGeom][solType_ctrl]  ->shape_funcs_current_elem(iqp, JacI_iqp, phi_ctrl, phi_ctrl_x, boost::none, space_dim);
    elem_all[qrule_i][ielGeom][solType_adj]  ->shape_funcs_current_elem(iqp, JacI_iqp, phi_adj, phi_adj_x, boost::none, space_dim);
    elem_all[qrule_i][ielGeom][solType_mu]  ->shape_funcs_current_elem(iqp, JacI_iqp, phi_mu, phi_mu_x, boost::none, space_dim);

   
//     // 	msh->_finiteElement[ielGeom][solType_u]   ->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), iqp, AbsDetJxWeight_iqp, phi_u, phi_u_x, boost::none);
//     msh->_finiteElement[ielGeom][solType_ctrl]->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), iqp, AbsDetJxWeight_iqp, phi_ctrl, phi_ctrl_x, boost::none);
//     msh->_finiteElement[ielGeom][solType_adj] ->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), iqp, AbsDetJxWeight_iqp, phi_adj, phi_adj_x, boost::none);
// 	msh->_finiteElement[ielGeom][solType_mu]  ->Jacobian(geom_element_iel.get_coords_at_dofs_3d(), iqp, AbsDetJxWeight_iqp, phi_mu, phi_mu_x, boost::none);

	
    sol_u_gss = 0.;
	sol_adj_gss = 0.;
	sol_ctrl_gss = 0.;
	sol_mu_gss = 0.;
    
	std::fill(sol_u_x_gss.begin(),sol_u_x_gss.end(), 0.);
	std::fill(sol_adj_x_gss.begin(), sol_adj_x_gss.end(), 0.);
	std::fill(sol_ctrl_x_gss.begin(), sol_ctrl_x_gss.end(), 0.);
	std::fill(sol_mu_x_gss.begin(), sol_mu_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < nDof_u; i++) {
	                                                sol_u_gss      += sol_u[i] * phi_u[i];
                   for (unsigned d = 0; d < dim; d++)   sol_u_x_gss[d] += sol_u[i] * phi_u_x[i * space_dim + d];
          }
	
	for (unsigned i = 0; i < nDof_adj; i++) {
	                                                sol_adj_gss      += sol_adj[i] * phi_adj[i];
                   for (unsigned d = 0; d < dim; d++)   sol_adj_x_gss[d] += sol_adj[i] * phi_adj_x[i * space_dim + d];
        }
	
	for (unsigned i = 0; i < nDof_ctrl; i++) {
	                                                sol_ctrl_gss      += sol_ctrl[i] * phi_ctrl[i];
                   for (unsigned d = 0; d < dim; d++)   sol_ctrl_x_gss[d] += sol_ctrl[i] * phi_ctrl_x[i * space_dim + d];
        }
        
    for (unsigned i = 0; i < nDof_mu; i++) {
	                                                sol_mu_gss      += sol_mu[i] * phi_mu[i];
		   for (unsigned d = 0; d < dim; d++)   sol_mu_x_gss[d] += sol_mu[i] * phi_mu_x[i * space_dim + d];
	}
        
//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
	      double laplace_rhs_du_adj_i = 0.; //
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_u )         laplace_rhs_du_adj_i             +=  (phi_u_x   [i * space_dim + kdim] * sol_adj_x_gss[kdim]);
	      }
	      
              double laplace_rhs_dctrl_ctrl_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_ctrl )         laplace_rhs_dctrl_ctrl_i      +=  (phi_ctrl_x   [i * space_dim + kdim] * sol_ctrl_x_gss[kdim]);
	      }
	      
	      double laplace_rhs_dctrl_adj_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_ctrl )         laplace_rhs_dctrl_adj_i       +=  (phi_ctrl_x   [i * space_dim + kdim] * sol_adj_x_gss[kdim]);
	      }
	      
	      double laplace_rhs_dadj_u_i = 0.;  //
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj )         laplace_rhs_dadj_u_i           +=  (phi_adj_x   [i * space_dim + kdim] * sol_u_x_gss[kdim]);
	      }
	      
	      double laplace_rhs_dadj_ctrl_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj )         laplace_rhs_dadj_ctrl_i        +=  (phi_adj_x   [i * space_dim + kdim] * sol_ctrl_x_gss[kdim]);
	      }
	      
              double laplace_rhs_du_u_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
                 if ( i < nDof_u )       laplace_rhs_du_u_i            +=  (phi_u_x  [i * space_dim + kdim] *sol_u_x_gss[kdim]);
	      }
	      
              double laplace_rhs_du_ctrl_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_u )          laplace_rhs_du_ctrl_i      +=  (phi_u_x   [i * space_dim + kdim] * sol_ctrl_x_gss[kdim]);
	      }
	      
              double laplace_rhs_dctrl_u_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_ctrl )         laplace_rhs_dctrl_u_i      +=  (phi_ctrl_x   [i * space_dim + kdim] * sol_u_x_gss[kdim]);
	      }
	      
//======================Residuals - BEGIN =======================
          // FIRST ROW
	  if (i < nDof_u)                      Res[0      + i] += - AbsDetJxWeight_iqp * (
          - laplace_rhs_du_adj_i 
     + target_flag * 
#if COST_FUNCTIONAL_TYPE == 0
   phi_u[i] * ( sol_u_gss + sol_ctrl_gss - u_des) 
#elif COST_FUNCTIONAL_TYPE == 1
    ( laplace_rhs_du_u_i + laplace_rhs_du_ctrl_i )
#endif
          );
          // SECOND ROW
	  if (i < nDof_ctrl)  {
	     if ( control_el_flag == 1)    {    Res[nDof_u + i] +=  /*(control_node_flag[i]) **/ - AbsDetJxWeight_iqp * ( 
													                                                  + OP_L2 * alpha * phi_ctrl[i] * sol_ctrl_gss
		                                                                                              + OP_H1 * beta * laplace_rhs_dctrl_ctrl_i
		                                                                                              - laplace_rhs_dctrl_adj_i /// @todo this is here because variations of restricted functions are also restricted
													                                                  /*+ ineq_flag * sol_mu_gss*/ 
                                                                                                      + target_flag *
                                                                                 #if COST_FUNCTIONAL_TYPE == 0
                                                                                                      phi_ctrl[i] * ( sol_u_gss + sol_ctrl_gss - u_des)
                                                                                  #elif COST_FUNCTIONAL_TYPE == 1
                                                                                                      ( laplace_rhs_dctrl_u_i + laplace_rhs_dctrl_ctrl_i )
                                                                                  #endif
                                                                                                                    ); 
             
        }
	      else if ( control_el_flag == 0) { Res[nDof_u + i] +=   (- penalty_outside_control_domain)  *  (1 - control_node_flag[i]) * (sol_ctrl[i] - 0.); }
	  }
          // THIRD ROW
          if (i < nDof_adj)        Res[nDof_u + nDof_ctrl + i] += /*-AbsDetJxWeight_iqp * phi_adj[i] * sol_adj_gss - 6.;*/- AbsDetJxWeight_iqp *  ( - laplace_rhs_dadj_u_i - laplace_rhs_dadj_ctrl_i ) ;

          // FOURTH ROW
       if (i < nDof_mu)  {
            if ( control_el_flag == 0) {  Res[nDof_u + nDof_ctrl + nDof_adj + i] +=   (- penalty_outside_control_domain) *  (1 - control_node_flag[i]) * (sol_mu[i] - 0.); }   //MU
       }
       
           
//====================== Residuals - END =======================
	      
//====================== Jacobians - BEGIN =======================
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
              //double laplace_mat_du_u = 0.;
              double laplace_mat_du_adj = 0.;
              double laplace_mat_dadj_u = 0.;
              double laplace_mat_dctrl_adj = 0.;
              double laplace_mat_dadj_ctrl = 0.;
              double laplace_mat_dctrl_ctrl = 0.;
              double laplace_mat_du_u = 0.;
              double laplace_mat_du_ctrl = 0.;
              double laplace_mat_dctrl_u = 0.;

              for (unsigned kdim = 0; kdim < dim; kdim++) {
              //if ( i < nDof_u && j < nDof_u )           laplace_mat_du_u           += (phi_u_x   [i * space_dim + kdim] * phi_u_x   [j * space_dim + kdim]);
              if ( i < nDof_u && j < nDof_adj )         laplace_mat_du_adj         += (phi_u_x   [i * space_dim + kdim] * phi_adj_x [j * space_dim + kdim]);
              if ( i < nDof_adj && j < nDof_u )         laplace_mat_dadj_u         += (phi_adj_x [i * space_dim + kdim] * phi_u_x   [j * space_dim + kdim]);  //equal to the previous
              if ( i < nDof_ctrl && j < nDof_adj )      laplace_mat_dctrl_adj      += (phi_ctrl_x[i * space_dim + kdim] * phi_adj_x [j * space_dim + kdim]);
              if ( i < nDof_adj && j < nDof_ctrl )      laplace_mat_dadj_ctrl      += (phi_adj_x [i * space_dim + kdim] * phi_ctrl_x[j * space_dim + kdim]);  //equal to the previous
              if ( i < nDof_ctrl && j < nDof_ctrl )     laplace_mat_dctrl_ctrl     += (phi_ctrl_x  [i * space_dim + kdim] * phi_ctrl_x  [j * space_dim + kdim]);
              
              if ( i < nDof_u && j < nDof_u )     laplace_mat_du_u           += (phi_u_x  [i * space_dim + kdim] * phi_u_x  [j * space_dim + kdim]);
              if ( i < nDof_u && j < nDof_ctrl )     laplace_mat_du_ctrl     += (phi_u_x  [i * space_dim + kdim] * phi_ctrl_x  [j * space_dim + kdim]);
              if ( i < nDof_ctrl && j < nDof_u )     laplace_mat_dctrl_u     += (phi_ctrl_x  [i * space_dim + kdim] * phi_u_x  [j * space_dim + kdim]);
	      }

              //============ delta_state row ============================
              //DIAG BLOCK delta_state - state
	      if ( i < nDof_u && j < nDof_u )       
		Jac[ (0 + i) * nDof_AllVars   +
		     (0 + j)                            ]  += AbsDetJxWeight_iqp * target_flag *
#if COST_FUNCTIONAL_TYPE == 0
		     phi_u[j] *  phi_u[i]
#elif COST_FUNCTIONAL_TYPE == 1
             laplace_mat_du_u 
#endif
              ;
              
	      // BLOCK  delta_state - control
              if ( i < nDof_u && j < nDof_ctrl )   
		Jac[ (0 + i) * nDof_AllVars   +
                     (nDof_u + j)                       ]  += AbsDetJxWeight_iqp * target_flag  * 
#if COST_FUNCTIONAL_TYPE == 0
                     phi_ctrl[j] *  phi_u[i]
#elif COST_FUNCTIONAL_TYPE == 1
                     laplace_mat_du_ctrl
#endif
	      ;
          
              // BLOCK  delta_state - adjoint
              if ( i < nDof_u && j < nDof_adj )  
		Jac[  (0 + i) * nDof_AllVars  +
                      (nDof_u + nDof_ctrl + j)          ]  += AbsDetJxWeight_iqp * (-1.) * laplace_mat_du_adj;
              
	      //=========== delta_control row ===========================     
	      if ( control_el_flag == 1)  {

	      //BLOCK delta_control - state
              if ( i < nDof_ctrl   && j < nDof_u   ) 
		Jac[ (nDof_u + i) * nDof_AllVars  +
		     (0 + j)                            ]  += 
		     ( control_node_flag[i]) * AbsDetJxWeight_iqp * target_flag * 
#if COST_FUNCTIONAL_TYPE == 0
		     phi_u[j] * phi_ctrl[i]
#elif COST_FUNCTIONAL_TYPE == 1
             laplace_mat_dctrl_u
#endif
             ;
             
	      //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl   )
		Jac[ (nDof_u + i) * nDof_AllVars +
		     (nDof_u + j)                       ]  += ( control_node_flag[i]) * AbsDetJxWeight_iqp * ( OP_H1 * beta * /*control_el_flag  **/ laplace_mat_dctrl_ctrl
		                                                                                + OP_L2 * alpha * /*control_el_flag **/ phi_ctrl[i] * phi_ctrl[j]
		                                                                                            + target_flag *
                                                                            #if COST_FUNCTIONAL_TYPE == 0
		                                                                                            phi_ctrl[i] * phi_ctrl[j] 
                                                                            #elif COST_FUNCTIONAL_TYPE == 1
                                                                                               laplace_mat_dctrl_ctrl
                                                                            #endif
                                                                                                     );
              
	      //BLOCK delta_control - adjoint
              if ( i < nDof_ctrl   && j < nDof_adj  ) 
		Jac[ (nDof_u + i) * nDof_AllVars  + 
		     (nDof_u + nDof_ctrl + j)           ]  += ( control_node_flag[i]) * AbsDetJxWeight_iqp * (-1.) * laplace_mat_dctrl_adj;
	      	      
	        }
	      
	      else if ( control_el_flag == 0)  {  
		
              //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl   &&  i==j ) {
		 Jac[ (nDof_u + i) * nDof_AllVars +
		      (nDof_u + j)                      ]  +=  penalty_outside_control_domain * (1 - control_node_flag[i]);
		}
	      
	      }
	      
		     
	      //=========== delta_adjoint row ===========================
              // BLOCK delta_adjoint - state	      
              if ( i < nDof_adj && j < nDof_u )   
		Jac[ (nDof_u + nDof_ctrl + i) * nDof_AllVars +
		     (0 + j)                            ]  += AbsDetJxWeight_iqp * (-1) * laplace_mat_dadj_u;   
	      
              // BLOCK delta_adjoint - control   
              if ( i < nDof_adj && j < nDof_ctrl )  
		Jac[ (nDof_u + nDof_ctrl + i)  * nDof_AllVars +
		     (nDof_u  + j)                      ]  += AbsDetJxWeight_iqp * (-1) * laplace_mat_dadj_ctrl; 
		     
		     
// 	      // BLOCK delta_adjoint - adjoint   
//               if ( i < nDof_adj && j < nDof_adj )  
// 		Jac[ (nDof_u + nDof_ctrl + i)  * nDof_AllVars +
// 		     (nDof_u + nDof_ctrl + j)                      ]  += AbsDetJxWeight_iqp * phi_adj[j] *  phi_adj[i]; 
    
	      
	      //============= delta_mu row ===============================
                if ( i < nDof_mu && j < nDof_mu && i==j )  {   
	        if ( control_el_flag == 0)  {  
		  Jac[   (nDof_u + nDof_ctrl + nDof_adj + i) * nDof_AllVars +  (nDof_u + nDof_ctrl + nDof_adj + j) ]  += penalty_outside_control_domain * (1 - control_node_flag[i]);    //MU
                }
	      }
	      
	      
            } // end phi_j loop
          } // endif assemble_matrix
//====================== Jacobians - END =======================

        } // end phi_i loop
        
      } // end gauss point loop

      
      
    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
   // std::cout << " ************* Element ************** " << iel << " **************************************** " << std::endl;     

// // //     if (control_el_flag == 0) {  //elements that should have zero control
// //          for (unsigned i_unk = 0; i_unk < n_unknowns; i_unk++) {
// //     std::cout << " ======= Row === " << i_unk << " =================================================== " << std::endl;     
// //         unsigned int row_block_offset = 0;
// // 	         for (unsigned k = 0; k < i_unk; k++) row_block_offset += Sol_n_el_dofs_Mat_vol[k];
// //          for (unsigned j_unk = 0; j_unk < n_unknowns; j_unk++) {
// //     std::cout << " ======= Column === " << j_unk << " ================== " << std::endl;     
// //         unsigned int column_block_offset = 0;
// // 	         for (unsigned k = 0; k < j_unk; k++) column_block_offset += Sol_n_el_dofs_Mat_vol[k];
// // 	  
// //          for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[i_unk]; i++) {
// // // 	      std::cout << Res[nDof_u + nDof_ctrl + nDof_adj + i ] << " " << std::endl;
// // 	   for (unsigned j = 0; j < Sol_n_el_dofs_Mat_vol[j_unk]; j++) {
// // 	      std::cout <<  " " << std::setfill(' ') << std::setw(10) << Jac[ (row_block_offset + i) * nDof_AllVars + ( column_block_offset + j) ] << " ";
// // 	    }
// // 	      std::cout << std::endl;
// // 	 }
// // 
// // 	 } //j_unk
// // 	} //i_unk
	 
	 
// // // 	}


 //========== end of integral-based part

    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap_AllVars);
      if (assembleMatrix)  JAC->add_matrix_blocked(Jac, l2GMap_AllVars, l2GMap_AllVars);

      
     if (print_algebra_local) {
         assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs_Mat_vol, 10, 5);
     }

       
  } //end element loop for each process
      
  //MU
  
  
// // // 
// // //   unsigned int global_ctrl_size = pdeSys->KKoffset[ctrl_index+1][iproc] - pdeSys->KKoffset[ctrl_index][iproc];
// // //   
// // //   std::vector<double>  one_times_mu(global_ctrl_size, 0.);
// // //   std::vector<int>    positions(global_ctrl_size);
// // // //  double position_mu_i;
// // //   for (unsigned i = 0; i < positions.size(); i++) {
// // //     positions[i] = pdeSys->KKoffset[ctrl_index][iproc] + i;
// // // //     position_mu_i = pdeSys->KKoffset[mu_index][iproc] + i;
// // // //     std::cout << position_mu_i << std::endl;
// // //     one_times_mu[i] = ineq_flag * 1. * (*sol->_Sol[solIndex_mu])(i/*position_mu_i*/) ;
// // //   }
// // //     RES->add_vector_blocked(one_times_mu, positions);
  std::vector<unsigned int> ctrl_index(1);  ctrl_index[0] = mlPdeSys->GetSolPdeIndex("control");
  std::vector<unsigned int>   mu_index(1);    mu_index[0] = mlPdeSys->GetSolPdeIndex("mu");
    
ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::add_one_times_mu_res_ctrl(iproc,
                               ineq_flag,
                               ctrl_index,
                               mu_index,
                               SolIndex,
                               sol,
                               mlPdeSys,
                               pdeSys,
                               RES);
    
    
    
     
RES->close();
if (assembleMatrix) JAC->close();  /// This is needed for the parallel, when splitting the add part from the insert part!!!
      // ***************** ADD PART - END  *******************
    
    
 


//   ***************** INSERT PART - BEGIN (must go AFTER the sum, clearly) *******************
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

// -------
   geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
      
   geom_element_iel.set_elem_center_3d(iel, solType_coords);
// -------
   
// -------
    el_dofs_unknowns_vol(sol, msh, pdeSys, iel,
                        SolFEType,
                        SolIndex,
                        SolPdeIndex,
                        Sol_n_el_dofs_Mat_vol,
                        sol_eldofs_Mat,
                        L2G_dofmap_Mat);
      
    //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ctrl:: square_or_cube:: lifting_internal< ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
 
    
    if (control_el_flag == 1) {

        
  ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::update_active_set_flag_for_current_nonlinear_iteration
  (msh,
   sol,
   iel,
   geom_element_iel.get_coords_at_dofs/*_3d*/(),
   sol_eldofs_Mat,
   Sol_n_el_dofs_Mat_vol,
   c_compl,
   pos_mu_in_mat,
   pos_ctrl_in_mat,
   solIndex_act_flag_sol,
   ctrl_lower,
   ctrl_upper,
   sol_actflag);
  
      


    ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::node_insertion(iel,
                   msh,
                   L2G_dofmap_Mat,
                   pos_mu_in_mat,
                   pos_ctrl_in_mat,
                   sol_eldofs_Mat,
                   Sol_n_el_dofs_Mat_vol,
                   sol_actflag,
                   ctrl_lower,
                   ctrl_upper,
                   ineq_flag,
                   c_compl,
                   JAC,
                   RES,
                   assembleMatrix
                   );
   

       }

 
    
    
 //============= delta_ctrl-delta_mu row ===============================
  if (assembleMatrix) { JAC->matrix_set_off_diagonal_values_blocked(L2G_dofmap_Mat[pos_mat_ctrl], L2G_dofmap_Mat[pos_mat_mu], ineq_flag * 1.); }
  
  
  }
//   ***************** INSERT PART - END (must go AFTER the sum, clearly) *******************
  
  RES->close();

  if (assembleMatrix) JAC->close();
  
  // ***************** END ASSEMBLY *******************
  
  
  print_global_residual_jacobian(print_algebra_global,
                                 ml_prob,
                                 mlPdeSys,
                                 pdeSys,
                                 RES,
                                 JAC,
                                 iproc,
                                 assembleMatrix);


  return;
}


};




class lifting_external: public femus::lifting_external  {

public: 
    
 static void compute_coordinates_bdry_one_face(std::vector< std::vector <double> > & coords_at_dofs_bdry, const int solType_coords, const unsigned int iel, const int jface, const Mesh * msh)  {
    
     const unsigned int dim = coords_at_dofs_bdry.size();
     
            unsigned nDofx_bdry    = msh->GetElementFaceDofNumber(iel,jface,solType_coords);

            for (unsigned idim = 0; idim < dim; idim++) {
                coords_at_dofs_bdry[idim].resize(nDofx_bdry);
            }
            
            for(unsigned i=0; i < nDofx_bdry; i++) {
                unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
                unsigned iDof = msh->GetSolutionDof(i_vol, iel, solType_coords);
                for(unsigned idim=0; idim<dim; idim++) {
                    coords_at_dofs_bdry[idim][i]=(*msh->GetTopology()->_Sol[idim])(iDof);
                }
            }

}
 
 
// //============ find interface boundary elements (now we do with coordinates, later we can do also with flag) =======================================
 static bool find_control_boundary_nodes(std::vector<unsigned int> & interface_node_flag, const std::vector<double> & elem_center_bdry, const unsigned int nDofu_bdry, const unsigned int iel, const int jface, const Mesh * msh) {
   
       enum begin_end {begin = 0, end};
       
      std::vector< std::vector<double> > coords_bdry_control_face(2);
   
                    for (int be = 0; be < 2; be++)  {
                       coords_bdry_control_face[be].resize(elem_center_bdry.size());
                        std::fill(coords_bdry_control_face[be].begin(), coords_bdry_control_face[be].end(), 0.);
                    }
                    
                    coords_bdry_control_face[begin][0] = 1.;
                    coords_bdry_control_face[end][0] = 1.;
                    coords_bdry_control_face[begin][1] = SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MIN;
                    coords_bdry_control_face[end][1] = SQUARE_OR_CUBE__CONTROL_FACE_TANG_COORD_MAX;
                    
     
  bool interface_elem_flag = false;
            const double my_eps = 1.e-6;
            
//             if (elem_center_bdry[0] > 1. - my_eps   && elem_center_bdry[0] < 1.   + my_eps  &&
//                     elem_center_bdry[1] > 0. - my_eps && elem_center_bdry[1] < 1. + my_eps)
            if (elem_center_bdry[0] > coords_bdry_control_face[begin][0] - my_eps   && elem_center_bdry[0] < coords_bdry_control_face[end][0] + my_eps  &&
                elem_center_bdry[1] > coords_bdry_control_face[begin][1] - my_eps   && elem_center_bdry[1] < coords_bdry_control_face[end][1] + my_eps)
            {
                
                std::cout << " bdry elem on interface with center " << "(" << elem_center_bdry[0] << "," << elem_center_bdry[1] << ")" << std::endl;
                
             interface_elem_flag = true;
                
                for (int i_bdry = 0; i_bdry < nDofu_bdry; i_bdry++)  {
                    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
                    interface_node_flag[i_vol] = 1;
                   }


               }

      return interface_elem_flag;
      
}


static void assemble_elliptic_dirichlet_control(MultiLevelProblem& ml_prob) {
    //  ml_prob is the global object from/to where get/set all the data

    //  level is the level of the PDE system to be assembled
    //  levelMax is the Maximum level of the MultiLevelProblem
    //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

    //  extract pointers to the several objects that we are going to use

    NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
    const unsigned level = mlPdeSys->GetLevelToAssemble();
    const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

    Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);

    MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
    Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

    LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
    SparseMatrix*             JAC = pdeSys->_KK;
    NumericVector*           RES = pdeSys->_RES;

    const unsigned  dim = msh->GetDimension();
    const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

    const unsigned    iproc = msh->processor_id();
  
    
    constexpr bool print_algebra_global = true;
    constexpr bool print_algebra_local = false;
  


//***************************************************
  CurrentElem < double > geom_element_iel(dim, msh); 
  
    const int solType_coords = CONTINUOUS_BIQUADRATIC;  //biquadratic

//************** geometry (at dofs) *************************************
    std::vector < std::vector < double > > coords_at_dofs(dim);
    std::vector < std::vector < double > > coords_at_dofs_bdry(dim);
    for (unsigned idim = 0; idim < dim; idim++) {
        coords_at_dofs[idim].reserve(max_size);
        coords_at_dofs_bdry[idim].reserve(max_size);
    }

//************** geometry (at quadrature points) *************************************
    std::vector < double > coord_at_qp(dim);


//************* shape functions (at dofs and quadrature points) **************************************
    double weight_qp = 0.;      // gauss point weight
    double weight_qp_bdry = 0.; // gauss point weight on the boundary

    std::vector < std::vector < double > > phi_fe_qp(NFE_FAMS);
    std::vector < std::vector < double > > phi_x_fe_qp(NFE_FAMS);
    std::vector < std::vector < double > > phi_xx_fe_qp(NFE_FAMS);

    for(int fe=0; fe < NFE_FAMS; fe++) {
        phi_fe_qp[fe].reserve(max_size);
        phi_x_fe_qp[fe].reserve(max_size*dim);
        phi_xx_fe_qp[fe].reserve(max_size*(3*(dim-1)));
    }

//************* bdry shape functions (at dofs and quadrature points) **************************************
    std::vector < std::vector < double > > phi_fe_qp_bdry(NFE_FAMS);
    std::vector < std::vector < double > > phi_x_fe_qp_bdry(NFE_FAMS);

    for(int fe=0; fe < NFE_FAMS; fe++) {
        phi_fe_qp_bdry[fe].reserve(max_size);
        phi_x_fe_qp_bdry[fe].reserve(max_size * dim);
    }

//********************* vol-at-bdry adjoint *******************
    std::vector <double> phi_adj_vol_at_bdry;
    phi_adj_vol_at_bdry.reserve(max_size);
    std::vector <double> phi_adj_x_vol_at_bdry;
    phi_adj_x_vol_at_bdry.reserve(max_size * dim);
    std::vector <double> sol_adj_x_vol_at_bdry_gss(dim);
//***************************************************
    
//********************* vol-at-bdry adjoint_ext *******************
    std::vector <double> phi_adj_ext_vol_at_bdry;
    phi_adj_ext_vol_at_bdry.reserve(max_size);
    std::vector <double> phi_adj_ext_x_vol_at_bdry;
    phi_adj_ext_x_vol_at_bdry.reserve(max_size * dim);
    std::vector <double> sol_adj_ext_x_vol_at_bdry_gss(dim);
//***************************************************

    const unsigned int n_components_ctrl = 1;
    const unsigned int first_loc_comp_ctrl = 0;

  //************** act flag ****************************   
    std::vector <unsigned int> solIndex_act_flag_sol(n_components_ctrl);

  ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::store_act_flag_in_old(mlPdeSys, ml_sol, sol, solIndex_act_flag_sol);
    


    //********* variables for ineq constraints *****************
     std::vector <std::vector < double/*int*/ > > sol_actflag(n_components_ctrl);    //flag for active set
     std::vector <std::vector < double > > ctrl_lower(n_components_ctrl);  
     std::vector <std::vector < double > > ctrl_upper(n_components_ctrl);  
  
     for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {
          sol_actflag[kdim].reserve(max_size);
           ctrl_lower[kdim].reserve(max_size);
           ctrl_upper[kdim].reserve(max_size);
     }
     
    const int ineq_flag = INEQ_FLAG;
    const double c_compl = C_COMPL;
    //***************************************************

//***************************************************
//********* WHOLE SET OF VARIABLES ******************
//***************************************************
    enum Pos_in_matrix {pos_mat_state = 0, pos_mat_ctrl, pos_mat_adj, pos_mat_adj_ext, pos_mat_mu}; //these are known at compile-time 
                    ///@todo these are the positions in the MlSol object or in the Matrix? I'd say the matrix, but we have to check where we use it...

                    
    assert(pos_mat_state   == mlPdeSys->GetSolPdeIndex("state"));
    assert(pos_mat_ctrl    == mlPdeSys->GetSolPdeIndex("control"));
    assert(pos_mat_adj     == mlPdeSys->GetSolPdeIndex("adjoint"));
    assert(pos_mat_adj_ext == mlPdeSys->GetSolPdeIndex("adjoint_ext"));
    assert(pos_mat_mu      == mlPdeSys->GetSolPdeIndex("mu"));
//***************************************************


 std::vector < unsigned >   pos_ctrl_in_mat(n_components_ctrl);  pos_ctrl_in_mat[0] = pos_mat_ctrl;
 std::vector < unsigned >   pos_mu_in_mat(n_components_ctrl);    pos_mu_in_mat[0] = pos_mat_mu;
    
    
    
    const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();

    enum Sol_pos {pos_state=0, pos_ctrl, pos_adj, pos_adj_ext, pos_mu}; //these are known at compile-time

    assert(pos_state   == mlPdeSys->GetSolPdeIndex("state"));
    assert(pos_ctrl    == mlPdeSys->GetSolPdeIndex("control"));
    assert(pos_adj     == mlPdeSys->GetSolPdeIndex("adjoint"));
    assert(pos_adj_ext == mlPdeSys->GetSolPdeIndex("adjoint_ext"));
    assert(pos_mu      == mlPdeSys->GetSolPdeIndex("mu"));


    std::vector < std::string > Solname(n_unknowns);
    Solname[0] = "state";
    Solname[1] = "control";
    Solname[2] = "adjoint";
    Solname[3] = "adjoint_ext";
    Solname[4] = "mu";

    std::vector < unsigned > SolPdeIndex(n_unknowns);
    std::vector < unsigned > SolIndex(n_unknowns);
    std::vector < unsigned > SolFEType(n_unknowns);


    for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
        SolPdeIndex[ivar] = mlPdeSys->GetSolPdeIndex(Solname[ivar].c_str());
        SolIndex[ivar]    = ml_sol->GetIndex        (Solname[ivar].c_str());
        SolFEType[ivar]   = ml_sol->GetSolutionType(SolIndex[ivar]);
    }

    std::vector < unsigned int > Sol_n_el_dofs_Mat_vol(n_unknowns);

//***************************************************
    //----------- quantities (at dof objects) ------------------------------
    std::vector < int >       L2G_dofmap_Mat_AllVars;
    L2G_dofmap_Mat_AllVars.reserve( n_unknowns * max_size );
    std::vector < std::vector < int > >     L2G_dofmap_Mat(n_unknowns);
    for(int i = 0; i < n_unknowns; i++) {
        L2G_dofmap_Mat[i].reserve(max_size);
    }
    
    std::vector < std::vector < double > >  sol_eldofs_Mat(n_unknowns);
    for(int k = 0; k < n_unknowns; k++) {        sol_eldofs_Mat[k].reserve(max_size);    }
    
    std::vector < double > Res;
    Res.reserve( n_unknowns * max_size);
    std::vector < double > Jac;
    Jac.reserve( n_unknowns * max_size * n_unknowns * max_size);

//***************************************************
    //------------ quantities (at quadrature points) ---------------------
    std::vector <double>        sol_qp(n_unknowns);
    std::vector < std::vector <double> > sol_grad_qp(n_unknowns);

    std::fill(sol_qp.begin(), sol_qp.end(), 0.);
    for (unsigned  k = 0; k < n_unknowns; k++) {
        sol_grad_qp[k].resize(dim);
        std::fill(sol_grad_qp[k].begin(), sol_grad_qp[k].end(), 0.);
    }


//********************* DATA ************************
    const double u_des = ctrl::square_or_cube :: cost_functional_without_regularization::DesiredTargetVec()[0];
    const double alpha = ALPHA_CTRL_VOL;
    const double beta  = BETA_CTRL_VOL;
    const double penalty_strong_ctrl = 1.e30;
    const double penalty_strong_u =    1.e30;
    const double penalty_interface = 1.e10;         //penalty for u = q and for weak continuity of adjoint Neumann
//***************************************************

    RES->zero();
    if (assembleMatrix)  JAC->zero();


    // element loop: each process loops only on the elements that owns
    for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    short unsigned ielGeom = geom_element_iel.geom_type();

    int group_flag         = msh->GetElementGroup(iel);
//    std::cout << " ======= grp_flag === " << group_flag << " ================== " << std::endl;
//     int face_no         = msh->GetElementFaceNumber(iel);
//     std::cout << " ======= face# === " << face_no << " ================== " << std::endl;

//******************** GEOMETRY *********************
        unsigned nDofx = msh->GetElementDofNumber(iel, solType_coords);    // number of coordinate element dofs
        for (int i = 0; i < dim; i++)  coords_at_dofs[i].resize(nDofx);
        // local storage of coordinates
        for (unsigned i = 0; i < nDofx; i++) {
            unsigned xDof  = msh->GetSolutionDof(i, iel, solType_coords);  // global to global mapping between coordinates node and coordinate dof

            for (unsigned jdim = 0; jdim < dim; jdim++) {
                coords_at_dofs[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
            }
        }

        // elem average point
        std::vector < double > elem_center(dim);
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

//****** set target domain flag *********************
        int target_flag = 0;
        target_flag = ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(elem_center);
//***************************************************

        //all vars###################################################################
        for (unsigned  k = 0; k < n_unknowns; k++) {
            unsigned  ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
            Sol_n_el_dofs_Mat_vol[k] = ndofs_unk;
            sol_eldofs_Mat[k].resize(ndofs_unk);
            L2G_dofmap_Mat[k].resize(ndofs_unk);
            for (unsigned i = 0; i < ndofs_unk; i++) {
                unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);                        // global to global mapping between solution node and solution dof // via local to global solution node
                sol_eldofs_Mat[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);                            // global extraction and local storage for the solution
                L2G_dofmap_Mat[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
            }
        }
        //all vars###################################################################

        ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::update_active_set_flag_for_current_nonlinear_iteration
         (msh, sol, iel, coords_at_dofs, sol_eldofs_Mat, Sol_n_el_dofs_Mat_vol, 
          c_compl, 
          pos_mu_in_mat,
          pos_ctrl_in_mat, 
          solIndex_act_flag_sol,
          ctrl_lower,
          ctrl_upper, 
          sol_actflag);


//******************** ALL VARS *********************
    unsigned int nDof_max          = ElementJacRes<double>::compute_max_n_dofs(Sol_n_el_dofs_Mat_vol);
    unsigned int sum_Sol_n_el_dofs = ElementJacRes<double>::compute_sum_n_dofs(Sol_n_el_dofs_Mat_vol);
    
        
        Res.resize(sum_Sol_n_el_dofs);
        std::fill(Res.begin(), Res.end(), 0.);
        Jac.resize(sum_Sol_n_el_dofs * sum_Sol_n_el_dofs);
        std::fill(Jac.begin(), Jac.end(), 0.);

        L2G_dofmap_Mat_AllVars.resize(0);
        for (unsigned  k = 0; k < n_unknowns; k++)     L2G_dofmap_Mat_AllVars.insert(L2G_dofmap_Mat_AllVars.end(),L2G_dofmap_Mat[k].begin(),L2G_dofmap_Mat[k].end());
//***************************************************

        const int n_faces = msh->GetElementFaceNumber(iel);
        
        //setting up control boundary region - BEGIN ***************************
        std::vector <bool> interface_elem_flag(n_faces);      for(unsigned j = 0; j < n_faces; j++) interface_elem_flag[j] = false;
        std::vector <unsigned int> is_dof_on_Gamma_c(Sol_n_el_dofs_Mat_vol[pos_state]);   ///@todo maybe it should be geometry-based
        std::fill(is_dof_on_Gamma_c.begin(), is_dof_on_Gamma_c.end(), 0);

     for(unsigned jface = 0; jface < n_faces; jface++) {
            
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

       //************** later try to avoid repeating this - BEGIN
          lifting_external::compute_coordinates_bdry_one_face(coords_at_dofs_bdry, solType_coords, iel, jface, msh);
          
            unsigned nDofu_bdry    = msh->GetElementFaceDofNumber(iel,jface,SolFEType[pos_state]);
            unsigned nDofctrl_bdry = msh->GetElementFaceDofNumber(iel,jface,SolFEType[pos_ctrl]);
            if (nDofu_bdry != nDofctrl_bdry) {
                std::cout << "State and control need to have the same FE space" << std::endl;
                abort();
            }
         //************** later try to avoid repeating this - END
                    
           interface_elem_flag[jface] = lifting_external::find_control_boundary_nodes(is_dof_on_Gamma_c, geom_element_iel.get_elem_center_bdry_3d(), nDofu_bdry, iel, jface, msh);
    
        }
       //setting up control boundary region - END ***************************
       
        
        
        for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {


            
     lifting_external::compute_coordinates_bdry_one_face(coords_at_dofs_bdry, solType_coords, iel, jface, msh);

            unsigned nDofu_bdry    = msh->GetElementFaceDofNumber(iel,jface,SolFEType[pos_state]);
            unsigned nDofctrl_bdry = msh->GetElementFaceDofNumber(iel,jface,SolFEType[pos_ctrl]);
            if (nDofu_bdry != nDofctrl_bdry) {
                std::cout << "State and control need to have the same FE space" << std::endl;
                abort();
            }
            
    

    if (interface_elem_flag[jface] == true) {

              std::vector<double> normal_qp(dim, 0.);
              
            const unsigned felt_bdry = msh->GetElementFaceType(iel, jface);
              
          for(unsigned ig_bdry=0; ig_bdry < msh->_finiteElement[felt_bdry][SolFEType[pos_ctrl]]->GetGaussPointNumber(); ig_bdry++) {

                    // *** get gauss point weight, test function and test function partial derivatives ***
                    for(int fe=0; fe < NFE_FAMS; fe++) {
                        msh->_finiteElement[felt_bdry][fe]->JacobianSur(coords_at_dofs_bdry, ig_bdry, weight_qp_bdry, phi_fe_qp_bdry[fe], phi_x_fe_qp_bdry[fe], normal_qp);
                    }
                    //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
                    msh->_finiteElement[felt_bdry][solType_coords]->JacobianSur(coords_at_dofs_bdry, ig_bdry, weight_qp_bdry, phi_fe_qp_bdry[solType_coords], phi_x_fe_qp_bdry[solType_coords], normal_qp);

                    if (ielGeom != QUAD) {
                        std::cout << "fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem not implemented" << std::endl;
                        abort();
                    }
                    
                    msh->_finiteElement[ielGeom][SolFEType[pos_adj]]->fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem(coords_at_dofs, coords_at_dofs_bdry, jface, ig_bdry, phi_adj_vol_at_bdry, phi_adj_x_vol_at_bdry);
                    msh->_finiteElement[ielGeom][SolFEType[pos_adj_ext]]->fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem(coords_at_dofs, coords_at_dofs_bdry, jface, ig_bdry, phi_adj_ext_vol_at_bdry, phi_adj_ext_x_vol_at_bdry);

                    
//=============== grad dot n for residual =========================================
//     compute gauss quantities on the boundary through VOLUME interpolation
                    std::fill(sol_adj_x_vol_at_bdry_gss.begin(), sol_adj_x_vol_at_bdry_gss.end(), 0.);
                    for (int iv = 0; iv < Sol_n_el_dofs_Mat_vol[pos_adj]; iv++)  {

                        for (int d = 0; d < dim; d++) {
                            sol_adj_x_vol_at_bdry_gss[d] += sol_eldofs_Mat[pos_adj][iv] * phi_adj_x_vol_at_bdry[iv * dim + d];//notice that the convention of the orders x y z is different from vol to bdry
                        }
                    }

                    double grad_adj_dot_n_res = 0.;
                    for (unsigned d = 0; d < dim; d++) {
                        grad_adj_dot_n_res += sol_adj_x_vol_at_bdry_gss[d] * normal_qp[d];
                    }
//=============== grad dot n  for residual =========================================

//=============== grad dot n for residual =========================================
//     compute gauss quantities on the boundary through VOLUME interpolation
                    std::fill(sol_adj_ext_x_vol_at_bdry_gss.begin(), sol_adj_ext_x_vol_at_bdry_gss.end(), 0.);
                    for (int iv = 0; iv < Sol_n_el_dofs_Mat_vol[pos_adj_ext]; iv++)  {

                        for (int d = 0; d < dim; d++) {
                            sol_adj_ext_x_vol_at_bdry_gss[d] += sol_eldofs_Mat[pos_adj_ext][iv] * phi_adj_ext_x_vol_at_bdry[iv * dim + d];//notice that the convention of the orders x y z is different from vol to bdry
                        }
                    }

                    double grad_adj_ext_dot_n_res = 0.;
                    for (unsigned d = 0; d < dim; d++) {
                        grad_adj_ext_dot_n_res += sol_adj_ext_x_vol_at_bdry_gss[d] * normal_qp[d];
                    }
//=============== grad dot n  for residual =========================================


                    for (int i_bdry = 0; i_bdry < nDofu_bdry; i_bdry++)  {
                        unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

//============ Bdry Residuals - BEGIN ==================

                        if ( group_flag == GROUP_INTERNAL ) Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_state, i_vol)  ]  +=  -  weight_qp_bdry *  ( grad_adj_dot_n_res * phi_fe_qp_bdry[SolFEType[pos_state]][i_bdry] );

                        
                        if ( group_flag == GROUP_EXTERNAL ) Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_ctrl, i_vol)  ]  +=  -  weight_qp_bdry *  ( grad_adj_ext_dot_n_res * phi_fe_qp_bdry[SolFEType[pos_ctrl]][i_bdry] );
  ///@todo check this term in control equation
                        
                        const unsigned gamma_c_u_minus_q_pos = pos_adj/*pos_state*/;
                        const unsigned gamma_c_Neum_adj_continuity_pos = pos_adj_ext/*pos_ctrl*/;
                        
//                         if ( group_flag == GROUP_INTERNAL ) Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_adj, i_vol)   ]  +=  -  penalty_interface * ( sol_eldofs_Mat[pos_state][i_vol] - sol_eldofs_Mat[pos_ctrl][i_vol] ) ;    // u = q
                        if ( group_flag == GROUP_INTERNAL ) Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, gamma_c_u_minus_q_pos, i_vol)  ] +=  -  _u_minus_q_strong * penalty_interface * ( sol_eldofs_Mat[pos_state][i_vol]);    // u
                        if ( group_flag == GROUP_EXTERNAL ) Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, gamma_c_u_minus_q_pos, i_vol)  ] +=  -  _u_minus_q_strong * penalty_interface * ( - sol_eldofs_Mat[pos_ctrl][i_vol]);    // - q
                        
                        if ( group_flag == GROUP_INTERNAL ) Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, gamma_c_Neum_adj_continuity_pos, i_vol)   ]  +=  - _neumann_adjoint_explicit * penalty_interface *  weight_qp_bdry * phi_fe_qp_bdry[SolFEType[pos_adj_ext]][i_bdry] * ( _neumann_adjoint_continuity_sign ) * grad_adj_dot_n_res ;
                        if ( group_flag == GROUP_EXTERNAL ) Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, gamma_c_Neum_adj_continuity_pos, i_vol)   ]  +=  - _neumann_adjoint_explicit * penalty_interface *  weight_qp_bdry * phi_fe_qp_bdry[SolFEType[pos_adj_ext]][i_bdry] * (  grad_adj_ext_dot_n_res ) ;
//============ Bdry Residuals - END ==================

//============ Bdry Jacobians on Bdry - BEGIN ==================
                        for(unsigned j_bdry=0; j_bdry < nDofu_bdry; j_bdry ++) {
                            unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);

//============ u = q =============================
                            if (i_vol == j_vol)  {
                                if ( group_flag == GROUP_INTERNAL ) Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, gamma_c_u_minus_q_pos, pos_state, i_vol, j_vol) ]  += _u_minus_q_strong * penalty_interface *  ( 1.);
                                if ( group_flag == GROUP_EXTERNAL ) Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, gamma_c_u_minus_q_pos, pos_ctrl, i_vol, j_vol) ]   += _u_minus_q_strong * penalty_interface *  (-1.);
                            }
//============ u = q =============================

                        } //end j_bdry
//============ Bdry Jacobians on Bdry - END ==================

//============ Bdry Jacobians on Volume - BEGIN ==================
//===================loop over j in the VOLUME (while i is in the boundary)
                        for(unsigned j=0; j < nDof_max; j ++) {

//=============== grad dot n  =========================================
                            double grad_adj_dot_n_mat = 0.;
                            for(unsigned d=0; d<dim; d++) {
                                grad_adj_dot_n_mat += phi_adj_x_vol_at_bdry[j * dim + d] * normal_qp[d];  //notice that the convention of the orders x y z is different from vol to bdry
                            }
//=============== grad dot n  =========================================

//=============== grad dot n  =========================================
                            double grad_adj_ext_dot_n_mat = 0.;
                            for(unsigned d=0; d<dim; d++) {
                                grad_adj_ext_dot_n_mat += phi_adj_ext_x_vol_at_bdry[j * dim + d] * normal_qp[d];  //notice that the convention of the orders x y z is different from vol to bdry
                            }
//=============== grad dot n  =========================================

                            if ( group_flag == GROUP_INTERNAL ) {
                                Jac[  assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_state, pos_adj, i_vol, j) ]  += weight_qp_bdry * grad_adj_dot_n_mat * phi_fe_qp_bdry[SolFEType[pos_state]][i_bdry];
                            }
                            
                            if ( group_flag == GROUP_EXTERNAL ) {  ///@todo check this term in control equation
                                Jac[  assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_ctrl, pos_adj_ext, i_vol, j) ]  += weight_qp_bdry * grad_adj_ext_dot_n_mat * phi_fe_qp_bdry[SolFEType[pos_ctrl]][i_bdry];
                            }
                            

                                if ( group_flag == GROUP_INTERNAL ) Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, gamma_c_Neum_adj_continuity_pos, pos_adj, i_vol, j) ]  += _neumann_adjoint_explicit * penalty_interface * weight_qp_bdry * phi_fe_qp_bdry[SolFEType[pos_adj_ext]][i_bdry] * ( _neumann_adjoint_continuity_sign ) * grad_adj_dot_n_mat;
                                if ( group_flag == GROUP_EXTERNAL ) Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, gamma_c_Neum_adj_continuity_pos, pos_adj_ext, i_vol, j) ]   += _neumann_adjoint_explicit * penalty_interface * weight_qp_bdry * phi_fe_qp_bdry[SolFEType[pos_adj_ext]][i_bdry] * ( 1.) * grad_adj_ext_dot_n_mat;


                        } //end j
//============ Bdry Jacobians on Volume - END ==================

                    } //end i_bdry

                }  //end ig_bdry loop


            }  //bdry interface elements

        }  //end boundary face loop


        
        
//             for (unsigned j = 0; j < is_dof_on_Gamma_c.size(); j++) {
//                 std::cout <<  " ** " << is_dof_on_Gamma_c[j] << " ";
//             }
//             
//            std::cout << std::endl; 
            

        for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {

            // *** get gauss point weight, test function and test function partial derivatives ***
            for(int fe=0; fe < NFE_FAMS; fe++) {
                msh->_finiteElement[ielGeom][fe]->Jacobian(coords_at_dofs,ig,weight_qp,phi_fe_qp[fe],phi_x_fe_qp[fe],phi_xx_fe_qp[fe]);
            }
            //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
            msh->_finiteElement[ielGeom][solType_coords]->Jacobian(coords_at_dofs,ig,weight_qp,phi_fe_qp[solType_coords],phi_x_fe_qp[solType_coords],phi_xx_fe_qp[solType_coords]);

//========= fill gauss value xyz ==================
            std::fill(coord_at_qp.begin(), coord_at_qp.end(), 0.);
            for (unsigned  d = 0; d < dim; d++) {
                for (unsigned i = 0; i < coords_at_dofs[d].size(); i++) {
                    coord_at_qp[d] += coords_at_dofs[d][i] * phi_fe_qp[solType_coords][i];
                }
            }
            //========= fill gauss value xyz ==================

//========= fill gauss value quantities ==================
            std::fill(sol_qp.begin(), sol_qp.end(), 0.);
            for (unsigned  k = 0; k < n_unknowns; k++) {
                std::fill(sol_grad_qp[k].begin(), sol_grad_qp[k].end(), 0.);
            }

            for (unsigned  k = 0; k < n_unknowns; k++) {
                for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[k]; i++) {
                    sol_qp[k]    += sol_eldofs_Mat[k][i] *   phi_fe_qp[SolFEType[k]][i];
                    for (unsigned d = 0; d < dim; d++)   sol_grad_qp[k][d] += sol_eldofs_Mat[k][i] * phi_x_fe_qp[SolFEType[k]][i * dim + d];
                }
            }
//========= fill gauss value quantities ==================





//==========FILLING WITH THE EQUATIONS ===========
            // *** phi_i loop ***
            for (unsigned i = 0; i < nDof_max; i++) {

                double laplace_rhs_du_adj_i = 0.;
                double laplace_rhs_dadj_u_i = 0.;
                double laplace_rhs_dctrl_adj_ext_i = 0.;
                double laplace_rhs_dadj_ext_ctrl_i = 0.;
                double laplace_rhs_dctrl_ctrl_i = 0.;

                for (unsigned kdim = 0; kdim < dim; kdim++) {
                    if ( i < Sol_n_el_dofs_Mat_vol[pos_state] )       laplace_rhs_du_adj_i          +=  (phi_x_fe_qp[SolFEType[pos_state]]   [i * dim + kdim] * sol_grad_qp[pos_adj][kdim]);
                    if ( i < Sol_n_el_dofs_Mat_vol[pos_adj] )         laplace_rhs_dadj_u_i          +=  (phi_x_fe_qp[SolFEType[pos_adj]]     [i * dim + kdim] * sol_grad_qp[pos_state][kdim]);
                    if ( i < Sol_n_el_dofs_Mat_vol[pos_ctrl] )        laplace_rhs_dctrl_adj_ext_i   +=  (phi_x_fe_qp[SolFEType[pos_ctrl]]    [i * dim + kdim] * sol_grad_qp[pos_adj_ext][kdim]);
                    if ( i < Sol_n_el_dofs_Mat_vol[pos_adj_ext] )     laplace_rhs_dadj_ext_ctrl_i   +=  (phi_x_fe_qp[SolFEType[pos_adj_ext]] [i * dim + kdim] * sol_grad_qp[pos_ctrl][kdim]);
                    if ( i < Sol_n_el_dofs_Mat_vol[pos_ctrl] )        laplace_rhs_dctrl_ctrl_i      +=  (phi_x_fe_qp[SolFEType[pos_ctrl]]    [i * dim + kdim] * sol_grad_qp[pos_ctrl][kdim]);
                }

//======================Volume Residuals - BEGIN =======================
//--- equations ---                
                if ( group_flag == GROUP_INTERNAL )     {
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol,pos_state,i) ]  += - weight_qp * (target_flag * phi_fe_qp[SolFEType[pos_state]][i] * ( sol_qp[pos_state] - u_des) - laplace_rhs_du_adj_i);
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol,pos_adj,i)]     += - weight_qp *  ( - laplace_rhs_dadj_u_i    - 0.) ;
                }

                else if ( group_flag == GROUP_EXTERNAL )  {
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol,pos_ctrl,i) ]   += - weight_qp * ( OP_L2 * alpha * phi_fe_qp[SolFEType[pos_ctrl]][i] * sol_qp[pos_ctrl]
                            + OP_H1 * beta * laplace_rhs_dctrl_ctrl_i
                            - laplace_rhs_dctrl_adj_ext_i /*+ 7000. * phi_fe_qp[SolFEType[pos_ctrl]][i]*/);
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol,pos_adj_ext,i)] += - weight_qp *  ( - laplace_rhs_dadj_ext_ctrl_i - 0.);
                }

//--- extensions to zero ---
            const double exclude_Dirichlet_Gamma_c_for_Adjoint =  (1 - is_dof_on_Gamma_c[i]);
            const double exclude_Dirichlet_Gamma_c_for_Adjoint_Ext =  (1 - is_dof_on_Gamma_c[i]);
            
                if ( group_flag == GROUP_INTERNAL )     {
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_adj_ext, i)] += - exclude_Dirichlet_Gamma_c_for_Adjoint_Ext * penalty_strong_ctrl * ( sol_eldofs_Mat[pos_adj_ext][i] - 0.);
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_ctrl, i)]    += - (1 - is_dof_on_Gamma_c[i]) * penalty_strong_ctrl * (sol_eldofs_Mat[pos_ctrl][i] - 0.);
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_mu, i)]      += - (1 - is_dof_on_Gamma_c[i]) * penalty_strong_ctrl * ( sol_eldofs_Mat[pos_mu][i] - 0.);
                }

                else if ( group_flag == GROUP_EXTERNAL )  {
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_adj, i)]     += - exclude_Dirichlet_Gamma_c_for_Adjoint * penalty_strong_u * (  sol_eldofs_Mat[pos_adj][i] - 0.);
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs_Mat_vol, pos_state, i) ]  += - (1 - is_dof_on_Gamma_c[i]) * penalty_strong_u * (sol_eldofs_Mat[pos_state][i] - 0.);
                }
//======================Volume Residuals - END =======================

                if (assembleMatrix) {

                    // *** phi_j loop ***
                    for (unsigned j = 0; j < nDof_max; j++) {

                        double laplace_mat_du_adj = 0.;
                        double laplace_mat_dadj_u = 0.;
                        double laplace_mat_dctrl_adj_ext = 0.;
                        double laplace_mat_dadj_ext_ctrl = 0.;
                        double laplace_mat_dctrl_ctrl = 0.;

                        for (unsigned kdim = 0; kdim < dim; kdim++) {
                            if ( i < Sol_n_el_dofs_Mat_vol[pos_state]   && j < Sol_n_el_dofs_Mat_vol[pos_adj] )      laplace_mat_du_adj        += (phi_x_fe_qp[SolFEType[pos_state]]  [i * dim + kdim] * phi_x_fe_qp[SolFEType[pos_adj]]    [j * dim + kdim]);
                            if ( i < Sol_n_el_dofs_Mat_vol[pos_adj]     && j < Sol_n_el_dofs_Mat_vol[pos_state] )    laplace_mat_dadj_u        += (phi_x_fe_qp[SolFEType[pos_adj]]    [i * dim + kdim] * phi_x_fe_qp[SolFEType[pos_state]]  [j * dim + kdim]);  //equal to the previous
                            if ( i < Sol_n_el_dofs_Mat_vol[pos_ctrl]    && j < Sol_n_el_dofs_Mat_vol[pos_adj_ext] )  laplace_mat_dctrl_adj_ext += (phi_x_fe_qp[SolFEType[pos_ctrl]]   [i * dim + kdim] * phi_x_fe_qp[SolFEType[pos_adj_ext]][j * dim + kdim]);
                            if ( i < Sol_n_el_dofs_Mat_vol[pos_adj_ext] && j < Sol_n_el_dofs_Mat_vol[pos_ctrl] )     laplace_mat_dadj_ext_ctrl += (phi_x_fe_qp[SolFEType[pos_adj_ext]][i * dim + kdim] * phi_x_fe_qp[SolFEType[pos_ctrl]]   [j * dim + kdim]);  //equal to the previous
                            if ( i < Sol_n_el_dofs_Mat_vol[pos_ctrl]    && j < Sol_n_el_dofs_Mat_vol[pos_ctrl] )     laplace_mat_dctrl_ctrl    += (phi_x_fe_qp[SolFEType[pos_ctrl]]   [i * dim + kdim] * phi_x_fe_qp[SolFEType[pos_ctrl]]   [j * dim + kdim]);
                        }


//--- equations ---                
                        if ( group_flag == GROUP_INTERNAL ) {

                            Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_state, pos_state, i, j) ] += weight_qp * target_flag * phi_fe_qp[SolFEType[pos_state]][j] * phi_fe_qp[SolFEType[pos_state]][i];
                            Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_state, pos_adj, i, j) ]   += weight_qp * (-1.) * laplace_mat_du_adj;

                            Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_adj, pos_state, i, j) ]   += weight_qp * (-1.) * laplace_mat_dadj_u;

                        }

                        else if ( group_flag == GROUP_EXTERNAL ) {

                            
                            Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_ctrl, pos_ctrl, i, j) ]     += weight_qp * ( 
                                 OP_L2 * alpha * phi_fe_qp[SolFEType[pos_ctrl]][i] * phi_fe_qp[SolFEType[pos_ctrl]][j]
                                + OP_H1 * beta * laplace_mat_dctrl_ctrl
                                    );

                            Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_ctrl, pos_adj_ext, i, j) ]  += weight_qp * (-1.) * laplace_mat_dctrl_adj_ext;

                            Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_adj_ext, pos_ctrl, i, j) ]  += weight_qp * (-1.) * laplace_mat_dadj_ext_ctrl;

                        }

                        

//--- extensions to zero ---                
                       if (  i == j ) {
                                
                        if ( group_flag == GROUP_INTERNAL ) {
                             Jac[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_adj_ext, pos_adj_ext, i, j)]  += exclude_Dirichlet_Gamma_c_for_Adjoint_Ext * penalty_strong_ctrl;
                             Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_ctrl, pos_ctrl, i, j) ]      += (1 - is_dof_on_Gamma_c[i]) * penalty_strong_ctrl;
                             Jac[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_mu, pos_mu, i, j) ]           += (1 - is_dof_on_Gamma_c[i]) * penalty_strong_ctrl;
                            }
                            
                        else if ( group_flag == GROUP_EXTERNAL ) {
                             Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_adj, pos_adj, i, j) ]      += exclude_Dirichlet_Gamma_c_for_Adjoint *  penalty_strong_u;
                             Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs_Mat_vol, sum_Sol_n_el_dofs, pos_state, pos_state, i, j) ]  += (1 - is_dof_on_Gamma_c[i]) *  penalty_strong_u;
                            }

                        }
                        
                        
                        
                        
                        
                        
                        
                    } // end phi_j loop

                } // endif assemble_matrix

            } // end phi_i loop

        } // end gauss point loop


//========== sum-based part ================================

        RES->add_vector_blocked(Res, L2G_dofmap_Mat_AllVars);
        if (assembleMatrix) JAC->add_matrix_blocked(Jac, L2G_dofmap_Mat_AllVars, L2G_dofmap_Mat_AllVars);

        
        
     if (print_algebra_local) {
    assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs_Mat_vol, 10, 5);
    assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs_Mat_vol, 10, 5);
     }
        
  } //end element loop for each process
       
       
       
     //this part is to add a bunch of 1's in the residual, so it is not based on the element loop
    unsigned int ctrl_index = mlPdeSys->GetSolPdeIndex("control");

    unsigned int global_ctrl_size = pdeSys->KKoffset[ctrl_index+1][iproc] - pdeSys->KKoffset[ctrl_index][iproc];

    std::vector<double>  one_times_mu(global_ctrl_size, 0.);
    std::vector<int>    positions(global_ctrl_size);

    for (unsigned i = 0; i < positions.size(); i++) {
        positions[i] = pdeSys->KKoffset[ctrl_index][iproc] + i;
        one_times_mu[i] = ineq_flag * 1. * (*sol->_Sol[SolIndex[pos_mu]])(i/*position_mu_i*/) ;
    }
    RES->add_vector_blocked(one_times_mu, positions);

    
    RES->close();
    if (assembleMatrix)   JAC->close();
      // ***************** ADD PART - END  *******************
      
        
        
        
//   ***************** INSERT PART - BEGIN (must go AFTER the sum, clearly) *******************
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

      

// -------
   geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
      
   geom_element_iel.set_elem_center_3d(iel, solType_coords);
// -------
   
// -------
    el_dofs_unknowns_vol(sol, msh, pdeSys, iel,
                        SolFEType,
                        SolIndex,
                        SolPdeIndex,
                        Sol_n_el_dofs_Mat_vol,
                        sol_eldofs_Mat,
                        L2G_dofmap_Mat);
      
    //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ctrl:: square_or_cube:: lifting_internal < ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
 
    
    if (control_el_flag == 1) {

        
  ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::update_active_set_flag_for_current_nonlinear_iteration
  (msh,
   sol,
   iel,
   geom_element_iel.get_coords_at_dofs/*_3d*/(),
   sol_eldofs_Mat,
   Sol_n_el_dofs_Mat_vol,
   c_compl,
   pos_mu_in_mat,
   pos_ctrl_in_mat,
   solIndex_act_flag_sol,
   ctrl_lower,
   ctrl_upper,
   sol_actflag);
  
      


    ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::node_insertion(iel,
                   msh,
                   L2G_dofmap_Mat,
                   pos_mu_in_mat,
                   pos_ctrl_in_mat,
                   sol_eldofs_Mat,
                   Sol_n_el_dofs_Mat_vol,
                   sol_actflag,
                   ctrl_lower,
                   ctrl_upper,
                   ineq_flag,
                   c_compl,
                   JAC,
                   RES,
                   assembleMatrix
                   );
   

       }
      
      
// // // //============= delta_mu row ===============================
// // //         std::vector<double> Res_mu (Sol_n_el_dofs[pos_mu]);
// // //         std::fill(Res_mu.begin(),Res_mu.end(), 0.);
// // //         for (unsigned i = 0; i < sol_actflag.size(); i++) {
// // //             if (sol_actflag[i] == 0) {  //inactive
// // //                 Res_mu [i] = - ineq_flag * ( 1. * sol_eldofs[pos_mu][i] - 0. );
// // //             }
// // //             else if (sol_actflag[i] == 1) {  //active_a
// // //                 Res_mu [i] = - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i] - c_compl * ctrl_lower[i]);
// // //             }
// // //             else if (sol_actflag[i] == 2) {  //active_b
// // //                 Res_mu [i]  = - ineq_flag * ( c_compl * sol_eldofs[pos_ctrl][i] - c_compl * ctrl_upper[i]);
// // //             }
// // //         }
// // // 
// // //      
// // //         RES->insert(Res_mu, L2G_dofmap_Mat[pos_mu]);
// // // 
// // //         //============= delta_mu-delta_ctrl row ===============================
// // //         for (unsigned i = 0; i < sol_actflag.size(); i++) if (sol_actflag[i] != 0 ) sol_actflag[i] = ineq_flag * c_compl;
// // // 
// // //         JAC->matrix_set_off_diagonal_values_blocked(L2G_dofmap_Mat[pos_mu], L2G_dofmap_Mat[pos_ctrl], sol_actflag);
// // // 
// // //         //============= delta_mu-delta_mu row ===============================
// // //         for (unsigned i = 0; i < sol_actflag.size(); i++) sol_actflag[i] =  ineq_flag * (1 - sol_actflag[i]/c_compl)  + (1-ineq_flag) * 1.;
// // // 
// // //         JAC->matrix_set_off_diagonal_values_blocked(L2G_dofmap_Mat[pos_mu], L2G_dofmap_Mat[pos_mu], sol_actflag );
// // // 
     
        //============= delta_ctrl-delta_mu row ===============================
        JAC->matrix_set_off_diagonal_values_blocked( L2G_dofmap_Mat[pos_ctrl], L2G_dofmap_Mat[pos_mu], ineq_flag * 1.);

     
    } //end element loop for each process    
 //   ***************** INSERT PART - END (must go AFTER the sum, clearly) *******************
   
  RES->close();
  if (assembleMatrix) JAC->close();

 
  
  print_global_residual_jacobian(print_algebra_global,
                                 ml_prob,
                                 mlPdeSys,
                                 pdeSys,
                                 RES,
                                 JAC,
                                 iproc,
                                 assembleMatrix);
  


  
    return;
}


};


} //end namespace


} //end namespace



#endif
 
