#ifndef OPT_SYSTEMS_NS_DIRICHLET_HPP
#define OPT_SYSTEMS_NS_DIRICHLET_HPP

#include "01_opt_system.hpp"

#include  "opt_systems_boundary_control_eqn_sobolev_integer.hpp"

#include "SparseMatrix.hpp"




namespace femus  {


namespace navier_stokes  {


class pure_boundary : public femus::pure_boundary {

public: 
   

 //Unknown definition  - BEGIN ==================
   enum pos_vector_quantities {pos_index_state = 0, pos_index_adj, pos_index_ctrl, pos_index_mu};
 
static const std::vector< unsigned >  provide_state_adj_ctrl_mu_offsets(const unsigned int dimension) {
     
     std::vector<unsigned>  opt_control_offsets(4);  
     
     opt_control_offsets[pos_index_state] = 0; 
     opt_control_offsets[pos_index_adj]   = dimension + 1; 
     
     opt_control_offsets[pos_index_mu]    = 2 * (dimension + 1); 
     
#if INEQ_FLAG != 0
  opt_control_offsets[pos_index_ctrl]  = opt_control_offsets[pos_index_mu] + dimension;
#else
  opt_control_offsets[pos_index_ctrl]  =  2 * (dimension + 1);
#endif
     
     
     return opt_control_offsets;
     
 } 
 
 
static const std::vector< Unknown >  provide_list_of_unknowns(const unsigned int dimension) {


 std::vector< Unknown >  unknowns( 3*(dimension+1) + dimension);
 
 
  std::vector< unsigned >  vector_offsets = provide_state_adj_ctrl_mu_offsets(dimension);


const int state_pos_begin   =  vector_offsets[pos_index_state];
  const int adj_pos_begin   =  vector_offsets[pos_index_adj];
 const int ctrl_pos_begin   =  vector_offsets[pos_index_ctrl];
   const int mu_pos_begin   =  vector_offsets[pos_index_mu];
   

  //--- names - BEGIN
                                        unknowns[state_pos_begin + 0]._name      = "u_0";
                                        unknowns[state_pos_begin + 1]._name      = "u_1";
  if (dimension == 3)                   unknowns[state_pos_begin + 2]._name      = "u_2";
                                unknowns[state_pos_begin + dimension]._name      = "u_p";
  
                        unknowns[adj_pos_begin + 0]._name      = "adj_0";
                        unknowns[adj_pos_begin + 1]._name      = "adj_1";
  if (dimension == 3)   unknowns[adj_pos_begin + 2]._name      = "adj_2";
                unknowns[adj_pos_begin + dimension]._name      = "adj_p";
  
                       unknowns[ctrl_pos_begin + 0]._name      = "ctrl_0";
                       unknowns[ctrl_pos_begin + 1]._name      = "ctrl_1";
  if (dimension == 3)  unknowns[ctrl_pos_begin + 2]._name      = "ctrl_2";
  
               unknowns[ctrl_pos_begin + dimension]._name      = "theta";

#if INEQ_FLAG != 0
                       unknowns[mu_pos_begin + 0]._name      = "mu_0";
                       unknowns[mu_pos_begin + 1]._name      = "mu_1";
  if (dimension == 3)  unknowns[mu_pos_begin + 2]._name      = "mu_2";
#endif  
 
   //--- names - END
 

  //--- fe family - BEGIN
                                        unknowns[state_pos_begin + 0]._fe_family      = LAGRANGE;
                                        unknowns[state_pos_begin + 1]._fe_family      = LAGRANGE;
  if (dimension == 3)                   unknowns[state_pos_begin + 2]._fe_family      = LAGRANGE;
                                unknowns[state_pos_begin + dimension]._fe_family      = LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/;
  
                        unknowns[adj_pos_begin + 0]._fe_family      =  LAGRANGE;
                        unknowns[adj_pos_begin + 1]._fe_family      =  LAGRANGE;
  if (dimension == 3)   unknowns[adj_pos_begin + 2]._fe_family      =  LAGRANGE;
                unknowns[adj_pos_begin + dimension]._fe_family      =  LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/;
  
                       unknowns[ctrl_pos_begin + 0]._fe_family      =   LAGRANGE;
                       unknowns[ctrl_pos_begin + 1]._fe_family      =   LAGRANGE;
  if (dimension == 3)  unknowns[ctrl_pos_begin + 2]._fe_family      =   LAGRANGE;
                                                                        
               unknowns[ctrl_pos_begin + dimension]._fe_family      =   DISCONTINUOUS_POLYNOMIAL;     //theta

#if INEQ_FLAG != 0
                       unknowns[mu_pos_begin + 0]._fe_family       =   LAGRANGE;
                       unknowns[mu_pos_begin + 1]._fe_family       =   LAGRANGE;
  if (dimension == 3)  unknowns[mu_pos_begin + 2]._fe_family       =   LAGRANGE;
#endif  
 
  //--- fe family - END


   //--- fe order - BEGIN
                                        unknowns[state_pos_begin + 0]._fe_order      = SECOND;
                                        unknowns[state_pos_begin + 1]._fe_order      = SECOND;
  if (dimension == 3)                   unknowns[state_pos_begin + 2]._fe_order      = SECOND;
                                unknowns[state_pos_begin + dimension]._fe_order      = FIRST;
  
                        unknowns[adj_pos_begin + 0]._fe_order          = SECOND;
                        unknowns[adj_pos_begin + 1]._fe_order          = SECOND;
  if (dimension == 3)   unknowns[adj_pos_begin + 2]._fe_order          = SECOND;
                unknowns[adj_pos_begin + dimension]._fe_order          = FIRST;
  
                       unknowns[ctrl_pos_begin + 0]._fe_order        = SECOND;
                       unknowns[ctrl_pos_begin + 1]._fe_order        = SECOND;
  if (dimension == 3)  unknowns[ctrl_pos_begin + 2]._fe_order        = SECOND;
                                                                        
               unknowns[ctrl_pos_begin + dimension]._fe_order      =   ZERO;     //theta

#if INEQ_FLAG != 0
                       unknowns[mu_pos_begin + 0]._fe_order          = SECOND;
                       unknowns[mu_pos_begin + 1]._fe_order          = SECOND;
  if (dimension == 3)  unknowns[mu_pos_begin + 2]._fe_order          = SECOND;
#endif  
 
  
  
  
  
   //--- fe order - END
     

  
  
     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._time_order = 0;
              unknowns[u]._is_pde_unknown = true;
              unknowns[u]._is_sparse = true;
              
     }
     
 
 for (unsigned int u = 0; u < dimension; u++) {
   unknowns[ctrl_pos_begin + u]._is_sparse = IS_CTRL_FRACTIONAL_SOBOLEV ? false: true;
     }
 
 
 unknowns[ctrl_pos_begin + dimension]._is_sparse = false;    //theta
 
 
   return unknowns;
     
}

 //Unknown definition  - END ==================

           

static void assemble_ns_dirichlet_control_pure_boundary(MultiLevelProblem& ml_prob) {

  // Main objects - BEGIN *******************************************
 //System
 NonLinearImplicitSystemWithPrimalDualActiveSetMethod * mlPdeSys  = & ml_prob.get_system< NonLinearImplicitSystemWithPrimalDualActiveSetMethod >("NSOpt");
  
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  //System, Equation
  const char* system_name            = mlPdeSys->name().c_str();
  LinearEquationSolver*  pdeSys	 = mlPdeSys->_LinSolver[level];   
  bool assembleMatrix = mlPdeSys->GetAssembleMatrix(); 
   
  SparseMatrix*	 JAC	= pdeSys->_KK;
  NumericVector* RES 	= pdeSys->_RES;
    

  constexpr bool print_algebra_global = true;
  constexpr bool print_algebra_local = false;

  //Mesh
  Mesh*		 msh    = ml_prob._ml_msh->GetLevel(level);
  elem*		 el	= msh->GetMeshElements();
  const unsigned dim 	= msh->GetDimension();
  unsigned dim2     = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned nel		= msh->GetNumberOfElements();
  unsigned igrid	= msh->GetLevel();
  unsigned iproc 	= msh->processor_id();
  unsigned   nprocs = msh->n_processors();
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));


  //Solution
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;
  Solution*	 sol  	         = ml_prob._ml_sol->GetSolutionLevel(level);
  
  // Main objects - END *******************************************


  // ======= Geometry at Dofs - BEGIN  =======
   unsigned coordXType = 2; /*CONTINUOUS_BIQUADRATIC*/// get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
   unsigned solType_coords = coordXType;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = /*dim*/  3  /*2*/    ;
 
  std::vector < std::vector < double> > coordX(dim);
  for(int i = 0; i < dim; i++) {   coordX[i].reserve(max_size);   }
  
  std::vector < std::vector < double> > coordX_bd(/*space_*/dim);
  for(int i = 0; i < dim; i++) { coordX_bd[i].reserve(max_size);  }
 
  const unsigned dim_bdry = dim - 1;
  // ======= Geometry at Dofs - END  =======


 // ======= Solutions, Unknowns - BEGIN =======

  std::vector< unsigned >  vector_offsets = provide_state_adj_ctrl_mu_offsets(dim);


const int state_pos_begin   =  vector_offsets[pos_index_state];
  const int adj_pos_begin   =  vector_offsets[pos_index_adj];
 const int ctrl_pos_begin   =  vector_offsets[pos_index_ctrl];
   const int mu_pos_begin   =  vector_offsets[pos_index_mu];
   


  
    const int press_type_pos   = dim;

  const int theta_index      = ctrl_pos_begin + dim;
  
  
    std::vector<std::string> ctrl_name;
    ctrl_name.resize(dim);
    ctrl_name[0] = "ctrl_0";
    ctrl_name[1] = "ctrl_1";
    if (dim == 3)  ctrl_name[2] = "ctrl_2";
 
  const int pos_mat_ctrl = ctrl_pos_begin;
  const int pos_sol_ctrl = ctrl_pos_begin;

    const unsigned int n_components_ctrl = dim;

  std::vector< Unknown > unknowns = provide_list_of_unknowns( dim );
  
  const int n_unknowns = unknowns.size();
    

  std::vector < std::string > Solname_Mat(n_unknowns);  
  std::vector < unsigned > SolPdeIndex(n_unknowns);
  std::vector < unsigned > SolIndex_Mat(n_unknowns);  
  std::vector < unsigned > SolFEType_Mat(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    Solname_Mat[ivar] = unknowns[ivar]._name;
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(Solname_Mat[ivar].c_str());
    SolIndex_Mat[ivar]	= ml_sol->GetIndex        (Solname_Mat[ivar].c_str());
    SolFEType_Mat[ivar]	= ml_sol->GetSolutionType(SolIndex_Mat[ivar]);
  }

  std::vector < unsigned > Sol_n_el_dofs_Mat_vol(n_unknowns);
  

  const double cost_functional_coeff = COST_FUNCTIONAL_COEFF;
  const double alpha = ALPHA_CTRL_BDRY;
  const double beta  = BETA_CTRL_BDRY;
 // ======= Solutions, Unknowns - END =======
  
     
  //=== Sol (quantities, not only unknowns) - BEGIN ========================================================
  
  const unsigned int n_quantities = ml_sol->GetSolutionSize();
  
 //***************************************************
    std::vector < std::string > Solname_quantities(n_quantities);
    
        for(unsigned ivar=0; ivar < Solname_quantities.size(); ivar++) {
            Solname_quantities[ivar] = ml_sol->GetSolName_from_index(ivar);
        }
 //***************************************************
 
    std::vector < unsigned > SolIndex_quantities(n_quantities);      //should have Sol order
    std::vector < unsigned > SolFEType_quantities(n_quantities);     //should have Sol order
    std::vector < unsigned > Sol_n_el_dofs_quantities(n_quantities); //should have Sol order
    
    for(unsigned ivar=0; ivar < n_quantities; ivar++) {
        SolIndex_quantities[ivar]    = ml_sol->GetIndex        (Solname_quantities[ivar].c_str());
        SolFEType_quantities[ivar]   = ml_sol->GetSolutionType(SolIndex_quantities[ivar]);
    }
    
    
#if INEQ_FLAG != 0
  //MU
  //************** variables for ineq constraints: act flag ****************************   
  std::vector<unsigned int> solIndex_act_flag_sol(n_components_ctrl); 
  
  ctrl::mixed_state_or_ctrl_inequality< femus:: ctrl:: square_or_cube::mixed_state_or_ctrl_inequality >::store_act_flag_in_old(mlPdeSys, ml_sol, sol, solIndex_act_flag_sol);
  //************** variables for ineq constraints: act flag ****************************   
    

  const int ineq_flag = INEQ_FLAG;
  const double c_compl = C_COMPL;
  std::vector < std::vector < double/*int*/ > > sol_actflag(n_components_ctrl);   
  std::vector < std::vector < double > > ctrl_lower(n_components_ctrl);           
  std::vector < std::vector < double > > ctrl_upper(n_components_ctrl);           
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
      sol_actflag[c].reserve(max_size);
      ctrl_lower[c].reserve(max_size);
      ctrl_upper[c].reserve(max_size);
      }

    
//MU
  std::vector<unsigned int> ctrl_index_in_mat(n_components_ctrl); 
  std::vector<unsigned int>   mu_index_in_mat(n_components_ctrl);    
      for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
           ctrl_index_in_mat[kdim] =  SolPdeIndex[ctrl_pos_begin + kdim];
             mu_index_in_mat[kdim] =  SolPdeIndex[mu_pos_begin + kdim];
      }
#endif
      
    
    
  //=== Sol (quantities, not only unknowns) - END ========================================================
  

  
  // ======= Solutions, Unknowns at dofs - BEGIN =======
  std::vector < std::vector < double > > Sol_eldofs_Mat(n_unknowns);
  std::vector < std::vector < double > > gradSol_eldofs_Mat(n_unknowns);
  
  for(int k = 0; k < n_unknowns; k++) {
    Sol_eldofs_Mat[k].reserve(max_size);
    gradSol_eldofs_Mat[k].reserve(max_size*dim);    
  }
  
  
// ****** Solutions, Unknowns at dofs, theta value from proc0 - BEGIN 
  double solTheta = femus::get_theta_value(msh->n_processors(), sol, SolIndex_Mat[theta_index]);
// ****** Solutions, Unknowns at dofs, theta value from proc0  - END

  // ======= Solutions, Unknowns at dofs - END =======


  // ======= Solutions, Unknowns at quadrature points - BEGIN =======
    std::vector < double > SolVAR_qp(n_unknowns);   //sol_V,P_gss_of_st,adj,ctrl_ie@quadraturepoints
    std::vector < std::vector < double > > gradSolVAR_qp(n_unknowns);
    for(int k = 0; k < n_unknowns; k++) {  gradSolVAR_qp[k].resize(dim_offset_grad /*space_dim*/);  }


  std::vector < std::vector < double > >  sol_adj_x_vol_at_bdry_gss(dim);
  for (int ldim =0; ldim < dim; ldim++) sol_adj_x_vol_at_bdry_gss[ldim].reserve(max_size);
  
  std::vector < double > grad_adj_dot_n_res_qp;
  std::vector < double > grad_adj_dot_n_jac_qp;
  grad_adj_dot_n_res_qp.reserve(max_size);
  grad_adj_dot_n_jac_qp.reserve(max_size);
  // ======= Solutions, Unknowns at quadrature points - END =======
      
 
  // ====== Geometry at Quadrature points - BEGIN ==============================================================================
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
    double detJac_iqp;
    double AbsDetJxWeight_iqp = 0.;

     std::vector < std::vector < double > >  JacI_iqp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_iqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_iqp_bdry.size(); d++) {   Jac_iqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp_bdry.size(); d++) { JacI_iqp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
    double AbsDetJxWeight_iqp_bdry = 0.;
    
    std::vector<double> normal_iqp(dim_offset_grad /*space_dim*/, 0.);
  // ====== Geometry at Quadrature points - END ==============================================================================
  
     
// ======= FE at Quadrature, all - BEGIN =======
      //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all;
  ml_prob.get_all_abstract_fe_multiple(elem_all);
  
  //==========================================================================================
  std::vector < std::vector < double > > phi_gss_fe(NFE_FAMS);
  std::vector < std::vector < double > > phi_x_gss_fe(NFE_FAMS);
  
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(max_size);
      phi_x_gss_fe[fe].reserve(max_size * dim_offset_grad);
   }
   
   
  //boundary adjoint & ctrl shape functions  
  std::vector < std::vector < double > > phi_bd_gss_fe(NFE_FAMS);
  std::vector < std::vector < double > > phi_x_bd_gss_fe(NFE_FAMS);

    for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_bd_gss_fe[fe].reserve(max_size);
      phi_x_bd_gss_fe[fe].reserve(max_size * dim_offset_grad);
  //bdry vol adj  evaluated at bdry points
    }
    
    
  //bdry vol adj  evaluated at bdry points
   std::vector < std::vector < double > > phi_vol_at_bdry_fe(NFE_FAMS);
   std::vector < std::vector < double > > phi_x_vol_at_bdry_fe(NFE_FAMS);

    for(int fe=0; fe < NFE_FAMS; fe++) {  
         phi_vol_at_bdry_fe[fe].reserve(max_size);
       phi_x_vol_at_bdry_fe[fe].reserve(max_size * dim_offset_grad);
    }
  //==========================================================================================

 //********************* bdry cont *******************
 //*************************************************** 
  std::vector <double> phi_ctrl_bdry;  
  std::vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(max_size);
  phi_ctrl_x_bdry.reserve(max_size * space_dim);
 //*************************************************** 
// ======= FE at Quadrature, all - END =======
    
  
  // ======= Equation, local - BEGIN =======
  std::vector < std::vector < int > > L2G_dofmap_Mat(n_unknowns); 
  std::vector < std::vector < double > > Res(n_unknowns);
  std::vector < std::vector < std::vector < double > > > Jac(n_unknowns);
  
  for(int i = 0; i < n_unknowns; i++) {     
    L2G_dofmap_Mat[i].reserve(max_size);
       Res[i].reserve(max_size);
  }
   
  if (assembleMatrix) {
    
    for(int i = 0; i < n_unknowns; i++) {
      Jac[i].resize(n_unknowns);    
      for(int j = 0; j < n_unknowns; j++) {
	Jac[i][j].reserve(max_size * max_size);	
      }
    }

  }
  
  std::vector < std::vector < double > > Jac_outer(dim);
  std::vector < double > Res_outer(1);

         for(int i = 0; i < dim; i++) {  Jac_outer[i].reserve(max_size); }
  // ======= Equation, local - END =======


  // ======= Equation, global - BEGIN =======
    RES->zero();
    if(assembleMatrix) JAC->zero();
  // ======= Equation, global - END =======

 
  // ======= Parameters - BEGIN ======= 
  const double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();

  //****** Control ********************************
 double penalty_outside_control_domain_boundary = PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY;       // penalty for zero control outside Gamma_c
 double penalty_dirichlet_bc_u_equal_q = PENALTY_DIRICHLET_BC_U_EQUAL_Q_BOUNDARY;         //penalty for u = q

 double theta_value_outside_fake_element = 0.;
 //**************************************

 // ======= Parameters - END =======
  
    
   
  // ======= Fractional - BEGIN =======
  const double s_frac = _s_frac;

  const double check_limits = 1.;//1./(1. - s_frac); // - s_frac;
   
  //--- quadrature rules -------------------
  constexpr unsigned qrule_i = QRULE_I;
  constexpr unsigned qrule_j = QRULE_J;
  constexpr unsigned qrule_k = QRULE_K;
  //----------------------
  // ======= Fractional - END =======

  
  
   ///@todo to avoid Petsc complaint about out-of-bounds allocation *******************************
//  MatSetOption(static_cast< PetscMatrix* >(JAC)->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
   //************************************
     
  
  
    
      // ***************** ADD PART - BEGIN  *******************

    if ( IS_CTRL_FRACTIONAL_SOBOLEV ) {
  
     ctrl::Gamma_control_equation_fractional_sobolev_differentiability_index<
                femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES, 
                femus::ctrl:: square_or_cube:: pure_boundary< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >
                >
::control_eqn_bdry(iproc,
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
                    AbsDetJxWeight_iqp_bdry,
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
                    //-----------
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
                    AbsDetJxWeight_iqp_bdry,
                     //-----------
                    phi_ctrl_bdry,
                    phi_ctrl_x_bdry, ///@todo this should be done for all components of ctrl
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
 
 
 
 
 
 
 
 
 
 
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

  // geometry *****************************
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element_iel.geom_type();
      
      geom_element_iel.set_elem_center_3d(iel, solType_coords);
  // geometry end *****************************
  
  // equation *****************************
    unsigned nDofsV = msh->GetElementDofNumber(iel, SolFEType_Mat[state_pos_begin]);
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType_Mat[state_pos_begin + press_type_pos]);
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel,SolFEType_Mat[adj_pos_begin]);
    unsigned nDofsPadj = msh->GetElementDofNumber(iel,SolFEType_Mat[adj_pos_begin + press_type_pos]);

    unsigned nDofsGctrl     = msh->GetElementDofNumber(iel, SolFEType_Mat[ctrl_pos_begin]);
    unsigned nDofsThetactrl = msh->GetElementDofNumber(iel, SolFEType_Mat[theta_index] );

    unsigned nDofsVP = dim * nDofsV + nDofsP; //same for state and adjoint
    unsigned nDofsVPctrl = dim * nDofsGctrl + nDofsThetactrl; //control
   
    unsigned nDofsVP_tot = 2 * nDofsVP + (nDofsVPctrl);
  // equation end *****************************
  
  //***** set target domain flag ********************************** 
   int target_flag = 0;
       target_flag = ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
   //***************************************       
   
 //************ set control flag - BEGIN *********************
        const bool does_iel_contain_Gamma_c = ctrl::square_or_cube :: pure_boundary< ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >  ::volume_elem_contains_a_Gamma_control_face( sol, msh, iel );
  int does_iel_contain_a_bdry_control_face = does_iel_contain_Gamma_c? 1 : 0;

 //************ initialize control node flag: for each Volume Elem, tell me if we have a Boundary Control dof *********************
    std::vector< std::vector<int> > control_node_flag_iel_jface(n_components_ctrl);
	    for(unsigned idim=0; idim < control_node_flag_iel_jface.size(); idim++) {
	          control_node_flag_iel_jface[idim].resize(nDofsGctrl);
    std::fill(control_node_flag_iel_jface[idim].begin(), control_node_flag_iel_jface[idim].end(), 0);
	    }
 //****************************- END *********************** 
  
   //###########################  fake_iel_flag - BEGIN ########################################  
    unsigned int fake_iel_flag = 0;
    unsigned int global_row_index_bdry_constr = pdeSys->KKoffset[SolPdeIndex[theta_index]][iproc];
  for (unsigned  k = 0; k < n_unknowns; k++) {
	unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType_Mat[k]);	//nDofs_V,P_of_st,adj,ctrl
	Sol_n_el_dofs_Mat_vol[k] = ndofs_unk;
	Sol_eldofs_Mat[k].resize(ndofs_unk);	//sol_V,P_of_st,adj,ctrl
	L2G_dofmap_Mat[k].resize(ndofs_unk); 
    for (unsigned i = 0; i < ndofs_unk; i++) {
	      unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType_Mat[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
	  Sol_eldofs_Mat[k][i] = (*sol->_Sol[SolIndex_Mat[k]])(solDof);      // global extraction and local storage for the solution
	  L2G_dofmap_Mat[k][i] = pdeSys->GetSystemDof(SolIndex_Mat[k], SolPdeIndex[k], i, iel);
 
    if (k == SolPdeIndex[theta_index] && L2G_dofmap_Mat[k][i] == global_row_index_bdry_constr) {       fake_iel_flag = iel;  }
    }
  }
  //############################ fake_iel_flag - END #######################################
    
  //************ fake theta flag: - BEGIN this flag tells me what degrees of freedom of the current element are fake for the theta variable  *********************
    std::vector<int>  bdry_int_constr_pos_vec(1, global_row_index_bdry_constr); /*KKoffset[SolPdeIndex[PADJ]][iproc]*/
    std::vector<int> fake_theta_flag(nDofsThetactrl, 0);
    for (unsigned i = 0; i < nDofsThetactrl; i++) {
      if ( L2G_dofmap_Mat[ SolPdeIndex[theta_index] ] [i] == bdry_int_constr_pos_vec[0]) { 	fake_theta_flag[i] = 1;       }
    }
 //************ fake theta flag - END *********************

// setting Jac and Res to zero  - BEGIN ******************************* 
    for(int ivar=0; ivar<n_unknowns; ivar++) {
              Res[SolPdeIndex[ivar]].resize(Sol_n_el_dofs_Mat_vol[ivar]);
     std::fill(Res[SolPdeIndex[ivar]].begin(), Res[SolPdeIndex[ivar]].end(), 0.);         
//       memset(&Res[SolPdeIndex[ivar]][0],0.,Sol_n_el_dofs_Mat_vol[ivar]*sizeof(double));
    }
   
    for(int ivar = 0; ivar < n_unknowns; ivar++) {
      for(int jvar = 0; jvar < n_unknowns; jvar++) {
	    if(assembleMatrix) {  //MISMATCH
                  Jac[ SolPdeIndex[ivar] ] [SolPdeIndex[jvar] ].resize( Sol_n_el_dofs_Mat_vol[ivar] * Sol_n_el_dofs_Mat_vol[jvar] );
     std::fill(Jac[ SolPdeIndex[ivar] ] [SolPdeIndex[jvar] ].begin(), Jac[ SolPdeIndex[ivar] ] [SolPdeIndex[jvar] ].end(), 0.);         
// 		  memset(&Jac[ SolPdeIndex[ivar] ] [SolPdeIndex[jvar] ][0], 0., Sol_n_el_dofs_Mat_vol[ivar]*Sol_n_el_dofs_Mat_vol[jvar]*sizeof(double));
           }
        }
     }
     
    for(int ivar = 0; ivar < dim; ivar++)     std::fill(Jac_outer[ivar].begin(), Jac_outer[ivar].end(), 0.); //did not use Jac_outer as Jac itself was placing the values as expected
    Res_outer[0] = 0.;
// setting Jac and Res to zero  - END ******************************* 

  
  
//======== BoundaryLoop - BEGIN  =====================================================================

  // Perform face loop over elements that contain some control face
  if (does_iel_contain_a_bdry_control_face == 1) {
	  
      double tau = 0.;
      std::vector <double> normal_iqp(dim_offset_grad, 0);
	       
	  // loop on faces of the current element

      for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
          
//-------          
       const unsigned ielGeom_bd = msh->GetElementFaceType(iel, jface);    
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();
//-------          

	    // look for boundary faces
            const int bdry_index = el->GetFaceElementIndex(iel, jface);
            
	    if( bdry_index < 0) {
	      unsigned int face_in_rectangle_domain = -( msh->GetMeshElements()->GetFaceElementIndex(iel, jface) + 1);

// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"u_0",tau,face,0.) && tau!=0.){
          if(  face_in_rectangle_domain == FACE_FOR_CONTROL) { //control face
              
//=================================================== 
		   //we use the dirichlet flag to say: if dirichlet == true, we set 1 on the diagonal. if dirichlet == false, we put the boundary equation
		  std::vector<bool> is_bc_for_control_dirichlet_on_jface(n_components_ctrl);
		  for(unsigned idim = 0; idim < is_bc_for_control_dirichlet_on_jface.size(); idim++) {
		      is_bc_for_control_dirichlet_on_jface[idim] = ml_sol->GetBdcFunctionMLProb()(& ml_prob, geom_element_iel.get_elem_center_bdry_3d(), ctrl_name[idim].c_str(), tau, face_in_rectangle_domain, 0.);
		  }
	  
	
//========= initialize gauss quantities on the boundary ============================================
	 std::vector < double >                      SolVAR_bd_qp(n_unknowns);
	 std::vector < std::vector < double > >       gradSolVAR_bd_qp(n_unknowns);
		for(int k=0; k<n_unknowns; k++) {  gradSolVAR_bd_qp[k].resize(dim_offset_grad /*space_dim*/);  }

//========= gauss_loop boundary===============================================================
		  for(unsigned iqp_bdry = 0; iqp_bdry < ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussPointsNumber(); iqp_bdry++) {
    
    elem_all[qrule_i][ielGeom_bd][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), iqp_bdry, Jac_iqp_bdry, JacI_iqp_bdry, detJac_iqp_bdry, space_dim);
	elem_all[qrule_i][ielGeom_bd][solType_coords]->compute_normal(Jac_iqp_bdry, normal_iqp);
    
    AbsDetJxWeight_iqp_bdry = detJac_iqp_bdry * ml_prob.GetQuadratureRule(ielGeom_bd).GetGaussWeightsPointer()[iqp_bdry];

    elem_all[qrule_i][ielGeom_bd][SolFEType_Mat[theta_index]] ->shape_funcs_current_elem(iqp_bdry, JacI_iqp_bdry,phi_bd_gss_fe[SolFEType_Mat[theta_index]],phi_x_bd_gss_fe[SolFEType_Mat[theta_index]], boost::none, space_dim);
    
    for (unsigned ldim = 0; ldim < dim; ldim++) {
    elem_all[qrule_i][ielGeom_bd][ SolFEType_Mat[ldim + ctrl_pos_begin]] ->shape_funcs_current_elem(iqp_bdry, JacI_iqp_bdry,
                                                                                                    phi_bd_gss_fe[   SolFEType_Mat[ldim + ctrl_pos_begin] ],
                                                                                                    phi_x_bd_gss_fe[ SolFEType_Mat[ldim + ctrl_pos_begin] ], boost::none , space_dim);
    }
                

    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv_vol_at_bdry_new(geom_element_iel.get_coords_at_dofs_3d(), iqp_bdry, jface, Jac_iqp/*not_needed_here*/, JacI_iqp, detJac_iqp/*not_needed_here*/, space_dim);
    
    for (unsigned ldim = 0; ldim < dim; ldim++) {
    elem_all[qrule_i][ielGeom][SolFEType_Mat[ldim + adj_pos_begin]]->shape_funcs_vol_at_bdry_current_elem(iqp_bdry, jface, JacI_iqp, 
                                                                                                          phi_vol_at_bdry_fe  [ SolFEType_Mat[ldim + adj_pos_begin] ],                                                                                                  phi_x_vol_at_bdry_fe[ SolFEType_Mat[ldim + adj_pos_begin] ], boost::none, space_dim);
    }
     

//========== compute gauss quantities on the boundary - BEGIN ===============================================

//=============== Boundary control - BEGIN ========================================= 
		    for (unsigned  kdim = 0; kdim < n_components_ctrl; kdim++) {
                
			    unsigned int ctrl_index = kdim + ctrl_pos_begin;
                
					  SolVAR_bd_qp[ SolPdeIndex[ctrl_index] ] = 0.;
			  for(unsigned ivar2 = 0; ivar2 < dim_offset_grad; ivar2++) {  gradSolVAR_bd_qp[ SolPdeIndex[ctrl_index] ][ivar2] = 0.; }
	  
	         const unsigned ndof_bdry = msh->GetElementFaceDofNumber(iel, jface, SolFEType_Mat[ctrl_index]);

			  for(int i_bd = 0; i_bd < ndof_bdry; i_bd++) {
		                  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
                          
                       SolVAR_bd_qp[SolPdeIndex[ctrl_index]]            += phi_bd_gss_fe  [ SolFEType_Mat[ctrl_index] ][i_bd]                  * Sol_eldofs_Mat[ SolPdeIndex[ ctrl_index ] ][i_vol];

			      for(unsigned ivar2 = 0; ivar2 < dim_offset_grad; ivar2++) { 
                      gradSolVAR_bd_qp[SolPdeIndex[ctrl_index]][ivar2]  += phi_x_bd_gss_fe[ SolFEType_Mat[ctrl_index] ][i_bd * dim_offset_grad + ivar2] * Sol_eldofs_Mat[ SolPdeIndex[ ctrl_index ] ][i_vol];
                        }
                        
			  }
			  
		    }
//=============== Boundary control - END ========================================= 
 
//=============== grad dot n - BEGIN ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
		for(unsigned ldim = 0; ldim < dim; ldim++) {   sol_adj_x_vol_at_bdry_gss[ldim].resize(dim_offset_grad); }
		grad_adj_dot_n_res_qp.resize(dim);
		grad_adj_dot_n_jac_qp.resize(dim);
		for(unsigned ldim=0; ldim<dim; ldim++) {   std::fill(sol_adj_x_vol_at_bdry_gss[ldim].begin(), sol_adj_x_vol_at_bdry_gss[ldim].end(), 0.);  }
		
		for(unsigned ldim=0; ldim<dim; ldim++) {  
		  for (int iv = 0; iv < nDofsVadj; iv++)  {
                     for (int d = 0; d < dim_offset_grad /*space_dim*/; d++) {
			   sol_adj_x_vol_at_bdry_gss[ldim][d] += Sol_eldofs_Mat[SolPdeIndex[ldim + adj_pos_begin]][iv] * phi_x_vol_at_bdry_fe[SolFEType_Mat[ldim + adj_pos_begin]][iv * dim_offset_grad /*space_dim*/ + d];//notice that the convention of the orders x y z is different from vol to bdry
                    }
		  }  
		      
		  grad_adj_dot_n_res_qp[ldim] = 0.;
		  for(unsigned d=0; d < dim_offset_grad /*space_dim*/; d++) {
		      grad_adj_dot_n_res_qp[ldim] += sol_adj_x_vol_at_bdry_gss[ldim][d] * normal_iqp[d];  
		  }
		}
		
//=============== grad dot n - END =========================================       
		  
//========== compute gauss quantities on the boundary - END ================================================

		
//=============== construct control node flag - BEGIN =========================================    
//this is all based on jface only!!!

	      /* (control_node_flag_iel_jface)       picks nodes on \Gamma_c
	         (1 - control_node_flag_iel_jface)   picks nodes on \Omega \setminus \Gamma_c
	       */
		    for(unsigned idim = 0; idim < n_components_ctrl; idim++) {
                
                unsigned int ctrl_index = idim + ctrl_pos_begin;

			if (is_bc_for_control_dirichlet_on_jface[idim] == false) {
                
	         const unsigned ndof_bdry = msh->GetElementFaceDofNumber(iel, jface, SolFEType_Mat[ctrl_index]);
             
		for(unsigned i_bdry = 0; i_bdry < ndof_bdry; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
                control_node_flag_iel_jface[idim][i_vol] = 1;
 			    }
		    }
		
        }
//=============== construct control node flag - END =========================================    
        

        const unsigned ndofs_bdry_max = msh->GetElementFaceDofNumber(iel, jface, solType_coords);
        const unsigned nDofsGctrl_bdry = msh->GetElementFaceDofNumber(iel, jface,  SolFEType_Mat[ctrl_pos_begin]);
        

            for(unsigned i_bdry = 0; i_bdry < ndofs_bdry_max; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

//============ Boundary-Boundary and Boundary-Volume Residuals - BEGIN ============================================================================================
		  
		      for (unsigned  kdim = 0; kdim < dim; kdim++) {
			
			    double lap_res_dctrl_ctrl_bd_kdim = 0.;
                 if (i_bdry < nDofsGctrl_bdry) {
			      for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
				  lap_res_dctrl_ctrl_bd_kdim += gradSolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim] *
				                           phi_x_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry * dim_offset_grad + jdim /*i_bdry + jdim*ndofs_bdry_max*/];
			      }//jdim
                 }
			
/*delta_state row */	   if (i_vol < nDofsV)      Res[kdim]                 [i_vol]  += - control_node_flag_iel_jface[kdim][i_vol] * penalty_dirichlet_bc_u_equal_q * (Sol_eldofs_Mat[SolPdeIndex[kdim + state_pos_begin]][i_vol] - Sol_eldofs_Mat[SolPdeIndex[kdim + ctrl_pos_begin]][i_vol]);	    //u-g
/*delta_adjoint row */     if (i_vol < nDofsVadj)   Res[kdim + adj_pos_begin] [i_vol]  += 0.;	   
/*delta_control row */     if (i_vol < nDofsGctrl)  Res[kdim + ctrl_pos_begin][i_vol]  += - control_node_flag_iel_jface[kdim][i_vol] * AbsDetJxWeight_iqp_bdry * (
                                                                                          _is_block_dctrl_ctrl_inside_main_big_assembly * alpha * SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_bd_gss_fe[SolFEType_Mat[kdim +  ctrl_pos_begin]][i_bdry]
                                                                                        + _is_block_dctrl_ctrl_inside_main_big_assembly * beta * lap_res_dctrl_ctrl_bd_kdim
                                                                                        - IRe * grad_adj_dot_n_res_qp[kdim]  * phi_bd_gss_fe[SolFEType_Mat[kdim +  ctrl_pos_begin]][i_bdry]
                                                                                        - /*(*sol->_Sol[SolIndex_Mat[theta_index]])(0)*/solTheta * phi_bd_gss_fe[SolFEType_Mat[kdim +  ctrl_pos_begin]][i_bdry] * normal_iqp[kdim]      //*sol->_Sol[SolIndex_Mat[theta_index]])(0) finds the global value from KKDof pos(72, 169,etc), SolVAReldof_theta gets the value in the boundary point which will be zero. Theta is just a const
                                                                                        );

		      }//kdim  

//============ Boundary-Boundary and Boundary-Volume Residuals - END ==================================================================================================

//============ Boundary-Boundary Jacobian i-BDRY/j-BDRY - BEGIN  ==================================================================================================
		      for(unsigned j_bdry = 0; j_bdry < ndofs_bdry_max; j_bdry ++) {
			  unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);

				std::vector < double > lap_jac_dctrl_ctrl_bd(dim, 0.);
                    if (i_bdry < nDofsGctrl_bdry && j_bdry < nDofsGctrl_bdry) {

			     for (unsigned kdim = 0; kdim < dim; kdim++) {
				for (unsigned ldim = 0; ldim < dim_offset_grad /*space_dim*/; ldim++) {
///@todo the error is in here
                    lap_jac_dctrl_ctrl_bd[kdim] += phi_x_bd_gss_fe[ SolFEType_Mat[kdim + ctrl_pos_begin] ][i_bdry * dim_offset_grad + ldim] * 
                                                   phi_x_bd_gss_fe[ SolFEType_Mat[kdim + ctrl_pos_begin] ][j_bdry * dim_offset_grad + ldim];
                                        }
                                    }
                    }
                    
                 
			  for (unsigned  kdim = 0; kdim < dim; kdim++) {
                  
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
			    if(i_vol < nDofsV && j_vol < nDofsV && i_vol == j_vol)       		          Jac[kdim][kdim][i_vol * nDofsV + j_vol]	  +=  (1.) * penalty_dirichlet_bc_u_equal_q * control_node_flag_iel_jface[kdim][i_vol];  //u
			 
//BLOCK delta_state - control------------------------------------------------------------------------------------
			    if(i_vol < nDofsV && j_vol < nDofsGctrl && i_vol == j_vol) 	Jac[kdim][kdim + ctrl_pos_begin][i_vol * nDofsGctrl + j_vol]  += (-1.) * penalty_dirichlet_bc_u_equal_q * control_node_flag_iel_jface[kdim][i_vol];  //-g

//DIAG BLOCK delta_control - control  --------------------------------------------------------------------------------------
			  if(i_vol < nDofsGctrl && j_vol < nDofsGctrl) {
				      Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i_vol * nDofsGctrl + j_vol] +=   control_node_flag_iel_jface[kdim][i_vol] * AbsDetJxWeight_iqp_bdry * (
                                                                                                       _is_block_dctrl_ctrl_inside_main_big_assembly * alpha * phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin] ][i_bdry] * phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin] ][j_bdry]
                                                                                                    + _is_block_dctrl_ctrl_inside_main_big_assembly * beta * lap_jac_dctrl_ctrl_bd[kdim]
                                                                                                 );
                  } 
                  
               }//kdim
			
                   
            }//end j_bdry loop
//============ Boundary-Boundary Jacobian i-BDRY/j-BDRY - END  ==================================================================================================
		    
//============ Boundary-Volume Jacobian i-BDRY/j-VOL - BEGIN  ==================================================================================================
//BLOCK delta_control - adjoint------------------------------------------------------------------------------------------------
//===================loop over j in the VOLUME (while i is in the boundary)	      
		for(unsigned j = 0; j < nDofsVadj; j ++) {
            
			//=============== grad dot n  =========================================    
		    for(unsigned ldim = 0; ldim < dim; ldim++) {
			    grad_adj_dot_n_jac_qp[ldim] = 0.;
			    for(unsigned d = 0; d < dim_offset_grad /*space_dim*/; d++) {
				grad_adj_dot_n_jac_qp[ldim] += phi_x_vol_at_bdry_fe[SolFEType_Mat[ldim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + d] * normal_iqp[d];  //notice that the convention of the orders x y z is different from vol to bdry
			    }
		    }
			  //=============== grad dot n  =========================================    

			  for (unsigned kdim = 0; kdim < dim; kdim++) {
				Jac[kdim + ctrl_pos_begin][kdim + adj_pos_begin][i_vol*nDofsVadj + j] += control_node_flag_iel_jface[kdim][i_vol] * (-1.) * (AbsDetJxWeight_iqp_bdry  * phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * IRe * grad_adj_dot_n_jac_qp[kdim]);    		      
			  }
		} // end j loop for volume 
//============ Boundary-Volume Jacobian i-BDRY/j-VOL - END  ==================================================================================================

		    }//end i_bdry loop


		    
//============ Boundary Integral Constraint Residual - BEGIN ============================================================================================
	  for (unsigned  kdim = 0; kdim < dim; kdim++) {
// 		for(unsigned i=0; i < nDofsThetactrl; i ++) { avoid because it is an element dof
/*delta_theta row */ 	/* Res[theta_index][i]*/ Res_outer[0] +=  /*fake_theta_flag[i] **/ AbsDetJxWeight_iqp_bdry * SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * normal_iqp[kdim] ;
// 		}  
	  }
		  
//============ Boundary Integral Constraint Residual - END ============================================================================================
		
//============ Boundary Integral Constraint Jacobian - BEGIN ============================================================================================
       		for(unsigned i_bdry = 0; i_bdry < ndofs_bdry_max; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
            
 		    for (unsigned  kdim = 0; kdim < dim; kdim++) { 
			  for(unsigned i =0; i < nDofsThetactrl; i ++) {
			    if(i_vol < nDofsGctrl) {
				double temp = AbsDetJxWeight_iqp_bdry * ( phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * normal_iqp[kdim]);
//ROW_BLOCK delta_theta - control -- loop over i in the VOLUME (while j(/i_vol) is in the boundary) -------------------------------------------------------------------------------------------------------------
			      Jac[theta_index][ctrl_pos_begin + kdim][i*nDofsGctrl + i_vol]     += - temp; /*AbsDetJxWeight_iqp_bdry * ( phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * normal_iqp[kdim])*/
//COLUMN_BLOCK delta_control - theta ---- loop over j in the VOLUME (while i(/i_vol) is in the boundary) ---------------------------------------------------------------------------------------------------
			      Jac[ctrl_pos_begin + kdim][theta_index][i_vol*nDofsThetactrl + i] += - control_node_flag_iel_jface[kdim][i_vol] /** phi_bd_gss_fe[SolFEType_Mat[theta_index]][i]*/*temp; /*AbsDetJxWeight_iqp_bdry * ( phi_bd_gss_fe[SolFEType_Mat[kdim + ctrl_pos_begin]][i_bdry] * normal_iqp[kdim]);*/
			    }//endif
			  }// i 
		    }//kdim           
            
            
            }		    
//============ Boundary Integral Constraint Jacobian - END ============================================================================================
		    
		    
                }  //end iqp_bdry loop
	  
             }    //end if control face
	     }  //end if boundary faces
      }  // loop over element faces //jface   
  } //end if control element flag

//======== BoundaryLoop - END ==== Boundary Residuals  and Jacobians ==================	
    
    ///@todo at this point, the node flag has been filled for ALL faces
    
 
 
 
//======================= VolumeLoop with Integration (and Zero boundary control outside Gamma_c) - BEGIN =====================================================    

for(unsigned iqp = 0; iqp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); iqp++) {
	
	// *** get Jacobian and test function and test function derivatives ***
       // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), iqp, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);
    AbsDetJxWeight_iqp = detJac_iqp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[iqp];
   
     for(int fe=0; fe < NFE_FAMS; fe++) {
    elem_all[qrule_i][ielGeom][fe]->shape_funcs_current_elem(iqp, JacI_iqp,phi_gss_fe[fe],phi_x_gss_fe[fe], boost::none , space_dim);
      }
 
 
//geometry eval at gauss points - BEGIN ********************************
 std::vector < double > coordX_gss(dim, 0.);
 	for(unsigned k = 0; k <  dim; k++) {
	  for(unsigned i = 0; i < geom_element_iel.get_coords_at_dofs_3d()[k].size(); i++) { 
         coordX_gss[k] += geom_element_iel.get_coords_at_dofs_3d()[k][i] * phi_gss_fe[ solType_coords ][i];
      }
    }
//geometry eval at gauss points - END  ********************************

//begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk < n_unknowns; unk++) {
	    SolVAR_qp[unk] = 0.;
	    for(unsigned ivar2 = 0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++) {
           gradSolVAR_qp[unk][ivar2] = 0.; 
	    }
	  
	    for(unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[unk]; i++) {
                        SolVAR_qp[unk] += phi_gss_fe[ SolFEType_Mat[unk] ][i]             * Sol_eldofs_Mat[SolPdeIndex[unk]][i];
		for(unsigned ivar2 = 0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++) {       
		    gradSolVAR_qp[unk][ivar2]  += phi_x_gss_fe[ SolFEType_Mat[unk] ][i * dim_offset_grad /*space_dim*/ + ivar2] * Sol_eldofs_Mat[SolPdeIndex[unk]][i]; 
		}
	    }//ndofsunk
	  
	} //unk 
 //end unknowns eval at gauss points ********************************
	
#if exact_sol_flag == 1
//======= computation of RHS (force and desired velocity) using MMS - BEGIN =============================================== 
//state values-------------------- //non-hom bdry
vector <double>  exact_stateVel(dim);
  mms_lid_driven::value_stateVel(coordX_gss, exact_stateVel);
vector <double>  exact_lap_stateVel(dim);
  mms_lid_driven::laplace_stateVel(coordX_gss, exact_lap_stateVel);
vector < std::vector < double > > exact_grad_stateVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_stateVel[k].resize(dim);
    std::fill(exact_grad_stateVel[k].begin(), exact_grad_stateVel[k].end(), 0.);
}
  mms_lid_driven::gradient_stateVel(coordX_gss,exact_grad_stateVel);

//adjoint values--------------------//hom bdry
vector <double>  exact_adjVel(dim);
  mms_lid_driven::value_adjVel(coordX_gss, exact_adjVel);
vector <double>  exact_lap_adjVel(dim);
  mms_lid_driven::laplace_adjVel(coordX_gss, exact_lap_adjVel);
vector < std::vector < double > > exact_grad_adjVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_adjVel[k].resize(dim);
    std::fill(exact_grad_adjVel[k].begin(), exact_grad_adjVel[k].end(), 0.);
}
  mms_lid_driven::gradient_adjVel(coordX_gss,exact_grad_adjVel);
vector <double> exact_grad_adjPress(dim);
  mms_lid_driven::gradient_adjPress(coordX_gss, exact_grad_adjPress);

//convection terms from delta_state-------------------------------------
vector <double>  exact_conv_u_nabla_u(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_u[k] += exact_grad_stateVel[k][i] * exact_stateVel[i] ; 
    }
}

//convection terms from delta_adjoint-------------------------
vector <double>  exact_conv_u_nabla_uadj(dim,0.);
vector <double>  exact_conv_nabla_uT_uadj(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_uadj[k] += exact_grad_adjVel[k][i] * exact_stateVel[i] ; 
    exact_conv_nabla_uT_uadj[k] += exact_grad_stateVel[i][k] * exact_adjVel[i];
    }
}

//force and desired velocity ---------------------------------------------
vector <double> exactForce(dim,0.);
vector <double> exactVel_d(dim,0.);
for (unsigned k = 0; k < dim; k++){
    exactForce[k] = - IRe * exact_lap_stateVel[k]  + advection_flag * exact_conv_u_nabla_u[k] + exact_grad_adjPress[k];
    exactVel_d[k] =   exact_stateVel[k] + (1./cost_functional_coeff) * (IRe * exact_lap_adjVel[k] - exact_grad_adjPress[k]) 
                    + (1./cost_functional_coeff) * advection_flag * (exact_conv_u_nabla_uadj[k] - exact_conv_nabla_uT_uadj[k]);
}
//======= computation of RHS (force and desired velocity) using MMS - END =============================================== 
#endif



//============ delta_state row - BEGIN ============================================================================================

  for (unsigned i = 0; i < nDofsV; i++) {
// FIRST ROW
	for (unsigned  kdim = 0; kdim < dim; kdim++) { // velocity block row 
                double lap_res_du_u		= 0.; 
                double adv_res_uold_nablauold 	= 0.;
	      for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		    lap_res_du_u 	       += gradSolVAR_qp[SolPdeIndex[kdim]][jdim] * phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + jdim];
          }
	      for (unsigned jdim = 0; jdim < dim; jdim++) {
		   adv_res_uold_nablauold  += SolVAR_qp[SolPdeIndex[jdim]]  * gradSolVAR_qp[SolPdeIndex[kdim]][jdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i];
	      }      
	      Res[kdim][i]   +=  (           
#if exact_sol_flag == 0
                                         + force[kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]
 #endif                                      
 #if exact_sol_flag == 1
                                       + exactForce[kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]
 #endif
                                       - IRe*lap_res_du_u 
                                       - advection_flag * adv_res_uold_nablauold 
                                       + SolVAR_qp[SolPdeIndex[press_type_pos]] * phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + kdim]) * AbsDetJxWeight_iqp; 
	}	    
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
	for (unsigned j = 0; j < nDofsV; j++) {
		      std::vector < double > lap_jac_du_u(dim,0.);
		      std::vector < double > adv_uold_nablaunew(dim,0.);
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		    lap_jac_du_u[kdim] += phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + jdim]*phi_x_gss_fe[SolFEType_Mat[kdim]][j * dim_offset_grad/*space_dim*/ + jdim];
            }
          }
        for (unsigned  kdim = 0; kdim < dim; kdim++) { 
         for (unsigned  jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
		    adv_uold_nablaunew[kdim] 	 += SolVAR_qp[SolPdeIndex[jdim]] * phi_x_gss_fe[ SolFEType_Mat[kdim] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i];
                }  //jdim
	      }
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
		Jac[kdim][kdim][i*nDofsV + j] += (   IRe * lap_jac_du_u[kdim] 
                                            + advection_flag * adv_uold_nablaunew[kdim] 		// c(u_old, u_new, delta_lambda)
                                            + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]	 // c(u_new,u_old,delta_lambda) diagonal blocks  ..... unew_nablauold
                                         ) * AbsDetJxWeight_iqp; 
              unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		Jac[kdim][off_kdim][i*nDofsV + j] += (	advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[off_kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][off_kdim] * phi_gss_fe[ SolFEType_Mat[kdim] ][i]	// c(u_new,u_old,delta_lambda) off-diagonal blocks  ..... unew_nablauold
                                             ) * AbsDetJxWeight_iqp;
	      }
	} //j_du_u loop
     
//BLOCK Pressure
      for (unsigned j = 0; j < nDofsP; j++) {
	    for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim][press_type_pos][i*nDofsP + j] += -( phi_gss_fe[SolFEType_Mat[press_type_pos]][j] * phi_x_gss_fe[SolFEType_Mat[kdim]][i * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
	    }
      }//j_press loop
   }//i_state loop

//DIV_state
  for (unsigned i = 0; i < nDofsP; i++) {
		    double div_u_du_qp =0.;
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      div_u_du_qp += gradSolVAR_qp[SolPdeIndex[kdim]][kdim] ;
      }
      Res[press_type_pos][i]  +=  ( (div_u_du_qp) * phi_gss_fe[SolFEType_Mat[press_type_pos]][i] ) * AbsDetJxWeight_iqp;
      for (unsigned j = 0; j < nDofsV; j++) {
	  for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      Jac[press_type_pos][kdim][i*nDofsV + j] += -( phi_gss_fe[SolFEType_Mat[press_type_pos]][i] * phi_x_gss_fe[SolFEType_Mat[kdim]][j * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
	  }
      } //j loop
   }//i_div_state
//============ delta_state row - END ============================================================================================


    
//============ delta_adjoint row - BEGIN =============================================================================================
  
  for (unsigned i = 0; i < nDofsVadj; i++) {
// SECOND ROW
     for (unsigned kdim = 0; kdim < dim; kdim++) { 
		    double lap_res_dadj_adj 			= 0.;
		    double adv_res_phiadj_nablauold_uadjold 	= 0.;
		    double adv_res_uold_nablaphiadj_uadjold 	= 0.;
	   for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		lap_res_dadj_adj 		             += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]][jdim]*phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim];
       }
	   for (unsigned jdim = 0; jdim < dim; jdim++) {
		adv_res_phiadj_nablauold_uadjold     += phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim]][kdim] 			* SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
		adv_res_uold_nablaphiadj_uadjold     += SolVAR_qp[SolPdeIndex[jdim]]		       * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
	   }
	  Res[kdim + adj_pos_begin][i] += ( 
#if exact_sol_flag == 0
                            - cost_functional_coeff * target_flag * ctrl::square_or_cube :: cost_functional_without_regularization::DesiredTargetVec()[kdim] 			      * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i]
 #endif                                      
 #if exact_sol_flag == 1
                            - cost_functional_coeff * target_flag * exactVel_d[kdim] 			      * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i]
 #endif
                                        + cost_functional_coeff * target_flag * SolVAR_qp[SolPdeIndex[kdim]] * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i]
                                        - IRe*lap_res_dadj_adj
                                        - advection_flag * adv_res_phiadj_nablauold_uadjold
                                        - advection_flag * adv_res_uold_nablaphiadj_uadjold
                                        + SolVAR_qp[SolPdeIndex[press_type_pos + adj_pos_begin]] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim]
                                        ) * AbsDetJxWeight_iqp;
      }
      
//BLOCK delta_adjoint - state------------------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsV; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][kdim][i*nDofsV + j] += ( - cost_functional_coeff * target_flag * phi_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i] * phi_gss_fe[SolFEType_Mat[kdim]][j] 
                                                             + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i]    * phi_x_gss_fe[ SolFEType_Mat[kdim] ][j*dim_offset_grad /*space_dim*/ + kdim] 		* SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u, u_new, lambda_old)  diagonal blocks  ......phiadj_nablaunew_uadjold 
                                                             + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim] ][j] 			* phi_x_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_new, delta_u, lambda_old) diagonal blocks  ......unew_nablaphiadj_uadjold
                                                            ) * AbsDetJxWeight_iqp;
              unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		Jac[kdim + adj_pos_begin][off_kdim][i*nDofsV + j] += (  advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i] * phi_x_gss_fe[ SolFEType_Mat[off_kdim] ][j*dim_offset_grad /*space_dim*/ + kdim]		      * SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u, u_new, lambda_old)  off-diagonal blocks  ......phiadj_nablaunew_uadjold 
                                                              + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType_Mat[off_kdim] ][j] 		  * phi_x_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	   //c(u_new, delta_u, lambda_old) off-diagonal blocks  ......unew_nablaphiadj_uadjold
                                                             ) * AbsDetJxWeight_iqp;
	  }            
     }//j_dadj_u loop


//DIAG BLOCK delta_adjoint - adjoint---------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsVadj; j++) {
		    std::vector < double > lap_jac_dadj_adj(dim,0.);
		    std::vector < double > adv_uold_nablaphiadj_uadjnew(dim, 0.);
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		  lap_jac_dadj_adj[kdim] += phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
            }
          }
       for (unsigned  kdim = 0; kdim < dim; kdim++) { 
	   for (unsigned jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
	     adv_uold_nablaphiadj_uadjnew[kdim]     += SolVAR_qp[SolPdeIndex[jdim]]  * phi_x_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][j] ;
	   }
	  }
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += (   IRe*lap_jac_dadj_adj[kdim] 
                                                                                + advection_flag * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][j]   //c(delta_u, u_old, lambda_new)  diagonal blocks  ......phiadj_nablauold_uadjnew  
                                                                                + advection_flag * adv_uold_nablaphiadj_uadjnew[kdim] 	//c(u_old, delta_u, lambda_new)
                                                                               ) * AbsDetJxWeight_iqp;
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		  Jac[kdim + adj_pos_begin][off_kdim + adj_pos_begin][i*nDofsVadj + j] += ( advection_flag * phi_gss_fe[ SolFEType_Mat[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim]][kdim] * phi_gss_fe[ SolFEType_Mat[off_kdim + adj_pos_begin] ][j]   //c(delta_u, u_old, lambda_new)  off-diagonal blocks  ......phiadj_nablauold_uadjnew   
                                                                                  ) * AbsDetJxWeight_iqp;
	  }
    } //j_dadj_adj loop
      
//BLOCK Pressure_adj
    for (unsigned j = 0; j < nDofsPadj; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][press_type_pos + adj_pos_begin][i*nDofsPadj + j] += -( phi_gss_fe[SolFEType_Mat[press_type_pos + adj_pos_begin]][j] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
	  }
    }//j_press_adj loop
  }//i_adj loop

//DIV_adj
  for (unsigned i = 0; i < nDofsPadj; i++) {
		double div_adj_dadj_qp = 0.;
      for (unsigned kdim = 0; kdim < dim; kdim++) {
	    div_adj_dadj_qp += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin ]][kdim] ;
      }
      Res[press_type_pos + adj_pos_begin][i] += ( (div_adj_dadj_qp) * phi_gss_fe[SolFEType_Mat[press_type_pos + adj_pos_begin]][i] ) * AbsDetJxWeight_iqp;
      for (unsigned j = 0; j < nDofsVadj; j++) {
        for (unsigned kdim = 0; kdim < dim; kdim++) {
            Jac[press_type_pos + adj_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += - ( phi_gss_fe[SolFEType_Mat[press_type_pos + adj_pos_begin]][i] * phi_x_gss_fe[SolFEType_Mat[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
        }
      }//j loop
  }//i_div_adj

//============ delta_adjoint row - END =============================================================================================

//============ delta_control row - BEGIN  ==================================================================================================
// delta_control
    for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {
        
         for (unsigned i = 0; i < nDofsGctrl; i++) {
       Res[kdim + ctrl_pos_begin][i] += - penalty_outside_control_domain_boundary * ( (1 - control_node_flag_iel_jface[kdim][i]) *
                             (  Sol_eldofs_Mat[SolPdeIndex[kdim + ctrl_pos_begin]][i] -  PENALTY_OUTSIDE_CONTROL_DOMAIN_BOUNDARY_VALUE_CONSISTENT_WITH_BOUNDARY_OF_BOUNDARY)  );              //enforce control zero outside the control boundary


// //DIAG BLOCK delta_control - control--------------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsGctrl; j++) {
	    if (i == j) {
		Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsGctrl + j] += penalty_outside_control_domain_boundary * (1 - control_node_flag_iel_jface[kdim][i]);              //enforce control zero outside the control boundary
                  } //end i==j
        }//j_dctrl_ctrl loop
     }//i_ctrl loop
  
      }  //kdim

//============ delta_control row - END  ==================================================================================================
 

 
//============ delta_mu row - BEGIN  ============================================================================================

#if INEQ_FLAG != 0
  //MU
//************ Residual, BEGIN *********************
  for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
          
  for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; i++) {
      
       Res[mu_pos_begin + kdim][i]  +=  (- penalty_outside_control_domain_boundary) *  (1 - control_node_flag_iel_jface[kdim][i]) * (Sol_eldofs_Mat[mu_pos_begin + kdim][i] - 0.);
      
     }
  }
//************ Residual, END *********************

  //MU
//************ Jacobian, BEGIN *********************
  for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
    for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; i++) {
      for (unsigned j = 0; j < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; j++) {
            if (i == j) {
               Jac[mu_pos_begin + kdim][mu_pos_begin + kdim][i * Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim] + j]  +=  penalty_outside_control_domain_boundary * (1 - control_node_flag_iel_jface[kdim][i]);
            }
         }
      }
  }
//************ Jacobian, END *********************
#endif

//============ delta_mu row - END  ============================================================================================
 
 
 
      }  // end quadrature point loop
      
//======================= VolumeLoop with Integration (and Zero boundary control outside Gamma_c) - END =====================================================    

    
    
//======================= Loop without Integration For Fake Dofs related to Compatibility Condition - BEGIN =====================================================    
    
        //============ delta_theta - theta row ==================================================================================================
  for (unsigned i = 0; i < nDofsThetactrl; i++) {
             /* if ( fake_theta_flag[i] != 1 ) */             Res[ theta_index ][i]    = - (1 - fake_theta_flag[i]) * ( theta_value_outside_fake_element - Sol_eldofs_Mat[SolPdeIndex[theta_index]][i]);  // Res_outer for the exact row (i.e. when fakeflag=1 , res =0(use Res_outer) and if not 1 this loop) and this is to take care of fake placement for the rest of dofs of theta values as 8
     for (unsigned j = 0; j < nDofsThetactrl; j++) {
			         if(i == j)  Jac[ theta_index ][ theta_index ][i*nDofsThetactrl + j] = (1 - fake_theta_flag[i]) * 1.; //likewise Jac_outer (actually Jac itself works in the correct placement) for bdry integral and this is for rest of dofs
             }//j_theta loop
        }//i_theta loop
   
 //============ delta_theta row ==================================================================================================
//======================= Loop without Integration For Fake Dofs related to Compatibility Condition - END =====================================================    



 //======================= From local to global - BEGIN =====================================================
    // FIRST ALL THE BLOCKS WITHOUT THETA ROW OR COLUMN - BEGIN 
    for(unsigned i_unk = 0; i_unk < theta_index; i_unk++) {
      RES->add_vector_blocked(Res[SolPdeIndex[i_unk]], L2G_dofmap_Mat[i_unk]);
        for(unsigned j_unk = 0; j_unk < theta_index; j_unk++) {
	  if(assembleMatrix) JAC->add_matrix_blocked( Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ], L2G_dofmap_Mat[i_unk], L2G_dofmap_Mat[j_unk]);
        }
    }
    // FIRST ALL THE BLOCKS WITHOUT THETA ROW OR COLUMN - END 
    
    // THEN THE BLOCKS WITH THETA ROW OR COLUMN - BEGIN
	/*delta_theta-theta*/    JAC->add_matrix_blocked( Jac[ SolPdeIndex[ theta_index ] ][ SolPdeIndex[ theta_index ] ], L2G_dofmap_Mat[ theta_index ], L2G_dofmap_Mat[ theta_index ]);
	    
     if (does_iel_contain_a_bdry_control_face == 1) {
	      for (unsigned kdim = 0; kdim < dim; kdim++) {
// // //                           /*delta_control*/       RES->add_vector_blocked(Res[SolPdeIndex[n_unknowns-2-kdim]], L2G_dofmap_Mat[n_unknowns-2-kdim]); ///@todo why was this here?!?
		if(assembleMatrix) {
                          /*delta_theta-control*/ JAC->add_matrix_blocked( Jac[ SolPdeIndex[ theta_index ] ][ SolPdeIndex[theta_index -1 -kdim] ], bdry_int_constr_pos_vec, L2G_dofmap_Mat[theta_index -1 -kdim]);
                          /*delta_control-theta*/ JAC->add_matrix_blocked( Jac[ SolPdeIndex[theta_index -1 -kdim] ][ SolPdeIndex[ theta_index ] ], L2G_dofmap_Mat[theta_index -1 -kdim], bdry_int_constr_pos_vec); 
		}
      }  //kdim
     }  //add control boundary element contributions
     
     
      if (does_iel_contain_a_bdry_control_face == 1) {
          /*delta_theta(bdry constr)*/         RES->add_vector_blocked(Res_outer, bdry_int_constr_pos_vec);
	  }
	  
     /* if (L2G_dofmap_Mat[n_unknowns-1][0] != bdry_int_constr_pos_vec[0]) */ /*delta_theta(fake)*/          RES->add_vector_blocked( Res[ SolPdeIndex[ theta_index ]],       L2G_dofmap_Mat[ theta_index ]);
    // THEN THE BLOCKS WITH THETA ROW OR COLUMN - END
	  

  
     if (print_algebra_local) {
         
         //extract only the ctrl range
         std::vector <unsigned>   Sol_n_el_dofs_Mat_vol_only_ctrl(1, 0);
	      for (unsigned kdim = 0; kdim < dim; kdim++) {
	       Sol_n_el_dofs_Mat_vol_only_ctrl[0] = Sol_n_el_dofs_Mat_vol[kdim + ctrl_pos_begin];
         assemble_jacobian<double,double>::print_element_residual(iel, Res[kdim + ctrl_pos_begin], Sol_n_el_dofs_Mat_vol_only_ctrl, 10, 5);
         assemble_jacobian<double,double>::print_element_jacobian(iel, Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin], Sol_n_el_dofs_Mat_vol_only_ctrl, 10, 5);
          }
    }
     
 //======================= From local to global - END =====================================================    
   
  } //end list of elements loop for each subdomain
  

#if INEQ_FLAG != 0
  //MU in res ctrl - BEGIN  ***********************************
ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::add_one_times_mu_res_ctrl(iproc,
                               ineq_flag,
                               ctrl_index_in_mat,
                               mu_index_in_mat,
                               SolIndex_Mat,
                               sol,
                               mlPdeSys,
                               pdeSys,
                               RES);
  //MU in res ctrl - END ***********************************
#endif
  
  
  
RES->close();
if (assembleMatrix) JAC->close();  /// This is needed for the parallel, when splitting the add part from the insert part!!!
      // ***************** ADD PART - END  *******************



#if INEQ_FLAG != 0
//   ***************** INSERT PART - BEGIN (must go AFTER the sum, clearly) *******************
    
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

	if ( ctrl:: square_or_cube:: pure_boundary< femus::ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >:: volume_elem_contains_a_Gamma_control_face( sol, msh, iel ) ) {


    	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);


        std::pair< int, unsigned int > pair_control_iface = femus::face_is_a_Gamma_control_face_of_some_index< femus::ctrl:: GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >(msh->GetMeshElements(), iel, iface);

        const int  iface_is_a_boundary_control  = pair_control_iface.first;

	    if( iface_is_a_boundary_control ) {

       ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::update_active_set_flag_for_current_nonlinear_iteration_bdry
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
 

  ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::node_insertion_bdry(iel, iface, 
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
#endif




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
  

  
}


 
    

// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv 
//i.e. ||prol_(u_(i-1)) - u_(i)|| = err(i) => err(i-1)/err(i) = 2^alpha ,implemented as log(err(i)/err(i+1))/log2

static double*  GetErrorNorm(const MultiLevelProblem & ml_prob, MultiLevelSolution* ml_sol, Solution* sol_coarser_prolongated) {
  
    static double ErrorNormArray[ pure_boundary_norms::no_of_l2_norms + pure_boundary_norms::no_of_h1_norms ];
    
  unsigned level = ml_sol->GetMLMesh()->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = ml_sol->GetMLMesh()->GetLevel(level);    // pointer to the mesh (level) object
  elem*     el  = msh->GetMeshElements();  // pointer to the elem object in msh (level)
  Solution* sol = ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

 // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  std::vector < std::vector < double > > coordX(dim);    // local coordinates

  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(max_size);
  }
   
  //geometry *******************************

 // solution variables *******************************************
  const int n_vars_state = dim+1;
  const int n_unknowns = 3*n_vars_state; //(2.*dim)+1; //state , adjoint of velocity terms and one pressure term
  const int press_type_pos = dim;
  const int state_pos_begin = 0;
  const int adj_pos_begin   = dim+1;
  const int ctrl_pos_begin   = 2*(dim+1);
  const int theta_index = press_type_pos + ctrl_pos_begin;
  
  std::vector < std::string > Solname_Mat(n_unknowns);  // const char Solname_Mat[4][8] = {"u_0","u_1","u_2","u_p"};
  Solname_Mat              [state_pos_begin+0] =                "u_0";
  Solname_Mat              [state_pos_begin+1] =                "u_1";
  if (dim == 3) Solname_Mat[state_pos_begin+2] =                "u_2";
  Solname_Mat              [state_pos_begin + press_type_pos] = "u_p";
  
  Solname_Mat              [adj_pos_begin + 0] =              "adj_0";
  Solname_Mat              [adj_pos_begin + 1] =              "adj_1";
  if (dim == 3) Solname_Mat[adj_pos_begin + 2] =              "adj_2";
  Solname_Mat              [adj_pos_begin + press_type_pos] = "adj_p";

  Solname_Mat              [ctrl_pos_begin + 0] =              "ctrl_0";
  Solname_Mat              [ctrl_pos_begin + 1] =              "ctrl_1";
  if (dim == 3) Solname_Mat[ctrl_pos_begin + 2] =              "ctrl_2";
  Solname_Mat              [ctrl_pos_begin + press_type_pos] = "theta";
  
  std::vector < unsigned > SolIndex_Mat(n_unknowns);  
  std::vector < unsigned > SolFEType_Mat(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolIndex_Mat[ivar]	= ml_sol->GetIndex        (Solname_Mat[ivar].c_str());
    SolFEType_Mat[ivar]	= ml_sol->GetSolutionType(SolIndex_Mat[ivar]);
  }

  std::vector < double > Sol_n_el_dofs(n_unknowns);
  
  //==========================================================================================
  // velocity ************************************
  std::vector < std::vector < double > > phi_gss_fe(NFE_FAMS);
  std::vector < std::vector < double > > phi_x_gss_fe(NFE_FAMS);
  std::vector < std::vector < double > > phi_xx_gss_fe(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(max_size);
      phi_x_gss_fe[fe].reserve(max_size*dim);
     phi_xx_gss_fe[fe].reserve(max_size*(3*(dim-1)));
   }
  
  //=================================================================================================
  
  // quadratures ********************************
  double AbsDetJxWeight_iqp;
  
  
  //----------- dofs ------------------------------
  std::vector < std::vector < double > > Sol_eldofs_Mat(n_unknowns);
  std::vector < std::vector < double > > gradSol_eldofs_Mat(n_unknowns);
  
  std::vector < std::vector < double > > SolVAR_coarser_prol_eldofs(n_unknowns);
  std::vector < std::vector < double > > gradSolVAR_coarser_prol_eldofs(n_unknowns);


  for(int k = 0; k < n_unknowns; k++) {
    Sol_eldofs_Mat[k].reserve(max_size);
    gradSol_eldofs_Mat[k].reserve(max_size*dim); 
    
    SolVAR_coarser_prol_eldofs[k].reserve(max_size);
    gradSolVAR_coarser_prol_eldofs[k].reserve(max_size*dim);    
  }

  //------------ at quadrature points ---------------------
  std::vector < double > SolVAR_qp(n_unknowns);
  std::vector < double > SolVAR_coarser_prol_qp(n_unknowns);
  std::vector < std::vector < double > > gradSolVAR_qp(n_unknowns);
  std::vector < std::vector < double > > gradSolVAR_coarser_prol_qp(n_unknowns);
  for(int k = 0; k < n_unknowns; k++) {
      gradSolVAR_qp[k].reserve(max_size);  
      gradSolVAR_coarser_prol_qp[k].reserve(max_size);  
  }
      
   std::vector < double > l2norm ( pure_boundary_norms::no_of_l2_norms ,0.);
  std::vector < double > seminorm ( pure_boundary_norms::no_of_h1_norms ,0.);

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    
  // geometry *****************************
    short unsigned ielGeom = msh->GetElementType(iel);
    
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    for (unsigned  k = 0; k < dim; k++) {       coordX[k].resize(nDofsX);    }
  
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->GetTopology()->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }
    
      // elem average point 
    std::vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
  for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofsX; i++) {
         elem_center[j] += coordX[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofsX; }
  //*************************************** 
  
  // geometry end *****************************
  
  
 // equation *****************************
    unsigned nDofsV = msh->GetElementDofNumber(iel, SolFEType_Mat[state_pos_begin]);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType_Mat[state_pos_begin + press_type_pos]);    // number of solution element dofs
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel,SolFEType_Mat[adj_pos_begin]);    // number of solution element dofs
    unsigned nDofsPadj = msh->GetElementDofNumber(iel,SolFEType_Mat[adj_pos_begin + press_type_pos]);    // number of solution element dofs

   unsigned nDofsGctrl = msh->GetElementDofNumber(iel,SolFEType_Mat[ctrl_pos_begin]);
    unsigned nDofsThetactrl = msh->GetElementDofNumber(iel,SolFEType_Mat[theta_index] );

    unsigned nDofsVP = dim * nDofsV + nDofsP; //same for state and adjoint
    unsigned nDofsVPctrl = dim * nDofsGctrl + nDofsThetactrl; //control
   
    unsigned nDofsVP_tot = 2*nDofsVP + (nDofsVPctrl);
  // equation end *****************************


   //STATE###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType_Mat[k]);
	Sol_n_el_dofs[k]=ndofs_unk;
       Sol_eldofs_Mat[k].resize(ndofs_unk);
       SolVAR_coarser_prol_eldofs[k].resize(ndofs_unk);
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType_Mat[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
       Sol_eldofs_Mat[k][i] = (*sol->_Sol[SolIndex_Mat[k]])(solDof);      // global extraction and local storage for the solution
       SolVAR_coarser_prol_eldofs[k][i] = (*sol_coarser_prolongated->_Sol[SolIndex_Mat[k]])(solDof);      // global extraction and local storage for the solution
      }
    }
  //CTRL###################################################################

 
      // ********************** Gauss point loop *******************************
      for(unsigned i_qp = 0; i_qp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); i_qp++) {
	
 
      for(int fe=0; fe < NFE_FAMS; fe++) {
	msh->_finiteElement[ielGeom][fe]->Jacobian(coordX, i_qp, AbsDetJxWeight_iqp,phi_gss_fe[fe],phi_x_gss_fe[fe],phi_xx_gss_fe[fe]);
      }
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
  	msh->_finiteElement[ielGeom][CONTINUOUS_BIQUADRATIC]->Jacobian(coordX,i_qp,AbsDetJxWeight_iqp,phi_gss_fe[CONTINUOUS_BIQUADRATIC],phi_x_gss_fe[CONTINUOUS_BIQUADRATIC],phi_xx_gss_fe[CONTINUOUS_BIQUADRATIC]);

 //begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk < n_unknowns; unk++) {
	  SolVAR_qp[unk] = 0.;
	  SolVAR_coarser_prol_qp[unk] = 0.;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_qp[unk][ivar2] = 0.; 
	    gradSolVAR_coarser_prol_qp[unk][ivar2] = 0.; 
	  }
    }
	  
	for(unsigned unk = 0; unk <  n_unknowns; unk++) {
	  for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
	    SolVAR_qp[unk] += phi_gss_fe[ SolFEType_Mat[unk] ][i] * Sol_eldofs_Mat[unk][i];
	    SolVAR_coarser_prol_qp[unk] += phi_gss_fe[ SolFEType_Mat[unk] ][i] * SolVAR_coarser_prol_eldofs[unk][i];
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType_Mat[unk] ][i*dim+ivar2] * Sol_eldofs_Mat[unk][i]; 
	      gradSolVAR_coarser_prol_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType_Mat[unk] ][i*dim+ivar2] * SolVAR_coarser_prol_eldofs[unk][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************


	for(unsigned unk = 0; unk < n_unknowns; unk++) {
        l2norm[unk] += ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * AbsDetJxWeight_iqp ; 
        
     }
    
    
	for(unsigned unk = 0; unk < dim; unk++) {
        for(int j = 0; j < dim; j++){
        seminorm[unk] += (gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * ( gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * AbsDetJxWeight_iqp ;
        seminorm[unk + dim] += (gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * ( gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * AbsDetJxWeight_iqp ;
        seminorm[unk + 2*dim] += (gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * ( gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * AbsDetJxWeight_iqp ;
        
        }
     }
    
    } // end gauss point loop
  } //end element loop for each process


    // add the norms of all processes
  NumericVector* norm_vec_inexact;
  norm_vec_inexact = NumericVector::build().release();
  norm_vec_inexact->init(msh->n_processors(), 1 , false, AUTOMATIC);

	for(unsigned unk = 0; unk <  pure_boundary_norms::no_of_l2_norms ; unk++) {
        norm_vec_inexact->set(iproc, l2norm[unk]);
        norm_vec_inexact->close();
        l2norm[unk] = norm_vec_inexact->l1_norm();
    }

	for(unsigned unk = 0; unk <  pure_boundary_norms::no_of_h1_norms ; unk++) {
        norm_vec_inexact->set(iproc, seminorm[unk]);
        norm_vec_inexact->close();
        seminorm[unk] = norm_vec_inexact->l1_norm();
    }


  delete norm_vec_inexact;
  
 
	for(unsigned unk = 0; unk <  pure_boundary_norms::no_of_l2_norms ; unk++) {
        ErrorNormArray[unk] = sqrt(l2norm[unk]);
    }
 	for(unsigned unk = 0; unk <  pure_boundary_norms::no_of_h1_norms ; unk++) {
        ErrorNormArray[unk +  pure_boundary_norms::no_of_l2_norms ] = sqrt(seminorm[unk]);
    }
  
   return ErrorNormArray;
  
  
}


    
    
    
};




class lifting_internal: public femus::lifting_internal  {

public: 
    
    
 //Unknown definition  - BEGIN ==================
   enum pos_vector_quantities {pos_index_state = 0, pos_index_adj, pos_index_ctrl, pos_index_mu};
 
static const std::vector< unsigned >  provide_state_adj_ctrl_mu_offsets(const unsigned int dimension) {
     
     std::vector<unsigned>  opt_control_offsets(4);  
     
     opt_control_offsets[pos_index_state] = 0; 
     opt_control_offsets[pos_index_adj]   = dimension + 1; 
     opt_control_offsets[pos_index_ctrl]    = 2 * (dimension + 1);
     opt_control_offsets[pos_index_mu]    = 3 * (dimension + 1); 

     return opt_control_offsets;
     
 } 
 




static  void  name_of_unknowns(std::vector< Unknown > & unknowns, const unsigned int dimension) {

 
  std::vector< unsigned >  vector_offsets = provide_state_adj_ctrl_mu_offsets(dimension);


const int state_pos_begin   =  vector_offsets[pos_index_state];
  const int adj_pos_begin   =  vector_offsets[pos_index_adj];
 const int ctrl_pos_begin   =  vector_offsets[pos_index_ctrl];
   const int mu_pos_begin   =  vector_offsets[pos_index_mu];
   
  
                        unknowns[state_pos_begin + 0]._name    = "u_0";
                        unknowns[state_pos_begin + 1]._name    = "u_1";
  if (dimension == 3)   unknowns[state_pos_begin + 2]._name    = "u_2";
                                unknowns[dimension]._name      = "p_u";
  
                        unknowns[adj_pos_begin + 0]._name      = "adj_0";
                        unknowns[adj_pos_begin + 1]._name      = "adj_1";
  if (dimension == 3)   unknowns[adj_pos_begin + 2]._name      = "adj_2";
                unknowns[adj_pos_begin + dimension]._name      = "p_adj";
  
                       unknowns[ctrl_pos_begin + 0]._name      = "ctrl_0";
                       unknowns[ctrl_pos_begin + 1]._name      = "ctrl_1";
  if (dimension == 3)  unknowns[ctrl_pos_begin + 2]._name      = "ctrl_2";
               unknowns[ctrl_pos_begin + dimension]._name      = "p_ctrl";

                       unknowns[mu_pos_begin + 0]._name      = "mu_0";
                       unknowns[mu_pos_begin + 1]._name      = "mu_1";
  if (dimension == 3)  unknowns[mu_pos_begin + 2]._name      = "mu_2";
 
     
}


 //Unknown definition  ==================
static  const std::vector< Unknown >  provide_list_of_unknowns(const unsigned int dimension) {
     
     
  std::vector< FEFamily > feFamily;
  std::vector< FEOrder >   feOrder;

                        feFamily.push_back(LAGRANGE);   //state
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/);
  
                        feFamily.push_back(LAGRANGE);   //adjoint
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/);
  
                        feFamily.push_back(LAGRANGE);   //control
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
 
                        feFamily.push_back(LAGRANGE);   //mu
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
 
  
  
                        feOrder.push_back(SECOND);
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);
                        feOrder.push_back(FIRST);
  
                        feOrder.push_back(SECOND);
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);
                        feOrder.push_back(FIRST);
  
                        feOrder.push_back(SECOND);
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);
                        feOrder.push_back(FIRST);
 
                        feOrder.push_back(SECOND);   //mu
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Unknown >  unknowns(feFamily.size());
 
    name_of_unknowns(unknowns, dimension);
  
     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._fe_family  = feFamily[u];
              unknowns[u]._fe_order   = feOrder[u];
              unknowns[u]._time_order = 0;
              unknowns[u]._is_pde_unknown = true;
              unknowns[u]._is_sparse = true;
              
     }
     
 
   return unknowns;
     
}


 //Unknown definition  - END ==================
  
    
    

static void assemble_ns_dirichlet_control_lifting_internal_AD(MultiLevelProblem& ml_prob) {

  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled
    
  
  const double cost_functional_coeff = COST_FUNCTIONAL_COEFF;
  const double alpha = ALPHA_CTRL_VOL;
  const double beta  = BETA_CTRL_VOL;

  
  std::cout << " ********************************  AD SYSTEM ******************************************** " << std::endl;
  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem * mlPdeSys   = & ml_prob.get_system<NonLinearImplicitSystem> ("NSOpt");   // pointer to the nonlinear implicit system named "NSOpt" 
   const unsigned level = mlPdeSys->GetLevelToAssemble();
 
  Mesh*          msh          	= ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         	= msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution*  ml_sol    = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        	= ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


   LinearEquationSolver*  pdeSys	 = mlPdeSys->_LinSolver[level];   
 SparseMatrix*    JAC         	= pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  unsigned coordXType = 2; /*CONTINUOUS_BIQUADRATIC*/// get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  unsigned solType_coords = coordXType;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = dim  /*3*/  /*2*/    ;
 
  std::vector < std::vector < double > > coordX(dim);    // local coordinates

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(max_size);
  }
  //geometry *******************************

//STATE######################################################################
  //velocity *******************************
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = ml_sol->GetIndex("u_0");
  solVIndex[1] = ml_sol->GetIndex("u_1");

  if (dim == 3) solVIndex[2] = ml_sol->GetIndex("u_2");

  unsigned solVType = ml_sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
  std::vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("u_0");
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("u_1");

  if (dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("u_2");
  
  std::vector < std::vector < adept::adouble > >  solV(dim);    // local solution
   std::vector < std::vector < adept::adouble > > aResV(dim);    // local redidual vector
   
 for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(max_size);
    aResV[k].reserve(max_size);
  }

  
  std::vector <double> phiV_gss;  // local test function
  std::vector <double> phiV_x_gss; // local test function first order partial derivatives
  std::vector <double> phiV_xx_gss; // local test function second order partial derivatives

  phiV_gss.reserve(max_size);
  phiV_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  phiV_xx_gss.reserve(max_size * dim2);
  

  //velocity *******************************

  //pressure *******************************
  unsigned solPIndex;
  solPIndex = ml_sol->GetIndex("p_u");    // get the position of "p_u" in the ml_sol object
  unsigned solPType = ml_sol->GetSolutionType(solPIndex);    // get the finite element type for "p_u"

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("p_u");    // get the position of "p_u" in the pdeSys object

  std::vector < adept::adouble >  solP; // local solution
  std::vector < adept::adouble > aResP; // local redidual vector
  
  solP.reserve(max_size);
  aResP.reserve(max_size);
  
  double* phiP_gss;
  //pressure *******************************
//STATE######################################################################
  
//ADJOINT######################################################################
  //velocity *******************************
  std::vector < unsigned > solVadjIndex(dim);
  solVadjIndex[0] = ml_sol->GetIndex("adj_0");
  solVadjIndex[1] = ml_sol->GetIndex("adj_1");

  if (dim == 3) solVadjIndex[2] = ml_sol->GetIndex("adj_2");

  unsigned solVadjType = ml_sol->GetSolutionType(solVadjIndex[0]);
 std::vector < unsigned > solVPdeadjIndex(dim);
  solVPdeadjIndex[0] = mlPdeSys->GetSolPdeIndex("adj_0");
  solVPdeadjIndex[1] = mlPdeSys->GetSolPdeIndex("adj_1");

  if (dim == 3) solVPdeadjIndex[2] = mlPdeSys->GetSolPdeIndex("adj_2");
  
  std::vector < std::vector < adept::adouble > >  solVadj(dim);
   std::vector < std::vector < adept::adouble > > aResVadj(dim);
   
 for (unsigned  k = 0; k < dim; k++) {
    solVadj[k].reserve(max_size);
    aResVadj[k].reserve(max_size);
  }

  
  std::vector <double> phiVadj_gss;  // local test function
  std::vector <double> phiVadj_x_gss; // local test function first order partial derivatives
  std::vector <double> phiVadj_xx_gss; // local test function second order partial derivatives

  phiVadj_gss.reserve(max_size);
  phiVadj_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  phiVadj_xx_gss.reserve(max_size * dim2);
  
  //velocity *******************************

  //pressure *******************************
  unsigned solPadjIndex;
  solPadjIndex = ml_sol->GetIndex("p_adj");
  unsigned solPadjType = ml_sol->GetSolutionType(solPadjIndex);    // get the finite element type for "p_adj"

  unsigned solPPdeadjIndex;
  solPPdeadjIndex = mlPdeSys->GetSolPdeIndex("p_adj");    // get the position of "p_adj" in the pdeSys object

  std::vector < adept::adouble >  solPadj; // local solution
  std::vector < adept::adouble > aResPadj; // local redidual vector
  
  solPadj.reserve(max_size);
  aResPadj.reserve(max_size);
  
  double* phiPadj_gss;
  //pressure *******************************
//ADJOINT######################################################################

  
//CONTROL######################################################################
  //velocity *******************************
  std::vector < unsigned > solVctrlIndex(dim);
  solVctrlIndex[0] = ml_sol->GetIndex("ctrl_0");
  solVctrlIndex[1] = ml_sol->GetIndex("ctrl_1");

  if (dim == 3) solVctrlIndex[2] = ml_sol->GetIndex("ctrl_2");

  unsigned solVctrlType = ml_sol->GetSolutionType(solVctrlIndex[0]);
 std::vector < unsigned > solVPdectrlIndex(dim);
 
  solVPdectrlIndex[0] = mlPdeSys->GetSolPdeIndex("ctrl_0");
  solVPdectrlIndex[1] = mlPdeSys->GetSolPdeIndex("ctrl_1");
  if (dim == 3) solVPdectrlIndex[2] = mlPdeSys->GetSolPdeIndex("ctrl_2");
  
  std::vector < std::vector < adept::adouble > >  solVctrl(dim);    // local solution
   std::vector < std::vector < adept::adouble > > aResVctrl(dim);    // local redidual vector
   
 for (unsigned  k = 0; k < dim; k++) {
    solVctrl[k].reserve(max_size);
    aResVctrl[k].reserve(max_size);
  }

  
  std::vector <double> phiVctrl_gss;  // local test function
  std::vector <double> phiVctrl_x_gss; // local test function first order partial derivatives
  std::vector <double> phiVctrl_xx_gss; // local test function second order partial derivatives

  phiVctrl_gss.reserve(max_size);
  phiVctrl_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  phiVctrl_xx_gss.reserve(max_size * dim2);
  

  //velocity *******************************

  //pressure *******************************
  unsigned solPctrlIndex;
  solPctrlIndex = ml_sol->GetIndex("p_ctrl");    // get the position of "p_ctrl" in the ml_sol object
  unsigned solPctrlType = ml_sol->GetSolutionType(solPctrlIndex);    // get the finite element type for "p_ctrl"

  unsigned solPPdectrlIndex;
  solPPdectrlIndex = mlPdeSys->GetSolPdeIndex("p_ctrl");    // get the position of "p_ctrl" in the pdeSys object

  std::vector < adept::adouble >  solPctrl; // local solution
  std::vector < adept::adouble > aResPctrl; // local redidual vector
  
  solPctrl.reserve(max_size);
  aResPctrl.reserve(max_size);
  
  double* phiPctrl_gss;
  //pressure *******************************
//CONTROL######################################################################


  
  //Nondimensional values ******************
  double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  //Nondimensional values ******************
  
  
  
  
  std::vector < int > L2G_dofmap_Mat; // local to global pdeSys dofs
  L2G_dofmap_Mat.reserve(3 *(dim + 1) *max_size);

  std::vector < double > Res; // local redidual vector
  Res.reserve(3 *(dim + 1) *max_size);

  std::vector < double > Jac;
  Jac.reserve(3* (dim + 1) *max_size * 3*(dim + 1) *max_size);

  JAC->zero(); // Set to zero all the entries of the Global Matrix

//*************************************************** 
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
   double AbsDetJxWeight_iqp;
   double detJac_iqp;

    //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
//*************************************************** 
 
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {


   // geometry *****************************
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element_iel.geom_type();
  // geometry end *****************************

  // equation
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsVP = dim * nDofsV + nDofsP;
   
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel, solVadjType);    // number of solution element dofs
    unsigned nDofsPadj = msh->GetElementDofNumber(iel, solPadjType);    // number of solution element dofs

     unsigned nDofsVctrl = msh->GetElementDofNumber(iel, solVctrlType);    // number of solution element dofs
    unsigned nDofsPctrl = msh->GetElementDofNumber(iel, solPctrlType);    // number of solution element dofs
    
    
     unsigned nDofsVP_tot = 3*nDofsVP;
     
     
    
    for (unsigned  k = 0; k < dim; k++)  {
      solV[k].resize(nDofsV);
      solVadj[k].resize(nDofsVadj);
      solVctrl[k].resize(nDofsVctrl);
    }
    solP.resize(nDofsP);
    solPadj.resize(nDofsPadj);
    solPctrl.resize(nDofsPctrl);
    

//element matrices and vectors
    // resize local arrays
    L2G_dofmap_Mat.resize(nDofsVP_tot);
    
//     Jac.resize(nDofsVP * nDofsVP);

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].resize(nDofsV);    //resize
      std::fill(aResV[k].begin(), aResV[k].end(), 0);    //set aRes to zero
    
      aResVadj[k].resize(nDofsVadj);    //resize
      std::fill(aResVadj[k].begin(), aResVadj[k].end(), 0);    //set aRes to zero
      
      aResVctrl[k].resize(nDofsVctrl);    //resize
      std::fill(aResVctrl[k].begin(), aResVctrl[k].end(), 0);    //set aRes to zero

    }

    aResP.resize(nDofsP);    //resize
    std::fill(aResP.begin(), aResP.end(), 0);    //set aRes to zero

     aResPadj.resize(nDofsPadj);    //resize
    std::fill(aResPadj.begin(), aResPadj.end(), 0);    //set aRes to zero
    
    aResPctrl.resize(nDofsPctrl);    //resize
    std::fill(aResPctrl.begin(), aResPctrl.end(), 0);    //set aRes to zero

  //*************************************** 
  
  //***** set target domain flag ********************************** 
  geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = femus::ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(geom_element_iel.get_elem_center_3d()/*elem_center*/);
//***************************************   
    
    
   //STATE###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof // via local to global solution node

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        L2G_dofmap_Mat[i + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);
      }
    }
    
    // pressure *************
    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // global to global mapping between solution node and solution dof // via local to global solution node
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      L2G_dofmap_Mat[i + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
//STATE###################################################################

//ADJ###################################################################
     // velocity ************
    for (unsigned i = 0; i < nDofsVadj; i++) {
      unsigned solVadjDof = msh->GetSolutionDof(i, iel, solVadjType);    // global to global mapping between solution node and solution dof // via local to global solution node

      for (unsigned  k = 0; k < dim; k++) {
        solVadj[k][i] = (*sol->_Sol[solVadjIndex[k]])(solVadjDof);      // global extraction and local storage for the solution
        L2G_dofmap_Mat[i + k * nDofsV +nDofsVP] = pdeSys->GetSystemDof(solVadjIndex[k], solVPdeadjIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    
    // pressure *************
    for (unsigned i = 0; i < nDofsPadj; i++) {
      unsigned solPadjDof = msh->GetSolutionDof(i, iel, solPadjType);    // global to global mapping between solution node and solution dof // via local to global solution node
      solPadj[i] = (*sol->_Sol[solPadjIndex])(solPadjDof);      // global extraction and local storage for the solution
      L2G_dofmap_Mat[i + dim * nDofsV +nDofsVP] = pdeSys->GetSystemDof(solPadjIndex, solPPdeadjIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
//ADJ###################################################################


//CTRL###################################################################
     // velocity ************
    for (unsigned i = 0; i < nDofsVctrl; i++) {
      unsigned solVctrlDof = msh->GetSolutionDof(i, iel, solVctrlType);    // global to global mapping between solution node and solution dof // via local to global solution node

      for (unsigned  k = 0; k < dim; k++) {
        solVctrl[k][i] = (*sol->_Sol[solVctrlIndex[k]])(solVctrlDof);      // global extraction and local storage for the solution
        L2G_dofmap_Mat[i + k * nDofsV + 2*nDofsVP] = pdeSys->GetSystemDof(solVctrlIndex[k], solVPdectrlIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    
    // pressure *************
    for (unsigned i = 0; i < nDofsPctrl; i++) {
      unsigned solPctrlDof = msh->GetSolutionDof(i, iel, solPctrlType);    // global to global mapping between solution node and solution dof // via local to global solution node
      solPctrl[i] = (*sol->_Sol[solPctrlIndex])(solPctrlDof);      // global extraction and local storage for the solution
      L2G_dofmap_Mat[i + dim * nDofsV + 2*nDofsVP] = pdeSys->GetSystemDof(solPctrlIndex, solPPdectrlIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
//CTRL###################################################################




      // start a new recording of all the operations involving adept::adouble variables
      s.new_recording();

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {

//STATE#############################################################################3	
        // *** get gauss point weight, test function and test function partial derivatives ***
   elem_all[ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), ig, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);
   
    AbsDetJxWeight_iqp = detJac_iqp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[ig];
   
    elem_all[ielGeom][solVType]->shape_funcs_current_elem(ig, JacI_iqp, phiV_gss, phiV_x_gss, phiV_xx_gss , space_dim);

        phiP_gss = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

	
        std::vector < adept::adouble > solV_gss(dim, 0);
        std::vector < std::vector < adept::adouble > > gradSolV_gss(dim);
        std::vector < double > coordX_gss(dim, 0.);

        for (unsigned  k = 0; k < dim; k++) {
          gradSolV_gss[k].resize(dim_offset_grad);
          std::fill(gradSolV_gss[k].begin(), gradSolV_gss[k].end(), 0);
        }

          for (unsigned  k = 0; k < dim; k++) {
        for (unsigned i = 0; i < geom_element_iel.get_coords_at_dofs_3d()[k].size(); i++) {
            coordX_gss[k] += geom_element_iel.get_coords_at_dofs_3d()[k][i] * phiV_gss[i];    ///@todo change phiV into phi_coords!!!
             }
          }
          
          for (unsigned  k = 0; k < dim; k++) {
        for (unsigned i = 0; i < solV[k].size(); i++) {
            
            solV_gss[k]   += solV[k][i] * phiV_gss[i];
        

          for (unsigned j = 0; j < dim_offset_grad; j++) {
              gradSolV_gss[k][j] += solV[k][i] * phiV_x_gss[i * dim_offset_grad + j] ;
            }
            
          }
        }

        adept::adouble solP_gss = 0;

        for (unsigned i = 0; i < nDofsP; i++) {
          solP_gss += phiP_gss[i] * solP[i];
        }


//STATE###############################################################################

//ADJOINT#############################################################################3	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[ielGeom][solVadjType]->shape_funcs_current_elem(ig, JacI_iqp, phiVadj_gss, phiVadj_x_gss, phiVadj_xx_gss , space_dim);

        phiPadj_gss = msh->_finiteElement[ielGeom][solPadjType]->GetPhi(ig);

	
        std::vector < adept::adouble > solVadj_gss(dim, 0);
        std::vector < std::vector < adept::adouble > > gradSolVadj_gss(dim);

        for (unsigned  k = 0; k < dim; k++) {
          gradSolVadj_gss[k].resize(dim_offset_grad /*space_dim*/);
          std::fill(gradSolVadj_gss[k].begin(), gradSolVadj_gss[k].end(), 0);
        }

        for (unsigned i = 0; i < nDofsVadj; i++) {
          for (unsigned  k = 0; k < dim; k++) {
            solVadj_gss[k] += phiVadj_gss[i] * solVadj[k][i];
          }

          for (unsigned j = 0; j < dim_offset_grad /*space_dim*/; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              gradSolVadj_gss[k][j] += phiVadj_x_gss[i * dim_offset_grad /*space_dim*/ + j] * solVadj[k][i];
            }
          }
        }

        adept::adouble solPadj_gss = 0;

        for (unsigned i = 0; i < nDofsPadj; i++) {
          solPadj_gss += phiPadj_gss[i] * solPadj[i];
        }

//ADJOINT###############################################################################



//CONTROL#############################################################################3	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[ielGeom][solVctrlType]->shape_funcs_current_elem(ig, JacI_iqp, phiVctrl_gss, phiVctrl_x_gss, phiVctrl_xx_gss , space_dim);

        phiPctrl_gss = msh->_finiteElement[ielGeom][solPctrlType]->GetPhi(ig);

	
        std::vector < adept::adouble > solVctrl_gss(dim, 0);
        std::vector < std::vector < adept::adouble > > gradSolVctrl_gss(dim);

        for (unsigned  k = 0; k < dim; k++) {
          gradSolVctrl_gss[k].resize(dim_offset_grad /*space_dim*/);
          std::fill(gradSolVctrl_gss[k].begin(), gradSolVctrl_gss[k].end(), 0);
        }

        for (unsigned i = 0; i < nDofsVctrl; i++) {
          for (unsigned  k = 0; k < dim; k++) {
            solVctrl_gss[k] += phiVctrl_gss[i] * solVctrl[k][i];
          }

          for (unsigned j = 0; j < dim_offset_grad /*space_dim*/; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              gradSolVctrl_gss[k][j] += phiVctrl_x_gss[i * dim_offset_grad /*space_dim*/ + j] * solVctrl[k][i];
            }
          }
        }

        adept::adouble solPctrl_gss = 0;

        for (unsigned i = 0; i < nDofsPctrl; i++) {
          solPctrl_gss += phiPctrl_gss[i] * solPctrl[i];
        }
//CONTROL###############################################################################


#if exact_sol_flag == 1
//computation of RHS (force and desired velocity) using MMS=============================================== 
//state values--------------------
vector <double>  exact_stateVel(dim, 0.);
   mms_state_control::value_stateVel(coordX_gss, exact_stateVel);
vector < std::vector < double > > exact_grad_stateVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_stateVel[k].resize(dim);
    std::fill(exact_grad_stateVel[k].begin(), exact_grad_stateVel[k].end(), 0.);
}
   mms_state_control::gradient_stateVel(coordX_gss,exact_grad_stateVel);
vector <double>  exact_lap_stateVel(dim, 0.);
   mms_state_control::laplace_stateVel(coordX_gss, exact_lap_stateVel);
vector <double> exact_grad_statePress(dim, 0.);
   mms_state_control::gradient_statePress(coordX_gss, exact_grad_statePress);

//control values-------------------------------
vector <double>  exact_ctrlVel(dim);
   mms_state_control::value_ctrlVel(coordX_gss, exact_ctrlVel);
vector < std::vector < double > > exact_grad_ctrlVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_ctrlVel[k].resize(dim);
    std::fill(exact_grad_ctrlVel[k].begin(), exact_grad_ctrlVel[k].end(), 0.);
}
   mms_state_control::gradient_ctrlVel(coordX_gss,exact_grad_ctrlVel);
vector <double>  exact_lap_ctrlVel(dim);
   mms_state_control::laplace_ctrlVel(coordX_gss, exact_lap_ctrlVel);

//convection terms from delta_state-------------------------------------
vector <double>  exact_conv_u_nabla_u(dim,0.);
vector <double>  exact_conv_u_nabla_uctrl(dim,0.);
vector <double>  exact_conv_uctrl_nabla_u(dim,0.);
vector <double>  exact_conv_uctrl_nabla_uctrl(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_u[k] += exact_grad_stateVel[k][i] * exact_stateVel[i] ; 
    exact_conv_u_nabla_uctrl[k] += exact_grad_ctrlVel[k][i] * exact_stateVel[i] ; 
    exact_conv_uctrl_nabla_u[k] += exact_grad_stateVel[k][i] * exact_ctrlVel[i] ; 
    exact_conv_uctrl_nabla_uctrl[k] += exact_grad_ctrlVel[k][i] * exact_ctrlVel[i] ; 
    }
}

//convection terms from delta_adjoint-------------------------
vector <double>  exact_conv_u_nabla_uadj(dim,0.);
vector <double>  exact_conv_nabla_uT_uadj(dim,0.);
vector <double>  exact_conv_nabla_uctrlT_uadj(dim,0.);
vector <double>  exact_conv_uctrl_nabla_uadj(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_uadj[k] += exact_grad_stateVel[k][i] * exact_stateVel[i] ; 
    exact_conv_nabla_uT_uadj[k] += exact_grad_stateVel[i][k] * exact_stateVel[i];
    exact_conv_nabla_uctrlT_uadj[k] += exact_grad_ctrlVel[i][k] * exact_stateVel[i];  
    exact_conv_uctrl_nabla_uadj[k] += exact_grad_stateVel[k][i] * exact_ctrlVel[i] ; 
    }
}

//force and desired velocity ---------------------------------------------
vector <double> exactForce(dim,0.);
vector <double> exactVel_d(dim,0.);
for (unsigned k = 0; k < dim; k++){
    exactForce[k] = - IRe * exact_lap_stateVel[k] - IRe * exact_lap_ctrlVel[k] 
                    + advection_flag * (exact_conv_u_nabla_u[k] + exact_conv_u_nabla_uctrl[k] + exact_conv_uctrl_nabla_u[k] + exact_conv_uctrl_nabla_uctrl[k]) 
                    + exact_grad_statePress[k];
    exactVel_d[k] =   exact_stateVel[k] + exact_ctrlVel[k] 
                    + (1./cost_functional_coeff) * ( IRe * exact_lap_stateVel[k] - exact_grad_statePress[k]) 
                    + (1./cost_functional_coeff) * advection_flag * (exact_conv_u_nabla_uadj[k] - exact_conv_nabla_uT_uadj[k] - exact_conv_nabla_uctrlT_uadj[k] + exact_conv_uctrl_nabla_uadj[k]);
}

//computation of RHS (force and desired velocity) using MMS=============================================== 
#endif



        // *** phiV_i loop ***
        for (unsigned i = 0; i < nDofsV; i++) {
          std::vector < adept::adouble > NSV_gss(dim, 0.);
	  std::vector < adept::adouble > NSVadj_gss(dim, 0.);
	  std::vector < adept::adouble > NSVctrl_gss(dim, 0.);
	  
          for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      
          for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { //focus on single partial derivative
	    
          
	      NSV_gss[kdim]   	 	+=  IRe*phiV_x_gss[i * dim_offset_grad /*space_dim*/ + jdim]*gradSolV_gss[kdim][jdim]; 
						    // /*deformation_tensor*//*IRe * phiV_x_gss[i * dim + jdim] * 
						    // (gradSolV_gss[kdim][jdim] + gradSolV_gss[jdim][kdim])*/;  //diffusion
	      NSV_gss[kdim] 		+= IRe*phiV_x_gss[i * dim_offset_grad /*space_dim*/ + jdim]*gradSolVctrl_gss[kdim][jdim];	 //delta_state-control
	      NSVadj_gss[kdim]   	+=  IRe*phiVadj_x_gss[i * dim_offset_grad /*space_dim*/ + jdim]*gradSolVadj_gss[kdim][jdim];  
	      NSVctrl_gss[kdim]   	+=   beta * phiVctrl_x_gss[i * dim_offset_grad /*space_dim*/ + jdim] * gradSolVctrl_gss[kdim][jdim];
	      NSVctrl_gss[kdim] 	+=   - IRe*phiVctrl_x_gss[i * dim_offset_grad /*space_dim*/ + jdim]*gradSolVadj_gss[kdim][jdim];  //nabla_delta_control-nabla_adjoint
	  }  //jdim loop

          for (unsigned jdim = 0; jdim < dim; jdim++) { //focus on single partial derivative
 
          NSV_gss[kdim]  		+=  advection_flag * phiV_gss[i] * (solV_gss[jdim] * gradSolV_gss[kdim][jdim]);                                     //advection (u_hat . \nabla) u_hat
	      NSV_gss[kdim] 		+=  advection_flag * phiV_gss[i] * (solV_gss[jdim] * gradSolVctrl_gss[kdim][jdim]);                                  //advection (u_hat . \nabla) u_0
	      NSV_gss[kdim] 		+=  advection_flag * phiV_gss[i] * (solVctrl_gss[jdim] * gradSolV_gss[kdim][jdim]);                                  //advection (u_0 . \nabla) u_hat 
	      NSV_gss[kdim] 		+=  advection_flag * phiV_gss[i] * (solVctrl_gss[jdim] * gradSolVctrl_gss[kdim][jdim]);                              //advection (u_0 . \nabla) u_0
	      
	      NSVadj_gss[kdim]		+=   advection_flag * phiVadj_gss[i] * gradSolV_gss[jdim][kdim] * solVadj_gss[jdim];           // c(delta u,u,lambda)
	      NSVadj_gss[kdim]		+=   advection_flag * solV_gss[jdim] * phiVadj_x_gss[i * dim + jdim] * solVadj_gss[kdim] ;     // c(u,delta u, lambda)
	      NSVadj_gss[kdim]		+=   advection_flag * phiVadj_gss[i] * gradSolVctrl_gss[jdim][kdim] * solVadj_gss[jdim];       // c(delta u,u0,lambda)
	      NSVadj_gss[kdim]		+=   advection_flag * solVctrl_gss[jdim] * phiVadj_x_gss[i * dim + jdim] * solVadj_gss[kdim] ;       // c(u0,delta u,lambda)
	      
	      NSVctrl_gss[kdim]		+= -  advection_flag * solV_gss[jdim] * phiVctrl_x_gss[i * dim + jdim] * solVadj_gss[kdim] ;     // c(u,delta u0, lambda)
	      NSVctrl_gss[kdim]		+= -  advection_flag * phiVctrl_gss[i] * gradSolV_gss[jdim][kdim] * solVadj_gss[jdim];           // c(delta u0,u,lambda)
	      NSVctrl_gss[kdim]		+= -  advection_flag * phiVctrl_gss[i] * gradSolVctrl_gss[jdim][kdim] * solVadj_gss[jdim];       // c(delta u0,u0,lambda)
	      NSVctrl_gss[kdim]		+= -  advection_flag * solVctrl_gss[jdim] * phiVctrl_x_gss[i * dim + jdim] * solVadj_gss[kdim] ;       // c(u0,delta u0,lambda)
						  
	  }  //jdim loop
	  
 #if exact_sol_flag == 0
              NSV_gss[kdim]     += - force[kdim] * phiV_gss[i];
	      NSVadj_gss[kdim] 		+=  + cost_functional_coeff* target_flag * femus::ctrl::square_or_cube :: cost_functional_without_regularization::DesiredTargetVec()[kdim] * phiVadj_gss[i];
  	      NSVctrl_gss[kdim]   	+=  - cost_functional_coeff* target_flag * femus::ctrl::square_or_cube :: cost_functional_without_regularization::DesiredTargetVec()[kdim] * phiVctrl_gss[i];
#endif
 #if exact_sol_flag == 1
              NSV_gss[kdim]     += - exactForce[kdim] * phiV_gss[i];
	      NSVadj_gss[kdim] 		+=  + cost_functional_coeff* target_flag * exactVel_d[kdim] * phiVadj_gss[i];
 	      NSVctrl_gss[kdim]   	+=  - cost_functional_coeff* target_flag * exactVel_d[kdim] * phiVctrl_gss[i];
#endif        
              
              
          NSVadj_gss[kdim]		+=  - cost_functional_coeff * target_flag * solV_gss[kdim]*phiVadj_gss[i]; //delta_adjoint-state
	      NSVadj_gss[kdim] 		+=  - cost_functional_coeff * target_flag * solVctrl_gss[kdim]*phiVadj_gss[i]; //delta_adjoint-control
	      NSVctrl_gss[kdim] 	+=    cost_functional_coeff * target_flag             * solV_gss[kdim]*phiVctrl_gss[i]; //delta_control-state
	      NSVctrl_gss[kdim]   	+=   (cost_functional_coeff * target_flag + alpha) * solVctrl_gss[kdim] * phiVctrl_gss[i] ;
           
            //velocity-pressure block
          NSV_gss[kdim] 	    += - solP_gss * phiV_x_gss[i * dim_offset_grad /*space_dim*/ + kdim];
	      NSVadj_gss[kdim] 	    += - solPadj_gss * phiVadj_x_gss[i * dim_offset_grad /*space_dim*/ + kdim];
          NSVctrl_gss[kdim] 	+= - solPctrl_gss * phiVctrl_x_gss[i * dim_offset_grad /*space_dim*/ + kdim];
	    
	  } //kdim loop


          for (unsigned  kdim = 0; kdim < dim; kdim++) { // - (b-Ax)
            aResV[kdim][i] 		+=   NSV_gss[kdim]     * AbsDetJxWeight_iqp;
	    aResVadj[kdim][i]   	+=   NSVadj_gss[kdim]  * AbsDetJxWeight_iqp;
            aResVctrl[kdim][i]    	+=   NSVctrl_gss[kdim] * AbsDetJxWeight_iqp;
							 	    
	  }
        } // end phiV_i loop

        // *** phiP_i loop ***
        for (unsigned i = 0; i < nDofsP; i++) {
          for (int kdim = 0; kdim < dim; kdim++) {
            aResP[i] 		+= - (gradSolV_gss[kdim][kdim]) * phiP_gss[i]  * AbsDetJxWeight_iqp;
	    aResPadj[i]  	+= - (gradSolVadj_gss[kdim][kdim]) * phiPadj_gss[i]  * AbsDetJxWeight_iqp;
	    aResPctrl[i]   	+= - (gradSolVctrl_gss[kdim][kdim]) * phiPctrl_gss[i]  * AbsDetJxWeight_iqp;
	    
          }
        } // end phiP_i loop
      } // end gauss point loop
      
              



    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsVP_tot);    //resize

    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
        Res[ i +  kdim * nDofsV ]            =  - aResV[kdim][i].value();
	Res[ i +  kdim * nDofsV + nDofsVP]   =  - aResVadj[kdim][i].value();
	Res[ i +  kdim * nDofsV + 2*nDofsVP] =  - aResVctrl[kdim][i].value();
      }
    }

    for (int i = 0; i < nDofsP; i++) {
      Res[ i + dim * nDofsV ]            = - aResP[i].value();
      Res[ i + dim * nDofsV + nDofsVP]   = - aResPadj[i].value();
      Res[ i + dim * nDofsV + 2*nDofsVP] = - aResPctrl[i].value();
    }
  

    RES->add_vector_blocked(Res, L2G_dofmap_Mat);

    //Extarct and store the Jacobian
       Jac.resize(nDofsVP_tot * nDofsVP_tot);

      // define the dependent variables
      for (unsigned  kdim = 0; kdim < dim; kdim++) { s.dependent(&aResV[kdim][0], nDofsV);}              s.dependent(&aResP[0], nDofsP);
      for (unsigned  kdim = 0; kdim < dim; kdim++) { s.dependent(&aResVadj[kdim][0], nDofsVadj); }       s.dependent(&aResPadj[0], nDofsPadj);
      for (unsigned  kdim = 0; kdim < dim; kdim++) { s.dependent(&aResVctrl[kdim][0], nDofsVctrl);  }    s.dependent(&aResPctrl[0], nDofsPctrl);
      
      // define the independent variables
      for (unsigned  kdim = 0; kdim < dim; kdim++) {  s.independent(&solV[kdim][0], nDofsV); }         s.independent(&solP[0], nDofsP);
      for (unsigned  kdim = 0; kdim < dim; kdim++) { s.independent(&solVadj[kdim][0], nDofsVadj); }    s.independent(&solPadj[0], nDofsPadj);
      for (unsigned  kdim = 0; kdim < dim; kdim++) { s.independent(&solVctrl[kdim][0], nDofsVctrl); }  s.independent(&solPctrl[0], nDofsPctrl);

      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0] , true);
      
      JAC->add_matrix_blocked(Jac, L2G_dofmap_Mat, L2G_dofmap_Mat);
 
      s.clear_independents();
      s.clear_dependents();
    
  } //end element loop for each process

  RES->close();
//   RES->print();

  JAC->close();

  // ***************** END ASSEMBLY *******************
}



static void assemble_ns_dirichlet_control_lifting_internal_nonAD(MultiLevelProblem& ml_prob) {
     
 std::cout << " ********************************  NON-AD SYSTEM ******************************************** " << std::endl;

  // ======= Main objects - BEGIN =======
 //System
 NonLinearImplicitSystemWithPrimalDualActiveSetMethod * mlPdeSys  = & ml_prob.get_system< NonLinearImplicitSystemWithPrimalDualActiveSetMethod >("NSOpt");
  
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  //System, Equation
  const char* system_name            = mlPdeSys->name().c_str();
  LinearEquationSolver*  pdeSys	 = mlPdeSys->_LinSolver[level];   
  bool assembleMatrix = mlPdeSys->GetAssembleMatrix(); 
   
  SparseMatrix*	 JAC	= pdeSys->_KK;
  NumericVector* RES 	= pdeSys->_RES;
    

  constexpr bool print_algebra_global = true;
  constexpr bool print_algebra_local = true;

  //Mesh
  Mesh*		 msh    = ml_prob._ml_msh->GetLevel(level);
  elem*		 el	= msh->GetMeshElements();
  const unsigned dim 	= msh->GetDimension();
  unsigned dim2     = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned nel		= msh->GetNumberOfElements();
  unsigned igrid	= msh->GetLevel();
  unsigned iproc 	= msh->processor_id();
  unsigned   nprocs = msh->n_processors();
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));


  //Solution
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;
  Solution*	 sol  	         = ml_prob._ml_sol->GetSolutionLevel(level);
  
  // ======= Main objects - END =======
  
  
  
  // ======= Geometry at Dofs - BEGIN  =======
  unsigned coordXType = 2; /*CONTINUOUS_BIQUADRATIC*/// get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  unsigned solType_coords = coordXType;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  std::vector < std::vector < double> > coordX(dim);
  for(int i = 0; i < dim; i++) { coordX[i].reserve(max_size);  }

  constexpr unsigned int space_dim = 3;
 
  const unsigned int dim_offset_grad = /*dim*/  3  /*2*/    ;
  // ======= Geometry at Dofs - END  =======
  
  // ======= Solutions, Unknowns - BEGIN =======

  std::vector< unsigned >  vector_offsets = navier_stokes::lifting_internal:: provide_state_adj_ctrl_mu_offsets(dim);


const int state_pos_begin   =  vector_offsets[navier_stokes::lifting_internal:: pos_index_state];
  const int adj_pos_begin   =  vector_offsets[navier_stokes::lifting_internal:: pos_index_adj];
 const int ctrl_pos_begin   =  vector_offsets[navier_stokes::lifting_internal:: pos_index_ctrl];
   const int mu_pos_begin   =  vector_offsets[navier_stokes::lifting_internal:: pos_index_mu];
   

  const int press_type_pos = state_pos_begin + dim;

  
  std::vector< Unknown > unknowns = navier_stokes::lifting_internal:: provide_list_of_unknowns( dim );
  
  const int n_unknowns = unknowns.size();
 
  
  std::vector < unsigned > SolPdeIndex(n_unknowns);
  std::vector < unsigned > SolIndex(n_unknowns);  
  std::vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(unknowns[ivar]._name.c_str());
    SolIndex[ivar]	= ml_sol->GetIndex        (unknowns[ivar]._name.c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

  std::vector < unsigned > Sol_n_el_dofs_Mat_vol(n_unknowns);
  
  
      const unsigned int n_components_ctrl = dim;
  double penalty_outside_control_domain = _lifting_internal_penalty_outside_control_domain;         ///@todo  this number affects convergence or not! // penalty for zero control outside
  // ======= Solutions, Unknowns - END =======

      
  // ======= Solutions, not Unknowns - BEGIN =======

  //MU
  //************** variables for ineq constraints: act flag ****************************   
  std::vector<unsigned int> solIndex_act_flag_sol(n_components_ctrl); 
  
  femus::ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::store_act_flag_in_old(mlPdeSys, ml_sol, sol, solIndex_act_flag_sol);
  //************** variables for ineq constraints: act flag ****************************   
    

  const int ineq_flag = INEQ_FLAG;
  const double c_compl = C_COMPL;
  std::vector < std::vector < double/*int*/ > > sol_actflag(n_components_ctrl);   
  std::vector < std::vector < double > > ctrl_lower(n_components_ctrl);           
  std::vector < std::vector < double > > ctrl_upper(n_components_ctrl);           
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
      sol_actflag[c].reserve(max_size);
      ctrl_lower[c].reserve(max_size);
      ctrl_upper[c].reserve(max_size);
      }

//MU
  std::vector<unsigned int> ctrl_index_in_mat(n_components_ctrl); 
  std::vector<unsigned int>   mu_index_in_mat(n_components_ctrl);    
      for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
           ctrl_index_in_mat[kdim] =  SolPdeIndex[ctrl_pos_begin + kdim];
             mu_index_in_mat[kdim] =  SolPdeIndex[mu_pos_begin + kdim];
      }
      
  // ======= Solutions, not Unknowns - END =======
 

  
  
  // ======= Solutions, Unknowns at dofs - BEGIN =======
  std::vector < std::vector < double > > Sol_eldofs_Mat(n_unknowns);
  std::vector < std::vector < double > > gradSol_eldofs_Mat(n_unknowns);
  
  for(int k=0; k<n_unknowns; k++) {
    Sol_eldofs_Mat[k].reserve(max_size);
    gradSol_eldofs_Mat[k].reserve(max_size*dim);    
  }
  // ======= Solutions, Unknowns at dofs - END =======

  // ======= Solutions, Unknowns at quadrature points - BEGIN =======
  std::vector < double > SolVAR_qp(n_unknowns);
    std::vector < std::vector < double > > gradSolVAR_qp(n_unknowns);
    for(int k=0; k<n_unknowns; k++) {  gradSolVAR_qp[k].resize(dim_offset_grad /*space_dim*/);  }
  // ======= Solutions, Unknowns at quadrature points - END =======
      
  //============ Geometry at Quadrature points - BEGIN ==============================================================================
     std::vector < std::vector < double > >  JacI_iqp(space_dim);
     std::vector < std::vector < double > >  Jac_iqp(dim);
    for (unsigned d = 0; d < Jac_iqp.size(); d++) {   Jac_iqp[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp.size(); d++) { JacI_iqp[d].resize(dim); }
    
    double AbsDetJxWeight_iqp = 0.;
    double detJac_iqp;
  //============ Geometry at Quadrature points - END ==============================================================================

   
// ======= FE at Quadrature, all - BEGIN =======
    //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
  
  //==========================================================================================
  std::vector < std::vector < double > > phi_gss_fe(NFE_FAMS);
  std::vector < std::vector < double > > phi_x_gss_fe(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(max_size);
      phi_x_gss_fe[fe].reserve(max_size * dim_offset_grad /*space_dim*/);
   }
// ======= FE at Quadrature, all - END =======


  
  // ======= Equation, local - BEGIN =======
  std::vector < std::vector < int > > L2G_dofmap_Mat(n_unknowns); 
  std::vector < std::vector < double > > Res(n_unknowns);
  std::vector < std::vector < std::vector < double > > > Jac(n_unknowns);
 
  for(int i = 0; i < n_unknowns; i++) {     
    L2G_dofmap_Mat[i].reserve(max_size);
      Res[i].reserve(max_size);
  }
   
  if(assembleMatrix) {
    for(int i = 0; i < n_unknowns; i++) {
      Jac[i].resize(n_unknowns);    
      for(int j = 0; j < n_unknowns; j++) {
	Jac[i][j].reserve(max_size*max_size);	
      }
    }
  }
  // ======= Equation, local - END =======
  
  // ======= Equation, global - BEGIN =======
   RES->zero();
    if(assembleMatrix) JAC->zero();
  // ======= Equation, global - END =======
  
    
    
  // ======= Parameters - BEGIN ======= 
  double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  
  const double cost_functional_coeff = COST_FUNCTIONAL_COEFF;
  const double alpha = ALPHA_CTRL_VOL;
  const double beta  = BETA_CTRL_VOL;
  // ======= Parameters - END =======

  
   
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

  // geometry *****************************
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
      const short unsigned ielGeom = geom_element_iel.geom_type();
  // geometry end *****************************

  
  // equation *****************************
    unsigned nDofsV = msh->GetElementDofNumber(iel, SolFEType[state_pos_begin]);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType[state_pos_begin + press_type_pos]);    // number of solution element dofs
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin]);    // number of solution element dofs
    unsigned nDofsPadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin + press_type_pos]);    // number of solution element dofs

    unsigned nDofsVctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin]);    // number of solution element dofs
    unsigned nDofsPctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin + press_type_pos] );    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + nDofsP;
    unsigned nDofsVP_tot = 3*nDofsVP;
  // equation end *****************************

    
    
    
    
  //***** set target domain flag ********************************** 
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = ctrl::square_or_cube :: cost_functional_without_regularization::ElementTargetFlag(geom_element_iel.get_elem_center_3d()/*elem_center*/);
   //***************************************       
   
   //###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	Sol_n_el_dofs_Mat_vol[k]=ndofs_unk;
       Sol_eldofs_Mat[k].resize(ndofs_unk);
       L2G_dofmap_Mat[k].resize(ndofs_unk); 
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
       Sol_eldofs_Mat[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
       L2G_dofmap_Mat[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    
    //###################################################################
    for(int ivar=0; ivar<n_unknowns; ivar++) {
      
      Res[SolPdeIndex[ivar]].resize(Sol_n_el_dofs_Mat_vol[ivar]);
      memset(&Res[SolPdeIndex[ivar]][0],0.,Sol_n_el_dofs_Mat_vol[ivar]*sizeof(double));
    }
   
    for(int ivar=0; ivar<n_unknowns; ivar++) {
      for(int jvar=0; jvar<n_unknowns; jvar++) {
      if(assembleMatrix){  //MISMATCH
	Jac[ SolPdeIndex[ivar] ][ SolPdeIndex[jvar] ].resize(Sol_n_el_dofs_Mat_vol[ivar]*Sol_n_el_dofs_Mat_vol[jvar]);
	memset(&Jac[SolPdeIndex[ivar]][SolPdeIndex[jvar]][0],0.,Sol_n_el_dofs_Mat_vol[ivar]*Sol_n_el_dofs_Mat_vol[jvar]*sizeof(double));
      }
    }
  }
  
    //=============================================================================

    


 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = ctrl:: square_or_cube:: lifting_internal< ctrl::GAMMA_CONTROL_LIST_OF_FACES_WITH_EXTREMES >::ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());

  std::vector< std::vector< int > > control_node_flag(n_components_ctrl);
       
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
              control_node_flag[c].resize(Sol_n_el_dofs_Mat_vol[ctrl_pos_begin + c]);
              std::fill(control_node_flag[c].begin(), control_node_flag[c].end(), 0);   
         }

  if (control_el_flag == 1) {
	  for (unsigned c = 0; c < n_components_ctrl; c++) {
      std::fill(control_node_flag[c].begin(), control_node_flag[c].end(), 1);
         }
  }
  //*************************************************** 
   
 ///@todo if I want to restrict the control lifting, I might just set to zero the dof values here! Also, I have to remove the equations for the unneeded dofs with a PENALTY!
 /// Well, this can also be done in the Initialization function I think... 
 /// The problem there is that it is only DOF-BASED, it doesn't receive the ELEMENT INFORMATION. We should change that

//    if (control_el_flag == 0) {
// 	  for (unsigned c = 0; c < n_components_ctrl; c++) {
//       std::fill(Sol_eldofs_Mat[ctrl_pos_begin + c].begin(), Sol_eldofs_Mat[ctrl_pos_begin + c].end(), 0.);
//          }
//    }

// it seems that with this the convergence for the control variables is worse...
 
 
   
      for(unsigned iqp = 0; iqp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); iqp++) {
	
    elem_all[ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), iqp, Jac_iqp, JacI_iqp, detJac_iqp, space_dim);
    AbsDetJxWeight_iqp = detJac_iqp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[iqp];
   
     for(int fe=0; fe < NFE_FAMS; fe++) {
    elem_all[ielGeom][fe]->shape_funcs_current_elem(iqp, JacI_iqp,phi_gss_fe[fe],phi_x_gss_fe[fe], boost::none, space_dim);
      }


 //begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk < n_unknowns; unk++) {
	  SolVAR_qp[unk] = 0.;
	  for(unsigned ivar2=0; ivar2<dim_offset_grad /*space_dim*/; ivar2++){ 
	    gradSolVAR_qp[unk][ivar2] = 0.; 
	  }
	  
	  for(unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[unk]; i++) {
	    SolVAR_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i] * Sol_eldofs_Mat[unk][i];
	    for(unsigned ivar2=0; ivar2 < dim_offset_grad /*space_dim*/; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim_offset_grad /*space_dim*/+ivar2] * Sol_eldofs_Mat[unk][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************
	
      std::vector < double > coordX_gss(dim, 0.);
 	for(unsigned k = 0; k <  dim; k++) {
	  for(unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[k]; i++) {
         coordX_gss[k] += coordX[k][i] * phi_gss_fe[ SolFEType[k] ][i];
      }
    }
	
//  // I x = 5 test ********************************
// 	for(unsigned i_unk = 0; i_unk<n_unknowns; i_unk++) { 
// 	    for(unsigned i_dof=0; i_dof < Sol_n_el_dofs_Mat_vol[i_unk]; i_dof++) {
// 		/*if ( i_unk!=0 && i_unk!=1 && i_unk!=2 && i_unk!=3 && i_unk!=4 && i_unk!=6 && i_unk!=7 )*/  Res[SolPdeIndex[i_unk]][i_dof] +=  (               0.* phi_gss_fe[SolFEType[i_unk]][i_dof] 
// 		                                    - SolVAR_qp[i_unk]*phi_gss_fe[SolFEType[i_unk]][i_dof] )*AbsDetJxWeight_iqp;
// 		  for(unsigned j_unk = 0; j_unk<n_unknowns; j_unk++) {
// 		  	for(unsigned j_dof=0; j_dof < Sol_n_el_dofs_Mat_vol[j_unk]; j_dof++) {
// 			  
// 		              if (i_unk==j_unk /*&& i_unk!=0 && i_unk!=1 && i_unk!=2 && i_unk!=3 && i_unk!=4 && i_unk!=6 && i_unk!=7*/)   {
// 				Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ][ i_dof*Sol_n_el_dofs_Mat_vol[i_unk] + j_dof ] += 
// 				        ( phi_gss_fe[SolFEType[i_unk]][i_dof]*phi_gss_fe[SolFEType[j_unk]][j_dof] )*AbsDetJxWeight_iqp;
// 			      }
// 			  
// 			} //j_dof
// 		  }  //j_unk
// 	    }  //i_dof
// 	}  //i_unk
//  // I x = 5 test ********************************
 

#if exact_sol_flag == 1
//computation of RHS (force and desired velocity) using MMS - BEGIN =============================================== 
//state values--------------------
vector <double>  exact_stateVel(dim, 0.);
   mms_state_control::value_stateVel(coordX_gss, exact_stateVel);
vector < std::vector < double > > exact_grad_stateVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_stateVel[k].resize(dim);
    std::fill(exact_grad_stateVel[k].begin(), exact_grad_stateVel[k].end(), 0.);
}
   mms_state_control::gradient_stateVel(coordX_gss,exact_grad_stateVel);
vector <double>  exact_lap_stateVel(dim, 0.);
   mms_state_control::laplace_stateVel(coordX_gss, exact_lap_stateVel);
vector <double> exact_grad_statePress(dim, 0.);
   mms_state_control::gradient_statePress(coordX_gss, exact_grad_statePress);

//control values-------------------------------
vector <double>  exact_ctrlVel(dim);
   mms_state_control::value_ctrlVel(coordX_gss, exact_ctrlVel);
vector < std::vector < double > > exact_grad_ctrlVel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_ctrlVel[k].resize(dim);
    std::fill(exact_grad_ctrlVel[k].begin(), exact_grad_ctrlVel[k].end(), 0.);
}
   mms_state_control::gradient_ctrlVel(coordX_gss,exact_grad_ctrlVel);
vector <double>  exact_lap_ctrlVel(dim);
   mms_state_control::laplace_ctrlVel(coordX_gss, exact_lap_ctrlVel);

//convection terms from delta_state-------------------------------------
vector <double>  exact_conv_u_nabla_u(dim,0.);
vector <double>  exact_conv_u_nabla_uctrl(dim,0.);
vector <double>  exact_conv_uctrl_nabla_u(dim,0.);
vector <double>  exact_conv_uctrl_nabla_uctrl(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_u[k] += exact_grad_stateVel[k][i] * exact_stateVel[i] ; 
    exact_conv_u_nabla_uctrl[k] += exact_grad_ctrlVel[k][i] * exact_stateVel[i] ; 
    exact_conv_uctrl_nabla_u[k] += exact_grad_stateVel[k][i] * exact_ctrlVel[i] ; 
    exact_conv_uctrl_nabla_uctrl[k] += exact_grad_ctrlVel[k][i] * exact_ctrlVel[i] ; 
    }
}

//convection terms from delta_adjoint-------------------------
vector <double>  exact_conv_u_nabla_uadj(dim,0.);
vector <double>  exact_conv_nabla_uT_uadj(dim,0.);
vector <double>  exact_conv_nabla_uctrlT_uadj(dim,0.);
vector <double>  exact_conv_uctrl_nabla_uadj(dim,0.);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_u_nabla_uadj[k] += exact_grad_stateVel[k][i] * exact_stateVel[i] ; 
    exact_conv_nabla_uT_uadj[k] += exact_grad_stateVel[i][k] * exact_stateVel[i];
    exact_conv_nabla_uctrlT_uadj[k] += exact_grad_ctrlVel[i][k] * exact_stateVel[i];  
    exact_conv_uctrl_nabla_uadj[k] += exact_grad_stateVel[k][i] * exact_ctrlVel[i] ; 
    }
}

//force and desired velocity ---------------------------------------------
vector <double> exactForce(dim,0.);
vector <double> exactVel_d(dim,0.);
for (unsigned k = 0; k < dim; k++){
    exactForce[k] = - IRe * exact_lap_stateVel[k] - IRe * exact_lap_ctrlVel[k] 
                    + advection_flag * (exact_conv_u_nabla_u[k] + exact_conv_u_nabla_uctrl[k] + exact_conv_uctrl_nabla_u[k] + exact_conv_uctrl_nabla_uctrl[k]) 
                    + exact_grad_statePress[k];
    exactVel_d[k] =   exact_stateVel[k] + exact_ctrlVel[k] 
                    + (1./cost_functional_coeff) * ( IRe * exact_lap_stateVel[k] - exact_grad_statePress[k]) 
                    + (1./cost_functional_coeff) * advection_flag * (exact_conv_u_nabla_uadj[k] - exact_conv_nabla_uT_uadj[k] - exact_conv_nabla_uctrlT_uadj[k] + exact_conv_uctrl_nabla_uadj[k]);
}

//computation of RHS (force and desired velocity) using MMS - END =============================================== 
#endif

 
 
 
//============ delta_state row - BEGIN  ============================================================================================

//************ Residual, Velocity, BEGIN *********************

    for (unsigned  kdim = 0; kdim < dim; kdim++) { // velocity block row 
          
         for (unsigned i = 0; i < nDofsV; i++) {

	              double lap_res_du_u_kdim_i 			= 0.; 
		      double lap_res_du_ctrl_kdim_i 			= 0.;
		      double adv_res_uold_nablauold_kdim_i 		= 0.;
		      double adv_res_uold_nablauctrlold_kdim_i 	= 0.;
		      double adv_res_uctrlold_nablauold_kdim_i 	= 0.;
		      double adv_res_uctrlold_nablauctrlold_kdim_i 	= 0.;
              
	      for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		    lap_res_du_u_kdim_i  		     += gradSolVAR_qp[SolPdeIndex[kdim]][jdim]		            * phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad + jdim];
		    lap_res_du_ctrl_kdim_i 		 += gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim] * phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad + jdim];
          }
	      for (unsigned jdim = 0; jdim < dim; jdim++) {
		   adv_res_uold_nablauold_kdim_i 	     += SolVAR_qp[SolPdeIndex[jdim]] 		  * gradSolVAR_qp[SolPdeIndex[kdim]][jdim]		    * phi_gss_fe[ SolFEType[kdim] ][i];
		  adv_res_uold_nablauctrlold_kdim_i 	 += SolVAR_qp[SolPdeIndex[jdim]] 		  * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
		  adv_res_uctrlold_nablauold_kdim_i	     += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]] * gradSolVAR_qp[SolPdeIndex[kdim]][jdim] 		    * phi_gss_fe[ SolFEType[kdim] ][i];
		  adv_res_uctrlold_nablauctrlold_kdim_i  += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
	      }
	      
	      
	      Res[kdim][i]   +=   AbsDetJxWeight_iqp * (         
#if exact_sol_flag == 0
                                         + force[kdim] * phi_gss_fe[ SolFEType[kdim] ][i]
 #endif                                      
 #if exact_sol_flag == 1
                                       + exactForce[kdim] * phi_gss_fe[ SolFEType[kdim] ][i]
 #endif
                                           - IRe*lap_res_du_u_kdim_i 
                                           - IRe*lap_res_du_ctrl_kdim_i
                                           - advection_flag * adv_res_uold_nablauold_kdim_i 
                                           - advection_flag * adv_res_uold_nablauctrlold_kdim_i
                                           - advection_flag * adv_res_uctrlold_nablauold_kdim_i
					                       - advection_flag * adv_res_uctrlold_nablauctrlold_kdim_i
					 );
	}
      
}
	
//************ Residual, Velocity, END *********************

//************ Jacobian, Velocity, BEGIN *********************

//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
              
  for (unsigned  kdim = 0; kdim < dim; kdim++) {
              
    for (unsigned i = 0; i < nDofsV; i++) {

	  for (unsigned j = 0; j < nDofsV; j++) {

    double lap_jac_du_u_kdim_i_j = 0.;
    double adv_uold_nablaunew_kdim_i_j = 0.;
    double adv_uctrlold_nablaunew_kdim_i_j = 0.;

        
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		    lap_jac_du_u_kdim_i_j += phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + jdim]	* 	phi_x_gss_fe[SolFEType[kdim]][ j * dim_offset_grad /*space_dim*/ + jdim ];
            }
          for (unsigned  jdim = 0; jdim < dim; jdim++) {
		    adv_uold_nablaunew_kdim_i_j     += SolVAR_qp[SolPdeIndex[jdim]]  * phi_x_gss_fe[ SolFEType[kdim] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
		    adv_uctrlold_nablaunew_kdim_i_j += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]] * phi_x_gss_fe[ SolFEType[kdim] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
            }  //jdim
                
              
        Jac[kdim][kdim][i * nDofsV + j] += (  IRe * lap_jac_du_u_kdim_i_j 
						    + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] 		  * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_hat_new,u_hat_old,delta_lambda) diagonal blocks  ..... unew_nablauold
						    + advection_flag * adv_uold_nablaunew_kdim_i_j 								 // c(u_hat_old, u_hat_new, delta_lambda)
						    + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][kdim] * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_hat_new,u_0_old,delta_lambda) diagonal blocks  ..... unew_nablauctrlold
						    + advection_flag * adv_uctrlold_nablaunew_kdim_i_j 						  	 // c(u_0_old, u_hat_new, delta_lambda)
						    ) * AbsDetJxWeight_iqp; 
						    
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		Jac[kdim][off_kdim][i * nDofsV + j] += ( +	advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][off_kdim] 		    * phi_gss_fe[ SolFEType[kdim] ][i]	// c(u_hat_new,u_hat_old,delta_lambda) off-diagonal blocks  ..... unew_nablauold
						      + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim] ][j] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][off_kdim] * phi_gss_fe[ SolFEType[kdim] ][i]	// c(u_hat_new,u_0_old,delta_lambda) off-diagonal blocks  ..... unew_nablauctrlold
						      ) * AbsDetJxWeight_iqp;
                              
	      
	} //j_du_u loop
	
   }//i_state

       
}
	      

//BLOCK delta_state - control------------------------------------------------------------------------------------
  for (unsigned  kdim = 0; kdim < dim; kdim++) {
              
   for (unsigned i = 0; i < nDofsV; i++) {
       
	for (unsigned j = 0; j < nDofsVctrl; j++) {
        
		       double  lap_jac_du_ctrl_kdim_i_j = 0.;
		       double  adv_uold_nablauctrlnew_kdim_i_j = 0.;
		       double  adv_uctrlold_nablauctrlnew_kdim_i_j = 0.;
              
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		    lap_jac_du_ctrl_kdim_i_j += phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + jdim] * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
            }
            
            
		  for (unsigned  jdim = 0; jdim < dim; jdim++) {  //diagonal blocks only
		    adv_uold_nablauctrlnew_kdim_i_j	 += SolVAR_qp[SolPdeIndex[jdim]] 		  * phi_x_gss_fe[ SolFEType[kdim +  ctrl_pos_begin] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
		    adv_uctrlold_nablauctrlnew_kdim_i_j += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]] * phi_x_gss_fe[ SolFEType[kdim +  ctrl_pos_begin] ][j * dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
                }  //jdim
	      


	      Jac[kdim][kdim + ctrl_pos_begin ][i*nDofsVctrl + j] += (+ IRe * lap_jac_du_ctrl_kdim_i_j 
									+ advection_flag * adv_uold_nablauctrlnew_kdim_i_j														 		 // c(u_hat_old, u_0_new, delta_lambda)
									+ advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim]		       * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_0_new,u_hat_old,delta_lambda) diagonal blocks  ..... uctrlnew_nablauold
									+ advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][kdim] * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_0_new,u_0_old,delta_lambda) diagonal blocks  ..... uctrlnew_nablauctrlold
									+ advection_flag * adv_uctrlold_nablauctrlnew_kdim_i_j															 // c(u_0_old, u_0_new, delta_lambda)
									) * AbsDetJxWeight_iqp;
									
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
		Jac[kdim][off_kdim + ctrl_pos_begin][i*nDofsVctrl + j] += (   +	advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j] * gradSolVAR_qp[SolPdeIndex[kdim]][off_kdim]		     * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_0_new,u_hat_old,delta_lambda) off-diagonal blocks  ..... uctrlnew_nablauold
									      + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][off_kdim] * phi_gss_fe[ SolFEType[kdim] ][i]	 // c(u_0_new,u_0_old,delta_lambda) off-diagonal blocks  ..... uctrlnew_nablauctrlold
									  ) * AbsDetJxWeight_iqp;
	} //j_du_ctrl loop

   }//i_state loop

       
}
   
   
      
//************ Jacobian, Velocity, END *********************


//************ Residual & Jacobian, Pressure, BEGIN *********************
   
//BLOCK Pressure
   for (unsigned  kdim = 0; kdim < dim; kdim++) {
         
     for (unsigned i = 0; i < nDofsV; i++) {
       	      Res[kdim][i]   +=   AbsDetJxWeight_iqp *  SolVAR_qp[SolPdeIndex[press_type_pos]] * phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + kdim];
              
      for (unsigned j = 0; j < nDofsP; j++) {
         Jac[kdim][press_type_pos][i * nDofsP + j] += -( phi_gss_fe[SolFEType[press_type_pos]][j] * phi_x_gss_fe[SolFEType[kdim]][i * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
         } //j_press loop
      
    }//i_state loop
      
  }

   
//DIV_state -----------
	    double div_u_du_qp = 0.;
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      div_u_du_qp += gradSolVAR_qp[SolPdeIndex[kdim]][kdim];
      }


  for (unsigned i = 0; i < nDofsP; i++) {  
      Res[press_type_pos][i]  +=  ( (div_u_du_qp) * phi_gss_fe[SolFEType[press_type_pos]][i] ) * AbsDetJxWeight_iqp;
   }
      
   for (unsigned  kdim = 0; kdim < dim; kdim++) {
     for (unsigned i = 0; i < nDofsP; i++) {
      for (unsigned j = 0; j < nDofsV; j++) {
	      Jac[press_type_pos][kdim][i * nDofsV + j] += - ( phi_gss_fe[SolFEType[press_type_pos]][i] * phi_x_gss_fe[SolFEType[kdim]][j * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
        } //j loop
	  }
	  
   }
//************ Residual & Jacobian, Pressure, END *********************


//============ delta_state row - END  ============================================================================================


    
//============ delta_adjoint row - BEGIN  =============================================================================================
  
  for (unsigned kdim = 0; kdim < dim; kdim++) { 
          
  for (unsigned i = 0; i < nDofsVadj; i++) {

		    double lap_res_dadj_adj_kdim_i 			= 0.;
		    double adv_res_phiadj_nablauold_uadjold_kdim_i 	= 0.;
		    double adv_res_uold_nablaphiadj_uadjold_kdim_i 	= 0.;
		    double adv_res_phiadj_nablauctrlold_uadjold_kdim_i = 0.;
		    double adv_res_uctrlold_nablaphiadj_uadjold_kdim_i = 0.;
            
        for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		   lap_res_dadj_adj_kdim_i    += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]][jdim]  *  phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim];
          }
          
       for (unsigned jdim = 0; jdim < dim; jdim++) {
		adv_res_phiadj_nablauold_uadjold_kdim_i     += phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim]][kdim] 			* SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
		adv_res_uold_nablaphiadj_uadjold_kdim_i     += SolVAR_qp[SolPdeIndex[jdim]]		       * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
		adv_res_phiadj_nablauctrlold_uadjold_kdim_i += phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim  + ctrl_pos_begin]][kdim]	* SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
		adv_res_uctrlold_nablaphiadj_uadjold_kdim_i += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]]  * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
	   }
	   
	   
	  Res[kdim + adj_pos_begin][i] += ( 
#if exact_sol_flag == 0
                            - cost_functional_coeff * target_flag * ctrl::square_or_cube :: cost_functional_without_regularization::DesiredTargetVec()[kdim] 			      * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
 #endif                                      
 #if exact_sol_flag == 1
                            - cost_functional_coeff * target_flag * exactVel_d[kdim] 			      * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
 #endif
					    + cost_functional_coeff * target_flag * SolVAR_qp[SolPdeIndex[kdim]] 		      * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
					    + cost_functional_coeff * target_flag * SolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
					    - IRe * lap_res_dadj_adj_kdim_i
					    - advection_flag * adv_res_phiadj_nablauold_uadjold_kdim_i
					    - advection_flag * adv_res_uold_nablaphiadj_uadjold_kdim_i
					    - advection_flag * adv_res_phiadj_nablauctrlold_uadjold_kdim_i
					    - advection_flag * adv_res_uctrlold_nablaphiadj_uadjold_kdim_i
					    ) * AbsDetJxWeight_iqp;
      }
  }
  
//BLOCK delta_adjoint - state------------------------------------------------------------------------------------------
  for (unsigned kdim = 0; kdim < dim; kdim++) {
          
   for (unsigned i = 0; i < nDofsVadj; i++) {
 
     for (unsigned j = 0; j < nDofsV; j++) {
         
	      Jac[kdim + adj_pos_begin][kdim][i*nDofsV + j] += ( - cost_functional_coeff * target_flag * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i] * phi_gss_fe[SolFEType[kdim]][j] 
								 + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i]    * phi_x_gss_fe[ SolFEType[kdim] ][j*dim_offset_grad /*space_dim*/ + kdim] 		* SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u, u_hat_new, lambda_old)  diagonal blocks  ......phiadj_nablaunew_uadjold 
								 + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim] ][j] 			* phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_hat_new, delta_u, lambda_old) diagonal blocks  ......unew_nablaphiadj_uadjold
								     ) * AbsDetJxWeight_iqp;
                                     
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
               
		Jac[kdim + adj_pos_begin][off_kdim][i * nDofsV + j] += (+ advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * phi_x_gss_fe[ SolFEType[off_kdim] ][j*dim_offset_grad /*space_dim*/ + kdim]		      * SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u, u_hat_new, lambda_old)  off-diagonal blocks  ......phiadj_nablaunew_uadjold 
								      + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim] ][j] 		  * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	   //c(u_hat_new, delta_u, lambda_old) off-diagonal blocks  ......unew_nablaphiadj_uadjold
								      ) * AbsDetJxWeight_iqp;
	    
	  }//j_dadj_u loop
     }
   }
   
   
//BLOCK delta_adjoint - control-----------------------------------------------------------------------------------------
  for (unsigned kdim = 0; kdim < dim; kdim++) {
      
   for (unsigned i = 0; i < nDofsVadj; i++) {
       
     for (unsigned j = 0; j < nDofsVctrl; j++) {
         
	     Jac[kdim + adj_pos_begin][kdim + ctrl_pos_begin][i*nDofsVctrl + j] += ( - cost_functional_coeff * target_flag * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j] 
										    + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i]  * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u, u_0_new, lambda_old)  diagonal blocks  ......phiadj_nablauctrlnew_uadjold 
										    + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j] * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(u_0_new, delta_u, lambda_old) diagonal blocks  ......uctrlnew_nablaphiadj_uadjold
										    ) * AbsDetJxWeight_iqp;
                                            
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
	     Jac[kdim + adj_pos_begin][off_kdim + ctrl_pos_begin][i*nDofsVctrl + j] += (+ advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i]      * phi_x_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u, u_0_new, lambda_old)  off-diagonal blocks  ......phiadj_nablauctrlnew_uadjold 
											+ advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j] * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	   //c(u_0_new, delta_u, lambda_old) off-diagonal blocks  ......uctrlnew_nablaphiadj_uadjold
											) * AbsDetJxWeight_iqp;
	    
	  }//j_dadj_ctrl loop
	  
     }
       
   }

//DIAG BLOCK delta_adjoint - adjoint---------------------------------------------------------------------------------
  for (unsigned  kdim = 0; kdim < dim; kdim++) {
              
   for (unsigned i = 0; i < nDofsVadj; i++) {
       
     for (unsigned j = 0; j < nDofsVadj; j++) {
         
		     double lap_jac_dadj_adj_kdim_i_j = 0.;
		     double adv_uold_nablaphiadj_uadjnew_kdim_i_j = 0.;
		     double adv_uctrlold_nablaphiadj_uadjnew_kdim_i_j = 0.;
            
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		  lap_jac_dadj_adj_kdim_i_j  +=  phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
            }

            for (unsigned jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
	     adv_uold_nablaphiadj_uadjnew_kdim_i_j    += SolVAR_qp[SolPdeIndex[jdim]] 		     * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j] ;
	     adv_uctrlold_nablaphiadj_uadjnew_kdim_i_j  += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]] * phi_x_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j] ;
	   }
	   

	   Jac[kdim + adj_pos_begin][kdim + adj_pos_begin][i * nDofsVadj + j] += ( + IRe*lap_jac_dadj_adj_kdim_i_j 
										    + advection_flag * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] 		  * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j]   //c(delta_u, u_hat_old, lambda_new)  diagonal blocks  ......phiadj_nablauold_uadjnew  
										    + advection_flag * adv_uold_nablaphiadj_uadjnew_kdim_i_j 	//c(u_hat_old, delta_uhat, lambda_new)
										    + advection_flag * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][kdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j]   //c(delta_u, u_0_old, lambda_new)  diagonal blocks  ......phiadj_nablauctrlold_uadjnew  
										    + advection_flag * adv_uctrlold_nablaphiadj_uadjnew_kdim_i_j   //c(u_0_old, delta_uhat, lambda_new)
											) * AbsDetJxWeight_iqp;
                                            
               unsigned int off_kdim = (kdim+1) % dim; //off-diagonal blocks
		  Jac[kdim + adj_pos_begin][off_kdim + adj_pos_begin][i * nDofsVadj + j] += (+ advection_flag * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim]][kdim]		    * phi_gss_fe[ SolFEType[off_kdim + adj_pos_begin] ][j]   //c(delta_u, u_hat_old, lambda_new)  off-diagonal blocks  ......phiadj_nablauold_uadjnew   
											   + advection_flag * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim + ctrl_pos_begin]][kdim] * phi_gss_fe[ SolFEType[off_kdim + adj_pos_begin] ][j]   //c(delta_u, u_0_old, lambda_new)  off-diagonal blocks  ......phiadj_nablauctrlold_uadjnew  
											  ) * AbsDetJxWeight_iqp;
	    
      } //j_dadj_adj loop
      
  }//i_adj loop
  
 }

  
//************ Residual & Jacobian, Pressure Adj, BEGIN *********************
//BLOCK Pressure_adj
  
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
  for (unsigned i = 0; i < nDofsVadj; i++) {
      
      Res[kdim + adj_pos_begin][i] +=  AbsDetJxWeight_iqp * (SolVAR_qp[SolPdeIndex[press_type_pos + adj_pos_begin]] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim]);
      
    for (unsigned j = 0; j < nDofsPadj; j++) {
	      Jac[kdim + adj_pos_begin][press_type_pos + adj_pos_begin][i*nDofsPadj + j] += - (    phi_gss_fe[SolFEType[press_type_pos + adj_pos_begin]][j] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
           }//j_press_adj loop
        }//i_adj loop
	  }

//DIV_adj
		double div_adj_dadj_qp = 0.;
      for (unsigned kdim = 0; kdim < dim; kdim++) {
	    div_adj_dadj_qp += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin ]][kdim] ;
      }
      
  for (unsigned i = 0; i < nDofsPadj; i++) {
      Res[press_type_pos + adj_pos_begin][i] += ( (div_adj_dadj_qp) * phi_gss_fe[ SolFEType[press_type_pos + adj_pos_begin] ][i] ) * AbsDetJxWeight_iqp;
  }//i_div_adj
      
    for (unsigned kdim = 0; kdim < dim; kdim++) {
     for (unsigned i = 0; i < nDofsPadj; i++) {
      for (unsigned j = 0; j < nDofsVadj; j++) {
	    Jac[press_type_pos + adj_pos_begin][kdim + adj_pos_begin][i * nDofsVadj + j] += - (   phi_gss_fe[SolFEType[press_type_pos + adj_pos_begin]][i] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
          }//j loop
        }//i_div_adj
	   }
//************ Residual & Jacobian, Pressure Adj, END *********************

//============ delta_adjoint row - END  =============================================================================================


//============ delta_control row - BEGIN ==================================================================================================


//************ Residual, BEGIN *********************

      if ( control_el_flag == 1)    {
          
      for (unsigned kdim = 0; kdim < dim; kdim++) {
          
         for (unsigned i = 0; i < nDofsVctrl; i++) {
      
		    double lap_res_dctrl_ctrl_kdim_i 			 = 0.;
		    double lap_res_dctrl_adj_kdim_i 			 = 0.;
		    double adv_res_phictrl_nablauold_uadjold_kdim_i 	 = 0.;
		    double adv_res_uold_nablaphictrl_uadjold_kdim_i 	 = 0.;
		    double adv_res_phictrl_nablauctrlold_uadjold_kdim_i = 0.;
		    double adv_res_uctrlold_nablaphictrl_uadjold_kdim_i = 0.;
            
     for (unsigned jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) {
		lap_res_dctrl_ctrl_kdim_i		      += gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim] 	*   phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim];
		lap_res_dctrl_adj_kdim_i 		      += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]][jdim] 	*   phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim];
         }
         
      for (unsigned jdim = 0; jdim < dim; jdim++) {
		adv_res_uold_nablaphictrl_uadjold_kdim_i     += SolVAR_qp[SolPdeIndex[jdim]]			 * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
		adv_res_phictrl_nablauold_uadjold_kdim_i     += phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim]][kdim]   			   * SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
		adv_res_phictrl_nablauctrlold_uadjold_kdim_i += phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[jdim  + ctrl_pos_begin]][kdim]	   * SolVAR_qp[SolPdeIndex[jdim + adj_pos_begin]];
		adv_res_uctrlold_nablaphictrl_uadjold_kdim_i += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]]   * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]  * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
      }
      
      
      Res[kdim + ctrl_pos_begin][i] +=  AbsDetJxWeight_iqp * (
#if exact_sol_flag == 0
                     + cost_functional_coeff * target_flag * ctrl::square_or_cube :: cost_functional_without_regularization::DesiredTargetVec()[kdim] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]
 #endif                                      
 #if exact_sol_flag == 1
                     + cost_functional_coeff * target_flag * exactVel_d[kdim] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]
 #endif
					 - cost_functional_coeff * target_flag * SolVAR_qp[SolPdeIndex[kdim]] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]
					 - cost_functional_coeff * target_flag * SolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]
					 - alpha * SolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]
					 - beta * lap_res_dctrl_ctrl_kdim_i
					 + IRe * lap_res_dctrl_adj_kdim_i
					+ advection_flag * adv_res_uold_nablaphictrl_uadjold_kdim_i
					+ advection_flag * adv_res_phictrl_nablauold_uadjold_kdim_i
					+ advection_flag * adv_res_phictrl_nablauctrlold_uadjold_kdim_i
					+ advection_flag * adv_res_uctrlold_nablaphictrl_uadjold_kdim_i				
					);
//       }
      
      
      }

          
    }
    
      }
      
    else if ( control_el_flag == 0) {
    
        for (unsigned kdim = 0; kdim < dim; kdim++) {
           for (unsigned i = 0; i < nDofsVctrl; i++) {

        Res[kdim + ctrl_pos_begin][i] +=       (- penalty_outside_control_domain)  *  (1 - control_node_flag[kdim][i]) * (Sol_eldofs_Mat[kdim + ctrl_pos_begin][i] - 0.);
                 }
            }
    }
  
//************ Residual, END *********************

//************ Jacobian, BEGIN *********************
      if ( control_el_flag == 1)    {
//BLOCK delta_control - state------------------------------------------------------------------------------------------------
   for (unsigned kdim = 0; kdim < dim; kdim++) {

     for (unsigned i = 0; i < nDofsVctrl; i++) {
    
      for (unsigned j = 0; j < nDofsV; j++) {
          
	      Jac[kdim + ctrl_pos_begin][kdim][i*nDofsV + j] +=  AbsDetJxWeight_iqp * (
                                 + cost_functional_coeff * target_flag * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i] * phi_gss_fe[SolFEType[kdim]][j] 
								 - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim] ][j] 		       * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_hat_new, delta_u0, lambda_old) diagonal blocks  ......unew_nablaphictrl_uadjold
								 - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i]  * phi_x_gss_fe[ SolFEType[kdim] ][j*dim_offset_grad /*space_dim*/ + kdim] 			* SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u0, u_hat_new, lambda_old)  diagonal blocks  ......phictrl_nablaunew_uadjold 
								    );
          
               unsigned int off_kdim = (kdim+1) % dim; //off-diagonal blocks
               
	      Jac[kdim + ctrl_pos_begin][off_kdim][i * nDofsV + j] +=  AbsDetJxWeight_iqp * ( - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim] ][j] 		   * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_hat_new, delta_u0, lambda_old) off-diagonal blocks  ......unew_nablaphictrl_uadjold
								      - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i] * phi_x_gss_fe[ SolFEType[off_kdim] ][j*dim_offset_grad /*space_dim*/ + kdim] 			* SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u0, u_hat_new, lambda_old)  off-diagonal blocks  ......phictrl_nablaunew_uadjold 
								    );
	    
        }//j_dctrl_u loop
     }
  }
      
//BLOCK delta_control - adjoint------------------------------------------------------------------------------------------------
 for (unsigned  kdim = 0; kdim < dim; kdim++) {
              
   for (unsigned i = 0; i < nDofsVctrl; i++) {
      for (unsigned j = 0; j < nDofsVadj; j++) {
          
		    double  lap_jac_dctrl_adj_kdim_i_j = 0.;
		    double  adv_uold_nablaphictrl_uadjnew_kdim_i_j = 0.;
		    double  adv_uctrlold_nablaphictrl_uadjnew_kdim_i_j = 0.;
            
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		lap_jac_dctrl_adj_kdim_i_j += phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]	*	phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
            }

            for (unsigned jdim = 0; jdim < dim; jdim++) { //diagonal blocks only
	     adv_uold_nablaphictrl_uadjnew_kdim_i_j     += SolVAR_qp[SolPdeIndex[jdim]] 		     * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j] ;
	     adv_uctrlold_nablaphictrl_uadjnew_kdim_i_j += SolVAR_qp[SolPdeIndex[jdim + ctrl_pos_begin]]* phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + jdim] * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j] ;
	   }
	   

	   Jac[kdim + ctrl_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += AbsDetJxWeight_iqp * ( - IRe * lap_jac_dctrl_adj_kdim_i_j 
										    - advection_flag * adv_uold_nablaphictrl_uadjnew_kdim_i_j 	//c(u_hat_old, delta_u0, lambda_new)
										    - advection_flag * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] 		    * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j]   //c(delta_u0, u_hat_old, lambda_new)  diagonal blocks  ......phictrl_nablauold_uadjnew  
										    - advection_flag * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][kdim]  * phi_gss_fe[ SolFEType[kdim + adj_pos_begin] ][j]   //c(delta_u0, u_0_old, lambda_new)  diagonal blocks  ......phictrl_nablauctrlold_uadjnew  
										    - advection_flag * adv_uctrlold_nablaphictrl_uadjnew_kdim_i_j 	//c(u_0_old, delta_u0, lambda_new)
										      );
                                              
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
               
	      Jac[kdim + ctrl_pos_begin][off_kdim + adj_pos_begin][i*nDofsVadj + j] += AbsDetJxWeight_iqp *  ( - advection_flag * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim]][kdim]  		     * phi_gss_fe[ SolFEType[off_kdim + adj_pos_begin] ][j]   //c(delta_u0, u_hat_old, lambda_new)  off-diagonal blocks  ......phictrl_nablauold_uadjnew  
										         - advection_flag * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i] * gradSolVAR_qp[SolPdeIndex[off_kdim + ctrl_pos_begin]][kdim]  * phi_gss_fe[ SolFEType[off_kdim + adj_pos_begin] ][j]   //c(delta_u0, u_0_old, lambda_new)  off-diagonal blocks  ......phictrl_nablauctrlold_uadjnew  
										      );
	    
	     } //j_dctrl_adj loop
      }
 }
 
//DIAG BLOCK delta_control - control--------------------------------------------------------------------------------------
 for (unsigned  kdim = 0; kdim < dim; kdim++) {
    for (unsigned i = 0; i < nDofsVctrl; i++) {
      for (unsigned j = 0; j < nDofsVctrl; j++) {
          
            double  lap_jac_dctrl_ctrl_kdim_i_j = 0.;
              
            for (unsigned  jdim = 0; jdim < dim_offset_grad /*space_dim*/; jdim++) { 
		lap_jac_dctrl_ctrl_kdim_i_j += phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + jdim]	*	phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j * dim_offset_grad /*space_dim*/ + jdim];
            }

            Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsVctrl + j] += (  + (cost_functional_coeff * target_flag + alpha) * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j]
											+  beta * lap_jac_dctrl_ctrl_kdim_i_j 
											- advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i]    * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]   //c(delta_u0, u_0_new, lambda_old)  diagonal blocks  ......phictrl_nablauctrlnew_uadjold 
											- advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][j] 	* phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_0_new, delta_u0, lambda_old) diagonal blocks  ......uctrlnew_nablaphictrl_uadjold
											) * AbsDetJxWeight_iqp;
                                            
               unsigned int off_kdim = (kdim+1)%dim; //off-diagonal blocks
               
	      Jac[kdim + ctrl_pos_begin][off_kdim + ctrl_pos_begin][i*nDofsVctrl + j] += ( - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i]     * phi_x_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j*dim_offset_grad /*space_dim*/ + kdim] * SolVAR_qp[SolPdeIndex[off_kdim + adj_pos_begin]]   //c(delta_u0, u_0_new, lambda_old)  off-diagonal blocks  ......phictrl_nablauctrlnew_uadjold 
											   - advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[off_kdim + ctrl_pos_begin] ][j] * phi_x_gss_fe[ SolFEType[kdim + ctrl_pos_begin] ][i*dim_offset_grad /*space_dim*/ + off_kdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]]	 //c(u_0_new, delta_u0, lambda_old) off-diagonal blocks  ......uctrlnew_nablaphictrl_uadjold
											) * AbsDetJxWeight_iqp;
	    
      }//j_dctrl_ctrl loop

    }//i_ctrl loop

  }
    
   
  }
      
    else if ( control_el_flag == 0) {
        
    for (unsigned  kdim = 0; kdim < dim; kdim++) {
     for (unsigned i = 0; i < nDofsVctrl; i++) {
      for (unsigned j = 0; j < nDofsVctrl; j++) {
          if (i == j) {
          Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsVctrl + j]  +=  penalty_outside_control_domain * (1 - control_node_flag[kdim][i]);
          }
      }
    }
  }
        
    } 
 //************ Jacobian, END *********************
   
    
    
//BLOCK Pressure_ctrl
  for (unsigned kdim = 0; kdim < dim; kdim++) {
          
for (unsigned i = 0; i < nDofsVctrl; i++) {      
            Res[kdim + ctrl_pos_begin][i] +=  AbsDetJxWeight_iqp * ( SolVAR_qp[SolPdeIndex[press_type_pos + ctrl_pos_begin]] * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim]);   //pressure

      for (unsigned j = 0; j < nDofsPctrl; j++) {
	      Jac[kdim + ctrl_pos_begin][press_type_pos + ctrl_pos_begin][i*nDofsPctrl + j] += -( phi_gss_fe[SolFEType[press_type_pos + ctrl_pos_begin]][j] * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
      }//j_press_ctrl
      
  }//i_ctrl loop
  
 }

  
//DIV_ctrl
		  double div_ctrl_dctrl_qp = 0.;
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
		div_ctrl_dctrl_qp += gradSolVAR_qp[ SolPdeIndex[ctrl_pos_begin + kdim] ][kdim] ;
	  }
	  
  for (unsigned i = 0; i < nDofsPctrl; i++) {
	  Res[press_type_pos + ctrl_pos_begin][i] += ( (div_ctrl_dctrl_qp) * phi_gss_fe[ SolFEType[ctrl_pos_begin + press_type_pos ] ][i] ) * AbsDetJxWeight_iqp;
  }//i_div_ctrl
  
      
 for (unsigned kdim = 0; kdim < dim; kdim++) {
   for (unsigned i = 0; i < nDofsPctrl; i++) {
	  for (unsigned j = 0; j < nDofsVctrl; j++) {
		Jac[press_type_pos + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsVctrl + j] += - ( phi_gss_fe[SolFEType[press_type_pos + ctrl_pos_begin]][i] * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j * dim_offset_grad /*space_dim*/ + kdim] ) * AbsDetJxWeight_iqp;
	  }//j loop
    }//i_div_ctrl
 }
 
//============ delta_control row - END ==================================================================================================
 
//============ delta_mu row - BEGIN  ============================================================================================
  //MU
//************ Residual, BEGIN *********************
  for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
          
  for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; i++) {
      
       Res[mu_pos_begin + kdim][i]  +=  (- penalty_outside_control_domain) *  (1 - control_node_flag[kdim][i]) * (Sol_eldofs_Mat[mu_pos_begin + kdim][i] - 0.);
      
     }
  }
//************ Residual, END *********************

  //MU
//************ Jacobian, BEGIN *********************
  for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
    for (unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; i++) {
      for (unsigned j = 0; j < Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim]; j++) {
            if (i == j) {
               Jac[mu_pos_begin + kdim][mu_pos_begin + kdim][i * Sol_n_el_dofs_Mat_vol[mu_pos_begin + kdim] + j]  +=  penalty_outside_control_domain * (1 - control_node_flag[kdim][i]);
            }
         }
      }
  }
//************ Jacobian, END *********************


//============ delta_mu row - END  ============================================================================================
 
 
 
      }  // end quadrature point loop
 

      //***************************************************************************************************************
      

    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned i_unk = 0; i_unk < n_unknowns; i_unk++) {
      RES->add_vector_blocked(Res[SolPdeIndex[i_unk]],L2G_dofmap_Mat[i_unk]);
        for(unsigned j_unk=0; j_unk < n_unknowns; j_unk++) {
	  if(assembleMatrix) JAC->add_matrix_blocked( Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ], L2G_dofmap_Mat[i_unk], L2G_dofmap_Mat[j_unk]);
        }
    }
 
   //--------------------------------------------------------------------------------------------------------  
  } //end list of elements loop for each subdomain
  
  
  //MU in res ctrl - BEGIN  ***********************************
femus::ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::add_one_times_mu_res_ctrl(iproc,
                               ineq_flag,
                               ctrl_index_in_mat,
                               mu_index_in_mat,
                               SolIndex,
                               sol,
                               mlPdeSys,
                               pdeSys,
                               RES);
  //MU in res ctrl - END ***********************************
    
    
    
     
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
                        Sol_eldofs_Mat,
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
   Sol_eldofs_Mat,
   Sol_n_el_dofs_Mat_vol,
   c_compl,
   mu_index_in_mat,
   ctrl_index_in_mat,
   solIndex_act_flag_sol,
   ctrl_lower,
   ctrl_upper,
   sol_actflag);
  
      


    ctrl::mixed_state_or_ctrl_inequality< femus::ctrl::square_or_cube::mixed_state_or_ctrl_inequality >::node_insertion(iel,
                   msh,
                   L2G_dofmap_Mat,
                   mu_index_in_mat,
                   ctrl_index_in_mat,
                    Sol_eldofs_Mat,
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
  if (assembleMatrix) { 
             for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) { 
JAC->matrix_set_off_diagonal_values_blocked(L2G_dofmap_Mat[ ctrl_index_in_mat[kdim] ], L2G_dofmap_Mat[ mu_index_in_mat[kdim] ], ineq_flag * 1.);
             }
            }
  
  
  }
//   ***************** INSERT PART - END (must go AFTER the sum, clearly) *******************
  
  
  
 RES->close();
 JAC->close();
  // ***************** END ASSEMBLY *******************

    
  
  print_global_residual_jacobian(print_algebra_global,
                                 ml_prob,
                                 mlPdeSys,
                                 pdeSys,
                                 RES,
                                 JAC,
                                 iproc,
                                 assembleMatrix);
  


}
 
    
 
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv 
//i.e. ||prol_(u_(i-1)) - u_(i)|| = err(i) => err(i-1)/err(i) = 2^alpha ,implemented as log(err(i)/err(i+1))/log2

static double*  GetErrorNorm(const MultiLevelProblem& ml_prob, MultiLevelSolution* ml_sol, Solution* sol_coarser_prolongated) {
  
    static double ErrorNormArray[ lifting_internal_norms::no_of_l2_norms + lifting_internal_norms::no_of_h1_norms ];
    
  unsigned level = ml_sol->GetMLMesh()->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = ml_sol->GetMLMesh()->GetLevel(level);    // pointer to the mesh (level) object
  elem*     el  = msh->GetMeshElements();  // pointer to the elem object in msh (level)
  Solution* sol = ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

 // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  std::vector < std::vector < double > > coordX(dim);    // local coordinates

  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(max_size);
  }
   
  //geometry *******************************

 // solution variables *******************************************
  const int n_vars = dim+1;
  const int n_unknowns = 3*n_vars; //(2.*dim)+1; //state , adjoint of velocity terms and one pressure term
  const int vel_type_pos = 0;
  const int press_type_pos = dim;
  const int state_pos_begin = 0;
  const int adj_pos_begin   = dim+1;
  const int ctrl_pos_begin   = 2*(dim+1);
  
  std::vector < std::string > Solname(n_unknowns);  // const char Solname[4][8] = {"u_0","u_1","u_2","p_u"};
  Solname              [state_pos_begin+0] =                "u_0";
  Solname              [state_pos_begin+1] =                "u_1";
  if (dim == 3) Solname[state_pos_begin+2] =                "u_2";
  Solname              [state_pos_begin + press_type_pos] = "p_u";
  
  Solname              [adj_pos_begin + 0] =              "adj_0";
  Solname              [adj_pos_begin + 1] =              "adj_1";
  if (dim == 3) Solname[adj_pos_begin + 2] =              "adj_2";
  Solname              [adj_pos_begin + press_type_pos] = "p_adj";

  Solname              [ctrl_pos_begin + 0] =              "ctrl_0";
  Solname              [ctrl_pos_begin + 1] =              "ctrl_1";
  if (dim == 3) Solname[ctrl_pos_begin + 2] =              "ctrl_2";
  Solname              [ctrl_pos_begin + press_type_pos] = "p_ctrl";
  
  std::vector < unsigned > SolIndex(n_unknowns);  
  std::vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

  std::vector < double > Sol_n_el_dofs_Mat_vol(n_unknowns);
  
  //==========================================================================================
  // velocity ************************************
  std::vector < std::vector < double > > phi_gss_fe(NFE_FAMS);
  std::vector < std::vector < double > > phi_x_gss_fe(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(max_size);
      phi_x_gss_fe[fe].reserve(max_size*dim);
   }
  
  //=================================================================================================
  
  // quadratures ********************************
  double AbsDetJxWeight_iqp;
  
  
  //----------- dofs ------------------------------
  std::vector < std::vector < double > > Sol_eldofs_Mat(n_unknowns);
  std::vector < std::vector < double > > gradSol_eldofs_Mat(n_unknowns);
  
  std::vector < std::vector < double > > SolVAR_coarser_prol_eldofs(n_unknowns);
  std::vector < std::vector < double > > gradSolVAR_coarser_prol_eldofs(n_unknowns);


  for(int k = 0; k < n_unknowns; k++) {
    Sol_eldofs_Mat[k].reserve(max_size);
    gradSol_eldofs_Mat[k].reserve(max_size*dim); 
    
    SolVAR_coarser_prol_eldofs[k].reserve(max_size);
    gradSolVAR_coarser_prol_eldofs[k].reserve(max_size*dim);    
  }

  //------------ at quadrature points ---------------------
  std::vector < double > SolVAR_qp(n_unknowns);
  std::vector < double > SolVAR_coarser_prol_qp(n_unknowns);
  std::vector < std::vector < double > > gradSolVAR_qp(n_unknowns);
  std::vector < std::vector < double > > gradSolVAR_coarser_prol_qp(n_unknowns);
  for(int k = 0; k < n_unknowns; k++) {
      gradSolVAR_qp[k].reserve(max_size);  
      gradSolVAR_coarser_prol_qp[k].reserve(max_size);  
  }
      
  std::vector < double > l2norm ( lifting_internal_norms::no_of_l2_norms ,0.);
  std::vector < double > seminorm ( lifting_internal_norms::no_of_h1_norms ,0.);

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    
  // geometry *****************************
    short unsigned ielGeom = msh->GetElementType(iel);
    
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    for (unsigned  k = 0; k < dim; k++) {       coordX[k].resize(nDofsX);    }
  
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->GetTopology()->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }
    
      // elem average point 
    std::vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
  for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofsX; i++) {
         elem_center[j] += coordX[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofsX; }
  //*************************************** 
  
  // geometry end *****************************
  
  
 // equation *****************************
    unsigned nDofsV = msh->GetElementDofNumber(iel, SolFEType[vel_type_pos]);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType[state_pos_begin + press_type_pos]);    // number of solution element dofs
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin]);    // number of solution element dofs
    unsigned nDofsPadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin + press_type_pos]);    // number of solution element dofs

    unsigned nDofsVctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin]);    // number of solution element dofs
    unsigned nDofsPctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin + press_type_pos] );    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + nDofsP;
    unsigned nDofsVP_tot = 3*nDofsVP;
  // equation end *****************************


   //STATE###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	Sol_n_el_dofs_Mat_vol[k]=ndofs_unk;
       Sol_eldofs_Mat[k].resize(ndofs_unk);
       SolVAR_coarser_prol_eldofs[k].resize(ndofs_unk);
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
       Sol_eldofs_Mat[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
       SolVAR_coarser_prol_eldofs[k][i] = (*sol_coarser_prolongated->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
      }
    }
  //CTRL###################################################################

 
      // ********************** Gauss point loop *******************************
      for(unsigned iqp = 0;iqp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); iqp++) {
	
 
      for(int fe=0; fe < NFE_FAMS; fe++) {
	msh->_finiteElement[ielGeom][fe]->Jacobian(coordX,iqp, AbsDetJxWeight_iqp,phi_gss_fe[fe],phi_x_gss_fe[fe], boost::none);
      }
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
  	msh->_finiteElement[ielGeom][CONTINUOUS_BIQUADRATIC]->Jacobian(coordX,iqp,AbsDetJxWeight_iqp,phi_gss_fe[CONTINUOUS_BIQUADRATIC],phi_x_gss_fe[CONTINUOUS_BIQUADRATIC], boost::none);

 //begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk < n_unknowns; unk++) {
	  SolVAR_qp[unk] = 0.;
	  SolVAR_coarser_prol_qp[unk] = 0.;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_qp[unk][ivar2] = 0.; 
	    gradSolVAR_coarser_prol_qp[unk][ivar2] = 0.; 
	  }
    }
	  
	for(unsigned unk = 0; unk <  n_unknowns; unk++) {
	  for(unsigned i = 0; i < Sol_n_el_dofs_Mat_vol[unk]; i++) {
	    SolVAR_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i] * Sol_eldofs_Mat[unk][i];
	    SolVAR_coarser_prol_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i] * SolVAR_coarser_prol_eldofs[unk][i];
//         std::cout << SolVAR_qp[unk] << " \t " << SolVAR_coarser_prol_qp[unk] << std::endl;
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim+ivar2] * Sol_eldofs_Mat[unk][i]; 
	      gradSolVAR_coarser_prol_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim+ivar2] * SolVAR_coarser_prol_eldofs[unk][i]; 
//         std::cout << gradSolVAR_qp[unk][ivar2] << " \t " << gradSolVAR_coarser_prol_qp[unk][ivar2] << std::endl;
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************


	for(unsigned unk = 0; unk < n_unknowns; unk++) {
        l2norm[unk] += ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * AbsDetJxWeight_iqp ; 
        
     }
    
    for(int  i = 0; i < dim; i++) {

        l2norm[n_unknowns + i] += ( ( SolVAR_qp[vel_type_pos + i] + SolVAR_qp[ctrl_pos_begin + i] ) - ( SolVAR_coarser_prol_qp[vel_type_pos + i]  + SolVAR_coarser_prol_qp[ctrl_pos_begin + i] ) ) * ( ( SolVAR_qp[vel_type_pos + i] + SolVAR_qp[ctrl_pos_begin + i] ) - ( SolVAR_coarser_prol_qp[vel_type_pos + i] + SolVAR_coarser_prol_qp[ctrl_pos_begin + i] ) )  * AbsDetJxWeight_iqp ;
    
    }

          
	for(unsigned unk = 0; unk < dim; unk++) {
        for(int j = 0; j < dim; j++){
        seminorm[unk] += (gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * ( gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * AbsDetJxWeight_iqp ;
        seminorm[unk + dim] += (gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * ( gradSolVAR_qp[unk + adj_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + adj_pos_begin][j] ) * AbsDetJxWeight_iqp ;
        seminorm[unk + 2*dim] += (gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * ( gradSolVAR_qp[unk + ctrl_pos_begin][j] - gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j] ) * AbsDetJxWeight_iqp ;
        seminorm[unk + 3*dim] += ((gradSolVAR_qp[unk][j]+gradSolVAR_qp[unk + ctrl_pos_begin][j]) - (gradSolVAR_coarser_prol_qp[unk][j]+gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j])) * ((gradSolVAR_qp[unk][j]+gradSolVAR_qp[unk + ctrl_pos_begin][j]) - (gradSolVAR_coarser_prol_qp[unk][j]+gradSolVAR_coarser_prol_qp[unk + ctrl_pos_begin][j])) * AbsDetJxWeight_iqp ;
        }
     }
    
           
    } // end gauss point loop
  } //end element loop for each process


    // add the norms of all processes
  NumericVector* norm_vec_inexact;
  norm_vec_inexact = NumericVector::build().release();
  norm_vec_inexact->init(msh->n_processors(), 1 , false, AUTOMATIC);

	for(unsigned unk = 0; unk <  lifting_internal_norms::no_of_l2_norms ; unk++) {
        norm_vec_inexact->set(iproc, l2norm[unk]);
        norm_vec_inexact->close();
        l2norm[unk] = norm_vec_inexact->l1_norm();
    }

	for(unsigned unk = 0; unk <  lifting_internal_norms::no_of_h1_norms ; unk++) {
        norm_vec_inexact->set(iproc, seminorm[unk]);
        norm_vec_inexact->close();
        seminorm[unk] = norm_vec_inexact->l1_norm();
    }

  delete norm_vec_inexact;
  
 
	for(unsigned unk = 0; unk <  lifting_internal_norms::no_of_l2_norms ; unk++) {
        ErrorNormArray[unk] = sqrt(l2norm[unk]);
    }
	for(unsigned unk = 0; unk <  lifting_internal_norms::no_of_h1_norms ; unk++) {
        ErrorNormArray[unk +  lifting_internal_norms::no_of_l2_norms ] = sqrt(seminorm[unk]);
    }
   
   return ErrorNormArray;
  
  
}

    
    
};

    
    
}


}





#endif
