#ifndef __00_COST_FUNCTIONAL_HPP__
#define __00_COST_FUNCTIONAL_HPP__



#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "FElemTypeEnum_list.hpp"
#include "Parallel.hpp"

#include "PetscVector.hpp"


#include "opt_common.hpp"

#include "00_cost_functional_without_regularization.hpp"





namespace femus {
    

 
//************** how to retrieve theta from proc0 ************************************* 
static const double get_theta_value(const unsigned int nprocs, const Solution * sol, const unsigned int sol_theta_index) {
    
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


namespace ctrl {



    
template < class LIST_OF_CTRL_FACES, class DOMAIN_CONTAINING_CTRL_FACES, class COST_WITHOUT_REG >
class cost_functional {


public:
  
    

  
  /** This function computes a functional with a volume part and a boundary part
     We pass a 2 Solution objects: the first for the cost functional, the second for the regularization
    */

static void compute_cost_functional_regularization_bdry(const MultiLevelProblem & ml_prob, 
                     const unsigned level, 
                     const unsigned iteration,
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars,
                     const unsigned alpha_in,
                     const unsigned beta_in,
                     const unsigned qrule_order_in
                    )  {
    
  std::cout << "=== Compute cost functional parts with boundary regularization =============" << std::endl;  
  
  
  
  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->GetMeshElements();

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //=============== Geometry ========================================
   unsigned solType_coords = CONTINUOUS_BIQUADRATIC;
 
  CurrentElem < double > geom_element_iel(dim, msh);
    
  constexpr unsigned int space_dim = 3;
  
  std::vector<double> normal(space_dim, 0.);
 //***************************************************

  //=============== Integration ========================================

 //***************************************************
  const double alpha = alpha_in;
  const double beta  = beta_in;
  
 //*************** state ***************************** 
 //***************************************************
  const unsigned n_components_state = state_vars.size();
  
  std::vector< unsigned > solIndex_u(n_components_state);
  std::vector< unsigned > solType_u(n_components_state);

  for (unsigned c = 0; c < solIndex_u.size(); c++) { 
     solIndex_u[c] = ml_sol->GetIndex( state_vars[c].c_str() );
     solType_u[c]  = ml_sol->GetSolutionType(solIndex_u[c]);
   }




  std::vector <double> phi_u;     phi_u.reserve(max_size);
  std::vector <double> phi_u_x;   phi_u_x.reserve(max_size * space_dim);
//   std::vector <double> phi_u_xx;  phi_u_xx.reserve(max_size * dim2);
 
  std::vector < double >  sol_u;
  sol_u.reserve(max_size);
  
  double u_gss = 0.;
  double u_x_gss = 0.;
 //*************************************************** 
 //***************************************************

  
 //************** cont *******************************
 //***************************************************
  const unsigned n_components_ctrl = ctrl_vars.size();
  
  std::vector <double> phi_ctrl_bdry;  
  std::vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(max_size);
  phi_ctrl_x_bdry.reserve(max_size * space_dim);

  
  std::vector< unsigned > solIndex_ctrl(n_components_ctrl);
  std::vector< unsigned > solType_ctrl(n_components_ctrl);
  
  for (unsigned c = 0; c < solIndex_ctrl.size(); c++) { 
  solIndex_ctrl[c] = ml_sol->GetIndex( ctrl_vars[c].c_str() );
  solType_ctrl[c] = ml_sol->GetSolutionType(solIndex_ctrl[c]);
  }
  

   std::vector < double >  sol_ctrl;   sol_ctrl.reserve(max_size);
 //***************************************************
 //*************************************************** 
  
  
 //************** desired ****************************
 //***************************************************
  std::vector <double> phi_udes;
  std::vector <double> phi_udes_x;

    phi_udes.reserve(max_size);
    phi_udes_x.reserve(max_size * space_dim);
 

  std::vector < double >  sol_udes;
  sol_udes.reserve(max_size);

  double udes_gss = 0.;
 //***************************************************
 //***************************************************

 //********** DATA *********************************** 
  double u_des = COST_WITHOUT_REG:: DesiredTargetVec()[0];
 //*************************************************** 
  
  double integral_target = 0.;
  double integral_alpha  = 0.;
  double integral_beta   = 0.;

  

 //*************************************************** 
  const unsigned qrule_i = qrule_order_in;
  
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
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element_iel.geom_type();

  //************* set target domain flag **************
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = COST_WITHOUT_REG::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
 //***************************************************

   
 //*********** state ********************************* 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u[0]);
    sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
      unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u[0]);
      sol_u[i] = (*sol->_Sol[solIndex_u[0] ])(solDof_u);
    }
 //*********** state ********************************* 


 //*********** cont ********************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl[0]);
    sol_ctrl    .resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl[0]);
      sol_ctrl[i] = (*sol->_Sol[ solIndex_ctrl[0] ])(solDof_ctrl);
    } 

 //*********** cont ********************************** 
 
 
 //*********** udes ********************************** 
    unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u[0]);
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
  
// // // 	if ( femus::ctrl::Gamma_control::volume_elem_contains_a_Gamma_control_face( geom_element_iel.get_elem_center_3d() ) ) {
	  
	       
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       const unsigned nve_bdry_ctrl = msh->GetElementFaceDofNumber(iel,iface,solType_ctrl[0]);
       
// ----------
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();
// ----------

          std::pair< int, unsigned int > pair_control_iface = femus::face_is_a_Gamma_control_face_of_some_index< LIST_OF_CTRL_FACES >(msh->GetMeshElements(), iel, iface);

          const int  iface_is_a_boundary_control  = pair_control_iface.first;

	    if( iface_is_a_boundary_control ) {
	
		//============ initialize gauss quantities on the boundary ==========================================
                double sol_ctrl_bdry_gss = 0.;
                std::vector<double> sol_ctrl_x_bdry_gss(space_dim);
		//============ initialize gauss quantities on the boundary ==========================================
		
		for(unsigned ig_bdry = 0; ig_bdry < ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussPointsNumber(); ig_bdry++) {
		  
    elem_all[qrule_i][ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_qp_bdry, JacI_qp_bdry, detJac_iqp_bdry, space_dim);
    weight_iqp_bdry = detJac_iqp_bdry * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];
    elem_all[qrule_i][ielGeom_bdry][solType_ctrl[0]] ->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_ctrl_bdry, phi_ctrl_x_bdry, boost::none, space_dim);

		  
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
	  
// // // 	} //end if control element flag

//=================== BOUNDARY PART - END ==================================================================================================  

  
  
   
//=================== VOLUME PART - BEGIN ==================================================================================================  

   for (unsigned ig = 0; ig < ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom).GetGaussPointsNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
    elem_all[qrule_i][ielGeom][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_3d(), ig, Jac_qp, JacI_qp, detJac_iqp, space_dim);
    weight_iqp = detJac_iqp * ml_prob.GetQuadratureRuleMultiple(qrule_i, ielGeom).GetGaussWeightsPointer()[ig];

    elem_all[qrule_i][ielGeom][solType_u[0]]                 ->shape_funcs_current_elem(ig, JacI_qp, phi_u, phi_u_x, boost::none, space_dim);
    elem_all[qrule_i][ielGeom][solType_u[0]/*solTypeTdes*/]  ->shape_funcs_current_elem(ig, JacI_qp, phi_udes, phi_udes_x, boost::none, space_dim);
    
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
  

  
  

static void compute_cost_functional_regularization_lifting_internal(const MultiLevelProblem & ml_prob, 
                     const unsigned level, 
                     const unsigned iteration,
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars,  
                    const unsigned alpha_in,
                     const unsigned beta_in,
                     const unsigned qrule_order_in

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
 unsigned solType_coords = CONTINUOUS_BIQUADRATIC;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
 //*************************************************** 

  
 //***************************************************  
  double alpha = alpha_in;
  double beta  = alpha_in;

 //******************** state ************************ 
 //*************************************************** 
  const unsigned n_components_state = state_vars.size();

  std::vector <double> phi_u;
  std::vector <double> phi_u_x;
  std::vector <double> phi_u_xx;

  phi_u.reserve(max_size);
  phi_u_x.reserve(max_size * dim);
  phi_u_xx.reserve(max_size * dim2);
  
 
  std::vector< unsigned > solIndex_u(n_components_state);
  std::vector< unsigned > solType_u(n_components_state);

  for (unsigned c = 0; c < solIndex_u.size(); c++) { 
     solIndex_u[c] = ml_sol->GetIndex( state_vars[c].c_str() );
     solType_u[c]  = ml_sol->GetSolutionType(solIndex_u[c]);
   }

  std::vector < double >  sol_u; // local solution
  sol_u.reserve(max_size);
  
  double u_gss = 0.;
double u_x_gss = 0.;
 //*************************************************** 

 //******************** control ********************** 
 //*************************************************** 
  const unsigned n_components_ctrl = ctrl_vars.size();

  std::vector <double> phi_ctrl;
  std::vector <double> phi_ctrl_x;
  std::vector <double> phi_ctrl_xx;

  phi_ctrl.reserve(max_size);
  phi_ctrl_x.reserve(max_size * dim);
  phi_ctrl_xx.reserve(max_size * dim2);
  
  std::vector< unsigned > solIndex_ctrl(n_components_ctrl);
  std::vector< unsigned > solType_ctrl(n_components_ctrl);
  
  for (unsigned c = 0; c < solIndex_ctrl.size(); c++) { 
  solIndex_ctrl[c] = ml_sol->GetIndex( ctrl_vars[c].c_str() );
  solType_ctrl[c] = ml_sol->GetSolutionType(solIndex_ctrl[c]);
  }
  
  std::vector < double >  sol_ctrl; // local solution
  sol_ctrl.reserve(max_size);
  
  double ctrl_gss = 0.;
  double ctrl_x_gss = 0.;
 //*************************************************** 
 //***************************************************  

  
 //********************* desired ********************* 
 //*************************************************** 
  std::vector <double> phi_udes;
  std::vector <double> phi_udes_x;
  std::vector <double> phi_udes_xx;

  phi_udes.reserve(max_size);
  phi_udes_x.reserve(max_size * dim);
  phi_udes_xx.reserve(max_size * dim2);
 
  
//  unsigned solIndex_udes;
//   solIndex_udes = ml_sol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solType_udes = ml_sol->GetSolutionType(solIndex_udes);    // get the finite element type for "state"

  std::vector < double >  sol_udes; // local solution
  sol_udes.reserve(max_size);

  double udes_gss = 0.;
 //*************************************************** 
 //*************************************************** 

 //********************* DATA ************************ 
  double u_des = COST_WITHOUT_REG::DesiredTargetVec()[0];
 //*************************************************** 
  
  double integral_target = 0.;
  double integral_alpha  = 0.;
  double integral_beta   = 0.;


 //***************************************************  
       //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < std::vector < /*const*/ elem_type_templ_base<double, double> *  > > > elem_all;
  ml_prob.get_all_abstract_fe_multiple(elem_all);
  
   const unsigned qrule_i = qrule_order_in;
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
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {
 
    geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    const short unsigned ielGeom = geom_element_iel.geom_type();

  //************* set target domain flag **************
   geom_element_iel.set_elem_center_3d(iel, solType_coords);

   int target_flag = 0;
   target_flag = COST_WITHOUT_REG::ElementTargetFlag( geom_element_iel.get_elem_center_3d() );
 //*************************************************** 

 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag = DOMAIN_CONTAINING_CTRL_FACES ::ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
 //*************************************************** 
   
 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u[0] );
        sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u[0] );
            sol_u[i] = (*sol->_Sol[ solIndex_u[0] ])(solDof_u);
    }
 //*************************************************** 


 //************** control **************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl[0]);
    sol_ctrl    .resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl[0]);
      sol_ctrl[i] = (*sol->_Sol[ solIndex_ctrl[0] ])(solDof_ctrl);
    } 
 //***************************************************  
 
 
 //**************** u_des **************************** 
    unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u[0] );
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

    
    elem_all[qrule_i][ielGeom][ solType_u[0] ]  ->shape_funcs_current_elem(iqp, JacI_iqp, phi_u, phi_u_x, boost::none, space_dim);
    elem_all[qrule_i][ielGeom][ solType_ctrl[0] ]  ->shape_funcs_current_elem(iqp, JacI_iqp, phi_ctrl, phi_ctrl_x, boost::none, space_dim);
    elem_all[qrule_i][ielGeom][ solType_u[0] ]  ->shape_funcs_current_elem(iqp, JacI_iqp, phi_udes, phi_udes_x, boost::none, space_dim);
    
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
  


static void compute_cost_functional_regularization_lifting_external(const MultiLevelProblem& ml_prob, 
                     const unsigned level, 
                     const unsigned iteration,
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars,  
                    const unsigned alpha_in,
                     const unsigned beta_in,
                     const unsigned qrule_order_in,
                     const unsigned flag_for_internal_group
)    {



    Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);            // pointer to the mesh (level) object

    MultiLevelSolution*    ml_sol = ml_prob._ml_sol;                             // pointer to the multilevel solution object
    Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

    const unsigned           dim = msh->GetDimension();                         // get the domain dimension of the problem
    const unsigned          dim2 = (3 * (dim - 1) + !(dim - 1));                // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
    const unsigned       max_size = static_cast< unsigned >(ceil(pow(3, dim)));  // conservative: based on line3, quad9, hex27

    const unsigned         iproc = msh->processor_id();                         // get the process_id (for parallel computation)

//***************************************************
    std::vector < std::vector < double > > coords_at_dofs(dim);   // local coordinates
    unsigned solType_coords  = CONTINUOUS_BIQUADRATIC;  // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
    for (unsigned i = 0; i < dim; i++) {
        coords_at_dofs[i].reserve(max_size);
    }
//***************************************************

//***************************************************
    double weight_qp; // gauss point weight

//***************************************************
    double alpha = alpha_in;
    double beta  = beta_in;

//******************** state ************************
//***************************************************
    std::vector <double> phi_u;    // local test function
    std::vector <double> phi_u_x;  // local test function first order partial derivatives
    std::vector <double> phi_u_xx; // local test function second order partial derivatives

    phi_u.reserve(max_size);
    phi_u_x.reserve(max_size * dim);
    phi_u_xx.reserve(max_size * dim2);


    unsigned solIndex_u;
    solIndex_u = ml_sol->GetIndex("state");                    // get the position of "state" in the ml_sol object
    unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);  // get the finite element type for "state"

    std::vector < double >  sol_u; // local solution
    sol_u.reserve(max_size);

    double u_gss = 0.;
//***************************************************
//***************************************************

//******************** control **********************
//***************************************************
    std::vector <double> phi_ctrl;    // local test function
    std::vector <double> phi_ctrl_x;  // local test function first order partial derivatives
    std::vector <double> phi_ctrl_xx; // local test function second order partial derivatives

    phi_ctrl.reserve(max_size);
    phi_ctrl_x.reserve(max_size * dim);
    phi_ctrl_xx.reserve(max_size * dim2);

    unsigned solIndex_ctrl;
    solIndex_ctrl = ml_sol->GetIndex("control");
    unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

    std::vector < double >  sol_ctrl; // local solution
    sol_ctrl.reserve(max_size);

    double ctrl_gss = 0.;
    double ctrl_x_gss = 0.;
//***************************************************
//***************************************************


//********************* desired *********************
//***************************************************
    std::vector <double> phi_udes;    // local test function
    std::vector <double> phi_udes_x;  // local test function first order partial derivatives
    std::vector <double> phi_udes_xx; // local test function second order partial derivatives

    phi_udes.reserve(max_size);
    phi_udes_x.reserve(max_size * dim);
    phi_udes_xx.reserve(max_size * dim2);

    std::vector < double >  sol_udes; // local solution
    sol_udes.reserve(max_size);

    double udes_gss = 0.;
//***************************************************
//***************************************************

//********************* DATA ************************
    double u_des = COST_WITHOUT_REG::DesiredTargetVec()[0];
//***************************************************

    double integral_target = 0.;
    double integral_alpha  = 0.;
    double integral_beta   = 0.;


    // element loop: each process loops only on the elements that owns
    for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

        int group_flag         = msh->GetElementGroup(iel);      // element group flag (Exterior = GROUP_EXTERNAL, Interior = GROUP_INTERNAL)
        short unsigned ielGeom = msh->GetElementType(iel);       // element geometry type

//***************** GEOMETRY ************************
        unsigned nDofx = msh->GetElementDofNumber(iel, solType_coords); // number of coordinate element dofs
        for (int i = 0; i < dim; i++)  coords_at_dofs[i].resize(nDofx);
        // local storage of coordinates
        for (unsigned i = 0; i < nDofx; i++) {
            unsigned xDof  = msh->GetSolutionDof(i, iel, solType_coords); // global to global mapping between coordinates node and coordinate dof

            for (unsigned jdim = 0; jdim < dim; jdim++) {
                coords_at_dofs[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);  // global extraction and local storage for the element coordinates
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

//************** set target domain flag *************
        int target_flag = 0;
        target_flag = COST_WITHOUT_REG::ElementTargetFlag(elem_center);
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
            integral_alpha  += (group_flag - flag_for_internal_group) * weight_qp * ctrl_gss * ctrl_gss;
            integral_beta   += (group_flag - flag_for_internal_group) * weight_qp * ctrl_x_gss * ctrl_x_gss;

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

  

  
  
  

  
  

static void compute_cost_functional_regularization_bdry_vec(const MultiLevelProblem& ml_prob, 
                     const unsigned level,
                     const unsigned iteration,
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars,  
                    const unsigned alpha_in,
                     const unsigned beta_in,
                     const unsigned cost_functional_coeff_in,
                     const unsigned qrule_order_in
 
) {

    
    
  const double cost_functional_coeff = cost_functional_coeff_in;
  const double alpha = alpha_in;
  const double beta  = beta_in;

  
  
  const Mesh*          msh          	= ml_prob._ml_msh->GetLevel(level); 
  elem*          el         	= msh->GetMeshElements();

  MultiLevelSolution*  ml_sol    = ml_prob._ml_sol;
  Solution*    sol        	= ml_prob._ml_sol->GetSolutionLevel(level);
  
  unsigned    iproc = msh->processor_id();
  
  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1)); 
  
  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  //geometry *******************************
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  unsigned solType_coords = coordXType;
 
  CurrentElem < double > geom_element_iel(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = /*dim*/  3  /*2*/    ;
 
  std::vector<double> normal_iqp(dim_offset_grad /*space_dim*/, 0.);

  std::vector < std::vector < double > > coordX(dim);    // local coordinates
  std::vector< std::vector < double> > coordX_bd(dim);

  for (unsigned  k = 0; k < dim; k++) { 
        coordX[k].reserve(max_size);
        coordX_bd[k].reserve(max_size); 
  }
  
  double AbsDetJxWeight_iqp = 0.;
  double AbsDetJxWeight_iqp_bdry = 0.;
  
  
  //geometry *******************************

//STATE######################################################################
  //velocity *******************************
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = ml_sol->GetIndex("u_0");
  solVIndex[1] = ml_sol->GetIndex("u_1");

  if (dim == 3) solVIndex[2] = ml_sol->GetIndex("u_2");

  unsigned solVType = ml_sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
  
  std::vector < std::vector < double > >  solV(dim);    // local solution
  std::vector <double >  V_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(max_size);
  }

  
  std::vector <double> phiV_gss;  // local test function
  std::vector <double> phiV_x_gss; // local test function first order partial derivatives

  phiV_gss.reserve(max_size);
  phiV_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
   
//STATE######################################################################
  

//CONTROL_@bdry######################################################################
  std::vector < unsigned > solVctrlIndex(dim);
  solVctrlIndex[0] = ml_sol->GetIndex("ctrl_0");
  solVctrlIndex[1] = ml_sol->GetIndex("ctrl_1");
  if (dim == 3) solVctrlIndex[2] = ml_sol->GetIndex("ctrl_2");

  unsigned solVctrlType = ml_sol->GetSolutionType(solVctrlIndex[0]);    // get the finite element type for "u"
  
  std::vector < std::vector < double > >  solVctrl(dim);    // local solution
  std::vector < double >   Vctrl_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solVctrl[k].reserve(max_size);
  }

  
  std::vector <double> phiVctrl_gss_bd;  // local test function
  std::vector <double> phiVctrl_x_gss_bd; // local test function first order partial derivatives

  phiVctrl_gss_bd.reserve(max_size);
  phiVctrl_x_gss_bd.reserve(max_size * dim_offset_grad );
  
//CONTROL_@bdry######################################################################
  
//Theta value ######################################################################
   const unsigned solThetaIndex = ml_sol->GetIndex("theta");
   const unsigned solThetaType = ml_sol->GetSolutionType(solThetaIndex);
   
//    double solTheta = (*sol->_Sol[solThetaIndex])(0)/*0.*/;
   //************** how to retrieve theta from proc0 ************************************* 
 double solTheta = femus::get_theta_value(msh->n_processors(), sol, solThetaIndex);
//*************************************************** 
// 		     solTheta = (*sol->_Sol[solThetaIndex])(0);
//Theta value ######################################################################


// Vel_desired##################################################################
  std::vector <double> phiVdes_gss;  // local test function
  std::vector <double> phiVdes_x_gss; // local test function first order partial derivatives

  phiVdes_gss.reserve(max_size);
  phiVdes_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);

//   std::vector< std::vector < double > >  solVdes(dim);    // local solution
  std::vector <double>  solVdes(dim,0.);
  std::vector<double> Vdes_gss(dim, 0.);  
  
//  for (unsigned  k = 0; k < dim; k++) {
//     solVdes[k].reserve(max_size);
//   }
//   


// Vel_desired##################################################################

  
std::vector<double> integral(dim);

double  integral_target_alpha = 0.;

double	integral_beta   = 0.;
double	integral_gamma  = 0.;

double integral_g_dot_n = 0.;
  

// double	integral_div_ctrl  = 0.;

//*************************************************** 
  //--- quadrature rules -------------------
  const unsigned qrule_i = qrule_order_in;
  
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

  
  

  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

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
   target_flag = COST_WITHOUT_REG::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
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
        solVdes[k]/*[i]*/ = COST_WITHOUT_REG::DesiredTargetVec()[k] /*(*sol->_Sol[solVIndex[k]])(solVdesDof)*/;      // global extraction and local storage for the solution
     }
//     }
 //DESIRED VEL###################################################################

 

//========BoundaryLoop=====================================================================

// // // 	if ( femus::ctrl::Gamma_control::volume_elem_contains_a_Gamma_control_face( geom_element_iel.get_elem_center_3d() ) ) {

        
    std::vector<double> normal_iqp(dim_offset_grad /*space_dim*/, 0.);
	  
    for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {

       const unsigned ielGeom_bd = msh->GetElementFaceType(iel, iface);
       
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
       geom_element_iel.set_elem_center_bdry_3d();

       std::pair< int, unsigned int > pair_control_iface = femus::face_is_a_Gamma_control_face_of_some_index< LIST_OF_CTRL_FACES >(msh->GetMeshElements(), iel, iface);

          const int iface_is_a_boundary_control  = pair_control_iface.first;

	    if( iface_is_a_boundary_control ) {

	  
//=================================================== 
		
//========= initialize gauss quantities on the boundary ============================================
    std::vector < double >   Vctrl_bd_qp(dim, 0.);    //  solution@bdry
    std::vector < std::vector < double > > gradVctrl_bd_qp(dim);
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
                 const unsigned ndof_bdry = msh->GetElementFaceDofNumber(iel, iface, solVctrlType);

		   for(int i_bd = 0; i_bd < ndof_bdry; i_bd++) {
		       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bd);
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

        
    }  // loop over element faces //iface
      
// // //   } //end if control element flag

  
  
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
  if (Parallel::get_rank() == 0 ) {
      
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




static void compute_cost_functional_regularization_lifting_internal_vec(
                     const MultiLevelProblem & ml_prob, 
                     const unsigned level,
                     const unsigned iteration,
                     const std::vector<std::string> state_vars,  
                     const std::vector<std::string> ctrl_vars,  
                    const unsigned alpha_in,
                     const unsigned beta_in,
                     const unsigned cost_functional_coeff_in
) {

    
    
  const double cost_functional_coeff = cost_functional_coeff_in;
  const double alpha = alpha_in;
  const double beta  = beta_in;

  
  Mesh*          msh          	= ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         	= msh->GetMeshElements();  // pointer to the elem object in msh (level)

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
 
  std::vector < std::vector < double > > coordX(dim);    // local coordinates

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(max_size);
  }
  
  double AbsDetJxWeight_iqp;
  
  
  //geometry *******************************

//STATE######################################################################
  std::vector < unsigned > solVIndex(dim);
  solVIndex[0] = ml_sol->GetIndex("u_0");
  solVIndex[1] = ml_sol->GetIndex("u_1");

  if (dim == 3) solVIndex[2] = ml_sol->GetIndex("u_2");

  unsigned solVType = ml_sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
  
  std::vector < std::vector < double > >  solV(dim);    // local solution
  std::vector <double >  V_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(max_size);
  }

  
  std::vector <double> phiV_gss;  // local test function
  std::vector <double> phiV_x_gss; // local test function first order partial derivatives

  phiV_gss.reserve(max_size);
  phiV_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  
//STATE######################################################################
  

//CONTROL######################################################################
  std::vector < unsigned > solVctrlIndex(dim);
  solVctrlIndex[0] = ml_sol->GetIndex("ctrl_0");
  solVctrlIndex[1] = ml_sol->GetIndex("ctrl_1");

  if (dim == 3) solVctrlIndex[2] = ml_sol->GetIndex("ctrl_2");

  unsigned solVctrlType = ml_sol->GetSolutionType(solVctrlIndex[0]);
  
  std::vector < std::vector < double > >  solVctrl(dim);    // local solution
  std::vector < double >   Vctrl_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solVctrl[k].reserve(max_size);
  }

  
  std::vector <double> phiVctrl_gss;  // local test function
  std::vector <double> phiVctrl_x_gss; // local test function first order partial derivatives

  phiVctrl_gss.reserve(max_size);
  phiVctrl_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);
  
//CONTROL######################################################################

// Vel_desired##################################################################
  std::vector <double> phiVdes_gss;  // local test function
  std::vector <double> phiVdes_x_gss; // local test function first order partial derivatives

  phiVdes_gss.reserve(max_size);
  phiVdes_x_gss.reserve(max_size * dim_offset_grad /*space_dim*/);

  std::vector <double>  solVdes(dim,0.);
  std::vector<double> Vdes_gss(dim, 0.);  
  
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
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

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
   target_flag = COST_WITHOUT_REG::ElementTargetFlag(geom_element_iel.get_elem_center_3d());
//***************************************       
    
 //***** set control flag ****************************
  int control_el_flag = 0;
  control_el_flag =  DOMAIN_CONTAINING_CTRL_FACES :: ControlDomainFlag_internal_restriction(geom_element_iel.get_elem_center_3d());
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
        solVdes[k]/*[i]*/ =COST_WITHOUT_REG::DesiredTargetVec()[k] /*(*sol->_Sol[solVIndex[k]])(solVdesDof)*/;      // global extraction and local storage for the solution
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

	
	  std::vector < std::vector < double > > gradVctrl_gss(dim);
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
  if (Parallel::get_rank() == 0 ) {
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


};



}


}



#endif

