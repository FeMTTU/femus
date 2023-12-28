#ifndef __femus_00_laplacian_with_all_dirichlet_bc_hpp__
#define __femus_00_laplacian_with_all_dirichlet_bc_hpp__
 

#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"
#include "CurrentElem.hpp"
 
using namespace femus;


/**
 * This function assemble the stiffness matrix KK and the residual vector Res
 * for the Newton iterative scheme
 *                  J(u0) w = Res(u_0) = f(x) - J u_0  ,
 *                  with u = u0 + w
 *                  J = \grad_u F
 * 
 * 
 **/


template < class system_type, class real_num, class real_num_mov >
void System_assemble_flexible_Laplacian_With_Manufactured_Sol(
                              const std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > > & elem_all,
                              const std::vector < std::vector < /*const*/ elem_type_templ_base<real_num_mov, real_num_mov> *  > > & elem_all_for_domain,
                              const std::vector<Gauss> & quad_rules,
                              system_type * mlPdeSys,
                              MultiLevelMesh * ml_mesh_in,
                              MultiLevelSolution * ml_sol_in,
                              const std::vector< Unknown > &  unknowns,
                              const std::vector< Math::Function< double > * > & exact_sol) {

    //  level is the level of the PDE system to be assembled
    //  levelMax is the Maximum level of the MultiLevelProblem


    const unsigned level = mlPdeSys->GetLevelToAssemble();
    const bool 			assembleMatrix 		    = mlPdeSys->GetAssembleMatrix();

    Mesh*                    msh = ml_mesh_in->GetLevel(level);

    MultiLevelSolution*   ml_sol = ml_sol_in;
    Solution*                sol = ml_sol->GetSolutionLevel(level);

    LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
    SparseMatrix*             KK = pdeSys->_KK;
    NumericVector*           RES = pdeSys->_RES;

    const unsigned  dim = msh->GetDimension();

    const unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
    
    

    RES->zero();
    if (assembleMatrix)   KK->zero();



    adept::Stack & stack = FemusInit::_adeptStack;  // call the adept stack object for potential use of AD

    const assemble_jacobian< real_num, double > * unk_assemble_jac/* = assemble_jacobian_base::build(dim)*/;


 //***************************************************  
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = 3;

  std::vector < std::vector < /*double*/ real_num_mov > >  JacI_qp(space_dim);
  std::vector < std::vector < /*double*/ real_num_mov > >  Jac_qp(dim);
    for (unsigned d = 0; d < dim; d++) { Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < space_dim; d++) { JacI_qp[d].resize(dim); }
    
    real_num_mov detJac_qp;
 //*************************************************** 

    //=============== Integration ========================================
    real_num_mov weight_qp;    // must be adept if the domain is moving

    //=============== Geometry - BEGIN ========================================
    unsigned xType = CONTINUOUS_BIQUADRATIC;

    CurrentElem < real_num_mov > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
            Phi < real_num_mov > geom_element_phi_dof_qp(dim_offset_grad/*dim*/);                   // must be adept if the domain is moving, otherwise double
    //=============== Geometry - END ========================================

    //=============== Unknowns - BEGIN ========================================

    const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();

    std::vector < UnknownLocal < real_num > >       unknowns_local(n_unknowns);                             //-- at dofs
    std::vector <          Phi < real_num/*real_num_mov*/ > >   unknowns_phi_dof_qp(n_unknowns, Phi< real_num/*real_num_mov*/ >(dim_offset_grad/*dim*/));   //-- at dofs and quadrature points ---------------

    for(int u = 0; u < n_unknowns; u++) {
        unknowns_local[u].initialize(dim_offset_grad/*dim*/, unknowns[u], ml_sol, mlPdeSys);
        assert(u == unknowns_local[u].pde_index());  //I would like this ivar order to coincide either with SolIndex or with SolPdeIndex, otherwise too many orders... it has to be  with SolPdeIndex, because SolIndex is related to the Solution object which can have more fields... This has to match the order of the unknowns[] argument
    }


    //===============  Unknowns, Elem matrix and Rhs - BEGIN ========================================
    ElementJacRes < real_num >       unk_element_jac_res(dim, unknowns_local);
    //===============  Unknowns, Elem matrix and Rhs - END ========================================
    
    std::vector < unsigned int > unk_num_elem_dofs_interface(n_unknowns); //to avoid recomputing offsets inside quadrature
    //=============== Unknowns - END ========================================


    //=============== Quantities that are not unknowns - BEGIN ========================================
    
     UnknownLocal < double >  sol_exact;
     sol_exact.initialize(dim_offset_grad, 
                          unknowns_local[0].name(), 
                          unknowns_local[0].fe_type(), 
                          unknowns_local[0].pde_index(), 
                          unknowns_local[0].sol_index());
     
    //=============== Quantities that are not unknowns - END ========================================
 
  
  
    // element loop: each process loops only on the elements that owns
    for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {


    //=============== Geometry - BEGIN ========================================
        geom_element.set_coords_at_dofs_and_geom_type(iel, xType);
        
        geom_element.set_elem_center_3d(iel, xType);
        
        const short unsigned ielGeom = geom_element.geom_type();
    //=============== Geometry - END ========================================

        
    //=============== Unknowns - BEGIN ========================================
        for (unsigned  u = 0; u < n_unknowns; u++) {
            unknowns_local[u].set_elem_dofs(iel, msh, sol);
        }
    //=============== Unknowns - END ========================================

    //=============== Quantities that are not unknowns - BEGIN ========================================
        //must be called AFTER the unknowns 
        sol_exact.set_elem_dofs(unknowns_local[0].num_elem_dofs(), geom_element, * exact_sol[0] );
    //=============== Quantities that are not unknowns - END ========================================

        

    //===============  Unknowns, Elem matrix and Rhs - BEGIN ========================================
        unk_element_jac_res.set_loc_to_glob_map(iel, msh, pdeSys);
    //===============  Unknowns, Elem matrix and Rhs - END ========================================


        unk_assemble_jac->prepare_before_integration_loop(stack);


//interface to avoid computation inside quadrature - BEGIN 
        unsigned sum_unk_num_elem_dofs_interface = 0;
        for (unsigned  u = 0; u < n_unknowns; u++) {
            unk_num_elem_dofs_interface[u]   = unknowns_local[u].num_elem_dofs();
            sum_unk_num_elem_dofs_interface += unk_num_elem_dofs_interface[u];
        }
//interface to avoid computation inside quadrature - END




        // *** Gauss point loop ***
        for (unsigned ig = 0; ig < quad_rules[ielGeom].GetGaussPointsNumber(); ig++) {

      elem_all/*_for_domain*/[ielGeom][xType]->JacJacInv(geom_element.get_coords_at_dofs_3d(), ig, Jac_qp, JacI_qp, detJac_qp, space_dim);
      weight_qp = detJac_qp * quad_rules[ielGeom].GetGaussWeightsPointer()[ig];
            
            // *** get gauss point weight, test function and test function partial derivatives ***
     for (unsigned  u = 0; u < n_unknowns; u++) {
         elem_all[ielGeom][unknowns_local[u].fe_type()]->shape_funcs_current_elem(ig, JacI_qp, unknowns_phi_dof_qp[u].phi(), unknowns_phi_dof_qp[u].phi_grad(), unknowns_phi_dof_qp[u].phi_hess(), space_dim);
     }

     elem_all_for_domain[ielGeom][xType]->shape_funcs_current_elem(ig, JacI_qp, geom_element_phi_dof_qp.phi(), geom_element_phi_dof_qp.phi_grad(), geom_element_phi_dof_qp.phi_hess(), space_dim);
     

            // evaluate the solution, the solution derivatives and the coordinates in the gauss point
            real_num solu_gss = 0.;
            std::vector < real_num > gradSolu_gss(dim_offset_grad, 0.);
            std::vector < real_num >   gradSolu_exact_gss(dim_offset_grad, 0.);

            for (unsigned  u = 0; u < n_unknowns; u++) {
                for (unsigned i = 0; i < unknowns_local[u].num_elem_dofs(); i++) {
                    solu_gss += unknowns_phi_dof_qp[u].phi(i) * unknowns_local[u].elem_dofs()[i];

                    for (unsigned jdim = 0; jdim < dim_offset_grad; jdim++) {
                        gradSolu_gss[jdim]       += unknowns_phi_dof_qp[u].phi_grad(i * dim_offset_grad + jdim) * unknowns_local[0].elem_dofs()[i];
                        gradSolu_exact_gss[jdim] += unknowns_phi_dof_qp[u].phi_grad(i * dim_offset_grad + jdim) * sol_exact.elem_dofs()[i];
                    }
                }
            }


            std::vector < double > x_gss(dim, 0.);
            for (unsigned i = 0; i < geom_element.get_coords_at_dofs()[0].size(); i++) {
                for (unsigned jdim = 0; jdim < x_gss.size(); jdim++) {
                    x_gss[jdim] += geom_element.get_coords_at_dofs(jdim,i) * geom_element_phi_dof_qp.phi(i);
                }
            }



            for (unsigned i = 0; i < unknowns_local[0].num_elem_dofs(); i++) {

                real_num laplace_weak = 0.;
                real_num laplace_weak_exact = 0.;

                for (unsigned jdim = 0; jdim < dim_offset_grad; jdim++) {
                    laplace_weak       +=  unknowns_phi_dof_qp[0].phi_grad(i * dim_offset_grad + jdim) * gradSolu_gss[jdim];
                    laplace_weak_exact +=  unknowns_phi_dof_qp[0].phi_grad(i * dim_offset_grad + jdim) * gradSolu_exact_gss[jdim];
                }


// // Mass(u) = 1
//                 double mass_exact = 1.;
//                 unk_element_jac_res.res()[i] += ( ( mass_exact - solu_gss ) * unknowns_phi_dof_qp[0].phi(i) ) * weight_qp;
//                 
// // Mass(u) = Mass(u_0)
//                 double mass_exact = exact_sol.value(x_gss);
//                 unk_element_jac_res.res()[i] += ( ( mass_exact - solu_gss ) * unknowns_phi_dof_qp[0].phi(i) ) * weight_qp;

// Helmholtz(u) = source : strong
//               double source_term = exact_sol.value(x_gss);
//         unk_element_jac_res.res()[i] += ( source_term * phi()[i] - solu_gss * phi()[i] - laplace_weak ) * weight_qp;

// Helmholtz(u) = Helmholtz(u_0) : strong
//                 double helmholtz_strong_exact = exact_sol.helmholtz(x_gss);
//                 unk_element_jac_res.res()[i] += (helmholtz_strong_exact * unknowns_phi_dof_qp[0].phi(i) - solu_gss * unknowns_phi_dof_qp[0].phi(i) - laplace_weak) * weight_qp;


// grad(u) grad(v)  =  - Laplace(u_0) v : strong
               double laplace_strong_exact = exact_sol[0]->laplacian(x_gss);
        unk_element_jac_res.res()[i] += ( - laplace_strong_exact * unknowns_phi_dof_qp[0].phi(i) - laplace_weak) * weight_qp;        //strong form of RHS and weak form of LHS

// // grad(u) grad(v)  =   grad(u_0) grad(v) : weak
//            unk_element_jac_res.res()[i] += (  laplace_weak_exact  - laplace_weak) * weight_qp;                  //weak form of RHS and weak form of LHS


// // Laplace(u) = source 
//               double source_term = 1.;
// //               double source_term;
// //               if (geom_element.get_elem_center_3d()[0] < 0.5 + 1.e-5) source_term = 1.;
// //               else source_term = 2.;
//         unk_element_jac_res.res()[i] += (  source_term * unknowns_phi_dof_qp[0].phi(i) - laplace_weak) * weight_qp;        //strong form of RHS and weak form of LHS




                unk_assemble_jac->compute_jacobian_inside_integration_loop(i, dim_offset_grad, unk_num_elem_dofs_interface, sum_unk_num_elem_dofs_interface, unknowns_local, unknowns_phi_dof_qp, weight_qp, unk_element_jac_res.jac());



            } // end phi_i loop


        } // end gauss point loop


        unk_assemble_jac->compute_jacobian_outside_integration_loop(stack, unknowns_local, unk_element_jac_res.res(), unk_element_jac_res.jac(), unk_element_jac_res.dof_map(), RES, KK);


    } //end element loop for each process


    RES->close();
    KK->close();

    // ***************** END ASSEMBLY *******************
}


namespace femus {
  
// here, the peculiarity of the application is not dealt with using virtuality, but with template specialization ---------------------------------

// template specialization for double
template < >
void assemble_jacobian< double, double >::prepare_before_integration_loop(adept::Stack& stack) const { }


// template specialization for double
template < >
void  assemble_jacobian < double, double > ::compute_jacobian_outside_integration_loop (
    adept::Stack & stack,
    const std::vector< std::vector< double > > & solu,
    const std::vector< double > & Res,
    std::vector< double > & Jac,
    const std::vector< int > & loc_to_glob_map_all_vars,
    NumericVector*           RES,
    SparseMatrix*             KK
)  const {

    RES->add_vector_blocked(Res, loc_to_glob_map_all_vars);
    KK->add_matrix_blocked(Jac, loc_to_glob_map_all_vars, loc_to_glob_map_all_vars);

}


// template specialization for double
template < >
void  assemble_jacobian < double, double > ::compute_jacobian_outside_integration_loop (
    adept::Stack & stack,
    const std::vector< UnknownLocal < double > > & unk_vec,
    const std::vector< double > & Res,
    std::vector< double > & Jac,
    const std::vector< int > & loc_to_glob_map_all_vars,
    NumericVector*           RES,
    SparseMatrix*             KK
)  const {

    RES->add_vector_blocked(Res, loc_to_glob_map_all_vars);
    KK->add_matrix_blocked(Jac, loc_to_glob_map_all_vars, loc_to_glob_map_all_vars);

}


// template specialization for double
// this is where you fill the jacobian in traditional (much faster) way
template < >
void  assemble_jacobian< double, double >::compute_jacobian_inside_integration_loop (
    const unsigned i,
    const unsigned dim,
    const std::vector < unsigned int > Sol_n_el_dofs,
    const unsigned int sum_Sol_n_el_dofs,
    const std::vector< UnknownLocal < double > > & unk_vec,
    const std::vector< Phi <double> > &  phi,
    const double weight_qp,
    std::vector< double > & Jac )  const {

        
  const unsigned int dim_offset_grad = 3;

  
    constexpr unsigned int pos_unk = 0/*unk_vec[0].SolPdeIndex*/;

// *** phi_j loop ***
    for (unsigned j = 0; j < Sol_n_el_dofs[ pos_unk ]; j++) {
        /*real_num*/double laplace_weak_jac = 0.;

        for (unsigned kdim = 0; kdim < dim; kdim++) {
            laplace_weak_jac += (phi[ pos_unk].phi_grad(i * dim_offset_grad + kdim) * phi[ pos_unk ].phi_grad(j * dim_offset_grad + kdim));
        }

        Jac[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_unk, pos_unk, i, j) ] += (laplace_weak_jac /*+   phi[ pos_unk ].phi(i) * phi[ pos_unk ].phi(j) */ ) * weight_qp;
    } // end phi_j loop


}


} //end namespace

#endif
