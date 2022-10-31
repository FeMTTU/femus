#ifndef __femus_applications_00_tutorial_hpp__
#define __femus_applications_00_tutorial_hpp__

#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"
#include "CurrentElem.hpp"

using namespace femus;


 
const std::vector< Unknown >  systems__provide_list_of_unknowns_lagrangian() {


    std::vector< FEFamily >     feFamily = {LAGRANGE, LAGRANGE, LAGRANGE};
    std::vector< FEOrder >       feOrder = {FIRST, SERENDIPITY, SECOND};
    std::vector< int >        time_order = {0, 0, 0};  //0 = steady, 2 = time-dependent
    std::vector< bool >   is_pde_unknown = {true, true, true};

    assert( feFamily.size() == feOrder.size());
    assert( feFamily.size() == is_pde_unknown.size());
    assert( feFamily.size() == time_order.size());

    std::vector< Unknown >  unknowns(feFamily.size());

    for (unsigned int fe = 0; fe < unknowns.size(); fe++) {

        std::ostringstream unk;
        unk << "u" << "_" << feFamily[fe] << "_" << feOrder[fe];
        unknowns[fe]._name           = unk.str();
        unknowns[fe]._fe_family      = feFamily[fe];
        unknowns[fe]._fe_order       = feOrder[fe];
        unknowns[fe]._time_order     = time_order[fe];
        unknowns[fe]._is_pde_unknown = is_pde_unknown[fe];

    }


    return unknowns;

}



const std::vector< Unknown >  systems__provide_list_of_unknowns_all_fe() {


    std::vector< FEFamily >     feFamily = {LAGRANGE, LAGRANGE, LAGRANGE, DISCONTINUOUS_POLYNOMIAL, DISCONTINUOUS_POLYNOMIAL};
    std::vector< FEOrder >       feOrder = {FIRST, SERENDIPITY, SECOND, ZERO, FIRST};
    std::vector< int >        time_order = {0, 0, 0, 0, 0};  //0 = steady, 2 = time-dependent
    std::vector< bool >   is_pde_unknown = {true, true, true, true, true};

    assert( feFamily.size() == feOrder.size());
    assert( feFamily.size() == is_pde_unknown.size());
    assert( feFamily.size() == time_order.size());

    std::vector< Unknown >  unknowns(feFamily.size());

    for (unsigned int fe = 0; fe < unknowns.size(); fe++) {

        std::ostringstream unk;
        unk << "u" << "_" << feFamily[fe] << "_" << feOrder[fe];
        unknowns[fe]._name           = unk.str();
        unknowns[fe]._fe_family      = feFamily[fe];
        unknowns[fe]._fe_order       = feOrder[fe];
        unknowns[fe]._time_order     = time_order[fe];
        unknowns[fe]._is_pde_unknown = is_pde_unknown[fe];

    }


    return unknowns;

}




template < class type = double >
class Zero : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return 0.;
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  =  0.;
        solGrad[1]  =  0.;

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return  0.; 
    }

      
};





namespace  Domain_square_01by01  {
    


template < class type = double >
class Function_NonZero_on_boundary_1 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return x[0] * x[0] * (1. - x[0]) + sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  =  - x[0] * x[0] +  (1. - x[0]) * 2. * x[0]  + pi * cos( pi * (x[0]) ) * sin( pi * (x[1]) );
        solGrad[1]  =                                              pi * sin( pi * (x[0]) ) * cos( pi * (x[1]) );

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return   - 2. * x[0] +  2. * (1. - x[0])  -  2. * x[0]  - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) ) - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return + sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  = pi * cos( pi * (x[0]) ) * sin( pi * (x[1]) );
        solGrad[1]  = pi * sin( pi * (x[0]) ) * cos( pi * (x[1]) );

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return -pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) ) - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};


//this solution shows SUPERCONVERGENCE for SERENDIPITY FE, and it is like SUPER PERFECT for BIQUADRATIC FE... it is because of the MESH!
template < class type = double >
class Function_Zero_on_boundary_2 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return  x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  = (1. - 2. * x[0]) *  x[1] * (1. - x[1]);
        solGrad[1]  = (1. - 2. * x[1]) *  x[0] * (1. - x[0]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return -2. * ( x[0] * (1. - x[0])  + x[1] * (1. - x[1]) );
    }



};


//this solution does not have SUPERCONVERGENCE even with the straight mesh
template < class type = double >
class Function_Zero_on_boundary_3 : public Math::Function< type > {

public:

// manufactured Laplacian =============
    type value(const std::vector < type >& x) const {
        
        return    x[0] *  x[0] * (1. - x[0]) * x[1] * (1. - x[1]) ;
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  =  x[0] * (1. - 2. * x[0]) *  x[1] * (1. - x[1]) +  x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
        solGrad[1]  =  x[0] * (1. - 2. * x[1]) *  x[0] * (1. - x[0]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return     x[0] *  (  -2. *  x[1] * (1. - x[1])   ) + 2. *  (1. - 2. * x[0]) *  x[1] * (1. - x[1])   +  x[0] * ( -2. *  x[0] * (1. - x[0])) ;
    }



};



} //end namespace


namespace  Domain_square_m05p05  {

 template < class type = double >
class Function_Zero_on_boundary_4 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
  return cos(pi * x[0]) * cos(pi * x[1]);
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
        solGrad[1]  = -pi * cos(pi * x[0]) * sin(pi * x[1]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return -pi * pi * cos(pi * x[0]) * cos(pi * x[1]) - pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};

   
    

std::string quad9_all_mesh_generation_methods_Structured(const unsigned method_flag, MultiLevelMesh & ml_mesh, const std::string fe_quad_rule) {
    

    
  double scalingFactor = 1.;

    std::string mesh_name("square");
    
    switch(method_flag) {
        case 0: {
            mesh_name += "_femus";
  // ======= Mesh, coarse gen, I: from function - BEGIN  ========================
  const unsigned int nsub_x = 2;
  const unsigned int nsub_y = 2;
  const unsigned int nsub_z = 0;
  const std::vector<double> xyz_min = {-0.5,-0.5,0.};
  const std::vector<double> xyz_max = { 0.5, 0.5,0.};
  const ElemType geom_elem_type = QUAD9/*TRI6*/;
  ml_mesh.GenerateCoarseBoxMesh(nsub_x, nsub_y, nsub_z, xyz_min[0], xyz_max[0], xyz_min[1], xyz_max[1], xyz_min[2], xyz_max[2], geom_elem_type, fe_quad_rule.c_str() );
  // ======= Mesh, coarse gen, I: from function - END ========================
 
            break;
        }
        
        case 1: {
            mesh_name += "_salome";
//   // ======= Mesh, coarse gen, III: from Salome - BEGIN  ========================
   ml_mesh.ReadCoarseMesh("./input/square_2x2_centered_at_origin.med", fe_quad_rule.c_str(), scalingFactor);
//   // ======= Mesh, coarse gen, III: from Salome - END ========================
            break;
        }
        case 2: {
            mesh_name += "_gambit";
//   // ======= Mesh, coarse gen, II: from Gambit - BEGIN  ========================
    ml_mesh.ReadCoarseMesh("./input/square_quad.neu", fe_quad_rule.c_str(), scalingFactor);
//     ml_mesh.ReadCoarseMesh("./input/square_tri.neu", fe_quad_rule.c_str(), scalingFactor);
//   //   ml_mesh.ReadCoarseMesh("./input/cube_tet.neu", fe_quad_rule.c_str(), scalingFactor);
  // ======= Mesh, coarse gen, II: from Gambit - END ========================
            break;
        }
        
        default: { abort(); }
    }
    
    
    return mesh_name;
    
}

    
}



namespace  Domain_L_shaped  {
    
  template < class type = double >
class Function_NonZero_on_boundary_2 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return x[0] * x[0] * (1. - x[0]) + sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  =  - x[0] * x[0] +  (1. - x[0]) * 2. * x[0]  + pi * cos( pi * (x[0]) ) * sin( pi * (x[1]) );
        solGrad[1]  =                                              pi * sin( pi * (x[0]) ) * cos( pi * (x[1]) );

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return   - 2. * x[0] +  2. * (1. - x[0])  -  2. * x[0]  - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) ) - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};
  
    
}



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
void System_assemble_flexible_Laplacian_With_Manufactured_Sol(const std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > > & elem_all,
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
    
    double detJac_qp;
 //*************************************************** 

    //=============== Integration ========================================
    real_num_mov weight_qp;    // must be adept if the domain is moving

    //=============== Geometry - BEGIN ========================================
    unsigned xType = BIQUADR_FE;

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
    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {


    //=============== Geometry - BEGIN ========================================
        geom_element.set_coords_at_dofs_and_geom_type(iel, xType);
        
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

      elem_all[ielGeom][xType]->JacJacInv(geom_element.get_coords_at_dofs_3d(), ig, Jac_qp, JacI_qp, detJac_qp, space_dim);
      weight_qp = detJac_qp * quad_rules[ielGeom].GetGaussWeightsPointer()[ig];
            
            // *** get gauss point weight, test function and test function partial derivatives ***
     for (unsigned  u = 0; u < n_unknowns; u++) {
         elem_all[ielGeom][unknowns_local[u].fe_type()]->shape_funcs_current_elem(ig, JacI_qp, unknowns_phi_dof_qp[u].phi(), unknowns_phi_dof_qp[u].phi_grad(), unknowns_phi_dof_qp[u].phi_hess(), space_dim);
     }

     elem_all[ielGeom][xType]->shape_funcs_current_elem(ig, JacI_qp, geom_element_phi_dof_qp.phi(), geom_element_phi_dof_qp.phi_grad(), geom_element_phi_dof_qp.phi_hess(), space_dim);
     

            // evaluate the solution, the solution derivatives and the coordinates in the gauss point
            real_num solu_gss = 0.;
            std::vector < real_num > gradSolu_gss(dim_offset_grad, 0.);
            std::vector < double >   gradSolu_exact_gss(dim_offset_grad, 0.);

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

                real_num laplace = 0.;
                real_num laplace_weak_exact = 0.;

                for (unsigned jdim = 0; jdim < dim_offset_grad; jdim++) {
                    laplace            +=  unknowns_phi_dof_qp[0].phi_grad(i * dim_offset_grad + jdim) * gradSolu_gss[jdim];
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
//         unk_element_jac_res.res()[i] += ( source_term * phi()[i] - solu_gss * phi()[i] - laplace ) * weight_qp;

// Helmholtz(u) = Helmholtz(u_0) : strong
//                 double helmholtz_strong_exact = exact_sol.helmholtz(x_gss);
//                 unk_element_jac_res.res()[i] += (helmholtz_strong_exact * unknowns_phi_dof_qp[0].phi(i) - solu_gss * unknowns_phi_dof_qp[0].phi(i) - laplace) * weight_qp;

// Laplace(u) = Laplace(u_0) : strong
//                double laplace_strong_exact = exact_sol[0]->laplacian(x_gss);
//         unk_element_jac_res.res()[i] += (- laplace_strong_exact * unknowns_phi_dof_qp[0].phi(i) - laplace) * weight_qp;        //strong form of RHS and weak form of LHS

// Laplace(u) = source 
              double source_term = 1.;
        unk_element_jac_res.res()[i] += (  source_term * unknowns_phi_dof_qp[0].phi(i) - laplace) * weight_qp;        //strong form of RHS and weak form of LHS


// grad(u) grad(v) = grad(u_0) grad(v) : weak
//            unk_element_jac_res.res()[i] += (laplace_weak_exact  - laplace) * weight_qp;                  //weak form of RHS and weak form of LHS



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
        /*real_num*/double laplace_jac = 0.;

        for (unsigned kdim = 0; kdim < dim; kdim++) {
            laplace_jac += (phi[ pos_unk].phi_grad(i * dim_offset_grad + kdim) * phi[ pos_unk ].phi_grad(j * dim_offset_grad + kdim));
        }

        Jac[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_unk, pos_unk, i, j) ] += (laplace_jac /*+   phi[ pos_unk ].phi(i) * phi[ pos_unk ].phi(j) */ ) * weight_qp;
    } // end phi_j loop


}


} //end namespace








#endif
