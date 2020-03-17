//compute convergence using Cauchy convergence test
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv
//

/** tutorial/Ex2
 * This example shows how to set and solve the weak form of the Poisson problem
 *                    $$ \Delta u = 1 \text{ on }\Omega, $$
 *          $$ u=0 \text{ on } \Gamma, $$
 * on a square domain $\Omega$ with boundary $\Gamma$;
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "LinearImplicitSystem.hpp"
#include "Files.hpp"
#include "FE_convergence.hpp"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"
#include "CurrentElem.hpp"

#include "adept.h"

// command to view matrices
// ./tutorial_ex2_c -mat_view ::ascii_info_detail
// ./tutorial_ex2_c -mat_view > matview_print_in_file.txt

using namespace femus;





template < class type = double >
class My_exact_solution : public Math::Function< type > {

public:

// manufactured Laplacian =============
    type value(const std::vector < type >& x) const {
        
        return sin( pi * (x[0]) ) * sin( pi * (x[1]) );
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



// constant =============
//   double value(const std::vector < double >& x) const {  return 1.; }
//
//
//  vector < double >  gradient(const std::vector < double >& x) const {
//
//     vector < double > solGrad(x.size());
//
//    for (int d = 0; d < x.size(); d++)   solGrad[d]  = 0.;
//
//   return solGrad;
// }
//
//
//  double laplacian(const std::vector < double >& x) const {  return 0.; }

  private: 
    
   static constexpr double pi = acos(-1.);
      
};


double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char * name) {

            double pi = acos(-1.);
   double value = 0.; /*sin( pi * (x[0]) ) * sin( pi * (x[1]) );*/ /*(x[0]) * (1. - x[0]) * (x[1]) * (1. - x[1]);*//*observe the interpolation convergence of this last one...*/

   return value;   

}


bool Solution_set_boundary_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {

    bool dirichlet = true; //dirichlet
    value = 0.;

    return dirichlet;
}


template < class system_type, class real_num, class real_num_mov >
void System_assemble_interface(MultiLevelProblem & ml_prob);

template < class system_type, class real_num, class real_num_mov >
void System_assemble_flexible(const std::vector < std::vector < const elem_type_templ_base<real_num, real_num_mov> *  > > & elem_all,
                              const std::vector<Gauss> & quad_rules,
                              MultiLevelMesh * ml_mesh,
                              MultiLevelSolution * ml_sol,
                              system_type * mlPdeSys,
                              const std::vector< Unknown > &  unknowns,
                              const Math::Function< double > & exact_sol);


//Unknown definition  ==================
const std::vector< Unknown >  provide_list_of_unknowns() {


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



// this is the class that I pass to the FE_convergence functions. Basically I pass a class, instead of passing a function pointer.
// notice that this is a class template, not a class with a function template
template < class real_num >
class My_main_single_level : public Main_single_level {

public:

    const MultiLevelSolution  run_on_single_level(const Files & files,
                                                  MultiLevelProblem & ml_prob,
                                                  const std::vector< Unknown > & unknowns,
                                                  const MultiLevelSolution::BoundaryFuncMLProb  SetBoundaryCondition_in,
                                                  const MultiLevelSolution::InitFuncMLProb      SetInitialCondition_in,
                                                  MultiLevelMesh & ml_mesh,
                                                  const unsigned i) const;

};






int main(int argc, char** args) {

    // ======= Init ==========================
    FemusInit mpinit(argc, args, MPI_COMM_WORLD);

    // ======= Files =========================
    Files files;
    files.CheckIODirectories();
    files.RedirectCout();

    // ======= Quad Rule ========================
    std::string quad_rule_order("seventh");
    /* this is the order of accuracy that is used in the gauss integration scheme
       In the future it is not going to be an argument of the mesh function   */
    
    // ======= Problem ========================
    MultiLevelProblem ml_prob;
    ml_prob.SetFilesHandler(&files);
    ml_prob.SetQuadratureRuleAllGeomElems(quad_rule_order);
    ml_prob.set_all_abstract_fe();


    // ======= Mesh ========================
    const ElemType geom_elem_type = QUAD9;
    const std::vector< unsigned int > nsub = {2,2,0};
    const std::vector< double >      xyz_min = {0.,0.,0.};
    const std::vector< double >      xyz_max = {1.,1.,0.};


    MultiLevelMesh ml_mesh;
    ml_mesh.GenerateCoarseBoxMesh(nsub[0], nsub[1], nsub[2], xyz_min[0], xyz_max[0], xyz_min[1], xyz_max[1], xyz_min[2], xyz_max[2], geom_elem_type, quad_rule_order.c_str());
//    std::string input_file = "Lshape_4.med";
//    std::string input_file = "Lshape.med";
//     std::string input_file = "interval.med";
//     std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
//     const std::string infile = mystream.str();
//   ml_mesh.ReadCoarseMesh(infile.c_str(),quad_rule_order.c_str(),1.);


    // ======= Unknowns ========================
    std::vector< Unknown > unknowns = provide_list_of_unknowns();



    // ======= Normal run (without convergence study) ========================
    My_main_single_level< /*adept::a*/double > my_main;
//     const unsigned int n_levels = 3;
//     my_main.run_on_single_level(files, ml_prob, unknowns, Solution_set_boundary_conditions, Solution_set_initial_conditions, ml_mesh, n_levels);

    // ======= Convergence study ========================

    // set total number of levels ================
    unsigned max_number_of_meshes = 6;

    if (ml_mesh.GetDimension() == 3) max_number_of_meshes = 4;

    //set coarse storage mesh (///@todo should write the copy constructor or "=" operator to copy the previous mesh) ==================
    MultiLevelMesh ml_mesh_all_levels;
    ml_mesh_all_levels.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,quad_rule_order.c_str());
//    ml_mesh_all_levels.ReadCoarseMesh(infile.c_str(),quad_rule_order.c_str(),1.);


    // convergence choices ================
    My_exact_solution<> exact_sol;         //provide exact solution, if available ==============
    const unsigned conv_order_flag = 0;    //Choose how to compute the convergence order ============== //0: incremental 1: absolute (with analytical sol)  2: absolute (with projection of finest sol)...
    const unsigned norm_flag = 1;          //Choose what norms to compute (//0 = only L2: //1 = L2 + H1) ==============


    // object ================
    FE_convergence<>  fe_convergence;

    fe_convergence.convergence_study(files, ml_prob, unknowns, Solution_set_boundary_conditions, Solution_set_initial_conditions, ml_mesh, ml_mesh_all_levels, max_number_of_meshes, norm_flag, conv_order_flag, my_main/*, & exact_sol*/);


    return 0;

}









template < class real_num >
const MultiLevelSolution  My_main_single_level< real_num >::run_on_single_level(const Files & files,
                                                                                MultiLevelProblem & ml_prob,
                                                                                const std::vector< Unknown > &  unknowns,
                                                                                const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition_in,
                                                                                const MultiLevelSolution::InitFuncMLProb SetInitialCondition_in,
                                                                                MultiLevelMesh & ml_mesh,
                                                                                const unsigned lev) const {


    //Mesh  ==================
    unsigned numberOfUniformLevels = lev + 1;
    unsigned numberOfSelectiveLevels = 0;
    ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
    ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);

    ml_mesh.PrintInfo();


    //Solution  ==================
    MultiLevelSolution ml_sol_single_level(&ml_mesh);

    ml_sol_single_level.SetWriter(VTK);
    ml_sol_single_level.GetWriter()->SetDebugOutput(true);

    // ======= Problem ========================
    ml_prob.SetMultiLevelMeshAndSolution(& ml_mesh,& ml_sol_single_level);
    
    ml_prob.get_systems_map().clear();  //at every lev we'll have a different map of systems

    //only one Problem,
    //only one Mesh per level,
    //only one Solution per level,
    // a separate System for each unknown (I want to do like this to see how to handle multiple coupled systems later on)

    for (unsigned int u = 0; u < unknowns.size(); u++) {

        // ======= Solution, II ==================
        ml_sol_single_level.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);
        ml_sol_single_level.Initialize(unknowns[u]._name.c_str(), SetInitialCondition_in, & ml_prob);
        ml_sol_single_level.AttachSetBoundaryConditionFunction(SetBoundaryCondition_in);
        ml_sol_single_level.GenerateBdc(unknowns[u]._name.c_str(), "Steady", & ml_prob);

// If you just want an interpolation study, without equation, just initialize every Solution with some function and comment out all the following System part        
        // ======= System ========================
        std::ostringstream sys_name;
        sys_name << unknowns[u]._name;  //give to each system the name of the unknown it solves for!

        LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > (sys_name.str());

        system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
        std::vector< Unknown > unknowns_vec(1);
        unknowns_vec[0] = unknowns[u]; //need to turn this into a vector
        system.set_unknown_list_for_assembly(unknowns_vec); //way to communicate to the assemble function, which doesn't belong to any class
        ml_prob.set_current_system_number(u);               //way to communicate to the assemble function, which doesn't belong to any class

        system.SetAssembleFunction(System_assemble_interface< LinearImplicitSystem, real_num, double >);

        // initialize and solve the system
        system.init();
        system.ClearVariablesToBeSolved();
        system.AddVariableToBeSolved("All");

//             system.SetDebugLinear(true);
//             system.SetMaxNumberOfLinearIterations(6);
//             system.SetAbsoluteLinearConvergenceTolerance(1.e-4);

        system.SetOuterSolver(GMRES);
        system.MGsolve();  //everything is stored into the Solution after this

        // ======= Print ========================
        std::vector < std::string > variablesToBePrinted;
        variablesToBePrinted.push_back(unknowns[u]._name);
        ml_sol_single_level.GetWriter()->Write(unknowns[u]._name, files.GetOutputPath(), "biquadratic", variablesToBePrinted, lev);

    }


    return ml_sol_single_level;
}



template < class system_type, class real_num, class real_num_mov >
void System_assemble_interface(MultiLevelProblem& ml_prob) {
    //  ml_prob is the global object from/to where get/set all the data

    // this is meant to be like a tiny addition to the main function, because we cannot pass these arguments through the function pointer
//what is funny is that this function is attached to a system, but I cannot retrieve the system to which it is attached unless I know the name or the number of it!
//all I have at hand is the MultiLevelProblem, which contains a Vector of Systems
// all I can do is put in the MultiLevelProblem a number that tells me what is the current system being solved


    My_exact_solution< double > exact_sol;  //this one I reproduce it here, otherwise I should pass it in the main to the MultiLevelProblem

    const unsigned current_system_number = ml_prob.get_current_system_number();

   //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);

    System_assemble_flexible< system_type, real_num, real_num_mov > (elem_all,
                                                                     ml_prob.GetQuadratureRuleAllGeomElems(),
                                                                     ml_prob._ml_msh,
                                                                     ml_prob._ml_sol,
                                                                     & ml_prob.get_system< system_type >(current_system_number),
                                                                     ml_prob.get_system< system_type >(current_system_number).get_unknown_list_for_assembly(),
                                                                     exact_sol);

}


/**
 * This function assemble the stiffness matrix KK and the residual vector Res
 * for the Newton iterative scheme
 *                  J(u0) w = Res(u_0) = f(x) - J u_0  ,
 *                  with u = u0 + w
 *                  J = \grad_u F
 **/

template < class system_type, class real_num, class real_num_mov >
void System_assemble_flexible(const std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > > & elem_all,
                              const std::vector<Gauss> & quad_rules,
                              MultiLevelMesh * ml_mesh_in,
                              MultiLevelSolution * ml_sol_in,
                              system_type * mlPdeSys,
                              const std::vector< Unknown > &  unknowns,
                              const Math::Function< double > & exact_sol) {

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

  std::vector < std::vector < double > >  JacI_qp(space_dim);
  std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < dim; d++) { Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < space_dim; d++) { JacI_qp[d].resize(dim); }
    
    double detJac_qp;
 //*************************************************** 

    //=============== Integration ========================================
    real_num_mov weight_qp;    // must be adept if the domain is moving

    //=============== Geometry ========================================
    unsigned xType = BIQUADR_FE;

    CurrentElem < real_num_mov > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
            Phi < real_num_mov > geom_element_phi_dof_qp(dim_offset_grad/*dim*/);                   // must be adept if the domain is moving, otherwise double

    //=============== Unknowns ========================================

    const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();

    std::vector < UnknownLocal < real_num > >       unknowns_local(n_unknowns);                             //-- at dofs
    std::vector <          Phi < real_num_mov > >   unknowns_phi_dof_qp(n_unknowns, Phi< real_num_mov >(dim_offset_grad/*dim*/));   //-- at dofs and quadrature points ---------------

    for(int u = 0; u < n_unknowns; u++) {
        unknowns_local[u].initialize(dim_offset_grad/*dim*/, unknowns[u], ml_sol, mlPdeSys);
        assert(u == unknowns_local[u].pde_index());  //I would like this ivar order to coincide either with SolIndex or with SolPdeIndex, otherwise too many orders... it has to be  with SolPdeIndex, because SolIndex is related to the Solution object which can have more fields... This has to match the order of the unknowns[] argument
    }


    //=============== Elem matrix ========================================
    ElementJacRes < real_num >       unk_element_jac_res(dim, unknowns_local);
    std::vector < unsigned int > unk_num_elem_dofs_interface(n_unknowns); //to avoid recomputing offsets inside quadrature


    //=============== Quantities that are not unknowns ========================================
    
     UnknownLocal < double >  sol_exact;
     sol_exact.initialize(dim_offset_grad, 
                          unknowns_local[0].name(), 
                          unknowns_local[0].fe_type(), 
                          unknowns_local[0].pde_index(), 
                          unknowns_local[0].sol_index());
     
 
  
  
    // element loop: each process loops only on the elements that owns
    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {


        geom_element.set_coords_at_dofs_and_geom_type(iel, xType);
        
        const short unsigned ielGeom = geom_element.geom_type();

        
        for (unsigned  u = 0; u < n_unknowns; u++) {
            unknowns_local[u].set_elem_dofs(iel, msh, sol);
        }

        //must be called AFTER the unknowns 
        sol_exact.set_elem_dofs(unknowns_local[0].num_elem_dofs(), geom_element, exact_sol);

        

        unk_element_jac_res.set_loc_to_glob_map(iel, msh, pdeSys);


        unk_assemble_jac->prepare_before_integration_loop(stack);


//interface to avoid computation inside quadrature
        unsigned sum_unk_num_elem_dofs_interface = 0;
        for (unsigned  u = 0; u < n_unknowns; u++) {
            unk_num_elem_dofs_interface[u]   = unknowns_local[u].num_elem_dofs();
            sum_unk_num_elem_dofs_interface += unk_num_elem_dofs_interface[u];
        }
//interface to avoid computation inside quadrature



        if (dim != 2) abort(); //only implemented in 2D now

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

// Helmholtz(u) = source - strong
//               double source_term = exact_sol.value(x_gss);
//         unk_element_jac_res.res()[i] += ( source_term * phi()[i] - solu_gss * phi()[i] - laplace ) * weight_qp;

// Helmholtz(u) = Helmholtz(u_0) - strong
                double helmholtz_strong_exact = exact_sol.helmholtz(x_gss);
                unk_element_jac_res.res()[i] += (helmholtz_strong_exact * unknowns_phi_dof_qp[0].phi(i) - solu_gss * unknowns_phi_dof_qp[0].phi(i) - laplace) * weight_qp;

// Laplace(u) = Laplace(u_0) - strong
//                double laplace_strong_exact = exact_sol.laplacian(x_gss);
//         unk_element_jac_res.res()[i] += (- laplace_strong_exact * phi()[i] - phi[i] * solu_gss - laplace) * weight_qp;        //strong form of RHS and weak form of LHS

// grad(u) grad(v) = grad(u_0) grad(v) - weak
//            unk_element_jac_res.res()[i] += (laplace_weak_exact - phi()[i] * solu_gss - laplace) * weight_qp;                  //weak form of RHS and weak form of LHS



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


    constexpr unsigned int pos_unk = 0/*unk_vec[0].SolPdeIndex*/;

// *** phi_j loop ***
    for (unsigned j = 0; j < Sol_n_el_dofs[ pos_unk ]; j++) {
        /*real_num*/double laplace_jac = 0.;

        for (unsigned kdim = 0; kdim < dim; kdim++) {
            laplace_jac += (phi[ pos_unk].phi_grad(i * dim + kdim) * phi[ pos_unk ].phi_grad(j * dim + kdim));
        }

        Jac[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_unk, pos_unk, i, j) ] += (laplace_jac + phi[ pos_unk ].phi(i) * phi[ pos_unk ].phi(j)) * weight_qp;
    } // end phi_j loop


}


} //end namespace

//-------------------------------------------------------------





///@todo: check bc for discontinuous FE
///@todo: compute error in L-\infty norm
///@todo: compute nonlinear convergence rate
///@todo: compute time convergence rate, pointwise and then in norms
///@todo: uncouple Gauss from Mesh
///@todo: make non-isoparametric Jacobian routines (abstract Jacobian) for all dims
///@todo: check solver and prec choices
///@todo: create a copy constructor/operator=() for the Mesh
///@todo: check FE convergence in 3D for various operators
///@todo: check face names in 3D
///@todo: test with linear compressible solid first
///@todo: check FE convergence with EXACT solution
///@todo: problems in the P0 or P1 should happen because of the boundary conditions not being enforced
