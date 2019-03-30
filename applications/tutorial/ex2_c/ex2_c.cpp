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
  double pi = acos(-1.);
  return cos(pi * x[0]) * cos(pi * x[1]);
}


 vector < type >  gradient(const std::vector < type >& x) const {
    
    vector < type > solGrad(x.size());
    
  double pi = acos(-1.);
  solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
  solGrad[1]  = -pi * cos(pi * x[0]) * sin(pi * x[1]);
  
  return solGrad;
}


 type laplacian(const std::vector < type >& x) const {
  double pi = acos(-1.);
  return -pi * pi * cos(pi * x[0]) * cos(pi * x[1]) - pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
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
 
 
  };

  


bool Solution_set_boundary_conditions(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
    
  bool dirichlet = true; //dirichlet
  value = 0;

  return dirichlet;
}


template < class system_type, class real_num, class real_num_mov = double >
void System_assemble_interface(MultiLevelProblem & ml_prob);

template < class system_type, class real_num, class real_num_mov = double >
void System_assemble_flexible(MultiLevelProblem& ml_prob, 
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
         
            std::ostringstream unk; unk << "u" << "_" << feFamily[fe] << "_" << feOrder[fe];
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
                                                   const std::vector< Unknown > & unknowns,  
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
  std::string fe_quad_rule("seventh");
 /* this is the order of accuracy that is used in the gauss integration scheme
    In the future it is not going to be an argument of the mesh function   */


   // ======= Mesh ========================
  const ElemType geom_elem_type = QUAD9;
  const std::vector< unsigned int > nsub = {2,2,0};
  const std::vector< double >      xyz_min = {-0.5,-0.5,0.};
  const std::vector< double >      xyz_max = { 0.5, 0.5,0.};


   MultiLevelMesh ml_mesh;
   ml_mesh.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());
//    std::string input_file = "Lshape_4.med";
//    std::string input_file = "Lshape.med";
//     std::string input_file = "interval.med";
//     std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
//     const std::string infile = mystream.str();
//   ml_mesh.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),1.);


   // ======= Unknowns ========================
   std::vector< Unknown > unknowns = provide_list_of_unknowns();
   


  // ======= Normal run ========================
    My_main_single_level< /*adept::a*/double > my_main;
//     const unsigned int n_levels = 3;
//      my_main.run_on_single_level(files, unknowns, ml_mesh, n_levels); //if you don't want the convergence study
    
   // ======= Convergence study ========================
    
   // set total number of levels ================  
   unsigned max_number_of_meshes = 6;
   
   if (ml_mesh.GetDimension() == 3) max_number_of_meshes = 4;
  
   //set coarse storage mesh (///@todo should write the copy constructor or "=" operator to copy the previous mesh) ==================
   MultiLevelMesh ml_mesh_all_levels;
   ml_mesh_all_levels.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());
//    ml_mesh_all_levels.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),1.);

 
   // convergence choices ================  
   My_exact_solution<> exact_sol;         //provide exact solution, if available ==============
   const unsigned conv_order_flag = 0;    //Choose how to compute the convergence order ============== //0: incremental 1: absolute (with analytical sol)  2: absolute (with projection of finest sol)...
   const unsigned norm_flag = 1;          //Choose what norms to compute (//0 = only L2: //1 = L2 + H1) ==============

   
   // object ================  
    FE_convergence<>  fe_convergence;
    
    fe_convergence.convergence_study(files, unknowns, Solution_set_boundary_conditions, ml_mesh, ml_mesh_all_levels, max_number_of_meshes, norm_flag, conv_order_flag, my_main);
         

  return 0;
  
}









template < class real_num > 
const MultiLevelSolution  My_main_single_level< real_num >::run_on_single_level(const Files & files,
                                                                                const std::vector< Unknown > &  unknowns,  
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
            MultiLevelProblem ml_prob(&ml_sol_single_level);
            
            ml_prob.SetFilesHandler(&files);
            
           //only one Mesh
           //only one Solution, 
           //only one Problem,
           // a separate System for each unknown (I want to do like this to see how to handle multiple coupled systems later on)
            
        for (unsigned int u = 0; u < unknowns.size(); u++) {
            
            // ======= Solution, II ==================
            ml_sol_single_level.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);
            ml_sol_single_level.Initialize(unknowns[u]._name.c_str());
            ml_sol_single_level.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
            ml_sol_single_level.GenerateBdc(unknowns[u]._name.c_str());

           // ======= System ========================
           std::ostringstream sys_name; sys_name << unknowns[u]._name;  //give to each system the name of the unknown it solves for!
           
            LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > (sys_name.str());

            system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
            std::vector< Unknown > unknowns_vec(1); unknowns_vec[0] = unknowns[u]; //need to turn this into a vector
            system.set_unknown_list_for_assembly(unknowns_vec); //way to communicate to the assemble function, which doesn't belong to any class
            ml_prob.set_current_system_number(u);               //way to communicate to the assemble function, which doesn't belong to any class
            
            system.SetAssembleFunction(System_assemble_interface< LinearImplicitSystem, real_num >);

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



template < class system_type, class real_num, class real_num_mov = double >
void System_assemble_interface(MultiLevelProblem& ml_prob) {
// this is meant to be like a tiny addition to the main function, because we cannot pass these arguments through the function pointer
//what is funny is that this function is attached to a system, but I cannot retrieve the system to which it is attached unless I know the name or the number of it!
//all I have at hand is the MultiLevelProblem, which contains a Vector of Systems
// all I can do is put in the MultiLevelProblem a number that tells me what is the current system being solved

    
  My_exact_solution< double > exact_sol;  //this one I reproduce it here, otherwise I should pass it in the main to the MultiLevelProblem
  
  const unsigned current_system_number = ml_prob.get_current_system_number();
  
  System_assemble_flexible< system_type, real_num > (ml_prob, 
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
void System_assemble_flexible(MultiLevelProblem& ml_prob, 
                              system_type * mlPdeSys,
                              const std::vector< Unknown > &  unknowns, 
                              const Math::Function< double > & exact_sol) {
    
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  //  extract pointers to the several objects that we are going to use

  const unsigned level = mlPdeSys->GetLevelToAssemble();
  bool 			assembleMatrix 		    = mlPdeSys->GetAssembleMatrix(); 

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;

  MultiLevelSolution*   ml_sol = ml_prob._ml_sol;  
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension(); 
  const unsigned max_size_elem_dofs = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  
  

  adept::Stack & stack = FemusInit::_adeptStack;  // call the adept stack object for potential use of AD

  const assemble_jacobian< real_num, double > * assemble_jac/* = assemble_jacobian_base::build(dim)*/;
  
  //=============== Integration ========================================
  real_num_mov weight;    // must be adept if the domain is moving
  
  //=============== Geometry ========================================
  unsigned xType = BIQUADR_FE;   

  CurrentElem< real_num_mov > element(dim, msh);            // must be adept if the domain is moving, otherwise double
  
  Phi< real_num_mov > phi_coords(dim);      // must be adept if the domain is moving, otherwise double
  
  //=============== Unknowns ========================================
  
  const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();
  if (n_unknowns > 1) { std::cout << "Only scalar variable now, haven't checked with vector PDE"; abort(); }
  
  std::vector < UnknownLocal  < real_num > > unk_loc(n_unknowns);
  
  for(int u = 0; u < n_unknowns; u++) {
      unk_loc[u].initialize(dim, unknowns[u], ml_sol, mlPdeSys);
      assert(u == unk_loc[u].SolPdeIndex);  //I would like this ivar order to coincide either with SolIndex or with SolPdeIndex, otherwise too many orders... it has to be  with SolPdeIndex, because SolIndex is related to the Solution object which can have more fields... This has to match the order of the unknowns[] argument
  }
  
  vector < unsigned int > Sol_n_el_dofs_interface(n_unknowns);

  vector < double >    solu_exact_at_dofs;  solu_exact_at_dofs.reserve(max_size_elem_dofs);



//-- at dofs and quadrature points --------------- 
  vector < Phi< real_num_mov > > phi_dof_qp(n_unknowns, Phi< real_num_mov >(dim));
  
//-- 
  ElementJacRes < real_num > element_jacres(dim, unk_loc);

  vector < double >   Jac;   Jac.reserve( n_unknowns * max_size_elem_dofs * n_unknowns * max_size_elem_dofs);
  vector < real_num > Res;   Res.reserve( n_unknowns * max_size_elem_dofs);
           vector < int >       loc_to_glob_map_all_vars;   loc_to_glob_map_all_vars.reserve( n_unknowns *max_size_elem_dofs);

  RES->zero();
  if (assembleMatrix)   KK->zero();


  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
      
    short unsigned ielGeom = msh->GetElementType(iel);
    
    element.set_coords_at_dofs(iel, xType);
    
    
    
    unsigned nDofu  = msh->GetElementDofNumber(iel, unk_loc[0].SolFEType);
    
    loc_to_glob_map_all_vars.resize(nDofu);
    Res.resize(nDofu);         std::fill(Res.begin(), Res.end(), 0.);
    Jac.resize(nDofu * nDofu);  std::fill(Jac.begin(), Jac.end(), 0.);
 
       // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
         loc_to_glob_map_all_vars[i] = pdeSys->GetSystemDof(unk_loc[0].SolIndex, unk_loc[0].SolPdeIndex, i, iel);
    }
    
    
    
    solu_exact_at_dofs.resize(nDofu);
    for (unsigned i = 0; i < nDofu; i++) {
            std::vector< double > x_at_node(dim,0.);
        for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = element.get_coords_at_dofs(jdim,i);
      solu_exact_at_dofs[i] = exact_sol.value(x_at_node);
    }
    

    
    for (unsigned  u = 0; u < n_unknowns; u++) { unk_loc[u].set_elem(iel, msh, sol); }

    unsigned sum_Sol_n_el_dofs_interface = 0;
    for (unsigned  u = 0; u < n_unknowns; u++) {
        Sol_n_el_dofs_interface[u]   = unk_loc[u].Sol_n_el_dofs;
        sum_Sol_n_el_dofs_interface += Sol_n_el_dofs_interface[u]; 
        
    }
    

    assemble_jac->prepare_before_integration_loop(stack);

    
    if (dim != 2) abort(); //only implemented in 2D now

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][unk_loc[0].SolFEType]->GetGaussPointNumber(); ig++) {
        
      // *** get gauss point weight, test function and test function partial derivatives ***
  for (unsigned  u = 0; u < n_unknowns; u++) {
      static_cast<const elem_type_2D*>( msh->_finiteElement[ielGeom][unk_loc[0].SolFEType] )
                                         ->Jacobian_type_non_isoparametric< double >( static_cast<const elem_type_2D*>( msh->_finiteElement[ielGeom][xType] ), element.get_coords_at_dofs(), ig, weight, phi_dof_qp[u].phi, phi_dof_qp[u].phi_x, phi_dof_qp[u].phi_xx);
  }
  
//       msh->_finiteElement[ielGeom][SolFEType[0]]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);
      msh->_finiteElement[ielGeom][xType]->Jacobian(element.get_coords_at_dofs(), ig, weight, phi_coords.phi, phi_coords.phi_x, phi_coords.phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
               real_num solu_gss = 0.;
      vector < real_num > gradSolu_gss(dim, 0.);
      vector < double > gradSolu_exact_gss(dim, 0.);

   for (unsigned  u = 0; u < n_unknowns; u++) {
     for (unsigned i = 0; i < unk_loc[u].Sol_n_el_dofs; i++) {
        solu_gss += phi_dof_qp[u].phi[i] * unk_loc[u].Sol_eldofs[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
                gradSolu_gss[jdim] += phi_dof_qp[u].phi_x[i * dim + jdim] * unk_loc[0].Sol_eldofs[i];
          gradSolu_exact_gss[jdim] += phi_dof_qp[u].phi_x[i * dim + jdim] * solu_exact_at_dofs[i];
        }
      }
   }

      
      vector < double > x_gss(dim, 0.);
      for (unsigned i = 0; i < element.get_coords_at_dofs()[0].size(); i++) {
        for (unsigned jdim = 0; jdim < dim; jdim++) {
          x_gss[jdim] += element.get_coords_at_dofs(jdim,i) * phi_coords.phi[i];
        }          
      }
      
      
      
      for (unsigned i = 0; i < nDofu; i++) {

        real_num laplace = 0.;
        real_num laplace_weak_exact = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          laplace            +=  phi_dof_qp[0].phi_x[i * dim + jdim] * gradSolu_gss[jdim];
          laplace_weak_exact +=  phi_dof_qp[0].phi_x[i * dim + jdim] * gradSolu_exact_gss[jdim];
        }
        

// arbitrary rhs
//               double source_term = exact_sol.value(x_gss);
//         Res[i] += ( source_term * phi[i] - solu_gss * phi[i] - laplace ) * weight;
        
// manufactured Helmholtz - strong
             double helmholtz_strong_exact = exact_sol.helmholtz(x_gss);
        Res[i] += (helmholtz_strong_exact * phi_dof_qp[0].phi[i] - solu_gss * phi_dof_qp[0].phi[i] - laplace) * weight;

// manufactured Laplacian - strong
//                double laplace_strong_exact = exact_sol.laplacian(x_gss);
//         Res[i] += (- laplace_strong_exact * phi[i] - phi[i] * solu_gss - laplace) * weight;        //strong form of RHS and weak form of LHS

// manufactured Laplacian - weak
//            Res[i] += (laplace_weak_exact - phi[i] * solu_gss - laplace) * weight;                  //weak form of RHS and weak form of LHS


        
        assemble_jac->compute_jacobian_inside_integration_loop(i, dim, Sol_n_el_dofs_interface, sum_Sol_n_el_dofs_interface, unk_loc, phi_dof_qp, weight, Jac);  //rethink of these arguments when you have more unknowns
        
      
        
      } // end phi_i loop
      
      
    } // end gauss point loop

    
 assemble_jac->compute_jacobian_outside_integration_loop(stack, unk_loc, Res, Jac, loc_to_glob_map_all_vars, RES, KK);
 
    
  } //end element loop for each process

  
  RES->close();
  KK->close();

  // ***************** END ASSEMBLY *******************
}


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
                                                         const double weight, 
                                                         std::vector< double > & Jac )  const {
                                                             
                                                             
     constexpr unsigned int pos_unk = 0/*unk_vec[0].SolPdeIndex*/;                         

// *** phi_j loop ***
        for (unsigned j = 0; j < Sol_n_el_dofs[ pos_unk ]; j++) {
          /*real_num*/double laplace_jac = 0.;

          for (unsigned kdim = 0; kdim < dim; kdim++) {
            laplace_jac += (phi[ pos_unk].phi_x[i * dim + kdim] * phi[ pos_unk ].phi_x[j * dim + kdim]);
          }

          Jac[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_unk, pos_unk, i, j) ] += (laplace_jac + phi[ pos_unk ].phi[i] * phi[ pos_unk ].phi[j]) * weight;
        } // end phi_j loop

        
}

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
