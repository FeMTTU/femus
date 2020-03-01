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
                              const std::vector< Math::Unknown > &  unknowns, 
                              const Math::Function< double > & exact_sol);


 //Unknown definition  ==================
 const std::vector< Math::Unknown >  provide_list_of_unknowns() {
     
     
  std::vector< FEFamily > feFamily = {LAGRANGE, LAGRANGE,  LAGRANGE, DISCONTINUOUS_POLYNOMIAL, DISCONTINUOUS_POLYNOMIAL};
  std::vector< FEOrder >   feOrder = {FIRST, SERENDIPITY ,SECOND,ZERO,FIRST};

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Math::Unknown >  unknowns(feFamily.size());
 
     for (unsigned int fe = 0; fe < unknowns.size(); fe++) {
         
            std::ostringstream unk; unk << "u" << "_" << feFamily[fe] << "_" << feOrder[fe];
              unknowns[fe]._name      = unk.str();
              unknowns[fe]._fe_family = feFamily[fe];
              unknowns[fe]._fe_order  = feOrder[fe];
              
     }
 
 
   return unknowns;
     
}



// this is the class that I pass to the FE_convergence functions. Basically I pass a class, instead of passing a function pointer.
// notice that this is a class template, not a class with a function template
template < class real_num > 
class My_main_single_level : public Main_single_level {
    
public:
    
const MultiLevelSolution  run_on_single_level(const Files & files, 
                                                   const std::vector< Math::Unknown > & unknowns,  
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
//   std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
//   const std::string infile = mystream.str();
//   ml_mesh.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),1.);


   // ======= Unknowns ========================
   std::vector< Math::Unknown > unknowns = provide_list_of_unknowns();
   

   // ======= Normal run ========================
    My_main_single_level< /*adept::a*/double > my_main;
//  const unsigned int n_levels = 3;
//  my_main.run_on_single_level(files, unknowns, ml_mesh, n_levels); if you don't want the convergence study
    
    
   // ======= Convergence study ========================
    
   // set total number of levels ================  
   unsigned max_number_of_meshes;

   if (nsub[2] == 0)   max_number_of_meshes = 6;
   else                max_number_of_meshes = 4;
  

   //set coarse storage mesh (should write the copy constructor or "=" operator to copy the previous mesh) ==================
   MultiLevelMesh ml_mesh_all_levels;
   ml_mesh_all_levels.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());
   //   ml_mesh_all_levels.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),1.);

 
   // convergence choices ================  
   My_exact_solution<> exact_sol;                                            //provide exact solution, if available ==============
   const unsigned conv_order_flag = 0;                                               //Choose how to compute the convergence order ============== //0: incremental 1: absolute (with analytical sol)  2: absolute (with projection of finest sol)...
   const unsigned norm_flag = 1;                                                     //Choose what norms to compute (//0 = only L2: //1 = L2 + H1) ==============

   
   // object ================  
    FE_convergence<>  fe_convergence;
    
    fe_convergence.convergence_study(files, unknowns, Solution_set_boundary_conditions, ml_mesh, ml_mesh_all_levels, max_number_of_meshes, norm_flag, conv_order_flag, my_main);
    

  return 0;
  
}









template < class real_num > 
const MultiLevelSolution  My_main_single_level< real_num >::run_on_single_level(const Files & files,
                                                                                const std::vector< Math::Unknown > &  unknowns,  
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

         for (unsigned int u = 0; u < unknowns.size(); u++) {
             
            ml_sol_single_level.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order);
            ml_sol_single_level.Initialize(unknowns[u]._name.c_str());
            ml_sol_single_level.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
            ml_sol_single_level.GenerateBdc(unknowns[u]._name.c_str());
      
            
           // ======= Problem ========================
            MultiLevelProblem ml_prob(&ml_sol_single_level);
            
            ml_prob.SetFilesHandler(&files);
      
      
           // ======= System ========================
            LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("Equation");

            system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());

            std::vector< Math::Unknown > unknowns_vec(1); unknowns_vec[0] = unknowns[u]; //need to turn this into a vector
            ml_prob.set_unknown_list_for_assembly(unknowns_vec); //way to communicate to the assemble function, which doesn't belong to any class
            
            system.SetAssembleFunction(System_assemble_interface< LinearImplicitSystem, real_num >);

            // initialize and solve the system
            system.init();
            system.ClearVariablesToBeSolved();
            system.AddVariableToBeSolved("All");

            ml_sol_single_level.SetWriter(VTK);
            ml_sol_single_level.GetWriter()->SetDebugOutput(true);
  
//             system.SetDebugLinear(true);
//             system.SetMaxNumberOfLinearIterations(6);
//             system.SetAbsoluteLinearConvergenceTolerance(1.e-4);
            
            system.SetOuterSolver(PREONLY);
            system.MGsolve();
      
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
    
   My_exact_solution< double > exact_sol;  //this one I reproduce it here, otherwise I should pass it in the main to the MultiLevelProblem

  if (ml_prob.n_systems() > 1) { std::cout << "Haven't tested it yet" << std::endl; abort(); }
  
  System_assemble_flexible< system_type, real_num > (ml_prob, & ml_prob.get_system< system_type >(0), ml_prob.get_unknown_list_for_assembly(), exact_sol);

}


/**
 * This function assemble the stiffness matrix KK and the residual vector Res
 * for the Newton iterative scheme
 *                  J(u0) w = Res(u_0) = f(x) - J u_0  ,
 *                  with u = u0 + w
 *                  J = \grad_u F
 **/

template < class system_type, class real_num, class real_num_mov = double >
void System_assemble_flexible(MultiLevelProblem& ml_prob, 
                              system_type * mlPdeSys,
                              const std::vector< Math::Unknown > &  unknowns, 
                              const Math::Function< double > & exact_sol) {
    
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  //  extract pointers to the several objects that we are going to use

  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;  
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension(); 
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned max_size_elem_dofs = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)


  const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();
  if (n_unknowns > 1) { std::cout << "Only scalar variable now, haven't checked with vector PDE"; abort(); }
  
  vector < std::string >  Solname(n_unknowns);     for(unsigned ivar=0; ivar < n_unknowns; ivar++) { Solname[ivar] = unknowns[ivar]._name; }
  vector < unsigned int > SolPdeIndex(n_unknowns);
  vector < unsigned int > SolIndex(n_unknowns);  
  vector < unsigned int > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(Solname[ivar].c_str());
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }
  
   vector < unsigned int > Sol_n_el_dofs(n_unknowns);

  
  //----------- at dofs ------------------------------
  vector < vector < real_num_mov > > x(dim);  unsigned xType = BIQUADR_FE;     // must be adept if the domain is moving, otherwise double
  for (unsigned i = 0; i < dim; i++)  x[i].reserve(max_size_elem_dofs);

  vector < vector < real_num > > SolVAR_eldofs(n_unknowns);
  for(int k = 0; k < n_unknowns; k++) {    SolVAR_eldofs[k].reserve(max_size_elem_dofs);  }
  
  vector < double >    solu_exact_at_dofs;  solu_exact_at_dofs.reserve(max_size_elem_dofs);


  //------------ at quadrature points ---------------------
  real_num_mov weight;    // must be adept if the domain is moving
//-----------------  
  vector < double > phi_coords;
  vector < real_num_mov > phi_coords_x;   // must be adept if the domain is moving
  vector < real_num_mov > phi_coords_xx;  // must be adept if the domain is moving

  phi_coords.reserve(max_size_elem_dofs);
  phi_coords_x.reserve(max_size_elem_dofs * dim);
  phi_coords_xx.reserve(max_size_elem_dofs * dim2);

//-----------------  
  vector < double > phi;
  vector < double > phi_x;
  vector < double > phi_xx;

  phi.reserve(max_size_elem_dofs);
  phi_x.reserve(max_size_elem_dofs * dim);
  phi_xx.reserve(max_size_elem_dofs * dim2);

  

  vector < int >       loc_to_glob_map;  loc_to_glob_map.reserve(max_size_elem_dofs);
  vector < real_num >  Res;                          Res.reserve(max_size_elem_dofs);  
  vector < double >    Jac;                          Jac.reserve(max_size_elem_dofs * max_size_elem_dofs);
  

  adept::Stack & stack = FemusInit::_adeptStack;  // call the adept stack object for potential use of AD

  const assemble_jacobian< real_num, double > assemble_jac;

  RES->zero();
  KK->zero();


  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, SolFEType[0]);
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);

    // resize local arrays
    loc_to_glob_map.resize(nDofu);

    solu_exact_at_dofs.resize(nDofu);

    for (int i = 0; i < dim; i++)    x[i].resize(nDofx);

    Res.resize(nDofu);         std::fill(Res.begin(), Res.end(), 0.);
    Jac.resize(nDofu * nDofu);  std::fill(Jac.begin(), Jac.end(), 0.);

    
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    } 
    
     
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
        std::vector< double > x_at_node(dim,0.);
        for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = x[jdim][i];
      solu_exact_at_dofs[i] = exact_sol.value(x_at_node);
         loc_to_glob_map[i] = pdeSys->GetSystemDof(SolIndex[0], SolPdeIndex[0], i, iel);
    }

      for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
       Sol_n_el_dofs[k] = ndofs_unk;
       SolVAR_eldofs[k].resize(ndofs_unk);
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);
       SolVAR_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);
//        JACDof[i + k *nDofsD]/*[k][i]*/ = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);
      }
    }
    
    unsigned sum_Sol_n_el_dofs = 0;
    for (unsigned  k = 0; k < n_unknowns; k++) { sum_Sol_n_el_dofs += Sol_n_el_dofs[k]; }



    assemble_jac.prepare_before_integration_loop(stack);

    
    if (dim != 2) abort(); //only implemented in 2D now

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][SolFEType[0]]->GetGaussPointNumber(); ig++) {
        
      // *** get gauss point weight, test function and test function partial derivatives ***
      static_cast<const elem_type_2D*>( msh->_finiteElement[ielGeom][SolFEType[0]] )
                                         ->Jacobian_type_non_isoparametric< double >( static_cast<const elem_type_2D*>( msh->_finiteElement[ielGeom][xType] ), x, ig, weight, phi, phi_x, phi_xx);
//       msh->_finiteElement[ielGeom][SolFEType[0]]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);
      msh->_finiteElement[ielGeom][xType]->Jacobian(x, ig, weight, phi_coords, phi_coords_x, phi_coords_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
               real_num solu_gss = 0.;
      vector < real_num > gradSolu_gss(dim, 0.);
      vector < double > gradSolu_exact_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * SolVAR_eldofs[0][i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
                gradSolu_gss[jdim] += phi_x[i * dim + jdim] * SolVAR_eldofs[0][i];
          gradSolu_exact_gss[jdim] += phi_x[i * dim + jdim] * solu_exact_at_dofs[i];
        }
      }

      
      vector < double > x_gss(dim, 0.);
      for (unsigned i = 0; i < nDofx; i++) {
        for (unsigned jdim = 0; jdim < dim; jdim++) {
          x_gss[jdim] += x[jdim][i] * phi_coords[i];
        }          
      }
      
      
      
      for (unsigned i = 0; i < nDofu; i++) {

        real_num laplace = 0.;
        real_num laplace_weak_exact = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          laplace            +=  phi_x[i * dim + jdim] * gradSolu_gss[jdim];
          laplace_weak_exact +=  phi_x[i * dim + jdim] * gradSolu_exact_gss[jdim];
        }
        

// arbitrary rhs
//               double source_term = exact_sol.value(x_gss);
//         Res[i] += ( source_term * phi[i] - solu_gss * phi[i] - laplace ) * weight;
        
// manufactured Helmholtz - strong
             double helmholtz_strong_exact = exact_sol.helmholtz(x_gss);
        Res[i] += (helmholtz_strong_exact * phi[i] - solu_gss * phi[i] - laplace) * weight;

// manufactured Laplacian - strong
//                double laplace_strong_exact = exact_sol.laplacian(x_gss);
//         Res[i] += (- laplace_strong_exact * phi[i] - phi[i] * solu_gss - laplace) * weight;        //strong form of RHS and weak form of LHS

// manufactured Laplacian - weak
//            Res[i] += (laplace_weak_exact - phi[i] * solu_gss - laplace) * weight;                  //weak form of RHS and weak form of LHS


        
        assemble_jac.compute_jacobian_inside_integration_loop(i, dim, Sol_n_el_dofs, sum_Sol_n_el_dofs, phi, phi_x, weight, Jac);  //rethink of these arguments when you have more unknowns
        
      
        
      } // end phi_i loop
      
      
    } // end gauss point loop

    
 assemble_jac.compute_jacobian_outside_integration_loop(stack, SolVAR_eldofs, Res, Jac, loc_to_glob_map, RES, KK);
 
    
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
                                               const std::vector< int > & loc_to_glob_map,
                                               NumericVector*           RES,
                                               SparseMatrix*             KK
                                                                   )  const {
    
    RES->add_vector_blocked(Res, loc_to_glob_map);
    KK->add_matrix_blocked(Jac, loc_to_glob_map, loc_to_glob_map);
    
}

 // template specialization for double
 // this is where you fill the jacobian in traditional (much faster) way
template < >
 void  assemble_jacobian< double, double >::compute_jacobian_inside_integration_loop (
                                                         const unsigned i,
                                                         const unsigned dim, 
                                                         const std::vector < unsigned int > Sol_n_el_dofs,
                                                         const unsigned int sum_Sol_n_el_dofs,
                                                         const std::vector< double > &  phi,
                                                         const std::vector< double > &  phi_x, 
                                                         const double weight, 
                                                         std::vector< double > & Jac )  const { 

// *** phi_j loop ***
        for (unsigned j = 0; j < Sol_n_el_dofs[0]; j++) {
          /*real_num*/double laplace_jac = 0.;

          for (unsigned kdim = 0; kdim < dim; kdim++) {
            laplace_jac += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]);
          }

          Jac[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, /*SolPdeIndex[0]*/ 0, /*SolPdeIndex[0]*/ 0, i, j) ] += (laplace_jac + phi[i] * phi[j]) * weight;
        } // end phi_j loop

        
}


}//end namespace
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
