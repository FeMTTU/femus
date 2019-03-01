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
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "LinearImplicitSystem.hpp"
#include "Files.hpp"
#include "Math.hpp"
#include "adept.h"

// command to view matrices
// ./tutorial_ex2_c -mat_view ::ascii_info_detail
// ./tutorial_ex2_c -mat_view > matview_print_in_file.txt

using namespace femus;


// template function: definition 
template < class real_num >
void prepare_before_integration_loop(adept::Stack& stack) { }
 

// template function: specialization
template < >
void prepare_before_integration_loop< adept::adouble  > (adept::Stack & stack) { 
    
  stack.new_recording();    // start a new recording of all the operations involving adept variables

}

// 
// //********************************
// template function: explicit instantiations: not even needed because the specialization above acts also as explicit instantiation
// template void prepare_before_integration_loop< adept::adouble >(adept::Stack& stack);
// // template class prepare_before_elem_loop< double >;


template < class real_num >
void  compute_jacobian_inside_integration_loop(const unsigned i,
                                               const unsigned dim,
                                               const unsigned nDofu,
                                               const std::vector< real_num > & phi,
                                               const std::vector< real_num > &  phi_x, 
                                               const real_num weight,
                                               std::vector< double > & Jac) { };
  

template < >
void  compute_jacobian_inside_integration_loop< double >(const unsigned i,
                                                         const unsigned dim, 
                                                         const unsigned nDofu, 
                                                         const std::vector< double > &  phi,
                                                         const std::vector< double > &  phi_x, 
                                                         const double weight, 
                                                         std::vector< double > & Jac) { 

// *** phi_j loop ***
        for (unsigned j = 0; j < nDofu; j++) {
          /*real_num*/double laplace_jac = 0.;

          for (unsigned kdim = 0; kdim < dim; kdim++) {
            laplace_jac += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]);
          }

          Jac[i * nDofu + j] += (laplace_jac + phi[i] * phi[j]) * weight;
        } // end phi_j loop

        
}


// //********************************
template < class real_num >
void  compute_jacobian_outside_integration_loop(adept::Stack & stack,
                                               const std::vector< real_num > & solu,
                                               const std::vector< real_num > & Res,
                                               std::vector< double > & Jac, 
                                               const std::vector< int > & loc_to_glob_map,
                                               NumericVector*           RES,
                                               SparseMatrix*             KK
                                                                   ) { }
                                               
                                               
template < >
void  compute_jacobian_outside_integration_loop < adept::adouble > (adept::Stack & stack,
                                                                    const std::vector< adept::adouble > & solu,
                                                                    const std::vector< adept::adouble > & Res,
                                                                    std::vector< double > & Jac,
                                                                    const std::vector< int > & loc_to_glob_map,
                                                                    NumericVector*           RES,
                                                                    SparseMatrix*             KK
                                                                   ) {
    
    //copy the value of the adept::adoube Res in double Res and store
    //all the calculations were done with adept variables
    
 ///convert to vector of double to send to the global matrix
  std::vector < double > Res_double(Res.size());

  for (int i = 0; i < Res_double.size(); i++) {
      Res_double[i] = - Res[i].value();
    }


    stack.dependent(  & Res[0], Res_double.size());      // define the dependent variables
    stack.independent(&solu[0],       solu.size());    // define the independent variables
    
    stack.jacobian(&Jac[0], true);    // get the jacobian matrix (ordered by row major)

    stack.clear_independents();
    stack.clear_dependents();    
    
    RES->add_vector_blocked(Res_double, loc_to_glob_map);
    KK->add_matrix_blocked(Jac, loc_to_glob_map, loc_to_glob_map);

}


template < >
void  compute_jacobian_outside_integration_loop < double > (adept::Stack & stack,
                                               const std::vector< double > & solu,
                                               const std::vector< double > & Res,
                                               std::vector< double > & Jac,
                                               const std::vector< int > & loc_to_glob_map,
                                               NumericVector*           RES,
                                               SparseMatrix*             KK
                                                                   ) {
    
    RES->add_vector_blocked(Res, loc_to_glob_map);
    KK->add_matrix_blocked(Jac, loc_to_glob_map, loc_to_glob_map);
    
}






 
template < class type >
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

  


bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
    
  bool dirichlet = true; //dirichlet
  value = 0;

  return dirichlet;
}


template < class real_num > void AssembleProblem_interface(MultiLevelProblem & ml_prob);

template < class real_num > void AssembleProblem_flexible(MultiLevelProblem & ml_prob,
                                                    const std::string system_name,
                                                    const std::string unknown,
                                                    const Math::Function< real_num/*type*/ > & exact_sol);


 //Unknown definition  ==================
 const std::vector< Math::Unknowns_definition >  provide_list_of_unknowns() {
     
     
  std::vector< FEFamily > feFamily = {LAGRANGE, LAGRANGE,  LAGRANGE, DISCONTINOUS_POLYNOMIAL, DISCONTINOUS_POLYNOMIAL};
  std::vector< FEOrder >   feOrder = {FIRST, SERENDIPITY ,SECOND,ZERO,FIRST};

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Math::Unknowns_definition >  unknowns(feFamily.size());
 
     for (unsigned int fe = 0; fe < unknowns.size(); fe++) {
         
            std::ostringstream unk; unk << "u" << "_" << feFamily[fe] << "_" << feOrder[fe];
              unknowns[fe]._name      = unk.str();
              unknowns[fe]._fe_family = feFamily[fe];
              unknowns[fe]._fe_order  = feOrder[fe];
              
     }
 
 
   return unknowns;
     
}



template < class real_num > 
class Main_single_level : public FE_convergence< real_num, double > {
    

const MultiLevelSolution  run_main_on_single_level(const Files & files, 
                                                   const std::vector< Math::Unknowns_definition > & unknowns,  
                                                   MultiLevelMesh & ml_mesh, 
                                                   const unsigned i);
  
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

   std::string input_file = "Lshape_4.med";
//    std::string input_file = "Lshape.med";
   MultiLevelMesh ml_mesh;
//   ml_mesh.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << input_file;
  const std::string infile = mystream.str();
  
  ml_mesh.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),1.);


// set total number of levels ================  
  unsigned max_number_of_meshes;

  if (nsub[2] == 0)   max_number_of_meshes = 3;
  else                max_number_of_meshes = 4;
  

 //set coarse storage mesh (should write the copy constructor or "=" operator to copy the previous mesh) ==================
  MultiLevelMesh ml_mesh_all_levels;
//   ml_mesh_all_levels.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());
  ml_mesh_all_levels.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),1.);

 
   My_exact_solution< double > exact_sol;                                            //provide exact solution, if available ==============
   const unsigned conv_order_flag = 0;                                               //Choose how to compute the convergence order ============== //0: incremental 1: absolute (with analytical sol)  2: absolute (with projection of finest sol)...
   const unsigned norm_flag = 1;                                                     //Choose what norms to compute (//0 = only L2: //1 = L2 + H1) ==============
   std::vector< Math::Unknowns_definition > unknowns = provide_list_of_unknowns();   //provide list of unknowns ==============

    
    Main_single_level< /*adept::a*/double >  fe_convergence;
    
    fe_convergence.convergence_study(files, unknowns, SetBoundaryCondition, ml_mesh, ml_mesh_all_levels, max_number_of_meshes, norm_flag, conv_order_flag);
    

  return 0;
  
}









template < class real_num > 
const MultiLevelSolution  Main_single_level< real_num >::run_main_on_single_level(const Files & files,
                                                                                const std::vector< Math::Unknowns_definition > &  unknowns,  
                                                                                MultiLevelMesh & ml_mesh,
                                                                                const unsigned i)  {
      
      
            //Mesh  ==================
            unsigned numberOfUniformLevels = i + 1;
            unsigned numberOfSelectiveLevels = 0;
            ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
            ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);

            ml_mesh.PrintInfo();
                  
      
           //Solution  ==================
            MultiLevelSolution ml_sol_single_level(&ml_mesh); 

         for (unsigned int u = 0; u < unknowns.size(); u++) {
             
            ml_sol_single_level.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order);
            ml_sol_single_level.Initialize(unknowns[u]._name.c_str());
            ml_sol_single_level.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
            ml_sol_single_level.GenerateBdc(unknowns[u]._name.c_str());
      
            
            
            // define the multilevel problem attach the ml_sol_single_level object to it
            MultiLevelProblem mlProb(&ml_sol_single_level);

            
            mlProb.SetFilesHandler(&files);
      
      
            // add system Poisson in mlProb as a Linear Implicit System
            LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("Equation");

            // add solution "u" to system
            system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());

            
            mlProb.set_current_unknown_assembly(unknowns[u]._name); //way to communicate to the assemble function, which doesn't belong to any class
            
            // attach the assembling function to system
            system.SetAssembleFunction(AssembleProblem_interface< real_num >);

            // initialize and solve the system
            system.init();
            system.ClearVariablesToBeSolved();
            system.AddVariableToBeSolved("All");

            ml_sol_single_level.SetWriter(VTK);
            ml_sol_single_level.GetWriter()->SetDebugOutput(true);
  
//             system.SetDebugLinear(true);
//             system.SetMaxNumberOfLinearIterations(6);
//             system.SetAbsoluteLinearConvergenceTolerance(1.e-4);
            
            system.MLsolve();
      
            // ======= Print ========================
            std::vector < std::string > variablesToBePrinted;
            variablesToBePrinted.push_back(unknowns[u]._name);
            ml_sol_single_level.GetWriter()->Write(unknowns[u]._name, files.GetOutputPath(), "biquadratic", variablesToBePrinted, i);  
     

         }
         

            return ml_sol_single_level;
}



template <class real_num >
void AssembleProblem_interface(MultiLevelProblem& ml_prob) {
// this is meant to be like a tiny addition to the main function, because we cannot pass these arguments through the function pointer
    
   My_exact_solution< real_num > exact_sol;
const std::string system_name = "Equation"; //I cannot get this from the system because there may be more than one

      AssembleProblem_flexible< real_num > (ml_prob, system_name, ml_prob.get_current_unknown_assembly(), exact_sol);

}


/**
 * This function assemble the stiffnes matrix Jac and the residual vector Res
 * such that
 *                  Jac w = RES = F - Jac u0,
 * and consequently
 *        u = u0 + w satisfies Jac u = F
 **/


/**
 * This function assemble the stiffnes matrix KK and the residual vector Res
 * Using automatic differentiation for Newton iterative scheme
 *                  J(u0) w =  - F(u0)  ,
 *                  with u = u0 + w
 *                  - F = f(x) - J u = Res
 *                  J = \grad_u F
 *
 * thus
 *                  J w = f(x) - J u0
 **/
template <class real_num>
void AssembleProblem_flexible(MultiLevelProblem& ml_prob, const std::string system_name, const std::string unknown, const Math::Function< real_num > & exact_sol) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

    

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> (system_name);   // pointer to the linear implicit system 
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


  //solution variable
  unsigned soluIndex = ml_sol->GetIndex(unknown.c_str());    // get the position of "u" in the ml_sol object
  unsigned soluType  = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"
  unsigned soluPdeIndex = mlPdeSys->GetSolPdeIndex(unknown.c_str());    // get the position of "u" in the pdeSys object
  if (soluPdeIndex > 0) { std::cout << "Only scalar variable now, haven't checked with vector PDE"; abort(); }

  
  real_num weight; // gauss point weight

      
  vector < real_num >  solu;  solu.reserve(max_size_elem_dofs);
  vector < real_num >  solu_exact_at_dofs;  solu_exact_at_dofs.reserve(max_size_elem_dofs);

  vector < vector < real_num > > x(dim);  unsigned xType = BIQUADR_FE;

  for (unsigned i = 0; i < dim; i++)  x[i].reserve(max_size_elem_dofs);

//-----------------  
  vector < real_num > phi_coords;
  vector < real_num > phi_coords_x;
  vector < real_num > phi_coords_xx;

  phi_coords.reserve(max_size_elem_dofs);
  phi_coords_x.reserve(max_size_elem_dofs * dim);
  phi_coords_xx.reserve(max_size_elem_dofs * dim2);

//-----------------  
  vector < real_num > phi;
  vector < real_num > phi_x;
  vector < real_num > phi_xx;

  phi.reserve(max_size_elem_dofs);
  phi_x.reserve(max_size_elem_dofs * dim);
  phi_xx.reserve(max_size_elem_dofs * dim2);

  

  vector < int > loc_to_glob_map;  loc_to_glob_map.reserve(max_size_elem_dofs);
  vector < real_num >  Res;    Res.reserve(max_size_elem_dofs);  
  vector < double > Jac;     Jac.reserve(max_size_elem_dofs * max_size_elem_dofs);  //this has to be double, not real_num
  

  adept::Stack & stack = FemusInit::_adeptStack;  // call the adept stack object for potential use of AD

  
  KK->zero();


  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);

    // resize local arrays
    loc_to_glob_map.resize(nDofu);
    solu.resize(nDofu);
    solu_exact_at_dofs.resize(nDofu);

    for (int i = 0; i < dim; i++)    x[i].resize(nDofx);


    Res.resize(nDofu);         std::fill(Res.begin(), Res.end(), 0);
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
        std::vector< real_num > x_at_node(dim,0.);
        for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = x[jdim][i];
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
                    solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      solu_exact_at_dofs[i] = exact_sol.value(x_at_node);
                  loc_to_glob_map[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }



    prepare_before_integration_loop< real_num >(stack);

    
    if (dim != 2) abort(); //only implemented in 2D now

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        
      // *** get gauss point weight, test function and test function partial derivatives ***
      static_cast<const elem_type_2D*>( msh->_finiteElement[ielGeom][soluType] )
                                         ->Jacobian_type_non_isoparametric< real_num >( static_cast<const elem_type_2D*>( msh->_finiteElement[ielGeom][xType] ), x, ig, weight, phi, phi_x, phi_xx);
//       msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);
      msh->_finiteElement[ielGeom][xType]->Jacobian(x, ig, weight, phi_coords, phi_coords_x, phi_coords_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
               real_num solu_gss = 0.;
      vector < real_num > gradSolu_gss(dim, 0.);
      vector < real_num > gradSolu_exact_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
                gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          gradSolu_exact_gss[jdim] += phi_x[i * dim + jdim] * solu_exact_at_dofs[i];
        }
      }

      
      vector < real_num > x_gss(dim, 0.);
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
//         Res[i] += ( source_term * phi[i] - phi[i] * solu_gss - laplace ) * weight;
        
// manufactured Helmholtz - strong
             real_num helmholtz_strong_exact = exact_sol.helmholtz(x_gss);
        Res[i] += (helmholtz_strong_exact * phi[i] - solu_gss * phi[i] - laplace) * weight;

// manufactured Laplacian - strong
//                double laplace_strong_exact = exact_sol.laplacian(x_gss);
//         Res[i] += (- laplace_strong_exact * phi[i] - phi[i] * solu_gss - laplace) * weight;        //strong form of RHS and weak form of LHS

// manufactured Laplacian - weak
//            Res[i] += (laplace_weak_exact - phi[i] * solu_gss - laplace) * weight;                  //weak form of RHS and weak form of LHS


        
        compute_jacobian_inside_integration_loop< real_num > (i, dim, nDofu, phi, phi_x, weight, Jac);
        
      
        
      } // end phi_i loop
      
      
    } // end gauss point loop

    
 compute_jacobian_outside_integration_loop < real_num > (stack, solu, Res, Jac, loc_to_glob_map, RES, KK);
 
    
  } //end element loop for each process

  
  RES->close();
  KK->close();

  // ***************** END ASSEMBLY *******************
}


///@todo: check bc for discontinuous FE
///@todo: compute error in L-\infty norm
///@todo: compute nonlinear convergence rate
///@todo: compute time convergence rate, pointwise and then in norms
///@todo: uncouple Gauss from Mesh
///@todo: make non-isoparametric Jacobian routines (abstract Jacobian)
///@todo: check solver and prec choices
///@todo: create a copy constructor/operator=() for the Mesh
