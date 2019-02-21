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

 




bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;

//   if (faceName == 2)
//     dirichlet = false;

  return dirichlet;
}


void AssemblePoissonProblem(MultiLevelProblem& ml_prob);
void AssemblePoissonProblem_AD(MultiLevelProblem& ml_prob);
void AssemblePoissonProblem_AD_flexible(MultiLevelProblem& ml_prob, const std::string system_name, const std::string unknown);

 

 //Unknown definition  ==================
 const std::vector< FE_convergence::Unknowns_definition >  provide_list_of_unknowns() {
     
     
  std::vector< FEFamily > feFamily = {LAGRANGE, LAGRANGE,  LAGRANGE/*, DISCONTINOUS_POLYNOMIAL*//*, DISCONTINOUS_POLYNOMIAL*/};
  std::vector< FEOrder >   feOrder = {FIRST, SERENDIPITY ,SECOND/*,ZERO, FIRST*/};

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< FE_convergence::Unknowns_definition >  unknowns(feFamily.size());
 
     for (unsigned int fe = 0; fe < unknowns.size(); fe++) {
         
            std::ostringstream unk; unk << "u" << "_" << feFamily[fe] << "_" << feOrder[fe];
              unknowns[fe]._name      = unk.str();
              unknowns[fe]._fe_family = feFamily[fe];
              unknowns[fe]._fe_order  = feOrder[fe];
              
     }
 
 
   return unknowns;
     
}


const MultiLevelSolution  run_main_on_single_level(const std::vector< FE_convergence::Unknowns_definition > & unknowns,  const Files & files, MultiLevelMesh & ml_mesh, const unsigned i);
  


int main(int argc, char** args) {

  // ======= Init ==========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Files =========================
  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();
 
  // ======= Quad Rule ========================
  std::string fe_quad_rule("seventh");
 /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    In the future it is not going to be an argument of the mesh function   */


  
  const ElemType geom_elem_type = QUAD9;
  const std::vector< unsigned int > nsub = {2,2,0};
  const std::vector< double >      xyz_min = {-0.5,-0.5,0.};
  const std::vector< double >      xyz_max = { 0.5, 0.5,0.};

   // ======= Mesh ========================
  MultiLevelMesh ml_mesh;
  ml_mesh.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());
//   ml_mesh.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);


  
   //Mesh all levels ==================
  unsigned max_number_of_meshes;

  if (nsub[2] == 0) {
    max_number_of_meshes = 6;
  } else {
    max_number_of_meshes = 4;
  }

  MultiLevelMesh ml_mesh_all_levels;
  ml_mesh_all_levels.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());
//   ml_mesh_all_levels.ReadCoarseMesh("./input/square_quad.neu", "seventh", scalingFactor);

  
  
    std::vector< FE_convergence::Unknowns_definition > unknowns = provide_list_of_unknowns();
  
    const unsigned norm_flag = 1; //0 only L2: //1 L2+H1
    vector < vector < vector < double > > > norms = FE_convergence::initialize_vector_of_norms( unknowns.size(), max_number_of_meshes,  norm_flag );

     MultiLevelSolution         ml_sol_all_levels = FE_convergence::prepare_convergence_study(unknowns, ml_mesh_all_levels, max_number_of_meshes, SetBoundaryCondition);
    
            
       for (int i = 0; i < max_number_of_meshes; i++) {   // loop on the mesh level
                  
            const MultiLevelSolution ml_sol_single_level  =   run_main_on_single_level(unknowns, files, ml_mesh, i);

        FE_convergence::prolongate_coarser_and_compute_error_norm_on_level( & ml_sol_single_level, & ml_sol_all_levels, unknowns, i, norm_flag, norms);

        
   }   //end h refinement
   
 
  
       FE_convergence::output_convergence_rate_all(unknowns, norms, norm_flag, max_number_of_meshes);
   

  return 0;
}








  const MultiLevelSolution  run_main_on_single_level(const std::vector< FE_convergence::Unknowns_definition > &  unknowns,  const Files & files, MultiLevelMesh & ml_mesh, const unsigned i)  {
      
      
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
            LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("Poisson");

            // add solution "u" to system
            system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());

            
            mlProb.set_current_unknown_assembly(unknowns[u]._name); //way to communicate to the assemble function, which doesn't belong to any class
            
            // attach the assembling function to system
            system.SetAssembleFunction(AssemblePoissonProblem_AD);

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
            variablesToBePrinted.push_back(unknowns[u]._name/*"All"*/);
            ml_sol_single_level.GetWriter()->Write(unknowns[u]._name, files.GetOutputPath(), "biquadratic", variablesToBePrinted, i);  
     

         }
         

            return ml_sol_single_level;
}




/**
 * This function assemble the stiffnes matrix Jac and the residual vector Res
 * such that
 *                  Jac w = RES = F - Jac u0,
 * and consequently
 *        u = u0 + w satisfies Jac u = F
 **/
void AssemblePoissonProblem(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Poisson");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex = 0; // ml_sol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex =  0; //mlPdeSys->GetSolPdeIndex("u");    // get the position of "u" in the pdeSys object

  vector < double >  solu; // local solution
  solu.reserve(maxSize);

  vector < double >  solu_exact_at_dofs; // local solution
  solu_exact_at_dofs.reserve(maxSize);

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);

  vector< double > Res; // local redidual vector
  Res.reserve(maxSize);

  vector< int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  vector < double > Jac;
  Jac.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    l2GMap.resize(nDofu);
    solu.resize(nDofu);
    solu_exact_at_dofs.resize(nDofu);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    Res.resize(nDofu);          std::fill(Res.begin(), Res.end(), 0.); 
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
        std::vector<double> x_at_node(dim,0.);
        for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = x[jdim][i];
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      solu_exact_at_dofs[i] = FE_convergence::GetExactSolutionValue(x_at_node);
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }



    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double solu_gss = 0;
      vector < double > gradSolu_gss(dim, 0.);
      vector < double > gradSolu_exact_gss(dim, 0.);
      vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
                gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          gradSolu_exact_gss[jdim] += phi_x[i * dim + jdim] * solu_exact_at_dofs[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        double laplace = 0.;
        double laplace_weak_exact = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          laplace   +=  phi_x[i * dim + jdim] * gradSolu_gss[jdim];
          laplace_weak_exact +=  phi_x[i * dim + jdim] * gradSolu_exact_gss[jdim];
        }

//         double srcTerm = - GetExactSolutionLaplace(x_gss);
//         Res[i] += (srcTerm * phi[i] - laplace) * weight;
        Res[i] += (laplace_weak_exact - laplace) * weight;

        // *** phi_j loop ***
        for (unsigned j = 0; j < nDofu; j++) {
          laplace = 0.;

          for (unsigned kdim = 0; kdim < dim; kdim++) {
            laplace += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]) * weight;
          }

          Jac[i * nDofu + j] += laplace;
        } // end phi_j loop

      } // end phi_i loop
    } // end gauss point loop


    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}


void AssemblePoissonProblem_AD(MultiLevelProblem& ml_prob) {
    
    
    AssemblePoissonProblem_AD_flexible(ml_prob,"Poisson", ml_prob.get_current_unknown_assembly());

}



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
void AssemblePoissonProblem_AD_flexible(MultiLevelProblem& ml_prob, const std::string system_name, const std::string unknown) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object


  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> (system_name);   // pointer to the linear implicit system 
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex = /*0;*/ml_sol->GetIndex(unknown.c_str());    // get the position of "u" in the ml_sol object
  unsigned soluType  = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex = mlPdeSys->GetSolPdeIndex(unknown.c_str());    // get the position of "u" in the pdeSys object
  if (soluPdeIndex > 0) { std::cout << "Only scalar variable now"; abort(); }
      
      
  vector < adept::adouble >  solu; // local solution
  solu.reserve(maxSize);
  
  vector < adept::adouble >  solu_exact_at_dofs; // local solution
  solu_exact_at_dofs.reserve(maxSize);

  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }

  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);

  vector< adept::adouble > aRes; // local redidual vector
  aRes.reserve(maxSize);

  vector< int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  vector< double > Res; // local redidual vector
  Res.reserve(maxSize);
  vector < double > Jac;
  Jac.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    l2GMap.resize(nDofu);
    solu.resize(nDofu);
    solu_exact_at_dofs.resize(nDofu);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    aRes.resize(nDofu);    //resize
    std::fill(aRes.begin(), aRes.end(), 0);    //set aRes to zero

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    } 
    
     
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
        std::vector<double> x_at_node(dim,0.);
        for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = x[jdim][i];
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
                    solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      solu_exact_at_dofs[i] = FE_convergence::GetExactSolutionValue(x_at_node);
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }



    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      adept::adouble solu_gss = 0;
      vector < adept::adouble > gradSolu_gss(dim, 0.);
      vector < adept::adouble > gradSolu_exact_gss(dim, 0.);
      vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
                gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          gradSolu_exact_gss[jdim] += phi_x[i * dim + jdim] * solu_exact_at_dofs[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        adept::adouble laplace = 0.;
        adept::adouble laplace_weak_exact = 0.;

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          laplace            +=  phi_x[i * dim + jdim] * gradSolu_gss[jdim];
          laplace_weak_exact +=  phi_x[i * dim + jdim] * gradSolu_exact_gss[jdim];
        }
        

//         double srcTerm = - GetExactSolutionLaplace(x_gss);
//         aRes[i] += (srcTerm * phi[i] - laplace) * weight;        //strong form of RHS nad weak form of LHS
        aRes[i] += (laplace_weak_exact - laplace) * weight;         //weak form of RHS nad weak form of LHS

      } // end phi_i loop
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    Res.resize(nDofu);    //resize

    for (int i = 0; i < nDofu; i++) {
      Res[i] = - aRes[i].value();
    }

    RES->add_vector_blocked(Res, l2GMap);



    // define the dependent variables
    s.dependent(&aRes[0], nDofu);

    // define the independent variables
    s.independent(&solu[0], nDofu);

    // get the jacobian matrix (ordered by row major )
    Jac.resize(nDofu * nDofu);    //resize
    s.jacobian(&Jac[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

    s.clear_independents();
    s.clear_dependents();

  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}





