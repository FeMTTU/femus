/** tutorial/Ex10
 * This example shows how to set and solve the weak form of the NonLinearPoisson problem
 *                  $\nabla \cdot (a(u)\nabla u)=f$,  $in\Omega= [-1,1]^{2} $ \\
 *                  $a(u)\nabla u \cdot n=g_{N}$, on $\partial\Omega_{\text{left}}$ \\
 *                  $u=0$,  on $\partial\Omega_{rest}$
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"

#include "adept.h"


#include "FE_convergence.hpp"

#include "../tutorial_common.hpp"



using namespace femus;



double GetExactSolutionValue(const std::vector < double >& x) {
  double pi = acos(-1.);
  return sin(pi * x[0]) * cos(0.5 * pi * x[1]); // u(x,y)=sin(pi*x)cos(pi/2*y)
};

void GetExactSolutionGradient(const std::vector < double >& x, std::vector < double >& solGrad) {
  solGrad.resize(2);  
  const double pi = acos(-1.);
  solGrad[0] =   pi*cos(pi*x[0])*cos(0.5*pi*x[1]);
  solGrad[1] =  -0.5 * sin(pi*x[0])*pi*sin(0.5*pi*x[1]);
  //solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
  //solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
};

double GetExactSolutionLaplace(const std::vector < double >& x) {
  double pi = acos(-1.);
  return (1./4.)*pi*pi*sin(pi*x[0])*cos((1./2.)*pi*x[1])*(15.*cos((1./2.)*pi*x[1])*cos((1./2.)*pi*x[1])*cos(pi*x[0])*cos(pi*x[0])
  -2.*cos(pi*x[0])*cos(pi*x[0])-7.*cos((1./2.)*pi*x[1])*cos((1./2.)*pi*x[1])-3.); // This is the source term f which must be provided as part of the pde.
};



bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceIndex, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0.;
  if(faceIndex == 1){
    dirichlet = false;
    double u = GetExactSolutionValue(x);
    std::vector < double > solGrad;
    GetExactSolutionGradient(x, solGrad);// This is gonna return "solGrad" with input "x". Carefully note that this "x" must be the coordinates of the nodes related to face=1 where x = (-1,y).
    value = -(1.+u*u) * solGrad[0]; // a(u)*u_x
  }

  return dirichlet;
}


void AssemblePoissonProblem_AD(MultiLevelProblem& ml_prob);



int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define MultiLevel object "mlMsh". 
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  
  const std::string relative_path_to_build_directory =  "../../../";
  const std::string mesh_file = relative_path_to_build_directory + Files::mesh_folder_path() + "01_gambit/2d/square/minus1-plus1_minus1-plus1/square_2x2_quad_Two_boundary_groups.neu"; 
  mlMsh.ReadCoarseMesh(mesh_file.c_str(), "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
    probably in future it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension(); // Domain dimension of the problem.
  unsigned maxNumberOfMeshes; // The number of mesh levels.

  if (dim == 2) {
    maxNumberOfMeshes = 7; 
  } else {
    maxNumberOfMeshes = 4;
  }

  std::vector < std::vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes);

  std::vector < std::vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i + 1; //We apply uniform refinement.
    unsigned numberOfSelectiveLevels = 0; // We may want to see the solution on some levels.
    mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

    // print mesh info
    mlMsh.PrintInfo();

    std::vector< Unknown > unknowns = systems__generate_list_of_scalar_unknowns_for_each_FE_family_lagrangian();


    l2Norm[i].resize( unknowns.size() );
    semiNorm[i].resize( unknowns.size() );

    for (unsigned u = 0; u < unknowns.size(); u++) {   // loop on the FE Order
      // define the multilevel solution and attach the mlMsh object to it
      MultiLevelSolution ml_sol(&mlMsh);

      // add variables to ml_sol
      ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);
      ml_sol.Initialize("All");

      // attach the boundary condition function and generate boundary data
      ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
      ml_sol.GenerateBdc( unknowns[u]._name.c_str() );

      // define the multilevel problem attach the ml_sol object to it
      MultiLevelProblem ml_prob(&ml_sol);

      ml_prob.get_systems_map().clear();
      ml_prob.set_current_system_number(0/*u*/);               //way to communicate to the assemble function, which doesn't belong to any class

      // add system Poisson in ml_prob as a Non Linear Implicit System
      NonLinearImplicitSystem& system = ml_prob.add_system < NonLinearImplicitSystem > ("NonLinearPoisson");

      // add solution "u" to system
      system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
      
        // set unknown list
        std::vector< Unknown > unknowns_vec(1);
        unknowns_vec[0] = unknowns[u]; //need to turn this into a vector

        system.set_unknown_list_for_assembly(unknowns_vec); //way to communicate to the assemble function, which doesn't belong to any class

        
      // attach the assembling function to system
      system.SetAssembleFunction(AssemblePoissonProblem_AD);

      // initilaize and solve the system
      system.init();
      system.SetOuterSolver(PREONLY);
      system.MGsolve();
    // ======= System - END ========================
      
      

      std::pair< double , double > norm = GetErrorNorm_L2_H1_with_analytical_sol(&ml_sol, unknowns_vec, GetExactSolutionValue, GetExactSolutionGradient);
      l2Norm[i][u]  = norm.first;
      semiNorm[i][u] = norm.second;
      // print solutions
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");
            
      VTKWriter vtkIO(&ml_sol);
      
      vtkIO.SetGraphVariable( unknowns[u]._name.c_str() );
      vtkIO.SetDebugOutput(true);
      vtkIO.Write(unknowns[u]._name.c_str(), Files::_application_output_directory, fe_fams_for_files[ FILES_CONTINUOUS_BIQUADRATIC ], variablesToBePrinted, i);

    }
  }

  // ======= L2 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < 3; j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < 3; j++) {
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= L2 - END  ========================

  // ======= H1 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < 3; j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < 3; j++) {
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= H1 - END ========================





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
void AssemblePoissonProblem_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object


  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
    const unsigned current_system_number = ml_prob.get_current_system_number();
    NonLinearImplicitSystem * mlPdeSys  = & ml_prob.get_system< NonLinearImplicitSystem >(current_system_number);

    // II) Unknowns of the System
  std::vector< Unknown >   unknowns = ml_prob.get_system< NonLinearImplicitSystem >(current_system_number).get_unknown_list_for_assembly();
  
  
  const unsigned level = mlPdeSys->GetLevelToAssemble(); // We have different level of meshes. we assemble the problem on the specified one.

  
  
  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->GetMeshElements();  // pointer to the elem object in msh (level)

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim))); // Return a value of unsigned // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = ml_sol->GetIndex(  unknowns[0]._name.c_str() );    // get the position of "u" in the ml_sol object
  unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex(  unknowns[0]._name.c_str() );    // get the position of "u" in the pdeSys object

  std::vector < adept::adouble >  solu; // local solution
  solu.reserve(maxSize);

  std::vector < std::vector < double > > x(dim);    // local coordinates. x is now dim x m matrix.
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned k = 0; k < dim; k++) { 
    x[k].reserve(maxSize); // dim x maxsize is reserved for x.  
  }

  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives
  
  double weight; // gauss point weight
  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim); // This is probably gradient but he is doing the life difficult for me!
  
  std::vector < adept::adouble > aRes; // local redidual vector
  aRes.reserve(maxSize);

  std::vector < int > l2GMap; // local to global mapping
  l2GMap.reserve(maxSize);
  std::vector < double > Res; // local redidual vector
  Res.reserve(maxSize);
  std::vector < double > Jac;
  Jac.reserve(maxSize * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix
  RES->zero(); // Set to zero all the entries of the Global Residual

  // element loop: each process loops only on the elements that owns
  // Adventure starts here!
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {
     
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);// number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    l2GMap.resize(nDofu);
    solu.resize(nDofu);

    for (int k = 0; k < dim; k++) {
      x[k].resize(nDofx); // 
    }

    aRes.resize(nDofu);    //resize
    std::fill(aRes.begin(), aRes.end(), 0);    //set aRes to zero

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      l2GMap[i] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        x[k][i] = (*msh->GetTopology()->_Sol[k])(xDof);      // global extraction and local storage for the element coordinates
      }
    }


    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();
           
    
  // ======= BOUNDARY - BEGIN ========================
    
    // *** Face Gauss point loop (boundary Integral) ***
    for ( unsigned jface = 0; jface < msh->GetElementFaceNumber ( iel ); jface++ ) {
      int faceIndex = el->GetBoundaryIndex(iel, jface);
      // look for boundary faces
      if ( faceIndex == 1 ) {  
        const unsigned faceGeom = msh->GetElementFaceType ( iel, jface );
        unsigned faceDofs = msh->GetElementFaceDofNumber (iel, jface, soluType);         
        std::vector < std::vector <  double> > faceCoordinates ( dim ); // A matrix holding the face coordinates rowwise.
        for ( int k = 0; k < dim; k++ ) {
          faceCoordinates[k].resize (faceDofs);
        }
        for ( unsigned i = 0; i < faceDofs; i++ ) {
          unsigned inode = msh->GetLocalFaceVertexIndex ( iel, jface, i ); // face-to-element local node mapping.
          for ( unsigned k = 0; k < dim; k++ ) {
            faceCoordinates[k][i] =  x[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
        for ( unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][soluType]->GetGaussPointNumber(); ig++ ) { 
            // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh. 
          std::vector < double> normal;
          msh->_finiteElement[faceGeom][soluType]->JacobianSur ( faceCoordinates, ig, weight, phi, phi_x, normal );
            
          std::vector < double > xg(dim,0.);
          for ( unsigned i = 0; i < faceDofs; i++ ) {
            for( unsigned k=0; k<dim; k++){
              xg[k] += phi[i] * faceCoordinates[k][i]; // xg(ig)= \sum_{i=0}^faceDofs phi[i](xig) facecoordinates[i]     
            }
          }
          double tau; // a(u)*grad_u\cdot normal
          SetBoundaryCondition( xg, "u", tau, faceIndex, 0. ); // return tau
          // *** phi_i loop ***
          for ( unsigned i = 0; i < faceDofs; i++ ) {
            unsigned inode = msh->GetLocalFaceVertexIndex ( iel, jface, i );
            aRes[inode] +=  phi[i] * tau * weight; 
          }        
        }
      }
    }   
  // ======= BOUNDARY - END ========================
  
  
  // ======= VOLUME - BEGIN ========================
  
    // *** Element Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      // Note that we dont compute the hat functions at gauss points.
      adept::adouble solu_gss = 0;
      std::vector < adept::adouble > gradSolu_gss(dim, 0.);
      std::vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i]; 

        for (unsigned k = 0; k < dim; k++) {
          gradSolu_gss[k] += phi_x[i * dim + k] * solu[i]; 
          x_gss[k] += x[k][i] * phi[i]; // We map the gausspoints from the reference element to physical element.
        }
      }

      // *** phi_i loop ***
      for (unsigned i = 0; i < nDofu; i++) {

        adept::adouble auGradUGradv = 0.;

        for (unsigned k = 0; k < dim; k++) {
          auGradUGradv   +=   phi_x[i * dim + k] * gradSolu_gss[k];
        }
        auGradUGradv *= (1. + solu_gss * solu_gss);
        double srcTerm = - GetExactSolutionLaplace(x_gss);
        aRes[i] += (srcTerm * phi[i] - auGradUGradv) * weight;

      } // end phi_i loop
    } // end gauss point loop

  // ======= VOLUME - END ========================


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

}

