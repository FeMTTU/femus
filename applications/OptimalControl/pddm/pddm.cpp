#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "LinearImplicitSystem.hpp"
#include "NumericVector.hpp"


#define NSUB_X  8
#define NSUB_Y  8


using namespace femus;


double InitialValueDS(const std::vector < double >& x) {
  return 0.;
}


bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = false;
  value = 0;
  
  const double tolerance = 1.e-5;
  
  if( x[0] < 0. + tolerance || x[0] > 1. - tolerance) {
      dirichlet = true;
        value = 0;
  }
  
  return dirichlet;
 }


void AssembleProblem(MultiLevelProblem& ml_prob);



int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

//   mlMsh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,"seventh");
  std::string infile("./input/Mesh_1.med");
  mlMsh.ReadCoarseMesh(infile.c_str(),"seventh",scalingFactor);
 /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
 
  unsigned numberOfUniformLevels = 4;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.EraseCoarseLevels(numberOfUniformLevels + numberOfSelectiveLevels - 1);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

    // ******* Print mesh *******
//   mlSol.SetWriter(VTK);  //   mlSol.GetWriter()->SetDebugOutput(true);
//   mlSol.GetWriter()->Write(files.GetOutputPath(), "biquadratic");
//   exit(0);
//   mlMsh.SetWriter(VTK);  //this doesn't work and should be removed, no application uses it
//   mlMsh.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic", meshToBePrinted);  
    // ******* End print mesh *******
  
  // add variables to mlSol
  mlSol.AddSolution("d_s", LAGRANGE,SECOND/*DISCONTINUOUS_POLYNOMIAL, ZERO*/);
  
  mlSol.Initialize("All");    // initialize all variables to zero
  
  mlSol.Initialize("d_s", InitialValueDS);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("d_s");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

// *************************************
 // this problem is defined on an open boundary mesh, and the boundary mesh can change 
 // as a function of the fracture propagation criterion.
 // Therefore, all the structures may need to be re-allocated after that
 // The workflow is:
  
 // Fill the dense matrix, and solve it
  
 // we may start from reading a mesh file 
  
 // check the propagation
  
 // if propagation occurs, re-dimensionalize all the arrays
// *************************************
  

  mlProb.SetFilesHandler(&files);
  
 // add system  in mlProb as a Linear Implicit System
  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("Frac");
 
  system.AddSolutionToSystemPDE("d_s");
 
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleProblem);

//   system.SetMaxNumberOfLinearIterations(2);
  // initialize and solve the system
  system.SetMgType(F_CYCLE/*F_CYCLE*//*M_CYCLE*/); //it doesn't matter if I use only 1 level
  system.SetOuterSolver(GMRES);
  system.init();
  
  system.MGsolve();
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(files.GetOutputPath(), "biquadratic", variablesToBePrinted);

  return 0;
}




void AssembleProblem(MultiLevelProblem& ml_prob) {

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("Frac");  
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             KK = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id();

 //***************************************************  
  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }
 //***************************************************   

 //******************** quadrature *******************************  
  double weight; 
  

 //********************* unknowns *********************** 
 //***************************************************  
  const int n_vars = mlPdeSys->GetSolPdeIndex().size();
  std::cout << "************" << n_vars << "************";
  const int solType_max = 2;  //biquadratic
  vector <double> phi_u;
  vector <double> phi_u_x; 
  vector <double> phi_u_xx;

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = mlSol->GetIndex("d_s"); 
  unsigned solFEType_u = mlSol->GetSolutionType(solIndex_u); 

  unsigned solPdeIndex_u;
  solPdeIndex_u = mlPdeSys->GetSolPdeIndex("d_s");

  vector < double >  sol_u; // local solution
  sol_u.reserve(maxSize);
  vector< int > l2GMap_u;
  l2GMap_u.reserve(maxSize);
 //***************************************************  
 //***************************************************  

  
 //***************************************************  
 //********* WHOLE SET OF VARIABLES ****************** 

  vector< int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_vars*maxSize);
  
  vector< double > Res; // local redidual vector
  Res.reserve(n_vars*maxSize);

  vector < double > Jac;
  Jac.reserve( n_vars*maxSize * n_vars*maxSize);
 //***************************************************  

  
 //********************* DATA ************************ 
  double penalty_strong = 10e+14;
 //***************************************************  
  
  
  if (assembleMatrix)  KK->zero();

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned kelGeom = msh->GetElementType(iel);

 //******************** GEOMETRY ********************* 
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);  

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);
      }
    }

   // elem average point 
    vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
    for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofx; i++) {
         elem_center[j] += x[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofx; }
 //***************************************************  
  
 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solFEType_u);
    sol_u    .resize(nDof_u);
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solFEType_u);
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
      l2GMap_u[i] = pdeSys->GetSystemDof(solIndex_u, solPdeIndex_u, i, iel);
    }
 //***************************************************  
 
 //******************** ALL VARS ********************* 
    unsigned nDof_AllVars = nDof_u; 
    int nDof_max    =  nDof_u;   // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    Res.resize(nDof_AllVars);
    std::fill(Res.begin(), Res.end(), 0.);

    Jac.resize(nDof_AllVars * nDof_AllVars);
    std::fill(Jac.begin(), Jac.end(), 0.);
    
    l2GMap_AllVars.resize(0);
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_u.begin(),l2GMap_u.end());
 //*************************************************** 
    
 //========= gauss value quantities ==================   
	double sol_u_gss = 0.;
	std::vector<double> sol_u_x_gss(dim);     std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
 //===================================================   

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
	msh->_finiteElement[kelGeom][solFEType_u]   ->Jacobian(x, ig, weight, phi_u, phi_u_x, phi_u_xx);
	
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < nDof_u; i++) {
	                                                sol_u_gss      += sol_u[i] * phi_u[i];
                   for (unsigned d = 0; d < dim; d++)   sol_u_x_gss[d] += sol_u[i] * phi_u_x[i * dim + d];
          }
	
//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
	      double laplace_res_du_u_i = 0.;
              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_u )         laplace_res_du_u_i             +=  (phi_u_x   [i * dim + kdim] * sol_u_x_gss[kdim]);
	      }
	      
//======================Residuals=======================
          // FIRST ROW
	  if (i < nDof_u)                      Res[0      + i] += - weight * ( phi_u[i] * (  -1. ) - laplace_res_du_u_i);
//======================Residuals=======================
	      
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {

              double laplace_mat_du_u = 0.;

              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_u && j < nDof_u )           laplace_mat_du_u           += (phi_u_x   [i * dim + kdim] * phi_u_x   [j * dim + kdim]);
	      }

              //============ delta_state row ============================
              //DIAG BLOCK delta_state - state
	      if ( i < nDof_u && j < nDof_u )       
		Jac[ (0 + i) * nDof_AllVars   +
		(0 + j)                         ]  += weight * laplace_mat_du_u;
              
            } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop


    RES->add_vector_blocked(Res, l2GMap_AllVars);

    if (assembleMatrix) {
      KK->add_matrix_blocked(Jac, l2GMap_AllVars, l2GMap_AllVars);
    }
    
  } //end element loop for each process

  RES->close();

  if (assembleMatrix) KK->close();
//   KK->print();
//   RES->print();

  // ***************** END ASSEMBLY *******************

  return;
}

