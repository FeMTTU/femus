#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "LinearImplicitSystem.hpp"
#include "NumericVector.hpp"

using namespace femus;

double InitialValueThom(const std::vector < double >& x) {
  return 0.;
}

double InitialValueThomAdj(const std::vector < double >& x) {
  return 0.;
}

double InitialValueTcont(const std::vector < double >& x) {
  return 0.;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;

  if (faceName == 2)
    dirichlet = false;

  return dirichlet;
}


void AssembleLiftRestrProblem(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  // read coarse level mesh and generate finers level meshes
  mlMsh.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 3;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("Thom", LAGRANGE, SECOND);
  mlSol.AddSolution("ThomAdj", LAGRANGE, SECOND);
  mlSol.AddSolution("Tcont", LAGRANGE, SECOND);

  mlSol.Initialize("All");    // initialize all varaibles to zero

  mlSol.Initialize("Thom", InitialValueThom);
  mlSol.Initialize("ThomAdj", InitialValueThomAdj);
  mlSol.Initialize("Tcont", InitialValueTcont);    // note that this initialization is the same as piecewise constant element
 
  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("Thom");
  mlSol.GenerateBdc("ThomAdj");
  mlSol.GenerateBdc("Tcont");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  
 // add system  in mlProb as a Linear Implicit System
  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("LiftRestr");
 
  system.AddSolutionToSystemPDE("Thom");  
  system.AddSolutionToSystemPDE("ThomAdj");  
//   system.AddSolutionToSystemPDE("Tcont");  
  
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleLiftRestrProblem);

  // initilaize and solve the system
  system.init();
  system.solve();
  
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("Thom");
  variablesToBePrinted.push_back("ThomAdj");
  variablesToBePrinted.push_back("Tcont");

  VTKWriter vtkIO(&mlSol);
  vtkIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  GMVWriter gmvIO(&mlSol);
  variablesToBePrinted.push_back("all");
  gmvIO.SetDebugOutput(false);
  gmvIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  return 0;
}



void AssembleLiftRestrProblem(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const unsigned levelMax = mlPdeSys->GetLevelMax();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

   //*************************** 
  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }
 //*************************** 

 //*************************** 
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);
 //*************************** 
  
  
 //*************************** 
  unsigned solIndexThom;
  solIndexThom = mlSol->GetIndex("Thom");    // get the position of "Thom" in the ml_sol object
  unsigned solTypeThom = mlSol->GetSolutionType(solIndexThom);    // get the finite element type for "Thom"

  unsigned solPdeIndexThom;
  solPdeIndexThom = mlPdeSys->GetSolPdeIndex("Thom");    // get the position of "Thom" in the pdeSys object

  vector < double >  solThom; // local solution
  solThom.reserve(maxSize);
  vector< int > l2GMap_Thom;
  l2GMap_Thom.reserve(maxSize);
 //*************************** 

  
 //*************************** 
  unsigned solIndexThomAdj;
  solIndexThomAdj = mlSol->GetIndex("ThomAdj");    // get the position of "Thom" in the ml_sol object
  unsigned solTypeThomAdj = mlSol->GetSolutionType(solIndexThomAdj);    // get the finite element type for "Thom"

  unsigned solPdeIndexThomAdj;
  solPdeIndexThomAdj = mlPdeSys->GetSolPdeIndex("ThomAdj");    // get the position of "Thom" in the pdeSys object

  vector < double >  solThomAdj; // local solution
  solThomAdj.reserve(maxSize);
 vector< int > l2GMap_ThomAdj;
  l2GMap_ThomAdj.reserve(maxSize);
  //*************************** 
  
  const int n_vars = 2;
  


 //*************************** 
  vector< int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_vars*maxSize);
  
  vector< double > Res; // local redidual vector
  Res.reserve(n_vars*maxSize);

  vector < double > Jac;
  Jac.reserve( n_vars*maxSize * n_vars*maxSize);
 //*************************** 

  if (assembleMatrix)
    KK->zero(); // Set to zero all the entries of the Global Matrix

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->IS_Mts2Gmt_elem_offset[iproc]; iel < msh->IS_Mts2Gmt_elem_offset[iproc + 1]; iel++) {

    unsigned kel = msh->IS_Mts2Gmt_elem[iel]; // mapping between paralell dof and mesh dof
    short unsigned kelGeom = el->GetElementType(kel);    // element geometry type
    unsigned nDofx = el->GetElementDofNumber(kel, xType);    // number of coordinate element dofs

    unsigned nDofThom     = el->GetElementDofNumber(kel, solTypeThom);    // number of solution element dofs
    unsigned nDofThomAdj  = el->GetElementDofNumber(kel, solTypeThomAdj);    // number of solution element dofs

    unsigned nDof_AllVars = nDofThom + nDofThomAdj; 
    
    // resize local arrays
    solThom    .resize(nDofThom);
    l2GMap_Thom.resize(nDofThom);

    solThomAdj    .resize(nDofThomAdj);
    l2GMap_ThomAdj.resize(nDofThomAdj);

    l2GMap_AllVars.resize(0);

 //*************************** 
    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    Res.resize(nDof_AllVars);    //resize
    std::fill(Res.begin(), Res.end(), 0);    //set Res to zero

    Jac.resize(nDof_AllVars * nDof_AllVars);    //resize
    std::fill(Jac.begin(), Jac.end(), 0);    //set Jac to zero
 //*************************** 

    // local storage of global mapping and solution
    for (unsigned i = 0; i < solThom.size(); i++) {
      unsigned iNode = el->GetMeshDof(kel, i, solTypeThom);    // local to global solution node
      unsigned solDofThom = msh->GetMetisDof(iNode, solTypeThom);    // global to global mapping between solution node and solution dof
      solThom[i] = (*sol->_Sol[solIndexThom])(solDofThom);      // global extraction and local storage for the solution
      l2GMap_Thom[i] = pdeSys->GetKKDof(solIndexThom, solPdeIndexThom, iNode);    // global to global mapping between solution node and pdeSys dof
    }

    for (unsigned i = 0; i < solThomAdj.size(); i++) {
      unsigned iNode = el->GetMeshDof(kel, i, solTypeThomAdj);    // local to global solution node
      unsigned solDofThomAdj = msh->GetMetisDof(iNode, solTypeThomAdj);    // global to global mapping between solution node and solution dof
      solThomAdj[i] = (*sol->_Sol[solIndexThomAdj])(solDofThomAdj);      // global extraction and local storage for the solution
      l2GMap_ThomAdj[i] = pdeSys->GetKKDof(solIndexThomAdj, solPdeIndexThomAdj, iNode);    // global to global mapping between solution node and pdeSys dof
    }
    
    
//**** dof composition all vars *********************** 
//           for(int i=0; i < n_vars;i++){
// 	l2GMap_AllVars.insert( l2GMap_AllVars.end(), dofsVAR[i].begin(), dofsVAR[i].end() );
//       }
//*************************** 
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_Thom.begin(),l2GMap_Thom.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_ThomAdj.begin(),l2GMap_ThomAdj.end());

 //*************************** 
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, xType);    // local to global coordinates node
      unsigned xDof  = msh->GetMetisDof(iNode, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_coordinate->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }
 //*************************** 

    if (level == levelMax || !el->GetRefinedElementIndex(kel)) {      // do not care about this if now (it is used for the AMR)

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solTypeThom]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[kelGeom][solTypeThom]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);

        // evaluate the solution, the solution derivatives and the coordinates in the gauss point
        double solThom_gss = 0;
        vector < double > gradSolu_gss(dim, 0.);
        vector < double > x_gss(dim, 0.);

        for (unsigned i = 0; i < nDofThom; i++) {
          solThom_gss += phi[i] * solThom[i];

          for (unsigned jdim = 0; jdim < dim; jdim++) {
            gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solThom[i];
            x_gss[jdim] += x[jdim][i] * phi[i];
          }
        }

        // *** phi_i loop ***
        for (unsigned i = 0; i < nDofThom; i++) {

          double laplace_rhs = 0.;

          for (unsigned jdim = 0; jdim < dim; jdim++) {
            laplace_rhs   +=  phi_x[i * dim + jdim] * gradSolu_gss[jdim];
          }

          double srcTerm = 10.;
          Res[i]            += (srcTerm * phi[i] - laplace_rhs) * weight;
          Res[nDofThom + i] += (srcTerm * phi[i] - laplace_rhs) * weight;

          if (assembleMatrix) {
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDofThom; j++) {
              double laplace_mat = 0.;

              for (unsigned kdim = 0; kdim < dim; kdim++) {
                laplace_mat += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]) * weight;
              }

              Jac[i * (2*nDofThom) + j] += laplace_mat;
              Jac[(2*nDofThom)*(nDofThom) + i * (2* nDofThom) +(nDofThom + j)] += laplace_mat;
            } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
      } // end gauss point loop
    } // endif single element not refined or fine grid loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store
    RES->add_vector_blocked(Res, l2GMap_AllVars);

    if (assembleMatrix) {
      //store K in the global matrix KK
      KK->add_matrix_blocked(Jac, l2GMap_AllVars, l2GMap_AllVars);
    }
  } //end element loop for each process

  RES->close();

  if (assembleMatrix) KK->close();

  // ***************** END ASSEMBLY *******************
}
