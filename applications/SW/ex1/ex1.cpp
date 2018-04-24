/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution varables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "PetscMatrix.hpp"

#include "LinearImplicitSystem.hpp"

#include "slepceps.h"

using namespace femus;

double InitalValueU(const std::vector < double >& x)
{
  return x[0] + x[1];
}

double InitalValueP(const std::vector < double >& x)
{
  return x[0];
}

double InitalValueT(const std::vector < double >& x)
{
  return x[1];
}


bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time)
{
  bool dirichlet = false; //dirichlet
  if( facename == 1 || facename == 2) dirichlet = true;
  value = 0.;
  return dirichlet;
}


void ETD(MultiLevelProblem& ml_prob, const unsigned& numberOfLayers);


int main(int argc, char** args)
{

  SlepcInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;



  unsigned numberOfUniformLevels = 3;
  unsigned numberOfSelectiveLevels = 0;

  mlMsh.GenerateCoarseBoxMesh(2, 0, 0, -1, 1, 0., 0., 0., 0., EDGE3, "seventh");

  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  unsigned NumberOfLayers = 2;

  for(unsigned i = 0; i < NumberOfLayers; i++) {
    char name[10];
    sprintf(name, "h%d", i);
    mlSol.AddSolution(name, LAGRANGE, FIRST);
    sprintf(name, "v%d", i);
    mlSol.AddSolution(name, LAGRANGE, FIRST);
  }

  mlSol.Initialize("All");    // initialize all varaibles to zero

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  MultiLevelProblem ml_prob(&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ("SW");
  for(unsigned i = 0; i < NumberOfLayers; i++) {
    char name[10];
    sprintf(name, "h%d", i);
    system.AddSolutionToSystemPDE(name);
    sprintf(name, "v%d", i);
    system.AddSolutionToSystemPDE(name);
  }
  system.init();
  
  ETD(ml_prob,NumberOfLayers);

  mlSol.SetWriter(VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "linear", print_vars);


  return 0;
}


void ETD(MultiLevelProblem& ml_prob, const unsigned& NLayers)
{

  adept::Stack& s = FemusInit::_adeptStack;

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("SW");   // pointer to the linear implicit system named "Poisson"

  unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1u;

  Mesh* msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solIndexh(NLayers);
  std::vector < unsigned > solPdeIndexh(NLayers);

  std::vector < unsigned > solIndexv(NLayers);
  std::vector < unsigned > solPdeIndexv(NLayers);

  vector< int > l2GMap; // local to global mapping

  for(unsigned i = 0; i < NLayers; i++) {
    char name[10];
    sprintf(name, "h%d", i);
    solIndexh[i] = mlSol->GetIndex(name); // get the position of "hi" in the sol object
    solPdeIndexh[i] = mlPdeSys->GetSolPdeIndex(name); // get the position of "hi" in the pdeSys object

    sprintf(name, "v%d", i);
    solIndexv[i] = mlSol->GetIndex(name); // get the position of "vi" in the sol object
    solPdeIndexv[i] = mlPdeSys->GetSolPdeIndex(name); // get the position of "vi" in the pdeSys object
  }

  unsigned solTypeh = mlSol->GetSolutionType(solIndexh[0]);    // get the finite element type for "hi"
  unsigned solTypev = mlSol->GetSolutionType(solIndexv[0]);    // get the finite element type for "vi"

  vector < double > x;    // local coordinates
  vector< vector < adept::adouble > > solh(NLayers);    // local coordinates
  vector< vector < adept::adouble > > solv(NLayers);    // local coordinates
  
  vector< vector < bool > > bdch(NLayers);    // local coordinates
  vector< vector < bool > > bdcv(NLayers);    // local coordinates

  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector < vector< adept::adouble > > aResh(NLayers);
  vector < vector< adept::adouble > > aResv(NLayers);

  KK->zero();
  RES->zero();

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofh  = msh->GetElementDofNumber(iel, solTypeh);    // number of solution element dofs
    unsigned nDofv  = msh->GetElementDofNumber(iel, solTypev);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    l2GMap.resize(NLayers * (nDofh + nDofv) );

    for(unsigned i = 0; i < NLayers; i++) {
      solh[i].resize(nDofh);
      solv[i].resize(nDofv);
      bdch[i].resize(nDofh);
      bdcv[i].resize(nDofh);
      
      aResh[i].resize(nDofh);    //resize
      std::fill(aResh[i].begin(), aResh[i].end(), 0);    //set aRes to zero
      aResv[i].resize(nDofv);    //resize
      std::fill(aResv[i].begin(), aResv[i].end(), 0);    //set aRes to zero
    }
    x.resize(nDofx);

    //local storage of global mapping and solution
    for(unsigned i = 0; i < nDofh; i++) {
      unsigned solDofh = msh->GetSolutionDof(i, iel, solTypeh);    // global to global mapping between solution node and solution dof
      for(unsigned j = 0; j < NLayers; j++) {
        solh[j][i] = (*sol->_Sol[solIndexh[j]])(solDofh);      // global extraction and local storage for the solution
	bdch[j][i] = ( (*sol->_Bdc[solIndexh[j]])(solDofh) < 1.5)? true:false;
        l2GMap[ j * (nDofh + nDofv) + i] = pdeSys->GetSystemDof(solIndexh[j], solPdeIndexh[j], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    for(unsigned i = 0; i < nDofv; i++) {
      unsigned solDofv = msh->GetSolutionDof(i, iel, solTypev);    // global to global mapping between solution node and solution dof
      for(unsigned j = 0; j < NLayers; j++) {
        solv[j][i] = (*sol->_Sol[solIndexv[j]])(solDofv);      // global extraction and local storage for the solution
	bdcv[j][i] = ( (*sol->_Bdc[solIndexv[j]])(solDofv) < 1.5)? true:false;
        l2GMap[ nDofh + j * (nDofh + nDofv) + i] = pdeSys->GetSystemDof(solIndexv[j], solPdeIndexv[j], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    for(unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
      x[i] = (*msh->_topology->_Sol[0])(xDof);      // global extraction and local storage for the element coordinates
    }

    s.new_recording();
    
    double dx = x[nDofx - 1] - x[0];
    
    for(unsigned k = 0; k < NLayers; k++){
      for (unsigned i = 0; i < nDofh; i++){
	if(!bdch[k][i]){
	  aResh[k][i] = dx;
	  for (unsigned j = 0; j < nDofh; j++){
	    double sign = ( i == j)? -1.:1;
	    aResh[k][i] += sign * solh[k][j]/(dx*dx);
	  }
	}
      }
      for (unsigned i = 0; i < nDofv; i++){
	if(!bdcv[k][i]){
	  aResv[k][i] = dx;
	  for (unsigned j = 0; j < nDofv; j++){
	    double sign = ( i == j)? -1.:1;
	    aResv[k][i] += sign * solv[k][j]/(dx*dx);
	  }
	}
      }       
    }
    
    vector< double > Res(NLayers * (nDofh + nDofv)); // local redidual vector
    
      
    unsigned counter = 0;
    for(unsigned k = 0; k < NLayers; k++){
      for(int i = 0; i < nDofh; i++) {
	Res[counter] =  aResh[k][i].value();
	counter++;
      }
      for(int i = 0; i < nDofv; i++) {
	Res[counter] =  aResv[k][i].value();
	counter++;
      }
    }

    RES->add_vector_blocked(Res, l2GMap);


    for(unsigned k = 0; k < NLayers; k++){
      // define the dependent variables
      s.dependent(&aResh[k][0], nDofh);
      s.dependent(&aResv[k][0], nDofv);
      
      // define the independent variables
      s.independent(&solh[k][0], nDofh);
      s.independent(&solv[k][0], nDofv);
    }
    
    // get the jacobian matrix (ordered by row major )
    vector < double > Jac(NLayers * (nDofh + nDofv) * NLayers * (nDofh + nDofv));
    s.jacobian(&Jac[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

    s.clear_independents();
    s.clear_dependents();
    
  }
  
  std::cout<<"BBBBBBBBBBBBBBBBBBBBBBBBBBB"<<std::endl;
  
  RES->close();
  KK->close();
  
//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//   PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//   PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(),viewer);
//   double a;
//   std::cin>>a;
}



