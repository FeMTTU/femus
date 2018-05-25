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
#include <slepcmfn.h>

using namespace femus;

double rho[10]={1025,1027,1028}; // kg/m^3

double rho1[10]={1029,1029,1029}; // kg/m^3

//double rho[10]={1025,1030,1035};
double mu = 1.5 * 1.0e-3; // pa s 

double aa = 6371000;  //radius of earth [m]
double H_0 = 2400;    // max depth of ocean [m]
double H_shelf = 100; //depth of continental shelf [m]

double z_c = aa;      //z coordinate of center point
double bb = 1250000;

double phi = 0.1;     //relative width of continental shelf: 10%

double dt = 100.; //= dx / maxWaveSpeed * 0.85;

const unsigned NumberOfLayers = 2;

const double hRest[3]={40,50,90};

//const double hRest[3]={90,90};

double InitalValueV(const std::vector < double >& x)
{
  return 0;
}


double InitalValueH0(const std::vector < double >& x)
{
  double b = hRest[0];
  if(NumberOfLayers == 1){
    double zz = sqrt(aa * aa - x[0] * x[0]); // z coordinate of points on sphere
    double dd = aa * acos((zz * z_c) / (aa * aa)); // distance to center point on sphere [m]
    double hh = 1 - dd * dd / (bb * bb);
    b = ( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );
  }
  return b + ( 2.* exp(-(x[0] / 100000.) * (x[0] / 100000.)));
}

double InitalValueH1(const std::vector < double >& x)
{
  
  if(NumberOfLayers == 2){
    double zz = sqrt(aa * aa - x[0] * x[0]); // z coordinate of points on sphere
    double dd = aa * acos((zz * z_c) / (aa * aa)); // distance to center point on sphere [m]
    double hh = 1 - dd * dd / (bb * bb);
    double b = ( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );
    return b - hRest[1];
  }
  return hRest[1];
}

double InitalValueH2(const std::vector < double >& x)
{
  double zz = sqrt(aa * aa - x[0] * x[0]); // z coordinate of points on sphere
  double dd = aa * acos((zz * z_c) / (aa * aa)); // distance to center point on sphere [m]
  double hh = 1 - dd * dd / (bb * bb);
  double b = ( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );
  return b - hRest[2];
}

double InitalValueT0(const std::vector < double >& x)
{
  return 20;
}

double InitalValueT1(const std::vector < double >& x)
{
  return 15;
}

double InitalValueT2(const std::vector < double >& x)
{
  return 5;
}


double InitalValueB(const std::vector < double >& x)
{
  double zz = sqrt(aa * aa - x[0] * x[0]); // z coordinate of points on sphere
  double dd = aa * acos((zz * z_c) / (aa * aa)); // distance to center point on sphere [m]
  double hh = 1 - dd * dd / (bb * bb);
  return ( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );
}


bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time)
{
  bool dirichlet = false; //dirichlet
  if( facename == 1 || facename == 2) dirichlet = true;
  value = 0.;
  return dirichlet;
}


void ETD(MultiLevelProblem& ml_prob);


int main(int argc, char** args)
{

  SlepcInitialize(&argc, &args, PETSC_NULL, PETSC_NULL);

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;



  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;

  unsigned nx = static_cast<unsigned>(floor(pow(2.,3) + 0.5));
  
  double length = 2. * 1465700.;

  mlMsh.GenerateCoarseBoxMesh(nx, 0, 0, -length / 2, length / 2, 0., 0., 0., 0., EDGE3, "seventh");

  //mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  for(unsigned i = 0; i < NumberOfLayers; i++) {
    char name[10];
    sprintf(name, "h%d", i);
    mlSol.AddSolution(name, DISCONTINOUS_POLYNOMIAL, ZERO);
    sprintf(name, "v%d", i);
    mlSol.AddSolution(name, LAGRANGE, FIRST);
    sprintf(name, "T%d", i);
    mlSol.AddSolution(name, DISCONTINOUS_POLYNOMIAL, ZERO);
    sprintf(name, "HT%d", i);
    mlSol.AddSolution(name, DISCONTINOUS_POLYNOMIAL, ZERO);
  }
  
  mlSol.AddSolution("b", DISCONTINOUS_POLYNOMIAL, ZERO, 1, false);

  mlSol.AddSolution("eta", DISCONTINOUS_POLYNOMIAL, ZERO, 1, false);

  mlSol.Initialize("All");
  
  
  mlSol.Initialize("h0",InitalValueH0);
  mlSol.Initialize("T0",InitalValueT0);
  if(NumberOfLayers > 1){
    mlSol.Initialize("h1",InitalValueH1);
    mlSol.Initialize("T1",InitalValueT1);
    if(NumberOfLayers > 2){
      mlSol.Initialize("h2",InitalValueH2);
      mlSol.Initialize("T2",InitalValueT2);
    }
  }
  
  for(unsigned i = 0; i < NumberOfLayers; i++) {
    char name[10];
    sprintf(name, "v%d", i);
    mlSol.Initialize(name, InitalValueV);
  }

  mlSol.Initialize("b", InitalValueB);
  
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
    sprintf(name, "HT%d", i);
    system.AddSolutionToSystemPDE(name);
  }
  system.init();

  mlSol.SetWriter(VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back("All");
  //mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "linear", print_vars, 0);

  unsigned numberOfTimeSteps = 2000;
  for(unsigned i = 0; i < numberOfTimeSteps; i++) {
    ETD(ml_prob);
    mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "linear", print_vars, (i + 1)/1);
  }
  return 0;
}


void ETD(MultiLevelProblem& ml_prob)
{

  const unsigned& NLayers = NumberOfLayers;

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
  NumericVector* EPS = pdeSys->_EPS; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solIndexh(NLayers);
  std::vector < unsigned > solPdeIndexh(NLayers);

  std::vector < unsigned > solIndexv(NLayers);
  std::vector < unsigned > solPdeIndexv(NLayers);
  
  std::vector < unsigned > solIndexHT(NLayers);
  std::vector < unsigned > solPdeIndexHT(NLayers);
  
  
   std::vector < unsigned > solIndexT(NLayers);

  vector< int > l2GMapRow; // local to global mapping
  vector< int > l2GMapColumn; // local to global mapping
  
  for(unsigned i = 0; i < NLayers; i++) {
    char name[10];
    sprintf(name, "h%d", i);
    solIndexh[i] = mlSol->GetIndex(name); // get the position of "hi" in the sol object
    solPdeIndexh[i] = mlPdeSys->GetSolPdeIndex(name); // get the position of "hi" in the pdeSys object

    sprintf(name, "v%d", i);
    solIndexv[i] = mlSol->GetIndex(name); // get the position of "vi" in the sol object
    solPdeIndexv[i] = mlPdeSys->GetSolPdeIndex(name); // get the position of "vi" in the pdeSys object
    
    sprintf(name, "HT%d", i);
    solIndexHT[i] = mlSol->GetIndex(name); // get the position of "Ti" in the sol object
    solPdeIndexHT[i] = mlPdeSys->GetSolPdeIndex(name); // get the position of "Ti" in the pdeSys object
    
    sprintf(name, "T%d", i);
    solIndexT[i] = mlSol->GetIndex(name); // get the position of "Ti" in the sol object
    
  }

  unsigned solTypeh = mlSol->GetSolutionType(solIndexh[0]);    // get the finite element type for "hi"
  unsigned solTypev = mlSol->GetSolutionType(solIndexv[0]);    // get the finite element type for "vi"
  unsigned solTypeHT = mlSol->GetSolutionType(solIndexHT[0]);    // get the finite element type for "Ti"
  
//   vector < double > x;    // local coordinates
//   vector < adept::adouble > solhm(NLayers); 
//   vector < adept::adouble > solh(NLayers);    // local coordinates
//   vector < adept::adouble > solhp(NLayers); 
//   vector < adept::adouble > solvm(NLayers);    // local coordinates
//   vector < adept::adouble > solv(NLayers);    // local coordinates
//   vector < adept::adouble > solvp(NLayers);    // local coordinates
//   vector < adept::adouble > solHTm(NLayers);    // local coordinates
//   vector < adept::adouble > solHT(NLayers);    // local coordinates
//   vector < adept::adouble > solHTp(NLayers);    // local coordinates
//       
//   unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
// 
//   vector< adept::adouble > aResh(NLayers);
//   vector< adept::adouble > aResv(NLayers);
//   vector< adept::adouble > aResHT(NLayers);
//   
//   KK->zero();
//   RES->zero();
// 
//   MatSetOption((static_cast<PetscMatrix*>(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//   
//   double maxWaveSpeed = 0.;
//   
//   double dx;
  
  for(unsigned k=0; k<NumberOfLayers; k++){
    for(unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++){
      double valueT = (*sol->_Sol[solIndexT[k]])(i);
      double valueH = (*sol->_Sol[solIndexh[k]])(i);
            
      double valueHT = valueT * valueH;
    
      sol->_Sol[solIndexHT[k]]->set(i, valueHT);
    }
    sol->_Sol[solIndexHT[k]]->close();
  }
    
  unsigned start = msh->_dofOffset[solTypeHT][iproc];
  unsigned end = msh->_dofOffset[solTypeHT][iproc + 1];
  for(unsigned i =  start; i <  end; i++){
    
    vector < double > x;    // local coordinates
    vector < adept::adouble > solhm(NLayers); 
    vector < adept::adouble > solh(NLayers);    // local coordinates
    vector < adept::adouble > solhp(NLayers); 
    vector < adept::adouble > solvm(NLayers);    // local coordinates
    vector < adept::adouble > solv(NLayers);    // local coordinates
    vector < adept::adouble > solvp(NLayers);    // local coordinates
    vector < adept::adouble > solHTm(NLayers);    // local coordinates
    vector < adept::adouble > solHT(NLayers);    // local coordinates
    vector < adept::adouble > solHTp(NLayers);    // local coordinates
      
    unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

    vector< adept::adouble > aResh(NLayers);
    vector< adept::adouble > aResv(NLayers);
    vector< adept::adouble > aResHT(NLayers);
   
    KK->zero();
    RES->zero();

    MatSetOption((static_cast<PetscMatrix*>(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  
    double maxWaveSpeed = 0.;
  
    double dx;
    
    unsigned bc1 = (i == start)? 0 : 1;
    unsigned bc2 = (i == end-1)? 0 : 1;
    
    l2GMapRow.resize(2 * NLayers);
    l2GMapColumn.resize( 2 * ( 2 + bc1 + bc2) * NLayers);
    
    std::fill(aResh.begin(), aResh.end(), 0);    //set aRes to zero
    std::fill(aResHT.begin(), aResHT.end(), 0);  //set aRes to zero
    
    x.resize(2);
    for(unsigned j = 0; j < NLayers; j++) {
      
      solh[j] = (*sol->_Sol[solIndexh[j]])(i);      
      solHT[j] = (*sol->_Sol[solIndexHT[j]])(i);  
      l2GMapRow[j] = pdeSys->GetSystemDof(solIndexh[j], solPdeIndexh[j], 0, i);   
      l2GMapRow[NLayers + j] = pdeSys->GetSystemDof(solIndexHT[j], solPdeIndexHT[j], 0, i);   
      
      l2GMapColumn[j] = pdeSys->GetSystemDof(solIndexh[j], solPdeIndexh[j], 0, i);   
      l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof(solIndexHT[j], solPdeIndexHT[j], 0, i);   
      
      solvm[j] = (*sol->_Sol[solIndexv[j]])(i); 
      solvp[j] = (*sol->_Sol[solIndexv[j]])(i+1); 

      l2GMapColumn[2 * NLayers + j] = pdeSys->GetSystemDof(solIndexv[j], solPdeIndexv[j], 0, i);   
      l2GMapColumn[3 * NLayers + j] = pdeSys->GetSystemDof(solIndexv[j], solPdeIndexv[j], 1, i);   
            
      if(i > start){ 
	solhm[j] = (*sol->_Sol[solIndexh[j]])(i-1);      
	solHTm[j] = (*sol->_Sol[solIndexHT[j]])(i-1); 
	
	l2GMapColumn[4 * NLayers + j] = pdeSys->GetSystemDof(solIndexh[j], solPdeIndexh[j], 0, i-1);   
	l2GMapColumn[5 * NLayers + j] = pdeSys->GetSystemDof(solIndexHT[j], solPdeIndexHT[j], 0, i-1);   
	
      }      
      
      if( i < end - 1){ 
	solhp[j] = (*sol->_Sol[solIndexh[j]])(i+1);      
	solHTp[j] = (*sol->_Sol[solIndexHT[j]])(i+1);      
	
	l2GMapColumn[(4 + 2 * bc1) * NLayers + j] = pdeSys->GetSystemDof(solIndexh[j], solPdeIndexh[j], 0, i+1);   
	l2GMapColumn[(5 + 2 * bc1) * NLayers + j] = pdeSys->GetSystemDof(solIndexHT[j], solPdeIndexHT[j], 0, i+1);   
      }
    }
    s.new_recording();
    
    for(unsigned j = 0; j < 2; j++) {
      unsigned xDof  = msh->GetSolutionDof(j, i, xType);    // global to global mapping between coordinates node and coordinate dof
      x[j] = (*msh->_topology->_Sol[0])(xDof);      // global extraction and local storage for the element coordinates
    }
    dx = x[1] - x[0];
       
    double xmid = 0.5 * (x[0] + x[1]);
    double zz = sqrt(aa * aa - xmid * xmid); // z coordinate of points on sphere
    double dd = aa * acos((zz * z_c) / (aa * aa)); // distance to center point on sphere [m]
    double hh = 1 - dd * dd / (bb * bb);
    double b = ( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) ); 
        
    double hTot = 0.;
    for(unsigned k = 0; k < NLayers; k++) {
      hTot += solh[k].value();
    }
     
    std::vector < double > hALE(NLayers, 0.); 
      
    hALE[0] = hRest[0] + (hTot - b);
    for(unsigned k = 1; k < NLayers - 1; k++){
      hALE[k] = hRest[k];
    }
    hALE[NLayers - 1] = b - hRest[NLayers - 1];
       
    std::vector < double > w(NLayers+1, 0.);
	  
    for(unsigned k = NLayers; k>1; k--){
      w[k-1] = w[k] -  solh[k-1].value() * (solvp[k-1].value() - solvm[k-1].value() )/dx - ( hALE[k-1] - solh[k-1].value()) / dt;//TODO
    }
      
      
    for(unsigned k = 0; k < NLayers; k++) {
      
      if( i > start ){
	aResh[k] += 0.5 * (solhm[k] + solh[k]) * solvm[k] / dx; 
      }
      if(i < end - 1){
	aResh[k] -= 0.5 * (solh[k] + solhp[k]) * solvp[k] / dx; 
      }
      aResh[k] += w[k+1] - w[k];	    
      
	
      if( i > start ){
        aResHT[k] += 0.5 * (solHTm[k] + solHT[k]) * solvm[k]  / dx; 
      }
      if(i < end - 1){
	aResHT[k] -= 0.5 * (solHT[k] + solHTp[k]) * solvp[k]  / dx; 
      }
	
      if(k<NLayers-1){
	aResHT[k] += w[k+1] * 0.5 * (solHT[k]/solh[k] + solHT[k+1]/solh[k+1]);
      }
      if( k > 0){
	aResHT[k] -= w[k] * 0.5 * (solHT[k-1]/solh[k-1] + solHT[k]/solh[k] );
      }		
    }
      
    vector< double > Res(NLayers * 2); // local redidual vector
    for(unsigned k = 0; k < NLayers; k++) {
      Res[k] =  aResh[k].value();
      Res[NLayers + k] =  aResHT[k].value();
    }

    RES->add_vector_blocked(Res, l2GMapRow);
  
    s.dependent(&aResh[0], NLayers);
    s.dependent(&aResHT[0], NLayers);

    // define the independent variables
    s.independent(&solh[0], NLayers);
    s.independent(&solHT[0], NLayers);
    s.independent(&solvm[0], NLayers);
    s.independent(&solvp[0], NLayers);
    if(i > start){
      s.independent(&solhm[0], NLayers);
      s.independent(&solHTm[0], NLayers);
    }
    if(i < end - 1){
      s.independent(&solhp[0], NLayers);
      s.independent(&solHTp[0], NLayers);
    }
    
    // get the jacobian matrix (ordered by row major )
    vector < double > Jac(NLayers * 2 * NLayers * 4 * (2 + bc1 + bc2) );
    s.jacobian(&Jac[0], true);
    
    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMapRow, l2GMapColumn);

    s.clear_independents();
    s.clear_dependents();   
         
  }
  
  
  
  start = msh->_dofOffset[solTypev][iproc] + 1;
  end = msh->_dofOffset[solTypev][iproc + 1] - 1;
  for(unsigned i =  start; i <  end; i++){
    
    vector < double > x;    // local coordinates
    vector < adept::adouble > solhm(NLayers); 
    vector < adept::adouble > solh(NLayers);    // local coordinates
    vector < adept::adouble > solhp(NLayers); 
    vector < adept::adouble > solvm(NLayers);    // local coordinates
    vector < adept::adouble > solv(NLayers);    // local coordinates
    vector < adept::adouble > solvp(NLayers);    // local coordinates
    vector < adept::adouble > solHTm(NLayers);    // local coordinates
    vector < adept::adouble > solHT(NLayers);    // local coordinates
    vector < adept::adouble > solHTp(NLayers);    // local coordinates
      
    unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

    vector< adept::adouble > aResh(NLayers);
    vector< adept::adouble > aResv(NLayers);
    vector< adept::adouble > aResHT(NLayers);
  
    KK->zero();
    RES->zero();

    MatSetOption((static_cast<PetscMatrix*>(KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  
    double maxWaveSpeed = 0.;
  
    double dx;
        
    l2GMapRow.resize(NLayers);
    l2GMapColumn.resize( 7 * NLayers);
    
    std::fill(aResv.begin(), aResv.end(), 0);    //set aRes to zero
        
    for(unsigned j = 0; j < NLayers; j++) {
      
      solv[j] = (*sol->_Sol[solIndexv[j]])(i); 
      l2GMapRow[j] = pdeSys->GetSystemDof(solIndexv[j], solPdeIndexv[j], 0, i);  
      l2GMapColumn[j] = pdeSys->GetSystemDof(solIndexv[j], solPdeIndexv[j], 0, i);   
      
      solvm[j] = (*sol->_Sol[solIndexv[j]])(i-1); 
      solvp[j] = (*sol->_Sol[solIndexv[j]])(i+1); 
      l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof(solIndexv[j], solPdeIndexv[j], 0, i-1);   
      l2GMapColumn[2 * NLayers + j] = pdeSys->GetSystemDof(solIndexv[j], solPdeIndexv[j], 1, i);   
      
               
      solhm[j] = (*sol->_Sol[solIndexh[j]])(i-1);      
      solHTm[j] = (*sol->_Sol[solIndexHT[j]])(i-1); 
	
      l2GMapColumn[3 * NLayers + j] = pdeSys->GetSystemDof(solIndexh[j], solPdeIndexh[j], 0, i-1);   
      l2GMapColumn[4 * NLayers + j] = pdeSys->GetSystemDof(solIndexHT[j], solPdeIndexHT[j], 0, i-1);   
	
      solhp[j]  = (*sol->_Sol[solIndexh[j]])(i);      
      solHTp[j] = (*sol->_Sol[solIndexHT[j]])(i);      
	
      l2GMapColumn[5 * NLayers + j] = pdeSys->GetSystemDof(solIndexh[j], solPdeIndexh[j], 0, i);   
      l2GMapColumn[6 * NLayers + j] = pdeSys->GetSystemDof(solIndexHT[j], solPdeIndexHT[j], 0, i);   
      
    }
    s.new_recording();
    
    std::vector <double> xm(2);
    std::vector <double> xp(2);
    for(unsigned j = 0; j < 2; j++) {
      unsigned xDofm  = msh->GetSolutionDof(j, i-1, xType);    // global to global mapping between coordinates node and coordinate dof
      xm[j] = (*msh->_topology->_Sol[0])(xDofm);      // global extraction and local storage for the element coordinates
      
      unsigned xDofp  = msh->GetSolutionDof(j, i, xType);    // global to global mapping between coordinates node and coordinate dof
      xp[j] = (*msh->_topology->_Sol[0])(xDofp);      // global extraction and local storage for the element coordinates
    }
    double dxm = xm[1] - xm[0];
    double dxp = xp[1] - xp[0];
       
    double xmidm = 0.5 * (xm[0] + xm[1]);
    double zzm = sqrt(aa * aa - xmidm * xmidm); // z coordinate of points on sphere
    double ddm = aa * acos((zzm * z_c) / (aa * aa)); // distance to center point on sphere [m]
    double hhm = 1 - ddm * ddm / (bb * bb);
    double bm = ( H_shelf + H_0 / 2 * (1 + tanh(hhm / phi)) ); 
    
    double xmidp = 0.5 * (xp[0] + xp[1]);
    double zzp = sqrt(aa * aa - xmidp * xmidp); // z coordinate of points on sphere
    double ddp = aa * acos((zzp * z_c) / (aa * aa)); // distance to center point on sphere [m]
    double hhp = 1 - ddp * ddp / (bb * bb);
    double bp = ( H_shelf + H_0 / 2 * (1 + tanh(hhp / phi)) ); 
        
    double hTotm = 0.;
    double hTotp = 0.;    
    
    double beta = 0.;
    
    std::vector < adept::adouble > Pm(NLayers);
    std::vector < adept::adouble > Pp(NLayers);
    for(unsigned k = 0; k < NLayers; k++) {
      hTotm += solhm[k].value();
      hTotp += solhp[k].value();
      adept::adouble rhokm = rho1[k] - beta * solHTm[k]/solhm[k];
      adept::adouble rhokp = rho1[k] - beta * solHTp[k]/solhp[k];
      Pm[k] = - rhokm * 9.81 * bm; // bottom topography
      Pp[k] = - rhokp * 9.81 * bp; // bottom topography
      for( unsigned j = 0; j < NLayers; j++){
	 adept::adouble rhojm = (j <= k) ? (rho1[j] - beta * solHTm[k]/solhm[k]) : rhokm;
	 adept::adouble rhojp = (j <= k) ? (rho1[j] - beta * solHTp[k]/solhp[k]) : rhokp;
	 Pm[k] += rhojm * 9.81 * solhm[k];
	 Pp[k] += rhojp * 9.81 * solhp[k];
      }
      Pm[k] /= 1024;
      Pp[k] /= 1024;
    }
     
    std::vector < double > hALEm(NLayers, 0.); 
    std::vector < double > hALEp(NLayers, 0.); 
      
    hALEm[0] = hRest[0] + (hTotm - bm);
    hALEp[0] = hRest[0] + (hTotp - bp);
    for(unsigned k = 1; k < NLayers - 1; k++){
      hALEm[k] = hRest[k];
      hALEp[k] = hRest[k];
      
    }
    hALEm[NLayers - 1] = bm - hRest[NLayers - 1];
    hALEp[NLayers - 1] = bp - hRest[NLayers - 1];
       
    std::vector < double > wm(NLayers+1, 0.);
    std::vector < double > wp(NLayers+1, 0.);
	  
    for(unsigned k = NLayers; k>1; k--){
      wm[k-1] = wm[k] -  solhm[k-1].value() * (solv[k-1].value() - solvm[k-1].value() )/dxm - ( hALEm[k-1] - solhm[k-1].value()) / dt;
      wp[k-1] = wp[k] -  solhp[k-1].value() * (solvp[k-1].value() - solv[k-1].value() )/dxp - ( hALEp[k-1] - solhp[k-1].value()) / dt;
    }
      
      
    for(unsigned k = 0; k < NLayers; k++) {
      adept::adouble vMidm = 0.5 * (solvm[k] + solv[k]);
      adept::adouble fvm = 0.5 * vMidm * vMidm + Pm[k];
      aResv[k] +=  fvm / dxp;
      
      adept::adouble vMidp = 0.5 * (solv[k] + solvp[k]);
      adept::adouble fvp = 0.5 * vMidp * vMidp + Pp[k];
      
      aResv[k] -=  fvp / dxp;
      
      if ( k > 0 ){
	aResv[k] -= ( wm[k] + wp[k]) * (solv[k-1]-solv[k]) / (solhm[k-1]+solhm[k]+solhp[k-1]+solhp[k] );
      }
      if (k < NLayers - 1) {
	aResv[k] -= ( wm[k+1] + wp[k+1]) * (solv[k]-solv[k+1]) / (solhm[k]+solhm[k+1]+solhp[k]+solhp[k+1] ); 
      }
    }
      
    vector< double > Res(NLayers); // local redidual vector
    for(unsigned k = 0; k < NLayers; k++) {
      Res[k] =  aResv[k].value();
    }

    RES->add_vector_blocked(Res, l2GMapRow);
  
    s.dependent(&aResv[0], NLayers);
    
    // define the independent variables
    s.independent(&solv[0], NLayers);
    s.independent(&solvm[0], NLayers);
    s.independent(&solvp[0], NLayers);
    s.independent(&solhm[0], NLayers);
    s.independent(&solHTm[0], NLayers);
    s.independent(&solhp[0], NLayers);
    s.independent(&solHTp[0], NLayers);
        
    // get the jacobian matrix (ordered by row major )
    vector < double > Jac(NLayers * NLayers * 7 );
    s.jacobian(&Jac[0], true);
    
    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMapRow, l2GMapColumn);

    s.clear_independents();
    s.clear_dependents();   
         
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  /*
  
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofh  = msh->GetElementDofNumber(iel, solTypeh);    // number of solution element dofs
    unsigned nDofv  = msh->GetElementDofNumber(iel, solTypev);    // number of solution element dofs
    unsigned nDofHT = msh->GetElementDofNumber(iel, solTypeHT);     // number of coordinate element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);        // number of coordinate element dofs

    unsigned nDofs = nDofh + nDofv + nDofHT;

    // resize local arrays
    l2GMap.resize(NLayers * (nDofs) );

    for(unsigned i = 0; i < NLayers; i++) {
      solh[i].resize(nDofh);
      solv[i].resize(nDofv);
      solHT[i].resize(nDofHT);
      solT[i].resize(nDofHT);
      bdch[i].resize(nDofh);
      bdcv[i].resize(nDofv);
      bdcHT[i].resize(nDofHT);
      aResh[i].resize(nDofh);    //resize
      std::fill(aResh[i].begin(), aResh[i].end(), 0);    //set aRes to zero
      aResv[i].resize(nDofv);    //resize
      std::fill(aResv[i].begin(), aResv[i].end(), 0);    //set aRes to zero
      aResHT[i].resize(nDofHT);    //resize
      std::fill(aResHT[i].begin(), aResHT[i].end(), 0);    //set aRes to zero
    }
    x.resize(nDofx);

    //local storage of global mapping and solution
    for(unsigned i = 0; i < nDofh; i++) {
      unsigned solDofh = msh->GetSolutionDof(i, iel, solTypeh);    // global to global mapping between solution node and solution dof
      for(unsigned j = 0; j < NLayers; j++) {
        solh[j][i] = (*sol->_Sol[solIndexh[j]])(solDofh);      // global extraction and local storage for the solution
        bdch[j][i] = ( (*sol->_Bdc[solIndexh[j]])(solDofh) < 1.5) ? true : false;
        l2GMap[ j * nDofs + i] = pdeSys->GetSystemDof(solIndexh[j], solPdeIndexh[j], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    for(unsigned i = 0; i < nDofv; i++) {
      unsigned solDofv = msh->GetSolutionDof(i, iel, solTypev);    // global to global mapping between solution node and solution dof
      for(unsigned j = 0; j < NLayers; j++) {
        solv[j][i] = (*sol->_Sol[solIndexv[j]])(solDofv);      // global extraction and local storage for the solution
        bdcv[j][i] = ( (*sol->_Bdc[solIndexv[j]])(solDofv) < 1.5) ? true : false;
        l2GMap[ j * nDofs + nDofh + i] = pdeSys->GetSystemDof(solIndexv[j], solPdeIndexv[j], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    for(unsigned i = 0; i < nDofHT; i++) {
      unsigned solDofHT = msh->GetSolutionDof(i, iel, solTypeHT);    // global to global mapping between solution node and solution dof
      unsigned solDofh = msh->GetSolutionDof(i, iel, solTypeh);
      for(unsigned j = 0; j < NLayers; j++) {
        solHT[j][i] = (*sol->_Sol[solIndexHT[j]])(solDofHT);   // global extraction and local storage for the solution
	solT[j][i] = (*sol->_Sol[solIndexT[j]])(solDofHT);   // global extraction and local storage for the solution
        bdcHT[j][i] = ( (*sol->_Bdc[solIndexHT[j]])(solDofHT) < 1.5) ? true : false;
        l2GMap[ j * nDofs + nDofh + nDofv + i] = pdeSys->GetSystemDof(solIndexHT[j], solPdeIndexHT[j], i, iel);    // global to global mapping between solution node and pdeSys dof
	//std::cout << solT[j][i] << " ";
      }
    }

    s.new_recording();
    
    for(unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof
      x[i] = (*msh->_topology->_Sol[0])(xDof);      // global extraction and local storage for the element coordinates
    }
    

    double xmid = 0.5 * (x[0] + x[1]);
    //std::cout << xmid << std::endl;
    double zz = sqrt(aa * aa - xmid * xmid); // z coordinate of points on sphere
    double dd = aa * acos((zz * z_c) / (aa * aa)); // distance to center point on sphere [m]
    double hh = 1 - dd * dd / (bb * bb);
    double b = ( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) ); 
    
   
   
      
    
    double hTot = 0.;
    
    double beta = 0.2;
    
     std::vector < adept::adouble > P(NLayers);
    for(unsigned k = 0; k < NLayers; k++) {
      hTot += solh[k][0].value();
      adept::adouble rhok = rho1[k] - beta * solT[k][0];
      P[k] = - rhok * 9.81 * b; // bottom topography
      for( unsigned j = 0; j < NLayers; j++){
	 adept::adouble rhoj = (j <= k) ? (rho1[j] - beta * solT[j][0]) : rhok;
	 P[k] += rhoj * 9.81 * solh[j][0];
      }
      P[k] /= 1024;
    }
       
    std::vector < double > hALE(NLayers, 0.); 
      
    hALE[0] = hRest[0] + (hTot - b);
    for(unsigned k = 1; k < NLayers - 1; k++){
      hALE[k] = hRest[k];
    }
    hALE[NLayers - 1] = b - hRest[NLayers - 1];
       
    std::vector < double > w(NLayers+1, 0.);
    
     dx = x[1] - x[0];
     
    for(unsigned k = NLayers; k>1; k--){
      w[k-1] = w[k] -  solh[k-1][0].value() * (solv[k-1][1].value() - solv[k-1][0].value() )/dx- ( hALE[k-1] - solh[k-1][0].value()) / dt;
      //std::cout << hALE[k-1] << " " << w[k-1] << " ";
    }
    //std::cout<<std::endl;
    
    std::vector < adept::adouble > zMid(NLayers);
    for(unsigned k = 0; k < NLayers; k++) {
      zMid[k] = -b + solh[k][0]/2;
      for(unsigned i = k+1; i < NLayers; i++) {
        zMid[k] += solh[i][0];
      }
    }
    
    double celerity = sqrt( 9.81 * hTot);
    for(unsigned i = 0; i<NLayers; i++){
      double vmid = 0.5 * ( solv[i][0].value() + solv[i][1].value() );
      maxWaveSpeed = ( maxWaveSpeed > fabs(vmid) + fabs(celerity) )?  maxWaveSpeed : fabs(vmid) + fabs(celerity);
    }
       

    for(unsigned k = 0; k < NLayers; k++) {
      if(!bdch[k][0]) {
        for (unsigned j = 0; j < nDofv; j++) {
          double sign = ( j == 0) ? 1. : -1;
          aResh[k][0] += sign * solh[k][0] * solv[k][j] / dx; 
	}
	aResh[k][0] += w[k+1] - w[k];	    
      }
      adept::adouble vMid = 0.5 * (solv[k][0] + solv[k][1]);
      adept::adouble fv = 0.5 * vMid * vMid + P[k];
      for (unsigned i = 0; i < nDofv; i++) {
        if(!bdcv[k][i]) {
          double sign = ( i == 0) ? -1. : 1;
          aResv[k][i] += sign * fv / dx;
	  if ( k > 0 ){
	    aResv[k][i] -= 0.5 * ( w[k]   * 0.5 * ( solv[k-1][i] - solv[k][i] ) / (0.5 * ( solh[k-1][0] + solh[k][0] ) ) );
	  }
	  if (k < NLayers - 1) {
	    aResv[k][i] -= 0.5 * ( w[k+1] * 0.5 * ( solv[k][i] - solv[k+1][i] ) / (0.5 * ( solh[k][0] + solh[k+1][0] ) ) );
	  }
	  //aResv[k][i] += sign * 9.81  * rho1[k] / 1024. * zMid[k] / dx;
        }
      }
      if(!bdcHT[k][0]) {
        for (unsigned j = 0; j < nDofv; j++) {
          double sign = ( j == 0) ? 1. : -1;
          aResHT[k][0] += sign * solv[k][j] * solHT[k][0] / dx; 
	}
	if(k<NLayers-1){
	  aResHT[k][0] += w[k+1] * 0.5 * (solHT[k][0]/solh[k][0] + solHT[k+1][0]/solh[k+1][0]);
	  //aResHT[k][0] += w[k+1] * 0.5 * (solT[k][0] + solT[k+1][0]);
	}
	if( k > 0){
	  aResHT[k][0] -= w[k] * 0.5 * (solHT[k-1][0]/solh[k-1][0] + solHT[k][0]/solh[k][0] );
	  //aResHT[k][0] -= w[k] * 0.5 * (solT[k-1][0] + solT[k][0]);
	}
      }
    } 

    vector< double > Res(NLayers * nDofs); // local redidual vector

    unsigned counter = 0;
    for(unsigned k = 0; k < NLayers; k++) {
      for(int i = 0; i < nDofh; i++) {
        Res[counter] =  aResh[k][i].value();
        counter++;
      }
      for(int i = 0; i < nDofv; i++) {
        Res[counter] =  aResv[k][i].value();
        counter++;
      }
      for(int i = 0; i < nDofHT; i++) {
        Res[counter] =  aResHT[k][i].value();
        counter++;
      }
    }

    RES->add_vector_blocked(Res, l2GMap);


    for(unsigned k = 0; k < NLayers; k++) {
      // define the dependent variables
      s.dependent(&aResh[k][0], nDofh);
      s.dependent(&aResv[k][0], nDofv);
      s.dependent(&aResHT[k][0], nDofHT);

      // define the independent variables
      s.independent(&solh[k][0], nDofh);
      s.independent(&solv[k][0], nDofv);
      s.independent(&solHT[k][0], nDofHT);
    }

    // get the jacobian matrix (ordered by row major )
    vector < double > Jac(NLayers * nDofs * NLayers * nDofs);
    s.jacobian(&Jac[0], true);

    //store K in the global matrix KK
    KK->add_matrix_blocked(Jac, l2GMap, l2GMap);

    s.clear_independents();
    s.clear_dependents();

  }*/

  RES->close();
  KK->close();

  PetscViewer    viewer;
  PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
  PetscObjectSetName((PetscObject)viewer,"FSI matrix");
  PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
  MatView((static_cast<PetscMatrix*>(KK))->mat(),viewer);
  double a;
  std::cin>>a;
//
  
//  abort();
  MFN mfn;
  Mat A = (static_cast<PetscMatrix*>(KK))->mat();
  FN f, f1, f2, f3 , f4;
  
  //std::cout << "dt = " << dt << " dx = "<< dx << " maxWaveSpeed = "<<maxWaveSpeed << std::endl;
  std::cout << "dt = " << dt << std::endl;
  
  //dt = 100.;

  Vec v = (static_cast< PetscVector* >(RES))->vec();
  Vec y = (static_cast< PetscVector* >(EPS))->vec();

  MFNCreate( PETSC_COMM_WORLD, &mfn );

  MFNSetOperator( mfn, A );
  MFNGetFN( mfn, &f );


//   FNCreate(PETSC_COMM_WORLD, &f1);
//   FNCreate(PETSC_COMM_WORLD, &f2);
//   FNCreate(PETSC_COMM_WORLD, &f3);
//   FNCreate(PETSC_COMM_WORLD, &f4);
// 
//   FNSetType(f1, FNEXP);
// 
//   FNSetType(f2, FNRATIONAL);
//   double coeff1[1] = { -1};
//   FNRationalSetNumerator(f2, 1, coeff1);
//   FNRationalSetDenominator(f2, 0, PETSC_NULL);
// 
//   FNSetType( f3, FNCOMBINE );
// 
//   FNCombineSetChildren(f3, FN_COMBINE_ADD, f1, f2);
// 
//   FNSetType(f4, FNRATIONAL);
//   double coeff2[2] = {1., 0.};
//   FNRationalSetNumerator(f4, 2, coeff2);
//   FNRationalSetDenominator(f4, 0, PETSC_NULL);
// 
//   FNSetType( f, FNCOMBINE );
// 
//   FNCombineSetChildren(f, FN_COMBINE_DIVIDE, f3, f4);

  FNPhiSetIndex(f,1);
  FNSetType( f, FNPHI );
// FNView(f,PETSC_VIEWER_STDOUT_WORLD);

  FNSetScale( f, dt, dt);
  MFNSetFromOptions( mfn );

  MFNSolve( mfn, v, y);
  MFNDestroy( &mfn );

//   FNDestroy(&f1);
//   FNDestroy(&f2);
//   FNDestroy(&f3);
//   FNDestroy(&f4);


  sol->UpdateSol(mlPdeSys->GetSolPdeIndex(), EPS, pdeSys->KKoffset); 

  unsigned solIndexeta = mlSol->GetIndex("eta");
  unsigned solIndexb = mlSol->GetIndex("b");
  sol->_Sol[solIndexeta]->zero();
  for(unsigned k=0;k<NumberOfLayers;k++){
    sol->_Sol[solIndexeta]->add(*sol->_Sol[solIndexh[k]]);
  }
  sol->_Sol[solIndexeta]->add(-1,*sol->_Sol[solIndexb]);

  for(unsigned k=0; k<NumberOfLayers; k++){
    for(unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++){
      double valueHT = (*sol->_Sol[solIndexHT[k]])(i);
      double valueH = (*sol->_Sol[solIndexh[k]])(i);
            
      double valueT = valueHT/valueH;
    
      sol->_Sol[solIndexT[k]]->set(i, valueT);
    }
    
    sol->_Sol[solIndexT[k]]->close();
    
  }
  
  
}





