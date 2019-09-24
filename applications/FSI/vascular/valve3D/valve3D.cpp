#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "FElemTypeEnum.hpp"
#include "Files.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "VTKWriter.hpp"
#include "../include/FSITimeDependentAssemblySupgNonConservative.hpp"
#include <cmath>
double scale = 1000.;

using namespace std;
using namespace femus;

double SetVariableTimeStep(const double time);

bool SetBoundaryConditionVeinValve(const std::vector < double >& x, const char name[],
                                   double &value, const int facename, const double time);

bool SetBoundaryConditionVeinValve2(const std::vector < double >& x, const char name[],
                                   double &value, const int facename, const double time);

void GetSolutionFluxes(MultiLevelSolution& mlSol, std::vector <double> &fluxes);

// void StoreOldDispcacement(MultiLevelSolution& mlSol);
  
//------------------------------------------------------------------------------------------------------------------

int main(int argc, char **args)
{

  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  //Files files;
  //files.CheckIODirectories();
  //files.RedirectCout();

  bool dimension2D = false;

  // ******* Extract the problem dimension and simulation identifier based on the inline input *******

  //std::string infile = "./input/valve_coarsemesh.neu";
  std::string infile = "./input/mindcraft_valve_1.neu";

  // ******* Set physics parameters *******
  double Lref, Uref, rhof, muf, rhos, ni, E, E1;

  Lref = 1.;
  Uref = 1.;

  rhof = 1060.;
  muf = 2.2 * 1.0e-3;
  rhos = 960;
  ni = 0.5;
  //E = 3.3 * 1.0e6; //vein young modulus
  E = 1.0 * 1.0e6;
  E1 = 0.2 * 1.0e6; //leaflet young modulus

  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solid;
  solid = Solid(par, E, ni, rhos, "Mooney-Rivlin");

  Solid solid1;
  solid1 = Solid(par, E1, ni, rhos, "Mooney-Rivlin");

  cout << "Solid properties: " << endl;
  cout << solid << endl;

  // Generate Fluid Object
  Fluid fluid(par, muf, rhof, "Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  // ******* Init multilevel mesh from mesh.neu file *******
  unsigned short numberOfUniformRefinedMeshes, numberOfAMRLevels;

  numberOfUniformRefinedMeshes = 1;
  numberOfAMRLevels = 0;

  std::cout << 0 << std::endl;

  MultiLevelMesh ml_msh(numberOfUniformRefinedMeshes + numberOfAMRLevels, numberOfUniformRefinedMeshes,
                        infile.c_str(), "fifth", Lref, NULL);

  //ml_msh.EraseCoarseLevels(numberOfUniformRefinedMeshes - 2);

  ml_msh.PrintInfo();

  // mark Solid nodes

  //ml_msh.MarkStructureNode();

  // ******* Init multilevel solution ******
  MultiLevelSolution ml_sol(&ml_msh);

  // ******* Add solution variables to multilevel solution and pair them *******

  ml_sol.AddSolution ( "DX", LAGRANGE, SECOND, 2 );
  ml_sol.AddSolution ( "DY", LAGRANGE, SECOND, 2 );
  if ( !dimension2D ) ml_sol.AddSolution ( "DZ", LAGRANGE, SECOND, 2 );

  ml_sol.AddSolution ( "U", LAGRANGE, SECOND, 2 );
  ml_sol.AddSolution ( "V", LAGRANGE, SECOND, 2 );
  if ( !dimension2D ) ml_sol.AddSolution ( "W", LAGRANGE, SECOND, 2 );

  // Pair each velocity variable with the corresponding displacement variable
  ml_sol.PairSolution ( "U", "DX" ); // Add this line
  ml_sol.PairSolution ( "V", "DY" ); // Add this line
  if ( !dimension2D ) ml_sol.PairSolution ( "W", "DZ" ); // Add this line

  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution ( "P", DISCONTINUOUS_POLYNOMIAL, FIRST, 2 );
  ml_sol.AssociatePropertyToSolution ( "P", "Pressure", false ); // Add this line

  ml_sol.AddSolution ( "lmbd", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false );

  ml_sol.AddSolution ( "Um", LAGRANGE, SECOND, 0, false );
  ml_sol.AddSolution ( "Vm", LAGRANGE, SECOND, 0, false );
  if ( !dimension2D ) ml_sol.AddSolution ( "Wm", LAGRANGE, SECOND, 0, false );

  
  // ******* Initialize solution *******
  ml_sol.Initialize("All");

  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionVeinValve2);

  // ******* Set boundary conditions *******
  ml_sol.GenerateBdc("DX", "Steady");
  ml_sol.GenerateBdc("DY", "Steady");
  if ( !dimension2D ) ml_sol.GenerateBdc("DZ", "Steady");

  ml_sol.GenerateBdc("U", "Steady");
  ml_sol.GenerateBdc("V", "Steady");
  if ( !dimension2D ) ml_sol.GenerateBdc("W", "Steady");

  ml_sol.GenerateBdc("P", "Steady");

  // ******* Define the FSI Multilevel Problem *******

  MultiLevelProblem ml_prob(&ml_sol);
  // Add fluid object
  ml_prob.parameters.set<Fluid> ("Fluid") = fluid;
  // Add Solid Object
  ml_prob.parameters.set<Solid> ("Solid") = solid;
  ml_prob.parameters.set<Solid> ("Solid1") = solid1;

  // ******* Add FSI system to the MultiLevel problem *******
  TransientMonolithicFSINonlinearImplicitSystem & system = ml_prob.add_system<TransientMonolithicFSINonlinearImplicitSystem> ("Fluid-Structure-Interaction");

  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");
  if ( !dimension2D ) system.AddSolutionToSystemPDE("DZ");

  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  if ( !dimension2D ) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P");

  // ******* System Fluid-Structure-Interaction Assembly *******
  system.SetAssembleFunction(FSITimeDependentAssemblySupgNew2);

  // ******* set MG-Solver *******
  system.SetMgType(F_CYCLE);

  system.SetNonLinearConvergenceTolerance(1.e-7);
  //system.SetResidualUpdateConvergenceTolerance ( 1.e-15 );
  system.SetMaxNumberOfNonLinearIterations(20);
  //system.SetMaxNumberOfResidualUpdatesForNonlinearIteration ( 4 );

  system.SetMaxNumberOfLinearIterations(6);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-13);

  system.SetNumberPreSmoothingStep(0);
  system.SetNumberPostSmoothingStep(2);

  // ******* Set Preconditioner *******

  system.SetLinearEquationSolverType(FEMuS_ASM);

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(RICHARDSON);
  //system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(MLU_PRECOND);

  system.SetTolerances(1.e-12, 1.e-20, 1.e+50, 20, 10);

  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  system.SetNumberOfSchurVariables(1);

  // ******* Set block size for the ASM smoothers *******
  system.SetElementBlockNumber(2);



  unsigned time_step_start = 1;

  //char restart_file_name[256] = "./save/valve2D_iteration30";
  char restart_file_name[256] = "";

  if (strcmp (restart_file_name, "") != 0) {
    ml_sol.LoadSolution(restart_file_name);
    time_step_start = 31;
    system.SetTime( (time_step_start - 1) * 1. / 32);
  }

  // ******* Print solution *******
  ml_sol.SetWriter(VTK);


  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  if ( !dimension2D ) mov_vars.push_back("DZ");

  ml_sol.GetWriter()->SetDebugOutput(true);

  ml_sol.GetWriter()->SetMovingMesh(mov_vars);
  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step_start - 1);


  // ******* Solve *******
  std::cout << std::endl;
  std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;

  // time loop parameter
  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
  const unsigned int n_timesteps =1024;
  
  int  iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  
  std::ofstream outf;
  if(iproc == 0) {
    outf.open("fluxes.txt");
    if(!outf) {
      std::cout << "Error in opening file DataPrint.txt";
      return 1;
    }
  }
  
  std::vector < double > Qtot(3,0.);   
  std::vector<double> fluxes(2,0.);

  for (unsigned time_step = time_step_start; time_step <= n_timesteps; time_step++) {

//     StoreOldDispcacement(ml_sol);
    system.CopySolutionToOldSolution();
    
    for (unsigned level = 0; level < numberOfUniformRefinedMeshes; level++) {
      SetLambdaNew(ml_sol, level , SECOND, ELASTICITY);
    }

    if (time_step > 1)
      system.SetMgType(V_CYCLE);


    system.MGsolve();
    
    StoreMeshVelocity(ml_prob);
    
    double dt = system.GetIntervalTime();
    
    Qtot[0] += 0.5 * dt * fluxes[0];
    Qtot[1] += 0.5 * dt * fluxes[1];
    
    GetSolutionFluxes(ml_sol,fluxes);
    
    Qtot[0] += 0.5 * dt * fluxes[0];
    Qtot[1] += 0.5 * dt * fluxes[1];
    Qtot[2] = Qtot[0] + Qtot[1];
    
    std::cout<< fluxes[0] <<" "<<fluxes[1] << Qtot[0] << " " << Qtot[1] << " " << Qtot[2] << std::endl;
    
    if(iproc == 0) {
      outf << time_step <<" "<< system.GetTime() <<" "<< fluxes[0] <<" "<<fluxes[1]<<" " << Qtot[0] << " " << Qtot[1] << " " << Qtot[2] << std::endl;
    }

    ml_sol.GetWriter()->SetMovingMesh(mov_vars);
    ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step);


    if ( time_step % 1 == 0) ml_sol.SaveSolution("valve3D", time_step);

  }

  if(iproc == 0) {
    outf.close();
  }
  
  // ******* Clear all systems *******
  ml_prob.clear();
  return 0;
}

/*
void StoreOldDispcacement(MultiLevelSolution& mlSol){
  
  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  
  //const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();
  const unsigned dim = msh->GetDimension();
  const char varname[6][4] = {"DX", "DY", "DZ", "DX2", "DY2", "DZ2"};

  vector <unsigned> indexVAR(2 * dim);
  vector <unsigned> indVAR(2 * dim);
  vector <unsigned> SolType(2 * dim);

  for (unsigned ivar = 0; ivar < 2 * dim; ivar++) {
      indVAR[ivar] = mlSol.GetIndex(&varname[ivar][0]);
  }
  
  for (unsigned ig=0; ig< level; ig++) {
    Solution* solution  = mlSol.GetSolutionLevel(ig);
    solution->CopySolutionToOldSolution();
    
    for(unsigned ivar = 0; ivar < dim; ivar++) {
    // Copy the old vector
      *(solution->_Sol[indVAR[dim + ivar]]) = *(solution->_SolOld[indVAR[ivar]]);
    }
  }
  
}*/


double SetVariableTimeStep(const double time)
{
  double dt = 1. / 32;
  double shiftedTime = time - floor(time);
  if (time > 1 && shiftedTime >= 0.125 && shiftedTime < 0.25) {
    dt = 1. / 32;
  }
  std::cout << " Shifted Time = " << shiftedTime << " dt = " << dt << std::endl;

  return dt;
}


//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionVeinValve(const std::vector < double >& x, const char name[], double &value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);
  double ramp = (time < 2) ? sin(PI / 2 * time / 2.) : 1.;

  if (!strcmp(name, "U")) {
    if (7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if (!strcmp(name, "V")) {
    if (1 == facename || 2 == facename || 6 == facename || 7 == facename) {
      test = 0;
      value = 0;
    }
  }
  if (!strcmp(name, "W")) {
    if (6 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if (!strcmp(name, "P")) {
    test = 0;
    value = 0.;
    if (1 == facename) {
      value = (0 + 5 * sin(2 * PI * time)) * ramp;      //+ 4.5
    }
    else if (2 == facename) {
      value = (0 - 5 * sin(2 * PI * time)) * ramp;      //- 4.5
    }
  }
  else if ( (!strcmp(name, "DX"))) {
    if (5 == facename || 7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if ( (!strcmp(name, "DY"))) {
    if (5 == facename || 6 == facename || 7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if ( (!strcmp(name, "DZ"))) {
    if (5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;

}

//----------------------------------------------------------------//

bool SetBoundaryConditionVeinValve2(const std::vector < double >& x, const char name[], double &value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);
  double ramp = (time < 2) ? sin(PI / 2 * time / 2.) : 1.;

  if (!strcmp(name, "U")) {
    if (5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if (!strcmp(name, "V")) {
    if (5 == facename || 7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if (!strcmp(name, "W")) {
    if (1 == facename || 2 == facename || 5 == facename || 6 == facename || 7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if (!strcmp(name, "P")) {
    test = 0;
    value = 0.;
    if (1 == facename) {
      value = (0 + 20 * sin(2 * PI * time)) * ramp;      //+ 5
    }
    else if (2 == facename) {
      value = (0 - 20 * sin(2 * PI * time)) * ramp;      //- 5
    }
  }
  else if ( (!strcmp(name, "DX"))) {
    if (5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if ( (!strcmp(name, "DY"))) {
    if (5 == facename || 7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if ( (!strcmp(name, "DZ"))) {
    if (5 == facename || 6 == facename || 7 == facename ) {
      test = 0;
      value = 0;
    }
  }

  return test;

}

//-------------------------------------------------------------------------//

void GetSolutionFluxes(MultiLevelSolution& mlSol, std::vector <double> &fluxes)
{

  int  iproc, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MyVector<double> qTop(1,0);
  qTop.stack();
  
  MyVector<double> qBottom(1,0);
  qBottom.stack();
  
  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);
  elem* myel =  msh->el;
  
  const unsigned dim = msh->GetDimension();
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  vector< vector < double> >  sol(dim);
  vector< vector < double> > x(dim);
 
  const char varname[6][3] = {"U", "V", "W","DX", "DY", "DZ"};
  vector <unsigned> indVar(2 * dim);
  unsigned solType;

  for (unsigned ivar = 0; ivar < dim; ivar++) {
    for (unsigned k = 0; k < 2; k++) {
      indVar[ivar + k * dim] = mlSol.GetIndex(&varname[ivar + k * 3][0]);
    }
  }
  solType = mlSol.GetSolutionType(&varname[0][0]);
    
  
   std::vector < double > phi;
   std::vector < double > gradphi;
   std::vector< double > xx(dim, 0.);
   double weight;
  
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    vector < double> normal(dim, 0);
    
    // loop on faces
    for (unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
      

      int faceNumber = myel->GetBoundaryIndex(iel, jface);
      // look for boundary faces
      if ( faceNumber == 1 || faceNumber ==2) {
       
        unsigned nve = msh->GetElementFaceDofNumber(iel, jface, solType);
        const unsigned felt = msh->GetElementFaceType(iel, jface);
	
	for (unsigned d = 0; d < dim; d++) {
	  x[d].resize(nve);
	  sol[d].resize(nve);
	}
	
        for (unsigned i = 0; i < nve; i++) {
          unsigned int ilocal = msh->GetLocalFaceVertexIndex(iel, jface, i);
          unsigned idof = msh->GetSolutionDof(ilocal, iel, 2);
          for (unsigned d = 0; d < dim; d++) {
            x[d][i] = (*msh->_topology->_Sol[d])(idof) + (*solution->_Sol[indVar[d+dim]])(idof);;
	    sol[d][i] = (*solution->_Sol[indVar[d]])(idof);;
          }
        }

        double flux = 0.;
        for (unsigned igs = 0; igs < msh->_finiteElement[felt][solType]->GetGaussPointNumber(); igs++) {
          msh->_finiteElement[felt][solType]->JacobianSur(x, igs, weight, phi, gradphi, normal);
          double value;
	  for (unsigned i = 0; i < nve; i++) {
	    value = 0.;
	    for (unsigned d = 0; d < dim; d++) {
	      value += normal[d] * sol[d][i];
	    }
	    value *= phi[i];
	  }
	  flux += value * weight;
	}
	if(faceNumber == 1) qBottom[iproc] += flux;
	else qTop[iproc] += flux;
      }     
    }
  }
  
  fluxes[0] = 0.; 
  fluxes[1] = 0.;
  for(int j = 0; j < nprocs; j++) {
    qBottom.broadcast(j);
    qTop.broadcast(j);
    fluxes[0] += qBottom[j]; 
    fluxes[1] += qTop[j]; 
    qBottom.clearBroadcast();
    qTop.clearBroadcast();
  } 
}

