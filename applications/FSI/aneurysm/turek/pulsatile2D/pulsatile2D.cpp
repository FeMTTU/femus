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
#include "../../include/FSITimeDependentAssembly.hpp"

double scale = 1000.;

using namespace std;
using namespace femus;

double SetVariableTimeStep(const double time);

bool SetBoundaryConditionTurek2D(const std::vector < double >& x, const char name[],
                                 double &value, const int facename, const double time);

void GetSolutionNorm(MultiLevelSolution& mlSol, const unsigned & group);
//------------------------------------------------------------------------------------------------------------------

int main(int argc, char **args)
{

  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  //Files files;
  //files.CheckIODirectories();
  //files.RedirectCout();

  // ******* Extract the problem dimension and simulation identifier based on the inline input *******


  // ******* Extract the mesh.neu file name based on the simulation identifier *******
//   std::string infile = "./input/aneurysm_Sara_5.neu";
  //std::string infile = "./input/Turek_porous_60micron.neu";
  //std::string infile = "./input/Turek_stents_60micron.neu";
  //std::string infile = "./input/Turek_11stents_60micron.neu";
  std::string infile = "./input/Turek.neu";

  // ******* Set physics parameters *******
  double Lref, Uref, rhof, muf, rhos, ni, E;

  Lref = 1.;
  Uref = 1.;

  rhof = 1035.;
  muf = 3.5 * 1.0e-3; //wrong=3.38*1.0e-4*rhof, note:3.38*1.0e-6*rhof=3.5*1.0e-3
  rhos = 1120;
  ni = 0.5;
  E = 120000 * 1.e0; //turek:120000*1.e6;

  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solid;
  solid = Solid(par, E, ni, rhos, "Mooney-Rivlin");

  cout << "Solid properties: " << endl;
  cout << solid << endl;

  // Generate Fluid Object
  Fluid fluid(par, muf, rhof, "Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  // ******* Init multilevel mesh from mesh.neu file *******
  unsigned short numberOfUniformRefinedMeshes, numberOfAMRLevels;

  numberOfUniformRefinedMeshes = 3;
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
  ml_sol.AddSolution("DX", LAGRANGE, SECOND, 2);
  ml_sol.AddSolution("DY", LAGRANGE, SECOND, 2);

  ml_sol.AddSolution("U", LAGRANGE, SECOND, 2);
  ml_sol.AddSolution("V", LAGRANGE, SECOND, 2);

  // Pair each velocity variable with the corresponding displacement variable
  ml_sol.PairSolution("U", "DX"); // Add this line
  ml_sol.PairSolution("V", "DY"); // Add this line

  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P", DISCONTINOUS_POLYNOMIAL, FIRST, 2);
  ml_sol.AssociatePropertyToSolution("P", "Pressure", false); // Add this line

  // ******* Initialize solution *******
  ml_sol.Initialize("All");

  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionTurek2D);

  // ******* Set boundary conditions *******
  ml_sol.GenerateBdc("DX", "Steady");
  ml_sol.GenerateBdc("DY", "Steady");

  ml_sol.GenerateBdc("U", "Time_dependent");
  ml_sol.GenerateBdc("V", "Steady");

  ml_sol.GenerateBdc("P", "Steady");


  // ******* Define the FSI Multilevel Problem *******

  MultiLevelProblem ml_prob(&ml_sol);
  // Add fluid object
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  // Add Solid Object
  ml_prob.parameters.set<Solid>("Solid") = solid;

  // ******* Add FSI system to the MultiLevel problem *******
  TransientMonolithicFSINonlinearImplicitSystem & system = ml_prob.add_system<TransientMonolithicFSINonlinearImplicitSystem> ("Fluid-Structure-Interaction");
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");

  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  system.AddSolutionToSystemPDE("P");

  // ******* System Fluid-Structure-Interaction Assembly *******
  system.SetAssembleFunction(FSITimeDependentAssembly);

  // ******* set MG-Solver *******
  system.SetMgType(F_CYCLE);

  system.SetNonLinearConvergenceTolerance(1.e-9);
  system.SetResidualUpdateConvergenceTolerance(1.e-15);
  system.SetMaxNumberOfNonLinearIterations(4);
  system.SetMaxNumberOfResidualUpdatesForNonlinearIteration(4);

  system.SetNumberPreSmoothingStep(0);
  system.SetNumberPostSmoothingStep(2);

  // ******* Set Preconditioner *******

  system.SetMgSmoother(ASM_SMOOTHER);

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(RICHARDSON);
  //system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-12, 1.e-20, 1.e+50, 20, 10);

  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  system.SetNumberOfSchurVariables(1);

  // ******* Set block size for the ASM smoothers *******
  system.SetElementBlockNumber(2);

  // ******* Print solution *******
  ml_sol.SetWriter(VTK);


  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  //mov_vars.push_back("DZ");
  ml_sol.GetWriter()->SetMovingMesh(mov_vars);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  ml_sol.GetWriter()->SetDebugOutput(true);
  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);


  // ******* Solve *******
  std::cout << std::endl;
  std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;
  
   // time loop parameter
  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
  const unsigned int n_timesteps = 500;
  
  for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {

    if( time_step > 0 )
      system.SetMgType(V_CYCLE);
    system.CopySolutionToOldSolution();
    system.MGsolve();
    GetSolutionNorm(ml_sol, 9);
    ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR,"biquadratic",print_vars, time_step+1);
  }
  
  // ******* Clear all systems *******
  ml_prob.clear();
  return 0;
}

double SetVariableTimeStep(const double time) {
  double dt = 1./32;
//   if( turek_FSI == 2 ){
//     if ( time < 9 ) dt = 0.05;
//     else dt = 0.025;
//   }
//   else if ( turek_FSI == 3 ){
//     //if	    ( time < 5. ) dt = 0.1;
//     //else
//     if ( time < 6. ) dt = 0.01;
//     else             dt = 0.01;
//   }
//   else if ( simulation == 3 ) dt=0.001;
//   else if ( simulation == 4 ) dt=0.1;
//   else if ( simulation == 5 ) dt=0.1;
//   else if ( simulation == 6 ) dt=0.1;
//   else if ( simulation == 7 ) dt=0.001;
//   else{
//     std::cout << "Warning this simulation case has not been considered yet for the time dependent case"<<std::endl;
//     abort();
//   }
  return dt;
}


//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionTurek2D(const std::vector < double >& x, const char name[], double &value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);
  if ( !strcmp(name, "U") ) {

    if (1 == facename) {
      double ramp = (time < 1) ? sin(PI/2 * time) : 1.;
      value = 0.05 * (x[1] * 1000 - 6) * ( x[1] * 1000 - 8)*(1.+ 0.75*sin(2.*PI*time)) * ramp; //inflow
    }
    else if ( 2 == facename || 5 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if ( !strcmp(name, "V") ) {
    if ( 2 == facename || 5 == facename ) {
      test = 0;
      value = 0.;
    }
  }
  else if (!strcmp(name, "P")) {
    test = 0;
    value = 0.;
  }
  else if (!strcmp(name, "DX") ) {
    //if(2 == facename || 4 == facename || 5 == facename || 6 == facename) {
    if (5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if (!strcmp(name, "DY") ) {
    //if(1 == facename || 3 == facename || 5 == facename || 6 == facename) {
    if ( 5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;

}


void GetSolutionNorm(MultiLevelSolution& mlSol, const unsigned & group)
{

  int  iproc, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  NumericVector* p2;
  NumericVector* v2;
  NumericVector* vol;
  NumericVector* vol0;
  p2 = NumericVector::build().release();
  v2 = NumericVector::build().release();
  vol = NumericVector::build().release();
  vol0 = NumericVector::build().release();

  if (nprocs == 1) {
    p2->init(nprocs, 1, false, SERIAL);
    v2->init(nprocs, 1, false, SERIAL);
    vol->init(nprocs, 1, false, SERIAL);
    vol0->init(nprocs, 1, false, SERIAL);
  }
  else {
    p2->init(nprocs, 1, false, PARALLEL);
    v2->init(nprocs, 1, false, PARALLEL);
    vol->init(nprocs, 1, false, PARALLEL);
    vol0->init(nprocs, 1, false, PARALLEL);
  }

  p2->zero();
  v2->zero();
  vol->zero();
  vol0->zero();

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);


  const unsigned dim = msh->GetDimension();


  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  vector< double > solP;
  vector< vector < double> >  solV(dim);
  vector< vector < double> > x0(dim);
  vector< vector < double> > x(dim);

  solP.reserve(max_size);
  for (unsigned d = 0; d < dim; d++) {
    solV[d].reserve(max_size);
    x0[d].reserve(max_size);
    x[d].reserve(max_size);
  }
  double weight;
  double weight0;

  vector <double> phiV;
  vector <double> gradphiV;
  vector <double> nablaphiV;

  double *phiP;

  phiV.reserve(max_size);
  gradphiV.reserve(max_size * dim);
  nablaphiV.reserve(max_size * (3 * (dim - 1) + !(dim - 1)));

  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol.GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol.GetIndex("V");    // get the position of "V" in the ml_sol object
  if (dim == 3) solVIndex[2] = mlSol.GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol.GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  vector < unsigned > solDIndex(dim);
  solDIndex[0] = mlSol.GetIndex("DX");    // get the position of "U" in the ml_sol object
  solDIndex[1] = mlSol.GetIndex("DY");    // get the position of "V" in the ml_sol object
  if (dim == 3) solDIndex[2] = mlSol.GetIndex("DZ");      // get the position of "V" in the ml_sol object

  unsigned solDType = mlSol.GetSolutionType(solDIndex[0]);  
     
  unsigned solPIndex;
  solPIndex = mlSol.GetIndex("P");
  unsigned solPType = mlSol.GetSolutionType(solPIndex);

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    if ( msh->GetElementGroup(iel) == group ) {
      short unsigned ielt = msh->GetElementType(iel);
      unsigned ndofV = msh->GetElementDofNumber(iel, solVType);
      unsigned ndofP = msh->GetElementDofNumber(iel, solPType);
      unsigned ndofD = msh->GetElementDofNumber(iel, solDType);
      // resize

      phiV.resize(ndofV);
      gradphiV.resize(ndofV * dim);
      nablaphiV.resize(ndofV * (3 * (dim - 1) + !(dim - 1)));

      solP.resize(ndofP);
      for (int d = 0; d < dim; d++) {
        solV[d].resize(ndofV);
	x0[d].resize(ndofD);
        x[d].resize(ndofD);
      }
      // get local to global mappings
      for (unsigned i = 0; i < ndofD; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solDType);
        for (unsigned d = 0; d < dim; d++) {
	  x0[d][i] = (*msh->_topology->_Sol[d])(idof);
	  
          x[d][i] = (*msh->_topology->_Sol[d])(idof) +
		    (*solution->_Sol[solDIndex[d]])(idof);
        }
      }
      
      for (unsigned i = 0; i < ndofV; i++) {
	unsigned idof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof
	for (unsigned  d = 0; d < dim; d++) {
	  solV[d][i] = (*solution->_Sol[solVIndex[d]])(idof);      // global extraction and local storage for the solution
	}
      }
      
     

      for (unsigned i = 0; i < ndofP; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solPType);
        solP[i] = (*solution->_Sol[solPIndex])(idof);
      }


      for (unsigned ig = 0; ig < mlSol._mlMesh->_finiteElement[ielt][solVType]->GetGaussPointNumber(); ig++) {
        // *** get Jacobian and test function and test function derivatives ***
	msh->_finiteElement[ielt][solVType]->Jacobian(x0, ig, weight0, phiV, gradphiV, nablaphiV);
        msh->_finiteElement[ielt][solVType]->Jacobian(x, ig, weight, phiV, gradphiV, nablaphiV);
        phiP = msh->_finiteElement[ielt][solPType]->GetPhi(ig);

	vol0->add(iproc, weight0);
        vol->add(iproc, weight);
      
        std::vector < double> SolV2(dim, 0.);
        for (unsigned i = 0; i < ndofV; i++) {
          for (unsigned d = 0; d < dim; d++) {
            SolV2[d] += solV[d][i] * phiV[i];
          }
        }

        double V2 = 0.;
        for (unsigned d = 0; d < dim; d++) {
          V2 += SolV2[d] * SolV2[d];
        }
        v2->add(iproc, V2 * weight);

        double P2 = 0;
        for (unsigned i = 0; i < ndofP; i++) {
          P2 += solP[i] * phiP[i];
        }
        P2 *= P2;
        p2->add(iproc, P2 * weight);
      }
    }
  }

  p2->close();
  v2->close();
  vol0->close();
  vol->close();

  double p2_l2 = p2->l1_norm();
  double v2_l2 = v2->l1_norm();
  double VOL0 = vol0->l1_norm();
  double VOL = vol->l1_norm();

  std::cout.precision(14);
  std::scientific;
  std::cout << " vol0 = " << VOL0 << std::endl;
  std::cout << " vol = " << VOL << std::endl;
  std::cout << " (vol-vol0)/vol0 = " << (VOL-VOL0) / VOL0 << std::endl;
  std::cout << " p_l2 norm / sqrt(vol) = " << sqrt(p2_l2/VOL)  << std::endl;
  std::cout << " v_l2 norm / sqrt(vol) = " << sqrt(v2_l2/VOL)  << std::endl;

  delete p2;
  delete v2;
  delete vol;

}