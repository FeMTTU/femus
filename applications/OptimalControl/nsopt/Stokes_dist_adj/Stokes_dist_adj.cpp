// Solving Navier-Stokes problem using automatic differentiation and/or Picards method
// boundary conditions were set in 2D as, no slip in left,right of the box and top to bottom gravity is enforced
// therefore, U=V=0 on left and right, U=0 on top and bottom, V is free 
//distributed control

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
// #include "MultiLevelMesh.hpp"
// #include "LinearImplicitSystem.hpp"
// #include "WriterEnum.hpp"

#include "Fluid.hpp"
#include "Parameter.hpp"

#include <iomanip>



using namespace femus;

 double force[3] = {10.,0.,0.};
 double Vel_desired[3] = {19.,0.,0.};
 double alpha = 1.;
 double beta = 1e-5;


bool SetBoundaryConditionOpt(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  //1: bottom  //2: right  //3: top  //4: left
  
  bool dirichlet = true;
   value = 0.;

//       if (facename == 4) {  //left
// 	  if (!strcmp(SolName, "V"))    { 
// 	      if (x[1] < 0.5 && x[1] > -0.5 ) value = 1.; } 
//       }
// 
//       
// //       if (!strcmp(SolName, "P"))  { value = 0.;      }
  
  
  
  
// LEFT ==========================  
      if (facename == 4) {
       if (!strcmp(SolName, "U"))    { dirichlet = false; }
  else if (!strcmp(SolName, "V"))    {      value = 0.; } 
  else if (!strcmp(SolName, "UADJ")) {  /*value = 0.;*/  dirichlet = false;   }
  else if (!strcmp(SolName, "VADJ")) { value = 0.;   }
	
      }
      
// RIGHT ==========================  
     if (facename == 2) {
       if (!strcmp(SolName, "U"))    {  dirichlet = false; }
  else if (!strcmp(SolName, "V"))    {   value = 0.;  } 
  else if (!strcmp(SolName, "UADJ")) {  /*value = 0.;*/  dirichlet = false; }
  else if (!strcmp(SolName, "VADJ")) {   value = 0.;  } 
  
      }
      
//       if (!strcmp(SolName, "P"))  { 
// 	 dirichlet = false;
//            if (facename == 4)  value = 1.; 
//            if (facename == 2)  value = 0.;
//    
//       }
      
  return dirichlet;
}







// void AssembleNavierStokesOpt(MultiLevelProblem &ml_prob);

void AssembleNavierStokesOpt_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );

// double ComputeIntegral_AD(MultiLevelProblem& ml_prob);

int main(int argc, char** args) {



  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
//   
//   std::string med_file = "RectFracWithGroup.med";
//   std::ostringstream mystream; 
//   mystream << "./" << DEFAULT_INPUTDIR << "/" << med_file;
//   const std::string infile = mystream.str();
//  
   //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;
 // *** apparently needed by non-AD assemble only **********************
  // add fluid material
  Parameter parameter(Lref,Uref);
  
  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,1,1,"Newtonian");
  std::cout << "Fluid properties: " << std::endl;
  std::cout << fluid << std::endl;
  
// *************************

//   MultiLevelMesh mlMsh;
//  mlMsh.ReadCoarseMesh(infile.c_str(),"seventh",Lref);
    mlMsh.GenerateCoarseBoxMesh(8,8,0,-0.5,0.5,-0.5,0.5,0.,0.,QUAD9,"seventh");
    
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 1; 
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  // state =====================  
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);
  mlSol.AddSolution("P", LAGRANGE, FIRST);
  // adjoint =====================  
  mlSol.AddSolution("UADJ", LAGRANGE, SECOND);
  mlSol.AddSolution("VADJ", LAGRANGE, SECOND);
  if (dim == 3) mlSol.AddSolution("WADJ", LAGRANGE, SECOND);
  mlSol.AddSolution("PADJ", LAGRANGE, FIRST);
  // control =====================  

  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionOpt);
  mlSol.GenerateBdc("All");
  

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  mlProb.parameters.set<Fluid>("Fluid") = fluid;

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system_opt    = mlProb.add_system < NonLinearImplicitSystem > ("NSOpt");

  // NS ===================
  system_opt.AddSolutionToSystemPDE("U");
  system_opt.AddSolutionToSystemPDE("V");
  if (dim == 3) system_opt.AddSolutionToSystemPDE("W");
  system_opt.AddSolutionToSystemPDE("P");
//   // NSADJ ===================
  system_opt.AddSolutionToSystemPDE("UADJ");
  system_opt.AddSolutionToSystemPDE("VADJ");
  if (dim == 3) system_opt.AddSolutionToSystemPDE("WADJ");
  system_opt.AddSolutionToSystemPDE("PADJ");
  
  // attach the assembling function to system
  system_opt.SetAssembleFunction(AssembleNavierStokesOpt_AD);
//   system_opt.SetAssembleFunction(AssembleNavierStokesOpt);
    
  // initilaize and solve the system
  system_opt.init();
  system_opt.MLsolve();

//     ComputeIntegral_AD(mlProb);
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  return 0;
}


void AssembleNavierStokesOpt_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NSOpt");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          	= ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         	= msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol    = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        	= ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys  = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    JAC         	= pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  vector < vector < double > > coordX(dim);    // local coordinates
//   vector< vector < double> > coordX_bd(dim);	//local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(maxSize);
//     coordX_bd[k].reserve(maxSize);
  }
  //geometry *******************************

//STATE######################################################################
  //velocity *******************************
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if (dim == 3) solVIndex[2] = mlSol->GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
 vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys->GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys->GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object

  if (dim == 3) solVPdeIndex[2] = mlPdeSys->GetSolPdeIndex("W");
  
  vector < vector < adept::adouble > >  solV(dim);    // local solution
   vector< vector < adept::adouble > > aResV(dim);    // local redidual vector
   
 for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    aResV[k].reserve(maxSize);
  }

  
  vector <double> phiV_gss;  // local test function
  vector <double> phiV_x_gss; // local test function first order partial derivatives
  vector <double> phiV_xx_gss; // local test function second order partial derivatives

  phiV_gss.reserve(maxSize);
  phiV_x_gss.reserve(maxSize * dim);
  phiV_xx_gss.reserve(maxSize * dim2);
  
  
//   vector < double > phiV_gss_bd;
//   vector < double > phiV_x_gss_bd;
//   phiV_gss_bd.reserve(maxSize);
//   phiV_x_gss_bd.reserve(maxSize*dim);
  //velocity *******************************

  //pressure *******************************
  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

  vector < adept::adouble >  solP; // local solution
  vector< adept::adouble > aResP; // local redidual vector
  
  solP.reserve(maxSize);
  aResP.reserve(maxSize);
  
  double* phiP_gss;
  //pressure *******************************
//STATE######################################################################
  
//ADJOINT######################################################################
  //velocity *******************************
  vector < unsigned > solVadjIndex(dim);
  solVadjIndex[0] = mlSol->GetIndex("UADJ");    // get the position of "U" in the ml_sol object
  solVadjIndex[1] = mlSol->GetIndex("VADJ");    // get the position of "V" in the ml_sol object

  if (dim == 3) solVadjIndex[2] = mlSol->GetIndex("WADJ");      // get the position of "V" in the ml_sol object

  unsigned solVadjType = mlSol->GetSolutionType(solVadjIndex[0]);    // get the finite element type for "u"
 vector < unsigned > solVPdeadjIndex(dim);
  solVPdeadjIndex[0] = mlPdeSys->GetSolPdeIndex("UADJ");    // get the position of "U" in the pdeSys object
  solVPdeadjIndex[1] = mlPdeSys->GetSolPdeIndex("VADJ");    // get the position of "V" in the pdeSys object

  if (dim == 3) solVPdeadjIndex[2] = mlPdeSys->GetSolPdeIndex("WADJ");
  
  vector < vector < adept::adouble > >  solVadj(dim);    // local solution
   vector< vector < adept::adouble > > aResVadj(dim);    // local redidual vector
   
 for (unsigned  k = 0; k < dim; k++) {
    solVadj[k].reserve(maxSize);
    aResVadj[k].reserve(maxSize);
  }

  
  vector <double> phiVadj_gss;  // local test function
  vector <double> phiVadj_x_gss; // local test function first order partial derivatives
  vector <double> phiVadj_xx_gss; // local test function second order partial derivatives //NOT USED, ONLY PASSED

  phiVadj_gss.reserve(maxSize);
  phiVadj_x_gss.reserve(maxSize * dim);
  phiVadj_xx_gss.reserve(maxSize * dim2);   // NOT USED, ONLY PASSED
  
  //velocity *******************************

  //pressure *******************************
  unsigned solPadjIndex;
  solPadjIndex = mlSol->GetIndex("PADJ");    // get the position of "P" in the ml_sol object
  unsigned solPadjType = mlSol->GetSolutionType(solPadjIndex);    // get the finite element type for "u"

  unsigned solPPdeadjIndex;
  solPPdeadjIndex = mlPdeSys->GetSolPdeIndex("PADJ");    // get the position of "P" in the pdeSys object

  vector < adept::adouble >  solPadj; // local solution
  vector< adept::adouble > aResPadj; // local redidual vector
  
  solPadj.reserve(maxSize);
  aResPadj.reserve(maxSize);
  
  double* phiPadj_gss;
  //pressure *******************************
//ADJOINT######################################################################
 
 

// //target_vel####################################################################
// 
//   vector <double> phiVTarget_gss;  // local test function
//   vector <double> phiVTarget_x_gss; // local test function first order partial derivatives
//   vector <double> phiVTarget_xx_gss; // local test function second order partial derivatives
// 
//   phiVTarget_gss.reserve(maxSize);
//   phiVTarget_x_gss.reserve(maxSize * dim);
//   phiVTarget_xx_gss.reserve(maxSize * dim2);
//   
//  
//    vector <double> solVTarget;
//    solVTarget.reserve(maxSize);
//    
//    vector <double> VelTarget_gss;
//    VelTarget_gss.reserve(maxSize);
// //target_vel####################################################################

   
   
  
  
  //Nondimensional values ******************
  double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  //Nondimensional values ******************
  
  double weight; // gauss point weight
//   double weight_bd;
  
  
  vector< int > JACDof; // local to global pdeSys dofs
  JACDof.reserve(2 *(dim + 1) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve(2 *(dim + 1) *maxSize);

  vector < double > Jac;
  Jac.reserve(2* (dim + 1) *maxSize * 2*(dim + 1) *maxSize);

  
  double target_vel[3];
  for (unsigned v=0; v < dim; v++) {     target_vel[v]=Vel_desired[v]; }
  
  
  
  
  
  JAC->zero(); // Set to zero all the entries of the Global Matrix

  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

// geometry
    short unsigned ielGeom = msh->GetElementType(iel);    // element geometry type

    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs
    
// equation
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsVP = dim * nDofsV + nDofsP;
    unsigned nDofsVP_tot = 2*nDofsVP;
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel, solVadjType);    // number of solution element dofs
    unsigned nDofsPadj = msh->GetElementDofNumber(iel, solPadjType);    // number of solution element dofs

    for (unsigned  k = 0; k < dim; k++) {       coordX[k].resize(nDofsX);    }
     
    for (unsigned  k = 0; k < dim; k++)  {
      solV[k].resize(nDofsV);
      solVadj[k].resize(nDofsVadj);
    }
    solP.resize(nDofsP);
    solPadj.resize(nDofsPadj);
    

//element matrices and vectors
    // resize local arrays
    JACDof.resize(nDofsVP_tot);
    
    
    

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].resize(nDofsV);    //resize
      std::fill(aResV[k].begin(), aResV[k].end(), 0);    //set aRes to zero
    
      aResVadj[k].resize(nDofsVadj);    //resize
      std::fill(aResVadj[k].begin(), aResVadj[k].end(), 0);    //set aRes to zero

    }

    aResP.resize(nDofsP);    //resize
    std::fill(aResP.begin(), aResP.end(), 0);    //set aRes to zero

     aResPadj.resize(nDofsPadj);    //resize
    std::fill(aResPadj.begin(), aResPadj.end(), 0);    //set aRes to zero

    // geometry ************
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }
    
  //STATE###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        JACDof[i + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    
    // pressure *************
    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // global to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      JACDof[i + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
//STATE###################################################################

//ADJ###################################################################
     // velocity ************
    for (unsigned i = 0; i < nDofsVadj; i++) {
      unsigned solVadjDof = msh->GetSolutionDof(i, iel, solVadjType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solVadj[k][i] = (*sol->_Sol[solVadjIndex[k]])(solVadjDof);      // global extraction and local storage for the solution
        JACDof[i + k * nDofsVadj +nDofsVP] = pdeSys->GetSystemDof(solVadjIndex[k], solVPdeadjIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    
    // pressure *************
    for (unsigned i = 0; i < nDofsPadj; i++) {
      unsigned solPadjDof = msh->GetSolutionDof(i, iel, solPadjType);    // global to global mapping between solution node and solution dof
      solPadj[i] = (*sol->_Sol[solPadjIndex])(solPadjDof);      // global extraction and local storage for the solution
      JACDof[i + dim * nDofsVadj +nDofsVP] = pdeSys->GetSystemDof(solPadjIndex, solPPdeadjIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
//ADJ###################################################################


// //target_vel####################################################################
//       unsigned nDofsVTarget = el->GetElementDofNumber(iel, solVType);    // number of solution element dofs
//       solVTarget.resize(nDofsVTarget);
//       for (unsigned i = 0; i < solVTarget.size(); i++) {
// 	  solVTarget[i] = target_vel[i];
//       }
// 
//       
//       VelTarget_gss.resize(nDofsVTarget);
//       std::fill(VelTarget_gss.begin(),VelTarget_gss.end(),0.);
// //target_vel####################################################################




      // start a new recording of all the operations involving adept::adouble variables
      s.new_recording();

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {

//STATE#############################################################################	
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV_gss, phiV_x_gss, phiV_xx_gss);
        phiP_gss = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

	
        vector < adept::adouble > solV_gss(dim, 0);
        vector < vector < adept::adouble > > gradSolV_gss(dim);

        for (unsigned  k = 0; k < dim; k++) {
          gradSolV_gss[k].resize(dim);
          std::fill(gradSolV_gss[k].begin(), gradSolV_gss[k].end(), 0);
        }

        for (unsigned i = 0; i < nDofsV; i++) {
          for (unsigned  k = 0; k < dim; k++) {
            solV_gss[k] += phiV_gss[i] * solV[k][i];
          }

          for (unsigned j = 0; j < dim; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              gradSolV_gss[k][j] += phiV_x_gss[i * dim + j] * solV[k][i];
            }
          }
        }

        adept::adouble solP_gss = 0;

        for (unsigned i = 0; i < nDofsP; i++) {
          solP_gss += phiP_gss[i] * solP[i];
        }

//STATE###############################################################################

//ADJOINT#############################################################################3	
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][solVadjType]->Jacobian(coordX, ig, weight, phiVadj_gss, phiVadj_x_gss, phiVadj_xx_gss);
        phiPadj_gss = msh->_finiteElement[ielGeom][solPadjType]->GetPhi(ig);

	
        vector < adept::adouble > solVadj_gss(dim, 0);
        vector < vector < adept::adouble > > gradSolVadj_gss(dim);

        for (unsigned  k = 0; k < dim; k++) {
          gradSolVadj_gss[k].resize(dim);
          std::fill(gradSolVadj_gss[k].begin(), gradSolVadj_gss[k].end(), 0);
        }

        for (unsigned i = 0; i < nDofsVadj; i++) {
          for (unsigned  k = 0; k < dim; k++) {
            solVadj_gss[k] += phiVadj_gss[i] * solVadj[k][i];
          }

          for (unsigned j = 0; j < dim; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              gradSolVadj_gss[k][j] += phiVadj_x_gss[i * dim + j] * solVadj[k][i];
            }
          }
        }

        adept::adouble solPadj_gss = 0;

        for (unsigned i = 0; i < nDofsPadj; i++) {
          solPadj_gss += phiPadj_gss[i] * solPadj[i];
        }

//ADJOINT###############################################################################


// //target_vel################################################################################
//         msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiVTarget_gss, phiVTarget_x_gss, phiVTarget_xx_gss);
// 
// 	
// 	 for (unsigned i =0; i < nDofsVTarget; i++) {
// 	   for (unsigned  k = 0; k < dim; k++) {
// 	   
// 	  VelTarget_gss[k] += solVTarget[i] * phiVTarget_gss[i];
// 	  }
// 	}
// //target_vel################################################################################
// 	  
	  

        // *** phiV_i loop ***
        for (unsigned i = 0; i < nDofsV; i++) {
          vector < adept::adouble > NSV_gss(dim, 0.);
	  vector < adept::adouble > NSVadj_gss(dim, 0.);
	  
          for (unsigned  kdim = 0; kdim < dim; kdim++) { // velocity block row
	      
          for (unsigned jdim = 0; jdim < dim; jdim++) { //focus on single partial derivative
//               NSV_gss[kdim]   +=  phiV_gss[i] * (solV_gss[jdim] * gradSolV_gss[kdim][jdim]);                                  //advection
              NSV_gss[kdim]      +=  IRe*phiV_x_gss   [i * dim + jdim]*gradSolV_gss   [kdim][jdim];/*deformation_tensor*//*IRe * phiV_x_gss[i * dim + jdim] * (gradSolV_gss[kdim][jdim] + gradSolV_gss[jdim][kdim])*/  //diffusion
	      NSVadj_gss[kdim]   +=  IRe*phiVadj_x_gss[i * dim + jdim]*gradSolVadj_gss[kdim][jdim]; 
 
	  }  //jdim loop
            
            //velocity-pressure block
          NSV_gss[kdim]    += - solP_gss    * phiV_x_gss   [i * dim + kdim];
	  NSVadj_gss[kdim] += - solPadj_gss * phiVadj_x_gss[i * dim + kdim];
          
	    
	  } //kdim loop

          for (unsigned  kdim = 0; kdim < dim; kdim++) {
            aResV[kdim][i]      +=   (force[kdim] * phiV_gss[i]    - NSV_gss[kdim]   ) * weight;
	    aResVadj[kdim][i]   += (alpha*(solV_gss[kdim]   - /*VelTarget_gss*/target_vel[kdim]) * phiVadj_gss[i] - NSVadj_gss[kdim]) * weight;  /*- solVadj_gss[kdim]*phiVadj_gss[i]*weight; //*/ 
          }
        } // end phiV_i loop

        // *** phiP_i loop ***
        for (unsigned i = 0; i < nDofsP; i++) {
          for (int kdim = 0; kdim < dim; kdim++) {
            aResP[i]    +=   - (gradSolV_gss[kdim][kdim]) * phiP_gss[i]  * weight;
	    aResPadj[i] += - (gradSolVadj_gss[kdim][kdim]) * phiPadj_gss[i]  * weight;    /*-  solPadj_gss*phiPadj_gss[i]*weight;   //*/
          }
        } // end phiP_i loop

      } // end gauss point loop
      
              



    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsVP_tot);    //resize
    std::fill(Res.begin(), Res.end(), 0.);

    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
        Res[ i +  kdim * nDofsV ] = - aResV[kdim][i].value();
	Res[ i +  kdim * nDofsV + nDofsVP] =  - aResVadj[kdim][i].value();
      }
    }

    for (int i = 0; i < nDofsP; i++) {
      Res[ i + dim * nDofsV ] = - aResP[i].value();
      Res[ i + dim * nDofsV + nDofsVP] = - aResPadj[i].value();
    }

    
//     //print Res vector 
//      for (unsigned  l = 0; l < nDofsVP_tot; l++) {
// 	    std::cout << " " << std::setfill(' ') << std::setw(10) << Res[l] << std::endl;
// 	      }

    RES->add_vector_blocked(Res, JACDof);

    //Extarct and store the Jacobian
       Jac.resize(nDofsVP_tot * nDofsVP_tot);
       std::fill(Jac.begin(), Jac.end(), 0.);
       
//        //print empty Jac matrix 
//        for (unsigned  l = 0; l < nDofsVP_tot; l++) {
//             for (unsigned  p = 0; p < nDofsVP_tot; p++) {
//             std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[p+ l*nDofsVP_tot];
// 	    }
// 	    std::cout << std::endl;
// 	  }
       
       
      // define the dependent variables
       	//state_aRes##############################
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
        s.dependent(&aResV[kdim][0], nDofsV);
      }
	s.dependent(&aResP[0], nDofsP);
      
	//adjoint_aRes###############################
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
	 s.dependent(&aResVadj[kdim][0], nDofsVadj);
        }
	 s.dependent(&aResPadj[0], nDofsPadj);

    
	 
	 
	 

      // define the independent variables
	//state_solVAR############################## 
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
        s.independent(&solV[kdim][0], nDofsV);
      }

	s.independent(&solP[0], nDofsP);
     
	//adjoint_solVAR##############################
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
	s.independent(&solVadj[kdim][0], nDofsVadj);
      }
	s.independent(&solPadj[0], nDofsPadj);

	
	
      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0] , true);
      
      
//       //print Jac matrix 
//           for (unsigned  l = 0; l < nDofsVP_tot; l++) {
//             for (unsigned  p = 0; p < nDofsVP_tot; p++) {
//             std::cout << " "  << std::setfill(' ') << std::setw(15) << Jac[p+ l*nDofsVP_tot] << " ";
// 	    }
// 	    std::cout << std::endl;
// 	  }

      
      JAC->add_matrix_blocked(Jac, JACDof, JACDof);

      s.clear_independents();
      s.clear_dependents();
    
  } //end element loop for each process

  RES->close();

  JAC->close();

  // ***************** END ASSEMBLY *******************
}


double ComputeIntegral_AD(MultiLevelProblem& ml_prob) {
    adept::Stack& s = FemusInit::_adeptStack;

   NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NSOpt");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          	= ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         	= msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol    = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        	= ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys  = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  vector < vector < double > > coordX(dim);    // local coordinates

  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(maxSize);
  }
  
  double weight;
  
  
  //geometry *******************************

//STATE######################################################################
  //velocity *******************************
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if (dim == 3) solVIndex[2] = mlSol->GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
  
  vector < vector < adept::adouble > >  solV(dim);    // local solution
// //   vector < vector < adept::adouble > >  V_gss(dim);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
//     V_gss[k] = 0;
  }

  
  vector <double> phiV_gss;  // local test function
  vector <double> phiV_x_gss; // local test function first order partial derivatives
  vector <double> phiV_xx_gss; // local test function second order partial derivatives

  phiV_gss.reserve(maxSize);
  phiV_x_gss.reserve(maxSize * dim);
  phiV_xx_gss.reserve(maxSize * dim2);
  
  
  //velocity *******************************
   

//STATE######################################################################
  

// Vel_desired##################################################################
  vector <double> phiVdes_gss;  // local test function
  vector <double> phiVdes_x_gss; // local test function first order partial derivatives
  vector <double> phiVdes_xx_gss; // local test function second order partial derivatives

  phiVdes_gss.reserve(maxSize);
  phiVdes_x_gss.reserve(maxSize * dim);
  phiVdes_xx_gss.reserve(maxSize * dim2);

  vector< vector < adept::adouble > >  solVdes(dim);    // local solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solVdes[k].reserve(maxSize);
  }
//   double* Vdes_gss [3] = Vel_desired/*= 0.*/;


// Vel_desired##################################################################



// vector<adept::adouble> integralval;
vector<adept::adouble> integral(dim);

double integralval = 0.;
//    JAC->zero(); // Set to zero all the entries of the Global Matrix

  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

// geometry
    short unsigned ielGeom = msh->GetElementType(iel);    // element geometry type

    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs
    
// equation
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
//     unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
//     unsigned nDofsVP = dim * nDofsV + nDofsP;
//     unsigned nDofsVP_tot = 2*nDofsVP;

    unsigned nDofsVdes = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs

    
//     unsigned nDofsVadj = msh->GetElementDofNumber(iel, solVadjType);    // number of solution element dofs
//     unsigned nDofsPadj = msh->GetElementDofNumber(iel, solPadjType);    // number of solution element dofs

    for (unsigned  k = 0; k < dim; k++) {       coordX[k].resize(nDofsX);    }
     
    for (unsigned  k = 0; k < dim; k++)  {
      solV[k].resize(nDofsV);
      solVdes[k].resize(nDofsVdes);
    }


    // geometry ************
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }
  //STATE###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
//         JACDof[i + k * nDofsV] = pdeSys->GetKKDof(solVIndex[k], solVPdeIndex[k], iNode);    // global to global mapping between solution node and pdeSys dof
      }
    }
    

//STATE###################################################################




  //DESIRED VEL###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVdesDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < solVdes.size() /*dim*/; k++) {
        solVdes[k][i] = Vel_desired[k] /*(*sol->_Sol[solVIndex[k]])(solVdesDof)*/;      // global extraction and local storage for the solution
//         JACDof[i + k * nDofsV] = pdeSys->GetKKDof(solVIndex[k], solVPdeIndex[k], iNode);    // global to global mapping between solution node and pdeSys dof
      }
    }
    

//DESIRED VEL###################################################################


      // start a new recording of all the operations involving adept::adouble variables
      s.new_recording();

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {

//STATE#############################################################################3	
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV_gss, phiV_x_gss, phiV_xx_gss);
//         phiP_gss = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);

	msh->_finiteElement[ielGeom][solVType  /*solVdes*/]->Jacobian(coordX, ig, weight, phiVdes_gss, phiVdes_x_gss, phiVdes_xx_gss);

	vector<adept::adouble> V_gss(dim) ;
	vector<adept::adouble> Vdes_gss(dim);

	 for (unsigned  k = 0; k < dim; k++) { 
	   V_gss[k]=0.;
	   Vdes_gss[k] = 0.;
	}
	 
	for (unsigned  k = 0; k < dim; k++) {
	    for (unsigned i = 0; i < nDofsV; i++) {
		V_gss[k] += solV[k][i] * phiV_gss[i];
		Vdes_gss[k] += solVdes[k][i] * phiVdes_gss[i];
		}
	}

      for (unsigned  k = 0; k < dim; k++) {
	  integral[k] += (V_gss[k] - Vdes_gss[k]) * (V_gss[k] - Vdes_gss[k]);
      }

//       integralval= sqrt((integral[0]*integral[0]) +(integral[1]*integral[1]))*weight;
// // 
    std::cout << "The value of the integral is " << std::setw(11) << std::setprecision(10) << ((alpha/2)*10.+ (beta/2)*1.)*weight/*integralval*/ << std::endl;

    return /*integralval*/ 10.*weight; 
	    }
    
  }
}

// nonAD is in the old PETSc, edit this for the new PETSc
// void AssembleNavierStokesOpt(MultiLevelProblem &ml_prob){
//      
//   //pointers
//   LinearImplicitSystem& mlPdeSys  = ml_prob.get_system<LinearImplicitSystem>("NSOPT");
//   const unsigned level = mlPdeSys.GetLevelToAssemble();
//   const unsigned  levelMax= mlPdeSys.GetLevelMax();
//   bool assembleMatrix = mlPdeSys.GetAssembleMatrix(); 
//    
//   Solution*	 sol  	         = ml_prob._ml_sol->GetSolutionLevel(level);
//   LinearEquationSolver*  pdeSys	 = mlPdeSys._LinSolver[level];   
//   const char* pdename            = mlPdeSys.name().c_str();
//   
//   MultiLevelSolution* mlSol = ml_prob._ml_sol;
//   
//   Mesh*		 msh    = ml_prob._ml_msh->GetLevel(level);
//   elem*		 el	= msh->el;
//   SparseMatrix*	 JAC	= pdeSys->_KK;
//   NumericVector* RES 	= pdeSys->_RES;
//     
//   //data
//   const unsigned dim 	= msh->GetDimension();
//   unsigned nel		= msh->GetNumberOfElements();
//   unsigned igrid	= msh->GetLevel();
//   unsigned iproc 	= msh->processor_id();
//  
//   const unsigned maxSize = static_cast< unsigned > (ceil(pow(3,dim)));
// 
//   // geometry *******************************************
//   vector< vector < double> > coordX(dim);	//local coordinates
//   vector< vector < double> > coordX_bd(dim);	//local coordinates
//   unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
//   for(int i=0;i<dim;i++) {   
//        coordX[i].reserve(maxSize); 
//     coordX_bd[i].reserve(maxSize); 
//   }
//   double normal_bd[3] = {0.,0.,0.};
//   // geometry *******************************************
// 
//   vector < vector < int > > node_pos_sol(NFE_FAMS); 
//   for(int i=0; i < NFE_FAMS; i++) { node_pos_sol[i].reserve(maxSize); }   
// //   vector < int > node_pos_sol_u; //2 
// //   vector < int > node_pos_sol_p; //0
//     // reserve
// //   node_pos_sol_u.reserve(maxSize);
// //   node_pos_sol_p.reserve( static_cast< unsigned > (ceil(pow(2,dim))));
//   
//   
//   // solution variables *******************************************
//   const int n_vars = dim+1;
//   const int n_unknowns = (2.*dim)+1;		//state , adjoint of velocity terms and one pressure term
//   const int vel_type_pos = 0;
//   const int adj_vel_type_pos = vel_type_pos;
//   const int press_type_pos = dim;
//   const int state_pos_begin = 0;
//   const int adj_pos_begin   = dim+1;
//   vector < std::string > Solname(n_unknowns);  // const char Solname[4][8] = {"U","V","W","P"};
//   Solname              [state_pos_begin+0] = "U";
//   Solname              [state_pos_begin+1] = "V";
//   if (dim == 3) Solname[state_pos_begin+2] = "W";
//   Solname              [state_pos_begin + press_type_pos] = "P";
//   
//   Solname              [adj_pos_begin + 0] =              "UADJ";
//   Solname              [adj_pos_begin + 1] =              "VADJ";
//   if (dim == 3) Solname[adj_pos_begin + 2] =              "WADJ";
// //   Solname              [adj_pos_begin + press_type_pos] = "PADJ";
//   
//   vector < unsigned > SolPdeIndex(n_unknowns);
//   vector < unsigned > SolIndex(n_unknowns);  
//   vector < unsigned > SolFEType(n_unknowns);  
// 
// 
//   for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
//     SolPdeIndex[ivar]	= mlPdeSys.GetSolPdeIndex(Solname[ivar].c_str());
//     SolIndex[ivar]	= mlSol->GetIndex        (Solname[ivar].c_str());
//     SolFEType[ivar]	= mlSol->GetSolutionType(SolIndex[ivar]);
//   }
// 
//   vector < double > Sol_n_el_dofs(n_unknowns);
//   
//   //==========================================================================================
//   // velocity ************************************
//   //-----------state------------------------------
//   vector < vector < double > > phi_gss_fe(NFE_FAMS);
//   vector < vector < double > > phi_x_gss_fe(NFE_FAMS);
//   vector < vector < double > > phi_xx_gss_fe(NFE_FAMS);
//  
//   for(int fe=0; fe < NFE_FAMS; fe++) {  
//         phi_gss_fe[fe].reserve(maxSize);
//       phi_x_gss_fe[fe].reserve(maxSize*dim);
//      phi_xx_gss_fe[fe].reserve(maxSize*(3*(dim-1)));
//    }
//    
//   vector < double > phiV_gss_bd;
//   vector < double > phiV_x_gss_bd;
//   phiV_gss_bd.reserve(maxSize);
//   phiV_x_gss_bd.reserve(maxSize*dim);
//   
//   //=================================================================================================
//   
//   // quadratures ********************************
//   double weight;
//   double weight_bd;
//   
//   
//   // equation ***********************************
//   vector < vector < int > > JACDof(n_unknowns); 
//   vector < vector < double > > Res(n_unknowns); /*was F*/
//   vector < vector < vector < double > > > Jac(n_unknowns); /*was B*/
//  
//   for(int i = 0; i < n_unknowns; i++) {     
//     JACDof[i].reserve(maxSize);
//       Res[i].reserve(maxSize);
//   }
//    
//   if(assembleMatrix){
//     for(int i = 0; i < n_unknowns; i++) {
//       Jac[i].resize(n_unknowns);    
//       for(int j = 0; j < n_unknowns; j++) {
// 	Jac[i][j].reserve(maxSize*maxSize);	
//       }
//     }
//   }
//   
//   //-----------state------------------------------
//   vector < double > SolVAR(n_unknowns);
//   vector < vector < double > > gradSolVAR(n_unknowns);
//   
//   for(int i=0; i<n_unknowns; i++) {     
//     gradSolVAR[i].resize(dim);    
//   }
//   
// 
//   double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
//   
//   // SUPG - not needed *****************************
// //   double ILambda	= 0; 
// //   bool penalty 		= true; 
//   
// //   double alpha = 0.;
// //   if(solPType == solVType && solVType == 0) // if pressure and velocity are both linear, we need stabilization 
// //   {
// //     alpha = 0.013333; 
// //   }
//   // SUPG - not needed *****************************
//  
//   
//   // Set to zeto all the entries of the matrix
//   if(assembleMatrix) JAC->zero();
//   
//   // ****************** element loop *******************
//  
//   for (int iel=msh->IS_Mts2Gmt_elem_offset[iproc]; iel < msh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {
// 
//   // geometry *****************************
//     unsigned kel = msh->IS_Mts2Gmt_elem[iel];
//     short unsigned kelGeom = el->GetElementType(kel);
//      unsigned nDofsX = el->GetElementDofNumber(kel, coordXType);    // number of coordinate element dofs
// 
//     for(int ivar=0; ivar<dim; ivar++) {
//       coordX[ivar].resize(nDofsX);
//     }
//    for( unsigned i=0;i<nDofsX;i++) {
//       unsigned inode = el->GetElementVertexIndex(kel,i)-1u;
//       unsigned inode_coord_metis = msh->GetMetisDof(inode,coordXType);
//       for(unsigned ivar = 0; ivar < dim; ivar++) {
// 	coordX[ivar][i] = (*msh->_coordinate->_Sol[ivar])(inode_coord_metis);
//       }
//     }
//   // geometry end *****************************
// 
//     for(unsigned unk = 0; unk < n_unknowns; unk++) {
//       Sol_n_el_dofs[unk] = el->GetElementDofNumber(kel,SolFEType[unk]);
//        JACDof[unk].resize(Sol_n_el_dofs[unk]);
//     }
// //     unsigned  nDof_st = Sol_n_el_dofs[state_pos_begin];
// //     unsigned nDof_adj = Sol_n_el_dofs[adj_pos_begin];
//     
// // fe data
//   for(int fe=0; fe < NFE_FAMS; fe++) {
//     unsigned n_el_dof_fe = el->GetElementDofNumber(kel,fe);
//     node_pos_sol[fe].resize(n_el_dof_fe); 
//     phi_gss_fe[fe].resize(n_el_dof_fe);
//     phi_x_gss_fe[fe].resize(n_el_dof_fe*dim);
//     phi_xx_gss_fe[fe].resize(n_el_dof_fe*(3*(dim-1)));
// 
//        for( unsigned i=0;i< n_el_dof_fe; i++) {
//         unsigned inode; //TODO 
//        if (fe == 2) { inode = el->GetElementVertexIndex(kel,i)-1u;
// 	node_pos_sol[fe][i] =  msh->GetMetisDof(inode,fe);
//           }
//   else if (fe == 0) { 
//     inode = (fe < dim)?(el->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
//     node_pos_sol[fe][i] = inode;
//         }
//       }  
//   }   
// 
//     
// //kkdof
//     for(unsigned unk = 0; unk < n_unknowns; unk++) {
//        for( unsigned i=0;i<Sol_n_el_dofs[unk]; i++) {
//         unsigned inode; //TODO 
//        if (SolFEType[unk] == 2)  { inode = el->GetElementVertexIndex(kel,i)-1u;}
//   else if (SolFEType[unk] == 0)  { inode = (SolFEType[unk] < dim)?(el->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);}
// 
//        JACDof[unk][i] = pdeSys->GetKKDof(SolIndex[unk],SolPdeIndex[unk],inode);
//          }
//       }
//        
//     for(int ivar=0; ivar<n_unknowns; ivar++) {
//       Res[SolPdeIndex[ivar]].resize(Sol_n_el_dofs[ivar]);
//       memset(&Res[SolPdeIndex[ivar]][0],0,Sol_n_el_dofs[ivar]*sizeof(double));
//     }
//    
//     for(int ivar=0; ivar<n_unknowns; ivar++) {
//       for(int jvar=0; jvar<n_unknowns; jvar++) {
//       if(assembleMatrix){  //MISMATCH
// 	Jac[ SolPdeIndex[ivar] ][ SolPdeIndex[jvar] ].resize(Sol_n_el_dofs[ivar]*Sol_n_el_dofs[jvar]);
// 	memset(&Jac[SolPdeIndex[ivar]][SolPdeIndex[jvar]][0],0,Sol_n_el_dofs[ivar]*Sol_n_el_dofs[jvar]*sizeof(double));
//       }
//     }
//   }
//   
// //   std::fill(Res.begin(),Res.end(),0.);
// //   std::fill(Jac.begin(),Jac.end(),0.);
//     //=============================================================================
// 
// // SUPG - not needed
// //     double hk = sqrt( (coordX[0][2] - coordX[0][0])*(coordX[0][2] - coordX[0][0]) + 
// //       (coordX[1][2] - coordX[1][0])*(coordX[1][2] - coordX[1][0]) );
//     
// 
//     
// 
// //     int nDof_max = nDof_st;
//     //add the comparison with nDof_adj!! 
//        
//     if(igrid==levelMax || !el->GetRefinedElementIndex(kel)) {
//       
//       // ********************** Gauss point loop *******************************
//       for(unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelGeom][SolFEType[vel_type_pos]]->GetGaussPointNumber(); ig++) {
// 	
// 	// *** get Jacobian and test function and test function derivatives ***
//       for(int fe=0; fe < NFE_FAMS; fe++) {
// 	ml_prob._ml_msh->_finiteElement[kelGeom][fe]->Jacobian(coordX,ig,weight,phi_gss_fe[fe],phi_x_gss_fe[fe],phi_xx_gss_fe[fe]);
//       }
//          //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
//   	ml_prob._ml_msh->_finiteElement[kelGeom][BIQUADR_FE]->Jacobian(coordX,ig,weight,phi_gss_fe[BIQUADR_FE],phi_x_gss_fe[BIQUADR_FE],phi_xx_gss_fe[BIQUADR_FE]);
// 
// 
//  //begin unknowns eval at gauss points ********************************
// 	for(unsigned unk = 0; unk < /*n_vars*/ n_unknowns; unk++) {
// 	  SolVAR[unk]=0;
// 	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
// 	    gradSolVAR[unk][ivar2]=0; 
// 	  }
// 	  unsigned SolIndex = mlSol->GetIndex       (Solname[unk].c_str());
// 	  unsigned SolType  = mlSol->GetSolutionType(Solname[unk].c_str());
// 	  
// 	  for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
// 	    double soli   = (*sol->_Sol[SolIndex])(node_pos_sol[ SolFEType[unk] ][i]);
// // 	    double soli = (*sol->_Sol[SolIndex])(msh->GetMetisDof(node_pos_sol[ SolFEType[press_type_pos] ][i],SolType));
// 	    SolVAR[unk] += phi_gss_fe[ SolFEType[unk] ][i]*soli;
// 	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
// 	      gradSolVAR[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim+ivar2]*soli; 
// 	    }
// 	  }
// 	  
// 	}  
//  //end unknowns eval at gauss points ********************************
// 	
// 	
// 	
// 	
// 	
//   //begin NS block row *********************************
//      for(unsigned ivar_block=0; ivar_block<dim; ivar_block++) {  //1st row blocks A B' 
// 	// *** phi_i loop ***
// 	for(unsigned i_u=0; i_u < Sol_n_el_dofs[vel_type_pos]; i_u++) { //1st row
// 	
// 	  //*************************************************
// // 	  double Lap_rhs_i=0;
// // 	    for(unsigned ivar2=0; ivar2<dim; ivar2++) { //RHS column Velocity values
// // 	      Lap_rhs_i += phi_x_gss_fe[SolFEType[vel_type_pos]][i_u*dim+ivar2]*gradSolVAR[ivar_block][ivar2];
// // 	    }
// // 	    
// // 	    Res[SolPdeIndex[ivar_block]][i_u] += ( -IRe*Lap_rhs_i + /*Picard iteration*/SolVAR[dim]*phi_x_gss_fe[SolFEType[vel_type_pos]][i_u*dim+ivar_block] + force[ivar_block] * phi_gss_fe[SolFEType[vel_type_pos]][i_u])*weight;
// 	  //***************************************************
// 	  Res[SolPdeIndex[ivar_block]][i_u] = 5.; // testing with identity
// 	   
// 	   
// 	   
// 	   // *** phi_j loop *** 
// 	   for(unsigned j_u=0; j_u < Sol_n_el_dofs[vel_type_pos]; j_u++) { // Matrix 4x4 block 1st row vel values of 3x3 block, especially A
// 
// 	     //**************************************************
// // 	    double Lap_ij=0;
// // 	      for(unsigned ivar_lap=0; ivar_lap<dim; ivar_lap++) {
// // 		Lap_ij  += phi_x_gss_fe[SolFEType[vel_type_pos]][i_u*dim+ivar_lap]*phi_x_gss_fe[SolFEType[vel_type_pos]][j_u*dim+ivar_lap];
// // 	      }
// // 
// // 		Jac[ SolPdeIndex[ivar_block] ][ SolPdeIndex[ivar_block] ][ i_u*Sol_n_el_dofs[vel_type_pos]+j_u ] += ( IRe*Lap_ij)*weight;
// 	      //**************************************************
// 	     Jac[ SolPdeIndex[ivar_block] ][ SolPdeIndex[ivar_block] ][ i_u*Sol_n_el_dofs[vel_type_pos]+j_u ] = 1.; // reserving Identity values
// 	      
// 	    }//end phij loop
// 	      
// 	      
// 	      
// 	      //************************************************************
// // 	    // *** phiP_j loop ***
// // 	      for(unsigned j_p = 0; j_p < Sol_n_el_dofs[press_type_pos]; j_p++){ // Matrix block 1st row's last col values, especially B' 
// // 		Jac[ SolPdeIndex[ivar_block] ][ SolPdeIndex[press_type_pos] ][ i_u*Sol_n_el_dofs[press_type_pos]+j_p ]  -=  phi_x_gss_fe[SolFEType[vel_type_pos]][i_u*dim+ivar_block]*phi_gss_fe[SolFEType[press_type_pos]][j_p]*weight;
// // 	      }//end phiP_j loop
// 	      //************************************************************
// 	      
// 	      
// 	    }  //end phii loop
// 	    
//         }  //end ivar_block
//    //end NS block row *********************************
// 
//    
//    
//    
//    
//    
//    //begin div u block row *********************************
//     for(unsigned ivar_block=0; ivar_block<1; ivar_block++) { // Matrix block 2nd row values, B and null
//       
// 	  //*******************************************************************
// // 	  double div = 0;
// // 	  for(unsigned ivar=0; ivar<dim; ivar++) {
// // 	    div += gradSolVAR[ivar][ivar];
// // 	  }
// 	  //********************************************************************
//       
// 	for(unsigned i_p=0; i_p < Sol_n_el_dofs[press_type_pos]; i_p++) { //RHS column Pressure values
// 
// 	    //************************************************************************
// // 	  //RESIDUALS B block ===========================
// // 	  Res[SolPdeIndex[press_type_pos]][i_p] += phi_gss_fe[SolFEType[press_type_pos]][i_p]*div*weight;
// // // 	  Res[SolPdeIndex[press_type_pos]][i_p] += (phiP_gss[i]*div + /*penalty*ILambda*phiP_gss[i]*SolVAR[dim]*/ 
// // // 	                             + 0.*((hk*hk)/(4.*IRe))*alpha*(GradSolP[0]*phiV_x_gss[i*dim + 0] + GradSolP[1]*phiV_x_gss[i*dim + 1]) )*weight; //REMOVED !!
// 	  //*********************************************************************
// 	  
// 	  Res[SolPdeIndex[press_type_pos]][i_p]=2.;
// 	  
// 	}  //end phiP_i loop
// 
// 	    // *** phi_j loop ***
//     for(unsigned jvar_block=0; jvar_block<dim; jvar_block++) {
//      for(unsigned i_p=0; i_p<Sol_n_el_dofs[press_type_pos]; i_p++) {
// 	 for(unsigned j_u = 0; j_u < Sol_n_el_dofs[press_type_pos/*vel_type_pos*/]; j_u++) { // Matrix block 2nd row values, especially B
// 		Jac[ SolPdeIndex[press_type_pos] ][ SolPdeIndex[press_type_pos/*jvar_block*/] ][ i_p*Sol_n_el_dofs[press_type_pos/*vel_type_pos*/]+j_u ] =1.;
// 	
// 		//***************************************************
// // 	   Jac[ SolPdeIndex[press_type_pos] ][ SolPdeIndex[jvar_block] ][ i_p*Sol_n_el_dofs[vel_type_pos]+j_u ] -= phi_gss_fe[SolFEType[press_type_pos]][i_p]*phi_x_gss_fe[SolFEType[vel_type_pos]][j_u*dim+jvar_block]*weight;
// 	        //********************************************************
// 	   
// 		}  //end phij loop
// 	     }//end phiP_i loop
// 	     
// // 	if(assembleMatrix * penalty){  //block nDofsP
// // 	  // *** phi_i loop ***
// // 	  for(unsigned i=0; i<nDofsP; i++){
// // 	    // *** phi_j loop ***
// // 	    for(unsigned j=0; j<nDofsP; j++){
// // 	      //Jac[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nDofsP+j]-= ILambda*phiP_gss[i]*phiP_gss[j]*weight;
// // 	      for(unsigned ivar=0; ivar<dim; ivar++) {
// // // 	        Jac[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nDofsP+j] -= ((hk*hk)/(4.*IRe))*alpha*(phiV_x_gss[i*dim + ivar]*phiV_x_gss[j*dim + ivar])*weight; //REMOVED !!
// // 		Jac[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nDofsP+j] -= 0.;
// // 	      }
// // 	    }
// // 	  }
// // 	}   //end if penalty	     
// 	     
//          } //end column u jvar 
//  
//        }  
//    //end div u block row *********************************
//    
//    double vel_desired[2] = {10.,0.};
// 
//    
//    
// // /*   
// // //     //begin NSADJ block row *********************************
// //      for(unsigned ivar_block=0; ivar_block<dim; ivar_block++) {  //3rd row blocks I 
// //        unsigned ivar_block_adj=ivar_block+adj_pos_begin;
// // 	// *** phi_i loop ***
// // 	for(unsigned i_u=0; i_u < Sol_n_el_dofs[adj_vel_type_pos]; i_u++) { //RHS 2nd group vel
// // 	    
// // 	    //*********************************************************************
// // // 	    double Lap_rhs_i=0;
// // // 	    for(unsigned ivar2=0; ivar2<dim; ivar2++) { //RHS column Velocity values
// // // 	      Lap_rhs_i += phi_x_gss_fe[SolFEType[adj_vel_type_pos]][i_u*dim+ivar2]*gradSolVAR[ivar_block_adj][ivar2];
// // // 	    }
// // 	    
// // // 	    Res[SolPdeIndex[ivar_block_adj]][i_u] += 0.*( -IRe*Lap_rhs_i + /*Picard iteration*/SolVAR[dim]*phi_x_gss_fe[SolFEType[adj_vel_type_pos]][i_u*dim+ivar_block] +(/*SolVAR[ivar_block]*/- vel_desired[ivar_block])* phi_gss_fe[SolFEType[adj_vel_type_pos]][i_u*dim+ivar_block] )*weight;
// // 
// // 
// // 	    
// // // // 	     Res[SolPdeIndex[ivar_block]][i_u] -=/* fRHS[ivar_block-adj_pos_begin] */0. *phi_gss_fe[SolFEType[ivar_block]][i_u] *weight;
// // 
// // 	  //*************************************************************************************
// // 	  Res[SolPdeIndex[ivar_block_adj]][i_u]=7.;
// // 	  
// // 	  
// // 	    // *** phi_j loop ***
// // 	    for(unsigned j_u=0; j_u < Sol_n_el_dofs[adj_vel_type_pos]; j_u++) { // Matrix 3x3 block 3rd row vel values, I
// // 
// // 		//***************************************************************************
// // // 	    double Lap_ij=0;
// // // 	      for(unsigned ivar_lap=0; ivar_lap<dim; ivar_lap++) {
// // // 		Lap_ij  += phi_gss_fe[SolFEType[adj_pos_begin]][i_u*dim+ivar_lap]*phi_gss_fe[SolFEType[adj_pos_begin]][j_u*dim+ivar_lap];
// // // 	      }
// // 
// // // 		Jac[ SolPdeIndex[ivar_block_adj] ][ SolPdeIndex[ivar_block_adj] ][ i_u*Sol_n_el_dofs[adj_vel_type_pos]+j_u ] += (IRe* Lap_ij)*weight;
// // // 		
// // // // 		Jac[ SolPdeIndex[ivar_block_adj] ][ SolPdeIndex[ivar_block] ][ i_u*Sol_n_el_dofs[adj_vel_type_pos]+j_u ] -=SolVAR[ivar_block]* phi_gss_fe[SolFEType[adj_vel_type_pos]][i_u*dim+ivar_block]*weight;
// // 
// // 		//***********************************************************************************
// // 		
// // 		Jac[ SolPdeIndex[ivar_block_adj] ][ SolPdeIndex[ivar_block_adj] ][ i_u*Sol_n_el_dofs[adj_vel_type_pos]+j_u ] = 1.;
// // 		
// // 	      }//end phij loop
// // 	      
// // 	      
// // 	      
// // 	     	    
// // 	  
// // // 	      	    // *** phiP_j loop ***
// // // 	      for(unsigned j_p = 0; j_p < Sol_n_el_dofs[press_type_pos]; j_p++){ // Matrix block 1st row's last col values, especially B' 
// // // // 		Jac[ SolPdeIndex[ivar_block_adj] ][ SolPdeIndex[adj_vel_type_pos] ][ i_u*Sol_n_el_dofs[adj_vel_type_pos]+j_p ]  -= SolVAR[ivar_block]* phi_gss_fe[SolFEType[adj_vel_type_pos]][i_u*dim+ivar_block]*weight;
// // // 
// // // 		Jac[ SolPdeIndex[ivar_block_adj] ][ SolPdeIndex[press_type_pos] ][ i_u*Sol_n_el_dofs[adj_vel_type_pos+press_type_pos]+j_p]  -=  phi_x_gss_fe[SolFEType[adj_vel_type_pos]][i_u*dim+ivar_block]*phi_gss_fe[SolFEType[press_type_pos]][j_p]*weight;
// // // 	      }//end phiP_j loop
// // 
// // 
// // 
// // 
// // 
// // 	      
// // 	    }  //end phii loop
// // 	    
// //         }  //end ivar_block
// //     //end NSADJ block row **********************************/
// 
//     //begin DIV LAMBDA block row *********************************
// //         for(unsigned ivar_block=0; ivar_block<1; ivar_block++) { // Matrix block 2nd row values, B and null
// //       
// // 	  double div = 0;
// // 	  for(unsigned ivar=0; ivar<dim; ivar++) {
// // 	    div += gradSolVAR[ivar][ivar];
// // 	  }
// // 	  
// // 	for(unsigned i_p=0; i_p < Sol_n_el_dofs[press_type_pos]; i_p++) { //RHS column Pressure values
// // 
// // 	  //RESIDUALS B block ===========================
// // 	  Res[SolPdeIndex[press_type_pos]][i_p] += phi_gss_fe[SolFEType[press_type_pos]][i_p]*div*weight;
// // // 	  Res[SolPdeIndex[press_type_pos]][i_p] += (phiP_gss[i]*div + /*penalty*ILambda*phiP_gss[i]*SolVAR[dim]*/ 
// // // 	                             + 0.*((hk*hk)/(4.*IRe))*alpha*(GradSolP[0]*phiV_x_gss[i*dim + 0] + GradSolP[1]*phiV_x_gss[i*dim + 1]) )*weight; //REMOVED !!
// // 	  	}  //end phiP_i loop
// // 
// // 	    // *** phi_j loop ***
// //     for(unsigned jvar_block=0; jvar_block<dim; jvar_block++) {
// //      for(unsigned i_p=0; i_p<Sol_n_el_dofs[press_type_pos]; i_p++) {
// // 	 for(unsigned j_u = 0; j_u < Sol_n_el_dofs[vel_type_pos]; j_u++) { // Matrix block 2nd row values, especially B
// // 		Jac[ SolPdeIndex[press_type_pos] ][ SolPdeIndex[jvar_block] ][ i_p*Sol_n_el_dofs[vel_type_pos]+j_u ] -= phi_gss_fe[SolFEType[press_type_pos]][i_p]*phi_x_gss_fe[SolFEType[vel_type_pos]][j_u*dim+jvar_block]*weight;
// // 	        }  //end phij loop
// // 	     }//end phiP_i loop
// // 	          
// //          } //end column u jvar 
// //  
// //        }  
// 
//    //end DIV LAMBDA block row *********************************
// 
//  
//  
//       }  // end gauss point loop
//       
//       
// //       //***************************boundary loop ************************************************************************
// //       unsigned nfaces = el->GetElementFaceNumber(kel);
// //       
// // 	std::vector< double > xx_bd(dim,0.);
// // 	vector < double > normal_bd(dim,0);   
// // 	// loop on faces
// // 
// // 	
// // 	for(unsigned jface=0; jface < nfaces; jface++) {
// // 	  // look for boundary faces
// // 	  if(el->GetFaceElementIndex(kel,jface)<0) {
// // 	    unsigned int face = -(msh->el->GetFaceElementIndex(kel,jface)+1);
// //              double tau=0.;
// // 	     
// // 	     bool flag_bd = mlSol->_SetBoundaryConditionFunction(xx_bd,"P",tau,face,0.);
// // 	    if( ! flag_bd /*&& tau!=0.*/){
// // 	      unsigned nDofsV_bd      = msh->el->GetElementFaceDofNumber(kel,jface,SolFEType[vel_type_pos]);
// // 	      unsigned nDofsX_bd = msh->el->GetElementFaceDofNumber(kel,jface,coordXType);
// //                 phiV_gss_bd.resize(nDofsV_bd);
// //                 phiV_x_gss_bd.resize(nDofsV_bd*dim);
// // 	      const unsigned felt = msh->el->GetElementFaceType(kel, jface);
// // 	      for(unsigned i=0; i < nDofsX_bd; i++) {
// // 		unsigned inode = msh->el->GetFaceVertexIndex(kel,jface,i)-1u;
// // 		unsigned inode_Metis = msh->GetMetisDof(inode,coordXType);
// // 		for(unsigned idim=0; idim<dim; idim++) {
// // 		  coordX_bd[idim][i] = (*msh->_coordinate->_Sol[idim])(inode_Metis);
// // 		}
// // 	      }
// // 	      for(unsigned igs=0; igs < msh->_finiteElement[felt][SolFEType[vel_type_pos]]->GetGaussPointNumber(); igs++) {
// // 		msh->_finiteElement[felt][SolFEType[vel_type_pos]]->JacobianSur(coordX_bd,igs,weight_bd,phiV_gss_bd,phiV_x_gss_bd,normal_bd);
// // 		// *** phi_i loop ***
// // 		for(unsigned i_bd=0; i_bd < nDofsV_bd; i_bd++) {
// // 		  unsigned int i_vol = msh->el->GetLocalFaceVertexIndex(kel, jface, i_bd);
// // 		for(unsigned idim=0; idim<dim; idim++) {
// // 		  Res[idim][i_vol]   += -1. * weight_bd *  tau *phiV_gss_bd[i_bd]* normal_bd[idim];
// //  		 }//velocity blocks
// // 		}
// // 	      }
// // 	    } //flag_bd
// // 	  }
// // 	}    
// 
//       //***************************************************************************************************************
//       
//       //--------------------------------------------------------------------------------------------------------
//     } // endif single element not refined or fine grid loop
//     //--------------------------------------------------------------------------------------------------------
// 
//     //Sum the local matrices/vectors into the Global Matrix/Vector
//     for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
//       RES->add_vector_blocked(Res[SolPdeIndex[ivar]],JACDof[ivar]);
//         for(unsigned jvar=0; jvar < n_unknowns; jvar++) {
// 	  if(assembleMatrix) JAC->add_matrix_blocked( Jac[ SolPdeIndex[ivar] ][ SolPdeIndex[jvar] ], JACDof[ivar], JACDof[jvar]);
//         }
//     }
//  
//    //--------------------------------------------------------------------------------------------------------  
//   } //end list of elements loop for each subdomain
//   
//   
//   JAC->close();
//   RES->close();
//   // ***************** END ASSEMBLY *******************
// }
