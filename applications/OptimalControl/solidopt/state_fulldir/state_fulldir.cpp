#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "LinearImplicitSystem.hpp"
#include "Parameter.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Files.hpp"
#include "Math.hpp"

#include "PetscMatrix.hpp"
#include <stdio.h>

#define NSUB_X  16
#define NSUB_Y  16

// #define MODEL "Linear_elastic"
// #define MODEL "Mooney-Rivlin" 
#define MODEL "Neo-Hookean"

using namespace femus;

  
double SetInitialCondition (const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {
         
           double value = 0.;

             if(!strcmp(name,"DX")) {
                 value = 0.;
             }
             else if(!strcmp(name,"DY")) {
                 value = 0.;
             }
             else if(!strcmp(name,"DZ")) {
                 value = 0.;
             }
             else if(!strcmp(name,"DP")) {
                 value = 0.;
             }
           
    
      return value;   
}  
  
  
  

bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  //1: bottom  //2: right  //3: top  //4: left
  
  bool dirichlet; 
  
      if (facename == 1) {
       if (!strcmp(SolName, "DX"))    { dirichlet = false/*true*/; value = 0.; }
  else if (!strcmp(SolName, "DY"))    { dirichlet = false/*true*/; value = 0.; } 
  	
      }

      if (facename == 2) {
       if (!strcmp(SolName, "DX"))    { dirichlet = true; value = 0.; }
  else if (!strcmp(SolName, "DY"))    { dirichlet = true; value = 0.; } 
  	
      }

      if (facename == 3) {
       if (!strcmp(SolName, "DX"))    { dirichlet = false/*true*/; value = 0.; }
  else if (!strcmp(SolName, "DY"))    { dirichlet = false/*true*/; value = 0.; } 
  	
      }

      if (facename == 4) {
       if (!strcmp(SolName, "DX"))    { dirichlet = true; value = 0.; }
  else if (!strcmp(SolName, "DY"))    { dirichlet = true; value = 0.; } 
  	
      }
      
  return dirichlet;
}




template <  class real_num, class real_num_mov >
void AssembleSolidMech_AD(MultiLevelProblem& ml_prob);





int main(int argc, char** args) {



  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  Files files;
        files.CheckIODirectories();
        files.RedirectCout();


  MultiLevelMesh ml_mesh;  // define multilevel mesh
  double scalingFactor = 1.;  // read coarse level mesh and generate finers level meshes
  
    //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;
 // *** apparently needed by non-AD assemble only **********************
  // add fluid material
  Parameter par(Lref,Uref);
  
 // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  double rhof = 1000;
  Fluid fluid(par,1,rhof,"Newtonian");
  std::cout << "Fluid properties: " << std::endl;
  std::cout << fluid << std::endl;

  
  // Generate Solid Object
  double E = 1500000;
  double ni = 0.5;
  double rhos = 1000;
  Solid solid;
  solid = Solid(par,E,ni,rhos,MODEL);

  std::cout << "Solid properties: " << std::endl;
  std::cout << solid << std::endl;

  
  ml_mesh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,"seventh");
//   /* "seventh" is the order of accuracy that is used in the gauss integration scheme
//      probably in the furure it is not going to be an argument of this function   */
  unsigned dimension = ml_mesh.GetDimension();

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);
  ml_mesh.PrintInfo();

  MultiLevelSolution ml_sol(&ml_mesh);

  // add variables to ml_sol
  ml_sol.AddSolution("DX",LAGRANGE,SECOND);
  ml_sol.AddSolution("DY",LAGRANGE,SECOND);
  if ( dimension == 3 ) ml_sol.AddSolution("DZ",LAGRANGE,SECOND);
  ml_sol.AddSolution("P",LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/,FIRST);

    // ======= Problem ========================
  MultiLevelProblem ml_prob(&ml_sol);  // define the multilevel problem attach the ml_sol object to it

  ml_prob.SetFilesHandler(&files);

  ml_sol.Initialize("All");
                        ml_sol.Initialize("DX", SetInitialCondition, &ml_prob);
                        ml_sol.Initialize("DY", SetInitialCondition, &ml_prob);
  if ( dimension == 3 ) ml_sol.Initialize("DZ", SetInitialCondition, &ml_prob);
                        ml_sol.Initialize("P",  SetInitialCondition, &ml_prob);

  // attach the boundary condition function and generate boundary data
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  ml_sol.GenerateBdc("All");


  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  ml_prob.parameters.set<Solid>("Solid") = solid;

  NonLinearImplicitSystem& system = ml_prob.add_system < NonLinearImplicitSystem > ("SolidMech");

  // add solution "u" to system
                      system.AddSolutionToSystemPDE("DX");
                      system.AddSolutionToSystemPDE("DY");
  if (dimension == 3) system.AddSolutionToSystemPDE("DZ");
                      system.AddSolutionToSystemPDE("P");
 

  // attach the assembling function to system
  system.SetAssembleFunction( AssembleSolidMech_AD< adept::adouble, adept::adouble >);

  // initilaize and solve the system
  system.init();
  
  // Solver and preconditioner
  system.SetOuterKSPSolver("preonly");
  system.SetSolverFineGrids(PREONLY);
  system.SetPreconditionerFineGrids(LU_PRECOND);

  //for Vanka   
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);
  system.SetDebugNonlinear(true);

//   system.SetMaxNumberOfNonLinearIterations(2);
//   system.SetNonLinearConvergenceTolerance(1.e-30);
//   system.SetDebugLinear(true);
//   system.SetMaxNumberOfLinearIterations(4);
//   system.SetAbsoluteLinearConvergenceTolerance(1.e-10);
 
  system.MGsolve();

  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  if ( dimension == 3 ) mov_vars.push_back("DZ");
  ml_sol.GetWriter()->SetMovingMesh(mov_vars);
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

 ml_sol.GetWriter()->Write(files.GetOutputPath(),"biquadratic", variablesToBePrinted);
 
    
  return 0;
}



template < class real_num, class real_num_mov >
void AssembleSolidMech_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  NonLinearImplicitSystem& 	mlPdeSys   	= ml_prob.get_system<NonLinearImplicitSystem> ("SolidMech");
  const unsigned 		level 		    = mlPdeSys.GetLevelToAssemble();
  bool 			assembleMatrix 		    = mlPdeSys.GetAssembleMatrix(); 
  const char* 			pdename         = mlPdeSys.name().c_str();

  Mesh*          		msh    		= ml_prob._ml_msh->GetLevel(level);
  elem*          		el     		= msh->el;

  MultiLevelSolution*  		ml_sol  = ml_prob._ml_sol;
  Solution*    			sol      	= ml_prob._ml_sol->GetSolutionLevel(level); 


  LinearEquationSolver* 	pdeSys  = mlPdeSys._LinSolver[level];
  SparseMatrix*    		JAC    		= pdeSys->_KK;
  NumericVector*   		RES         = pdeSys->_RES;


  const unsigned 	 dim = msh->GetDimension();
  const unsigned    dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned   nel   = msh->GetNumberOfElements();
  const unsigned igrid   = msh->GetLevel();
  const unsigned  iproc  = msh->processor_id();
  // reserve memory for the local standard vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

 
  // geometry (at dofs) *******************************************
  vector< vector < real_num_mov > >   coords(dim);
  vector  < vector  < double > > coords_hat(dim);
  unsigned coordsType = BIQUADR_FE; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  for(int i=0;i<dim;i++) {   
       coords[i].reserve(maxSize); 
       coords_hat[i].reserve(maxSize); 
 }
  // geometry *******************************************

  
 // solution variables *******************************************
   const unsigned int n_unknowns = mlPdeSys.GetSolPdeIndex().size();
//      enum Sol_pos {pos_dx = 0, pos_dy, pos_p};  //these are known at compile-time
//      enum Sol_pos {pos_dx = 0, pos_dy, pos_dz, pos_p};
//   constexpr unsigned int pos_dx = 0;  
//   constexpr unsigned int pos_dy = 1;  
//   constexpr unsigned int pos_dp = 2;  
//   if (dim == 3) {std::cout << "Uncomment the 3d enum and recompile"; abort(); }
   
  constexpr int disp_type_pos = 0;     //known at compile time
  const     int press_type_pos = dim;      //known at run time
  constexpr int state_pos_begin = 0;   //known at compile time


  vector < std::string > Solname(n_unknowns);
  Solname              [state_pos_begin+0] =                "DX";
  Solname              [state_pos_begin+1] =                "DY";
  if (dim == 3) Solname[state_pos_begin+2] =                "DZ";
  Solname              [state_pos_begin + press_type_pos] = "P";
  
  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys.GetSolPdeIndex(Solname[ivar].c_str());
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

   vector < unsigned int > Sol_n_el_dofs(n_unknowns);
  
  //==========================================================================================
  vector < vector < double > > phi_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_hat_gss_fe(NFE_FAMS);
  vector < vector < real_num_mov > > phi_x_gss_fe(NFE_FAMS);   //the derivatives depend on the Jacobian, which depends on the unknown, so we have to be adept here  //some of these should be real_num, some real_num_mov...
  vector < vector < real_num_mov > > phi_xx_gss_fe(NFE_FAMS);  //the derivatives depend on the Jacobian, which depends on the unknown, so we have to be adept here
  vector < vector < double > > phi_x_hat_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_xx_hat_gss_fe(NFE_FAMS);

  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(maxSize);
    phi_hat_gss_fe[fe].reserve(maxSize);
      phi_x_gss_fe[fe].reserve(maxSize*dim);
  phi_x_hat_gss_fe[fe].reserve(maxSize*dim);
     phi_xx_gss_fe[fe].reserve(maxSize*(3*(dim-1)));
 phi_xx_hat_gss_fe[fe].reserve(maxSize*(3*(dim-1)));
   }
   
  //=================================================================================================
  
  // quadratures ********************************
  real_num_mov weight = 0.;
    double weight_hat = 0.;
 
  //----------- at dofs ------------------------------
  // equation ***********************************
  vector < int > JACDof;   JACDof.reserve( n_unknowns *maxSize);
  vector < double > Res;   Res.reserve( n_unknowns *maxSize);
  vector < double > Jac;   Jac.reserve( n_unknowns *maxSize * n_unknowns *maxSize);

  vector < vector < real_num > > aResVAR(n_unknowns);   //real num are those that are adept despite of the moving domain
  vector < vector < real_num > > SolVAR_eldofs(n_unknowns);
  vector < vector < real_num > > gradSolVAR_eldofs(n_unknowns);
   

  for(int k=0; k<n_unknowns; k++) {
    SolVAR_eldofs[k].reserve(maxSize);
    gradSolVAR_eldofs[k].reserve(maxSize*dim);
    aResVAR[k].reserve(maxSize);
  }

  //------------ at quadrature points ---------------------
    vector < real_num > SolVAR_qp(n_unknowns);
    vector < vector < real_num > > gradSolVAR_qp(n_unknowns);
    
    vector < vector < real_num > > gradSolVAR_hat_qp(n_unknowns);
    for(int k=0; k<n_unknowns; k++) { 
	gradSolVAR_qp[k].resize(dim); 
 	gradSolVAR_hat_qp[k].resize(dim); 
    }
      

   // ------------------------------------------------------------------------
    // Physical parameters
    const double rhof	 	= ml_prob.parameters.get < Fluid>("Fluid").get_density();
    const double mu_lame 	= ml_prob.parameters.get < Solid>("Solid").get_lame_shear_modulus();
    const double lambda_lame 	= ml_prob.parameters.get < Solid>("Solid").get_lame_lambda();
    const double mus		= mu_lame / rhof;
    const double lambda	= lambda_lame / rhof;
    const int    solid_model	= ml_prob.parameters.get < Solid>("Solid").get_physical_model();

    const bool incompressible = (0.5  ==  ml_prob.parameters.get < Solid>("Solid").get_poisson_coeff()) ? 1 : 0;
    const bool penalty = ml_prob.parameters.get < Solid>("Solid").get_if_penalty();

    // gravity
    double _gravity[3] = {0., 1., 0.};
    // -----------------------------------------------------------------
 
    RES->zero();
  if (assembleMatrix) JAC->zero();
  
  adept::Stack & stack = FemusInit::_adeptStack;  // call the adept stack object for potential use of AD
  const assemble_jacobian< real_num, real_num_mov > assemble_jac;


  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

  // geometry *****************************
    short unsigned ielGeom = msh->GetElementType(iel);

      unsigned nDofsX = msh->GetElementDofNumber(iel, coordsType);    // number of coordinate element dofs
    
    for(int ivar=0; ivar<dim; ivar++) {
      coords[ivar].resize(nDofsX);
      coords_hat[ivar].resize(nDofsX);
    }
    
   for( unsigned i=0;i<nDofsX;i++) {
      unsigned coordsDof  = msh->GetSolutionDof(i, iel, coordsType);    // global to global mapping between coordinates node and coordinate dof // via local to global solution node
      for(unsigned ivar = 0; ivar < dim; ivar++) {
          //Fixed coordinates (Reference frame)
	coords_hat[ivar][i] = (*msh->_topology->_Sol[ivar])(coordsDof);
      }
    }

  //***************************************  
  
  // geometry end *****************************

  // equation *****************************
    unsigned nDofsD = msh->GetElementDofNumber(iel, SolFEType[disp_type_pos]);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType[state_pos_begin + press_type_pos]);    // number of solution element dofs
    
    unsigned nDofsDP = dim * nDofsD + nDofsP;
  // equation end *****************************      

  //resize ###################################################################
    JACDof.resize(nDofsDP);
       Res.resize(nDofsDP);            std::fill(Res.begin(), Res.end(), 0.);
       Jac.resize(nDofsDP * nDofsDP);  std::fill(Jac.begin(), Jac.end(), 0.);
  //resize ###################################################################

   //STATE###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
       Sol_n_el_dofs[k] = ndofs_unk;
       SolVAR_eldofs[k].resize(ndofs_unk);
             aResVAR[k].resize(ndofs_unk);
      std::fill(aResVAR[k].begin(), aResVAR[k].end(), 0.);
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);
       SolVAR_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);
       JACDof[i + k *nDofsD]/*[k][i]*/ = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);
      }
    }
  //STATE###################################################################
  

    
    //Moving coordinates (Moving frame)
      for (unsigned idim = 0; idim < dim; idim++) {
        for (int j = 0; j < nDofsD; j++) {
          coords[idim][j] = coords_hat[idim][j] + SolVAR_eldofs[SolIndex[idim]][j];
        }
      }

     // start a new recording of all the operations involving real_num_mov variables
   if( assembleMatrix ) stack.new_recording();
  
    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][SolFEType[disp_type_pos]/*solDType*/]->GetGaussPointNumber(); ig++) {

	// *** get Jacobian and test function and test function derivatives ***
      for(int fe=0; fe < NFE_FAMS; fe++) {
	ml_prob._ml_msh->_finiteElement[ielGeom][fe]->Jacobian(coords,    ig,weight,    phi_gss_fe[fe],    phi_x_gss_fe[fe],    phi_xx_gss_fe[fe]);
	ml_prob._ml_msh->_finiteElement[ielGeom][fe]->Jacobian(coords_hat,ig,weight_hat,phi_hat_gss_fe[fe],phi_x_hat_gss_fe[fe],phi_xx_hat_gss_fe[fe]);
      }
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
  	ml_prob._ml_msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coords,ig,weight,phi_gss_fe[BIQUADR_FE],phi_x_gss_fe[BIQUADR_FE],phi_xx_gss_fe[BIQUADR_FE]);
  	ml_prob._ml_msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coords_hat,ig,weight_hat,phi_hat_gss_fe[BIQUADR_FE],phi_x_hat_gss_fe[BIQUADR_FE],phi_xx_hat_gss_fe[BIQUADR_FE]);


  //begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk <  n_unknowns; unk++) {
	  SolVAR_qp[unk] = 0.;

	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_qp[unk][ivar2] = 0.; 
	    gradSolVAR_hat_qp[unk][ivar2] = 0.; 
	  }
	  
	  for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
	    SolVAR_qp[unk]    += phi_gss_fe[ SolFEType[unk] ][i] * SolVAR_eldofs[unk][i];
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim+ivar2] * SolVAR_eldofs[unk][i]; 
	      gradSolVAR_hat_qp[unk][ivar2] += phi_x_hat_gss_fe[ SolFEType[unk] ][i*dim+ivar2] * SolVAR_eldofs[unk][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************

// // //   // I x = 5 test ********************************
// // // 	for(unsigned i_unk=0; i_unk<n_unknowns; i_unk++) { 
// // // 	    for(unsigned i_dof=0; i_dof < Sol_n_el_dofs[i_unk]; i_dof++) {
// // // 		/*Res[ i_dof +  i_unk * nDofsD ]*/aResVAR[i_unk][i_dof] +=  (   5.* phi_gss_fe[SolFEType[i_unk]][i_dof] - SolVAR_qp[i_unk]*phi_gss_fe[SolFEType[i_unk]][i_dof] )*weight;
// // // // std::cout << aResVAR[i_unk][i_dof] << "----" << std::endl;
// // // 		// 		  for(unsigned j_unk=dim; j_unk<n_unknowns; j_unk++) {
// // // // 		  	for(unsigned j_dof=0; j_dof < Sol_n_el_dofs[j_unk]; j_dof++) {
// // // // 			  
// // // // 		              if (i_unk == j_unk )   {
// // // // 				Jac[i_dof*Sol_n_el_dofs[i_unk] + j_dof i +  k * nDofsD][ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ][ i_dof*Sol_n_el_dofs[i_unk] + j_dof ] += 
// // // // 				        ( phi_gss_fe[SolFEType[i_unk]][i_dof]*phi_gss_fe[SolFEType[j_unk]][j_dof] )*weight;
// // // // 			      }
// // // // 			  
// // // // 			} //j_dof
// // // // 		  }  //j_unk
// // // 	    }  //i_dof
// // // 	}  //i_unk
// // //  // I x = 5 test ********************************

 

 //*******************************************************************************************************
   vector < vector < real_num_mov > > Cauchy(3); for (int i = 0; i < Cauchy.size(); i++) Cauchy[i].resize(3);
   real_num_mov J_hat;
   real_num_mov trace_e;

    Cauchy = Solid::get_Cauchy_stress_tensor< real_num_mov >(solid_model, mus, lambda, dim, press_type_pos, gradSolVAR_qp, gradSolVAR_hat_qp, SolVAR_qp, SolPdeIndex, J_hat, trace_e);

    

              //BEGIN residual Solid Momentum in moving domain
          for (unsigned i = 0; i < nDofsD; i++) {

              real_num_mov Cauchy_direction[3] = {0., 0., 0.};

              for (int idim = 0.; idim < dim; idim++) {
                for (int jdim = 0.; jdim < dim; jdim++) {
                  Cauchy_direction[idim] += phi_x_gss_fe[SolFEType[idim]][i * dim + jdim] * Cauchy[idim][jdim];
                }
              }

              for (int idim = 0; idim < dim; idim++) {
                aResVAR[SolPdeIndex[idim]][i] += ( +  phi_gss_fe[SolFEType[idim]][i] * _gravity[idim] - Cauchy_direction[idim]) * weight;
              }

            }
              //END residual Solid Momentum in moving domain
              

              //BEGIN residual solid mass balance
            for (unsigned i = 0; i < nDofsP; i++) {
                
                real_num_mov div_displ = 0.;
                 for (int idim = 0; idim < dim; idim++) div_displ += gradSolVAR_qp[SolPdeIndex[idim]][idim];
                              
              aResVAR[SolPdeIndex[press_type_pos]][i] += phi_gss_fe[SolFEType[press_type_pos]][i] * Solid::get_mass_balance< real_num_mov >(solid_model, penalty, incompressible, lambda, weight, weight_hat, div_displ, J_hat, SolVAR_qp, SolPdeIndex, press_type_pos);
                
            }
              //END residual solid mass balance


    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------

    //copy the value of the adept::adoube aRes in double Res and store them in RES

     Res.resize(nDofsDP);    //resize
   for (unsigned  k = 0; k < n_unknowns; k++)  {
       for (int i = 0; i < Sol_n_el_dofs[k]; i++) Res[ i +  k * nDofsD ] = - aResVAR[k][i].value();
    }
    
    RES->add_vector_blocked(Res, JACDof);

 if (assembleMatrix) {   //Extract and store the Jacobian
     Jac.resize(nDofsDP * nDofsDP);
   // define the dependent, independent variables
    for (unsigned  k = 0; k < n_unknowns; k++) {   
	stack.dependent(&aResVAR[k][0], Sol_n_el_dofs[k]);
    }
    
    for (unsigned  k = 0; k < n_unknowns; k++) {   
	stack.independent(&SolVAR_eldofs[k][0], Sol_n_el_dofs[k]); 
    }
   
    // get the and store jacobian matrix (row-major)
    stack.jacobian(&Jac[0] , true);

    assemble_jacobian<real_num,real_num_mov>::print_element_residual(iel,Res,Sol_n_el_dofs,9,5);
    assemble_jacobian<real_num,real_num_mov>::print_element_jacobian(iel,Jac,Sol_n_el_dofs,9,5);


   JAC->add_matrix_blocked(Jac, JACDof, JACDof);

    stack.clear_independents();
    stack.clear_dependents();
 }  //end assemble matrix
    
    
  } //end element loop for each process


 if (assembleMatrix) JAC->close();
                     RES->close();
  
  
  
    //print JAC and RES to files
//     const unsigned nonlin_iter = mlPdeSys.GetNonlinearIt();
//     assemble_jacobian< real_num_mov,double >::print_global_jacobian(assembleMatrix, ml_prob, JAC, nonlin_iter);
//     assemble_jacobian< real_num_mov,double >::print_global_residual(ml_prob, RES, nonlin_iter);

  
}


