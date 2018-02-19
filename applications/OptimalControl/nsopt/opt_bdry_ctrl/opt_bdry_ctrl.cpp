// Solving Navier-Stokes problem using automatic differentiation and/or Picards method
// boundary conditions were set in 2D as, no slip in left,right of the box and top to bottom gravity is enforced
// therefore, U=V=0 on left and right, U=0 on top and bottom, V is free 



#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "LinearImplicitSystem.hpp"

#include "Fluid.hpp"
#include "Parameter.hpp"
#include "Files.hpp"
#include <stdio.h>


//*********************** Sets Number of subdivisions in X and Y direction *****************************************

#define NSUB_X  2
#define NSUB_Y  2

#define CTRL_FACE_IDX  3


using namespace femus;

 double force[3] = {0./*1.*/,0.,0.};
 double Vel_desired[3] = {0.125,0.,0.};
 double alpha_val = 1.;
 double beta_val = 1.;
 double gamma_val = 1.;
 double penalty_outside_control_boundary = 1.e50;       // penalty for zero control outside Gamma_c
 double penalty_ctrl = 1.e10;         //penalty for u=q
 
  
 int ElementTargetFlag(const std::vector<double> & elem_center) {

 //***** set target domain flag ********************************** 
  int target_flag = 0;
  
   if ( elem_center[0] > 0.   - 1.e-5  &&  elem_center[0] < 0.5   + 1.e-5  && 
        elem_center[1] > 0.25 - 1.e-5  &&  elem_center[1] < 0.75  + 1.e-5
  ) {
     
     target_flag = 1;
     
  }
  
     return target_flag;

}

 
bool SetBoundaryConditionOpt(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  //1: bottom  //2: right  //3: top  //4: left
  
  bool dirichlet = true;
   value = 0.;
  
                if (!strcmp(SolName, "GX"))       { dirichlet = false; }
                if (!strcmp(SolName, "GY"))       { dirichlet = false; }
                if (!strcmp(SolName, "GZ"))       { dirichlet = false; }
                if (!strcmp(SolName, "THETA"))    { dirichlet = false; }
                
// TOP ==========================  
      if (facename == CTRL_FACE_IDX) {
       if (!strcmp(SolName, "U"))    { value = 70.; } //lid - driven
  else if (!strcmp(SolName, "V"))    { value = 0.; } 
  	
      }
      
  return dirichlet;

}


//*********************** Find volume elements that contain a Control Face element *********************************

int ControlDomainFlag(const std::vector<double> & elem_center) {

 //***** set control domain flag ***** 
  double mesh_size = 1./NSUB_Y;
  int control_el_flag = 0;
   if ( elem_center[1] >  1. - mesh_size ) { control_el_flag = 1; }

     return control_el_flag;
}



void AssembleNavierStokesOpt(MultiLevelProblem& ml_prob);    

// double ComputeIntegral(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {
  
  FemusInit mpinit(argc, args, MPI_COMM_WORLD); 	// init Petsc-MPI communicator
  
    // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
	files.RedirectCout();
  
  MultiLevelMesh mlMsh;			// define multilevel mesh
  double scalingFactor = 1.;		// read coarse level mesh and generate finers level meshes

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
    mlMsh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,"seventh");
    
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
  mlSol.AddSolution("U", LAGRANGE, FIRST/*SECOND*/);
  mlSol.AddSolution("V", LAGRANGE, FIRST/*SECOND*/);
  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, FIRST/*SECOND*/);
  mlSol.AddSolution("P", LAGRANGE, FIRST);
  // adjoint =====================  
  mlSol.AddSolution("UADJ", LAGRANGE, FIRST/*SECOND*/);
  mlSol.AddSolution("VADJ", LAGRANGE, FIRST/*SECOND*/);
  if (dim == 3) mlSol.AddSolution("WADJ", LAGRANGE, FIRST/*SECOND*/);
  mlSol.AddSolution("PADJ", LAGRANGE, FIRST);
  // boundary condition =====================
  mlSol.AddSolution("GX", LAGRANGE, FIRST/*SECOND*/);
  mlSol.AddSolution("GY", LAGRANGE, FIRST/*SECOND*/);
  if (dim == 3) mlSol.AddSolution("GZ", LAGRANGE, FIRST/*SECOND*/);
  mlSol.AddSolution("THETA", DISCONTINOUS_POLYNOMIAL, ZERO/*SECOND*/);
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
//   LinearImplicitSystem& system_opt    = mlProb.add_system < LinearImplicitSystem > ("NSOpt");

  // ST ===================
  system_opt.AddSolutionToSystemPDE("U");
  system_opt.AddSolutionToSystemPDE("V");
  if (dim == 3) system_opt.AddSolutionToSystemPDE("W");
  system_opt.AddSolutionToSystemPDE("P");
//   ADJ ===================
  system_opt.AddSolutionToSystemPDE("UADJ");
  system_opt.AddSolutionToSystemPDE("VADJ");
  if (dim == 3) system_opt.AddSolutionToSystemPDE("WADJ");
  system_opt.AddSolutionToSystemPDE("PADJ");
  // BD ===================
  system_opt.AddSolutionToSystemPDE("GX");
  system_opt.AddSolutionToSystemPDE("GY");
  if (dim == 3)  system_opt.AddSolutionToSystemPDE("GZ");
  system_opt.AddSolutionToSystemPDE("THETA");
 
  
  // attach the assembling function to system
   system_opt.SetAssembleFunction(AssembleNavierStokesOpt);
    
  // initilaize and solve the system
  system_opt.init();
  system_opt.ClearVariablesToBeSolved();
  system_opt.AddVariableToBeSolved("All");
  
  system_opt.SetMaxNumberOfNonLinearIterations(1);
  system_opt.MLsolve();

//     ComputeIntegral(mlProb);
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/,"biquadratic", variablesToBePrinted);
 
  //Destroy all the new systems
  mlProb.clear();
 

  return 0;
}


// nonAD is in the old PETSc, edit this for the new PETSc
void AssembleNavierStokesOpt(MultiLevelProblem& ml_prob){
     
  //pointers
  NonLinearImplicitSystem& mlPdeSys  = ml_prob.get_system<NonLinearImplicitSystem>("NSOpt");
  const unsigned level = mlPdeSys.GetLevelToAssemble();

  bool assembleMatrix = mlPdeSys.GetAssembleMatrix(); 
   
  Solution*	 sol  	         = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  pdeSys	 = mlPdeSys._LinSolver[level];   
  const char* pdename            = mlPdeSys.name().c_str();
  
  MultiLevelSolution* mlSol = ml_prob._ml_sol;
  
  Mesh*		 msh    = ml_prob._ml_msh->GetLevel(level);
  elem*		 el	= msh->el;
  SparseMatrix*	 JAC	= pdeSys->_KK;
  NumericVector* RES 	= pdeSys->_RES;
    
  //data
  const unsigned dim 	= msh->GetDimension();
  unsigned nel		= msh->GetNumberOfElements();
  unsigned igrid	= msh->GetLevel();
  unsigned iproc 	= msh->processor_id();
 
  const unsigned maxSize = static_cast< unsigned > (ceil(pow(3,dim)));

  // geometry *******************************************
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  vector< vector < double> > coordX(dim);
  vector< vector < double> > coordX_bd(dim);
  for(int i=0;i<dim;i++) {   
       coordX[i].reserve(maxSize); 
       coordX_bd[i].reserve(maxSize); 
  }
  // geometry *******************************************

  
  // solution variables *******************************************
  
    std::vector<std::string> ctrl_name;
    ctrl_name.resize(dim);
    ctrl_name[0] = "GX";
    ctrl_name[1] = "GY";
     if (dim == 3)  ctrl_name[2] = "GZ";
 
  
  const int n_vars = dim+1;
  const int n_unknowns = 3*n_vars; //(2.*dim)+1; //state , adjoint of velocity terms and one pressure term
  const int vel_type_pos = 0;
  const int press_type_pos = dim;
  const int adj_vel_type_pos = vel_type_pos;
  const int state_pos_begin = 0;
  const int adj_pos_begin   = dim+1;
  const int ctrl_pos_begin   = 2*(dim+1);
  
  vector < std::string > Solname(n_unknowns);  // const char Solname[4][8] = {"U","V","W","P"};
  Solname              [state_pos_begin+0] =                "U";
  Solname              [state_pos_begin+1] =                "V";
  if (dim == 3) Solname[state_pos_begin+2] =                "W";
  Solname              [state_pos_begin + press_type_pos] = "P";

  Solname              [adj_pos_begin + 0] =              "UADJ";
  Solname              [adj_pos_begin + 1] =              "VADJ";
  if (dim == 3) Solname[adj_pos_begin + 2] =              "WADJ";
  Solname              [adj_pos_begin + press_type_pos] = "PADJ";

  Solname              [ctrl_pos_begin + 0] =                  ctrl_name[0];
  Solname              [ctrl_pos_begin + 1] =                  ctrl_name[1];
  if (dim == 3) Solname[ctrl_pos_begin + 2] =                  ctrl_name[2];
  Solname              [ctrl_pos_begin + press_type_pos] = "THETA";
  
  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys.GetSolPdeIndex(Solname[ivar].c_str());
    SolIndex[ivar]	= mlSol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= mlSol->GetSolutionType(SolIndex[ivar]);
  }

  vector < double > Sol_n_el_dofs(n_unknowns);
  
  //==========================================================================================
  // velocity ************************************
  //-----------state------------------------------
  vector < vector < double > > phi_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_xx_gss_fe(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(maxSize);
      phi_x_gss_fe[fe].reserve(maxSize*dim);
     phi_xx_gss_fe[fe].reserve(maxSize*(3*(dim-1)));
   }
   
   
  //boundary adjoint & ctrl shape functions  
  vector < vector < double > > phi_bd_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_bd_gss_fe(NFE_FAMS);
  
  //bdry vol adj  evaluated at bdry points
   vector < vector < double > > phi_vol_at_bdry_fe(NFE_FAMS);
    vector < vector < double > > phi_x_vol_at_bdry_fe(NFE_FAMS);

    for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_bd_gss_fe[fe].reserve(maxSize);
      phi_x_bd_gss_fe[fe].reserve(maxSize*dim);
  //bdry vol adj  evaluated at bdry points
         phi_vol_at_bdry_fe[fe].reserve(maxSize);
       phi_x_vol_at_bdry_fe[fe].reserve(maxSize*dim);

      
    }
   

  vector < vector < double > >  sol_adj_x_vol_at_bdry_gss(dim);
  for (int ldim =0; ldim < dim; ldim++) sol_adj_x_vol_at_bdry_gss[ldim].reserve(dim);
  vector < double > grad_dot_n_adj_res;
  for (int ldim =0; ldim < dim; ldim++) grad_dot_n_adj_res.reserve(dim);
  vector < double > grad_adj_dot_n;
  for (int ldim =0; ldim < dim; ldim++) grad_adj_dot_n.reserve(dim);
  
  //=================================================================================================
  
  // quadratures ********************************
  double weight = 0.;
  double weight_bd = 0.;
  
  
  // equation ***********************************
  vector < vector < int > > JACDof(n_unknowns); 
  vector < vector < double > > Res(n_unknowns); /*was F*/
  vector < vector < vector < double > > > Jac(n_unknowns); /*was B*/
 
  for(int i = 0; i < n_unknowns; i++) {     
    JACDof[i].reserve(maxSize);
      Res[i].reserve(maxSize);
  }
   
  if(assembleMatrix){
    for(int i = 0; i < n_unknowns; i++) {
      Jac[i].resize(n_unknowns);    
      for(int j = 0; j < n_unknowns; j++) {
	Jac[i][j].reserve(maxSize*maxSize);	
      }
    }
  }
  
  //----------- dofs ------------------------------
  vector < vector < double > > SolVAR_eldofs(n_unknowns); //sol_V,P_of_st,adj,ctrl
  vector < vector < double > > gradSolVAR_eldofs(n_unknowns);
  
  for(int k=0; k<n_unknowns; k++) {
    SolVAR_eldofs[k].reserve(maxSize);
    gradSolVAR_eldofs[k].reserve(maxSize*dim);    
  }

  //------------ at quadrature points ---------------------
  vector < double > SolVAR_qp(n_unknowns);   //sol_V,P_gss_of_st,adj,ctrl_ie@quadraturepoints
    vector < vector < double > > gradSolVAR_qp(n_unknowns);
    for(int k=0; k<n_unknowns; k++) {  gradSolVAR_qp[k].resize(dim);  }
      
    
  
//   double penalty_outside_control_boundary = 1.e50;       // penalty for zero control outside Gamma_c

  double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();

  // Set to zero all the global structures
   RES->zero();
    if(assembleMatrix) JAC->zero();
  
  // ****************** element loop *******************
 
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    
  // geometry *****************************
    short unsigned ielGeom = msh->GetElementType(iel);

   unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs
    
    for(int ivar=0; ivar<dim; ivar++) {
      coordX[ivar].resize(nDofsX);
    }
    
   for( unsigned i=0;i<nDofsX;i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof // via local to global solution node
      for(unsigned ivar = 0; ivar < dim; ivar++) {
	coordX[ivar][i] = (*msh->_topology->_Sol[ivar])(coordXDof);
      }
    }

     // elem average point 
    vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
    for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofsX; i++) {
         elem_center[j] += coordX[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofsX; }
  //***************************************  
  
  // geometry end *****************************
  
  // equation *****************************
    unsigned nDofsV = msh->GetElementDofNumber(iel, SolFEType[vel_type_pos]);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType[state_pos_begin + press_type_pos]);    // number of solution element dofs
    
    unsigned nDofsVadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin]);    // number of solution element dofs
    unsigned nDofsPadj = msh->GetElementDofNumber(iel,SolFEType[adj_pos_begin + press_type_pos]);    // number of solution element dofs

    unsigned nDofsGctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin]);    // number of solution element dofs
    unsigned nDofsThetactrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin + press_type_pos] );    // number of solution element dofs

    unsigned nDofsVP = dim * nDofsV + nDofsP;
    unsigned nDofsVP_tot = 3*nDofsVP;
  // equation end *****************************
  
  //***** set target domain flag ********************************** 
   int target_flag = 0;
   target_flag = ElementTargetFlag(elem_center);
   //***************************************       
   
 //************ set control flag *********************
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag(elem_center);
  std::vector< std::vector<int> > control_node_flag(dim);
	    for(unsigned idim=0; idim<dim; idim++) {
	      control_node_flag[idim].resize(nDofsGctrl);
    std::fill(control_node_flag[idim].begin(), control_node_flag[idim].end(), 0);
	    }
 //*************************************************** 
  
   //STATE###################################################################  
   unsigned int fake_iel_flag = 0;
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);	//nDofs_V,P_of_st,adj,ctrl
	Sol_n_el_dofs[k]=ndofs_unk;
       SolVAR_eldofs[k].resize(ndofs_unk);	//sol_V,P_of_st,adj,ctrl
       JACDof[k].resize(ndofs_unk); 
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
       SolVAR_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
       JACDof[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof

    if (k == SolPdeIndex[ctrl_pos_begin + press_type_pos] && JACDof[k][i] == 72) {
      fake_iel_flag = iel;
      std::cout << "ffffffffffffffffff " << i << " " << fake_iel_flag << std::endl;
    }
    }
    }
  //CTRL###################################################################
    
   //************ set fake theta flag *********************
   unsigned int bdry_int_constr_pos = 72; /*KKoffset[SolPdeIndex[PADJ]][iproc]*/
  std::vector<int> fake_theta_flag(nDofsThetactrl,0);
    for (unsigned i = 0; i < nDofsThetactrl; i++) {
      if ( JACDof[ SolPdeIndex[ctrl_pos_begin + press_type_pos] ] [i] == bdry_int_constr_pos) { 
	fake_theta_flag[i] = 1;
      }
       std::cout << "##################################################################################################### " << i << " " << JACDof[ SolPdeIndex[ctrl_pos_begin + press_type_pos] ] [i] << " " << fake_theta_flag[i] <<  std::endl; 
   }
    
   
 //************ end set fake theta flag *********************

// setting Jac and Res to zero ******************************* 
    for(int ivar=0; ivar<n_unknowns; ivar++) {
      
      Res[SolPdeIndex[ivar]].resize(Sol_n_el_dofs[ivar]);
      memset(&Res[SolPdeIndex[ivar]][0],0.,Sol_n_el_dofs[ivar]*sizeof(double));
    }
   
    for(int ivar=0; ivar<n_unknowns; ivar++) {
      for(int jvar=0; jvar<n_unknowns; jvar++) {
      if(assembleMatrix){  //MISMATCH
	Jac[ SolPdeIndex[ivar] ][ SolPdeIndex[jvar] ].resize(Sol_n_el_dofs[ivar]*Sol_n_el_dofs[jvar]);
	memset(&Jac[SolPdeIndex[ivar]][SolPdeIndex[jvar]][0],0.,Sol_n_el_dofs[ivar]*Sol_n_el_dofs[jvar]*sizeof(double));
      }
    }
  }
 // setting Jac and Res to zero ******************************* 
  
  
  
  //============ boundary integral constraint : delta_theta - control (here I first loop over boundary columns (with i_bdry/i_vol) and then over volume rows (with i) )============================================================================================
 if (iel == fake_iel_flag)   {

   //  for(unsigned i=0; i < nDofsThetactrl; i ++) {
// // if (fake_theta_flag[i] == 1) {
//     for (unsigned  kdim = 0; kdim < dim; kdim++) { 
// // 		    if(j_vol < nDofsGctrl) {
//     Jac[ctrl_pos_begin + press_type_pos][ctrl_pos_begin + kdim][i*nDofsGctrl + i_vol] +=  /*fake_theta_flag[i]  **/ weight_bd * ( phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry] * normal[kdim]);
//  std::cout <<"  jac for ctrl  " <<  Jac[ctrl_pos_begin + press_type_pos][ctrl_pos_begin + kdim][i*nDofsGctrl + i_vol]  <<  " phi " << phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry]<< " fake "<< fake_theta_flag[i] <<  " " << i_vol  << std::endl;
//             }
// //     } 
// // 		    }
// 	    }
 }
 
 //============ boundary integral constraint ============================================================================================

  
 
//========BoundaryLoop=====================================================================

  // Perform face loop over elements that contain some control face
  if (control_el_flag == 1) {
	  
    double tau=0.;
    vector<double> normal(dim,0);
	  
    for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
	std::vector < double > xyz_bdc(3,0.);  //not being used, because the boundaries are identified by the face numbers
	// look for boundary faces
	if(el->GetFaceElementIndex(iel,jface) < 0) {
	   unsigned int face = -( msh->el->GetFaceElementIndex(iel,jface)+1);

	   if(  face == CTRL_FACE_IDX) { //control face
//=================================================== 
		//we use the dirichlet flag to say: if dirichlet == true, we set 1 on the diagonal. if dirichlet == false, we put the boundary equation
	    std::vector<bool> dir_bool(dim);
	    for(unsigned idim=0; idim<dim; idim++) {
		dir_bool[idim] = false; //mlSol->GetBdcFunction()(xyz_bdc,ctrl_name[idim].c_str(),tau,face,0.);
	    }
	  
//=================================================== 
	    unsigned nve_bd = msh->GetElementFaceDofNumber(iel,jface, SolFEType[ctrl_pos_begin] ); //AAAAAAAAAAAAAAAAA
	    const unsigned felt_bd = msh->GetElementFaceType(iel, jface);    
	    for(unsigned i=0; i < nve_bd; i++) {
		unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
                unsigned iDof = msh->GetSolutionDof(i_vol, iel, coordXType);
		for(unsigned idim=0; idim<dim; idim++) {
		    coordX_bd[idim][i]=(*msh->_topology->_Sol[idim])(iDof);
		}
	    }
		
//========= initialize gauss quantities on the boundary ============================================
	    vector < double > SolVAR_bd_qp(n_unknowns);
	    vector < vector < double > > gradSolVAR_bd_qp(n_unknowns);
	    for(int k=0; k<n_unknowns; k++) {  gradSolVAR_bd_qp[k].resize(dim);  }

//========= gauss_loop boundary===============================================================
	    for(unsigned ig_bd=0; ig_bd < msh->_finiteElement[felt_bd][SolFEType[ctrl_pos_begin]]->GetGaussPointNumber(); ig_bd++) {
// 		for(int fe=0; fe < NFE_FAMS; fe++) {
// 		     ml_prob._ml_msh->_finiteElement[felt_bd][fe]->JacobianSur(coordX_bd,ig_bd,weight_bd,phi_bd_gss_fe[fe],phi_x_bd_gss_fe[fe],normal);
// 		}
		ml_prob._ml_msh->_finiteElement[felt_bd][SolFEType[ctrl_pos_begin]]->JacobianSur(coordX_bd,ig_bd,weight_bd,phi_bd_gss_fe[SolFEType[ctrl_pos_begin]],phi_x_bd_gss_fe[SolFEType[ctrl_pos_begin]],normal);
		ml_prob._ml_msh->_finiteElement[ielGeom][SolFEType[adj_pos_begin]]->ShapeAtBoundary(coordX,ig_bd,phi_vol_at_bdry_fe[SolFEType[adj_pos_begin]],phi_x_vol_at_bdry_fe[SolFEType[adj_pos_begin]]);
		  
//========== temporary soln for surface gradient on a face parallel to the X axis ===================
		double dx_dxi = 0.;
		const elem_type_1D* myeltype = static_cast<const elem_type_1D*>(msh->_finiteElement[felt_bd][SolFEType[ctrl_pos_begin]]);
		const double* myptr = myeltype->GetDPhiDXi(ig_bd);
		for (int inode = 0; inode < nve_bd; inode++) {
		      dx_dxi += myptr[inode] * coordX_bd[0][inode];
		}  
		for (int inode = 0; inode < nve_bd; inode++) {
                     for (int d = 0; d < dim; d++) {
                          if (d == 0 ) phi_x_bd_gss_fe[SolFEType[ctrl_pos_begin]][inode + d*nve_bd] = myptr[inode]* (1./ dx_dxi);
                          else  phi_x_bd_gss_fe[SolFEType[ctrl_pos_begin]][inode + d*nve_bd] = 0.;
                     }
                }
//========== temporary soln for surface gradient on a face parallel to the X axis ===================
		  
//========== compute gauss quantities on the boundary ===============================================
// 	for(unsigned unk = 0; unk < /*n_vars*/ n_unknowns; unk++) {
// 	    unsigned ndofsface_unk = msh->GetElementFaceDofNumber(iel,jface, SolFEType[unk] );	//nDofs_V,P_of_st,adj,ctrl
		/*unk*/ SolVAR_bd_qp[/*unk*/SolPdeIndex[ctrl_pos_begin]] = 0.;
		   for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
		      gradSolVAR_bd_qp[/*unk*/SolPdeIndex[ctrl_pos_begin]][ivar2] = 0.; 
		   }
	  
		   for(int i_bd = 0; i_bd < /*ndofsface_unk*/ nve_bd; i_bd++/*unsigned i = 0; i < Sol_n_el_dofs[unk]; i++*/) {
		       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
		        SolVAR_bd_qp[/*unk*/SolPdeIndex[ctrl_pos_begin]] += phi_bd_gss_fe[ SolFEType[/*unk*/ctrl_pos_begin] ][i_bd] * SolVAR_eldofs[/*unk*/SolPdeIndex[ctrl_pos_begin]][i_vol];
		       for(unsigned ivar2=0; ivar2<dim; ivar2++) {
			   gradSolVAR_bd_qp[/*unk*/SolPdeIndex[ctrl_pos_begin]][ivar2] += phi_x_bd_gss_fe[ SolFEType[ctrl_pos_begin] ][i_bd + ivar2 * nve_bd] * SolVAR_eldofs[/*unk*/SolPdeIndex[ctrl_pos_begin]][i_vol]; 
		       }
		   }
// 	}
 //end unknowns eval at gauss points ********************************
		      
// // // //=============== grad dot n for residual ========================================= 
// // // //     compute gauss quantities on the boundary through VOLUME interpolation
// // // 		for(unsigned ldim=0; ldim<dim; ldim++) {
// // // 		    std::fill(sol_adj_x_vol_at_bdry_gss[ldim].begin(), sol_adj_x_vol_at_bdry_gss[ldim].end(), 0.);
// // // 		  for (int iv = 0; iv < nDofsVadj; iv++)  {
// // //                      for (int d = 0; d < dim; d++) {
// // // 			   sol_adj_x_vol_at_bdry_gss[ldim][d] += SolVAR_eldofs[SolPdeIndex[adj_pos_begin]][iv] * phi_x_vol_at_bdry_fe[SolFEType[adj_pos_begin]][iv * dim + d];//notice that the convention of the orders x y z is different from vol to bdry
// // // 		     }
// // // 		  }  
// // // 		      
// // // 		  grad_dot_n_adj_res[ldim] = 0.;
// // // 		  for(unsigned d=0; d<dim; d++) {
// // // 		      grad_dot_n_adj_res[ldim] += sol_adj_x_vol_at_bdry_gss[ldim][d]*normal[d];  
// // // 		  }
// // // 		}
// // // //=============== grad dot n  for residual =========================================       
		  
//========== compute gauss quantities on the boundary ================================================

		
		
  // *** phi_i loop ***
		for(unsigned i_bdry=0; i_bdry < nve_bd; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
                 
//=============== construct control node flag field on the go  =========================================    
	      /* (control_node_flag[i])       picks nodes on \Gamma_c
	         (1 - control_node_flag[i])   picks nodes on \Omega \setminus \Gamma_c
	       */
		    for(unsigned idim=0; idim<dim; idim++) {
			if (dir_bool[idim] == false) { 
// 		std::cout << " found boundary control nodes ==== " << std::endl;
			    for(unsigned k=0; k < control_node_flag[idim].size(); k++) {
				  control_node_flag[idim][i_vol] = 1;
			    }
			}
		    }
//=============== construct control node flag field on the go  =========================================    
		
		  
//Boundary Residuals  and Jacobians ==================	

//============ ============================================================================================

		      for (unsigned  kdim = 0; kdim < dim; kdim++) {
			
			    double lap_res_dctrl_ctrl_bd = 0.;
			    double dctrl_dot_n_res = 0.;
			    for (unsigned jdim = 0; jdim < dim; jdim++) {
				  lap_res_dctrl_ctrl_bd += gradSolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim]*phi_x_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry + jdim*nve_bd];
				  dctrl_dot_n_res +=  phi_x_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry + jdim*nve_bd] * normal[jdim];
			    }
			
/*delta_state row */	    if(i_vol<nDofsV)     Res[kdim]                 [i_vol]  += - control_node_flag[kdim][i_vol] * penalty_ctrl * (SolVAR_eldofs[SolPdeIndex[kdim + state_pos_begin]][i_vol] - SolVAR_eldofs[SolPdeIndex[kdim + ctrl_pos_begin]][i_vol]);	    //u-g
/*delta_adjoint row */     if(i_vol<nDofsVadj)   Res[kdim + adj_pos_begin] [i_vol]  += 0.;	   
/*delta_control row */     if(i_vol<nDofsGctrl)  Res[kdim + ctrl_pos_begin][i_vol]  += - control_node_flag[kdim][i_vol] * weight_bd * (
                                                                        beta_val* SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_bd_gss_fe[SolFEType[kdim +  ctrl_pos_begin]][i_bdry]
								     + gamma_val* lap_res_dctrl_ctrl_bd
// 								     - grad_dot_n_adj_res[kdim] *  phi_bd_gss_fe[SolFEType[kdim +  ctrl_pos_begin]][i_bdry]
// 								     + SolVAR_bd_qp[SolPdeIndex[ctrl_pos_begin + press_type_pos]] /** phi_x_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry + kdim*nve_bd] */* normal[kdim]  /*dctrl_dot_n_res*/     
								  );	    
//  std::cout << " res of ctrl "  << Res[kdim + ctrl_pos_begin][i_vol] ;
//  std::cout<< std::endl;

// /*delta_theta row */ 	if( i_vol < nDofsThetactrl ) Res[/*72*/press_type_pos + ctrl_pos_begin][i_vol] += - control_node_flag[kdim][i_vol] * weight_bd * SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * normal[kdim] /*ctrl_dot_n_res*/;
			
		  }  

//============  ==================================================================================================
		  
	    for(unsigned j_bdry=0; j_bdry < nve_bd; j_bdry ++) {
 		      unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);
// std::cout << "%%%%%%  ibd  is " << i_bdry  <<  " ivol value " << i_vol << " jbd  is " << j_bdry  << " jvol value " << j_vol << " w8bd " <<  weight_bd<< std::endl;
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
		      if(i_vol < nDofsV && j_vol < nDofsV && i_vol == j_vol) {
			
			 for (unsigned  kdim = 0; kdim < dim; kdim++) { 
			      Jac[kdim][kdim][i_vol*nDofsV + j_vol]	+= penalty_ctrl * control_node_flag[kdim][i_vol];  //u
			 }
		      }
//BLOCK delta_state - control------------------------------------------------------------------------------------
		      if(i_vol < nDofsV && j_vol < nDofsGctrl && i_vol == j_vol) {
			  for (unsigned  kdim = 0; kdim < dim; kdim++) { 
				Jac[kdim][kdim + ctrl_pos_begin][i_vol*nDofsV + j_vol]	+= (-1.) * penalty_ctrl *control_node_flag[kdim][i_vol];  //-g
			  }
		      }		      


//DIAG BLOCK delta_control - control  --------------------------------------------------------------------------------------
		    if(i_vol < nDofsGctrl && j_vol < nDofsGctrl) {

// // 			  std::cout << " EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ************** " << iel << " " << jface << " " << j_bdry << phi_bd_gss_fe[0][j_bdry] << std::endl;
// // 
			  double lap_jac_dctrl_ctrl_bd = 0.;
			  for (unsigned ldim = 0; ldim < dim; ldim++) lap_jac_dctrl_ctrl_bd += phi_x_bd_gss_fe[SolFEType[ldim + ctrl_pos_begin]][i_bdry + ldim*nve_bd]*phi_x_bd_gss_fe[SolFEType[ldim + ctrl_pos_begin]][j_bdry + ldim*nve_bd];

			  for (unsigned kdim = 0; kdim < dim; kdim++) {
			      Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i_vol*nDofsGctrl + j_vol] += control_node_flag[kdim][i_vol] * weight_bd *(
				                                                                                beta_val * phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin] ][i_bdry] * phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin] ][j_bdry]
                                                                                                             + gamma_val * lap_jac_dctrl_ctrl_bd
													      );
//  std::cout <<"  jac for ctrl  " <<  Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i_vol*nDofsGctrl + j_vol] << std::endl;
			  }
                   }
                   
// std::cout << " &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& " << std::endl;
                   


                   
//BLOCK delta_control - theta -------------------------------------------------------------------------------------------------------
// // // 		  if( i_vol < nDofsGctrl && j_vol == 72glob ) { /*first column of delta_theta*/
// // // 			for (unsigned kdim = 0; kdim < dim; kdim++) {
// // // 			      Jac[kdim + ctrl_pos_begin][press_type_pos + ctrl_pos_begin][i_vol*nDofsPctrl + j_vol] +=  phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry] * phi_bd_gss_fe[SolFEType[press_type_pos + ctrl_pos_begin]][j_bdry] *normal[kdim]  * weight_bd;
// // // 			}
// // // 		      }
//ROW_BLOCK delta_theta - control ---------------------------------------------------------------------------------------------------------------
// // // 		  if(i_vol == 72glob && j_vol < nDofsGctrl){ /*first row of delta_theta*/
// // // 		      for (unsigned  kdim = 0; kdim < dim; kdim++) {  
// // // // 			    double ctrl_dot_n_res = 0.;
// // // // 			    for (unsigned jdim = 0; jdim < dim; jdim++) {
// // // // 				  ctrl_dot_n_res += SolVAR_bd_qp[SolPdeIndex[jdim + ctrl_pos_begin]] * normal[jdim];
// // // // 			    }
// // // 		      }
// // // 		      for(unsigned j_bdry=0; j_bdry < nve_bd; j_bdry ++) {
// // // 			  unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);
// // // 			      for (unsigned kdim = 0; kdim < dim; kdim++) {
// // // 				  Jac[press_type_pos + ctrl_pos_begin][kdim + ctrl_pos_begin][i_vol*nDofsGctrl + j_vol] +=  control_node_flag[kdim][i_vol] * weight_bd * phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]] [j_bdry]* normal[kdim];
// // // 			      }
// // // 		      }
// // // 		  }

                         }//end j_bdry loop
		    
// //BLOCK delta_control - adjoint------------------------------------------------------------------------------------------------
// //===================loop over j in the VOLUME (while i is in the boundary)	      
// 		for(unsigned j=0; j < nDofsVadj; j ++) {
// 			//=============== grad dot n  =========================================    
// 		    for(unsigned ldim=0; ldim<dim; ldim++) {
// 			    grad_adj_dot_n[ldim] = 0.;
// 			    for(unsigned d=0; d<dim; d++) {
// 				grad_adj_dot_n[ldim] += phi_x_vol_at_bdry_fe[SolFEType[adj_pos_begin]][j * dim + d]*normal[d];  //notice that the convention of the orders x y z is different from vol to bdry
// 			    }
// 		    }
// 			
// 			  //=============== grad dot n  =========================================    
// 
// // 			  std::cout << " gradadjdotn " << grad_adj_dot_n << std::endl;
//   
// 			  for (unsigned kdim = 0; kdim < dim; kdim++) {
// 				Jac[kdim + ctrl_pos_begin][kdim + adj_pos_begin][i_vol*nDofsGctrl + j] += control_node_flag[kdim][i_vol] * (-1) * (weight_bd  * phi_bd_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i_bdry]* grad_adj_dot_n[kdim]);    		      
// 			  }
// 		} // end j loop for volume 
		  
		    }//end i_bdry loop

                }  //end ig_bdry loop
	  
             }    //end if control face
	 }  //end if boundary faces
      }  // loop over element faces //jface   
  } //end if control element flag

//End Boundary Residuals  and Jacobians ==================	
    
    
//=======================VolumeLoop=====================================================    
// ********************** Gauss point loop *******************************
 for(unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[ielGeom][SolFEType[vel_type_pos]]->GetGaussPointNumber(); ig++) {
	
	// *** get Jacobian and test function and test function derivatives ***
      for(int fe=0; fe < NFE_FAMS; fe++) {
	ml_prob._ml_msh->_finiteElement[ielGeom][fe]->Jacobian(coordX,ig,weight,phi_gss_fe[fe],phi_x_gss_fe[fe],phi_xx_gss_fe[fe]);
      }
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
  	ml_prob._ml_msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coordX,ig,weight,phi_gss_fe[BIQUADR_FE],phi_x_gss_fe[BIQUADR_FE],phi_xx_gss_fe[BIQUADR_FE]);


 //begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk < /*n_vars*/ n_unknowns; unk++) {
	  SolVAR_qp[unk] = 0.;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_qp[unk][ivar2] = 0.; 
	  }
	  
	  for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
	    SolVAR_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i] * SolVAR_eldofs[unk][i];
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim+ivar2] * SolVAR_eldofs[unk][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************
	
	
//  // I x = 5 test ********************************
// 	for(unsigned i_unk=dim; i_unk<n_unknowns; i_unk++) { 
// 	    for(unsigned i_dof=0; i_dof < Sol_n_el_dofs[i_unk]; i_dof++) {
// 		/*if ( i_unk!=0 && i_unk!=1 && i_unk!=2 && i_unk!=3 && i_unk!=4 && i_unk!=6 && i_unk!=7 )*/  Res[SolPdeIndex[i_unk]][i_dof] +=  (               0.* phi_gss_fe[SolFEType[i_unk]][i_dof] 
// 		                                    - SolVAR_qp[i_unk]*phi_gss_fe[SolFEType[i_unk]][i_dof] )*weight;
// 		  for(unsigned j_unk=dim; j_unk<n_unknowns; j_unk++) {
// 		  	for(unsigned j_dof=0; j_dof < Sol_n_el_dofs[j_unk]; j_dof++) {
// 			  
// 		              if (i_unk==j_unk /*&& i_unk!=0 && i_unk!=1 && i_unk!=2 && i_unk!=3 && i_unk!=4 && i_unk!=6 && i_unk!=7*/)   {
// 				Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ][ i_dof*Sol_n_el_dofs[i_unk] + j_dof ] += 
// 				        ( phi_gss_fe[SolFEType[i_unk]][i_dof]*phi_gss_fe[SolFEType[j_unk]][j_dof] )*weight;
// 			      }
// 			  
// 			} //j_dof
// 		  }  //j_unk
// 	    }  //i_dof
// 	}  //i_unk
//  // I x = 5 test ********************************
 
 
 
//============ delta_state row ============================================================================================

  for (unsigned i = 0; i < nDofsV; i++) {
// FIRST ROW
	for (unsigned  kdim = 0; kdim < dim; kdim++) { // velocity block row 
	              double lap_res_du_u = 0.; 
// 		      double adv_res_uold_uold = 0.;
	      for (unsigned jdim = 0; jdim < dim; jdim++) {
		    lap_res_du_u += gradSolVAR_qp[SolPdeIndex[kdim]][jdim]*phi_x_gss_fe[SolFEType[kdim]][i * dim + jdim];
// 		   adv_res_uold_uold += SolVAR_qp[SolPdeIndex[jdim]] * gradSolVAR_qp[SolPdeIndex[kdim]][jdim];
	      }      
	      Res[kdim][i]   +=  (         + force[kdim] * phi_gss_fe[SolFEType[kdim]][i]
                                           - IRe*lap_res_du_u 
//                                            - adv_res_uold_uold * phi_gss_fe[ SolFEType[kdim] ][i]
					    + SolVAR_qp[SolPdeIndex[press_type_pos]] * phi_x_gss_fe[SolFEType[kdim]][i * dim + kdim]) * weight; 
	}	    
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
	for (unsigned j = 0; j < nDofsV; j++) {
		      double lap_jac_du_u = 0.;
// 		      double adv_unew_uold = 0.;
// 		      double adv_uold_unew = 0.;
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
		    lap_jac_du_u += phi_x_gss_fe[SolFEType[kdim]][i * dim + kdim]*phi_x_gss_fe[SolFEType[kdim]][j * dim + kdim];
// 		    adv_uold_unew += SolVAR_qp[SolPdeIndex[kdim]]*phi_x_gss_fe[ SolFEType[kdim] ][j * dim + kdim] * phi_gss_fe[ SolFEType[kdim] ][i];
// 		    adv_unew_uold += phi_gss_fe[ SolFEType[kdim] ][i] * gradSolVAR_qp[SolPdeIndex[kdim]][kdim] * phi_gss_fe[ SolFEType[kdim] ][j];
	      }
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
		Jac[kdim][kdim][i*nDofsV + j] += (   IRe*lap_jac_du_u 
						    /*+ adv_uold_unew 
						    + adv_unew_uold*/) * weight; 
	      }
	} //j_du_u loop
     
//BLOCK Pressure
      for (unsigned j = 0; j < nDofsP; j++) {
	    for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim][press_type_pos][i*nDofsP + j] += -( phi_gss_fe[SolFEType[press_type_pos]][j] * phi_x_gss_fe[SolFEType[kdim]][i * dim + kdim] ) * weight;
	    }
      }//j_press loop
   }//i_state loop

//DIV_state
  for (unsigned i = 0; i < nDofsP; i++) {
		    double div_u_du_qp =0.;
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      div_u_du_qp += gradSolVAR_qp[SolPdeIndex[kdim]][kdim] ;
      }
      Res[press_type_pos][i]  +=  ( (div_u_du_qp) * phi_gss_fe[SolFEType[press_type_pos]][i] ) * weight;
      for (unsigned j = 0; j < nDofsV; j++) {
	  for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      Jac[press_type_pos][kdim][i*nDofsV + j] += -( phi_gss_fe[SolFEType[press_type_pos]][i] * phi_x_gss_fe[SolFEType[kdim]][j * dim + kdim] ) * weight;
	  }
      } //j loop
   }//i_div_state
    //============ delta_state row ============================================================================================


    
//============ delta_adjoint row =============================================================================================
  
  for (unsigned i = 0; i < nDofsVadj; i++) {
// SECOND ROW
		    double res_dadj_u = 0.;
     for (unsigned kdim = 0; kdim < dim; kdim++) { 
		    double lap_res_dadj_adj = 0.;
// 		  double adv_res_dadj_uold = 0.;
// 		  double adv_res_uold_dadj = 0.;
		  res_dadj_u += SolVAR_qp[SolPdeIndex[kdim]]*phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i];
	   for (unsigned jdim = 0; jdim < dim; jdim++) {
		lap_res_dadj_adj += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]][jdim]*phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim + jdim];
// 		   adv_res_dadj_uold += phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i] * gradSolVAR_qp[SolPdeIndex[kdim]][jdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
// 		  adv_res_uold_dadj += SolVAR_qp[SolPdeIndex[kdim]] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim + jdim] * SolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]];
	   }
	  Res[kdim + adj_pos_begin][i] += (  - alpha_val*target_flag*Vel_desired[kdim] * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
					     + alpha_val*target_flag*(res_dadj_u) 					    
					    - IRe*lap_res_dadj_adj
					    + SolVAR_qp[SolPdeIndex[press_type_pos + adj_pos_begin]] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim + kdim]
					    ) * weight;
      }
      
//BLOCK delta_adjoint - state------------------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsV; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][kdim][i*nDofsVadj + j] += ( -alpha_val*target_flag*phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]*phi_gss_fe[SolFEType[kdim]][j] ) * weight;
	  }
     }//j_dadj_u loop


//DIAG BLOCK delta_adjoint - adjoint---------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsVadj; j++) {
		    double lap_jac_dadj_adj = 0.;
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
		  lap_jac_dadj_adj += phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim + kdim]*phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][j * dim + kdim];
	  }
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += ( IRe*lap_jac_dadj_adj ) * weight;
	  }
      } //j_dadj_adj loop
      
//BLOCK Pressure_adj
    for (unsigned j = 0; j < nDofsPadj; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][press_type_pos + adj_pos_begin][i*nDofsPadj + j] += -( phi_gss_fe[SolFEType[press_type_pos + adj_pos_begin]][j] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim + kdim] ) * weight;
	  }
    }//j_press_adj loop
  }//i_adj loop

//DIV_adj
  for (unsigned i = 0; i < nDofsPadj; i++) {
		double div_adj_dadj_qp = 0.;
      for (unsigned kdim = 0; kdim < dim; kdim++) {
	    div_adj_dadj_qp += gradSolVAR_qp[SolFEType[kdim + adj_pos_begin ]][kdim] ;
      }
      Res[press_type_pos + adj_pos_begin][i] += ( (div_adj_dadj_qp) * phi_gss_fe[SolFEType[press_type_pos + adj_pos_begin]][i] ) * weight;
      for (unsigned j = 0; j < nDofsVadj; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	    Jac[press_type_pos + adj_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += - ( phi_gss_fe[SolFEType[press_type_pos + adj_pos_begin]][i] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][j * dim + kdim] ) * weight;
	  }
      }//j loop
  }//i_div_adj

      //============ delta_adjoint row =============================================================================================

      //============ delta_control row ==================================================================================================
// THIRD ROW
 for (unsigned i = 0; i < nDofsGctrl; i++) {
      for (unsigned kdim = 0; kdim < dim; kdim++) { 
       Res[kdim + ctrl_pos_begin][i] += - penalty_outside_control_boundary * ( (1 - control_node_flag[kdim][i]) * (  SolVAR_eldofs[kdim + ctrl_pos_begin][i] - 0.)  );
      }

// // //BLOCK delta_control - adjoint--------------------------------------------------------------------------------------
//      for (unsigned j = 0; j < nDofsVadj; j++) {
// 	    if (i==j && i!=2 && i!=3) {
// 	  for (unsigned kdim = 0; kdim < dim; kdim++) {
// 		Jac[kdim + ctrl_pos_begin][kdim + adj_pos_begin][i*nDofsGctrl + j] = (1 - control_node_flag[kdim][i]) * 1.;  /*penalty_outside_control_boundary;*/
// 	      }
// 	    } //end i==j
//       }//j_dctrl_ctrl loop

// //DIAG BLOCK delta_control - control--------------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsGctrl; j++) {
	    if (i==j) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
		Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsGctrl + j] += penalty_outside_control_boundary *(1 - control_node_flag[kdim][i]);
//  std::cout <<"  jac for ctrl  in vol loop " <<  Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsGctrl + j] << " for " << i*nDofsGctrl + j << std::endl;
	      }
	    } //end i==j
      }//j_dctrl_ctrl loop
  }//i_ctrl loop

//  if (control_el_flag == 1) {
// // 		    for(unsigned idim=0; idim<dim; idim++) {
// //     std::fill(control_node_flag[idim].begin(), control_node_flag[idim].end(), 0);
// // 	    }
//   
//     for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
// 	if(el->GetFaceElementIndex(iel,jface) < 0) {
// 	   unsigned int face = -( msh->el->GetFaceElementIndex(iel,jface)+1);
// 	    if(face==CTRL_FACE_IDX){
// 	    unsigned nve_bd = msh->GetElementFaceDofNumber(iel,jface, SolFEType[ctrl_pos_begin] ); //AAAAAAAAAAAAAAAAA
// 		for(unsigned i_bdry=0; i_bdry < nve_bd; i_bdry++) {
// 		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
// // 		    for(unsigned idim=0; idim<dim; idim++) {
// // 				  control_node_flag[idim][i_vol] = 0;
// // 			    }
//  for (unsigned i = 0; i < nDofsGctrl; i++) {
//       for (unsigned kdim = 0; kdim < dim; kdim++) { 
//      /* if(i!=i_vol)*/ Res[kdim + ctrl_pos_begin][i] += - penalty_outside_control_boundary * ( (1 - control_node_flag[kdim][i]) * (  SolVAR_eldofs[kdim + ctrl_pos_begin][i] - 0.)  );
//       }
// 
// 
//       for (unsigned j = 0; j < nDofsGctrl; j++) {
// 	    if (i==j && i!=i_vol /*&& i!=2 && i!=3*/) {
// 	  for (unsigned kdim = 0; kdim < dim; kdim++) {
// 		/* if(i!=i_vol)*/ Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsGctrl + j] = (1 - control_node_flag[kdim][i]) * 1.;  /*penalty_outside_control_boundary;*/
// //  std::cout <<"  jac for ctrl  in vol loop " <<  Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsGctrl + j] << " for " << i*nDofsGctrl + j << std::endl;
// 	      }
// 	    } //end i==j
//       }//j_dctrl_ctrl loop
//   }//i_ctrl loop
//  			}
// 		    }
// 
// 	}
//       }
// }
 //============ delta_control row ==================================================================================================
 
 //============ delta_theta - theta row ==================================================================================================
  for (unsigned i = 0; i < nDofsThetactrl; i++) {
      for (unsigned j = 0; j < nDofsThetactrl; j++) {
			  if(i==j) {
				  Jac[press_type_pos + ctrl_pos_begin ][press_type_pos + ctrl_pos_begin ][i*nDofsThetactrl + j] = (1 - fake_theta_flag[i]) * 1.;
			  }
       }//j_theta loop
    }//i_theta loop

 //============ delta_theta row ==================================================================================================
 
      }  // end gauss point loop
      
    
    
    
 //--------------------PRINTING------------------------------------------------------------------------------------
 // Add the local Matrix/Vector into the global Matrix/Vector
  std::cout << " -------------------------Element = " << iel << " ----------------------Jac -------------------------- " << std::endl;      
  for(unsigned i_unk=0; i_unk < n_unknowns; i_unk++) {
    unsigned ndofs_unk_i = msh->GetElementDofNumber(iel, SolFEType[i_unk]);	//nDofs_V,P_of_st,adj,ctrl//Sol_n_el_dofs[k]
    std::cout << " ======== Row ==== " << i_unk << " Unk_i ===================================== " << std::endl;      
    for(unsigned j_unk=0; j_unk < n_unknowns; j_unk++) {
	unsigned ndofs_unk_j = msh->GetElementDofNumber(iel, SolFEType[j_unk]);	//nDofs_V,P_of_st,adj,ctrl
	std::cout << " ======= Column ==== " << j_unk << " Unk_j ================== " << std::endl;      
	for (unsigned i = 0; i < ndofs_unk_i; i++) {
// // // 	      std::cout << " " << std::setfill(' ') << std::setw(10) << Res[SolPdeIndex[i_unk]][ i ];
// // // 	      std::cout << std::endl;
            for (unsigned j = 0; j < ndofs_unk_j; j++) {
	    std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ][ i*ndofs_unk_i + j ];
	    }  //j end
	    std::cout << std::endl;
	 } //i end
     } //j_unk end
   } //i_unk end
 //--------------------PRINTING------------------------------------------------------------------------------------

      //***************************************************************************************************************

    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned i_unk=0; i_unk < n_unknowns; i_unk++) {
      RES->add_vector_blocked(Res[SolPdeIndex[i_unk]],JACDof[i_unk]);
        for(unsigned j_unk=0; j_unk < n_unknowns; j_unk++) {
	  if(assembleMatrix) JAC->add_matrix_blocked( Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ], JACDof[i_unk], JACDof[j_unk]);
        }
    }
    
//      if (control_el_flag == 1) {
//        std::vector<int>  bdry_int_constr_pos_vec(1,bdry_int_constr_pos);
// 	  for (unsigned kdim = 0; kdim < dim; kdim++) { RES->add_vector_blocked(Res[SolPdeIndex[n_unknowns-2-kdim]],JACDof[n_unknowns-2-kdim]); }
//                                                         RES->add_vector_blocked(Res[SolPdeIndex[n_unknowns-1]],JACDof[n_unknowns-1]);
// 	  if(assembleMatrix) {
// 	    for (unsigned kdim = 0; kdim < dim; kdim++) {
// 	    JAC->add_matrix_blocked( Jac[ SolPdeIndex[n_unknowns-1] ][ SolPdeIndex[n_unknowns-2-kdim] ], bdry_int_constr_pos_vec, JACDof[n_unknowns-2-kdim]);
// 	    JAC->add_matrix_blocked( Jac[ SolPdeIndex[n_unknowns-2-kdim] ][ SolPdeIndex[n_unknowns-1] ], JACDof[n_unknowns-2-kdim], bdry_int_constr_pos_vec);
// 	    }
// 	 }
//      }  //add control boundary element contributions
    
 
   //--------------------------------------------------------------------------------------------------------  
  } //end list of elements loop for each subdomain
  
  
  JAC->close();
  RES->close();
  
  JAC->print();
  // ***************** END ASSEMBLY *******************
}



double ComputeIntegral(MultiLevelProblem& ml_prob) {


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
  vector< vector < double> > coordX_bd(dim);

  for (unsigned  k = 0; k < dim; k++) { 
    coordX[k].reserve(maxSize);
       coordX_bd[k].reserve(maxSize); 
  }
  
  double weight = 0.;
  double weight_bd = 0.;
  
  
  //geometry *******************************

//STATE######################################################################
  //velocity *******************************
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if (dim == 3) solVIndex[2] = mlSol->GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"
  
  vector < vector < double > >  solV(dim);    // local solution
  vector <double >  V_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
  }

  
  vector <double> phiV_gss;  // local test function
  vector <double> phiV_x_gss; // local test function first order partial derivatives
  vector <double> phiV_xx_gss; // local test function second order partial derivatives

  phiV_gss.reserve(maxSize);
  phiV_x_gss.reserve(maxSize * dim);
  phiV_xx_gss.reserve(maxSize * dim2);
  
  
  //velocity *******************************
   

//STATE######################################################################
  

//CONTROL_@bdry######################################################################
  vector < unsigned > solVctrlIndex(dim);
  solVctrlIndex[0] = mlSol->GetIndex("GX");    // get the position of "U" in the ml_sol object
  solVctrlIndex[1] = mlSol->GetIndex("GY");    // get the position of "V" in the ml_sol object
  if (dim == 3) solVctrlIndex[2] = mlSol->GetIndex("GZ");      // get the position of "V" in the ml_sol object

  unsigned solVctrlType = mlSol->GetSolutionType(solVctrlIndex[0]);    // get the finite element type for "u"
  
  vector < vector < double > >  solVctrl(dim);    // local solution
  vector < double >   Vctrl_gss(dim, 0.);    //  solution
   
 for (unsigned  k = 0; k < dim; k++) {
    solVctrl[k].reserve(maxSize);
  }

  
  vector <double> phiVctrl_gss_bd;  // local test function
  vector <double> phiVctrl_x_gss_bd; // local test function first order partial derivatives
  vector <double> phiVctrl_xx_gss_bd; // local test function second order partial derivatives

  phiVctrl_gss_bd.reserve(maxSize);
  phiVctrl_x_gss_bd.reserve(maxSize * dim);
  phiVctrl_xx_gss_bd.reserve(maxSize * dim2);
  
  
   

//CONTROL_@bdry######################################################################

// Vel_desired##################################################################
  vector <double> phiVdes_gss;  // local test function
  vector <double> phiVdes_x_gss; // local test function first order partial derivatives
  vector <double> phiVdes_xx_gss; // local test function second order partial derivatives

  phiVdes_gss.reserve(maxSize);
  phiVdes_x_gss.reserve(maxSize * dim);
  phiVdes_xx_gss.reserve(maxSize * dim2);

//   vector< vector < double > >  solVdes(dim);    // local solution
  vector <double>  solVdes(dim,0.);
  vector<double> Vdes_gss(dim, 0.);  
  
//  for (unsigned  k = 0; k < dim; k++) {
//     solVdes[k].reserve(maxSize);
//   }
//   
//   double* Vdes_gss [3] = Vel_desired/*= 0.*/;


// Vel_desired##################################################################

  
vector<double> integral(dim);

double  integral_target_alpha = 0.;

double	integral_beta   = 0.;
double	integral_gamma  = 0.;
  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

// geometry
    short unsigned ielGeom = msh->GetElementType(iel);    // element geometry type

    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs
    
// equation
    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
//     unsigned nDofsVdes = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsVctrl = msh->GetElementDofNumber(iel, solVctrlType);    // number of solution element dofs
    
    for (unsigned  k = 0; k < dim; k++) {       coordX[k].resize(nDofsX);    }
     
    for (unsigned  k = 0; k < dim; k++)  {
      solV[k].resize(nDofsV);
      solVctrl[k].resize(nDofsVctrl);
//       solVdes[k].resize(nDofsVdes);
    }


    // geometry ************
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }
    
     // elem average point 
    vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
  for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofsX; i++) {
         elem_center[j] += coordX[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofsX; }
  //*************************************** 
  
  //***** set target domain flag ********************************** 
   int target_flag = 0;
   target_flag = ElementTargetFlag(elem_center);
//***************************************       
    
    
 //STATE###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
      }
    }
//STATE###################################################################

//CONTROL###################################################################  
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVctrlDof = msh->GetSolutionDof(i, iel, solVctrlType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solVctrl[k][i] = (*sol->_Sol[solVctrlIndex[k]])(solVctrlDof);      // global extraction and local storage for the solution
      }
    }
//CONTROL###################################################################




  //DESIRED VEL###################################################################  
    // velocity ************
//     for (unsigned i = 0; i < nDofsV; i++) {
//       unsigned solVdesDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < solVdes.size() /*dim*/; k++) {
        solVdes[k]/*[i]*/ = Vel_desired[k] /*(*sol->_Sol[solVIndex[k]])(solVdesDof)*/;      // global extraction and local storage for the solution
      }
//     }
 //DESIRED VEL###################################################################

 
 //************ set control flag *********************
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag(elem_center);
  std::vector< std::vector<int> > control_node_flag(dim);
	    for(unsigned idim=0; idim<dim; idim++) {
	      control_node_flag[idim].resize(nDofsVctrl);
   /*if (control_el_flag == 0)*/ std::fill(control_node_flag[idim].begin(), control_node_flag[idim].end(), 0);
	    }
 //*************************************************** 

//========BoundaryLoop=====================================================================

  // Perform face loop over elements that contain some control face
  if (control_el_flag == 1) {
	  
    double tau=0.;
    vector<double> normal(dim,0);
	  
    for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
	std::vector < double > xyz_bdc(3,0.);  //not being used, because the boundaries are identified by the face numbers
	// look for boundary faces
	if(el->GetFaceElementIndex(iel,jface) < 0) {
	   unsigned int face = -( msh->el->GetFaceElementIndex(iel,jface)+1);

	   if(  face == CTRL_FACE_IDX) { //control face
// //=================================================== 
// 		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
// 	    std::vector<bool> dir_bool; dir_bool.resize(dim);
// 	    for(unsigned idim=0; idim<dim; idim++) {
// 		dir_bool[idim] = mlSol->GetBdcFunction()(xyz_bdc,ctrl_name[idim].c_str(),tau,face,0.);
// 	    }
	  
//=================================================== 
	    unsigned nve_bd = msh->GetElementFaceDofNumber(iel,jface, solVctrlType ); //AAAAAAAAAAAAAAAAA
	    const unsigned felt_bd = msh->GetElementFaceType(iel, jface);    
	    for(unsigned i=0; i < nve_bd; i++) {
		unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
                unsigned iDof = msh->GetSolutionDof(i_vol, iel, coordXType);
		for(unsigned idim=0; idim<dim; idim++) {
		    coordX_bd[idim][i]=(*msh->_topology->_Sol[idim])(iDof);
		}
	    }
		
//========= initialize gauss quantities on the boundary ============================================
    vector < double >   Vctrl_bd_qp(dim, 0.);    //  solution@bdry
    vector < vector < double > > gradVctrl_bd_qp(dim);
      for (unsigned  k = 0; k < dim; k++) {
          gradVctrl_bd_qp[k].resize(dim);
          std::fill(gradVctrl_bd_qp[k].begin(), gradVctrl_bd_qp[k].end(), 0);
        }

//========= gauss_loop boundary===============================================================
	    for(unsigned ig_bd=0; ig_bd < msh->_finiteElement[felt_bd][solVctrlType]->GetGaussPointNumber(); ig_bd++) {
		ml_prob._ml_msh->_finiteElement[felt_bd][solVctrlType]->JacobianSur(coordX_bd,ig_bd,weight_bd,phiVctrl_gss_bd,phiVctrl_x_gss_bd,normal);
		  
//========== temporary soln for surface gradient on a face parallel to the X axis ===================
		double dx_dxi = 0.;
		const elem_type_1D* myeltype = static_cast<const elem_type_1D*>(msh->_finiteElement[felt_bd][solVctrlType]);
		const double* myptr = myeltype->GetDPhiDXi(ig_bd);
		for (int inode = 0; inode < nve_bd; inode++) {
		      dx_dxi += myptr[inode] * coordX_bd[0][inode];
		}  
		for (int inode = 0; inode < nve_bd; inode++) {
                     for (int d = 0; d < dim; d++) {
                          if (d == 0 ) phiVctrl_x_gss_bd[inode + d*nve_bd] = myptr[inode]* (1./ dx_dxi);
                          else  phiVctrl_x_gss_bd[inode + d*nve_bd] = 0.;
                     }
                }
//========== temporary soln for surface gradient on a face parallel to the X axis ===================
		  
//========== compute gauss quantities on the boundary ===============================================
      for (unsigned i = 0; i < nDofsVctrl; i++) {
	for (unsigned  k = 0; k < dim; k++) {
		   for(int i_bd = 0; i_bd < nve_bd; i_bd++) {
		       unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
		       Vctrl_bd_qp[k] += phiVctrl_gss_bd[i_bd] * solVctrl[k][i_vol];
		       for(unsigned ivar2=0; ivar2<dim; ivar2++) {
			   gradVctrl_bd_qp[k][ivar2] += phiVctrl_x_gss_bd[i_bd + ivar2 * nve_bd] * solVctrl[k][i_vol]; 
		       }
		   }
	}
      }
 //end unknowns eval at gauss points ********************************
		      
		  
//========== compute gauss quantities on the boundary ================================================
      for (unsigned  k = 0; k < dim; k++) {
	 integral_beta	+= ((beta_val/2)*(Vctrl_bd_qp[k])*(Vctrl_bd_qp[k])*weight);
// 	integral_gamma	  += ((gamma_val/2)*(gradVctrl_bd_qp[k][0])*(gradVctrl_bd_qp[k][0])*weight);
      }
      for (unsigned  k = 0; k < dim; k++) {
	for (unsigned  j = 0; j < dim; j++) {	
		integral_gamma	  += ((gamma_val/2)*(gradVctrl_bd_qp[k][j])*(gradVctrl_bd_qp[k][j])*weight);
	}
      }


                }  //end ig_bdry loop
	  
             }    //end if control face
	 }  //end if boundary faces
      }  // loop over element faces //jface   
  } //end if control element flag

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {

//STATE#############################################################################	
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV_gss, phiV_x_gss, phiV_xx_gss);
	
	msh->_finiteElement[ielGeom][solVType  /*solVdes*/]->Jacobian(coordX, ig, weight, phiVdes_gss, phiVdes_x_gss, phiVdes_xx_gss);

    
      for (unsigned i = 0; i < nDofsV; i++) {
	 for (unsigned  k = 0; k < dim; k++) {
	   	V_gss[k] += solV[k][i] * phiV_gss[i];
		Vdes_gss[k] += solVdes[k]/*[i]*/ * phiVdes_gss[i];
		}
      }
	

      for (unsigned  k = 0; k < dim; k++) {
	 integral_target_alpha+=((alpha_val* target_flag/2 ) * (V_gss[k]  - Vdes_gss[k]) * (V_gss[k]  - Vdes_gss[k])*weight);
      }

      }// end gauss point loop
    } //end element loop  

    std::cout << "The value of the integral of target is " << std::setw(11) << std::setprecision(10) <<  integral_target_alpha << std::endl;
    std::cout << "The value of the integral of beta is " << std::setw(11) << std::setprecision(10) <<  integral_beta << std::endl;
    std::cout << "The value of the integral of gamma is " << std::setw(11) << std::setprecision(10) <<  integral_gamma << std::endl; 
    
    
    return integral_target_alpha + integral_beta /*+ integral_gamma*/ ; 
	  
  
}

