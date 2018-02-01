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

#define NSUB_X  32
#define NSUB_Y  32

#define CTRL_FACE_IDX  3


using namespace femus;

 double force[3] = {0./*1.*/,0.,0.};
 double Vel_desired[3] = {0.125,0.,0.};
 double alpha_val = 1.;
 double beta_val = 1.;
 double gamma_val = 1.;
 
  
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

double ComputeIntegral(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {
  
  FemusInit mpinit(argc, args, MPI_COMM_WORLD); 	// init Petsc-MPI communicator
  
//     // ======= Files ========================
//   Files files; 
//         files.CheckIODirectories();
// 	files.RedirectCout();
  
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
  
  
  

	
  char ordertobeprinted[30];
  int mix = sprintf(ordertobeprinted, "biquadratic alpha = %e beta = %e gamma = %e" , alpha_val,beta_val,gamma_val);
    // ==================================================
  

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
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);
  mlSol.AddSolution("P", LAGRANGE, FIRST);
  // adjoint =====================  
  mlSol.AddSolution("UADJ", LAGRANGE, SECOND);
  mlSol.AddSolution("VADJ", LAGRANGE, SECOND);
  if (dim == 3) mlSol.AddSolution("WADJ", LAGRANGE, SECOND);
  mlSol.AddSolution("PADJ", LAGRANGE, FIRST);
  // boundary condition =====================
  mlSol.AddSolution("GX", LAGRANGE, SECOND);
  mlSol.AddSolution("GY", LAGRANGE, SECOND);
  if (dim == 3) mlSol.AddSolution("GZ", LAGRANGE, SECOND);
  mlSol.AddSolution("THETA", LAGRANGE, SECOND);
  // control ===================== 
  
  
  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionOpt);
  mlSol.GenerateBdc("All");
  

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  mlProb.parameters.set<Fluid>("Fluid") = fluid;

  // add system Poisson in mlProb as a Linear Implicit System
//   NonLinearImplicitSystem& system_opt    = mlProb.add_system < NonLinearImplicitSystem > ("NSOpt");
  LinearImplicitSystem& system_opt    = mlProb.add_system < LinearImplicitSystem > ("NSOpt");

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
  system_opt.MLsolve();

    ComputeIntegral(mlProb);
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, ordertobeprinted /*"biquadratic"*/, variablesToBePrinted);
//   vtkIO.Write(files.GetOutputPath(),ordertobeprinted , variablesToBePrinted);

  

  return 0;
}


// nonAD is in the old PETSc, edit this for the new PETSc
void AssembleNavierStokesOpt(MultiLevelProblem& ml_prob){
     
  //pointers
  LinearImplicitSystem& mlPdeSys  = ml_prob.get_system<LinearImplicitSystem>("NSOpt");
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
   
   
  //boundary adjoint shape functions  
   
//   vector < double > phiV_bd_gss;
//   vector < double > phiV_x_bd_gss;
//   phiV_bd_gss.reserve(maxSize);
//   phiV_x_bd_gss.reserve(maxSize*dim);
  vector < vector < double > > phi_bd_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_bd_gss_fe(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_bd_gss_fe[fe].reserve(maxSize);
      phi_x_bd_gss_fe[fe].reserve(maxSize*dim);
   }
   
   vector <double> phi_adj_vol_at_bdry;  // local test function
  vector <double> phi_adj_x_vol_at_bdry; // local test function first order partial derivatives
  phi_adj_vol_at_bdry.reserve(maxSize);
  phi_adj_x_vol_at_bdry.reserve(maxSize * dim);
  vector <double> sol_adj_x_vol_at_bdry_gss;
  sol_adj_x_vol_at_bdry_gss.reserve(dim);
  
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
  vector < vector < double > > SolVAR_eldofs(n_unknowns);
  vector < vector < double > > gradSolVAR_eldofs(n_unknowns);
  
  for(int k=0; k<n_unknowns; k++) {
    SolVAR_eldofs[k].reserve(maxSize);
    gradSolVAR_eldofs[k].reserve(maxSize*dim);    
  }

  //------------ at quadrature points ---------------------
  vector < double > SolVAR_qp(n_unknowns);
    vector < vector < double > > gradSolVAR_qp(n_unknowns);
    for(int k=0; k<n_unknowns; k++) {  gradSolVAR_qp[k].resize(dim);  }
      
    
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

    unsigned nDofsVctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin]);    // number of solution element dofs
    unsigned nDofsPctrl = msh->GetElementDofNumber(iel,SolFEType[ctrl_pos_begin + press_type_pos] );    // number of solution element dofs

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
  std::vector<int> control_node_flag(nDofsX,0);
//   if (control_el_flag == 0) std::fill(control_node_flag.begin(), control_node_flag.end(), 0);
 //*************************************************** 
  
   //STATE###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	Sol_n_el_dofs[k]=ndofs_unk;
       SolVAR_eldofs[k].resize(ndofs_unk);
       JACDof[k].resize(ndofs_unk); 
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
       SolVAR_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
       JACDof[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
  //CTRL###################################################################

       
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
  
    //=============================================================================

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
		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
          std::vector<bool> dir_bool; dir_bool.resize(dim);
		  for(unsigned idim=0; idim<dim; idim++) {
		    dir_bool[idim] = mlSol->GetBdcFunction()(xyz_bdc,ctrl_name[idim].c_str(),tau,face,0.);
		  }
	  
 //=================================================== 
		
		unsigned nve_bdry = msh->GetElementFaceDofNumber(iel,jface, SolFEType[ctrl_pos_begin] ); //AAAAAAAAAAAAAAAAA
		const unsigned felt_bdry = msh->GetElementFaceType(iel, jface);    
		for(unsigned i=0; i < nve_bdry; i++) {
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
//========= initialize gauss quantities on the boundary ============================================
	
	
  for(unsigned ig_bd=0; ig_bd < msh->_finiteElement[felt_bdry][SolFEType[ctrl_pos_begin]]->GetGaussPointNumber(); ig_bd++) {
	for(int fe=0; fe < NFE_FAMS; fe++) {
	ml_prob._ml_msh->_finiteElement[felt_bdry][fe]->JacobianSur(coordX_bd,ig_bd,weight_bd,phi_bd_gss_fe[fe],phi_x_bd_gss_fe[fe],normal);
	}
	msh->_finiteElement[ielGeom][SolFEType[adj_pos_begin]]->ShapeAtBoundary(coordX,ig_bd,phi_adj_vol_at_bdry,phi_adj_x_vol_at_bdry);

		  
//========== temporary soln for surface gradient on a face parallel to the X axis ===================
		  double dx_dxi = 0.;
		 const elem_type_1D * myeltype = static_cast<const elem_type_1D*>(msh->_finiteElement[felt_bdry][SolFEType[ctrl_pos_begin]]);
		 const double * myptr = myeltype->GetDPhiDXi(ig_bd);
		      for (int inode = 0; inode < nve_bdry; inode++) dx_dxi += myptr[inode] * coordX_bd[0][inode];
  
		      for (int inode = 0; inode < nve_bdry; inode++) {
                            for (int d = 0; d < dim; d++) {
                              if (d==0 ) phi_x_bd_gss_fe[ctrl_pos_begin][inode + d*nve_bdry] = myptr[inode]* (1./ dx_dxi);
                              else  phi_x_bd_gss_fe[ctrl_pos_begin][inode + d*nve_bdry] = 0.;
                         }
                     }
//========== temporary soln for surface gradient on a face parallel to the X axis ===================
		  
//========== compute gauss quantities on the boundary ===============================================
	for(unsigned unk = 0; unk < n_unknowns; unk++) {
	  SolVAR_bd_qp[unk] = 0.;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_bd_qp[unk][ivar2] = 0.; 
	  }
	  
	  for(int i_bd = 0; i_bd < nve_bdry; i_bd++/*unsigned i = 0; i < Sol_n_el_dofs[unk]; i++*/) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
	    SolVAR_bd_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i_bd] * SolVAR_eldofs[unk][i_vol];
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR_bd_qp[unk][ivar2] += phi_x_bd_gss_fe[ SolFEType[unk] ][i_bd + ivar2 * nve_bdry] * SolVAR_eldofs[unk][i_vol]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************
		      
//=============== grad dot n for residual ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
           std::fill(sol_adj_x_vol_at_bdry_gss.begin(), sol_adj_x_vol_at_bdry_gss.end(), 0.);
		      for (int iv = 0; iv < nDofsVadj; iv++)  {
                            for (int d = 0; d < dim; d++) {
			      sol_adj_x_vol_at_bdry_gss[d] += SolVAR_eldofs[adj_pos_begin][iv] * phi_adj_x_vol_at_bdry[iv * dim + d];//notice that the convention of the orders x y z is different from vol to bdry
			    }
		      }  
		      
    double grad_dot_n_adj_res = 0.;
        for(unsigned d=0; d<dim; d++) {
	  grad_dot_n_adj_res += sol_adj_x_vol_at_bdry_gss[d]*normal[d];  
	}
//=============== grad dot n  for residual =========================================       
		  
//========== compute gauss quantities on the boundary ================================================

		  // *** phi_i loop ***
		  for(unsigned i_bdry=0; i_bdry < nve_bdry; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

                 double lap_rhs_dctrl_ctrl_bdry_gss_i = 0.;
                 for (unsigned d = 0; d < dim; d++) {
                       if ( i_vol < nDofsVctrl )  lap_rhs_dctrl_ctrl_bdry_gss_i +=  phi_x_bd_gss_fe[SolPdeIndex[ctrl_pos_begin]][i_bdry + d*nve_bdry] * gradSolVAR_bd_qp[ctrl_pos_begin][d];
                 }
                 
//=============== construct control node flag field on the go  =========================================    
	      /* (control_node_flag[i])       picks nodes on \Gamma_c
	         (1 - control_node_flag[i])   picks nodes on \Omega \setminus \Gamma_c
	       */
		  for(unsigned idim=0; idim<dim; idim++) {
	      if (dir_bool[idim] == false) { 
// 		std::cout << " found boundary control nodes ==== " << std::endl;
			for(unsigned k=0; k<control_node_flag.size(); k++) {
				  control_node_flag[i_vol] = 1;
			}
              }
		  }
//=============== construct control node flag field on the go  =========================================    
		vector <double>  SolU; SolU.resize(nDofsV);
		for (unsigned i = 0; i< SolU.size(); i++){
		SolU[i] = SolVAR_qp[SolPdeIndex[state_pos_begin]]; }
		vector <double>  Solg; Solg.resize(nDofsVctrl);
		for (unsigned i = 0; i< Solg.size(); i++){
		Solg[i] = SolVAR_qp[SolPdeIndex[ctrl_pos_begin]]; }
		
		  
//Boundary Residuals  and Jacobians ==================	
//============ delta_state row ============================================================================================
//   for (unsigned i_vol = 0; i_vol < nDofsV; i_vol++) {
	for (unsigned  kdim = 0; kdim < dim; kdim++) { // velocity block row 
	if(i_vol<nDofsV) Res[kdim][i_vol]  += - control_node_flag[i_vol] * (SolU[i_vol] - Solg[i_vol]);	    //u-g
	}		  
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
		    for(unsigned j_bdry=0; j_bdry < nve_bdry; j_bdry ++) {
		         unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);
//  for (unsigned j_vol = 0; j_vol < nDofsV; j_vol++) {
     for (unsigned  kdim = 0; kdim < dim; kdim++) { 
	if(i_vol==j_vol && i_vol<nDofsV) Jac[kdim][kdim][i_vol*nDofsV + j_vol]	+= control_node_flag[i_vol];  
      }
//BLOCK delta_state - control------------------------------------------------------------------------------------
     for (unsigned  kdim = 0; kdim < dim; kdim++) { 
	if(i_vol==j_vol && i_vol<nDofsV && j_vol<nDofsVctrl) Jac[kdim][kdim + ctrl_pos_begin][i_vol*nDofsV + j_vol]	+= (-1.) * control_node_flag[i_vol];  
      }

		      
		    }  //end j loop
      
//Boundary Residuals  and Jacobians ==================	

				}//end i loop		  
                	    }  //end ig_bdry loop
	  
                        }    //end if control face
	            }  //end if boundary faces
	  	  }  // loop over element faces   
               } //end if control element flag
    
    
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
		      double adv_res_uold_uold = 0.;
		      double adv_res_uold_uctrlold = 0.;

	      for (unsigned jdim = 0; jdim < dim; jdim++) {
		    lap_res_du_u += gradSolVAR_qp[SolPdeIndex[kdim]][jdim]*phi_x_gss_fe[SolFEType[kdim]][i * dim + jdim];
		   adv_res_uold_uold += SolVAR_qp[jdim] * gradSolVAR_qp[kdim][jdim];
		  adv_res_uold_uctrlold += SolVAR_qp[jdim] * gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim];

	      }      
	      Res[kdim][i]   +=  (         + force[kdim] * phi_gss_fe[SolFEType[kdim]][i]
                                           - IRe*lap_res_du_u 
//                                            - adv_res_uold_uold * phi_gss_fe[ SolFEType[kdim] ][i]
					    + SolVAR_qp[SolPdeIndex[press_type_pos]] * phi_x_gss_fe[SolFEType[kdim]][i * dim + kdim]) * weight; 
	}	    
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
	for (unsigned j = 0; j < nDofsV; j++) {
		      double lap_jac_du_u = 0.;
		      double adv_unew_uold = 0.;
		      double adv_uold_unew = 0.;
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
		    lap_jac_du_u += phi_x_gss_fe[SolFEType[kdim]][i * dim + kdim]*phi_x_gss_fe[SolFEType[kdim]][j * dim + kdim];
		    adv_uold_unew += SolVAR_qp[SolIndex[kdim]]*phi_x_gss_fe[ SolFEType[kdim] ][j * dim + kdim] * phi_gss_fe[ SolFEType[kdim] ][i];
		    adv_unew_uold += phi_gss_fe[ SolFEType[kdim] ][i] * gradSolVAR_qp[SolIndex[kdim]][kdim] * phi_gss_fe[ SolFEType[kdim] ][j];
	      }
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
		Jac[kdim][kdim][i*nDofsV + j] += (   IRe*lap_jac_du_u 
						    /*+ adv_uold_unew 
						    + adv_unew_uold*/) * weight; 
	      }
	} //j_du_u loop

//BLOCK delta_state - control------------------------------------------------------------------------------------
// 	for (unsigned j = 0; j < nDofsVctrl; j++) {
// 		      double lap_jac_du_ctrl = 0.;
// 	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
// 		    lap_jac_du_ctrl += phi_x_gss_fe[SolFEType[kdim]][i * dim + kdim]*phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j * dim + kdim];
// 	      }
// 	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
// 		Jac[kdim][kdim + ctrl_pos_begin ][i*nDofsVctrl + j] += ( IRe*lap_jac_du_ctrl ) * weight;
// 	      }
// 	} //j_du_ctrl loop

     
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
     for (unsigned kdim = 0; kdim < dim; kdim++) { 
		    double lap_res_dadj_adj = 0.;
		    double lap_res_dadj_u = 0.;
	   for (unsigned jdim = 0; jdim < dim; jdim++) {
		lap_res_dadj_adj += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]][jdim]*phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim + jdim];
		  lap_res_dadj_u += SolVAR_qp[SolPdeIndex[kdim]]*phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i];
	   }
	  Res[kdim + adj_pos_begin][i] += (   alpha_val*target_flag*(lap_res_dadj_u) 
					    - alpha_val*target_flag*Vel_desired[kdim] * phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]
					    - IRe*lap_res_dadj_adj
					    + SolVAR_qp[SolPdeIndex[press_type_pos + adj_pos_begin]] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim + kdim]) * weight;
      }
      
//BLOCK delta_adjoint - state------------------------------------------------------------------------------------------
     for (unsigned j = 0; j < nDofsV; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + adj_pos_begin][kdim][i*nDofsVadj + j] += ( -alpha_val*target_flag*phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]*phi_gss_fe[SolFEType[kdim]][j] ) * weight;
	  }
     }//j_dadj_u loop

//BLOCK delta_adjoint - control-----------------------------------------------------------------------------------------
//      for (unsigned j = 0; j < nDofsVctrl; j++) {
// 	  for (unsigned kdim = 0; kdim < dim; kdim++) {
// 	     Jac[kdim + adj_pos_begin][kdim + ctrl_pos_begin][i*nDofsVadj + j] += ( -alpha_val*target_flag*phi_gss_fe[SolFEType[kdim + adj_pos_begin]][i]*phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j] ) * weight;
// 	  }
//      }//j_dadj_ctrl loop

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
	    Jac[press_type_pos + adj_pos_begin][kdim + adj_pos_begin][i*nDofsVadj + j] += - ( phi_gss_fe[SolFEType[press_type_pos + adj_pos_begin]][i] * phi_x_gss_fe[SolFEType[kdim + adj_pos_begin]][i * dim + kdim] ) * weight;
	  }
      }//j loop
  }//i_div_adj

      //============ delta_adjoint row =============================================================================================


//============ delta_control row ==================================================================================================
// THIRD ROW
  for (unsigned i = 0; i < nDofsVctrl; i++) {
      for (unsigned kdim = 0; kdim < dim; kdim++) { 
		    double lap_res_dctrl_ctrl = 0.;
		    double lap_res_dctrl_adj = 0.;
      for (unsigned jdim = 0; jdim < dim; jdim++) {
	  lap_res_dctrl_ctrl += gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim]*phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim + jdim];
	   lap_res_dctrl_adj += gradSolVAR_qp[SolPdeIndex[kdim + adj_pos_begin]][jdim]*phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim + jdim];
      }
      Res[kdim + ctrl_pos_begin][i] += (+ beta_val* SolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]
					+ gamma_val * lap_res_dctrl_ctrl
					) * weight;
      }

//BLOCK delta_control - state------------------------------------------------------------------------------------------------
      for (unsigned j = 0; j < nDofsV; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + ctrl_pos_begin][kdim][i*nDofsVctrl + j] += 0.;
	  }
      }//j_dctrl_u loop
      
//BLOCK delta_control - adjoint------------------------------------------------------------------------------------------------
      for (unsigned j = 0; j < nDofsVadj; j++) {
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + ctrl_pos_begin][kdim + adj_pos_begin][i*nDofsVctrl + j] += 0.;
	  }
      }//j_dctrl_adj loop

//DIAG BLOCK delta_control - control--------------------------------------------------------------------------------------
      for (unsigned j = 0; j < nDofsVctrl; j++) {
		      double lap_jac_dctrl_ctrl = 0.;
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
		lap_jac_dctrl_ctrl += phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim + kdim]*phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j * dim + kdim];
	  }
	  for (unsigned kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsVctrl + j] += (  - beta_val * phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i]*phi_gss_fe[SolFEType[kdim + ctrl_pos_begin]][j]
											- gamma_val*lap_jac_dctrl_ctrl ) * weight;
	  }
      }//j_dctrl_ctrl loop

// // //BLOCK Pressure_ctrl
// //       for (unsigned j = 0; j < nDofsPctrl; j++) {
// // 	  for (unsigned kdim = 0; kdim < dim; kdim++) {
// // 	      Jac[kdim + ctrl_pos_begin][press_type_pos + ctrl_pos_begin][i*nDofsPctrl + j] += -( phi_gss_fe[SolFEType[press_type_pos + ctrl_pos_begin]][j] * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim + kdim] ) * weight;
// // 	  }
// //       }//j_press_ctrl
  }//i_ctrl loop

// // //DIV_ctrl
// //   for (unsigned i = 0; i < nDofsPctrl; i++) {
// // 		  double div_ctrl_dctrl_qp = 0.;
// // 	  for (unsigned kdim = 0; kdim < dim; kdim++) {
// // 		div_ctrl_dctrl_qp += gradSolVAR_qp[SolPdeIndex[kdim + ctrl_pos_begin]][kdim] ;
// // 	  }
// // 	  Res[press_type_pos + ctrl_pos_begin][i] += ( (div_ctrl_dctrl_qp) * phi_gss_fe[SolFEType[press_type_pos + ctrl_pos_begin]][i] ) * weight;
// // 	  for (unsigned j = 0; j < nDofsVctrl; j++) {
// // 	      for (unsigned kdim = 0; kdim < dim; kdim++) {
// // 		Jac[press_type_pos + ctrl_pos_begin][kdim + ctrl_pos_begin][i*nDofsVctrl + j] += - ( phi_gss_fe[SolFEType[press_type_pos + ctrl_pos_begin]][i] * phi_x_gss_fe[SolFEType[kdim + ctrl_pos_begin]][i * dim + kdim] ) * weight;
// // 	      }
// // 	  }//j loop
// //   }//i_div_ctrl
 
   //============ delta_control row ==================================================================================================
 
      }  // end gauss point loop
 

      //***************************************************************************************************************

    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned i_unk=0; i_unk < n_unknowns; i_unk++) {
      RES->add_vector_blocked(Res[SolPdeIndex[i_unk]],JACDof[i_unk]);
        for(unsigned j_unk=0; j_unk < n_unknowns; j_unk++) {
	  if(assembleMatrix) JAC->add_matrix_blocked( Jac[ SolPdeIndex[i_unk] ][ SolPdeIndex[j_unk] ], JACDof[i_unk], JACDof[j_unk]);
        }
    }
 
   //--------------------------------------------------------------------------------------------------------  
  } //end list of elements loop for each subdomain
  
  
  JAC->close();
  RES->close();
  // ***************** END ASSEMBLY *******************
}



double ComputeIntegral(MultiLevelProblem& ml_prob) {


   LinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<LinearImplicitSystem> ("NSOpt");   // pointer to the linear implicit system named "Poisson"
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
  

//CONTROL######################################################################
  //velocity *******************************
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

  
  vector <double> phiVctrl_gss;  // local test function
  vector <double> phiVctrl_x_gss; // local test function first order partial derivatives
  vector <double> phiVctrl_xx_gss; // local test function second order partial derivatives

  phiVctrl_gss.reserve(maxSize);
  phiVctrl_x_gss.reserve(maxSize * dim);
  phiVctrl_xx_gss.reserve(maxSize * dim2);
  
  
  //velocity *******************************
   

//CONTROL######################################################################

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


      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {

//STATE#############################################################################	
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV_gss, phiV_x_gss, phiV_xx_gss);
	
	msh->_finiteElement[ielGeom][solVctrlType]->Jacobian(coordX, ig, weight, phiVctrl_gss, phiVctrl_x_gss, phiVctrl_xx_gss);

	msh->_finiteElement[ielGeom][solVType  /*solVdes*/]->Jacobian(coordX, ig, weight, phiVdes_gss, phiVdes_x_gss, phiVdes_xx_gss);

	  vector < vector < double > > gradVctrl_gss(dim);
      for (unsigned  k = 0; k < dim; k++) {
          gradVctrl_gss[k].resize(dim);
          std::fill(gradVctrl_gss[k].begin(), gradVctrl_gss[k].end(), 0);
        }
	
//     for (unsigned  k = 0; k < dim; k++) {
//       V_gss[k]       = 0.;
//       Vdes_gss[k]    = 0.;
//        Vctrl_gss[k]  = 0.;
//     }
    
      for (unsigned i = 0; i < nDofsV; i++) {
	 for (unsigned  k = 0; k < dim; k++) {
	   	V_gss[k] += solV[k][i] * phiV_gss[i];
		Vdes_gss[k] += solVdes[k]/*[i]*/ * phiVdes_gss[i];
		}
      }
	
      for (unsigned i = 0; i < nDofsVctrl; i++) {
	 for (unsigned  k = 0; k < dim; k++) {
	   Vctrl_gss[k] += solVctrl[k][i] * phiVctrl_gss[i];
	 }
     for (unsigned j = 0; j < dim; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              gradVctrl_gss[k][j] += phiVctrl_x_gss[i * dim + j] * solVctrl[k][i];
            }
          }
      }
          
          
	
      for (unsigned  k = 0; k < dim; k++) {
// 	for (unsigned  j = 0; j < dim; k++) {
	 integral_target_alpha/* integral[k]*/ +=((alpha_val* target_flag/2 ) * (V_gss[k] + Vctrl_gss[k] - Vdes_gss[k]) * (V_gss[k] + Vctrl_gss[k] - Vdes_gss[k])*weight)
// 					  + ((beta_val/2)*(Vctrl_gss[k])*(Vctrl_gss[k])*weight)
					 /* + ((gamma_val/2)*(gradVctrl_gss[k][j])*(gradVctrl_gss[k][j])*weight)*/;
	 integral_beta	+= ((beta_val/2)*(Vctrl_gss[k])*(Vctrl_gss[k])*weight);
// 	}
      }
      for (unsigned  k = 0; k < dim; k++) {
	for (unsigned  j = 0; j < dim; j++) {	
		integral_gamma	  += ((gamma_val/2)*(gradVctrl_gss[k][j])*(gradVctrl_gss[k][j])*weight);
	}
      }
      
      
//       integralval= sqrt(((integral[0]*integral[0]) +(integral[1]*integral[1]))/**weight*/);
// // 

      }// end gauss point loop
    } //end element loop  

    std::cout << "The value of the integral of target is " << std::setw(11) << std::setprecision(10) <<  integral_target_alpha << std::endl;
    std::cout << "The value of the integral of beta is " << std::setw(11) << std::setprecision(10) <<  integral_beta << std::endl;
    std::cout << "The value of the integral of gamma is " << std::setw(11) << std::setprecision(10) <<  integral_gamma << std::endl; 
    
    
    return integral_target_alpha + integral_beta /*+ integral_gamma*/ ; 
	  
  
}

