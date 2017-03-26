#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "LinearImplicitSystem.hpp"
#include "NumericVector.hpp"

using namespace femus;

#define ALPHA_CTRL 1.
#define BETA_CTRL  1.e-3
#define GAMMA_CTRL 1.e-3

int ElementTargetFlag(const std::vector<double> & elem_center) {

 //***** set target domain flag ********************************** 
  int target_flag = 1; //set 0 to 1 to get the entire domain
  
   if ( elem_center[0] < 0.5 + (1./16. + 1./64.)  + 1.e-5  && elem_center[0] > 0.5 - (1./16. + 1./64.) - 1.e-5  && 
        elem_center[1] < 0.5 + (1./16. + 1./64.)  + 1.e-5  && elem_center[1] > 0.5 - (1./16. + 1./64.) - 1.e-5 
  ) {
     
     target_flag = 1;
     
  }
  
     return target_flag;

}

// find volume elements that contain a control face element
int ControlDomainFlag(const std::vector<double> & elem_center) {

 //***** set flag ********************************** 

  double mesh_size = 1./32.;
  int control_el_flag = 0;
   if ( elem_center[1] >  1. - mesh_size ) { control_el_flag = 1; }

     return control_el_flag;

}


double DesiredTarget() {
 
  return 1.;
}

double InitialValueContReg(const std::vector < double >& x) {
  return ControlDomainFlag(x);
}

double InitialValueTargReg(const std::vector < double >& x) {
  return ElementTargetFlag(x);
}

double InitialValueState(const std::vector < double >& x) {
  return 0.;
}

double InitialValueAdjoint(const std::vector < double >& x) {
  return 0.;
}

double InitialValueControl(const std::vector < double >& x) {
  return 0.;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = true; //dirichlet
  value = 0;

  if(!strcmp(name,"control")) {
  if (faceName == 3)
    dirichlet = false;
  }
  
  return dirichlet;
}


double ComputeIntegral(MultiLevelProblem& ml_prob);

void AssembleOptSys(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  mlMsh.GenerateCoarseBoxMesh(32,32,0,0.,1.,0.,1.,0.,0.,QUAD9,"seventh");
   //1: bottom  //2: right  //3: top  //4: left
  
 /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("state", LAGRANGE, SECOND);
  mlSol.AddSolution("control", LAGRANGE, SECOND);
  mlSol.AddSolution("adjoint", LAGRANGE, SECOND);
  mlSol.AddSolution("TargReg",  DISCONTINOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
  mlSol.AddSolution("ContReg",  DISCONTINOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field

  
  mlSol.Initialize("All");    // initialize all varaibles to zero

  mlSol.Initialize("state", InitialValueState);
  mlSol.Initialize("control", InitialValueControl);
  mlSol.Initialize("adjoint", InitialValueAdjoint);
  mlSol.Initialize("TargReg", InitialValueTargReg);
  mlSol.Initialize("ContReg", InitialValueContReg);

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("state");
  mlSol.GenerateBdc("control");
  mlSol.GenerateBdc("adjoint");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  
 // add system  in mlProb as a Linear Implicit System
  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("LiftRestr");
 
  system.AddSolutionToSystemPDE("state");  
  system.AddSolutionToSystemPDE("control");  
  system.AddSolutionToSystemPDE("adjoint");  
  
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleOptSys);

  // initilaize and solve the system
  system.init();
  system.MGsolve();
  
  ComputeIntegral(mlProb);
 
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("state");
  variablesToBePrinted.push_back("control");
  variablesToBePrinted.push_back("adjoint");
  variablesToBePrinted.push_back("TargReg");
  variablesToBePrinted.push_back("ContReg");

  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  return 0;
}


void AssembleOptSys(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
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
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  vector < vector < double > > x(dim);
  vector < vector < double> >  x_bdry(dim);
  for (unsigned i = 0; i < dim; i++) {
         x[i].reserve(maxSize);
	 x_bdry[i].reserve(maxSize);
  }

 //*************************** 

 //*************************** 
  double weight = 0.; // gauss point weight
  double weight_bdry = 0.; // gauss point weight on the boundary


 //******** state ******************* 
 //*************************** 
  vector <double> phi_u;  // local test function
  vector <double> phi_u_x; // local test function first order partial derivatives
  vector <double> phi_u_xx; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
  unsigned solIndex_u    = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u     = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"
  unsigned solPdeIndex_u = mlPdeSys->GetSolPdeIndex("state");    // get the position of "state" in the pdeSys object

  vector < double >  sol_u;     sol_u.reserve(maxSize);
  vector< int > l2GMap_u;    l2GMap_u.reserve(maxSize);
 //*************************** 
 //*************************** 

  
 //************ adjoint **********
 //*******************************
  vector <double> phi_adj;  // local test function
  vector <double> phi_adj_x; // local test function first order partial derivatives
  vector <double> phi_adj_xx; // local test function second order partial derivatives

  phi_adj.reserve(maxSize);
  phi_adj_x.reserve(maxSize * dim);
  phi_adj_xx.reserve(maxSize * dim2);
 
  unsigned solIndex_adj    = mlSol->GetIndex("adjoint");    // get the position of "state" in the ml_sol object
  unsigned solType_adj     = mlSol->GetSolutionType(solIndex_adj);    // get the finite element type for "state"
  unsigned solPdeIndex_adj = mlPdeSys->GetSolPdeIndex("adjoint");    // get the position of "state" in the pdeSys object

  vector < double >  sol_adj;   sol_adj.reserve(maxSize);
  vector < int > l2GMap_adj; l2GMap_adj.reserve(maxSize);
  //*************************** 
 //*************************** 

  
 //************ bdry cont *********
 //******************************** 
  vector <double> phi_ctrl_bdry;  
  vector <double> phi_ctrl_x_bdry; 

  phi_ctrl_bdry.reserve(maxSize);
  phi_ctrl_x_bdry.reserve(maxSize * dim);
  
  unsigned solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

  unsigned solPdeIndex_ctrl = mlPdeSys->GetSolPdeIndex("control");

  std::string ctrl_name = "control";
  
 vector < double >  sol_ctrl;   sol_ctrl.reserve(maxSize);
 vector< int > l2GMap_ctrl;   l2GMap_ctrl.reserve(maxSize);
  //*************************** 
 //*************************** 
  

 //*************************** 
  //********* WHOLE SET OF VARIABLES ****************** 
  const int solType_max = 2;  //biquadratic

  const int n_vars = 3;
 
  vector< int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_vars*maxSize);
  
  vector< double > Res; // local redidual vector
  Res.reserve(n_vars*maxSize);

  vector < double > Jac;
  Jac.reserve( n_vars*maxSize * n_vars*maxSize);
 //*************************** 

  
 //********** DATA ***************** 
  double u_des = DesiredTarget();
  double alpha = ALPHA_CTRL;
  double beta  = BETA_CTRL;
  double gamma = GAMMA_CTRL;
  double penalty_strong = 1.e10;
 //*************************** 
  
  
  if (assembleMatrix)  KK->zero();

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned kelGeom = msh->GetElementType(iel);    // element geometry type

 //********* GEOMETRY ******************
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof // via local to global solution node

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
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
  //*************************************** 
  
  //***** set target domain flag ********************************** 
   int target_flag = 0;
   target_flag = ElementTargetFlag(elem_center);
  //*************************************** 
   
  //***** set control flag ********************************** 
  int control_el_flag = 0;
        control_el_flag = ControlDomainFlag(elem_center);
  std::vector<int> control_node_flag(nDofx,0);
//   if (control_el_flag == 0) std::fill(control_node_flag.begin(), control_node_flag.end(), 0);
  //*************************************** 
    
 //*********** state **********************
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
    sol_u   .resize(nDof_u);
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
      unsigned solDof_u = msh->GetSolutionDof(i,iel, solType_u);    // global to global mapping between solution node and solution dof
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);      // global extraction and local storage for the solution
      l2GMap_u[i] = pdeSys->GetSystemDof(solIndex_u, solPdeIndex_u, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
 //*********** state **************************** 

 //*********** bdry cont *************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);    // number of solution element dofs
    sol_ctrl    .resize(nDof_ctrl);
    l2GMap_ctrl.resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);    // global to global mapping between solution node and solution dof
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);      // global extraction and local storage for the solution
      l2GMap_ctrl[i] = pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);    // global to global mapping between solution node and pdeSys dof
    } 
 //*********************************************** 

 //*********** adjoint **************************
    unsigned nDof_adj  = msh->GetElementDofNumber(iel, solType_adj);    // number of solution element dofs
    sol_adj    .resize(nDof_adj);
    l2GMap_adj.resize(nDof_adj);
    for (unsigned i = 0; i < sol_adj.size(); i++) {
      unsigned solDof_adj = msh->GetSolutionDof(i, iel, solType_adj);    // global to global mapping between solution node and solution dof
      sol_adj[i] = (*sol->_Sol[solIndex_adj])(solDof_adj);      // global extraction and local storage for the solution
      l2GMap_adj[i] = pdeSys->GetSystemDof(solIndex_adj, solPdeIndex_adj, i, iel);    // global to global mapping between solution node and pdeSys dof
    } 
 //*********** adjoint **************************** 

 
 //********** ALL VARS ***************** 
    unsigned nDof_AllVars = nDof_u + nDof_adj + nDof_ctrl; 
    int nDof_max    =  nDof_u;   // AAAAAAAAAAAAAAAAAAAAAAAAAAA TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    if(nDof_adj > nDof_max) 
    {
      nDof_max = nDof_adj;
      }
    
    if(nDof_ctrl > nDof_max)
    {
      nDof_max = nDof_ctrl;
    }
    
    
    Res.resize(nDof_AllVars);
    std::fill(Res.begin(), Res.end(), 0.);

    Jac.resize(nDof_AllVars * nDof_AllVars);
    std::fill(Jac.begin(), Jac.end(), 0.);
    
    l2GMap_AllVars.resize(0);
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_u.begin(),l2GMap_u.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_ctrl.begin(),l2GMap_ctrl.end());
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_adj.begin(),l2GMap_adj.end());
 //*************************** 

    
 //===========================   

	// Perform face loop over elements that contain some control face
	if (control_el_flag == 1) {
	  
	  double tau=0.;
	  vector<double> normal(dim,0);
	       
	  // loop on faces of the current element

	  for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {
            std::vector < double > xx(3,0.);  //not being used, because the boundaries are identified by the face numbers
	    // look for boundary faces
	    if(el->GetFaceElementIndex(iel,jface) < 0) {
	      unsigned int face = -( msh->el->GetFaceElementIndex(iel,jface)+1);
	      
		
// 	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
	      if(  face == 3) { //control face

// ==================================
		//we use the dirichlet flag to say: if dirichlet = true, we set 1 on the diagonal. if dirichlet = false, we put the boundary equation
	      bool  dir_bool = mlSol->GetBdcFunction()(xx,ctrl_name.c_str(),tau,face,0.);

// ==================================

		
		unsigned nve = msh->GetElementFaceDofNumber(iel,jface,solType_ctrl);
		const unsigned felt = msh->GetElementFaceType(iel, jface);    
		for(unsigned i=0; i<nve; i++) {
		  unsigned int ilocal = msh->GetLocalFaceVertexIndex(iel, jface, i);
                  unsigned iDof = msh->GetSolutionDof(ilocal, iel, xType);
		  for(unsigned idim=0; idim<dim; idim++) {
		      x_bdry[idim][i]=(*msh->_topology->_Sol[idim])(iDof);
		  }
		}
		for(unsigned igs=0; igs < msh->_finiteElement[felt][solType_ctrl]->GetGaussPointNumber(); igs++) {
		  msh->_finiteElement[felt][solType_ctrl]->JacobianSur(x_bdry,igs,weight_bdry,phi_ctrl_bdry,phi_ctrl_x_bdry,normal);
		  //phi1 =msh->_finiteElement[felt][SolType2]->GetPhi(igs);
		  // *** phi_i loop ***
		  for(unsigned i=0; i<nve; i++) {
		    unsigned int ilocal = msh->GetLocalFaceVertexIndex(iel, jface, i);

	      if (dir_bool == false) { 
		std::cout << " found boundary control nodes ==== " << std::endl;
			for(unsigned k=0; k<control_node_flag.size(); k++) {
				  control_node_flag[ilocal] = 1;
			}
              }

		    

		    for(unsigned j=0; j<nve; j++) {
		    unsigned int jlocal = msh->GetLocalFaceVertexIndex(iel, jface, j);


// block delta_adjoint/control ========
// 		   if ( i < nDof_adj    && j < nDof_ctrl   &&  i==j)   Jac[ 
// 			(nDof_u + nDof_ctrl + i) * (nDof_u + nDof_ctrl + nDof_adj)  +
// 		        (nDof_u + j)             ]  += - (1 - control_node_flag[i]) * penalty_strong;      
		    
// block delta_control / control ========
   
//                  Jac[  
// 		    (nDof_u + ilocal) * (nDof_u + nDof_ctrl + nDof_adj) +
// 	            (nDof_u + jlocal) ] 
// 			+=   control_node_flag[ilocal] * weight_bdry* (alpha*phi_ctrl_bdry[ilocal]*phi_ctrl_bdry[jlocal]);
		    
		    double grad_bdry = 0.;
		      for (unsigned d = 0; d < dim; d++) {   grad_bdry += phi_ctrl_x_bdry[i+d*nve] * phi_ctrl_x_bdry[j+d*nve];    }
// 	         Jac[
// 		    (nDof_u + ilocal) * (nDof_u + nDof_ctrl + nDof_adj) +
// 	            (nDof_u + jlocal) ] 
// 	                += control_node_flag[ilocal] * weight_bdry * beta * grad_bdry;
		   }
				  
		  }
		}
	      }
	      
	    }
	  }    
	  
	} //end if control element flag
	
	else { //here we set the diagonal to 1 and the rhs to 0
	  
	  
	}
    
 //===========================   

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
        //  ==== state 
	msh->_finiteElement[kelGeom][solType_u]  ->Jacobian(x, ig, weight, phi_u, phi_u_x, phi_u_xx);
        //  ==== adj 
        msh->_finiteElement[kelGeom][solType_adj]->Jacobian(x, ig, weight, phi_adj, phi_adj_x, phi_adj_xx);
          
	
       //FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {

          // FIRST ROW - adj
          if (i < nDof_u)      Res[0                  + i] += weight * ( alpha * target_flag * u_des * phi_u[i] );
  
          // SECOND ROW - ctrl
           if (i < nDof_ctrl)  Res[nDof_u             + i] += penalty_strong * ( (1 - control_node_flag[i]) * x[1][i] + control_node_flag[i]*19. );//weight * (0.);
	      
	  // THIRD ROW - state
	  if (i < nDof_adj)    Res[nDof_u + nDof_ctrl + i] += weight * (0.) ;
	      
	      
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
              double laplace_mat_u = 0.;
              double laplace_mat_adj = 0.;

              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_adj && j < nDof_u )         laplace_mat_u      += (phi_adj_x   [i * dim + kdim] * phi_u_x   [j * dim + kdim]);
              if ( i < nDof_u   && j < nDof_adj )   laplace_mat_adj        += (phi_u_x     [i * dim + kdim] * phi_adj_x [j * dim + kdim]);
		
	      }

              //adjoint row ==================
              // BLOCK delta_state / state	      
              if ( i < nDof_u && j < nDof_u )   Jac[ i * (nDof_u + nDof_ctrl + nDof_adj) +
		                                         (0 + j)                      ]       += weight * 1. * target_flag *  phi_u[i] * phi_u[j];   
//               //BLOCK delta_state / adjoint
              if ( i < nDof_u && j < nDof_adj )   Jac[ i * (nDof_u + nDof_ctrl + nDof_adj) +
							 (nDof_u + nDof_ctrl + j)                 ]  += weight * laplace_mat_adj;
	      
	      
              //control row ==================
             //enforce control zero outside the control boundary
	      if ( i < nDof_ctrl && j < nDof_ctrl && i==j)
		Jac[    
		(nDof_u + i) * (nDof_u + nDof_ctrl + nDof_adj)  +
		(nDof_u + j)                                ]  += penalty_strong * ( (1 - control_node_flag[i]) + control_node_flag[i] )/*weight * phi_adj[i]*phi_adj[j]*/;
	      
              //state row ==================
              // BLOCK delta_adjoint / state
// 	      if ( i < nDof_adj && j < nDof_u  &&  i==j) 
// 		Jac[    
// 		(nDof_u + nDof_ctrl + i) * (nDof_u + nDof_ctrl + nDof_adj)  +
// 		(0 + j)                                ]  += (1 - control_node_flag[i])*penalty_strong;
	      
              // BLOCK delta_adjoint / state
	      if ( i < nDof_adj && j < nDof_u )
		Jac[    
		(nDof_u + nDof_ctrl + i) * (nDof_u + nDof_ctrl + nDof_adj)  +
		(0 + j)                                ]  += weight * laplace_mat_u;
	      
              // BLOCK delta_adjoint / adjoint
	      if ( i < nDof_adj && j < nDof_adj )
		Jac[    
		(nDof_u + nDof_ctrl + i) * (nDof_u + nDof_ctrl + nDof_adj)  +
		(nDof_u + nDof_ctrl + j)                                ]  += weight * phi_adj[i]*phi_adj[j];
	      

// 	      std::cout << Jac[ i * (nDof_u + nDof_ctrl + nDof_adj) +j ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + i) * (nDof_u + nDof_ctrl + nDof_adj) +j ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + nDof_ctrl + i) * (nDof_u + nDof_ctrl + nDof_adj) +j ] << " " << std::endl;
// 	      std::cout << Jac[ i * (nDof_u + nDof_ctrl + nDof_adj) + (nDof_u + j) ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + i) * (nDof_u + nDof_ctrl + nDof_adj) + (nDof_u + j) ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + nDof_ctrl + i) * (nDof_u + nDof_ctrl + nDof_adj) + (nDof_u + j) ] << " " << std::endl;
// 	      std::cout << Jac[ i * (nDof_u + nDof_ctrl + nDof_adj) + (nDof_u + nDof_ctrl + j) ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + i) * (nDof_u + nDof_ctrl + nDof_adj) + (nDof_u + nDof_ctrl + j) ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + nDof_ctrl + i) * (nDof_u + nDof_ctrl + nDof_adj) + (nDof_u + nDof_ctrl + j) ] << " " << std::endl;
	      
	    } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop

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

  return;
}


double ComputeIntegral(MultiLevelProblem& ml_prob)    {
  
  
  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*                     el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

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
  double weight; // gauss point weight
  

 //******** state ******************* 
 //*************************** 
  vector <double> phi_u;  // local test function
  vector <double> phi_u_x; // local test function first order partial derivatives
  vector <double> phi_u_xx; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * dim);
  phi_u_xx.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

  vector < double >  sol_u; // local solution
  sol_u.reserve(maxSize);
  
  double u_gss = 0.;
 //*************************** 
 //*************************** 

  
 //************ adjoint *************** 
 //*************************** 
  vector <double> phi_udes;  // local test function
  vector <double> phi_udes_x; // local test function first order partial derivatives
  vector <double> phi_udes_xx; // local test function second order partial derivatives

    phi_udes.reserve(maxSize);
    phi_udes_x.reserve(maxSize * dim);
    phi_udes_xx.reserve(maxSize * dim2);
 
  
//  unsigned solIndexTdes;
//   solIndexTdes = mlSol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solTypeTdes = mlSol->GetSolutionType(solIndexTdes);    // get the finite element type for "state"

  vector < double >  sol_udes; // local solution
  sol_udes.reserve(maxSize);
  vector< int > l2GMap_Tdes;
  l2GMap_Tdes.reserve(maxSize);
  double udes_gss = 0.;
  //*************************** 
 //*************************** 

 //************ cont *************** 

  unsigned solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

  vector< int > l2GMap_ctrl;
    l2GMap_ctrl.reserve(maxSize);
  
  double ctrl_gss = 0.;
  //*************************** 
 //*************************** 
  

 //*************************** 
  //********* WHOLE SET OF VARIABLES ****************** 
  const int solType_max = 2;  //biquadratic

  const int n_vars = 3;
 
  vector< int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_vars*maxSize);
  
  vector< double > Res; // local redidual vector
  Res.reserve(n_vars*maxSize);

  vector < double > Jac;
  Jac.reserve( n_vars*maxSize * n_vars*maxSize);
 //*************************** 

  
 //********** DATA ***************** 
  double u_des = DesiredTarget();
  //*************************** 
  
  double integral = 0.;

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned kelGeom = msh->GetElementType(iel);    // element geometry type
    
 //********* GEOMETRY ****************** 
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
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
  //*************************************** 
  
  //***** set target domain flag ********************************** 
   int target_flag = 0;
   target_flag = ElementTargetFlag(elem_center);
  //*************************************** 

   
 //*********** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
    sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
      unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);    // global to global mapping between solution node and solution dof
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);      // global extraction and local storage for the solution
    }
 //*********** state **************************** 


 //*********** cont **************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);    // number of solution element dofs
 //*********** Tcont **************************** 
 
 
 //*********** udes **************************** 
    unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
        sol_udes    .resize(nDof_udes);
    for (unsigned i = 0; i < sol_udes.size(); i++) {
            sol_udes[i] = u_des;  //dof value
    } 
 //*********** udes **************************** 

 
 //********** ALL VARS ***************** 
    int nDof_max    =  nDof_u;   // AAAAAAAAAAAAAAAAAAAAAAAAAAA TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    if(nDof_udes > nDof_max) 
    {
      nDof_max = nDof_udes;
      }
    
    if(nDof_ctrl > nDof_max)
    {
      nDof_max = nDof_ctrl;
    }
    
  //*************************** 

   
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
        //  ==== state 
	msh->_finiteElement[kelGeom][solType_u]   ->Jacobian(x, ig, weight, phi_u, phi_u_x, phi_u_xx);
        //  ==== adjoint
        msh->_finiteElement[kelGeom][solType_u/*solTypeTdes*/]->Jacobian(x, ig, weight, phi_udes, phi_udes_x, phi_udes_xx);

	u_gss = 0.;  for (unsigned i = 0; i < nDof_u; i++) u_gss += sol_u[i] * phi_u[i];		
	udes_gss  = 0.; for (unsigned i = 0; i < nDof_udes; i++)  udes_gss  += sol_udes[i]  * phi_udes[i];  

               integral += target_flag * weight * (u_gss  - udes_gss) * (u_gss - udes_gss);
	  
      } // end gauss point loop
      
  } //end element loop

  std::cout << "The value of the integral is " << std::setw(11) << std::setprecision(10) << integral << std::endl;
  
return integral;
  
}
  
  

