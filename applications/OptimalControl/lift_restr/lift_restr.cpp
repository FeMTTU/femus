#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "LinearImplicitSystem.hpp"
#include "NumericVector.hpp"

using namespace femus;

#define NSUB_X  16
#define NSUB_Y  16
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

int ControlDomainFlag(const std::vector<double> & elem_center) {

 //***** set target domain flag ********************************** 
  // flag = 1: we are in the lifting nonzero domain
  int control_el_flag = 0;
   if ( elem_center[1] >  0.7 ) { control_el_flag = 1; }

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

void AssembleLiftRestrProblem(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  mlMsh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,"seventh");
 /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("state", LAGRANGE, FIRST);
  mlSol.AddSolution("control", LAGRANGE, FIRST);
  mlSol.AddSolution("adjoint", LAGRANGE, FIRST);
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
  system.SetAssembleFunction(AssembleLiftRestrProblem);

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

    // ******* Print solution *******
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

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
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object

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
  double weight; // gauss point weight
  

 //******** state ******************* 
 //*************************** 
  vector <double> phi_u;  // local test function
  vector <double> phi_x_u; // local test function first order partial derivatives
  vector <double> phi_xx_u; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_x_u.reserve(maxSize * dim);
  phi_xx_u.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

  unsigned solPdeIndexThom;
  solPdeIndexThom = mlPdeSys->GetSolPdeIndex("state");    // get the position of "state" in the pdeSys object

  vector < double >  sol_u; // local solution
  sol_u.reserve(maxSize);
  vector< int > l2GMap_u;
  l2GMap_u.reserve(maxSize);
 //*************************** 
 //*************************** 

  
 //************ control *************** 
 //*************************** 
  vector <double> phi_ctrl;  // local test function
  vector <double> phi_x_ctrl; // local test function first order partial derivatives
  vector <double> phi_xx_ctrl; // local test function second order partial derivatives

  phi_ctrl.reserve(maxSize);
  phi_x_ctrl.reserve(maxSize * dim);
  phi_xx_ctrl.reserve(maxSize * dim2);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

  unsigned solPdeIndex_ctrl;
  solPdeIndex_ctrl = mlPdeSys->GetSolPdeIndex("control");

  vector < double >  sol_ctrl; // local solution
  sol_ctrl.reserve(maxSize);
 vector< int > l2GMap_ctrl;
  l2GMap_ctrl.reserve(maxSize);
  //*************************** 
 //*************************** 
  
  
 //************ adjoint *************** 
 //*************************** 
   vector <double> phi_adj;  // local test function
  vector <double> phi_x_adj; // local test function first order partial derivatives
  vector <double> phi_xx_adj; // local test function second order partial derivatives

  phi_adj.reserve(maxSize);
    phi_x_adj.reserve(maxSize * dim);
  phi_xx_adj.reserve(maxSize * dim2);
 
  
 unsigned solIndex_adj;
    solIndex_adj = mlSol->GetIndex("adjoint");    // get the position of "state" in the ml_sol object
  unsigned solType_adj = mlSol->GetSolutionType(solIndex_adj);    // get the finite element type for "state"

  unsigned solPdeIndex_adj;
  solPdeIndex_adj = mlPdeSys->GetSolPdeIndex("adjoint");    // get the position of "state" in the pdeSys object

  vector < double >  sol_adj; // local solution
    sol_adj.reserve(maxSize);
 vector< int > l2GMap_adj;
    l2GMap_adj.reserve(maxSize);
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
  double T_des = DesiredTarget();
  double alpha = ALPHA_CTRL;
  double beta  = BETA_CTRL;
  double gamma = GAMMA_CTRL;
  double penalty_strong = 10e+14;
 //*************************** 
  
  
  if (assembleMatrix)  KK->zero();

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned kelGeom = msh->GetElementType(iel);    // element geometry type

 //********* GEOMETRY ****************** 
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);  // global to global mapping between coordinates node and coordinate dof

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
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);  // global to global mapping between solution node and solution dof
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);            // global extraction and local storage for the solution
      l2GMap_u[i] = pdeSys->GetSystemDof(solIndex_u, solPdeIndexThom, i, iel);  // global to global mapping between solution node and pdeSys dof
    }
 //*********** **************************** 

 //*********** control **************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);    // number of solution element dofs
    sol_ctrl    .resize(nDof_ctrl);
    l2GMap_ctrl.resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);    // global to global mapping between solution node and solution dof
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);      // global extraction and local storage for the solution
      l2GMap_ctrl[i] = pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);    // global to global mapping between solution node and pdeSys dof
    } 
 //*************************************** 
 

 //*********** adj **************************** 
    unsigned nDof_adj  = msh->GetElementDofNumber(iel, solType_adj);    // number of solution element dofs
        sol_adj    .resize(nDof_adj);
        l2GMap_adj.resize(nDof_adj);
    for (unsigned i = 0; i < sol_adj.size(); i++) {
      unsigned solDof_adj = msh->GetSolutionDof(i, iel, solType_adj);   // global to global mapping between solution node and solution dof
      sol_adj[i] = (*sol->_Sol[solIndex_adj])(solDof_adj);      // global extraction and local storage for the solution
      l2GMap_adj[i] = pdeSys->GetSystemDof(solIndex_adj, solPdeIndex_adj, i, iel);   // global to global mapping between solution node and pdeSys dof
    } 
 //*************************************** 

 //********** ALL VARS ***************** 
    unsigned nDof_AllVars = nDof_u + nDof_ctrl + nDof_adj; 
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
    
 //***** set control flag ********************************** 
  int control_el_flag = 0;
  control_el_flag = ControlDomainFlag(elem_center);
  std::vector<int> control_node_flag(nDof_ctrl,0);
  if (control_el_flag == 1) std::fill(control_node_flag.begin(), control_node_flag.end(), 1);
  //*************************************** 

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
	msh->_finiteElement[kelGeom][solType_u]   ->Jacobian(x, ig, weight, phi_u, phi_x_u, phi_xx_u);
        msh->_finiteElement[kelGeom][solType_ctrl]->Jacobian(x, ig, weight, phi_ctrl, phi_x_ctrl, phi_xx_ctrl);
        msh->_finiteElement[kelGeom][solType_adj] ->Jacobian(x, ig, weight, phi_adj, phi_x_adj, phi_xx_adj);
	
        //FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {

          // FIRST ROW
	  if (i < nDof_u)    Res[0                   + i] += weight * ( alpha * target_flag * T_des * phi_u[i] );
          // SECOND ROW
	  if (i < nDof_ctrl)  {
	      if ( control_el_flag == 1)       Res[nDof_u + i] += /* (1 - control_node_flag[i]) **/ weight * ( alpha * target_flag * T_des * phi_ctrl[i] );
	      else if ( control_el_flag == 0)  Res[nDof_u + i] +=  /*control_node_flag[i] **/ penalty_strong * 0.; 
	  }
          // THIRD ROW
          if (i < nDof_adj) Res[nDof_u + nDof_ctrl + i] += weight * (0.) ;
	      
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
              double laplace_mat_du_u = 0.;
              double laplace_mat_du_adj = 0.;
              double laplace_mat_dadj_u = 0.;
              double laplace_mat_dctrl_adj = 0.;
              double laplace_mat_dadj_ctrl = 0.;
              double laplace_mat_dctrl_ctrl = 0.;

              for (unsigned kdim = 0; kdim < dim; kdim++) {
              if ( i < nDof_u && j < nDof_u )           laplace_mat_du_u           += (phi_x_u   [i * dim + kdim] * phi_x_u   [j * dim + kdim]);
              if ( i < nDof_u && j < nDof_adj )         laplace_mat_du_adj         += (phi_x_u   [i * dim + kdim] * phi_x_adj [j * dim + kdim]);
              if ( i < nDof_adj && j < nDof_u )         laplace_mat_dadj_u         += (phi_x_adj [i * dim + kdim] * phi_x_u   [j * dim + kdim]);  //equal to the previous
              if ( i < nDof_ctrl && j < nDof_adj )      laplace_mat_dctrl_adj      += (phi_x_ctrl[i * dim + kdim] * phi_x_adj [j * dim + kdim]);
              if ( i < nDof_adj && j < nDof_ctrl )      laplace_mat_dadj_ctrl      += (phi_x_adj [i * dim + kdim] * phi_x_ctrl[j * dim + kdim]);  //equal to the previous
              if ( i < nDof_ctrl   && j < nDof_ctrl )   laplace_mat_dctrl_ctrl     += (phi_x_ctrl  [i * dim + kdim] * phi_x_ctrl  [j * dim + kdim]);
	      }

              //first row ==================
              //DIAG BLOCK delta_state - state
	      if ( i < nDof_u && j < nDof_u )       Jac[    (0 + i) * nDof_AllVars    +
		                                            (0 + j)                      ]  += weight * alpha * target_flag * phi_u[j] *  phi_u[i];
              // BLOCK  delta_state - control
              if ( i < nDof_u && j < nDof_ctrl )   Jac[  (0 + i) * nDof_AllVars   +
                                                            (nDof_u + j)                 ]  += weight * alpha * target_flag  * phi_ctrl[j] *  phi_u[i];
	      
              // BLOCK  delta_state - adjoint
              if ( i < nDof_u && j < nDof_adj )   Jac[  (0 + i) * nDof_AllVars   +
                                                        (nDof_u + nDof_ctrl + j)         ]  += weight * laplace_mat_du_adj;
              //second row ==================	      
	      if ( control_el_flag == 1)  {
	      
              //BLOCK delta_control - control
              if ( i < nDof_ctrl   && j < nDof_ctrl   ) Jac[ (nDof_u + i) * nDof_AllVars +
		                                             (nDof_u + j)                     ]  += /*(1 - control_node_flag[i]) **/ weight * ( gamma * control_el_flag  * laplace_mat_dctrl_ctrl + beta * control_el_flag * phi_ctrl[i] * phi_ctrl[j] + alpha  * target_flag  * phi_ctrl[i] * phi_ctrl[j]);
              //BLOCK delta_control - state
              if ( i < nDof_ctrl   && j < nDof_u   ) Jac[ (nDof_u + i) * nDof_AllVars  +
								(0 + j)                                           ]  += /*(1 - control_node_flag[i]) **/ weight * alpha * target_flag * phi_u[j] * phi_ctrl[i];
	      //BLOCK delta_control - adjoint
              if ( i < nDof_ctrl   && j < nDof_adj  ) Jac[ (nDof_u + i) * nDof_AllVars  + 
		                                           (nDof_u + nDof_ctrl + j) ]  +=  /*(1 - control_node_flag[i]) **/ weight * laplace_mat_dctrl_adj;
	      }
	      
	      else if ( control_el_flag == 0)  {  
		
              //BLOCK delta_control - control
               if ( i < nDof_ctrl   && j < nDof_ctrl &&  i==j ) {
		Jac[ (nDof_u + i) * nDof_AllVars +
		     (nDof_u + j)                     ] += /*control_node_flag[i] **/ penalty_strong;
		}
	      
	   }
	      
            //third row ==================
              // BLOCK delta_adjoint - state	      
              if ( i < nDof_adj && j < nDof_u )   Jac[    (nDof_u + nDof_ctrl + i) * nDof_AllVars +
		                                          (0 + j)                      ]                       += weight * laplace_mat_dadj_u;   
	      
              // BLOCK delta_adjoint - control   
              if ( i < nDof_adj && j < nDof_ctrl )   Jac[    (nDof_u + nDof_ctrl + i)  * nDof_AllVars +
		                                             (nDof_u  + j)                      ]  += weight * laplace_mat_dadj_ctrl; 

	      
            } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop

      
      
    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector
	if (control_el_flag == 0) {  //elements that should have zero control
	  
         for (unsigned i = 0; i < nDof_max; i++) {
            for (unsigned j = 0; j < nDof_max; j++) {
// 	      std::cout << Jac[ i * nDof_AllVars +j ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + i) * nDof_AllVars +j ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + nDof_ctrl + i) * nDof_AllVars +j ] << " " << std::endl;
// 	      std::cout << Jac[ i * nDof_AllVars + (nDof_u + j) ] << " " << std::endl;
// 	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + j) ];
	      std::cout << " " << std::setfill(' ') << std::setw(10) << Jac[ (nDof_u + nDof_adj + i) * nDof_AllVars + (nDof_u + nDof_adj + j) ];
// 	      std::cout << Jac[ (nDof_u + nDof_ctrl + i) * nDof_AllVars + (nDof_u + j) ] << " " << std::endl;
// 	      std::cout << Jac[ i * nDof_AllVars + (nDof_u + nDof_ctrl + j) ] << " " << std::endl;
//       std::cout << Jac[ (nDof_u + i) * nDof_AllVars + (nDof_u + nDof_ctrl + j) ] << " " << std::endl;
// 	      std::cout << Jac[ (nDof_u + nDof_ctrl + i) * nDof_AllVars + (nDof_u + nDof_ctrl + j) ] << " " << std::endl;
	    }
	      std::cout << std::endl;
	 }
	 
	}
    std::cout << " ========== " << iel << " ================== " << std::endl;     
    
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
  
 //*************************** 
  double alpha = ALPHA_CTRL;
  double beta  = BETA_CTRL;
  double gamma = GAMMA_CTRL;

 //******** state ******************* 
 //*************************** 
  vector <double> phi_u;  // local test function
  vector <double> phi_x_u; // local test function first order partial derivatives
  vector <double> phi_xx_u; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_x_u.reserve(maxSize * dim);
  phi_xx_u.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u;
    solIndex_u = mlSol->GetIndex("state");    // get the position of "state" in the ml_sol object
  unsigned solType_u = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "state"

  vector < double >  sol_u; // local solution
    sol_u.reserve(maxSize);
  
  double Thom_gss = 0.;
 //*************************** 
 //*************************** 

 //************ control *************** 
 //*************************** 
  vector <double> phi_ctrl;  // local test function
  vector <double> phi_x_ctrl; // local test function first order partial derivatives
  vector <double> phi_xx_ctrl; // local test function second order partial derivatives

  phi_ctrl.reserve(maxSize);
  phi_x_ctrl.reserve(maxSize * dim);
  phi_xx_ctrl.reserve(maxSize * dim2);
  
  unsigned solIndex_ctrl;
  solIndex_ctrl = mlSol->GetIndex("control");
  unsigned solType_ctrl = mlSol->GetSolutionType(solIndex_ctrl);

  vector < double >  sol_ctrl; // local solution
  sol_ctrl.reserve(maxSize);
 vector< int > l2GMap_ctrl;
  l2GMap_ctrl.reserve(maxSize);
  
  double Tcont_gss = 0.;
  double Tcontgrad_gss = 0.;
  //*************************** 
 //*************************** 

  
 //************ desired *************** 
 //*************************** 
   vector <double> phi_Tdes;  // local test function
  vector <double> phi_x_Tdes; // local test function first order partial derivatives
  vector <double> phi_xx_Tdes; // local test function second order partial derivatives

    phi_Tdes.reserve(maxSize);
    phi_x_Tdes.reserve(maxSize * dim);
    phi_xx_Tdes.reserve(maxSize * dim2);
 
  
//  unsigned solIndexTdes;
//   solIndexTdes = mlSol->GetIndex("Tdes");    // get the position of "state" in the ml_sol object
//   unsigned solTypeTdes = mlSol->GetSolutionType(solIndexTdes);    // get the finite element type for "state"

  vector < double >  solTdes; // local solution
  solTdes.reserve(maxSize);
 vector< int > l2GMap_Tdes;
    l2GMap_Tdes.reserve(maxSize);
  double Tdes_gss = 0.;
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
  double T_des = DesiredTarget();
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
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

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
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);  // global to global mapping between solution node and solution dof
            sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);            // global extraction and local storage for the solution
    }
 //*************************************** 


 //*********** control **************************** 
    unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);    // number of solution element dofs
    sol_ctrl    .resize(nDof_ctrl);
    for (unsigned i = 0; i < sol_ctrl.size(); i++) {
      unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl);    // global to global mapping between solution node and solution dof
      sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);      // global extraction and local storage for the solution
    } 
 //*************************************** 
 
 
 //*********** u_des **************************** 
    unsigned nDofTdes  = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
    solTdes    .resize(nDofTdes);
    for (unsigned i = 0; i < solTdes.size(); i++) {
      solTdes[i] = T_des;  //dof value
    } 
 //*************************************** 

 
 //********** ALL VARS ***************** 
    int nDof_max    =  nDof_u;   // AAAAAAAAAAAAAAAAAAAAAAAAAAA TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    if(nDofTdes > nDof_max) 
    {
      nDof_max = nDofTdes;
      }
    
    if(nDof_ctrl > nDof_max)
    {
      nDof_max = nDof_ctrl;
    }
    
  //*************************** 
   
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
        //  ==== Thom 
	msh->_finiteElement[kelGeom][solType_u]   ->Jacobian(x, ig, weight, phi_u, phi_x_u, phi_xx_u);
        //  ==== ThomAdj 
        msh->_finiteElement[kelGeom][solType_u/*solTypeTdes*/]->Jacobian(x, ig, weight, phi_Tdes, phi_x_Tdes, phi_xx_Tdes);
        //  ==== ThomCont 
        msh->_finiteElement[kelGeom][solType_ctrl]  ->Jacobian(x, ig, weight, phi_ctrl, phi_x_ctrl, phi_xx_ctrl);

	Thom_gss = 0.;  for (unsigned i = 0; i < nDof_u; i++) Thom_gss += sol_u[i] * phi_u[i];		
	Tcont_gss = 0.; for (unsigned i = 0; i < nDof_ctrl; i++) Tcont_gss += sol_ctrl[i] * phi_ctrl[i];  
	Tdes_gss  = 0.; for (unsigned i = 0; i < nDofTdes; i++)  Tdes_gss  += solTdes[i]  * phi_Tdes[i]; 
    Tcontgrad_gss  = 0.; for (unsigned i = 0; i < nDof_ctrl; i++)  
    {
       for (unsigned idim = 0; idim < dim; idim ++) Tcontgrad_gss  += sol_ctrl[i] * phi_x_ctrl[i + idim * nDof_ctrl];
    }

               integral += target_flag * (alpha * weight * (Thom_gss +  Tcont_gss - Tdes_gss) * (Thom_gss +  Tcont_gss - Tdes_gss)
                                        + beta * weight * Tcont_gss * Tcont_gss 
                                        + gamma * weight * Tcontgrad_gss * Tcontgrad_gss);
	  
      } // end gauss point loop
  } //end element loop

  std::cout << "The value of the integral is " << std::setw(11) << std::setprecision(10) << integral << std::endl;
  
return integral;
  
}
  
  

