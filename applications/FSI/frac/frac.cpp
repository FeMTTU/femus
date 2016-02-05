
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

using namespace femus;

	double force[3] = {1.,0.,0.};



bool SetBoundaryCondition(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  
  bool dirichlet = true; //dirichlet
  value = 0.;
  
      if (facename == 2) {
       if (!strcmp(SolName, "V")) {   value = 0.;  } 
  else if (!strcmp(SolName, "U")) {  dirichlet = false; }
      }
      
      if (facename == 4) {
       if (!strcmp(SolName, "V")) {   value = 0.;  } 
  else if (!strcmp(SolName, "U")) {  dirichlet = false; }
      }
      
  return dirichlet;
}

void AssembleNavierStokes(MultiLevelProblem &ml_prob);

void AssembleNavierStokes_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );


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
 
//   MultiLevelMesh mlMsh;
//  mlMsh.ReadCoarseMesh(infile.c_str(),"seventh",Lref);
    mlMsh.GenerateCoarseBoxMesh( 8,8,0,-0.5,0.5,-0.5,0.5,0.,0.,QUAD9,"seventh");
    
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
     probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 3; 
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);

  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);

  mlSol.AddSolution("P", LAGRANGE, FIRST);
  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

// *** apparently needed by non-AD assemble only **********************
  // add fluid material
  Parameter parameter(Lref,Uref);
  
  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,0.001,1,"Newtonian",0.001,1.);
  std::cout << "Fluid properties: " << std::endl;
  std::cout << fluid << std::endl;
  
  mlProb.parameters.set<Fluid>("Fluid") = fluid;
// *************************

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  if (dim == 3) system.AddSolutionToSystemPDE("W");

  system.AddSolutionToSystemPDE("P");

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleNavierStokes_AD);

  // initilaize and solve the system
  system.init();
  system.solve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  return 0;
}


void AssembleNavierStokes_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const unsigned levelMax = mlPdeSys->GetLevelMax();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  //geometry *******************************
  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) {
    coordX[k].reserve(maxSize);
  }
  //geometry *******************************

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

  double weight; // gauss point weight

  vector< int > KKDof; // local to global pdeSys dofs
  KKDof.reserve((dim + 1) *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve((dim + 1) *maxSize);

  vector < double > Jac;
  Jac.reserve((dim + 1) *maxSize * (dim + 1) *maxSize);

  if (assembleMatrix)   KK->zero(); // Set to zero all the entries of the Global Matrix

  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->IS_Mts2Gmt_elem_offset[iproc]; iel < msh->IS_Mts2Gmt_elem_offset[iproc + 1]; iel++) {

    unsigned kel = msh->IS_Mts2Gmt_elem[iel]; // mapping between paralell dof and mesh dof
    short unsigned kelGeom = el->GetElementType(kel);    // element geometry type

    unsigned nDofsX = el->GetElementDofNumber(kel, coordXType);    // number of coordinate element dofs
    
    unsigned nDofsV = el->GetElementDofNumber(kel, solVType);    // number of solution element dofs
    unsigned nDofsP = el->GetElementDofNumber(kel, solPType);    // number of solution element dofs
    unsigned nDofsVP = dim * nDofsV + nDofsP;

    for (unsigned  k = 0; k < dim; k++) {       coordX[k].resize(nDofsX);    }
     
    for (unsigned  k = 0; k < dim; k++)   solV[k].resize(nDofsV);
    solP.resize(nDofsP);

//element matrices and vectors
    // resize local arrays
    KKDof.resize(nDofsVP);
    
    Jac.resize(nDofsVP * nDofsVP);

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].resize(nDofsV);    //resize
      std::fill(aResV[k].begin(), aResV[k].end(), 0);    //set aRes to zero
    }

    aResP.resize(nDofsP);    //resize
    std::fill(aResP.begin(), aResP.end(), 0);    //set aRes to zero

    
    // geometry ************
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, coordXType);    // local to global coordinates node
      unsigned coordXDof  = msh->GetMetisDof(iNode, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_coordinate->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }
    
    // velocity ************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, solVType);    // local to global solution node
      unsigned solVDof = msh->GetMetisDof(iNode, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        KKDof[i + k * nDofsV] = pdeSys->GetKKDof(solVIndex[k], solVPdeIndex[k], iNode);    // global to global mapping between solution node and pdeSys dof
      }
    }
    
    // pressure *************
    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, solPType);    // local to global solution node
      unsigned solPDof = msh->GetMetisDof(iNode, solPType);    // global to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      KKDof[i + dim * nDofsV] = pdeSys->GetKKDof(solPIndex, solPPdeIndex, iNode);    // global to global mapping between solution node and pdeSys dof
    }

    
    if (level == levelMax || !el->GetRefinedElementIndex(kel)) {      // do not care about this if now (it is used for the AMR)
      // start a new recording of all the operations involving adept::adouble variables
      s.new_recording();

      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solVType]->GetGaussPointNumber(); ig++) {
        // *** get gauss point weight, test function and test function partial derivatives ***
        msh->_finiteElement[kelGeom][solVType]->Jacobian(coordX, ig, weight, phiV_gss, phiV_x_gss, phiV_xx_gss);
        phiP_gss = msh->_finiteElement[kelGeom][solPType]->GetPhi(ig);

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

        double nu = 1.;
	
        // *** phiV_i loop ***
        for (unsigned i = 0; i < nDofsV; i++) {
          vector < adept::adouble > NSV_gss(dim, 0.);

          for (unsigned j = 0; j < dim; j++) {
            for (unsigned  k = 0; k < dim; k++) {
              NSV_gss[k]   +=  nu * phiV_x_gss[i * dim + j] * (gradSolV_gss[k][j] + gradSolV_gss[j][k]);  //diffusion
              NSV_gss[k]   +=  phiV_gss[i] * (solV_gss[j] * gradSolV_gss[k][j]);                                  //advection
	      NSV_gss[k]   +=  force[k] * phiV_gss[i] ;
            }
          }

          for (unsigned  k = 0; k < dim; k++) {
            NSV_gss[k] += - solP_gss * phiV_x_gss[i * dim + k];
          }

          for (unsigned  k = 0; k < dim; k++) {
            aResV[k][i] += - NSV_gss[k] * weight;
          }
        } // end phiV_i loop

        // *** phiP_i loop ***
        for (unsigned i = 0; i < nDofsP; i++) {
          for (int k = 0; k < dim; k++) {
            aResP[i] += - (gradSolV_gss[k][k]) * phiP_gss[i]  * weight;
          }
        } // end phiP_i loop

      } // end gauss point loop
      
    } // endif single element not refined or fine grid loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsVP);    //resize

    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        Res[ i +  k * nDofsV ] = -aResV[k][i].value();
      }
    }

    for (int i = 0; i < nDofsP; i++) {
      Res[ i + dim * nDofsV ] = -aResP[i].value();
    }

    RES->add_vector_blocked(Res, KKDof);

    //Extarct and store the Jacobian
    if (assembleMatrix) {
      // define the dependent variables

      for (unsigned  k = 0; k < dim; k++) {
        s.dependent(&aResV[k][0], nDofsV);
      }

      s.dependent(&aResP[0], nDofsP);

      // define the independent variables
      for (unsigned  k = 0; k < dim; k++) {
        s.independent(&solV[k][0], nDofsV);
      }

      s.independent(&solP[0], nDofsP);

      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0] , true);
      KK->add_matrix_blocked(Jac, KKDof, KKDof);

      s.clear_independents();
      s.clear_dependents();
    }
  } //end element loop for each process

  RES->close();

  if (assembleMatrix) KK->close();

  // ***************** END ASSEMBLY *******************
}




void AssembleNavierStokes(MultiLevelProblem &ml_prob){
     
  //pointers
  LinearImplicitSystem& my_lin_impl_sys    = ml_prob.get_system<LinearImplicitSystem>("NS");
  const unsigned level = my_lin_impl_sys.GetLevelToAssemble();
  const unsigned gridn = my_lin_impl_sys.GetLevelMax();
  bool assemble_matrix = my_lin_impl_sys.GetAssembleMatrix(); 
   
  Solution*	 mysolution  	             = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  mylsyspde	     = my_lin_impl_sys._LinSolver[level];   
  const char* pdename                        = my_lin_impl_sys.name().c_str();
  
  MultiLevelSolution* ml_sol=ml_prob._ml_sol;
  
  
  Mesh*		 mymsh    	= ml_prob._ml_msh->GetLevel(level);
  elem*		 myel		= mymsh->el;
  SparseMatrix*	 myKK		= mylsyspde->_KK;
  NumericVector* myRES 		= mylsyspde->_RES;
    
  //data
  const unsigned dim = mymsh->GetDimension();
  unsigned nel= mymsh->GetNumberOfElements();
  unsigned igrid= mymsh->GetLevel();
  unsigned iproc = mymsh->processor_id();
  double ILambda= 0; 
  double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  bool penalty = true; 
  const bool symm_mat = false;
  
  // solution and coordinate variables
  const char Solname[4][2] = {"U","V","W","P"};
  vector < unsigned > SolPdeIndex(dim+1);
  vector < unsigned > SolIndex(dim+1);  

  vector< vector < double> > coordinates(dim);
  
  for(unsigned ivar=0; ivar<dim; ivar++) {
    SolPdeIndex[ivar]=my_lin_impl_sys.GetSolPdeIndex(&Solname[ivar][0]);
    SolIndex[ivar]=ml_sol->GetIndex(&Solname[ivar][0]);
  }
  SolPdeIndex[dim]=my_lin_impl_sys.GetSolPdeIndex(&Solname[3][0]);
  SolIndex[dim]=ml_sol->GetIndex(&Solname[3][0]);       
  //solution order
  unsigned order_ind_vel = ml_sol->GetSolutionType(SolIndex[0]);
  unsigned order_ind_p = ml_sol->GetSolutionType(SolIndex[dim]);
  
  double alpha = 0.;
  if(order_ind_p == order_ind_vel && order_ind_vel == 0) // if pressure and velocity are both linear, we need stabilization 
  {
    alpha = 0.013333; 
  }
  
  // declare 
  vector < int > metis_node2; 
  vector < int > node1;
  vector< vector< int > > KK_dof(dim+1); 
  vector <double> phi2;
  vector <double> gradphi2;
  vector <double> nablaphi2;
  const double *phi1;
  double Weight2;
  double normal[3];
  vector< vector< double > > F(dim+1);
  vector< vector< vector< double > > > B(dim+1); 
  
  // reserve
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  metis_node2.reserve(max_size);
  node1.reserve( static_cast< unsigned > (ceil(pow(2,dim))));
  for(int i=0;i<dim;i++) {
    coordinates[i].reserve(max_size);
  }
  phi2.reserve(max_size);
  gradphi2.reserve(max_size*dim);
  nablaphi2.reserve(max_size*(3*(dim-1)));	
  for(int i=0;i<dim;i++) {
    KK_dof[i].reserve(max_size);
  }
   
  for(int i=0;i<dim+1;i++) F[i].reserve(max_size);
    
  if(assemble_matrix){
    for(int i=0;i<dim+1;i++){
      B[i].resize(dim+1);
      for(int j=0;j<dim+1;j++){
	B[i][j].reserve(max_size*max_size);
      }
    }
  }
    
  vector < double > SolVAR(dim+1);
  vector < vector < double > > gradSolVAR(dim);
  for(int i=0;i<dim;i++) {
    gradSolVAR[i].resize(dim);  
  }
  
  // Set to zeto all the entries of the matrix
  if(assemble_matrix) myKK->zero();
  
  // *** element loop ***
 
  for (int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel = mymsh->IS_Mts2Gmt_elem[iel];
    short unsigned kelt=myel->GetElementType(kel);
    unsigned nve2=myel->GetElementDofNumber(kel,order_ind_vel);
    unsigned nve1=myel->GetElementDofNumber(kel,order_ind_p);
    
    //set to zero all the entries of the FE matrices
    metis_node2.resize(nve2);
    node1.resize(nve1);
    phi2.resize(nve2);
    gradphi2.resize(nve2*dim);
    nablaphi2.resize(nve2*(3*(dim-1)));
    for(int ivar=0; ivar<dim; ivar++) {
      coordinates[ivar].resize(nve2);
      KK_dof[ivar].resize(nve2);
      
      F[SolPdeIndex[ivar]].resize(nve2);
      memset(&F[SolPdeIndex[ivar]][0],0,nve2*sizeof(double));
      
      if(assemble_matrix){
	B[SolPdeIndex[ivar]][SolPdeIndex[ivar]].resize(nve2*nve2);
	B[SolPdeIndex[ivar]][SolPdeIndex[dim]].resize(nve2*nve1);
	B[SolPdeIndex[dim]][SolPdeIndex[ivar]].resize(nve1*nve2);
	memset(&B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][0],0,nve2*nve2*sizeof(double));
	memset(&B[SolPdeIndex[ivar]][SolPdeIndex[dim]][0],0,nve2*nve1*sizeof(double));
	memset(&B[SolPdeIndex[dim]][SolPdeIndex[ivar]][0],0,nve1*nve2*sizeof(double));
      }
    }
    KK_dof[dim].resize(nve1);
    F[SolPdeIndex[dim]].resize(nve1);
    memset(&F[SolPdeIndex[dim]][0],0,nve1*sizeof(double));
      
      
    if(assemble_matrix*penalty){
      B[SolPdeIndex[dim]][SolPdeIndex[dim]].resize(nve1*nve1,0.);
      memset(&B[SolPdeIndex[dim]][SolPdeIndex[dim]][0],0,nve1*nve1*sizeof(double));
    }
    
    for( unsigned i=0;i<nve2;i++){
      unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;
      unsigned inode_coord_metis=mymsh->GetMetisDof(inode,2);
      metis_node2[i]=mymsh->GetMetisDof(inode,order_ind_vel);
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*mymsh->_coordinate->_Sol[ivar])(inode_coord_metis);
	KK_dof[ivar][i]=mylsyspde->GetKKDof(SolIndex[ivar],SolPdeIndex[ivar],inode);
      }
    }
    
    double hk = sqrt( (coordinates[0][2] - coordinates[0][0])*(coordinates[0][2] - coordinates[0][0]) + 
      (coordinates[1][2] - coordinates[1][0])*(coordinates[1][2] - coordinates[1][0]) );
    
    for(unsigned i=0;i<nve1;i++) {
      unsigned inode=(order_ind_p<dim)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
      node1[i]=inode;
      KK_dof[dim][i]=mylsyspde->GetKKDof(SolIndex[dim],SolPdeIndex[dim],inode);
    }
   
    if(igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind_vel]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	ml_prob._ml_msh->_finiteElement[kelt][order_ind_vel]->Jacobian(coordinates,ig,Weight2,phi2,gradphi2,nablaphi2);
	phi1=ml_prob._ml_msh->_finiteElement[kelt][order_ind_p]->GetPhi(ig);

	double GradSolP[3] = {0.,0.,0.};
	//velocity variable
	for(unsigned ivar=0; ivar<dim; ivar++) {
	  SolVAR[ivar]=0;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR[ivar][ivar2]=0; 
	  }
	  unsigned SolIndex=ml_sol->GetIndex(&Solname[ivar][0]);
	  unsigned SolType=ml_sol->GetSolutionType(&Solname[ivar][0]);
	  for(unsigned i=0; i<nve2; i++) {
	    double soli = (*mysolution->_Sol[SolIndex])(metis_node2[i]);
	    SolVAR[ivar]+=phi2[i]*soli;
	    for(unsigned ivar2=0; ivar2<dim; ivar2++){
	      gradSolVAR[ivar][ivar2] += gradphi2[i*dim+ivar2]*soli; 
	    }
	  }
	}
	//pressure variable
	SolVAR[dim]=0;
	unsigned SolIndex=ml_sol->GetIndex(&Solname[3][0]);
	unsigned SolType=ml_sol->GetSolutionType(&Solname[3][0]);
	for(unsigned i=0; i<nve1; i++){
	  unsigned sol_dof = mymsh->GetMetisDof(node1[i],SolType);
	  double soli = (*mysolution->_Sol[SolIndex])(sol_dof);
	  SolVAR[dim]+=phi1[i]*soli;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){
	    GradSolP[ivar2] += gradphi2[i*dim+ivar2]*soli;
	  }
	}

	
	// *** phi_i loop ***
	for(unsigned i=0; i<nve2; i++){
	
	  //BEGIN RESIDUALS A block ===========================
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    double Lap_rhs=0;
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      Lap_rhs += gradphi2[i*dim+ivar2]*gradSolVAR[ivar][ivar2];
	    }
	    F[SolPdeIndex[ivar]][i]+= (-IRe*Lap_rhs + SolVAR[dim]*gradphi2[i*dim+ivar] + force[ivar] * phi2[i])*Weight2;
	  }
	  //END RESIDUALS A block ===========================
	  
	  if(assemble_matrix){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      double Lap=0;
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		// Laplacian
		Lap  += gradphi2[i*dim+ivar]*gradphi2[j*dim+ivar]*Weight2;
	      }

	      for(unsigned ivar=0; ivar<dim; ivar++) {    
		B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][i*nve2+j] += IRe*Lap;
	      }
  	    } //end phij loop
	    
	    // *** phi1_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[ivar]][SolPdeIndex[dim]][i*nve1+j] -= gradphi2[i*dim+ivar]*phi1[j]*Weight2;
	      }
	    } //end phi1_j loop
	  } // endif assemble_matrix
	} //end phii loop
  

	// *** phi1_i loop ***
	for(unsigned i=0; i<nve1; i++){
	  //BEGIN RESIDUALS B block ===========================
	  double div = 0;
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    div += gradSolVAR[ivar][ivar];
	  }
	  F[SolPdeIndex[dim]][i]+= (phi1[i]*div + /*penalty*ILambda*phi1[i]*SolVAR[dim]*/ 
	                             + 0.*((hk*hk)/(4.*IRe))*alpha*(GradSolP[0]*gradphi2[i*dim + 0] + GradSolP[1]*gradphi2[i*dim + 1]) )*Weight2;
				     
           
	  //END RESIDUALS  B block ===========================
	  
	  if(assemble_matrix){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[dim]][SolPdeIndex[ivar]][i*nve2+j]-= phi1[i]*gradphi2[j*dim+ivar]*Weight2;
	      }
	    }  //end phij loop
	  } // endif assemble_matrix
	}  //end phi1_i loop
	
	if(assemble_matrix * penalty){  //block nve1 nve1
	  // *** phi_i loop ***
	  for(unsigned i=0; i<nve1; i++){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      //B[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nve1+j]-= ILambda*phi1[i]*phi1[j]*Weight2;
	      for(unsigned ivar=0; ivar<dim; ivar++) {
	        B[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nve1+j] -= ((hk*hk)/(4.*IRe))*alpha*(gradphi2[i*dim + ivar]*gradphi2[j*dim + ivar])*Weight2; 
	      }
	    }
	  }
	}   //end if penalty
      }  // end gauss point loop
      
      //--------------------------------------------------------------------------------------------------------
      // Boundary Integral --> to be added
      //number of faces for each type of element
//       if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
//      
// 	unsigned nfaces = myel->GetElementFaceNumber(kel);
// 
// 	// loop on faces
// 	for(unsigned jface=0;jface<nfaces;jface++){ 
// 	  
// 	  // look for boundary faces
// 	  if(myel->GetFaceElementIndex(kel,jface)<0){
// 	    for(unsigned ivar=0; ivar<dim; ivar++) {
// 	      ml_prob.ComputeBdIntegral(pdename, &Solname[ivar][0], kel, jface, level, ivar);
// 	    }
// 	  }
// 	}	
//       }
      //--------------------------------------------------------------------------------------------------------
    } // endif single element not refined or fine grid loop
    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned ivar=0; ivar<dim; ivar++) {
      myRES->add_vector_blocked(F[SolPdeIndex[ivar]],KK_dof[ivar]);
      if(assemble_matrix){
	myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[ivar]],KK_dof[ivar],KK_dof[ivar]);  
	myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[dim]],KK_dof[ivar],KK_dof[dim]);
	myKK->add_matrix_blocked(B[SolPdeIndex[dim]][SolPdeIndex[ivar]],KK_dof[dim],KK_dof[ivar]);
      }
    }
    //Penalty
    if(assemble_matrix*penalty) myKK->add_matrix_blocked(B[SolPdeIndex[dim]][SolPdeIndex[dim]],KK_dof[dim],KK_dof[dim]);
    myRES->add_vector_blocked(F[SolPdeIndex[dim]],KK_dof[dim]);
    //--------------------------------------------------------------------------------------------------------  
  } //end list of elements loop for each subdomain
  
  
  if(assemble_matrix) myKK->close();
  myRES->close();
  // ***************** END ASSEMBLY *******************
}
