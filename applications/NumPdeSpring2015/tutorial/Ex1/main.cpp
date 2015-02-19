 #include "MultiLevelProblem.hpp"

#include "NumericVector.hpp"
// #include "Fluid.hpp"
// #include "Parameter.hpp"

// #include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "LinearImplicitSystem.hpp"
// #include "FElemTypeEnum.hpp"
// #include "Files.hpp"

// using std::cout;
// using std::endl;


#include "FemusInit.hpp"
using namespace femus;

double InitalValueU(const double &x, const double &y, const double &z){
  return x+y;
}

double InitalValueP(const double &x, const double &y, const double &z){
  return x;
}

double InitalValueT(const double &x, const double &y, const double &z){
  return y;
}

bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], double &value, const int facename, const double time) {
  bool dirichlet=1; //dirichlet
  value=0.;
  return dirichlet;
}

void AssemplePoissonProblem(MultiLevelProblem &ml_prob, unsigned level, const unsigned &levelMax, const bool &assembleMatrix);


int main(int argc, char **args) {
  
  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
 
  // read coarse level mesh and generate finers level meshes
  MultiLevelMesh mlMsh;
  double scalingFactor=1.; 
  mlMsh.ReadCoarseMesh("./input/square.neu","seventh",scalingFactor); 
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfLevels=1;
  unsigned numberOfSelectiveLevels=0;
  mlMsh.RefineMesh(numberOfLevels , numberOfLevels + numberOfSelectiveLevels, NULL);
  
  mlMsh.EraseCoarseLevels(numberOfLevels-1);
  
  mlMsh.PrintInfo();
  
  // define and initialize variables
  MultiLevelSolution mlSol(&mlMsh);
  
  mlSol.AddSolution("u",LAGRANGE, FIRST);
  
  //   mlSol.AddSolution("V",LAGRANGE, SERENDIPITY);
  //   mlSol.AddSolution("W",LAGRANGE, SECOND);
  //   mlSol.AddSolution("P",DISCONTINOUS_POLYNOMIAL, ZERO);
  //   mlSol.AddSolution("T",DISCONTINOUS_POLYNOMIAL, FIRST);
  
  //  mlSol.Initialize("All"); 
  
  mlSol.Initialize("u", InitalValueU);
  //   mlSol.Initialize("P", InitalValueP);
  //   mlSol.Initialize("T", InitalValueT); 
  // //note that this initialization is the same as piecewise constant element
  
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("u");
  
  MultiLevelProblem ml_prob(&mlMsh,&mlSol);  
  LinearImplicitSystem & system = ml_prob.add_system<LinearImplicitSystem>("Poisson");
  system.AddSolutionToSystemPDE("u");
   
  // Set MG Options
  system.SetAssembleFunction(AssemplePoissonProblem);  
  system.init();
  system.solve();
     
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("u");
//   variablesToBePrinted.push_back("P");
//   variablesToBePrinted.push_back("T");

  VTKWriter vtkIO(mlSol);
  vtkIO.write_system_solutions(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  GMVWriter gmvIO(mlSol);
  variablesToBePrinted.push_back("all");
  gmvIO.SetDebugOutput(false);
  gmvIO.write_system_solutions(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);
   
  return 0;
}

void AssemplePoissonProblem(MultiLevelProblem &ml_prob, unsigned level, const unsigned &levelMax, const bool &assembleMatrix) {
  // ml_prob is the global object from/to where get/set all the data 
  // level is the level of the PDE system to be assembled  
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled
    
  //pointers and references
  Mesh*         	msh	       	= ml_prob._ml_msh->GetLevel(level); // pointer to the mesh (level) object 
  elem*         	el	       	= msh->el;  // pointer to the elem object in msh (level) 				
  
  MultiLevelSolution* 	mlSol       	= ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* 		sol       	= ml_prob._ml_sol->GetSolutionLevel(level); // pointer to the solution (level) object
  
  LinearImplicitSystem& mlPdeSys 	= ml_prob.get_system<LinearImplicitSystem>("Poisson"); // reference to the linear implicit system named "Poisson" 
  LinearEquationSolver* pdeSys      	= mlPdeSys._LinSolver[level]; // pointer to the equation (level) object 
  SparseMatrix*  	KK	       	= pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* 	RES	       	= pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)
  
  const unsigned	dim	= msh->GetDimension(); // get the domain dimension of the problem
  unsigned 		iproc	= msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned solIndex;
  solIndex = mlSol->GetIndex("u"); // get the position of "u" in the ml_sol object
  unsigned solType = mlSol->GetSolutionType(solIndex); // get the finite element type for "u"  
    
  unsigned solPdeIndex;
  solPdeIndex = mlPdeSys.GetSolPdeIndex("u"); // get the position of "u" in the pdeSys object
    
  vector < double >  soli; // local solution
  
  vector < vector < double > > x(dim); // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > KKDof; // local to global solPde dofs
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight
  vector< double > Res; // local redidual vector
  vector< double > K; // local stiffness matrix
    
  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned > (ceil(pow(3,dim))); // conservative: based on line3, quad9, hex27
 
  
  soli.reserve(maxSize);
  for(unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);
  KKDof.reserve(maxSize);
  phi.reserve(maxSize);
  phi_x.reserve(maxSize*dim);
  unsigned dim2=(3*(dim-1)+!(dim-1));
  phi_xx.reserve(maxSize*dim2);
  Res.reserve(maxSize);
  K.reserve(maxSize*maxSize);
 
  double nu = 1.; //diffusive constant
  
  if( assembleMatrix ) 
    KK->zero(); // Set to zero all the entries of the Global Matrix
  
  // element loop: each process loops only on the elements that owns 
  for (int iel=msh->IS_Mts2Gmt_elem_offset[iproc]; iel < msh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {
    
    unsigned kel = msh->IS_Mts2Gmt_elem[iel]; // mapping between paralell dof and mesh dof 
    short unsigned kelGeom = el->GetElementType( kel ); // element geometry type
    unsigned nDofs  = el->GetElementDofNumber( kel, solType); // number of solution element dofs
    unsigned nDofs2 = el->GetElementDofNumber( kel, xType); // number of coordinate element dofs
    
    // resize local arrays
    KKDof.resize(nDofs);
    soli.resize(nDofs);
    for(int i=0; i<dim; i++) {
      x[i].resize(nDofs2);
    }
    Res.resize(nDofs); //resize
    memset(&Res[0], 0 , nDofs*sizeof(double) ); //set Res to zero
    if(assembleMatrix) {
      K.resize(nDofs*nDofs); //resize
      memset(&K[0], 0, nDofs*nDofs*sizeof(double) ); //set K to zero
    }
         
    // local storage of global mapping and solution
    for( unsigned i=0; i<nDofs; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, solType);  // local to global solution node   
      unsigned solDof = msh->GetMetisDof(iNode, solType); // global to global mapping between solution node and solution dof
      soli[i] = (*sol->_Sol[solIndex])(solDof); // global extraction and local storage for the solution
      KKDof[i] = pdeSys->GetKKDof(solIndex, solPdeIndex, iNode); // global to global mapping between solution node and pdeSys dof
    }
        
     // local storage of coordinates 
    for( unsigned i=0; i<nDofs2; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, xType);  // local to global coordinates node     
      unsigned xDof  = msh->GetMetisDof(iNode, xType); // global to global mapping between coordinates node and coordinate dof
      for(unsigned jdim=0; jdim<dim; jdim++) {
	x[jdim][i] = (*msh->_coordinate->_Sol[jdim])(xDof); // global extraction and local storage for the element coordinates
      }
    }    
        
    if( level == levelMax || !el->GetRefinedElementIndex(kel)) { // do not care now of this if now
      // *** Gauss point loop ***
      for(unsigned ig=0; ig < msh->_finiteElement[kelGeom][solType]->GetGaussPointNumber(); ig++) {
	// *** get gauss weight, test function and test function derivatives ***
	msh->_finiteElement[kelGeom][solType]->Jacobian(x,ig,weight,phi,phi_x,phi_xx);
	
	// evaluate the gauss point solution, solution derivatives and coordinates
	double solG = 0;
	vector < double > solG_x(dim,0.);
	vector < double > xG(dim,0.);
	
	for(unsigned i=0; i<nDofs; i++) {
	  solG+=phi[i]*soli[i];
	  for(unsigned jdim=0; jdim<dim; jdim++) {
	    solG_x[jdim] += phi_x[i*dim+jdim]*soli[i];
	    xG[jdim] += x[jdim][i]*phi[i]; 
	  }
	}
        
        // *** phi_i loop ***
	for(unsigned i=0; i<nDofs; i++) {
	 
	  double laplace = 0.;
	  for(unsigned jdim=0; jdim<dim; jdim++) {
	    laplace   +=  phi_x[i*dim+jdim]*solG_x[jdim];
	  }
	  double srcTerm = 1.;  
	  Res[i]+= (srcTerm*phi[i] - nu*laplace) * weight;
	
	  if( assembleMatrix ) {
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nDofs; j++) {
	      laplace = 0.;
	      for(unsigned kdim=0; kdim<dim; kdim++) {
		laplace += ( phi_x[i*dim + kdim] * phi_x[j*dim + kdim] )*weight;
	      }
	      K[i*nDofs+j] += nu * laplace;
	    } // end phi_j loop
	  } // endif assemble_matrix
	} // end phi_i loop
      } // end gauss point loop
    } // endif single element not refined or fine grid loop

    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the global Matrix/vector

    RES->add_vector_blocked(Res,KKDof);
    if(assembleMatrix) KK->add_matrix_blocked(K,KKDof,KKDof);
  } //end element loop for each process

  RES->close();
  if( assembleMatrix ) KK->close();
  // ***************** END ASSEMBLY *******************

  
}


