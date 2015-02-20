/** tutorial/Ex2
 * This example shows how to set and solve the weak form of the Poisson problem 
 *                    $$ \Delta u = 1 \text{ on }\Omega, $$
 * 		      $$ u=0 \text{ on } \Gamma, $$
 * on a square domain $\Omega$ with boundary $\Gamma$;
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "LinearImplicitSystem.hpp"

using namespace femus;

bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], double &value, const int facename, const double time) {
  bool dirichlet=1; //dirichlet
  value=0.;
  return dirichlet;
}

void AssemblePoissonProblem(MultiLevelProblem &ml_prob, unsigned level, const unsigned &levelMax, const bool &assembleMatrix);


int main(int argc, char **args) {
  
  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
 
  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor=1.; 
  mlMsh.ReadCoarseMesh("./input/square.neu","seventh",scalingFactor); 
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels=2;
  unsigned numberOfSelectiveLevels=0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  
  // erase all the coarse mesh levels 
  mlMsh.EraseCoarseLevels(numberOfUniformLevels-1);
  
  // print mesh info
  mlMsh.PrintInfo();
  
  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);
  
  // add variables to mlSol
  mlSol.AddSolution("u",LAGRANGE, FIRST);
  mlSol.Initialize("All");
  
  // attach the boundary condition function and generate boundary data 
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("u");
  
  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);  
  
  // add system Poisson in mlProb as a Linear Implicit System
  LinearImplicitSystem & system = mlProb.add_system < LinearImplicitSystem > ("Poisson");
  
  // add solution "u" to system
  system.AddSolutionToSystemPDE("u");
   
  // attach the assembling function to system
  system.SetAssembleFunction(AssemblePoissonProblem);  
  
  // initilaize and solve the system 
  system.init();
  system.solve();
     
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(mlSol);
  vtkIO.write_system_solutions(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  GMVWriter gmvIO(mlSol);
  gmvIO.SetDebugOutput(true);
  gmvIO.write_system_solutions(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);
   
  return 0;
}


/**
 * This function assemble the stiffnes matrix KK and the residual vector Res
 * such that        
 *                  KK w = RES = F - KK u0, 
 * and consequently 
 * 		    u = u0 + w satisfies KK u = F
 **/

void AssemblePoissonProblem(MultiLevelProblem &ml_prob, unsigned level, const unsigned &levelMax, const bool &assembleMatrix) {
  //  ml_prob is the global object from/to where get/set all the data 
  //  level is the level of the PDE system to be assembled  
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled
      
  //  extract pointers to the several objects that we are going to use 
  Mesh*         	msh	       	= ml_prob._ml_msh->GetLevel(level); // pointer to the mesh (level) object 
  elem*         	el	       	= msh->el;  // pointer to the elem object in msh (level) 				
  
  MultiLevelSolution* 	mlSol       	= ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* 		sol       	= ml_prob._ml_sol->GetSolutionLevel(level); // pointer to the solution (level) object
  
  LinearImplicitSystem* mlPdeSys 	= &ml_prob.get_system<LinearImplicitSystem>("Poisson"); // pointer to the linear implicit system named "Poisson" 
  LinearEquationSolver* pdeSys      	= mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object 
  SparseMatrix*  	KK	       	= pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* 	RES	       	= pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)
  
  const unsigned	dim	= msh->GetDimension(); // get the domain dimension of the problem
  unsigned 		iproc	= msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex;
  soluIndex = mlSol->GetIndex("u"); // get the position of "u" in the ml_sol object
  unsigned soluType = mlSol->GetSolutionType(soluIndex); // get the finite element type for "u"  
    
  unsigned soluPdeIndex;
  soluPdeIndex = mlPdeSys->GetSolPdeIndex("u"); // get the position of "u" in the pdeSys object
    
  vector < double >  solu; // local solution
  
  vector < vector < double > > x(dim); // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  vector< int > KKDof; // local to global pdeSys dofs
  vector <double> phi;  // local test function
  vector <double> phi_x; // local test function first order partial derivatives
  vector <double> phi_xx; // local test function second order partial derivatives
  double weight; // gauss point weight
  
  vector< double > Res; // local redidual vector
  vector< double > K; // local stiffness matrix
    
  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned > (ceil(pow(3,dim))); // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);
  for(unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);
  KKDof.reserve(maxSize);
  phi.reserve(maxSize);
  phi_x.reserve(maxSize*dim);
  unsigned dim2=(3*(dim-1)+!(dim-1)); // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  phi_xx.reserve(maxSize*dim2);
  Res.reserve(maxSize);
  K.reserve(maxSize*maxSize);
   
  if( assembleMatrix ) 
    KK->zero(); // Set to zero all the entries of the Global Matrix
  
  // element loop: each process loops only on the elements that owns 
  for (int iel=msh->IS_Mts2Gmt_elem_offset[iproc]; iel < msh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {
    
    unsigned kel = msh->IS_Mts2Gmt_elem[iel]; // mapping between paralell dof and mesh dof 
    short unsigned kelGeom = el->GetElementType( kel ); // element geometry type
    unsigned nDofs  = el->GetElementDofNumber( kel, soluType); // number of solution element dofs
    unsigned nDofs2 = el->GetElementDofNumber( kel, xType); // number of coordinate element dofs
    
    // resize local arrays
    KKDof.resize(nDofs);
    solu.resize(nDofs);
    for(int i=0; i<dim; i++) {
      x[i].resize(nDofs2);
    }
    Res.resize(nDofs); //resize
    std::fill(Res.begin(), Res.end(), 0); //set Res to zero
    if(assembleMatrix) {
      K.resize(nDofs*nDofs); //resize
      std::fill(K.begin(), K.end(), 0); //set K to zero
    }
         
    // local storage of global mapping and solution
    for( unsigned i=0; i<nDofs; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, soluType);  // local to global solution node   
      unsigned solDof = msh->GetMetisDof(iNode, soluType); // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof); // global extraction and local storage for the solution
      KKDof[i] = pdeSys->GetKKDof(soluIndex, soluPdeIndex, iNode); // global to global mapping between solution node and pdeSys dof
    }
        
     // local storage of coordinates 
    for( unsigned i=0; i<nDofs2; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, xType);  // local to global coordinates node     
      unsigned xDof  = msh->GetMetisDof(iNode, xType); // global to global mapping between coordinates node and coordinate dof
      for(unsigned jdim=0; jdim<dim; jdim++) {
	x[jdim][i] = (*msh->_coordinate->_Sol[jdim])(xDof); // global extraction and local storage for the element coordinates
      }
    }    
        
    if( level == levelMax || !el->GetRefinedElementIndex(kel)) { // do not care about this if now (it is used for the AMR)
      // *** Gauss point loop ***
      for(unsigned ig=0; ig < msh->_finiteElement[kelGeom][soluType]->GetGaussPointNumber(); ig++) {
	// *** get gauss point weight, test function and test function partial derivatives ***
	msh->_finiteElement[kelGeom][soluType]->Jacobian(x,ig,weight,phi,phi_x,phi_xx);
	
	// evaluate the solution, the solution derivatives and the coordinates in the gauss point
	double soluGauss = 0;
	vector < double > soluGauss_x(dim,0.);
	vector < double > xGauss(dim,0.);
	
	for(unsigned i=0; i<nDofs; i++) {
	  soluGauss+=phi[i]*solu[i];
	  for(unsigned jdim=0; jdim<dim; jdim++) {
	    soluGauss_x[jdim] += phi_x[i*dim+jdim]*solu[i];
	    xGauss[jdim] += x[jdim][i]*phi[i]; 
	  }
	}
        
        // *** phi_i loop ***
	for(unsigned i=0; i<nDofs; i++) {
	 
	  double laplace = 0.;
	  for(unsigned jdim=0; jdim<dim; jdim++) {
	    laplace   +=  phi_x[i*dim+jdim]*soluGauss_x[jdim];
	  }
	  double srcTerm = 1.;  
	  Res[i]+= (srcTerm*phi[i] - laplace) * weight;
	
	  if( assembleMatrix ) {
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nDofs; j++) {
	      laplace = 0.;
	      for(unsigned kdim=0; kdim<dim; kdim++) {
		laplace += ( phi_x[i*dim + kdim] * phi_x[j*dim + kdim] )*weight;
	      }
	      K[i*nDofs+j] += laplace;
	    } // end phi_j loop
	  } // endif assemble_matrix
	} // end phi_i loop
      } // end gauss point loop
    } // endif single element not refined or fine grid loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    RES->add_vector_blocked(Res,KKDof);
    if(assembleMatrix) KK->add_matrix_blocked(K,KKDof,KKDof);
  } //end element loop for each process

  RES->close();
  if( assembleMatrix ) KK->close();
  // ***************** END ASSEMBLY *******************

  
}


