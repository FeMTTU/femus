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

bool AssemplePoissonProblem(const double &x, const double &y, const double &z,const char name[], double &value, const int facename, const double time) {
  bool dirichlet=1; //dirichlet
  value=0.;
  return dirichlet;
}

void SetPoissonProblem(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assemble_matrix);

int main(int argc, char **args) {
  
  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
 
  // read coarse level mesh and generate finers level meshes
  MultiLevelMesh mlMsh;
  double scalingFactor=1.; 
  mlMsh.ReadCoarseMesh("./input/square.neu","seventh",scalingFactor); 
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfLevels=4;
  unsigned numberOfSelectiveLevels=0;
  mlMsh.RefineMesh(numberOfLevels , numberOfLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();
  
  // define and initialize variables
  MultiLevelSolution mlSol(&mlMsh);
  
  mlSol.AddSolution("U",LAGRANGE, FIRST);
  mlSol.AddSolution("V",LAGRANGE, SERENDIPITY);
  mlSol.AddSolution("W",LAGRANGE, SECOND);
  mlSol.AddSolution("P",DISCONTINOUS_POLYNOMIAL, ZERO);
  mlSol.AddSolution("T",DISCONTINOUS_POLYNOMIAL, FIRST);
  
  mlSol.Initialize("All"); 
  
  mlSol.Initialize("U", InitalValueU);
  mlSol.Initialize("P", InitalValueP);
  mlSol.Initialize("T", InitalValueT); 
  //note that this initialization is the same as piecewise constant element
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("U");
  variablesToBePrinted.push_back("P");
  variablesToBePrinted.push_back("T");

  VTKWriter vtkIO(mlSol);
  vtkIO.write_system_solutions(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  GMVWriter gmvIO(mlSol);
  variablesToBePrinted.push_back("all");
  gmvIO.SetDebugOutput(false);
  gmvIO.write_system_solutions(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);
   
  return 1;
}

void AssemplePoissonProblem(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assemble_matrix) {

  //pointers and references
  Mesh*         	msh	       	= ml_prob._ml_msh->GetLevel(level); // pointer to the mesh (level) object 
  elem*         	el	       	= msh->el;  // pointer to the elem object in msh (level) 				
  
  MultiLevelSolution* 	ml_sol       	= ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* 		solution       	= ml_prob._ml_sol->GetSolutionLevel(level); // pointer to the solution (level) object
  
  LinearImplicitSystem& ml_linsyspde 	= ml_prob.get_system<LinearImplicitSystem>("Poisson"); // reference to the linear implicit system named "Poisson" 
  LinearEquationSolver* lsyspde      	= ml_linsyspde._LinSolver[level]; // pointer to the equation (level) object 
  SparseMatrix*  	KK	       	= lsyspde->_KK;  // pointer to the stifness matrix object in equation (level)
  NumericVector* 	RES	       	= lsyspde->_RES; // pointer to the residual vector object in equation (level)
  
  const unsigned	dim	= msh->GetDimension(); // get the dimension of the mesh
  unsigned 		iproc	= msh->processor_id(); // get the process_id for parallel computation 

  //solution variable
  unsigned SolIndex;
  SolIndex = ml_sol->GetIndex("u"); // get the position of "u" in the ml_sol object
  unsigned solType = ml_sol->GetSolutionType(SolIndex); // get the finite element type of "u"  
    
  unsigned SolPdeIndex;
  SolPdeIndex = ml_linsyspde.GetSolPdeIndex("u"); // get the position of "u" in the lin_impl_sys object
       
  //coordinates
  vector< vector < double> > coordinates(dim);

  // declare
  vector< int > metis_node;
  vector< int > KK_dof;
  vector <double> phi;
  vector <double> gradphi;
  vector <double> nablaphi;
  double weight;
  vector< double > F;
  vector< double > B;
  vector<double> normal(3.0);
  double src_term = 0.;
  vector<double> xyzt(4,0.);
//   ParsedFunction* bdcfunc = NULL;

  // reserve
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  metis_node.reserve(max_size);
  KK_dof.reserve(max_size);
  for(int i=0; i<dim; i++)
    coordinates[i].reserve(max_size);
  phi.reserve(max_size);
  gradphi.reserve(max_size*dim);
  unsigned nabla_dim=(3*(dim-1)+!(dim-1));
  nablaphi.reserve(max_size*nabla_dim);
  F.reserve(max_size);
  B.reserve(max_size*max_size);
 
    
  // Set to zeto all the entries of the Global Matrix
  if(assemble_matrix) 
    KK->zero();

  
  // *** element loop ***
  for (int iel=msh->IS_Mts2Gmt_elem_offset[iproc]; iel < msh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel = msh->IS_Mts2Gmt_elem[iel];
    short unsigned kelt=el->GetElementType(kel);
    unsigned nve=el->GetElementDofNumber( kel, solType);
    unsigned nve2=el->GetElementDofNumber(kel,2);
    // resize
    metis_node.resize(nve);
    KK_dof.resize(nve);
    
    for(int i=0; i<dim; i++) {
      coordinates[i].resize(nve);
    }

    // set to zero all the entries of the FE matrices
    F.resize(nve);
    memset(&F[0],0,nve*sizeof(double));
    if(assemble_matrix) {
      B.resize(nve*nve);
      memset(&B[0],0,nve*nve*sizeof(double));
    }

    // get local to global mappings
    for( unsigned i=0; i<nve2; i++) {
      unsigned inode=el->GetMeshDof(kel,i, solType);
      unsigned inode_coord_metis=msh->GetMetisDof(inode,2);
     
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*msh->_coordinate->_Sol[ivar])(inode_coord_metis);
      }
      if(i<nve){
	metis_node[i]=msh->GetMetisDof(inode, solType);
	KK_dof[i]=lsyspde->GetKKDof(SolIndex,SolPdeIndex,inode);
      }
    }
    
   
      
    if(level==gridn || !el->GetRefinedElementIndex(kel)) {
      // *** Gauss point loop ***
                 
      // Supg stabilization tau evaluation	
      double V[3]={0.,0.,0.}; 
      double nu=1.;
      if(dim==1){
	V[0]=1.; 
	nu=0.01;
      }
      else if(dim==2){
// 	nu=0.0001;
// 	V[0]=sqrt(2)/2;
// 	V[1]=-sqrt(2)/2;
	nu=1.;
 	V[0]=0.;
 	V[1]=0.;
      }
      else if(dim==3){
	nu=1.;
	V[0]=V[1]=V[2]=0.;
      }
      
      double barNu=0.;
      double vL2Norm2=0.;
      for(int i=0;i<dim;i++){
	vL2Norm2 += V[i]*V[i];
	unsigned ip = referenceElementDirection[kelt][i][1];
	unsigned im = referenceElementDirection[kelt][i][0];
	double VxiHxi=0.;
	for(int j=0;j<dim;j++){
	  VxiHxi += (coordinates[j][ip]-coordinates[j][im]) * V[j];
	}	
	double PeXi=VxiHxi/(2.*nu);		
	double barXi = ( fabs( PeXi ) < 1.0e-10) ? 0. : 1./tanh(PeXi)-1./PeXi;
	barNu += barXi * VxiHxi /2.;
      }
      double supgTau = ( vL2Norm2 > 1.0e-15 ) ? barNu/vL2Norm2 : 0.;
      // End Stabilization stabilization tau evaluation
            
      for(unsigned ig=0; ig < msh->_finiteElement[kelt][solType]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	msh->_finiteElement[kelt][solType]->Jacobian(coordinates,ig,weight,phi,gradphi,nablaphi);
	//current solution
	double SolT=0;
	vector < double > gradSolT(dim,0.);
	vector < double > NablaSolT(dim,0.);
	
	xyzt.assign(4,0.);
	unsigned SolType=ml_sol->GetSolutionType("u");
	for(unsigned i=0; i<nve; i++) {
	  double soli = (*solution->_Sol[SolIndex])(metis_node[i]);
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    xyzt[ivar] += coordinates[ivar][i]*phi[i]; 
	  }
	  SolT+=phi[i]*soli;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	    gradSolT[ivar2] += gradphi[i*dim+ivar2]*soli;
	    NablaSolT[ivar2] += nablaphi[i*nabla_dim+ivar2]*soli;
	  }
	}
        // *** phi_i loop ***
	for(unsigned i=0; i<nve; i++) {
	  //BEGIN RESIDUALS A block ===========================
	  double advRhs=0.;
	  double lapRhs=0.;
	  double resRhs=0.;
	  double supgPhi=0.;  
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    lapRhs   +=  nu*gradphi[i*dim+ivar]*gradSolT[ivar];
	    advRhs   +=  V[ivar]*gradSolT[ivar]*phi[i];
	    resRhs   += -nu*NablaSolT[ivar] + V[ivar]*gradSolT[ivar];
	    supgPhi  += (V[ivar] * gradphi[i*dim+ivar] + nu*nablaphi[i*nabla_dim + ivar] )* supgTau;
	  }
	  
	  src_term = 1.;
            	  
	  F[i]+= (  src_term* phi[i] -lapRhs - advRhs   
		  +(src_term - resRhs) * supgPhi )*weight;
		    
	  //END RESIDUALS A block ===========================
	  if(assemble_matrix) {
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve; j++) {
	      double lap=0;
	      double adv=0;
	            
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		lap += nu*( gradphi[i*dim +ivar] * gradphi[j*dim+ivar] 
			   -nablaphi[j*nabla_dim + ivar] * supgPhi )*weight;
		adv += V[ivar]*gradphi[j*dim+ivar]* (phi[i] + supgPhi)*weight;
	      }
	      B[i*nve+j] += lap + adv ;
	    } // end phij loop
	  } // end phii loop
	} // endif assemble_matrix
      } // end gauss point loop
    } // endif single element not refined or fine grid loop

    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the global Matrix/vector

    RES->add_vector_blocked(F,KK_dof);
    if(assemble_matrix) KK->add_matrix_blocked(B,KK_dof,KK_dof);
  } //end list of elements loop for each subdomain

  RES->close();
  if(assemble_matrix) KK->close();

  // ***************** END ASSEMBLY *******************

}


