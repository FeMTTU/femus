
#include <sstream>
#include "FemusDefault.hpp"
#include "FemusInit.hpp"
#include "MultiLevelMesh.hpp"
#include "MultiLevelProblem.hpp"
#include "LinearImplicitSystem.hpp"
#include "NumericVector.hpp"
#include "WriterEnum.hpp"
#include "VTKWriter.hpp"



using namespace femus;

// So far, we need to be careful 
// - how we call the groups in Salome (Group_N_M syntax, with N>1)
// - how we export MED by selecting all groups and meshes
// Moreover, we can see the groups in ParaVIS, but only not simultaneously

void AssembleFracFSI(MultiLevelProblem& ml_prob);


double InitialValue_u(const std::vector < double >& x) {
  return 0.;
}


bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

  bool dirichlet = true; //dirichlet
  value = 0;
  
//   if(!strcmp(name,"u")) {
//   if (faceName == 3)
//     dirichlet = false;
//   }
  
  return dirichlet;
}


int main(int argc,char **args) {

  FemusInit init(argc,args,MPI_COMM_WORLD);
  
  std::string med_file = "RectFracWithGroup.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << med_file;
  const std::string infile = mystream.str();
 
  //Adimensional
  double Lref = 1.;
  
  MultiLevelMesh ml_msh;
  ml_msh.ReadCoarseMesh(infile.c_str(),"fifth",Lref);
   
  
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  ml_msh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  ml_msh.PrintInfo();

  ml_msh.SetWriter(XDMF);
  ml_msh.GetWriter()->write(DEFAULT_OUTPUTDIR,"biquadratic");
  ml_msh.SetWriter(GMV);
  ml_msh.GetWriter()->write(DEFAULT_OUTPUTDIR,"biquadratic");
  ml_msh.SetWriter(VTK);
  ml_msh.GetWriter()->write(DEFAULT_OUTPUTDIR,"biquadratic");
  
  
  // define the multilevel solution and attach the ml_msh object to it
  MultiLevelSolution mlSol(&ml_msh);

  // add variables to mlSol
  mlSol.AddSolution("u", LAGRANGE, SECOND);
 
  mlSol.Initialize("All");    // initialize all varaibles to zero

  mlSol.Initialize("u", InitialValue_u);
 
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition); 
  mlSol.GenerateBdc("u"); 
  
 // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  
 // add system  in mlProb as a Linear Implicit System
  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("FracFSI");
 
  system.AddSolutionToSystemPDE("u");  
  
   // attach the assembling function to system
  system.SetAssembleFunction(AssembleFracFSI);

  // initilaize and solve the system
  system.init();
  system.solve();
 
   // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("u");
 
  VTKWriter vtkIO(&mlSol);
  vtkIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  


  return 0;
}





void AssembleFracFSI(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data

  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use

  LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ("FracFSI");   // pointer to the linear implicit system named "LiftRestr"
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const unsigned levelMax = mlPdeSys->GetLevelMax();
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
  vector < vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  for (unsigned i = 0; i < dim; i++) {
    x[i].reserve(maxSize);
  }
 //*************************** 

 //*************************** 
  double weight; // gauss point weight
  

 //******** Thom ******************* 
 //*************************** 
  vector <double> phi_u;  // local test function
  vector <double> phi_x_u; // local test function first order partial derivatives
  vector <double> phi_xx_u; // local test function second order partial derivatives

  phi_u.reserve(maxSize);
  phi_x_u.reserve(maxSize * dim);
  phi_xx_u.reserve(maxSize * dim2);
  
 
  unsigned solIndex_u;
  solIndex_u = mlSol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned solType_u = mlSol->GetSolutionType(solIndex_u);    // get the finite element type for "u"

  unsigned solPdeIndex_u;
  solPdeIndex_u = mlPdeSys->GetSolPdeIndex("u");    // get the position of "Thom" in the pdeSys object

  vector < double >  sol_u; // local solution
  sol_u.reserve(maxSize);
  vector< int > l2GMap_u;
  l2GMap_u.reserve(maxSize);
 //*************************** 
 //*************************** 

   

 //*************************** 
  //********* WHOLE SET OF VARIABLES ****************** 
  const int solType_max = 2;  //biquadratic

  const int n_vars = 1;
 
  vector< int > l2GMap_AllVars; // local to global mapping
  l2GMap_AllVars.reserve(n_vars*maxSize);
  
  vector< double > Res; // local redidual vector
  Res.reserve(n_vars*maxSize);

  vector < double > Jac;
  Jac.reserve( n_vars*maxSize * n_vars*maxSize);
 //*************************** 

  
 //********** DATA ***************** 
  double T_des = 100.;
  double alpha = 10.e5;
  double beta  = 1.;
  double gamma = 1.;
  
 //*************************** 
  
  
  if (assembleMatrix)  KK->zero();

    
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->IS_Mts2Gmt_elem_offset[iproc]; iel < msh->IS_Mts2Gmt_elem_offset[iproc + 1]; iel++) {

    unsigned kel = msh->IS_Mts2Gmt_elem[iel]; // mapping between paralell dof and mesh dof
    short unsigned kelGeom = el->GetElementType(kel);    // element geometry type

 //********* GEOMETRY ****************** 
    unsigned nDofx = el->GetElementDofNumber(kel, xType);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  x[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned iNode = el->GetMeshDof(kel, i, xType);    // local to global coordinates node
      unsigned xDof  = msh->GetMetisDof(iNode, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_coordinate->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }
 //*************************************** 
    
 //*********** Thom **************************** 
    unsigned nDof_u     = el->GetElementDofNumber(kel, solType_u);    // number of solution element dofs
    sol_u   .resize(nDof_u);
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
      unsigned iNode = el->GetMeshDof(kel, i, solType_u);    // local to global solution node
      unsigned solDof_u = msh->GetMetisDof(iNode, solType_u);    // global to global mapping between solution node and solution dof
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);      // global extraction and local storage for the solution
      l2GMap_u[i] = pdeSys->GetKKDof(solIndex_u, solPdeIndex_u, iNode);    // global to global mapping between solution node and pdeSys dof
    }
 //*********** Thom **************************** 

 
 //********** ALL VARS ***************** 
    unsigned nDof_AllVars = nDof_u; 
    const int nDof_max    = nDof_u;   // AAAAAAAAAAAAAAAAAAAAAAAAAAA TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    Res.resize(nDof_AllVars);
    std::fill(Res.begin(), Res.end(), 0.);

    Jac.resize(nDof_AllVars * nDof_AllVars);
    std::fill(Jac.begin(), Jac.end(), 0.);
    
    l2GMap_AllVars.resize(0);
    l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_u.begin(),l2GMap_u.end());
 //*************************** 


    if (level == levelMax || !el->GetRefinedElementIndex(kel)) {      // do not care about this if now (it is used for the AMR)
   
      // *** Gauss point loop ***
      for (unsigned ig = 0; ig < msh->_finiteElement[kelGeom][solType_max]->GetGaussPointNumber(); ig++) {
	
        // *** get gauss point weight, test function and test function partial derivatives ***
        //  ==== Thom 
	msh->_finiteElement[kelGeom][solType_u]   ->Jacobian(x, ig, weight, phi_u, phi_x_u, phi_xx_u);
	
        //FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {

          double srcTerm = 10.;
	  
          // FIRST ROW
	  Res[0                      + i] += weight * (srcTerm*phi_u[i]) ;

          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {
              double laplace_mat_u = 0.;

              for (unsigned kdim = 0; kdim < dim; kdim++) {
                 laplace_mat_u        += (phi_x_u   [i * dim + kdim] * phi_x_u   [j * dim + kdim]);
	      }

              //DIAG BLOCK
	        Jac[    0    * (nDof_u)    +  i    * (nDof_u) + (0 + j) ]  += weight * laplace_mat_u;
	      
            } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop
    } // endif single element not refined or fine grid loop

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
}