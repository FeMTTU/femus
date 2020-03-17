// /** started from file Ex6.cpp

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NumericVector.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "LinearImplicitSystem.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "Files.hpp"

#define FACE_FOR_CONTROL  1


#include   "../../nsopt_params.hpp"

#define exact_sol_flag 0 // 1 = if we want to use manufactured solution; 0 = if we use regular convention
#define compute_conv_flag 0 // 1 = if we want to compute the convergence and error ; 0 =  no error computation

#define NO_OF_NORMS 5 // for L2 norm of U,V,P and H1 norm of U,V

using namespace femus;


bool SetBoundaryConditionBox(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  //1: bottom  //2: right  //3: top  //4: left
 
  bool dirichlet = true;
  value = 0.;

#if exact_sol_flag == 0
// b.c. for lid-driven cavity problem, wall u_top = 1 = shear_force, v_top = 0 and u=v=0 on other 3 walls ; rhs_f = body_force = {0,0}
// TOP ==========================  
      if (facename == 3) {
       if (!strcmp(SolName, "U"))    { dirichlet = false; /*value =  1.;*/ } //lid - driven
  else if (!strcmp(SolName, "V"))    { dirichlet = false; /*value =  0.;*/} 
  	
      }
#endif

#if exact_sol_flag == 1
  //b.c. for manufactured lid driven cavity
// TOP ==========================  
   double pi = acos(-1.);
      if (facename == 3) {
       if (!strcmp(SolName, "U"))    { value =  sin(pi* x[0]) * sin(pi* x[0]) * cos(pi* x[1]) - sin(pi* x[0]) * sin(pi* x[0]); } //lid - driven
  else if (!strcmp(SolName, "V"))    { value = - sin(2. * pi * x[0]) * sin(pi* x[1]) + pi * x[1] * sin(2. * pi * x[0]);} 
  	
      }
#endif
  
      
  return dirichlet;
}


void AssembleNS_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );

void AssembleNS_nonAD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );

double*  GetErrorNorm(const MultiLevelProblem& ml_prob, MultiLevelSolution* mlSol, Solution* sol_coarser_prolongated);
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv 
//i.e. ||prol_(u_(i-1)) - u_(i)|| = err(i) => err(i-1)/err(i) = 2^alpha ,implemented as log(err(i)/err(i+1))/log2

void output_convergence_rate( double norm_i, double norm_ip1, std::string norm_name, unsigned maxNumberOfMeshes, int loop_i );


int main(int argc, char** args) {



  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  Files files;
        files.CheckIODirectories();
        files.RedirectCout();

    // ======= Quad Rule ========================
    std::string fe_quad_rule("seventh");

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  MultiLevelMesh mlMsh_all_levels;
 // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
   
  
    //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;
 // *** apparently needed by non-AD assemble only **********************
  // add fluid material
  Parameter parameter(Lref,Uref);
  
 
  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,1,FLUID_DENSITY,"Newtonian");
  std::cout << "Fluid properties: " << std::endl;
  std::cout << fluid << std::endl;

  
  mlMsh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,fe_quad_rule.c_str());
  mlMsh_all_levels.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,fe_quad_rule.c_str());
//   mlMsh.ReadCoarseMesh("./input/cube_hex.neu", fe_quad_rule.c_str(), scalingFactor);
//   //mlMsh.ReadCoarseMesh ( "./input/square_quad.neu", "seventh", scalingFactor );
  
  unsigned dim = mlMsh.GetDimension();
  unsigned maxNumberOfMeshes;

  if (dim == 2) {
    maxNumberOfMeshes = 1;
  } else {
    maxNumberOfMeshes = 4;
  }

 
    double comp_conv[maxNumberOfMeshes][NO_OF_NORMS];

        unsigned numberOfUniformLevels_finest = maxNumberOfMeshes;
        mlMsh_all_levels.RefineMesh(numberOfUniformLevels_finest, numberOfUniformLevels_finest, NULL);
//      mlMsh_all_levels.EraseCoarseLevels(numberOfUniformLevels - 2);  // need to keep at least two levels to send u_(i-1) projected(prolongated) into next refinement
        
        //store the fine solution  ==================
            MultiLevelSolution * mlSol_all_levels;
            mlSol_all_levels = new MultiLevelSolution (& mlMsh_all_levels);  //with the declaration outside and a "new" inside it persists outside the loop scopes
         // add variables to mlSol_all_levels
        // state =====================  
            mlSol_all_levels->AddSolution("U", LAGRANGE, SECOND);
            mlSol_all_levels->AddSolution("V", LAGRANGE, SECOND);
            if (dim == 3) mlSol_all_levels->AddSolution("W", LAGRANGE, SECOND);
            mlSol_all_levels->AddSolution("P", LAGRANGE, FIRST);
            mlSol_all_levels->Initialize("All");
            mlSol_all_levels->AttachSetBoundaryConditionFunction(SetBoundaryConditionBox);
            mlSol_all_levels->GenerateBdc("All");

         for (int i = 0; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

  unsigned numberOfUniformLevels = i + 1;
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

// #if PRESS == 1
  mlSol.AddSolution("P", LAGRANGE, FIRST);
// #endif

  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionBox);
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  mlProb.parameters.set<Fluid>("Fluid") = fluid;
  mlProb.SetQuadratureRuleAllGeomElems(fe_quad_rule);
  mlProb.SetFilesHandler(&files);

  // add system NS_fulldir in mlProb as a NonLinear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  if (dim == 3) system.AddSolutionToSystemPDE("W");

// #if PRESS == 1
  system.AddSolutionToSystemPDE("P");
// #endif

  // attach the assembling function to system
  system.SetAssembleFunction(AssembleNS_AD);
//   system.SetAssembleFunction(AssembleNS_nonAD);

  
   // initilaize and solve the system
  system.init();
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

//print the nonlinear iterations  
  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);
  system.SetDebugNonlinear(true);

  system.SetOuterSolver(PREONLY);
  system.MGsolve();
//   system.MGsolve();

  system.compute_convergence_rate();

  
    if ( i > 0 ) {
        
//prolongation of coarser  
      mlSol_all_levels->RefineSolution(i);
      Solution* sol_coarser_prolongated = mlSol_all_levels->GetSolutionLevel(i);
  
  
      double* norm = GetErrorNorm(mlProb,&mlSol,sol_coarser_prolongated);
    
      for(int j = 0; j < NO_OF_NORMS; j++)       comp_conv[i-1][j] = norm[j];
 
    }

    
//store the last computed solution
// 
       const unsigned level_index_current = 0;
      //@todo there is a duplicate function in MLSol: GetSolutionLevel() and GetLevel()
       const unsigned n_vars = mlSol.GetSolutionLevel(level_index_current)->_Sol.size();
       
        for(unsigned short j = 0; j < n_vars; j++) {  
               *(mlSol_all_levels->GetLevel(i)->_Sol[j]) = *(mlSol.GetSolutionLevel(level_index_current)->_Sol[j]);
        }
        
   
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

 mlSol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/,"biquadratic", variablesToBePrinted, i);
 
  //Destroy all the new systems
//   mlProb.clear();
 }

//   delete mlSol_all_levels; 

#if compute_conv_flag == 1
std::vector< std::string > norm_names = {"L2-NORM of U","L2-NORM of V", "L2-NORM of P" , "H1-Norm of U" , "H1-Norm of V"};

   for(int j = 0; j <  NO_OF_NORMS; j++)  {
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << norm_names[j] << " ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\t\t" << norm_names[j] << "\t\t\t\torder of convergence\n"; 
   for(int i = 0; i <  maxNumberOfMeshes - 1; i++){
       output_convergence_rate(comp_conv[i][j], comp_conv[i + 1][j], norm_names[j], maxNumberOfMeshes , i );
    }
  }
#endif
  
  return 0;
}


void output_convergence_rate( double norm_i, double norm_ip1, std::string norm_name, unsigned maxNumberOfMeshes , int loop_i) {

    std::cout << loop_i + 1 << "\t\t" <<  std::setw(11) << std::setprecision(10) << norm_i << "\t\t\t\t" ;
  
    if (loop_i < maxNumberOfMeshes/*norm.size()*/ - 2) {
      std::cout << std::setprecision(3) << log( norm_i/ norm_ip1 ) / log(2.) << std::endl;
    }
  
}



//manufactured solution for lid-driven---------------------------------------------
void value_Vel(const std::vector < double >& x, vector < double >& val_Vel) {
  double pi = acos(-1.);
  val_Vel[0] =   sin(pi* x[0]) * sin(pi* x[0]) * cos(pi* x[1]) - sin(pi* x[0]) * sin(pi* x[0]);
  val_Vel[1] = - sin(2. * pi * x[0]) * sin(pi* x[1]) + pi * x[1] * sin(2. * pi * x[0]);
 };
 
 
void gradient_Vel(const std::vector < double >& x, vector < vector < double > >& grad_Vel) {
  double pi = acos(-1.);
  grad_Vel[0][0]  =   pi * sin(2. * pi * x[0]) * cos(pi* x[1]) - pi * sin(2. * pi * x[0]);
  grad_Vel[0][1]  = - pi * sin(pi* x[0]) * sin(pi* x[0]) *  sin(pi * x[1]); 
  grad_Vel[1][0]  = - 2. * pi * cos(2. * pi * x[0]) * sin(pi* x[1]) + 2. * pi * pi * x[1] * cos(2. * pi * x[0]);   
  grad_Vel[1][1]  = - pi * sin(2. * pi * x[0]) * cos(pi * x[1]) + pi * sin(2. * pi * x[0]); 
 };

  
void laplace_Vel(const std::vector < double >& x, vector < double >& lap_Vel) {
  double pi = acos(-1.);
  lap_Vel[0] = - 2. * pi * pi * cos(2. * pi * x[0]) - 0.5 * pi * pi * cos(pi * x[1]) + 2.5 * pi * pi * cos(2. * pi* x[0]) * cos(pi* x[1]);
  lap_Vel[1] = - 4. * pi * pi * pi * x[1] * sin(2. * pi * x[0]) + 5. * pi * pi * sin(2. * pi * x[0]) * sin(pi * x[1]);
};

double value_Press(const std::vector < double >& x) {
  double pi = acos(-1.);
  return sin(2. * pi * x[0]) * sin(2. * pi * x[1]); //p
 };

 void gradient_Press(const std::vector < double >& x, vector < double >& grad_statePress) {
  double pi = acos(-1.);
  grad_statePress[0]  =   2. * pi * cos(2. * pi * x[0]) * sin(2. * pi * x[1]); 
  grad_statePress[1]  =   2. * pi * sin(2. * pi * x[0]) * cos(2. * pi * x[1]);
 };

//manufactured solution for lid-driven---------------------------------------------



void AssembleNS_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem& mlPdeSys   = ml_prob.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the nonlinear implicit system named "NS"
  const unsigned level = mlPdeSys.GetLevelToAssemble();
  bool assembleMatrix = mlPdeSys.GetAssembleMatrix(); 

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys._LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

// total number of variables
  unsigned n_unknowns = dim + 1;

// #if PRESS == 1
//   n_unknowns += 1;
// #endif  

  //solution variable
  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if (dim == 3) solVIndex[2] = mlSol->GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"


  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys.GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys.GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object

  if (dim == 3) solVPdeIndex[2] = mlPdeSys.GetSolPdeIndex("W");

// #if PRESS == 1
  unsigned solPIndex;
  solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"
  
  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys.GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object
// #endif

//   const int n_unknowns = n_vars;  //state velocity terms and one pressure term
  const int vel_type_pos = 0;
  const int press_type_pos = dim;
  const int state_pos_begin = 0;
  
  vector < std::string > Solname(n_unknowns);  // const char Solname[4][8] = {"U","V","W","P"};
  Solname              [state_pos_begin+0] =                "U";
  Solname              [state_pos_begin+1] =                "V";
  if (dim == 3) Solname[state_pos_begin+2] =                "W";
// #if PRESS == 1
  Solname              [state_pos_begin + press_type_pos] = "P";
// #endif  
  
  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar]	= mlPdeSys.GetSolPdeIndex(Solname[ivar].c_str());
    SolIndex[ivar]	= mlSol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= mlSol->GetSolutionType(SolIndex[ivar]);
  }

  

  vector < vector < adept::adouble > >  solV(dim);    // local solution
  vector< vector < adept::adouble > > aResV(dim);    // local redidual vector

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    aResV[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

// #if PRESS == 1
  vector < adept::adouble >  solP; // local solution
  vector< adept::adouble > aResP; // local redidual vector
  solP.reserve(maxSize);
  aResP.reserve(maxSize);
// #endif


  double weight; // gauss point weight
  
  vector <double> phiV;  // local test function
  vector <double> phiV_x; // local test function first order partial derivatives
  vector <double> phiV_xx; // local test function second order partial derivatives

  phiV.reserve(maxSize);
  phiV_x.reserve(maxSize * dim);
  phiV_xx.reserve(maxSize * dim2);

// #if PRESS == 1
  double* phiP;
// #endif

  //Nondimensional values ******************
  double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  //Nondimensional values ******************
  
  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve( n_unknowns *maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve( n_unknowns *maxSize);

  vector < double > Jac;
  Jac.reserve( n_unknowns *maxSize * n_unknowns *maxSize);

  RES->zero(); // Set to zero all the entries of the Global Matrix
  if (assembleMatrix) KK->zero(); // Set to zero all the entries of the Global Matrix

// #if PRESS == 1
//     for (unsigned  k = 0; k < n_unknowns; k++) {
//   std::cout << "******************" << std::endl;
//         sol->_Sol[ SolIndex[k]]->print();
//     }
// #endif  
  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
// #if PRESS == 1
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
// #endif
    
    unsigned nDofsVP = dim * nDofsV + nDofsP;

// #if PRESS == 1
//     nDofsVP += nDofsP;
// #endif

    // resize local arrays
    sysDof.resize(nDofsVP);

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }

// #if PRESS == 1
    solP.resize(nDofsP);
// #endif

    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].resize(nDofsV);    //resize
      std::fill(aResV[k].begin(), aResV[k].end(), 0);    //set aRes to zero
    }

// #if PRESS == 1
    aResP.resize(nDofsP);    //resize
    std::fill(aResP.begin(), aResP.end(), 0);    //set aRes to zero
// #endif

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[i + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

// #if PRESS == 1
    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // global to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[i + dim * nDofsV] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
// #endif

    // local storage of coordinates
    for (unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);
// #if PRESS == 1
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);
// #endif
      
      vector < adept::adouble > solV_gss(dim, 0);
      vector < vector < adept::adouble > > gradSolV_gss(dim);
      vector < double > coordX_gss(dim, 0.);

      for (unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].resize(dim);
        std::fill(gradSolV_gss[k].begin(), gradSolV_gss[k].end(), 0);
      }

      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += phiV[i] * solV[k][i];    
          coordX_gss[k] += coordX[k][i] * phiV[i];
        }

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += phiV_x[i * dim + j] * solV[k][i];
          }
        }
      }

// #if PRESS == 1
      adept::adouble solP_gss = 0;
      for (unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }
// #endif



//computation of RHS force using MMS=============================================== 
vector <double>  exact_Vel(dim,0.);
value_Vel(coordX_gss,exact_Vel);
vector < vector < double > > exact_grad_Vel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_Vel[k].resize(dim);
    std::fill(exact_grad_Vel[k].begin(), exact_grad_Vel[k].end(), 0.);
}
gradient_Vel(coordX_gss,exact_grad_Vel);
vector <double>  exact_lap_Vel(dim,0.);
laplace_Vel(coordX_gss, exact_lap_Vel);
vector <double>  exact_conv_Vel(dim,0.);
vector <double> exact_grad_Press(dim,0.);
gradient_Press(coordX_gss, exact_grad_Press);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_Vel[k] += exact_grad_Vel[k][i] * exact_Vel[i] ; 
    }
}


vector <double> exactForce(dim,0.);
for (unsigned k = 0; k < dim; k++){
    exactForce[k] =  - IRe * exact_lap_Vel[k] + advection_flag * exact_conv_Vel[k] + exact_grad_Press[k] ;
}
//computation of RHS force using MMS=============================================== 


        // *** phiV_i loop ***
      for (unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > NSV(dim, 0.);

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            NSV[k]   +=  IRe * phiV_x[i * dim + j] * (gradSolV_gss[k][j]/* + gradSolV_gss[j][k]*/)
                        + advection_flag * phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]);
          }
        }

// #if PRESS == 1
        for (unsigned  k = 0; k < dim; k++) {
          NSV[k] += - solP_gss * phiV_x[i * dim + k];
        }
// #endif

        for (unsigned  k = 0; k < dim; k++) {
#if exact_sol_flag == 0
           aResV[k][i] += ( + NSV[k] - force[k] * phiV[i]  ) * weight;
#endif
#if exact_sol_flag == 1
          aResV[k][i] += ( + NSV[k] - exactForce[k] * phiV[i]  ) * weight;
#endif
       }
      } // end phiV_i loop

// #if PRESS == 1
      // *** phiP_i loop ***
      for (unsigned i = 0; i < nDofsP; i++) {
        for (int k = 0; k < dim; k++) {
          aResP[i] += - (gradSolV_gss[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop
// #endif

    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsVP);    //resize

    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        Res[ i +  k * nDofsV ] = - aResV[k][i].value();
      }
    }

// #if PRESS == 1
    for (int i = 0; i < nDofsP; i++) {
      Res[ i + dim * nDofsV ] = - aResP[i].value();
    }
// #endif

    RES->add_vector_blocked(Res, sysDof);

 if (assembleMatrix){   //Extarct and store the Jacobian

    Jac.resize(nDofsVP * nDofsVP);
    // define the dependent variables

    for (unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResV[k][0], nDofsV);
    }

// #if PRESS == 1
    s.dependent(&aResP[0], nDofsP);
// #endif
    
    // define the independent variables
    for (unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofsV);
    }

// #if PRESS == 1
    s.independent(&solP[0], nDofsP);
// #endif
    
    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0] , true);
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    s.clear_independents();
    s.clear_dependents();
 }  //end assemble matrix
    
    
  } //end element loop for each process


 if (assembleMatrix){   //Extarct and store the Jacobian
  KK->close();
  }
 
  RES->close();
//   RES->print();
//    std::cout << "solution iterate RESC" << std::endl;
//   pdeSys->_RESC->print();
// 
//   std::cout << "solution iterate EPS" << std::endl;
//   pdeSys->_EPS->print();
//   std::cout << "solution iterate EPSC" << std::endl;
//   pdeSys->_EPSC->print();
  // ***************** END ASSEMBLY *******************
}


void AssembleNS_nonAD(MultiLevelProblem& ml_prob){
     
  //pointers
  NonLinearImplicitSystem& mlPdeSys  = ml_prob.get_system<NonLinearImplicitSystem>("NS");
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
  vector< vector < double> > coordX(dim);	//local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  for(int i=0;i<dim;i++) {   
       coordX[i].reserve(maxSize); 
  }
  // geometry *******************************************

 
  
  // solution variables *******************************************
  int n_vars = dim + 1;
// #if PRESS == 1
//   n_vars += 1;
// #endif  
  const int n_unknowns = n_vars;  //state velocity terms and one pressure term
  const int vel_type_pos = 0;
  const int press_type_pos = dim;
  const int state_pos_begin = 0;
  
  vector < std::string > Solname(n_unknowns);  // const char Solname[4][8] = {"U","V","W","P"};
  Solname              [state_pos_begin+0] =                "U";
  Solname              [state_pos_begin+1] =                "V";
  if (dim == 3) Solname[state_pos_begin+2] =                "W";
// #if PRESS == 1
  Solname              [state_pos_begin + press_type_pos] = "P";
// #endif  
  
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
  vector < vector < double > > phi_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_xx_gss_fe(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(maxSize);
      phi_x_gss_fe[fe].reserve(maxSize*dim);
     phi_xx_gss_fe[fe].reserve(maxSize*(3*(dim-1)));
   }
  
  //=================================================================================================
  
  // quadratures ********************************
  double weight;
  
  
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
  
// #if PRESS == 1
//     for (unsigned  k = 0; k < n_unknowns; k++) {
//   std::cout << "******************" << std::endl;
//   sol->_Sol[SolIndex[k]]->print();
//     }
// #endif  
  
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
// #if PRESS == 1
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType[state_pos_begin + press_type_pos]);    // number of solution element dofs
// #endif
    
    unsigned nDofsVP = dim * nDofsV + nDofsP;
// #if PRESS == 1
//     nDofsVP += nDofsP;
// #endif
    // equation end *****************************
  
  
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

   
      // ********************** Gauss point loop *******************************
      for(unsigned ig=0;ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {
	
 
      for(int fe=0; fe < NFE_FAMS; fe++) {
	ml_prob._ml_msh->_finiteElement[ielGeom][fe]->Jacobian(coordX,ig,weight,phi_gss_fe[fe],phi_x_gss_fe[fe],phi_xx_gss_fe[fe]);
      }
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
  	ml_prob._ml_msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coordX,ig,weight,phi_gss_fe[BIQUADR_FE],phi_x_gss_fe[BIQUADR_FE],phi_xx_gss_fe[BIQUADR_FE]);

      vector < double > coordX_gss(dim, 0.);
 	for(unsigned k = 0; k <  dim; k++) {
	  for(unsigned i = 0; i < Sol_n_el_dofs[k]; i++) {
         coordX_gss[k] += coordX[k][i] * phi_gss_fe[ SolFEType[k] ][i];
      }
    }

    
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

// // //  //residuals and Jac------------------------------------------------------------------------------------------------
// // // //==========FILLING WITH THE EQUATIONS =========================================================================================================
// // // for(unsigned i_unk=0; i_unk<n_unknowns; i_unk++) { 
// // //     for (unsigned i = 0; i < Sol_n_el_dofs[i_unk]; i++) {
// // // 	    double div_u_du_qp = 0.;
// // // 	    double lap_res_du_u = 0.; 
// // // 	    double p_div_du_qp = 0.;
// // // // 	    double adv_rhs = 0.;
// // // 	    
// // // 	for (unsigned jdim = 0; jdim < dim; jdim++) {
// // // 	  if ( i_unk==0 || i_unk==1 ){	      lap_res_du_u  +=  gradSolVAR_qp[i_unk][jdim]*phi_x_gss_fe[ SolFEType[i_unk] ][i * dim + jdim];
// // // // 						adv_rhs	+= SolVAR_qp[jdim] * gradSolVAR_qp[i_unk][j_dim];
// // // 	  }
// // // #if PRESS == 1
// // //                 p_div_du_qp += SolVAR_qp[SolIndex[press_type_pos]] * phi_x_gss_fe[ SolFEType[i_unk] ][i * dim + jdim];
// // // //div--------------------------
// // // 	  	div_u_du_qp += gradSolVAR_qp[SolPdeIndex[jdim]][jdim] ;  //kdims are with jdims  
// // // #endif
// // // 	}//jdim 
// // // 
// // //         
// // // //======================Residuals===================================================================================================================
// // //     // FIRST ROW
// // // 	  if (i_unk==0 || i_unk==1)       	   	  	 Res[i_unk][i]  -=  (   + force[i_unk] * phi_gss_fe[ SolFEType[i_unk] ][i]
// // // 										    - IRe*lap_res_du_u 
// // // // 										    - adv_rhs * phi_gss_fe[ SolFEType[i_unk] ][i]
// // // #if PRESS == 1
// // // 										    + /*p_div_du_qp*/SolVAR_qp[SolIndex[press_type_pos]] * phi_x_gss_fe[ SolFEType[i_unk] ][i * dim + i_unk] 
// // // #endif
// // // 	  ) * weight; 
// // //   
// // // #if PRESS == 1
// // // 	  if (i_unk==2)       	       				 Res[i_unk][i]  -=  ( (div_u_du_qp) * phi_gss_fe[ SolFEType[i_unk] ][i] ) * weight;
// // // #endif
// // //  
// // //     
// // // // // //       }//kdim_Res
// // // 
// // // 	 
// // // 	 
// // // //======================Jacobian========================================================================================================================
// // // 	      
// // //    if (assembleMatrix) {
// // //     for(unsigned j_unk=0; j_unk< n_unknowns; j_unk++) { 
// // // 	for (unsigned j = 0; j < Sol_n_el_dofs[j_unk]; j++) {
// // // 	            double lap_jac_du_u = 0.;
// // // 	      
// // // 		for (unsigned kdim = 0; kdim < dim; kdim++) {
// // // 		  if ( i_unk==j_unk && (i_unk==0 ||i_unk==1) ) 	lap_jac_du_u += phi_x_gss_fe[ SolFEType[i_unk] ][i * dim + kdim]*phi_x_gss_fe[ SolFEType[j_unk] ][j * dim + kdim];
// // // 		}//kdim
// // // 	     
// // //     //============ delta_state row ============================================================================================
// // //        //DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
// // //     // FIRST ROW
// // // 		  if ( i_unk==j_unk && (i_unk==0 ||i_unk==1))          Jac[i_unk][j_unk][i*nDofsV + j] -=  (  IRe*lap_jac_du_u ) * weight; 
// // // #if PRESS == 1
// // // 		  if ((i_unk==0 ||i_unk==1) && j_unk==2)               Jac[i_unk][j_unk][i*nDofsP + j] -= -( phi_gss_fe[ SolFEType[j_unk] ][j] * phi_x_gss_fe[ SolFEType[i_unk] ][i * dim + i_unk] ) * weight;
// // // 		  if ( i_unk==2  && (j_unk==0 ||j_unk==1))             Jac[i_unk][j_unk][i*nDofsV + j] -= -( phi_gss_fe[ SolFEType[i_unk] ][i] * phi_x_gss_fe[ SolFEType[j_unk] ][j * dim + j_unk] ) * weight;
// // // #endif
// // // 	} //end j loop
// // //     } //end j_unk loop
// // //   } // endif assemble_matrix
// // // 
// // //     } // end i loop
// // // } // end i_unk loop



//computation of RHS force using MMS=============================================== 
vector <double>  exact_Vel(dim,0.);
value_Vel(coordX_gss,exact_Vel);
vector < vector < double > > exact_grad_Vel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_Vel[k].resize(dim);
    std::fill(exact_grad_Vel[k].begin(), exact_grad_Vel[k].end(), 0.);
}
gradient_Vel(coordX_gss,exact_grad_Vel);
vector <double>  exact_lap_Vel(dim,0.);
laplace_Vel(coordX_gss, exact_lap_Vel);
vector <double>  exact_conv_Vel(dim,0.);
vector <double> exact_grad_Press(dim,0.);
gradient_Press(coordX_gss, exact_grad_Press);

for (unsigned k = 0; k < dim; k++){
    for (unsigned i = 0; i < dim; i++){
    exact_conv_Vel[k] += exact_grad_Vel[k][i] * exact_Vel[i] ; 
    }
}


vector <double> exactForce(dim,0.);
for (unsigned k = 0; k < dim; k++){
    exactForce[k] =  - IRe * exact_lap_Vel[k] + advection_flag * exact_conv_Vel[k] + exact_grad_Press[k] ;
}
//computation of RHS force using MMS=============================================== 

 
//good old method for filling residuals and Jac  
//============ delta_state row ============================================================================================

  for (unsigned i = 0; i < nDofsV; i++) {
// FIRST ROW
	for (unsigned  kdim = 0; kdim < dim; kdim++) { // velocity block row 
	              double lap_res_du_u = 0.; 
		      double adv_res = 0.;
	      for (unsigned jdim = 0; jdim < dim; jdim++) {
		    lap_res_du_u += gradSolVAR_qp[SolPdeIndex[kdim]][jdim]*phi_x_gss_fe[ SolFEType[kdim] ][i * dim + jdim];
			adv_res	+= SolVAR_qp[jdim] * gradSolVAR_qp[kdim][jdim];
	      }      
	      Res[kdim][i]   +=  (         
#if exact_sol_flag == 0
                                         + force[kdim] * phi_gss_fe[ SolFEType[kdim] ][i]
 #endif                                      
 #if exact_sol_flag == 1
                                       + exactForce[kdim] * phi_gss_fe[ SolFEType[kdim] ][i]
 #endif
                                          - IRe*lap_res_du_u 
                                           - advection_flag * adv_res * phi_gss_fe[ SolFEType[kdim] ][i]
					    + SolVAR_qp[SolPdeIndex[press_type_pos]] * phi_x_gss_fe[ SolFEType[kdim] ][i * dim + kdim]) * weight; 
	}	    
//DIAG BLOCK delta_state - state--------------------------------------------------------------------------------
	for (unsigned j = 0; j < nDofsV; j++) {
		      double lap_jac_du_u = 0.;
		      vector < double > adv_uold_nablaunew(dim,0.);
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
		    lap_jac_du_u += phi_x_gss_fe[ SolFEType[kdim] ][i * dim + kdim]*phi_x_gss_fe[ SolFEType[kdim] ][j * dim + kdim];
		for (unsigned  jdim = 0; jdim < dim; jdim++) { 
		    adv_uold_nablaunew[kdim] += SolVAR_qp[SolIndex[jdim]]*phi_x_gss_fe[ SolFEType[kdim] ][j * dim + jdim] * phi_gss_fe[ SolFEType[kdim] ][i];
                }  //jdim
	      }
	      
	      for (unsigned  kdim = 0; kdim < dim; kdim++) { 
		Jac[kdim][kdim][i*nDofsV + j] += (   IRe*lap_jac_du_u 
						    + advection_flag * adv_uold_nablaunew[kdim] //this part goes only in diagonal blocks
						    + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim] ][j] * gradSolVAR_qp[SolIndex[kdim]][kdim] * phi_gss_fe[ SolFEType[kdim] ][i] // diagonal blocks of adv_unew_nablauold
						    ) * weight; 

               unsigned int kdim_p1 = (kdim+1)%dim; //off-diagonal blocks of adv_unew_nablauold
                Jac[kdim][kdim_p1][i*nDofsV + j] +=  + advection_flag * (1 - advection_Picard) * phi_gss_fe[ SolFEType[kdim_p1] ][j] * gradSolVAR_qp[SolIndex[kdim]][kdim_p1] * phi_gss_fe[ SolFEType[kdim] ][i]  * weight;
						    
	      }
	} //j_du_u loop


     
//BLOCK Pressure
      for (unsigned j = 0; j < nDofsP; j++) {
	    for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      Jac[kdim][press_type_pos][i*nDofsP + j] += - ( phi_gss_fe[ SolFEType[press_type_pos] ][j] * phi_x_gss_fe[ SolFEType[kdim] ][i * dim + kdim] ) * weight;
	    }
      }//j_press loop
   }//i_state loop

//DIV_state
  for (unsigned i = 0; i < nDofsP; i++) {
		    double div_u_qp =0.;
      for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      div_u_qp += gradSolVAR_qp[SolPdeIndex[kdim]][kdim] ;
      }
      Res[press_type_pos][i]  +=  ( (div_u_qp) * phi_gss_fe[ SolFEType[press_type_pos] ][i] ) * weight;
      for (unsigned j = 0; j < nDofsV; j++) {
	  for (unsigned  kdim = 0; kdim < dim; kdim++) {
	      Jac[press_type_pos][kdim][i*nDofsV + j] += - ( phi_gss_fe[ SolFEType[press_type_pos] ][i] * phi_x_gss_fe[ SolFEType[kdim] ][j * dim + kdim] ) * weight;
	  }
      } //j loop
   }//i_div_state
    //============ delta_state row ============================================================================================

 
 
 
 
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
//    if(mlPdeSys.GetNonlinearIt() == 0){
//      std::ostringstream mat_out; mat_out << "matrix_non_ad" << mlPdeSys.GetNonlinearIt()  << ".txt";
//   JAC->print_matlab(mat_out.str(),"ascii");}
  RES->close();
//   RES->print();
//   std::cout << "solution iterate RESC" << std::endl;
//   pdeSys->_RESC->print();
// 
//   std::cout << "solution iterate EPS" << std::endl;
//   pdeSys->_EPS->print();
//   std::cout << "solution iterate EPSC" << std::endl;
//   pdeSys->_EPSC->print();
// ***************** END ASSEMBLY *******************
}



double*  GetErrorNorm(const MultiLevelProblem& ml_prob, MultiLevelSolution* mlSol, Solution* sol_coarser_prolongated) {
  
    static double ErrorNormArray[NO_OF_NORMS];
    
  unsigned level = mlSol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = mlSol->_mlMesh->GetLevel(level);    // pointer to the mesh (level) object
  elem*     el  = msh->el;  // pointer to the elem object in msh (level)
  Solution* sol = mlSol->GetSolutionLevel(level);    // pointer to the solution (level) object

  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)
  
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
   
  //geometry *******************************

 // solution variables *******************************************
  int n_vars = dim + 1;
// #if PRESS == 1
//   n_vars += 1;
// #endif  
  const int n_unknowns = n_vars;  //state velocity terms and one pressure term
  const int vel_type_pos = 0;
  const int press_type_pos = dim;
  const int state_pos_begin = 0;
  
  vector < std::string > Solname(n_unknowns);  // const char Solname[4][8] = {"U","V","W","P"};
  Solname              [state_pos_begin+0] =                "U";
  Solname              [state_pos_begin+1] =                "V";
  if (dim == 3) Solname[state_pos_begin+2] =                "W";
// #if PRESS == 1
  Solname              [state_pos_begin + press_type_pos] = "P";
// #endif  
  
  vector < unsigned > SolIndex(n_unknowns);  
  vector < unsigned > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolIndex[ivar]	= mlSol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= mlSol->GetSolutionType(SolIndex[ivar]);
  }

  vector < double > Sol_n_el_dofs(n_unknowns);
  
  //==========================================================================================
  // velocity ************************************
  vector < vector < double > > phi_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_x_gss_fe(NFE_FAMS);
  vector < vector < double > > phi_xx_gss_fe(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_gss_fe[fe].reserve(maxSize);
      phi_x_gss_fe[fe].reserve(maxSize*dim);
     phi_xx_gss_fe[fe].reserve(maxSize*(3*(dim-1)));
   }
  
  //=================================================================================================
  
  // quadratures ********************************
  double weight;
  
 
  //----------- dofs ------------------------------
  vector < vector < double > > SolVAR_eldofs(n_unknowns);
  vector < vector < double > > gradSolVAR_eldofs(n_unknowns);
  
  vector < vector < double > > SolVAR_coarser_prol_eldofs(n_unknowns);
  vector < vector < double > > gradSolVAR_coarser_prol_eldofs(n_unknowns);


  for(int k = 0; k < n_unknowns; k++) {
    SolVAR_eldofs[k].reserve(maxSize);
    gradSolVAR_eldofs[k].reserve(maxSize*dim); 
    
    SolVAR_coarser_prol_eldofs[k].reserve(maxSize);
    gradSolVAR_coarser_prol_eldofs[k].reserve(maxSize*dim);    
  }

  //------------ at quadrature points ---------------------
  vector < double > SolVAR_qp(n_unknowns);
  vector < double > SolVAR_coarser_prol_qp(n_unknowns);
  vector < vector < double > > gradSolVAR_qp(n_unknowns);
  vector < vector < double > > gradSolVAR_coarser_prol_qp(n_unknowns);
  for(int k = 0; k < n_unknowns; k++) {
      gradSolVAR_qp[k].reserve(maxSize);  
      gradSolVAR_coarser_prol_qp[k].reserve(maxSize);  
  }
      
  vector  < double > l2norm (NO_OF_NORMS,0.);

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    
  // geometry *****************************
    short unsigned ielGeom = msh->GetElementType(iel);
    
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    for (unsigned  k = 0; k < dim; k++) {       coordX[k].resize(nDofsX);    }
  
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
  
  // geometry end *****************************
  
  
  // equation *****************************
    unsigned nDofsV = msh->GetElementDofNumber(iel, SolFEType[vel_type_pos]);    // number of solution element dofs
// #if PRESS == 1
    unsigned nDofsP = msh->GetElementDofNumber(iel, SolFEType[state_pos_begin + press_type_pos]);    // number of solution element dofs
// #endif
    
    unsigned nDofsVP = dim * nDofsV + nDofsP;
// #if PRESS == 1
//     nDofsVP += nDofsP;
// #endif
    // equation end *****************************


   //STATE###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	Sol_n_el_dofs[k]=ndofs_unk;
       SolVAR_eldofs[k].resize(ndofs_unk);
       SolVAR_coarser_prol_eldofs[k].resize(ndofs_unk);
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);    // global to global mapping between solution node and solution dof // via local to global solution node
       SolVAR_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
       SolVAR_coarser_prol_eldofs[k][i] = (*sol_coarser_prolongated->_Sol[SolIndex[k]])(solDof);      // global extraction and local storage for the solution
      }
    }
  //CTRL###################################################################

 
      // ********************** Gauss point loop *******************************
      for(unsigned ig=0;ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {
	
 
      for(int fe=0; fe < NFE_FAMS; fe++) {
	msh->_finiteElement[ielGeom][fe]->Jacobian(coordX,ig,weight,phi_gss_fe[fe],phi_x_gss_fe[fe],phi_xx_gss_fe[fe]);
      }
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
  	msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coordX,ig,weight,phi_gss_fe[BIQUADR_FE],phi_x_gss_fe[BIQUADR_FE],phi_xx_gss_fe[BIQUADR_FE]);

 //begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk <  n_unknowns; unk++) {
	  SolVAR_qp[unk] = 0.;
	  SolVAR_coarser_prol_qp[unk] = 0.;
      gradSolVAR_qp[unk].resize(dim);  
      gradSolVAR_coarser_prol_qp[unk].resize(dim);  
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_qp[unk][ivar2] = 0.; 
	    gradSolVAR_coarser_prol_qp[unk][ivar2] = 0.; 
	  }
    }
	  
      vector < double > coordX_gss(dim, 0.);
 	for(unsigned k = 0; k <  dim; k++) {
	  for(unsigned i = 0; i < Sol_n_el_dofs[k]; i++) {
         coordX_gss[k] += coordX[k][i] * phi_gss_fe[ SolFEType[k] ][i];
      }
    }
	for(unsigned unk = 0; unk <  n_unknowns; unk++) {
	  for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
	    SolVAR_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i] * SolVAR_eldofs[unk][i];
	    SolVAR_coarser_prol_qp[unk] += phi_gss_fe[ SolFEType[unk] ][i] * SolVAR_coarser_prol_eldofs[unk][i];
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim+ivar2] * SolVAR_eldofs[unk][i]; 
	      gradSolVAR_coarser_prol_qp[unk][ivar2] += phi_x_gss_fe[ SolFEType[unk] ][i*dim+ivar2] * SolVAR_coarser_prol_eldofs[unk][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************

#if exact_sol_flag == 1
// exact solution error norm ========================================================
vector <double>  exact_Vel(dim,0.);
value_Vel(coordX_gss,exact_Vel);
double exact_Press = value_Press(coordX_gss);
vector < vector < double > > exact_grad_Vel(dim);
for (unsigned k = 0; k < dim; k++){ 
    exact_grad_Vel[k].resize(dim);
    std::fill(exact_grad_Vel[k].begin(), exact_grad_Vel[k].end(), 0.);
}
gradient_Vel(coordX_gss,exact_grad_Vel);

        for(unsigned unk = 0; unk <  dim; unk++) {
                l2norm[unk] += ( SolVAR_qp[unk] - exact_Vel[unk] ) * ( SolVAR_qp[unk] - exact_Vel[unk] ) * weight ; 
        }
                l2norm[dim] += ( SolVAR_qp[dim] - exact_Press ) * ( SolVAR_qp[dim] - exact_Press ) * weight ; 
        for(unsigned unk = 0; unk <  dim; unk++) {
            for(int j = 0; j < dim; j++){
                l2norm[n_unknowns + unk] += (gradSolVAR_qp[unk][j] - exact_grad_Vel[unk][j] ) * ( gradSolVAR_qp[unk][j] - exact_grad_Vel[unk][j] ) * weight ;
        }
 } //seminorm
// exact solution error norm ========================================================
#endif   

#if exact_sol_flag == 0
    for(unsigned unk = 0; unk <  n_unknowns; unk++) {
        l2norm[unk] += ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * ( SolVAR_qp[unk] - SolVAR_coarser_prol_qp[unk] ) * weight ; 
    } //l2norm
     
    for(unsigned unk = 0; unk <  dim; unk++) {
        for(int j = 0; j < dim; j++){
    l2norm[n_unknowns + unk] += (gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * ( gradSolVAR_qp[unk][j] - gradSolVAR_coarser_prol_qp[unk][j] ) * weight ;
        }
    } //seminorm
#endif   
     
    } // end gauss point loop
  } //end element loop for each process


    // add the norms of all processes
  NumericVector* norm_vec_inexact;
  norm_vec_inexact = NumericVector::build().release();
  norm_vec_inexact->init(msh->n_processors(), 1 , false, AUTOMATIC);

	for(unsigned unk = 0; unk < NO_OF_NORMS; unk++) {
        norm_vec_inexact->set(iproc, l2norm[unk]);
        norm_vec_inexact->close();
        l2norm[unk] = norm_vec_inexact->l1_norm();
    }


  delete norm_vec_inexact;
  
 
	for(unsigned unk = 0; unk < NO_OF_NORMS; unk++) {
        ErrorNormArray[unk] = sqrt(l2norm[unk]);
    }
   
   return ErrorNormArray;
  
  
}

