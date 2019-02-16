#include "MultiLevelProblem.hpp"
#include "TransientSystem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "FElemTypeEnum.hpp"
#include "Files.hpp"




using std::cout;
using std::endl;
using namespace femus;

void AssembleMatrixRes(MultiLevelProblem &ml_prob);

//-------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const std::vector < double >& x,const char name[],
			  double &value, const int FaceName, const double time){
  bool test = 1; //Dirichlet
  value = 0.;

  if(!strcmp(name,"u")) {
    if (1==FaceName) { //inflow
      test=1;
     value=0.;
    }
    else if(2==FaceName ){  //outflow
     test=1;
     value=0.;
    }
    else if(3==FaceName ){  // no-slip fluid wall
      test=1;
      value=0.;
    }
    else if(4==FaceName ){  // no-slip solid wall
      test=1;
      value=0.;
    }
  }


  return test;
}


double SetInitialCondition (const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {
         
           double value = 0.;

             if(!strcmp(name,"u")) {
                 value = sin(M_PI * x[0]) * sin(M_PI * x[1]);
             }

           
      return value;   
}


int main(int argc,char **args) {

    // ======= Init ========================
  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
    // ======= Files ========================
  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();

  // ======= Quad Rule ===================
  std::string fe_quad_rule("seventh");

  // ======= Mesh ========================
  //Nondimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;
  const unsigned int nsub_x = 16;
  const unsigned int nsub_y = 16;
  const unsigned int nsub_z = 0;
  const std::vector<double> xyz_min = {0.,0.,0.};
  const std::vector<double> xyz_max = {1.,1.,0.};
  const ElemType geom_elem_type = QUAD9;

  MultiLevelMesh ml_msh;
  ml_msh.GenerateCoarseBoxMesh(nsub_x,nsub_y,nsub_z,xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  ml_msh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  ml_msh.PrintInfo();

  // ======= Solution ========================
  MultiLevelSolution ml_sol(&ml_msh);  // define the multilevel solution and attach the mlMsh object to it
  const unsigned int time_dep_flag = 2;
  ml_sol.AddSolution("u", LAGRANGE, FIRST, time_dep_flag);

    // ======= Problem ========================
  MultiLevelProblem ml_prob(&ml_sol);  // define the multilevel problem attach the ml_sol object to it

  ml_prob.SetFilesHandler(&files);
  
    // ======= Initial values ========================
  ml_sol.Initialize("All");    // initialize all variables to zero
  ml_sol.Initialize("u",   SetInitialCondition, &ml_prob);

    // ======= Boundary Conditions ========================
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);  // attach the boundary condition function and generate boundary data

//   ml_sol.GenerateBdc("All");  //this would do it also for the non-equation-related variables
  ml_sol.GenerateBdc("u"); //"Time_dependent");

    // ======= Parameters ========================
  Parameter parameter(Lref,Uref);


  // ======= System ========================
  TransientNonlinearImplicitSystem & system = ml_prob.add_system<TransientNonlinearImplicitSystem> ("Timedep");
  system.AddSolutionToSystemPDE("u");

  system.init();  //does it have to stay here or later?

  system.SetAssembleFunction(AssembleMatrixRes);
  system.SetMaxNumberOfLinearIterations(1);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-8);
  system.SetMgType(V_CYCLE);
  system.SetMaxNumberOfNonLinearIterations(15);

  //**************
   ml_sol.SetWriter(VTK);   //need to move this here for the DebugNonlinear function
//    ml_sol.GetWriter()->SetDebugOutput(true);
  //**************
  
  // time loop parameter
  system.SetIntervalTime(0.001);
  const unsigned int n_timesteps = 100;
  const unsigned int write_interval = 1;

    // ======= Time Loop ========================
  for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {

    // ======= Final Print ========================
    if ( !(time_step%write_interval) ) {

  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  ml_sol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/, "biquadratic", variablesToBePrinted,time_step);    // print solutions

    }
    
    // ======= Solve ========================
    std::cout << std::endl;
    std::cout << " *********** Timedep ************  " << std::endl;
    ml_prob.get_system("Timedep").MLsolve();

    //update Solution
    ml_prob.get_system<TransientNonlinearImplicitSystem>("Timedep").CopySolutionToOldSolution();


  } //end loop timestep


  // Destroy all the new systems
  ml_prob.clear();


  return 0;
}




//------------------------------------------------------------------------------------------------------------
void AssembleMatrixRes(MultiLevelProblem &ml_prob){

  TransientNonlinearImplicitSystem* mlPdeSys = & ml_prob.get_system<TransientNonlinearImplicitSystem>("Timedep");
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  MultiLevelSolution *ml_sol			      = ml_prob._ml_sol;
  Solution*	 sol  	                      = ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  pdeSys	              = mlPdeSys->_LinSolver[level];
  const char* pdename                                 = mlPdeSys->name().c_str();

  Mesh*		 msh    	   = ml_prob._ml_msh->GetLevel(level);
  elem*		 myel		   = msh->el;
  SparseMatrix*	 KK	 	   = pdeSys->_KK;
  NumericVector* RES 		   = pdeSys->_RES;

  //data
  const unsigned dim = msh->GetDimension();
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  unsigned nel = msh->GetNumberOfElements();
  unsigned igrid = msh->GetLevel();
  unsigned iproc = msh->processor_id();
  
  //time dep data
  double dt = mlPdeSys->GetIntervalTime();
  double theta = 0.5;

//************** geometry (at dofs and quadrature points) *************************************  
  vector < vector < double > > coords_at_dofs(dim);
  unsigned coords_fe_type = BIQUADR_FE; // get the finite element type for "x", it is always 2 (LAGRANGE BIQUADRATIC)
  for (unsigned i = 0; i < coords_at_dofs.size(); i++)    coords_at_dofs[i].reserve(max_size);

  vector < double > coord_at_qp(dim);
  
  
 //************* shape functions (at dofs and quadrature points) **************************************  
  const int solType_max = BIQUADR_FE;  //biquadratic

  double weight_qp; // gauss point weight
  
  vector < vector < double > > phi_fe_qp(NFE_FAMS);
  vector < vector < double > > phi_x_fe_qp(NFE_FAMS);
  vector < vector < double > > phi_xx_fe_qp(NFE_FAMS);
 
  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_fe_qp[fe].reserve(max_size);
      phi_x_fe_qp[fe].reserve(max_size*dim);
     phi_xx_fe_qp[fe].reserve(max_size*(3*(dim-1)));
   }

  
 //***************************************************  
 //********* WHOLE SET OF VARIABLES ****************** 
 //***************************************************  
  const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();

  vector < std::string > Solname(n_unknowns);     Solname[0] = "u";
  vector < unsigned > SolPdeIndex(n_unknowns);
  vector < unsigned > SolIndex(n_unknowns);
  vector < unsigned int > SolFEType(n_unknowns);    //FEtype of each MultilevelSolution       
  vector < unsigned int > Sol_n_el_dofs(n_unknowns); //number of element dofs

  std::fill(Sol_n_el_dofs.begin(), Sol_n_el_dofs.end(), 0);

  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolPdeIndex[ivar] = mlPdeSys->GetSolPdeIndex(  Solname[ivar].c_str() );
       SolIndex[ivar] = ml_sol->GetIndex         (  Solname[ivar].c_str() );
      SolFEType[ivar] = ml_sol->GetSolutionType  ( SolIndex[ivar]);
  }
  
  //------------ quantities (at quadrature points) ---------------------
            vector<double>        sol_qp(n_unknowns);
            vector<double>    sol_old_qp(n_unknowns);
    vector< vector<double> > sol_grad_qp(n_unknowns);
    vector< vector<double> > sol_old_grad_qp(n_unknowns);
    
      std::fill(sol_qp.begin(), sol_qp.end(), 0.);
      std::fill(sol_old_qp.begin(), sol_old_qp.end(), 0.);
    for (unsigned  k = 0; k < n_unknowns; k++) {
        sol_grad_qp[k].resize(dim);
        std::fill(sol_grad_qp[k].begin(), sol_grad_qp[k].end(), 0.);
        sol_old_grad_qp[k].resize(dim);
        std::fill(sol_old_grad_qp[k].begin(), sol_old_grad_qp[k].end(), 0.);
    }

  //----------- quantities (at dof objects) ------------------------------
  vector < vector < double > >     sol_eldofs(n_unknowns);
  vector < vector < double > >     sol_old_eldofs(n_unknowns);
  for(int k=0; k<n_unknowns; k++) { sol_eldofs[k].reserve(max_size);
                                sol_old_eldofs[k].reserve(max_size);
   }
  
  //******** linear system *******************************************  
  vector < vector < int > > L2G_dofmap(n_unknowns);     for(int i = 0; i < n_unknowns; i++) { L2G_dofmap[i].reserve(max_size); }
            vector< int >   L2G_dofmap_AllVars; L2G_dofmap_AllVars.reserve( n_unknowns*max_size );
          
  vector< vector< double > > F(n_unknowns);
  vector< vector< vector< double > > > B(n_unknowns);

   for(int i = 0; i < n_unknowns; i++) F[i].reserve(max_size);

   for(int i = 0; i < n_unknowns; i++) {
    B[i].resize(n_unknowns);
    for(int j = 0; j < n_unknowns; j++) {
      B[i][j].reserve(max_size*max_size);
    }
  }

  
  // Set to zero all the entries of the matrix
  RES->zero();
   KK->zero();

  // *** element loop ***

  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc+1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);    // element geometry type
    
     //******************** GEOMETRY ********************* 
    unsigned nDofx = msh->GetElementDofNumber(iel, coords_fe_type);    // number of coordinate element dofs
    for (int i = 0; i < dim; i++)  coords_at_dofs[i].resize(nDofx);
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, coords_fe_type);  // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        coords_at_dofs[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }
 //***************************************************  
    
   //all vars###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned  ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
	   Sol_n_el_dofs[k] = ndofs_unk;
          sol_eldofs[k].resize(ndofs_unk);
      sol_old_eldofs[k].resize(ndofs_unk);
          L2G_dofmap[k].resize(ndofs_unk); 
    for (unsigned i = 0; i < ndofs_unk; i++) {
           unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);                        // global to global mapping between solution node and solution dof // via local to global solution node
           sol_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);                            // global extraction and local storage for the solution
       sol_old_eldofs[k][i] = (*sol->_SolOld[SolIndex[k]])(solDof);                            // global extraction and local storage for the solution
           L2G_dofmap[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    
    
    unsigned nDof_AllVars = 0;
    for (unsigned  k = 0; k < n_unknowns; k++) { nDof_AllVars += Sol_n_el_dofs[k]; }
 // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    int nDof_max    =  0;   
      for (unsigned  k = 0; k < n_unknowns; k++)     {
          if(Sol_n_el_dofs[k] > nDof_max)    nDof_max = Sol_n_el_dofs[k];
      }
    

     for(int k = 0; k < n_unknowns; k++) {
     L2G_dofmap[k].resize(Sol_n_el_dofs[k]);

      F[SolPdeIndex[k]].resize(Sol_n_el_dofs[k]);
      memset(&F[SolPdeIndex[k]][0], 0., Sol_n_el_dofs[k]*sizeof(double));


      B[SolPdeIndex[k]][SolPdeIndex[k]].resize(Sol_n_el_dofs[k] * Sol_n_el_dofs[k]);
      memset(&B[SolPdeIndex[k]][SolPdeIndex[k]][0], 0., Sol_n_el_dofs[k] * Sol_n_el_dofs[k] * sizeof(double));
    }

  //all vars###################################################################    


      // *** quadrature loop ***
      for(unsigned ig = 0; ig < ml_prob._ml_msh->_finiteElement[ielGeom][solType_max]->GetGaussPointNumber(); ig++) {
          
      // *** get gauss point weight, test function and test function partial derivatives ***
      for(int fe=0; fe < NFE_FAMS; fe++) {
         msh->_finiteElement[ielGeom][fe]->Jacobian(coords_at_dofs,ig,weight_qp,phi_fe_qp[fe],phi_x_fe_qp[fe],phi_xx_fe_qp[fe]);
      }
   //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
         msh->_finiteElement[ielGeom][coords_fe_type]->Jacobian(coords_at_dofs,ig,weight_qp,phi_fe_qp[coords_fe_type],phi_x_fe_qp[coords_fe_type],phi_xx_fe_qp[coords_fe_type]);
         

 //========= fill gauss value quantities ==================   
   std::fill(sol_qp.begin(), sol_qp.end(), 0.);
   std::fill(sol_old_qp.begin(), sol_old_qp.end(), 0.);
   for (unsigned  k = 0; k < n_unknowns; k++) { std::fill(sol_grad_qp[k].begin(), sol_grad_qp[k].end(), 0.); 
                                                std::fill(sol_old_grad_qp[k].begin(), sol_old_grad_qp[k].end(), 0.);
                                            }
    
    for (unsigned  k = 0; k < n_unknowns; k++) {
	for (unsigned i = 0; i < Sol_n_el_dofs[k]; i++) {
	                                                         sol_qp[k]    +=     sol_eldofs[k][i] *   phi_fe_qp[SolFEType[k]][i];
	                                                     sol_old_qp[k]    += sol_old_eldofs[k][i] *   phi_fe_qp[SolFEType[k]][i];
                   for (unsigned d = 0; d < dim; d++) {      sol_grad_qp[k][d] +=     sol_eldofs[k][i] * phi_x_fe_qp[SolFEType[k]][i * dim + d]; 
                                                         sol_old_grad_qp[k][d] += sol_old_eldofs[k][i] * phi_x_fe_qp[SolFEType[k]][i * dim + d]; 
                                     }
       }        
    }
 //========= fill gauss value quantities ==================
         

	// *** phi_i loop ***
	for(unsigned i=0; i<nDof_max; i++){

	  //BEGIN RESIDUALS A block ===========================
	    double Lap_rhs = 0.;
	    double Lap_old_rhs = 0.;
	    for(unsigned d = 0;  d < dim;  d++) {
	      Lap_rhs 	  += phi_x_fe_qp[SolFEType[0]][i * dim + d] *     sol_grad_qp[0][d];
	      Lap_old_rhs += phi_x_fe_qp[SolFEType[0]][i * dim + d] * sol_old_grad_qp[0][d];
	    }
	    F[SolPdeIndex[0]][i] += weight_qp * ( 
                    -       theta  * dt * Lap_rhs                           // Laplacian
					- (1. - theta) * dt * Lap_old_rhs                  // Laplacian
					- (sol_qp[0] - sol_old_qp[0]) * phi_fe_qp[ SolFEType[0] ][i]        // acceleration
            );
	  //END RESIDUALS A block ===========================

	  {
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nDof_max; j++) {
            
          double Mass = phi_fe_qp[ SolFEType[0] ][i]*phi_fe_qp[ SolFEType[0] ][j];

          double Lap_mat = 0.;
          for(unsigned d = 0; d < dim; d++) {
	        Lap_mat  += phi_x_fe_qp[SolFEType[0]][i * dim + d] * phi_x_fe_qp[SolFEType[0]][j * dim + d];
	      }

		B[SolPdeIndex[0]][SolPdeIndex[0]][i * Sol_n_el_dofs[0] + j] += weight_qp * (Lap_mat + Mass);
	      
  	    } //end phij loop

	  } // endif assemble_matrix
	} //end phii loop


      }  // end gauss point loop


//--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
      RES->add_vector_blocked(F[SolPdeIndex[ivar]],L2G_dofmap[ivar]);
      KK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[ivar]],L2G_dofmap[ivar],L2G_dofmap[ivar]);
    }
    //--------------------------------------------------------------------------------------------------------
  } //end list of elements loop for each subdomain


  KK->close();
  RES->close();
  // ***************** END ASSEMBLY *******************
}


