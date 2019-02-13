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
                 value = 0.;
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

  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,0.001,1.,"Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  // add fluid material
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;


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
  system.SetIntervalTime(0.1);
  const unsigned int n_timesteps = 20;
  const unsigned int write_interval = 1;

    // ======= Time Loop ========================
  for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {

    // ======= Solve ========================
    std::cout << std::endl;
    std::cout << " *********** Timedep ************  " << std::endl;
    ml_prob.get_system("Timedep").MLsolve();

    //update Solution
    ml_prob.get_system<TransientNonlinearImplicitSystem>("Timedep").CopySolutionToOldSolution();

    // ======= Final Print ========================
    if ( !(time_step%write_interval) ) {

    // ======= Final Print ========================
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
  ml_sol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/, "biquadratic", variablesToBePrinted,time_step);    // print solutions

    }

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
  Solution*	 mysolution  	                      = ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  mylsyspde	              = mlPdeSys->_LinSolver[level];
  const char* pdename                                 = mlPdeSys->name().c_str();

  Mesh*		 msh    	   = ml_prob._ml_msh->GetLevel(level);
  elem*		 myel		   = msh->el;
  SparseMatrix*	 KK	 	   = mylsyspde->_KK;
  NumericVector* RES 		   = mylsyspde->_RES;

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
  

  // declare
  vector< int > metis_node2;
  vector< vector< int > > KK_dof(n_unknowns);
//   vector<double> phi2;
  vector<double> gradphi2;
  vector<double> nablaphi2;
  double Weight2;
  double normal[3];
  vector< vector< double > > F(n_unknowns);
  vector< vector< vector< double > > > B(n_unknowns);

  metis_node2.reserve(max_size);
  gradphi2.reserve(max_size*dim);
  nablaphi2.reserve(max_size*(3*(dim-1)));
  for(int i=0;i<dim;i++) {
    KK_dof[i].reserve(max_size);
  }

  for(int i = 0; i < n_unknowns; i++) F[i].reserve(max_size);


  for(int i = 0; i < n_unknowns; i++) {
    B[i].resize(n_unknowns);
    for(int j = 0; j < n_unknowns; j++) {
      B[i][j].reserve(max_size*max_size);
    }
  }


  vector < double > Sol(1);
  vector < vector < double > > gradSol(dim);
  for(int i = 0; i < dim; i++) {
    gradSol[i].resize(dim);
  }

  vector < double > SolOld(n_unknowns);
  vector < vector < double > > gradSolOld(dim);
  for(int i = 0; i < dim; i++) {
    gradSolOld[i].resize(dim);
  }
  //vector < double > AccSol(dim);

  // Set to zeto all the entries of the matrix
  KK->zero();
  RES->zero();

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
    
    
    
    
    unsigned nve2 = msh->GetElementDofNumber(iel,SolFEType[0]);

    //set to zero all the entries of the FE matrices
    metis_node2.resize(nve2);
    gradphi2.resize(nve2*dim);
    nablaphi2.resize(nve2*(3*(dim-1)));


    
    for(int ivar=0; ivar<n_unknowns; ivar++) {
      KK_dof[ivar].resize(nve2);

      F[SolPdeIndex[ivar]].resize(nve2);
      memset(&F[SolPdeIndex[ivar]][0],0,nve2*sizeof(double));


      B[SolPdeIndex[ivar]][SolPdeIndex[ivar]].resize(nve2*nve2);
      memset(&B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][0],0,nve2*nve2*sizeof(double));

    }




    for( unsigned i=0;i<nve2;i++) {
      unsigned inode_metis=msh->GetSolutionDof(i, iel, 2);
      metis_node2[i]=inode_metis;
      for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
           KK_dof[ivar][i]=mylsyspde->GetSystemDof(SolIndex[ivar],SolPdeIndex[ivar],i, iel);
        }
    }



      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[ielGeom][solType_max]->GetGaussPointNumber(); ig++) {
          
      // *** get gauss point weight, test function and test function partial derivatives ***
      for(int fe=0; fe < NFE_FAMS; fe++) {
         msh->_finiteElement[ielGeom][fe]->Jacobian(coords_at_dofs,ig,weight_qp,phi_fe_qp[fe],phi_x_fe_qp[fe],phi_xx_fe_qp[fe]);
      }
   //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
         msh->_finiteElement[ielGeom][coords_fe_type]->Jacobian(coords_at_dofs,ig,weight_qp,phi_fe_qp[coords_fe_type],phi_x_fe_qp[coords_fe_type],phi_xx_fe_qp[coords_fe_type]);
          
	//velocity variable
	for(unsigned ivar=0; ivar<n_unknowns; ivar++) {
	  Sol[ivar]=0.;
	  SolOld[ivar]=0.;
	  //AccSol[ivar]=0.;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++){
	    gradSol[ivar][ivar2]=0;
	    gradSolOld[ivar][ivar2]=0;
	  }
	  unsigned SolIndex=ml_sol->GetIndex(&Solname[ivar][0]);
	  unsigned SolType =ml_sol->GetSolutionType(&Solname[ivar][0]);
	  for(unsigned i=0; i<nve2; i++) {
	    double soli    = (*mysolution->_Sol[SolIndex])(metis_node2[i]);
	    double sololdi = (*mysolution->_SolOld[SolIndex])(metis_node2[i]);
	    Sol[ivar]    += phi_fe_qp[ SolFEType[0] ][i]*soli;
	    SolOld[ivar] += phi_fe_qp[ SolFEType[0] ][i]*sololdi;
	    for(unsigned ivar2=0; ivar2<dim; ivar2++){
	      gradSol[ivar][ivar2]    += gradphi2[i*dim+ivar2]*soli;
	      gradSolOld[ivar][ivar2] += gradphi2[i*dim+ivar2]*sololdi;
	    }
	  }
	}


	// *** phi_i loop ***
	for(unsigned i=0; i<nve2; i++){

	  //BEGIN RESIDUALS A block ===========================
	  for(unsigned ivar=0; ivar<n_unknowns; ivar++) {
	    double Lap_rhs=0;
	    double Lap_old_rhs=0;
	    for(unsigned ivar2=0; ivar2< dim; ivar2++) {
	      Lap_rhs 	  += gradphi2[i*dim+ivar2]*gradSol[ivar][ivar2];
	      Lap_old_rhs += gradphi2[i*dim+ivar2]*gradSolOld[ivar][ivar2];
	    }
	    F[SolPdeIndex[ivar]][i] += ( 
                    -theta*dt*Lap_rhs                           // Laplacian
					-(1.-theta)*dt*Lap_old_rhs                  // Laplacian
					+dt*Sol[dim]*gradphi2[i*dim+ivar]                  // pressure
					-(Sol[ivar] - SolOld[ivar])*phi_fe_qp[ SolFEType[0] ][i]        // acceleration
	    )*Weight2;
	  }
	  //END RESIDUALS A block ===========================

	  {
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      double Lap = 0.;
          double Mass = phi_fe_qp[ SolFEType[0] ][i]*phi_fe_qp[ SolFEType[0] ][j]*Weight2;
          
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		// Laplacian
		Lap  += gradphi2[i*dim+ivar]*gradphi2[j*dim+ivar]*Weight2;
		// advection term I
	      }

	      for(unsigned ivar=0; ivar< n_unknowns; ivar++) {
		B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][i*nve2+j] += theta*dt*Lap + Mass;
	      }
  	    } //end phij loop

	  } // endif assemble_matrix
	} //end phii loop


      }  // end gauss point loop


//--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
      RES->add_vector_blocked(F[SolPdeIndex[ivar]],KK_dof[ivar]);
      KK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[ivar]],KK_dof[ivar],KK_dof[ivar]);
    }
    //--------------------------------------------------------------------------------------------------------
  } //end list of elements loop for each subdomain


  KK->close();
  RES->close();
  // ***************** END ASSEMBLY *******************
}


