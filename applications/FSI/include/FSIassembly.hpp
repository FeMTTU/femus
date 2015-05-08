#ifndef __femus_include_FSIassembly_hpp__
#define __femus_include_FSIassembly_hpp__

#include "TransientSystem.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"

namespace femus {

void AssembleMatrixResFSI(MultiLevelProblem &ml_prob) {
    
  clock_t AssemblyTime=0;
  clock_t start_time, end_time;
  PetscErrorCode ierr;
  
  //pointers and references
  
  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem>("Fluid-Structure-Interaction");
  const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();
  const unsigned gridn = my_nnlin_impl_sys.GetLevelMax();
  bool assemble_matrix = my_nnlin_impl_sys.GetAssembleMatrix(); 
    
  MultiLevelSolution*	 ml_sol	                      = ml_prob._ml_sol;
  Solution*	 mysolution  	                      = ml_sol->GetSolutionLevel(level);
  LinearEquationSolver*  mylsyspde	              = my_nnlin_impl_sys._LinSolver[level];   
  Mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(level);
  elem		*myel		=  mymsh->el;
  SparseMatrix	*myKK		=  mylsyspde->_KK;
  NumericVector *myRES		=  mylsyspde->_RES;
  vector <int>	&myKKIndex	=  mylsyspde->KKIndex;
  
  
  const unsigned dim = mymsh->GetDimension();
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  
  // local objects
  vector<double> SolVAR(3*dim+1);
  vector<double> SolOldVAR(3*dim+1);
  vector<vector<double> > GradSolVAR(2*dim);
  vector<vector<double> > GradSolOldVAR(2*dim);
  for(int i=0;i<2*dim;i++){
    GradSolVAR[i].resize(dim);
    GradSolOldVAR[i].resize(dim);
  }
  vector<vector<double> > GradSolhatVAR(dim);
  vector<vector<double> > GradSolOldhatVAR(dim);
  for(int i=0;i<dim;i++){
    GradSolhatVAR[i].resize(dim);
    GradSolOldhatVAR[i].resize(dim);
  }
    
  vector <int> metis_node1;
  vector <int> metis_node2;
  vector <bool> solidmark;
  
  vector <double > phi;
  vector <double > phi_hat;
  vector <double > phi_old;
  
  vector <double> gradphi;
  vector <double> gradphi_hat;
  vector <double> gradphi_old;
  
  vector <double > nablaphi;
  
  metis_node1.reserve(max_size);
  metis_node2.reserve(max_size);
  solidmark.reserve(max_size);
  phi.reserve(max_size);
  phi_hat.reserve(max_size);
  phi_old.reserve(max_size);
  gradphi.reserve(max_size*dim);
  gradphi_hat.reserve(max_size*dim);
  gradphi_old.reserve(max_size*dim);
  
  nablaphi.reserve(max_size*(3*(dim-1)+!(dim-1)));
  
  const double *phi1;
    
  double Weight=0.;
  double Weight_nojac=0.;
  double Weight_hat=0.;
  double Weight_old=0.;
  
  vector <vector < double> > vx(dim);
  vector <vector < double> > vx_hat(dim);
  vector <vector < double> > vx_old(dim);
  
  for(int i=0;i<dim;i++){
    vx[i].reserve(max_size);
    vx_hat[i].reserve(max_size);
    vx_old[i].reserve(max_size);
  }
   
  vector< vector< double > > Rhs(2*dim+1);
  vector< vector< vector< double > > > B(2*dim+1); 
  for(int i=0;i<2*dim+1;i++){
    B[i].resize(2*dim+1);
  }
  vector< vector< int > > dofsVAR(2*dim+1); 
  
  // algorithm parameters
  double eps_pen 	= 1.;  // previous 1
  bool   newton		= 0;
  
  // ------------------------------------------------------------------------
  // Physical parameters
  double rhof	 	= ml_prob.parameters.get<Fluid>("Fluid").get_density();            // ml_prob._fluid->get_density();
  double rhos 		= ml_prob.parameters.get<Solid>("Solid").get_density();            // ml_prob._solid->get_density();
  double mu_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_shear_modulus(); // ml_prob._solid->get_lame_shear_modulus();
  double lambda_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_lambda();        // ml_prob._solid->get_lame_lambda();
  double mus		= mu_lame/rhof;
  double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();   // ml_prob._fluid->get_IReynolds_number();
  double lambda		= lambda_lame / rhof;
  double betafsi	= rhos / rhof;
  double betans		= 1.;
  int    solid_model	= ml_prob.parameters.get<Solid>("Solid").get_physical_model();     // ml_prob._solid->get_physical_model();

  //physical quantity
  double Jnp1_hat;
  double Jn_hat;
  double I_bleft;
  double I_e;
  double Cauchy[3][3];
  double Cauchy_old[3][3]; 
  double tg_stiff_matrix[3][3];
  //initialization C tensor: Saint-Venaint Kirchoff model : solid_model==1;
  const double Id2th[3][3]= {{ 1.,0.,0.}, { 0.,1.,0.}, { 0.,0.,1.}};
  double C_mat[3][3][3][3];
  for (int I=0; I<3; ++I) {
    for (int J=0; J<3; ++J) {
      for (int K=0; K<3; ++K) {
        for (int L=0; L<3; ++L) {
          C_mat[I][J][K][L] = 2.*mus*Id2th[I][K]*Id2th[J][L];
        }
      }
    }
  }

  // ale map
  double _lambda_map=0.;
  double _mu_ale[3] = {1.,1.,1.};
 
  // gravity
  double _gravity[3]={0.,0.,0.};

  // -----------------------------------------------------------------

  // time discretization algorithm paramters
  double dt =  my_nnlin_impl_sys.GetIntervalTime(); 
  
  const double gamma = 0.5;
  const double gammaratio = (1.-gamma)/gamma;
  double _theta_ns=1.0;
  
  // space discretization parameters
  unsigned order_ind2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U")); 
  unsigned order_ind1 = ml_sol->GetSolutionType(ml_sol->GetIndex("P"));  

  // mesh and procs
  unsigned nel    = mymsh->GetNumberOfElements();
  unsigned igrid  = mymsh->GetLevel();
  unsigned iproc  = mymsh->processor_id();

  //----------------------------------------------------------------------------------
  //variable-name handling
  const char varname[10][3] = {"DX","DY","DZ","U","V","W","P","AX","AY","AZ"};
  //const char coordname[3][2] = {"X","Y","Z"};
  vector <unsigned> indexVAR(2*dim+1);
  //vector <unsigned> indCOORD(dim);
  vector <unsigned> indVAR(3*dim+1);  
  vector <unsigned> SolType(3*dim+1);  /// @todo is this really 3*d + 1 or 2*d+1 would be enough?
  
  for(unsigned ivar=0; ivar<dim; ivar++) {
    //indCOORD[ivar]=ml_prob.GetIndex(&coordname[ivar][0]);
    indVAR[ivar]=ml_sol->GetIndex(&varname[ivar][0]);
    indVAR[ivar+dim]=ml_sol->GetIndex(&varname[ivar+3][0]);
    indVAR[ivar+2*dim+1]=ml_sol->GetIndex(&varname[ivar+7][0]);
    
    SolType[ivar]=ml_sol->GetSolutionType(&varname[ivar][0]);
    SolType[ivar+dim]=ml_sol->GetSolutionType(&varname[ivar+3][0]);
    SolType[ivar+2*dim+1]=ml_sol->GetSolutionType(&varname[ivar+7][0]);
    indexVAR[ivar]=my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
    indexVAR[ivar+dim]=my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar+3][0]);
  }
  indexVAR[2*dim]=my_nnlin_impl_sys.GetSolPdeIndex(&varname[6][0]);
  indVAR[2*dim]=ml_sol->GetIndex(&varname[6][0]);
  SolType[2*dim]=ml_sol->GetSolutionType(&varname[6][0]);
  //----------------------------------------------------------------------------------
  
  start_time=clock();
    
  myKK->zero();
  
  /// *** element loop ***
  for(int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel        = mymsh->IS_Mts2Gmt_elem[iel]; 
    short unsigned kelt = myel->GetElementType(kel);
    unsigned nve        = myel->GetElementDofNumber(kel,order_ind2);
    unsigned nve1       = myel->GetElementDofNumber(kel,order_ind1);
    int flag_mat        = myel->GetElementMaterial(kel);

    //*******************************************************************************************************
    
    //initialization of everything is in common fluid and solid
    
    //Rhs
    for(int i=0; i<2*dim; i++) {
      dofsVAR[i].resize(nve);	
      Rhs[indexVAR[i]].resize(nve);
      memset(&Rhs[indexVAR[i]][0],0,nve*sizeof(double));
    }
    dofsVAR[2*dim].resize(nve1);
    Rhs[indexVAR[2*dim]].resize(nve1);
    memset(&Rhs[indexVAR[2*dim]][0],0,nve1*sizeof(double));
      
    //Kinematic relation (solid) and ALE Map (fluid)
    for(int i=0; i<dim; i++) {
      B[indexVAR[i]][indexVAR[i]].resize(nve*nve);
      memset(&B[indexVAR[i]][indexVAR[i]][0],0,nve*nve*sizeof(double));
    }

    //Stiffness Matrix (solid)
    for(int i=0; i<dim; i++) {
      for(int j=0; j<dim; j++) {
	B[indexVAR[dim+i]][indexVAR[j]].resize(nve*nve);
	memset(&B[indexVAR[dim+i]][indexVAR[j]][0],0,nve*nve*sizeof(double));
      }
    }
      
    //Mass Matrix (solid and fluid) and Diffusion Matrix (fluid only) 
    for(int i=0; i<dim; i++) {
      B[indexVAR[dim+i]][indexVAR[dim+i]].resize(nve*nve);
      memset(&B[indexVAR[dim+i]][indexVAR[dim+i]][0],0,nve*nve*sizeof(double));
    }

    //Pressure gradient (fluid and solid)
    for(int i=0; i<dim; i++) {
      B[indexVAR[dim+i]][indexVAR[2*dim]].resize(nve*nve1);
      memset(&B[indexVAR[dim+i]][indexVAR[2*dim]][0],0,nve*nve1*sizeof(double));
    }

    // Pressure Mass Matrix
    B[indexVAR[2*dim]][indexVAR[2*dim]].resize(nve1*nve1);
    memset(&B[indexVAR[2*dim]][indexVAR[2*dim]][0],0,nve1*nve1*sizeof(double));
    
    if (flag_mat==2) { //initialization for fluid only
      // Fluid Continuity Matrix: divergence of the velocity
      for(int i=0; i<dim; i++) {
	B[indexVAR[2*dim]][indexVAR[dim+i]].resize(nve1*nve);
	memset(&B[indexVAR[2*dim]][indexVAR[dim+i]][0],0,nve1*nve*sizeof(double));
      }
    }
    else{ // initialization for solid only
      // Kinematic relation
      for(int i=0; i<dim; i++) {
	B[indexVAR[i]][indexVAR[dim+i]].resize(nve*nve);
	memset(&B[indexVAR[i]][indexVAR[dim+i]][0],0,nve*nve*sizeof(double));
      }
      // Solid Continuity Matrix: divergence of the displacemnet
      for(int i=0; i<dim; i++) {
	B[indexVAR[2*dim]][indexVAR[i]].resize(nve1*nve);
	memset(&B[indexVAR[2*dim]][indexVAR[i]][0],0,nve1*nve*sizeof(double));
      }
    }
    
    // ----------------------------------------------------------------------------------------
    // coordinates, displacement, velocity dofs
    
    metis_node2.resize(nve);
    metis_node1.resize(nve1);
    solidmark.resize(nve);
    phi.resize(nve);
    phi_hat.resize(nve);
    phi_old.resize(nve);
    gradphi.resize(nve*dim);
    gradphi_hat.resize(nve*dim);
    gradphi_old.resize(nve*dim);
        
    nablaphi.resize(nve*(3*(dim-1)+!(dim-1)));
    
    for(int i=0;i<dim;i++){
      vx[i].resize(nve);
      vx_old[i].resize(nve);
      vx_hat[i].resize(nve);
    }
    
    for (unsigned i=0;i<nve;i++) {
      // gambit nodes
      unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;
      // dof metis
      unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
      metis_node2[i]=inode_Metis;
      
      //unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
      // flag to know if the node "inode" lays on the fluid-solid interface
      solidmark[i]=myel->GetNodeRegion(inode); // to check
      for(int j=0; j<dim; j++) {
	//Updated coordinates (Moving frame)
        vx[j][i]= (*mymsh->_coordinate->_Sol[j])(inode_Metis) + (*mysolution->_Sol[indVAR[j]])(inode_Metis);
	//Old coordinates (Moving frame)
        vx_old[j][i]= (*mymsh->_coordinate->_Sol[j])(inode_Metis) + (*mysolution->_SolOld[indVAR[j]])(inode_Metis);
	//Fixed coordinates (Reference frame)
	vx_hat[j][i]= (*mymsh->_coordinate->_Sol[j])(inode_Metis);  
	// displacement dofs
	dofsVAR[j][i]= mylsyspde->GetKKDof(indVAR[j],indexVAR[j],inode); 
	// velocity dofs
	dofsVAR[j+dim][i]= mylsyspde->GetKKDof(indVAR[j+dim],indexVAR[j+dim],inode);   
      }
    }

    // pressure dofs
    for (unsigned i=0;i<nve1;i++) {
      unsigned inode=(order_ind1<3)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
      metis_node1[i]=mymsh->GetMetisDof(inode,SolType[2*dim]);
      dofsVAR[2*dim][i]=mylsyspde->GetKKDof(indVAR[2*dim],indexVAR[2*dim],inode);
    }
    // ----------------------------------------------------------------------------------------
       
    if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
      /// *** Gauss point loop ***
      for (unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->GetGaussPointNumber(); ig++) {

	// *** get Jacobian and test function and test function derivatives in the moving frame***
	ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->Jacobian(vx,ig,Weight,phi,gradphi,nablaphi);
	ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->Jacobian(vx_old,ig,Weight_old,phi_old,gradphi_old,nablaphi);
	ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->Jacobian(vx_hat,ig,Weight_hat,phi_hat,gradphi_hat,nablaphi);
	phi1=ml_prob._ml_msh->_finiteElement[kelt][order_ind1]->GetPhi(ig);
	if (flag_mat==2) Weight_nojac = ml_prob._ml_msh->_finiteElement[kelt][order_ind2]->GetGaussWeight(ig);

	// ---------------------------------------------------------------------------
	// displacement and velocity
	for(int i=0; i<2*dim; i++){
	  SolVAR[i]=0.;
	  SolOldVAR[i]=0.;
	  for(int j=0; j<dim; j++) {
	    GradSolVAR[i][j]=0.;
	    GradSolOldVAR[i][j]=0.;
	    if(i<dim){
	      GradSolhatVAR[i][j]=0.;
	      GradSolOldhatVAR[i][j]=0.;
	    }
	  }
	    
	  for (unsigned inode=0; inode<nve; inode++) {
	    unsigned sol_dof = metis_node2[inode];
	      
	    double soli = (*mysolution->_Sol[indVAR[i]])(sol_dof);
	    SolVAR[i]+=phi[inode]*soli;
	      
	    double soli_old = (*mysolution->_SolOld[indVAR[i]])(sol_dof);
	    SolOldVAR[i]+=phi[inode]*soli_old;
	      
	    for(int j=0; j<dim; j++) {
	      GradSolVAR[i][j]+=gradphi[inode*dim+j]*soli;
	      GradSolOldVAR[i][j]+=gradphi[inode*dim+j]*soli_old;
	      if(i<dim){ 
		GradSolhatVAR[i][j]   +=gradphi_hat[inode*dim+j]*soli;
		GradSolOldhatVAR[i][j]+=gradphi_hat[inode*dim+j]*soli_old;
	      }
	    }	      
	  }
	}  
		  
	// pressure
	SolVAR[2*dim]=0.;
	for (unsigned inode=0; inode<nve1; inode++) {
	  double soli = (*mysolution->_Sol[indVAR[2*dim]])(metis_node1[inode]);
	  SolVAR[2*dim]+=phi1[inode]*soli;
	}
  
	// acceleration (solid only)  
	if(flag_mat!=2){
	  for(int i=2*dim+1; i<3*dim+1; i++){
	    SolVAR[i]=0.;
	    for (unsigned inode=0; inode<nve; inode++) {
	      double soli = (*mysolution->_Sol[indVAR[i]])(metis_node2[inode]);
	      SolVAR[i]+=phi[i]*soli;  
	    }
	  }
	}
  
	// ---------------------------------------------------------------------------
	 	  
	//BEGIN FLUID ASSEMBLY ============
	  
	if(flag_mat==2){
	  
          //divergence of the velocity
	  double div_vel=0.;
	  double div_w=0.;
	  for(int i=0; i<dim; i++) {
	    div_vel+=GradSolVAR[dim+i][i];
	    div_w+=(GradSolVAR[i][i]-GradSolOldVAR[i][i])*(1./dt);
	  }
	  
          { // Laplace operator + adection operator + Mass operator

            const double *gradfi=&gradphi[0];
            const double *fi=&phi[0];

            // *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {

              //BEGIN RESIDUALS A block ===========================
 	      double LapmapVAR[3] = {0., 0., 0.};
 	      for(int idim=0; idim<dim; idim++) {
 		for(int idim2=0; idim2<dim; idim2++) {
 	          LapmapVAR[idim] += dt*( _mu_ale[idim2]*gradphi_hat[i*dim+idim2]*GradSolhatVAR[idim][idim2] );
 		}
 	      }
      
	      // Residual ALE equations
 	      for(int idim=0; idim<dim; idim++) {
 	        Rhs[indexVAR[idim]][i]+=(!solidmark[i])*(-LapmapVAR[idim]*Weight_nojac);
 	      }
	      
	      //-------------------------------------------------------------------------------
	      
	      double LapvelVAR[3]={0.,0.,0.};
	      double AdvaleVAR[3]={0.,0.,0.};
	      for(int idim=0.; idim<dim; idim++) {
		for(int idim2=0.; idim2<dim; idim2++) {
		  LapvelVAR[idim]+=gradphi[i*dim+idim2]*GradSolVAR[dim+idim][idim2];
		  AdvaleVAR[idim]+=((SolVAR[dim+idim2]*dt - (SolVAR[idim2]-SolOldVAR[idim2]))*GradSolVAR[dim+idim][idim2])*phi[i];;
		}
	      }
	      // Residual Momentum equations 
	      for(int idim=0; idim<dim; idim++) {
		Rhs[indexVAR[dim+idim]][i]+= (
					      -phi[i]*betans*SolVAR[dim+idim]*Weight
					      +phi[i]*betans*SolOldVAR[dim+idim]*Weight_old
					      -AdvaleVAR[idim]*0.5*(Weight+Weight_old)
					      -0.5*dt*div_vel*SolVAR[dim+idim]*phi[i]*0.5*(Weight+Weight_old)
					      +dt*div_w*SolVAR[dim+idim]*phi[i]*0.5*(Weight+Weight_old)
					      -dt*IRe*LapvelVAR[idim]*Weight
					      +dt*SolVAR[2*dim]*gradphi[i*dim+idim]*Weight
					      );
	      }
	      //END RESIDUALS A block ===========================
	      
              const double *gradfj=&gradphi[0];
              const double *fj=&phi[0];
              //  *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim,fj++) {

                //Laplacian
                double Lap=0.;
		for(int idim=0; idim<dim; idim++) {
		  Lap+=(*(gradfi+idim))*(*(gradfj+idim));
		}
                double LapXweight = Lap*Weight;

                //advection term I
		double Adv1=0.;
		for(int idim=0; idim<dim; idim++) {
		  Adv1+= ((SolVAR[dim+idim]*dt - (SolVAR[idim]-SolOldVAR[idim]))*(*(gradfj+idim))*(*(fi)))*0.5*(Weight+Weight_old);
		}
				
		double div_stab = 0.5*div_vel*((*fi))*((*fj))*0.5*(Weight+Weight_old);
		
		double div_ale = div_w*((*fi))*((*fj))*0.5*(Weight+Weight_old);

                double Mass = ((*fi))*((*fj))*Weight;

		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[dim+idim]][i*nve+j] += dt*IRe*LapXweight + Adv1 + dt*div_stab - dt*div_ale + betans*Mass;
		}
		
//                 for(int idim=0; idim<dim; idim++) {
// 		  for(int idim2=0; idim2<dim; idim2++) {
// 		    B[indexVAR[dim+idim]][indexVAR[idim2]][i*nve+j] += betans*SolVAR[dim+idim]*((*fi))*(*(gradfi+idim2))*Weight;
// 		  }
// 		}

		double Lap_ale=0.;
		for(int idim=0; idim<dim; idim++) {
		  Lap_ale+=_mu_ale[idim]*(*(gradfi+idim))*(*(gradfj+idim));
		}
		  
		// Laplacian ALE map
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[0+idim]][indexVAR[idim]][i*nve+j] += (!solidmark[i])*dt*Lap_ale*Weight_nojac;  
		}

              } // end phi_j loop
            } // end phi loop
          } // end A + Bt

          { //Gradient of Pressure operator

            const double *gradfi=&gradphi[0];
            const double *fi=phi1;

            /// *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {

              const double *fj=phi1;
              /// *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
                for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[2*dim]][i*nve1+j] -= dt*((*(gradfi+idim))*(*fj))*Weight;
		}
              } // end phi_j loop
            } // end phi_i loop
          } // End Bt
        
          { // Divergence of the Velocity operator
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {

              //BEGIN RESIDUALS B block ===========================
              Rhs[indexVAR[2*dim]][i] += -(-((*fi))*div_vel)*Weight;
              //END RESIDUALS  B block ===========================

              const double *gradfj=&gradphi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[2*dim]][indexVAR[idim+dim]][i*nve+j] -= ((*fi)*(*(gradfj+idim)))*Weight;
		}
              }
            }
          }
	}   
	//END FLUID ASSEMBLY ============
	//*******************************************************************************************************
	//BEGIN SOLID ASSEMBLY ============
	  
	else{
	  
	  //------------------------------------------------------------------------------------------------------------
          if (solid_model==0) {
	    double e[3][3];
	    double e_old[3][3];
	    //computation of the stress tensor
	    for(int i=0;i<dim;i++){
	      for(int j=0;j<dim;j++){
		e[i][j]=0.5*(GradSolhatVAR[i][j]+GradSolhatVAR[j][i]);
		e_old[i][j]=0.5*(GradSolOldhatVAR[i][j]+GradSolOldhatVAR[j][i]);
	      }
	    }
	    I_e=0;
	    for(int i=0;i<dim;i++){
	      I_e += e[i][i];
	    }
	    
	    for (int i=0; i<dim; i++) {
              for (int j=0; j<dim; j++) {
                //incompressible
                Cauchy[i][j] = 2*mus*e[i][j];
		Cauchy_old[i][j]=2*mus*e_old[i][j];
              }
            }
          }

          else if (solid_model==1) {
	    double F[3][3]={{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
	    double b_left[3][3];      
	    
	    for(int i=0;i<dim;i++){
	      for(int j=0;j<dim;j++){
		F[i][j]+=GradSolhatVAR[i][j];
	      }
	    }
	    
	    Jnp1_hat =  F[0][0]*F[1][1]*F[2][2] + F[0][1]*F[1][2]*F[2][0] + F[0][2]*F[1][0]*F[2][1]
		      - F[2][0]*F[1][1]*F[0][2] - F[2][1]*F[1][2]*F[0][0] - F[2][2]*F[1][0]*F[0][1];		      
	   
	    // computation of the the three deformation tensor b
	    for (int I=0; I<3; ++I) {
	      for (int J=0; J<3; ++J) {
		b_left[I][J]=0.;
		for (int K=0; K<3; ++K) {
		  //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
		  b_left[I][J] += F[I][K]*F[J][K];
		}
		Cauchy[I][J] = (mus/Jnp1_hat)*(b_left[I][J] - Id2th[I][J]);
	      }
	    }

	    I_bleft = b_left[0][0] + b_left[1][1] + b_left[2][2];
	    
	    //compressible case
	    //             for (int ii=0; ii<3; ++ii) {
	    //               for (int jj=0; jj<3; ++jj) {
	    //                 //Cauchy stress tensor
	    //                 Cauchy[ii][jj] = (_mus/Jnp1_hat)*(b_left[ii][jj] - Id2th[ii][jj]) + (_lambda/Jnp1_hat)*log(Jnp1_hat)*Id2th[ii][jj];
	    //                 for (int k=0; k<3; ++k) {
	    //                   for (int l=0; l<3; ++l) {
	    //                   }
	    //                 }
	    //               }
	    //             }
	    
	    //for the incompressible(nearly incompressible) case
	    for (int ii=0; ii<3; ++ii) {
	      for (int jj=0; jj<3; ++jj) {
		for (int kk=0; kk<3; ++kk) {
		  for (int ll=0; ll<3; ++ll) {
		    C_mat[ii][jj][kk][ll] = 2.*mus*pow(Jnp1_hat,-1.6666666666666)*(
										  0.333333333333*I_bleft*Id2th[ii][kk]*Id2th[jj][ll]              //1/3*I_c*i
										  // 	                        +0.111111111111*I_C*Id2th[i][j]*Id2th[k][l]             //1/9*I_b*IxI
										  // 				-0.333333333333*b_left[i][j]*Id2th[k][l]                //-1/3*b*I
										  // 				-0.333333333333*Id2th[i][j]*b_left[k][l]                //-1/3*b*I
										  )
		      -SolVAR[2*dim]*(Id2th[ii][jj]*Id2th[kk][ll]-2.*Id2th[ii][kk]*Id2th[jj][ll] );  // -p(IxI-2i)
		  }
		}
	      }
	    }
	    
	    //Old deformation gradient
	    double F_old[3][3]={{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
	    for(int i=0;i<dim;i++){
	      for(int j=0;j<dim;j++){
		F_old[i][j]+=GradSolhatVAR[i][j];
	      }
	    }
            
            Jn_hat =  F_old[0][0]*F_old[1][1]*F_old[2][2] + F_old[0][1]*F_old[1][2]*F_old[2][0] + F_old[0][2]*F_old[1][0]*F_old[2][1]
		    - F_old[2][0]*F_old[1][1]*F_old[0][2] - F_old[2][1]*F_old[1][2]*F_old[0][0] - F_old[2][2]*F_old[1][0]*F_old[0][1] ;

            // computation of the the three deformation tensor b
            for (int I=0; I<3; ++I) {
              for (int J=0; J<3; ++J) {
                b_left[I][J]=0.;
                for (int k=0; k<3; ++k) {
                  //left Cauchy-Green deformation tensor or F_oldinger tensor (b = F_old*F_old^T)
                  b_left[I][J] += F_old[I][k]*F_old[J][k];
                }
                Cauchy_old[I][J] = (mus/Jn_hat)*(b_left[I][J] - Id2th[I][J]);
              }
            }
          }
          //----------------------------------------------------------------------------------------------------------------------------

          /////////////

          ///Mass + Stiffness operator
          {
            const double *gradfi=&gradphi[0];
            const double *fi=&phi[0];

            /// *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {


              //BEGIN RESIDUALS A + Bt block ===========================
	      
	      // Residual ALE equations
 	      for(int idim=0; idim<dim; idim++) {
	        Rhs[indexVAR[idim]][i] += (
					   -phi[i]*eps_pen*(
							    +SolVAR[idim] - SolOldVAR[idim]
							    -dt*gamma*SolVAR[dim+idim]
							    -dt*(1.-gamma)*SolOldVAR[dim+idim]
							    )
					   )*Weight_hat;
              }
              
              double CauchyDIR[3]={0.,0.,0.};
	      for(int idim=0.; idim<dim; idim++) {
		for(int idim2=0.; idim2<dim; idim2++) {
		  CauchyDIR[idim]+= gradphi[i*dim+idim2]*Cauchy[idim][idim2];
		}
	      }

              // Residual Momentum equations
              for(int idim=0; idim<dim; idim++) {
	        Rhs[indexVAR[dim+idim]][i] += (
					       phi[i]*dt*_gravity[idim]*Weight_hat
					       -phi[i]*betafsi*(SolVAR[dim+idim] - SolOldVAR[dim+idim])*Weight_hat*(1./gamma)
					       +phi[i]*betafsi*SolVAR[2*dim+1+idim]*Weight_hat*gammaratio*dt
					       -dt*CauchyDIR[idim]*Weight
					       +dt*SolVAR[2*dim]*gradphi[i*dim+idim]*Weight
					       );

              }
              
              //---------------------------------------------------------------------------------------------------------------------------------

              //END RESIDUALS A + Bt block ===========================

              const double *gradfj=&gradphi[0];
              const double *fj=&phi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim,fj++) {
                /// Mass term
                // (v_n+1,psi)
                //gamma = 1 Backward Euler
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[dim+idim]][i*nve+j] += betafsi*(*(fi))*(*(fj))*Weight_hat*(1./gamma);
		}

                //Da collaudare il 3D
                for (int icount=0; icount<dim; ++icount) {
                  for (int jcount=0; jcount<dim; ++jcount) {
                    tg_stiff_matrix[icount][jcount] = 0.;
                    for (int kcount=0; kcount<dim; ++kcount) {
                      for (int lcount=0; lcount<dim; ++lcount) {
                        tg_stiff_matrix[icount][jcount] += (*(gradfi+kcount))*0.25*(
										    C_mat[icount][kcount][jcount][lcount]+C_mat[icount][kcount][lcount][jcount]
										    +C_mat[kcount][icount][jcount][lcount]+C_mat[kcount][icount][lcount][jcount]
										    )*(*(gradfj+lcount));
                      }
                    }
                  }
                }
                
                //geometric tangent stiffness matrix
                double geom_tg_stiff_matrx = 0.;
                for(int kcount=0; kcount<dim; ++kcount) {
                  for(int lcount=0; lcount<dim; ++lcount) {
                    geom_tg_stiff_matrx += (*(gradfi+kcount))*Cauchy[kcount][lcount]*(*(gradfj+lcount));
                  }
                }

                /// Stiffness operator -- Elasticity equation (Linear or not)
                for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[0+idim]][i*nve+j] += dt*geom_tg_stiff_matrx*Weight;
 		  for(int idim2=0; idim2<dim; idim2++) {
 		    B[indexVAR[dim+idim]][indexVAR[0+idim2]][i*nve+j] += dt*tg_stiff_matrix[0+idim][0+idim2]*Weight;
 		  }
 		}

                /// Kinematic equation v = du/dt
                //   -_theta*dt*(v_n+1,eta)
		for(int idim=0; idim<dim; idim++) {
		  //   -_theta*dt*(v_n+1,eta)
		  B[indexVAR[0+idim]][indexVAR[dim+idim]][i*nve+j] -= eps_pen*gamma*dt*(*(fi))*(*(fj))*Weight_hat;
		  // (u_n+1,eta)
		  B[indexVAR[0+idim]][indexVAR[0+idim]][i*nve+j] += eps_pen*(*(fi))*(*(fj))*Weight_hat;
		}
	      }
            }
          }
          ////////////
          { ///Gradient of Pressure
            const double *gradfi=&gradphi[0];
            // *** phi_i loop ***
            for (unsigned i=0; i<nve; i++,gradfi+=dim) {
              const double *fj=phi1;
              // *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[dim+idim]][indexVAR[2*dim]][i*nve1+j] -= dt*((*(gradfi+idim))*(*fj))*Weight;
		}
              }
            }
          }
          ////////////
          { ///Divergence of the Displacement
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {

              //BEGIN RESIDUALS B block ===========================

              if (solid_model==0) {
                Rhs[indexVAR[2*dim]][i] += -(-((*fi))*(I_e + (1./lambda)*SolVAR[2*dim] ) )*Weight_hat;
              }
              else if (solid_model==1) {
                Rhs[indexVAR[2*dim]][i] += -(-((*fi))*( log(Jnp1_hat)/Jnp1_hat + (1./lambda)*SolVAR[2*dim] ) )*Weight_hat;
              }

              //END RESIDUALS B block ===========================

              const double *gradfj=&gradphi[0];
              // *** phi_j loop ***
              for (unsigned j=0; j<nve; j++,gradfj+=dim) {
		for(int idim=0; idim<dim; idim++) {
		  B[indexVAR[2*dim]][indexVAR[idim]][i*nve+j] -= ((*fi)*(*(gradfj+idim)))*Weight;
		}
              }
            }
          }
          //  /////////////
          {  ///Pressure Mass term
            const double *fi=phi1;
            // *** phi_i loop ***
            for (unsigned i=0; i<nve1; i++,fi++) {
              const double *fj=phi1;
              // *** phi_j loop ***
              for (unsigned j=0; j<nve1; j++,fj++) {
                B[indexVAR[2*dim]][indexVAR[2*dim]][i*nve1+j] -= (1./lambda)*((*fi)*(*fj))*Weight_hat;
              }
            }
          }  //end pressure mass term
	  //---------------------------------------------------------------------------------------------------------------------------------
	}  
	//END SOLID ASSEMBLY ============
      }
    }

    //BEGIN local to global assembly 
    // ALE mapping
    for(int i=0; i<dim; i++) {
      myRES->add_vector_blocked(Rhs[indexVAR[i]],dofsVAR[i]);
      myKK ->add_matrix_blocked(B[indexVAR[i]][indexVAR[i]],dofsVAR[i],dofsVAR[i]);  
      if(flag_mat!=2){ //Solid only
	myKK->add_matrix_blocked(B[indexVAR[i]][indexVAR[i+dim]],dofsVAR[i],dofsVAR[i+dim]);
      }
    }
    
    // Momentum equation
    for(int i=0; i<dim; i++) {
      myRES->add_vector_blocked(Rhs[indexVAR[dim+i]],dofsVAR[dim+i]);
      myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[dim+i]],dofsVAR[dim+i],dofsVAR[dim+i]);
      myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[2*dim]],dofsVAR[dim+i],dofsVAR[2*dim]);
      for(int j=0; j<dim; j++) {
	myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[j]],dofsVAR[dim+i],dofsVAR[j]);
      }
    }
   
    //P-continuity equation
    myRES->add_vector_blocked(Rhs[indexVAR[2*dim]],dofsVAR[2*dim]);
    for(int i=0; i<dim; i++) {
      if(flag_mat==2){ //Fluid only
	myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[dim+i]],dofsVAR[2*dim],dofsVAR[dim+i]);
      }
      else{ //Solid only
	myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[i]],dofsVAR[2*dim],dofsVAR[i]);  
      } 
    }
    myKK->add_matrix_blocked(B[indexVAR[2*dim]][indexVAR[2*dim]],dofsVAR[2*dim],dofsVAR[2*dim]);
    //END local to global assembly
   
  } //end list of elements loop

  // close residual vector and matrix

  myKK->close();
  myRES->close();
  
  // *************************************
  end_time=clock();
  AssemblyTime+=(end_time-start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

}

} //end namespace femus

#endif 