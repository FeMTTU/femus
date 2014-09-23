#ifndef __IncompFSIassembly_hpp__
#define __IncompFSIassembly_hpp__

#include "MonolithicFSINonLinearImplicitSystem.hpp" 

namespace femus {
 
  void IncompressibleFSIAssembly(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix) {
    
    clock_t AssemblyTime=0;
    clock_t start_time, end_time;
  
    //pointers and references
    MultiLevelSolution*	 ml_sol	                      = ml_prob._ml_sol;
    Solution*	 mysolution  	                      = ml_sol->GetSolutionLevel(level);
    MonolithicFSINonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<MonolithicFSINonLinearImplicitSystem>("Fluid-Structure-Interaction");
    LinearEquationSolver*  myLinEqSolver	              = my_nnlin_impl_sys._LinSolver[level];   
    mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(level);
    elem		*myel		=  mymsh->el;
    SparseMatrix	*myKK		=  myLinEqSolver->_KK;
    NumericVector 	*myRES		=  myLinEqSolver->_RES;
    vector <int>	&myKKIndex	=  myLinEqSolver->KKIndex;
    
    const unsigned dim = mymsh->GetDimension();
    const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  
    // local objects
    vector<double> SolVAR(2*dim+1);
    vector<vector<double> > GradSolVAR(2*dim);
    for(int i=0;i<2*dim;i++){
      GradSolVAR[i].resize(dim);
    }
    vector<vector<double> > GradSolhatVAR(dim);
    for(int i=0;i<dim;i++){
      GradSolhatVAR[i].resize(dim);
    }
    
    vector <int> metis_node1;
    vector <int> metis_node2;
    vector <bool> solidmark;
  
    vector <double > phi;
    vector <double > phi_hat;
  
    vector <double> gradphi;
    vector <double> gradphi_hat;
  
    metis_node1.reserve(max_size);
    metis_node2.reserve(max_size);
    solidmark.reserve(max_size);
    phi.reserve(max_size);
    phi_hat.reserve(max_size);
    gradphi.reserve(max_size*dim);
    gradphi_hat.reserve(max_size*dim);
  
    const double *phi1;
    
    double Weight=0.;
    double Weight_nojac=0.;
    double Weight_hat=0.;
  
    vector <vector < double> > vx(dim);
    vector <vector < double> > vx_hat(dim);
  
    for(int i=0;i<dim;i++){
      vx[i].reserve(max_size);
      vx_hat[i].reserve(max_size);
    }
   
    vector< vector< double > > Rhs(2*dim+1);
    vector< vector< vector< double > > > B(2*dim+1); 
    for(int i=0;i<2*dim+1;i++){
      B[i].resize(2*dim+1);
    }
    vector< vector< int > > dofsVAR(2*dim+1); 
  
    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof	 	= ml_prob.parameters.get<Fluid>("Fluid").get_density();            
    double rhos 		= ml_prob.parameters.get<Solid>("Solid").get_density();            
    double mu_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_shear_modulus(); 
    double lambda_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_lambda();        
    double mus		= mu_lame/rhof;
    double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();   
    double lambda		= lambda_lame / rhof;
    double betafsi	= rhos / rhof;
    double betans		= 1.;
    int    solid_model	= ml_prob.parameters.get<Solid>("Solid").get_physical_model();     

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
  
    // newton algorithm
    bool nwtn_alg = true;

    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned order_ind2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));  
    unsigned end_ind2   = mymsh->GetEndIndex(order_ind2);

    unsigned order_ind1 = ml_sol->GetSolutionType(ml_sol->GetIndex("P"));  
    unsigned end_ind1   = mymsh->GetEndIndex(order_ind1);

    // mesh and procs
    unsigned nel    = mymsh->GetElementNumber();
    unsigned igrid  = mymsh->GetGridNumber();
    unsigned iproc  = mymsh->processor_id();

    //----------------------------------------------------------------------------------
    //variable-name handling
    const char varname[7][3] = {"DX","DY","DZ","U","V","W","P"};
    vector <unsigned> indexVAR(2*dim+1);
    vector <unsigned> indVAR(2*dim+1);  
    vector <unsigned> SolType(2*dim+1);  
  
    for(unsigned ivar=0; ivar<dim; ivar++) {
      indVAR[ivar]=ml_sol->GetIndex(&varname[ivar][0]);
      indVAR[ivar+dim]=ml_sol->GetIndex(&varname[ivar+3][0]);
      SolType[ivar]=ml_sol->GetSolutionType(&varname[ivar][0]);
      SolType[ivar+dim]=ml_sol->GetSolutionType(&varname[ivar+3][0]);
      indexVAR[ivar]=my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);
      indexVAR[ivar+dim]=my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar+3][0]);
    }
    indexVAR[2*dim]=my_nnlin_impl_sys.GetSolPdeIndex(&varname[6][0]);
    indVAR[2*dim]=ml_sol->GetIndex(&varname[6][0]);
    SolType[2*dim]=ml_sol->GetSolutionType(&varname[6][0]);
    //----------------------------------------------------------------------------------
      
    int nprocs=mymsh->n_processors();
    NumericVector *area_elem_first;
    area_elem_first = NumericVector::build().release();
   
    if(nprocs==1) { 
      area_elem_first->init(nprocs,1,false,SERIAL);     
    } 
    else { 
     area_elem_first->init(nprocs,1,false,PARALLEL); 	   	   
    }
    area_elem_first->zero();
    double rapresentative_area=1.;
    
    start_time=clock();
    
    myKK->zero();
    
    /// *** element loop ***
    for(int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

	unsigned kel        = mymsh->IS_Mts2Gmt_elem[iel]; 
	short unsigned kelt = myel->GetElementType(kel);
	unsigned nve        = myel->GetElementDofNumber(kel,end_ind2);
	unsigned nve1       = myel->GetElementDofNumber(kel,end_ind1);
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
	
	    // 	B[indexVAR[dim+i]][indexVAR[dim+j]].resize(nve*nve);
	    // 	memset(&B[indexVAR[dim+i]][indexVAR[dim+j]][0],0,nve*nve*sizeof(double));
	  }
	}
      
	//Mass Matrix (solid and fluid) and Diffusion Matrix (fluid only) 
	for(int i=0; i<dim; i++) {
	  B[indexVAR[dim+i]][indexVAR[dim+i]].resize(nve*nve);
	  memset(&B[indexVAR[dim+i]][indexVAR[dim+i]][0],0,nve*nve*sizeof(double));
	  if(nwtn_alg== true) {
	    for(int jdim=1; jdim<dim; jdim++) {
	      B[indexVAR[dim+i]][indexVAR[dim+(i+jdim)%dim]].resize(nve*nve);
	      memset(&B[indexVAR[dim+i]][indexVAR[dim+(i+jdim)%dim]][0],0,nve*nve*sizeof(double));
	    }
	  }
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
	gradphi.resize(nve*dim);
	gradphi_hat.resize(nve*dim);
        
	for(int i=0;i<dim;i++){
	  vx[i].resize(nve);
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
	    //Fixed coordinates (Reference frame)
	    vx_hat[j][i]= (*mymsh->_coordinate->_Sol[j])(inode_Metis);  
	    // displacement dofs
	    dofsVAR[j][i]= myLinEqSolver->GetKKDof(indVAR[j],indexVAR[j],inode); 
	    // velocity dofs
	    dofsVAR[j+dim][i]= myLinEqSolver->GetKKDof(indVAR[j+dim],indexVAR[j+dim],inode);   
	  }
	}

	// pressure dofs
	for (unsigned i=0;i<nve1;i++) {
	  unsigned inode=(order_ind1<3)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
	  metis_node1[i]=mymsh->GetMetisDof(inode,SolType[2*dim]);
	  dofsVAR[2*dim][i]=myLinEqSolver->GetKKDof(indVAR[2*dim],indexVAR[2*dim],inode);
	}
	// ----------------------------------------------------------------------------------------
       
	if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
	  /// *** Gauss point loop ***
	  double area=1.;
	  for (unsigned ig=0;ig < ml_prob._ml_msh->_type_elem[kelt][order_ind2]->GetGaussPointNumber(); ig++) {

	    // *** get Jacobian and test function and test function derivatives in the moving frame***
	    (ml_prob._ml_msh->_type_elem[kelt][order_ind2]->*(ml_prob._ml_msh->_type_elem[kelt][order_ind2])->Jacobian_ptr)(vx,ig,Weight,phi,gradphi);
	    (ml_prob._ml_msh->_type_elem[kelt][order_ind2]->*(ml_prob._ml_msh->_type_elem[kelt][order_ind2])->Jacobian_ptr)(vx_hat,ig,Weight_hat,phi_hat,gradphi_hat);
	    phi1=ml_prob._ml_msh->_type_elem[kelt][order_ind1]->GetPhi(ig);
	    if (flag_mat==2) {
	      if(ig==0){
		double GaussWeight = ml_prob._ml_msh->_type_elem[kelt][order_ind2]->GetGaussWeight(ig);
		area=Weight_hat/GaussWeight;
		if(iel==mymsh->IS_Mts2Gmt_elem_offset[iproc]){
		  area_elem_first->add(mymsh->processor_id(),area);
		  area_elem_first->close();
		  rapresentative_area=area_elem_first->l1_norm()/nprocs;
		}
	      }
	      Weight_nojac = Weight_hat/area*rapresentative_area;
	    }
	    // ---------------------------------------------------------------------------
	    // displacement and velocity
	    for(int i=0; i<2*dim; i++){
	      SolVAR[i]=0.;
	      for(int j=0; j<dim; j++) {
		GradSolVAR[i][j]=0.;
		if(i<dim){
		  GradSolhatVAR[i][j]=0.;
		}
	      }
	    
	      for (unsigned inode=0; inode<nve; inode++) {
		unsigned sol_dof = metis_node2[inode];
	      
		double soli = (*mysolution->_Sol[indVAR[i]])(sol_dof);
		SolVAR[i]+=phi[inode]*soli;
	      
		for(int j=0; j<dim; j++) {
		  GradSolVAR[i][j]+=gradphi[inode*dim+j]*soli;
		  if(i<dim){ 
		    GradSolhatVAR[i][j]   +=gradphi_hat[inode*dim+j]*soli;
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
  
	    // ---------------------------------------------------------------------------
	 	  
	    //BEGIN FLUID ASSEMBLY ============
	  
	    if(flag_mat==2){
	  
	      { // Laplace operator + adection operator + Mass operator

		const double *gradfi=&gradphi[0];
		const double *gradfi_hat=&gradphi_hat[0];
		const double *fi=&phi[0];

		// *** phi_i loop ***
		for (unsigned i=0; i<nve; i++,gradfi+=dim,gradfi_hat+=dim,fi++) {

		  //BEGIN RESIDUALS A + Bt block ===========================
 	      
		  // begin redidual Laplacian ALE map in the reference domain
		  double LapmapVAR[3] = {0., 0., 0.};
		  for(int idim=0; idim<dim; idim++) {
		    for(int jdim=0; jdim<dim; jdim++) {
		      LapmapVAR[idim] += _mu_ale[jdim]*(GradSolhatVAR[idim][jdim]*gradphi_hat[i*dim+jdim]) ;
		    }
		  }
		  for(int idim=0; idim<dim; idim++) {
		    Rhs[indexVAR[idim]][i]+=(!solidmark[i])*(-LapmapVAR[idim]*Weight_nojac);
		  }
		  // end redidual Laplacian ALE map

		  // begin redidual Navier-Stokes in the moving domain   
		  double LapvelVAR[3]={0.,0.,0.};
		  double AdvaleVAR[3]={0.,0.,0.};
		  for(int idim=0.; idim<dim; idim++) {
		    for(int jdim=0.; jdim<dim; jdim++) {
		      LapvelVAR[idim]+=GradSolVAR[dim+idim][jdim]*gradphi[i*dim+jdim];
		      AdvaleVAR[idim]+=SolVAR[dim+jdim]*GradSolVAR[dim+idim][jdim]*phi[i];
		    }
		  }
		  // Residual Momentum equations 
		  for(int idim=0; idim<dim; idim++) {
		    Rhs[indexVAR[dim+idim]][i]+= (-AdvaleVAR[idim]			 // advection term	
						  -IRe*LapvelVAR[idim]		 // viscous dissipation
						  +SolVAR[2*dim]*gradphi[i*dim+idim] // pressure gradient
						  )*Weight;
		    // End redidual Navier-Stokes    
		  }
		  //END RESIDUALS A + Bt block ===========================
	      
		  //BEGIN A block ===========================
		  const double *gradfj=&gradphi[0];
		  const double *gradfj_hat=&gradphi_hat[0];
		  const double *fj=&phi[0];
		  //  *** phi_j loop ***
		  for (unsigned j=0; j<nve; j++,gradfj+=dim,gradfj_hat+=dim,fj++) {

		    // begin Laplacian ALE map
		    double Lap_ale=0.;
		    for(int jdim=0; jdim<dim; jdim++) {
		      Lap_ale+=_mu_ale[jdim]*(*(gradfi_hat+jdim))*(*(gradfj_hat+jdim));
		    }	
		    for(int idim=0; idim<dim; idim++) {
		      B[indexVAR[idim]][indexVAR[idim]][i*nve+j] += (!solidmark[i])*Lap_ale*Weight_nojac;
		    }
		    // end Laplacian ALE map
		
		    //Laplacian
		    double Lap=0.;
		    for(int jdim=0; jdim<dim; jdim++) {
		      Lap+=(*(gradfi+jdim))*(*(gradfj+jdim))*Weight;
		    }
               
		    //advection term I
		    double Adv1 = 0.;
		    for(int jdim=0; jdim<dim; jdim++) {
		      Adv1+= SolVAR[dim+jdim]*(*(gradfj+jdim))*(*(fi))*Weight;
		    }
		
		    //advection term II
		    double Adv2 = ((*fi))*((*fj))*Weight;
				
		    for(int idim=0; idim<dim; idim++) {
		      B[indexVAR[dim+idim]][indexVAR[dim+idim]][i*nve+j] += (IRe*Lap + Adv1);
		      // Advection term II
		      if(nwtn_alg== true) {
			B[indexVAR[dim+idim]][indexVAR[dim+idim]][i*nve+j]  += Adv2*GradSolVAR[dim+idim][idim];
			for(unsigned jdim=1; jdim<dim; jdim++) {
			  B[indexVAR[dim+idim]][indexVAR[dim + (idim+jdim)%dim]][i*nve+j] += Adv2*GradSolVAR[dim+idim][(idim+jdim)%dim];
			}
		      }
		    }
		  } // end phi_j loop
		  //END A block ===========================
		} // end phi loop
	      } // end A 
	 
	 
	      //BEGIN Bt block ===========================
	      { //Gradient of Pressure operator

		const double *gradfi=&gradphi[0];
		const double *fi=phi1;
		/// *** phi_i loop ***
		  for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {
		    const double *fj=phi1;
		    /// *** phi_j loop ***
		      for (unsigned j=0; j<nve1; j++,fj++) {
			for(int idim=0; idim<dim; idim++) {
			  B[indexVAR[dim+idim]][indexVAR[2*dim]][i*nve1+j] -= ((*(gradfi+idim))*(*fj))*Weight;
			}
		      } // end phi_j loop
		  } // end phi_i loop
	      } // End Bt
	      //END Bt block ===========================
          
	      //BEGIN B block ===========================
	      {  	    
		//divergence of the velocity
		double div_vel=0.;
		for(int i=0; i<dim; i++) {
		  div_vel+=GradSolVAR[dim+i][i];
		}
		// Divergence of the Velocity operator
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
	      //END B block ===========================
	    }   
	    //END FLUID ASSEMBLY ============
	    
	    //*******************************************************************************************************
	    
	    //BEGIN SOLID ASSEMBLY ============
	    else{
	      //------------------------------------------------------------------------------------------------------------
	      if (solid_model==0) { //incompressible Neo-Hookean material
		double e[3][3];
		//computation of the stress tensor
		for(int i=0;i<dim;i++){
		  for(int j=0;j<dim;j++){
		    e[i][j]=0.5*(GradSolhatVAR[i][j]+GradSolhatVAR[j][i]);
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
		    Cauchy[I][J] = mus*(b_left[I][J] - Id2th[I][J]);
		  }
		}
		// C_Mat: newton-raphson Jacobi for Cauchy, to be used as gradfj_jj*C_mat[jdim][jj][idim][ii]*gradfi_ii
		C_mat[0][0][0][0]= 2.*F[0][0];	C_mat[0][0][0][1]= 1.*F[1][0]; 	C_mat[0][0][0][2]= 1.*F[2][0];
		C_mat[0][0][1][0]= 1.*F[1][0]; 	C_mat[0][0][1][1]= 0.;	 	C_mat[0][0][1][2]= 0.;
		C_mat[0][0][2][0]= 1.*F[2][0]; 	C_mat[0][0][2][1]= 0.;	  	C_mat[0][0][2][2]= 0.;
	    
		C_mat[0][1][0][0]= 2.*F[0][1]; 	C_mat[0][1][0][1]= 1.*F[1][1]; 	C_mat[0][1][0][2]= 1.*F[2][1];
		C_mat[0][1][1][0]= 1.*F[1][1]; 	C_mat[0][1][1][1]= 0.;	  	C_mat[0][1][1][2]= 0.;
		C_mat[0][1][2][0]= 1.*F[2][1]; 	C_mat[0][1][2][1]= 0.;	  	C_mat[0][1][2][2]= 0.;
	    
		C_mat[0][2][0][0]= 2.*F[0][2]; 	C_mat[0][2][0][1]= 1.*F[1][2]; 	C_mat[0][2][0][2]= 1.*F[2][2];
		C_mat[0][2][1][0]= 1.*F[1][2]; 	C_mat[0][2][1][1]= 0.;	  	C_mat[0][2][1][2]= 0.;
		C_mat[0][2][2][0]= 1.*F[2][2]; 	C_mat[0][2][2][1]= 0.;	  	C_mat[0][2][2][2]= 0.;
	    
		C_mat[1][0][0][0]= 0.;	  	C_mat[1][0][0][1]= 1.*F[0][0]; 	C_mat[1][0][0][2]= 0.;
		C_mat[1][0][1][0]= 1.*F[0][0]; 	C_mat[1][0][1][1]= 2.*F[1][0]; 	C_mat[1][0][1][2]= 1.*F[2][0];
		C_mat[1][0][2][0]= 0.;	   	C_mat[1][0][2][1]= 1.*F[2][0]; 	C_mat[1][0][2][2]= 0.;
	    	    
		C_mat[1][1][0][0]= 0.;	   	C_mat[1][1][0][1]= 1.*F[0][1]; 	C_mat[1][1][0][2]= 0.;
		C_mat[1][1][1][0]= 1.*F[0][1]; 	C_mat[1][1][1][1]= 2.*F[1][1]; 	C_mat[1][1][1][2]= 1.*F[2][1];
		C_mat[1][1][2][0]= 0.;	   	C_mat[1][1][2][1]= 1.*F[2][1]; 	C_mat[1][1][2][2]= 0.;
	    
		C_mat[1][2][0][0]= 0.;	   	C_mat[1][2][0][1]= 1.*F[0][2]; 	C_mat[1][2][0][2]= 0.;
		C_mat[1][2][1][0]= 1.*F[0][2]; 	C_mat[1][2][1][1]= 2.*F[1][2]; 	C_mat[1][2][1][2]= 1.*F[2][2];
		C_mat[1][2][2][0]= 0.;	   	C_mat[1][2][2][1]= 1.*F[2][2]; 	C_mat[1][2][2][2]= 0.;
	    
		C_mat[2][0][0][0]= 0.;	   	C_mat[2][0][0][1]= 0.;	  	C_mat[2][0][0][2]= 1.*F[0][0];
		C_mat[2][0][1][0]= 0.;	   	C_mat[2][0][1][1]= 0.;	  	C_mat[2][0][1][2]= 1.*F[1][0];
		C_mat[2][0][2][0]= 1.*F[0][0]; 	C_mat[2][0][2][1]= 1.*F[1][0]; 	C_mat[2][0][2][2]= 2.*F[2][0];
	    
		C_mat[2][1][0][0]= 0.;	  	C_mat[2][1][0][1]= 0.;    	C_mat[2][1][0][2]= 1.*F[0][1];
		C_mat[2][1][1][0]= 0.;	  	C_mat[2][1][1][1]= 0.; 	  	C_mat[2][1][1][2]= 1.*F[1][1];
		C_mat[2][1][2][0]= 1.*F[0][1]; 	C_mat[2][1][2][1]= 1.*F[1][1];	C_mat[2][1][2][2]= 2.*F[2][1];
	    
		C_mat[2][2][0][0]= 0.;	   	C_mat[2][2][0][1]= 0.;    	C_mat[2][2][0][2]= 1.*F[0][2];
		C_mat[2][2][1][0]= 0.;	   	C_mat[2][2][1][1]= 0.; 	  	C_mat[2][2][1][2]= 1.*F[1][2];
		C_mat[2][2][2][0]= 1.*F[0][2]; 	C_mat[2][2][2][1]= 1.*F[1][2]; 	C_mat[2][2][2][2]= 2.*F[2][2];
	    
		for (int ii=0; ii<3; ++ii) {
		  for (int jj=0; jj<3; ++jj) {
		    for (int kk=0; kk<3; ++kk) {
		      for (int ll=0; ll<3; ++ll) {
			C_mat[ii][jj][kk][ll]*=mus;
		      }
		    }
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
	      
		      // Residual ALE equations in the reference domain (remeber that phi=phi_hat)
		      for(int idim=0; idim<dim; idim++) {
			Rhs[indexVAR[idim]][i] += (-phi[i]*(-SolVAR[dim+idim] ))*Weight_hat;
		      }
              
		      double CauchyDIR[3]={0.,0.,0.};
		      for(int idim=0.; idim<dim; idim++) {
			for(int jdim=0.; jdim<dim; jdim++) {
			  CauchyDIR[idim]+= gradphi[i*dim+jdim]*Cauchy[idim][jdim];
			}
		      }

		      // Residual Momentum equations
		      for(int idim=0; idim<dim; idim++) {
			Rhs[indexVAR[dim+idim]][i] += (
						       -CauchyDIR[idim]*Weight
						       +SolVAR[2*dim]*gradphi[i*dim+idim]*Weight
						       );
		      }
              
		      //---------------------------------------------------------------------------------------------------------------------------------

		      //END RESIDUALS A + Bt block ===========================

		      const double *gradfj=&gradphi[0];
		      const double *gradfj_hat=&gradphi_hat[0];
		      const double *fj=&phi[0];
		      // *** phi_j loop ***
		      for (unsigned j=0; j<nve; j++,gradfj+=dim,gradfj_hat+=dim,fj++) {
			
			/// Kinematic equation v = du/dt --> In the steady state we write \deltau^n+1 - \deltav^n+1 = v - 0
			for(int idim=0; idim<dim; idim++) {
			  // -(v_n+1,eta)
			  B[indexVAR[idim]][indexVAR[dim+idim]][i*nve+j] -= (*(fi))*(*(fj))*Weight_hat;
			  //  (u_n+1,eta)
			  B[indexVAR[idim]][indexVAR[idim]][i*nve+j] 	+=  (*(fi))*(*(fj))*Weight_hat;
			}
			
			
			/// Stiffness operator -- Elasticity equation (Linear or not)
			for(int idim=0; idim<dim; idim++) {
			  for(int jdim=0; jdim<dim; jdim++) {
			    double tg_stiff_matrix = 0.;
			    for (int ii=0; ii<dim; ++ii) {
			      for (int jj=0; jj<dim; ++jj) {
				tg_stiff_matrix += (*(gradfj_hat+jj))*(C_mat[jdim][jj][idim][ii])*(*(gradfi+ii));
			      }
			    }
			    B[indexVAR[dim+idim]][indexVAR[jdim]][i*nve+j] += tg_stiff_matrix*Weight;
			  }
			}
			
		      }
		    }
		}
		// ********************************************
		{ ///Gradient of Pressure
		  const double *gradfi=&gradphi[0];
		  // *** phi_i loop ***
		  for (unsigned i=0; i<nve; i++,gradfi+=dim) {
		    const double *fj=phi1;
		    // *** phi_j loop ***
		    for (unsigned j=0; j<nve1; j++,fj++) {
		      for(int idim=0; idim<dim; idim++) {
			B[indexVAR[dim+idim]][indexVAR[2*dim]][i*nve1+j] -= ((*(gradfi+idim))*(*fj))*Weight;
		      }
		    }
		  }
		}
		////////////
		{           
		  const double *fi=phi1;
		  // *** phi_i loop ***
		  for (unsigned i=0; i<nve1; i++,fi++) {
	    
		    //BEGIN RESIDUALS B block ===========================

		    if (solid_model==0) {
		      Rhs[indexVAR[2*dim]][i] += -(-((*fi))*(I_e + (1./lambda)*SolVAR[2*dim] ) )*Weight_hat;
		    }
		    else if (solid_model==1) {
		      Rhs[indexVAR[2*dim]][i] +=(*fi)*(Jnp1_hat-1.)*Weight_hat;
		    }

		    //END RESIDUALS B block ===========================

		    const double *gradfj=&gradphi[0];
		    // *** phi_j loop ***
		    for (unsigned j=0; j<nve; j++,gradfj+=dim) {
		      for(int idim=0; idim<dim; idim++) {
			B[indexVAR[2*dim]][indexVAR[idim]][i*nve+j] -= (*fi)*(*(gradfj+idim))*Weight;
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
		      B[indexVAR[2*dim]][indexVAR[2*dim]][i*nve1+j] -= 0.;
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
	  if(nwtn_alg== true){
	    for(unsigned jdim=1; jdim<dim; jdim++) {
	      myKK->add_matrix_blocked(B[indexVAR[dim+i]][indexVAR[dim+(i+jdim)%dim]],dofsVAR[dim+i],dofsVAR[dim+(i+jdim)%dim]);  
	    }
	  }
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
  
      delete area_elem_first;
       
      // *************************************
      end_time=clock();
      AssemblyTime+=(end_time-start_time);
      // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

  }  
    
}

#endif
