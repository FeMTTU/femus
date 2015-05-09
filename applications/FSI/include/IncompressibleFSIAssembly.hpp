#ifndef __femus_include_IncompressibleFSIAssembly_hpp__
#define __femus_include_IncompressibleFSIAssembly_hpp__

#include "MonolithicFSINonLinearImplicitSystem.hpp" 
#include "adept.h"


namespace femus {
    
  void IncompressibleFSIAssemblyAD_DD(MultiLevelProblem &ml_prob) {
       
    clock_t AssemblyTime=0;
    clock_t start_time, end_time;
  
    adept::Stack & s = FemusInit::_adeptStack;
    //    static adept::Stack s; 
    
    //pointers and references
    
    MonolithicFSINonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<MonolithicFSINonLinearImplicitSystem>("Fluid-Structure-Interaction");
    const unsigned level = my_nnlin_impl_sys.GetLevelToAssemble();
    const unsigned gridn = my_nnlin_impl_sys.GetLevelMax();
    bool assemble_matrix = my_nnlin_impl_sys.GetAssembleMatrix(); 
    
    MultiLevelSolution*	 ml_sol	                      = ml_prob._ml_sol;
    Solution*	 mysolution  	                      = ml_sol->GetSolutionLevel(level);
    LinearEquationSolver*  myLinEqSolver	              = my_nnlin_impl_sys._LinSolver[level];   
    Mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(level);
    elem		*myel		=  mymsh->el;
    SparseMatrix	*myKK		=  myLinEqSolver->_KK;
    NumericVector 	*myRES		=  myLinEqSolver->_RES;
        
    const unsigned dim = mymsh->GetDimension();
    const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  
    // local objects
    vector<adept::adouble> SolVAR(2*dim+1);
    vector<vector<adept::adouble> > GradSolVAR(2*dim);
    vector<vector<adept::adouble> > GradSolhatVAR(2*dim);
    
    vector<vector<adept::adouble> > NablaSolVAR(2*dim);
    vector<vector<adept::adouble> > NablaSolhatVAR(2*dim);
    
    for(int i=0;i<2*dim;i++){
      GradSolVAR[i].resize(dim);
      GradSolhatVAR[i].resize(dim);
      
      NablaSolVAR[i].resize(3*(dim-1));
      NablaSolhatVAR[i].resize(3*(dim-1));
    }
   
    vector <bool> solidmark;
    vector <double > phi;
    vector <double > phi_hat;
    vector <adept::adouble> gradphi;
    vector <double> gradphi_hat;
    vector <adept::adouble> nablaphi;
    vector <double> nablaphi_hat;
       
    phi.reserve(max_size);
    solidmark.reserve(max_size);
    phi_hat.reserve(max_size);
    gradphi.reserve(max_size*dim);
    gradphi_hat.reserve(max_size*dim);
    nablaphi.reserve(max_size*3*(dim-1));
    nablaphi_hat.reserve(max_size*3*(dim-1));
    
    //     const double *phi;
   
    
    const double *phi1;
    
    adept::adouble Weight=0.;
    double Weight_nojac=0.;
    double Weight_hat=0.;
  
    vector <vector < adept::adouble> > vx(dim);
    vector <vector < adept::adouble> > vx_face(dim);
    vector <vector < double> > vx_hat(dim);
  
    for(int i=0;i<dim;i++){
      vx[i].reserve(max_size);
      vx_face[i].resize(9);
      vx_hat[i].reserve(max_size);
    }
   
    vector< vector< adept::adouble > > Soli(2*dim+1);
    vector< vector< int > > dofsVAR(2*dim+1); 
    for(int i=0;i<2*dim+1;i++){
      Soli[i].reserve(max_size);
      dofsVAR[i].reserve(max_size);
    }
    
    vector< vector< double > > Rhs(2*dim+1);
    vector< vector< adept::adouble > > aRhs(2*dim+1);
    for(int i=0;i<2*dim+1;i++){
      aRhs[i].reserve(max_size);
      Rhs[i].reserve(max_size);
    }     
    
    
    vector < int > dofsAll;
    dofsAll.reserve(max_size*(2+dim+1));
        
    vector < double > KKloc;
    KKloc.reserve(dim*max_size*(2*dim+1)*dim*max_size*(2*dim+1));
        
    vector < double > Jac;
    Jac.reserve(dim*max_size*(2*dim+1)*dim*max_size*(2*dim+1));
    
    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof	 	= ml_prob.parameters.get<Fluid>("Fluid").get_density();             
    double mu_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_shear_modulus(); 
    double lambda_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_lambda();        
    double mus		= mu_lame/rhof;
    double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();   
    double lambda	= lambda_lame / rhof;
    double betans	= 1.;
    int    solid_model	= ml_prob.parameters.get<Solid>("Solid").get_physical_model();     
       
    bool incompressible=( 0.5 == ml_prob.parameters.get<Solid>("Solid").get_poisson_coeff() )?1:0;
    const bool penalty = ml_prob.parameters.get<Solid>("Solid").get_if_penalty();
    
    // gravity
    double _gravity[3]={0.,0.,0.};
     
    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));  
    unsigned SolType1 = ml_sol->GetSolutionType(ml_sol->GetIndex("P"));  

    // mesh and procs
    unsigned nel    = mymsh->GetNumberOfElements();
    unsigned igrid  = mymsh->GetLevel();
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
    
    // *** element loop ***
    for(int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

      unsigned kel        = mymsh->IS_Mts2Gmt_elem[iel]; 
      short unsigned kelt = myel->GetElementType(kel);
      unsigned nve        = myel->GetElementDofNumber(kel,SolType2);
      unsigned nve1       = myel->GetElementDofNumber(kel,SolType1);
      int flag_mat        = myel->GetElementMaterial(kel);

      // *******************************************************************************************************
    
      //initialization of everything is in common fluid and solid
    
      //Rhs
      for(int i=0; i<2*dim; i++) {
	dofsVAR[i].resize(nve);
	Soli[indexVAR[i]].resize(nve);
	aRhs[indexVAR[i]].resize(nve);
	Rhs[indexVAR[i]].resize(nve);
      }
      dofsVAR[2*dim].resize(nve1);
      Soli[indexVAR[2*dim]].resize(nve1);
      aRhs[indexVAR[2*dim]].resize(nve1);
      Rhs[indexVAR[2*dim]].resize(nve1);
      
      dofsAll.resize(0);
      
      KKloc.resize((2*dim*nve+nve1)*(2*dim*nve+nve1));
      Jac.resize((2*dim*nve+nve1)*(2*dim*nve+nve1));
      
      // ----------------------------------------------------------------------------------------
      // coordinates, solutions, displacement, velocity dofs
    
      solidmark.resize(nve);
        
      for(int i=0;i<dim;i++){
	vx[i].resize(nve);
	vx_hat[i].resize(nve);
      }
    
      for (unsigned i=0;i<nve;i++) {
	unsigned inode=myel->GetMeshDof(kel,i,SolType2);
	unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
	// flag to know if the node "inode" lays on the fluid-solid interface
	solidmark[i]=myel->GetNodeRegion(inode); // to check
	for(int j=0; j<dim; j++) {
	  Soli[indexVAR[j]][i]     =  (*mysolution->_Sol[indVAR[j]])(inode_Metis);
	  Soli[indexVAR[j+dim]][i] =  (*mysolution->_Sol[indVAR[j+dim]])(inode_Metis);
	    	    
	  aRhs[indexVAR[j]][i]     = 0.;
	  aRhs[indexVAR[j+dim]][i] = 0.;
	   
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
	unsigned inode=myel->GetMeshDof(kel,i,SolType1);
	unsigned inode_Metis =mymsh->GetMetisDof(inode,SolType[2*dim]);
	dofsVAR[2*dim][i]=myLinEqSolver->GetKKDof(indVAR[2*dim],indexVAR[2*dim],inode);
	Soli[indexVAR[2*dim]][i] = (*mysolution->_Sol[indVAR[2*dim]])(inode_Metis);
	aRhs[indexVAR[2*dim]][i] = 0.;
      }
      
      // build dof ccomposition             
      for(int idim=0;idim<2*dim;idim++){
	dofsAll.insert( dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end() );
      }
      dofsAll.insert( dofsAll.end(), dofsVAR[2*dim].begin(), dofsVAR[2*dim].end() );
 
      if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {  
	
	s.new_recording();
	
	for (unsigned idim=0;idim<dim;idim++) {
	  for(int j=0; j<nve; j++) {
	    vx[idim][j]= vx_hat[idim][j] + Soli[indexVAR[idim]][j];
	  }
	}
	
	// Boundary integral
	{
	  double tau=0.;
	  vector<adept::adouble> normal(dim,0);
	       
	  // loop on faces
	  for(unsigned jface=0; jface<myel->GetElementFaceNumber(kel); jface++) {
            std::vector < double > xx(3,0.);
	    // look for boundary faces
	    if(myel->GetFaceElementIndex(kel,jface)<0) {
	      unsigned int face = -(mymsh->el->GetFaceElementIndex(kel,jface)+1);	      
	      if( !ml_sol->_SetBoundaryConditionFunction(xx,"U",tau,face,0.) && tau!=0.){
		unsigned nve = mymsh->el->GetElementFaceDofNumber(kel,jface,SolType2);
		const unsigned felt = mymsh->el->GetElementFaceType(kel, jface);  		  		  
		for(unsigned i=0; i<nve; i++) {
		  unsigned inode=mymsh->el->GetFaceVertexIndex(kel,jface,i)-1u;
		  unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
		  unsigned int ilocal = mymsh->el->GetLocalFaceVertexIndex(kel, jface, i);
		  for(unsigned idim=0; idim<dim; idim++) {
		    vx_face[idim][i]=(*mymsh->_coordinate->_Sol[idim])(inode_Metis) + Soli[indexVAR[idim]][ilocal];
		  }
		}
		for(unsigned igs=0; igs < mymsh->_finiteElement[felt][SolType2]->GetGaussPointNumber(); igs++) {
		  mymsh->_finiteElement[felt][SolType2]->JacobianSur(vx_face,igs,Weight,phi,gradphi,normal);
		  //phi1 =mymsh->_finiteElement[felt][SolType2]->GetPhi(igs);
		  // *** phi_i loop ***
		  for(unsigned i=0; i<nve; i++) {
		    adept::adouble value = - phi[i]*tau/rhof*Weight;
		    unsigned int ilocal = mymsh->el->GetLocalFaceVertexIndex(kel, jface, i);
		    
		    for(unsigned idim=0; idim<dim; idim++) {
		      if((!solidmark[ilocal])){
			aRhs[indexVAR[dim+idim]][ilocal]   += value*normal[idim];
		      }
		      else { //if interface node it goes to solid
			aRhs[indexVAR[idim]][ilocal]   += value*normal[idim];
		      }
		    }	    
		  }
		}
	      }
	    }
	  }    
	}
	  	  
	// *** Gauss point loop ***
	double area=1.;
	adept::adouble supg_tau;
	for (unsigned ig=0;ig < mymsh->_finiteElement[kelt][SolType2]->GetGaussPointNumber(); ig++) {
	  // *** get Jacobian and test function and test function derivatives in the moving frame***
	  mymsh->_finiteElement[kelt][SolType2]->Jacobian(vx,ig,Weight,phi,gradphi,nablaphi);
	  mymsh->_finiteElement[kelt][SolType2]->Jacobian(vx_hat,ig,Weight_hat,phi_hat,gradphi_hat,nablaphi_hat);
	  phi1=mymsh->_finiteElement[kelt][SolType1]->GetPhi(ig);
	  
	  if (flag_mat==2 || iel == mymsh->IS_Mts2Gmt_elem_offset[iproc]) {
	    if(ig==0){
	      double GaussWeight = mymsh->_finiteElement[kelt][SolType2]->GetGaussWeight(ig);
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
	      GradSolhatVAR[i][j]=0.;
	    }
	    for(int j=0; j<3*(dim-1); j++) {
	      NablaSolVAR[i][j]=0.;
	      NablaSolhatVAR[i][j]=0.;
	    }
	    for (unsigned inode=0; inode<nve; inode++) {
	      SolVAR[i]+=phi[inode]*Soli[indexVAR[i]][inode];
	      for(int j=0; j<dim; j++) {
		GradSolVAR[i][j]+=gradphi[inode*dim+j]*Soli[indexVAR[i]][inode];
		GradSolhatVAR[i][j] += gradphi_hat[inode*dim+j]*Soli[indexVAR[i]][inode];
	      }	 
	      for(int j=0; j<3*(dim-1); j++) {
		NablaSolVAR[i][j]+=nablaphi[inode*3*(dim-1)+j]*Soli[indexVAR[i]][inode];
		NablaSolhatVAR[i][j]+=nablaphi_hat[inode*3*(dim-1)+j]*Soli[indexVAR[i]][inode];
	      }
	    }
	  } 
	  // pressure
	  SolVAR[2*dim]=0.;
	  for (unsigned inode=0; inode<nve1; inode++) {
	    adept::adouble soli = Soli[indexVAR[2*dim]][inode];
	    SolVAR[2*dim]+=phi1[inode]*soli;
	  }
	  // ---------------------------------------------------------------------------
	  //BEGIN FLUID ASSEMBLY ============
	  if(flag_mat==2){
	    //BEGIN ALE + Momentum (Navier-Stokes)
	    { 
	      for (unsigned i=0; i<nve; i++){
		
		//BEGIN redidual Laplacian ALE map in the reference domain
		adept::adouble LapmapVAR[3] = {0., 0., 0.};
		for(int idim=0; idim<dim; idim++) {
		  for(int jdim=0; jdim<dim; jdim++) {
		    LapmapVAR[idim] += (GradSolhatVAR[idim][jdim]*gradphi_hat[i*dim+jdim]) ;
		  }
		}
		for(int idim=0; idim<dim; idim++) {
		  aRhs[indexVAR[idim]][i]+=(!solidmark[i])*(-LapmapVAR[idim]*Weight_nojac);
		}
		//END redidual Laplacian ALE map in reference domain

		//BEGIN redidual Navier-Stokes in moving domain   
		adept::adouble LapvelVAR[3]={0.,0.,0.};
		adept::adouble AdvaleVAR[3]={0.,0.,0.};		
		
		for(int idim=0.; idim<dim; idim++) {
		  for(int jdim=0.; jdim<dim; jdim++) {
		    LapvelVAR[idim]     += GradSolVAR[dim+idim][jdim]*gradphi[i*dim+jdim];
		    AdvaleVAR[idim]	+= SolVAR[dim+jdim]*GradSolVAR[dim+idim][jdim]*phi[i];   
		  }
		}
		
		for(int idim=0; idim<dim; idim++) {
		  adept::adouble value = (-AdvaleVAR[idim]      	     // advection term	
					  -IRe*LapvelVAR[idim]	   	     // viscous dissipation
					  +SolVAR[2*dim]*gradphi[i*dim+idim] // pressure gradient
					  )*Weight;
		  if((!solidmark[i])){
		    aRhs[indexVAR[dim+idim]][i]+=value; 
		  }
		  else{
		    aRhs[indexVAR[idim]][i]+= value;
		  }
		  //END redidual Navier-Stokes in moving domain    
		}
	      } 
	    } 
	    //END ALE + Momentum (Navier-Stokes)      
	    
	    //BEGIN continuity block 
	    {  	    
	      adept::adouble div_vel=0.;
	      for(int i=0; i<dim; i++) {
		div_vel+=GradSolVAR[dim+i][i];
	      }
	      for (unsigned i=0; i<nve1; i++) {
		aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*div_vel)*Weight;
	      }
	    }
	    //END continuity block ===========================
	  }   
	  //END FLUID ASSEMBLY ============
	    
	  //*******************************************************************************************************
	    
	  //BEGIN SOLID ASSEMBLY ============
	  else{	      
	    //BEGIN build Chauchy Stress in moving domain
	    //physical quantity
	    adept::adouble J_hat;
	    adept::adouble I_e;
	    adept::adouble Cauchy[3][3];
	    adept::adouble Id2th[3][3]= {{ 1.,0.,0.}, { 0.,1.,0.}, { 0.,0.,1.}};
	    
	    adept::adouble I1_B=0.;
	    adept::adouble I2_B=0.;
	    
	    if (solid_model==0) { // Saint-Venant
	      adept::adouble e[3][3];
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
		  Cauchy[i][j] = 2*mus*e[i][j]-2*mus*I_e*SolVAR[2*dim]*Id2th[i][j];
		  //+(penalty)*lambda*I_e*Id2th[i][j];
		}
	      }
	    }
	  
	    else { // hyperelastic non linear material
	      adept::adouble F[3][3]={{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
	      adept::adouble B[3][3];      
	    
	      for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
		  F[i][j]+=GradSolhatVAR[i][j];
		}
	      }
		
	      J_hat =   F[0][0]*F[1][1]*F[2][2] + F[0][1]*F[1][2]*F[2][0] + F[0][2]*F[1][0]*F[2][1]
		      - F[2][0]*F[1][1]*F[0][2] - F[2][1]*F[1][2]*F[0][0] - F[2][2]*F[1][0]*F[0][1];
	      for (int I=0; I<3; ++I) {
		for (int J=0; J<3; ++J) {
		  B[I][J]=0.;
		  for (int K=0; K<3; ++K) {
		    //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
		    B[I][J] += F[I][K]*F[J][K];
		  }
		}
	      }	  
	      if( solid_model <=4 ){ // Neo-Hookean
		I1_B = B[0][0] + B[1][1] + B[2][2];
		
	    	for (int I=0; I<3; ++I) {
		  for (int J=0; J<3; ++J) {
		    if	    ( 1 == solid_model ) Cauchy[I][J] = mus*B[I][J] 
							       -mus*I1_B*SolVAR[2*dim]*Id2th[I][J]; 	//Wood-Bonet J_hat  =1;
		    else if ( 2 == solid_model ) Cauchy[I][J] = mus/J_hat*B[I][J] 
							       -mus/J_hat*SolVAR[2*dim]*Id2th[I][J];    //Wood-Bonet J_hat !=1;
		    else if ( 3 == solid_model ) Cauchy[I][J] = mus*(B[I][J] - Id2th[I][J])/J_hat
							      + lambda/J_hat*log(J_hat)*Id2th[I][J]; 	//Wood-Bonet penalty
		    else if ( 4 == solid_model ) Cauchy[I][J] = mus*(B[I][J] - I1_B*Id2th[I][J]/3.)/pow(J_hat,5./3.)
						              + lambda*(J_hat-1.)*Id2th[I][J];  	  //Allan-Bower
		    
		  }
		}
	      }
	      else if ( 5 == solid_model ){ //Mooney-Rivlin
		adept::adouble detB=   B[0][0]* ( B[1][1]*B[2][2] - B[2][1]*B[1][2] ) 
				     - B[0][1]* ( B[2][2]*B[1][0] - B[1][2]*B[2][0] ) 
				     + B[0][2]* ( B[1][0]*B[2][1] - B[2][0]*B[1][1] );
		adept::adouble invdetB=1./detB;		      
		adept::adouble invB[3][3];
		
		invB[0][0] =  (B[1][1]*B[2][2]-B[2][1]*B[1][2])*invdetB;
		invB[1][0] = -(B[0][1]*B[2][2]-B[0][2]*B[2][1])*invdetB;
		invB[2][0] =  (B[0][1]*B[1][2]-B[0][2]*B[1][1])*invdetB;
		invB[0][1] = -(B[1][0]*B[2][2]-B[1][2]*B[2][0])*invdetB;
		invB[1][1] =  (B[0][0]*B[2][2]-B[0][2]*B[2][0])*invdetB;
		invB[2][1] = -(B[0][0]*B[1][2]-B[1][0]*B[0][2])*invdetB;
		invB[0][2] =  (B[1][0]*B[2][1]-B[2][0]*B[1][1])*invdetB;
		invB[1][2] = -(B[0][0]*B[2][1]-B[2][0]*B[0][1])*invdetB;
		invB[2][2] =  (B[0][0]*B[1][1]-B[1][0]*B[0][1])*invdetB;
		
		I1_B = B[0][0] + B[1][1] + B[2][2];
		I2_B = B[0][0] * B[1][1] + B[1][1] * B[2][2] + B[2][2]*B[0][0]
		      -B[0][1] * B[1][0] - B[1][2] * B[2][1] - B[2][0]*B[0][2];
		
		for (int I=0; I<3; ++I) {
		  for (int J=0; J<3; ++J) {
		    Cauchy[I][J] = mus*(2.*B[I][J] - invB[I][J] )/3.
				  -mus*(2.*I1_B - I2_B )/3.*SolVAR[2*dim]*Id2th[I][J];
		    ;
		  }
		}
		
	      }	  		
	    }
	    //END build Chauchy Stress in moving domain
	    
	    //BEGIN v=0 + Momentum (Solid)      
	    {
	      for (unsigned i=0; i<nve; i++) {

		//BEGIN redidual v=0 in fixed domain 
		for(int idim=0; idim<dim; idim++) {
		  aRhs[indexVAR[dim+idim]][i] += (-phi[i]*(-SolVAR[dim+idim]))*Weight_hat;
		}
                //END redidual v=0 in fixed domain 
                
                //BEGIN redidual Solid Momentum in moving domain 
		adept::adouble CauchyDIR[3]={0.,0.,0.};
		for(int idim=0.; idim<dim; idim++) {
		  for(int jdim=0.; jdim<dim; jdim++) {
		    CauchyDIR[idim]+= gradphi[i*dim+jdim]*Cauchy[idim][jdim];
		  }
		}
		for(int idim=0; idim<dim; idim++) {
		  aRhs[indexVAR[idim]][i] += (phi[i]*_gravity[idim] - CauchyDIR[idim])*Weight;			
		} 
		//END redidual Solid Momentum in moving domain 
	      }
	    }
	    //END v=0 + Momentum (Solid)   
	    
	    //BEGIN continuity block
	    {       
	        for (unsigned i=0; i<nve1; i++) {
		if(!penalty){		    
		  if (0 == solid_model) {
		    aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*(I_e + (!incompressible)/lambda*SolVAR[2*dim] ) )*Weight_hat;
		  }
		  else if (1 == solid_model || 5 == solid_model) {
		    aRhs[indexVAR[2*dim]][i] +=phi1[i]*(J_hat-1. + (!incompressible)/lambda*SolVAR[2*dim])*Weight_hat;
		  }
		  else if (2 == solid_model){
		    aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*( log(J_hat)/J_hat + (!incompressible)/lambda*SolVAR[2*dim] ) )*Weight_hat;
		  }
		}
		else if (3 == solid_model || 4 == solid_model){ // pressure = 0 in the solid
		  aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*( SolVAR[2*dim] ) )*Weight_hat;
		}
	      }
	    }
	    //END continuity block
	  }  
	  //END SOLID ASSEMBLY ============
	}
      }
	
      //BEGIN local to global assembly 	
      //copy adouble aRhs into double Rhs
      for (unsigned i=0;i<2*dim;i++) {
	for(int j=0; j<nve; j++) {
	  Rhs[indexVAR[i]][j] = aRhs[indexVAR[i]][j].value();
	}
      }
      for (unsigned j=0;j<nve1;j++) {
	Rhs[indexVAR[2*dim]][j] = aRhs[indexVAR[2*dim]][j].value();
      }	
      for(int i=0; i<2*dim+1; i++) {
	myRES->add_vector_blocked(Rhs[indexVAR[i]],dofsVAR[i]);
      }     

      //Store equations
      for(int i=0; i<2*dim; i++) {  
	s.dependent(&aRhs[indexVAR[i]][0], nve);
	s.independent(&Soli[indexVAR[i]][0], nve); 
      }
      s.dependent(&aRhs[indexVAR[2*dim]][0], nve1);
      s.independent(&Soli[indexVAR[2*dim]][0], nve1);   
      s.jacobian(&Jac[0]);	
      unsigned nveAll=(2*dim*nve+nve1);
      for (int inode=0;inode<nveAll;inode++){
	for (int jnode=0;jnode<nveAll;jnode++){
	   KKloc[inode*nveAll+jnode]=-Jac[jnode*nveAll+inode];
	}
      }
      myKK->add_matrix_blocked(KKloc,dofsAll,dofsAll);
      s.clear_independents();
      s.clear_dependents();
       
      //END local to global assembly
   
    } //end list of elements loop
   
    myKK->close();
    myRES->close();
  
    delete area_elem_first; 
	  
    // *************************************
    end_time=clock();
    AssemblyTime+=(end_time-start_time);
    // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

  }
  
  /*
  
  void IncompressibleFSIAssemblyAD(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assemble_matrix) {
       
    clock_t AssemblyTime=0;
    clock_t start_time, end_time;
  
    static adept::Stack s; 
    
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
    vector<adept::adouble> SolVAR(2*dim+1);
    vector<vector<adept::adouble> > GradSolVAR(2*dim);
    for(int i=0;i<2*dim;i++){
      GradSolVAR[i].resize(dim);
    }
    vector<vector<adept::adouble> > GradSolhatVAR(dim);
    for(int i=0;i<dim;i++){
      GradSolhatVAR[i].resize(dim);
    }
  
    vector <bool> solidmark;
    vector <double > phi_hat;
    vector <double > phi;
    vector <adept::adouble> gradphi;
    vector <adept::adouble> nablaphi;
    vector <double> gradphi_hat;
    vector <double> nablaphi_hat;
    
    solidmark.reserve(max_size);
    phi_hat.reserve(max_size);
    phi.reserve(max_size);
    gradphi.reserve(max_size*dim);
    nablaphi.reserve(max_size*3*(dim-1));
    gradphi_hat.reserve(max_size*dim);
    nablaphi_hat.reserve(max_size*3*(dim-1));
    
    const double *phi1;
    
    adept::adouble Weight=0.;
    double Weight_nojac=0.;
    double Weight_hat=0.;
  
    vector <vector < adept::adouble> > vx(dim);
    vector <vector < adept::adouble> > vx_face(dim);
    vector <vector < double> > vx_hat(dim);
  
    for(int i=0;i<dim;i++){
      vx[i].reserve(max_size);
      vx_face[i].resize(9);
      vx_hat[i].reserve(max_size);
    }
   
    vector< vector< adept::adouble > > Soli(2*dim+1);
    vector< vector< int > > dofsVAR(2*dim+1); 
    for(int i=0;i<2*dim+1;i++){
      Soli[i].reserve(max_size);
      dofsVAR[i].reserve(max_size);
    }
    
    vector< vector< double > > Rhs(2*dim+1);
    vector< vector< adept::adouble > > aRhs(2*dim+1);
    for(int i=0;i<2*dim+1;i++){
      aRhs[i].reserve(max_size);
      Rhs[i].reserve(max_size);
    }     
    
    vector < int > dofsVel;
    vector < int > dofsAll;
    dofsAll.reserve(max_size*(2+dim+1));
    dofsVel.reserve(max_size*dim);
    
    vector < double > B11;
    vector < double > B22;
    vector < double > B12;
    vector < double > Bmass;
    Bmass.reserve(max_size*max_size);
    B11.reserve(max_size*max_size);
    B22.reserve(max_size*max_size);
    B12.reserve(max_size*max_size);
    vector < double > Bm;
    Bm.reserve(max_size*dim*max_size*(2*dim+1));
        
    vector < double > Jac(max_size*dim*max_size*(2*dim+1));
    
    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof	 	= ml_prob.parameters.get<Fluid>("Fluid").get_density();             
    double mu_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_shear_modulus(); 
    double lambda_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_lambda();        
    double mus		= mu_lame/rhof;
    double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();   
    double lambda	= lambda_lame / rhof;
    double betans	= 1.;
    int    solid_model	= ml_prob.parameters.get<Solid>("Solid").get_physical_model();     
       
    bool incompressible=( 0.5 == ml_prob.parameters.get<Solid>("Solid").get_poisson_coeff() )?1:0;
    const bool penalty = ml_prob.parameters.get<Solid>("Solid").get_if_penalty();
    const bool mass_penalty = ml_prob.parameters.get<Solid>("Solid").get_if_mass_penalty();
        
    // gravity
    double _gravity[3]={0.,0.,0.};
     
    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));  
    unsigned end_ind2   = mymsh->GetEndIndex(SolType2);

    unsigned SolType1 = ml_sol->GetSolutionType(ml_sol->GetIndex("P"));  
    unsigned end_ind1   = mymsh->GetEndIndex(SolType1);

    // mesh and procs
    unsigned nel    = mymsh->GetElementNumber();
    unsigned igrid  = mymsh->GetLevel();
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
    
    // *** element loop ***
    for(int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

      unsigned kel        = mymsh->IS_Mts2Gmt_elem[iel]; 
      short unsigned kelt = myel->GetElementType(kel);
      unsigned nve        = myel->GetElementDofNumber(kel,end_ind2);
      unsigned nve1       = myel->GetElementDofNumber(kel,end_ind1);
      int flag_mat        = myel->GetElementMaterial(kel);

      // *******************************************************************************************************
    
      //initialization of everything is in common fluid and solid
    
      //Rhs
      for(int i=0; i<2*dim; i++) {
	dofsVAR[i].resize(nve);
	Soli[indexVAR[i]].resize(nve);
	aRhs[indexVAR[i]].resize(nve);
	Rhs[indexVAR[i]].resize(nve);
      }
      dofsVAR[2*dim].resize(nve1);
      Soli[indexVAR[2*dim]].resize(nve1);
      aRhs[indexVAR[2*dim]].resize(nve1);
      Rhs[indexVAR[2*dim]].resize(nve1);
      
      dofsAll.resize(0);
      dofsVel.resize(0);
      
      B11.resize(nve1*nve1);
      B22.resize(nve*nve);
      B12.resize(nve1*nve);
      Bm.resize(nve*dim*(2*dim*nve+nve1));
    
      Bmass.resize(nve*nve);
      memset(&Bmass[0],0,nve*nve*sizeof(double));
      
      // ----------------------------------------------------------------------------------------
      // coordinates, solutions, displacement, velocity dofs
    
      solidmark.resize(nve);
      phi_hat.resize(nve);
      phi.resize(nve);
      gradphi.resize(nve*dim);
      nablaphi.resize(nve*3*(dim-1));
      gradphi_hat.resize(nve*dim);
      nablaphi_hat.resize(nve*3*(dim-1));
      
      
      for(int i=0;i<dim;i++){
	vx[i].resize(nve);
	vx_hat[i].resize(nve);
      }
    
      for (unsigned i=0;i<nve;i++) {
	unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;
	unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
	// flag to know if the node "inode" lays on the fluid-solid interface
	solidmark[i]=myel->GetNodeRegion(inode); // to check
	for(int j=0; j<dim; j++) {
	  Soli[indexVAR[j]][i]     =  (*mysolution->_Sol[indVAR[j]])(inode_Metis);
	  Soli[indexVAR[j+dim]][i] =  (*mysolution->_Sol[indVAR[j+dim]])(inode_Metis);
	    	    
	  aRhs[indexVAR[j]][i]     = 0.;
	  aRhs[indexVAR[j+dim]][i] = 0.;
	   
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
	unsigned inode=(SolType1<3)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
	unsigned inode_Metis =mymsh->GetMetisDof(inode,SolType[2*dim]);
	dofsVAR[2*dim][i]=myLinEqSolver->GetKKDof(indVAR[2*dim],indexVAR[2*dim],inode);
	Soli[indexVAR[2*dim]][i] = (*mysolution->_Sol[indVAR[2*dim]])(inode_Metis);
	aRhs[indexVAR[2*dim]][i] = 0.;
      }
      
      // build dof ccomposition             
      for(int idim=0;idim<dim;idim++){
	dofsAll.insert( dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end() );
	dofsVel.insert( dofsVel.end(), dofsVAR[idim+dim].begin(), dofsVAR[idim+dim].end() );
      }
      dofsAll.insert( dofsAll.end(), dofsVel.begin(), dofsVel.end() );
      dofsAll.insert( dofsAll.end(), dofsVAR[2*dim].begin(), dofsVAR[2*dim].end() );
 
      if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {  
	
	s.new_recording();
	
	for (unsigned idim=0;idim<dim;idim++) {
	  for(int j=0; j<nve; j++) {
	    vx[idim][j]= vx_hat[idim][j] + Soli[indexVAR[idim]][j];
	  }
	}
	
	// Boundary integral
	{
	  double tau=0.;
	  vector<adept::adouble> normal(dim,0);
	       
	  // loop on faces
	  for(unsigned jface=0; jface<myel->GetElementFaceNumber(kel); jface++) {
		
	    // look for boundary faces
	    if(myel->GetFaceElementIndex(kel,jface)<0) {
	      unsigned int face = -(mymsh->el->GetFaceElementIndex(kel,jface)+1);	      
	      if( !ml_sol->_SetBoundaryConditionFunction(0.,0.,0.,"U",tau,face,0.) && tau!=0.){
		unsigned nve = mymsh->el->GetElementFaceDofNumber(kel,jface,SolType2);
		const unsigned felt = mymsh->el->GetElementFaceType(kel, jface); 
		  		  		  
		for(unsigned i=0; i<nve; i++) {
		  unsigned inode=mymsh->el->GetFaceVertexIndex(kel,jface,i)-1u;
		  unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
		  unsigned int ilocal = mymsh->el->GetLocalFaceVertexIndex(kel, jface, i);
		  for(unsigned idim=0; idim<dim; idim++) {
		    vx_face[idim][i]=(*mymsh->_coordinate->_Sol[idim])(inode_Metis) + Soli[indexVAR[idim]][ilocal];
		  }
		}
		for(unsigned igs=0; igs < ml_prob._ml_msh->_finiteElement[felt][SolType2]->GetGaussPointNumber(); igs++) {
		  (ml_prob._ml_msh->_finiteElement[felt][SolType2]->*(ml_prob._ml_msh->_finiteElement[felt][SolType2])->Jacobian_sur_AD_ptr)(vx_face,igs,Weight,gradphi,normal);
		  phi1 =ml_prob._ml_msh->_finiteElement[felt][SolType2]->GetPhi(igs);
		  // *** phi_i loop ***
		  for(unsigned i=0; i<nve; i++) {
		    adept::adouble value = - phi1[i]*tau/rhof*Weight;
		    unsigned int ilocal = mymsh->el->GetLocalFaceVertexIndex(kel, jface, i);
		    
		    for(unsigned idim=0; idim<dim; idim++) {
		      aRhs[indexVAR[dim+idim]][ilocal]   += value*normal[idim];
		    }
		    
		    
		    
		  }
		}
	      }
	    }
	  }    
	}
	  	  
	// *** Gauss point loop ***
	double area=1.;
	for (unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][SolType2]->GetGaussPointNumber(); ig++) {
	  // *** get Jacobian and test function and test function derivatives in the moving frame***
	  ml_prob._ml_msh->_finiteElement[kelt][SolType2]->Jacobian_AD(vx,ig,Weight,phi,gradphi,nablaphi);
	  ml_prob._ml_msh->_finiteElement[kelt][SolType2]->Jacobian(vx_hat,ig,Weight_hat,phi_hat,gradphi_hat,nablaphi_hat);
	  //phi =ml_prob._ml_msh->_finiteElement[kelt][SolType2]->GetPhi(ig);
	  phi1=ml_prob._ml_msh->_finiteElement[kelt][SolType1]->GetPhi(ig);
	  
	  if (flag_mat==2 || iel == mymsh->IS_Mts2Gmt_elem_offset[iproc]) {
	    if(ig==0){
	      double GaussWeight = ml_prob._ml_msh->_finiteElement[kelt][SolType2]->GetGaussWeight(ig);
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
	      SolVAR[i]+=phi[inode]*Soli[indexVAR[i]][inode];
	      for(int j=0; j<dim; j++) {
		GradSolVAR[i][j]+=gradphi[inode*dim+j]*Soli[indexVAR[i]][inode];
		if(i<dim){ 
		  GradSolhatVAR[i][j] += gradphi_hat[inode*dim+j]*Soli[indexVAR[i]][inode];
		}
	      }	      
	    }
	  } 
	  // pressure
	  SolVAR[2*dim]=0.;
	  for (unsigned inode=0; inode<nve1; inode++) {
	    adept::adouble soli = Soli[indexVAR[2*dim]][inode];
	    SolVAR[2*dim]+=phi1[inode]*soli;
	  }
	  // ---------------------------------------------------------------------------
	  //BEGIN FLUID ASSEMBLY ============
	  if(flag_mat==2){
	    //BEGIN ALE + Momentum (Navier-Stokes)
	    { 
	      for (unsigned i=0; i<nve; i++){
		
		//BEGIN redidual Laplacian ALE map in the reference domain
		adept::adouble LapmapVAR[3] = {0., 0., 0.};
		for(int idim=0; idim<dim; idim++) {
		  for(int jdim=0; jdim<dim; jdim++) {
		    LapmapVAR[idim] += (GradSolhatVAR[idim][jdim]*gradphi_hat[i*dim+jdim]) ;
		  }
		}
		for(int idim=0; idim<dim; idim++) {
		  aRhs[indexVAR[idim]][i]+=(!solidmark[i])*(-LapmapVAR[idim]*Weight_nojac);
		}
		//END redidual Laplacian ALE map in reference domain

		//BEGIN redidual Navier-Stokes in moving domain   
		adept::adouble LapvelVAR[3]={0.,0.,0.};
		adept::adouble AdvaleVAR[3]={0.,0.,0.};
		for(int idim=0.; idim<dim; idim++) {
		  for(int jdim=0.; jdim<dim; jdim++) {
		    LapvelVAR[idim]+=GradSolVAR[dim+idim][jdim]*gradphi[i*dim+jdim];
		    AdvaleVAR[idim]+=SolVAR[dim+jdim]*GradSolVAR[dim+idim][jdim]*phi[i];
		  }
		}
		for(int idim=0; idim<dim; idim++) {
		  aRhs[indexVAR[dim+idim]][i]+= (-AdvaleVAR[idim]      	      	      // advection term	
						 -IRe*LapvelVAR[idim]	 	      // viscous dissipation
						 +SolVAR[2*dim]*gradphi[i*dim+idim]   // pressure gradient
						 )*Weight;
		  //END redidual Navier-Stokes in moving domain    
		}
	      } 
	    } 
	    //END ALE + Momentum (Navier-Stokes)      
	    
	    //BEGIN continuity block 
	    {  	    
	      adept::adouble div_vel=0.;
	      for(int i=0; i<dim; i++) {
		div_vel+=GradSolVAR[dim+i][i];
	      }
	      for (unsigned i=0; i<nve1; i++) {
		aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*div_vel)*Weight;
	      }
	    }
	    //END continuity block ===========================
	  }   
	  //END FLUID ASSEMBLY ============
	    
	  //*******************************************************************************************************
	    
	  //BEGIN SOLID ASSEMBLY ============
	  else{	      
	    //BEGIN build Chauchy Stress in moving domain
	    //physical quantity
	    adept::adouble J_hat;
	    adept::adouble I_e;
	    adept::adouble Cauchy[3][3];
	    adept::adouble Id2th[3][3]= {{ 1.,0.,0.}, { 0.,1.,0.}, { 0.,0.,1.}};
	    
	    adept::adouble I1_B=0.;
	    adept::adouble I2_B=0.;
	    
	    if (solid_model==0) { // Saint-Venant
	      adept::adouble e[3][3];
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
		  Cauchy[i][j] = 2*mus*e[i][j]-2*mus*I_e*SolVAR[2*dim]*Id2th[i][j];
		  //+(penalty)*lambda*I_e*Id2th[i][j];
		}
	      }
	    }
	  
	    else { // hyperelastic non linear material
	      adept::adouble F[3][3]={{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
	      adept::adouble B[3][3];      
	    
	      for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
		  F[i][j]+=GradSolhatVAR[i][j];
		}
	      }
		
	      J_hat =   F[0][0]*F[1][1]*F[2][2] + F[0][1]*F[1][2]*F[2][0] + F[0][2]*F[1][0]*F[2][1]
		      - F[2][0]*F[1][1]*F[0][2] - F[2][1]*F[1][2]*F[0][0] - F[2][2]*F[1][0]*F[0][1];
	      for (int I=0; I<3; ++I) {
		for (int J=0; J<3; ++J) {
		  B[I][J]=0.;
		  for (int K=0; K<3; ++K) {
		    //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
		    B[I][J] += F[I][K]*F[J][K];
		  }
		}
	      }	  
	      if( solid_model <=4 ){ // Neo-Hookean
		I1_B = B[0][0] + B[1][1] + B[2][2];
		
	    	for (int I=0; I<3; ++I) {
		  for (int J=0; J<3; ++J) {
		    if	    ( 1 == solid_model ) Cauchy[I][J] = mus*B[I][J] 
							       -mus*I1_B*SolVAR[2*dim]*Id2th[I][J]; 	//Wood-Bonet J_hat  =1;
		    else if ( 2 == solid_model ) Cauchy[I][J] = mus/J_hat*B[I][J] 
							       -mus/J_hat*SolVAR[2*dim]*Id2th[I][J];    //Wood-Bonet J_hat !=1;
		    else if ( 3 == solid_model ) Cauchy[I][J] = mus*(B[I][J] - Id2th[I][J])/J_hat
							      + lambda/J_hat*log(J_hat)*Id2th[I][J]; 	//Wood-Bonet penalty
		    else if ( 4 == solid_model ) Cauchy[I][J] = mus*(B[I][J] - I1_B*Id2th[I][J]/3.)/pow(J_hat,5./3.)
						              + lambda*(J_hat-1.)*Id2th[I][J];  	  //Allan-Bower
		    
		  }
		}
	      }
	      else if ( 5 == solid_model ){ //Mooney-Rivlin
		adept::adouble detB=   B[0][0]* ( B[1][1]*B[2][2] - B[2][1]*B[1][2] ) 
				     - B[0][1]* ( B[2][2]*B[1][0] - B[1][2]*B[2][0] ) 
				     + B[0][2]* ( B[1][0]*B[2][1] - B[2][0]*B[1][1] );
		adept::adouble invdetB=1./detB;		      
		adept::adouble invB[3][3];
		
		invB[0][0] =  (B[1][1]*B[2][2]-B[2][1]*B[1][2])*invdetB;
		invB[1][0] = -(B[0][1]*B[2][2]-B[0][2]*B[2][1])*invdetB;
		invB[2][0] =  (B[0][1]*B[1][2]-B[0][2]*B[1][1])*invdetB;
		invB[0][1] = -(B[1][0]*B[2][2]-B[1][2]*B[2][0])*invdetB;
		invB[1][1] =  (B[0][0]*B[2][2]-B[0][2]*B[2][0])*invdetB;
		invB[2][1] = -(B[0][0]*B[1][2]-B[1][0]*B[0][2])*invdetB;
		invB[0][2] =  (B[1][0]*B[2][1]-B[2][0]*B[1][1])*invdetB;
		invB[1][2] = -(B[0][0]*B[2][1]-B[2][0]*B[0][1])*invdetB;
		invB[2][2] =  (B[0][0]*B[1][1]-B[1][0]*B[0][1])*invdetB;
		
		I1_B = B[0][0] + B[1][1] + B[2][2];
		I2_B = B[0][0] * B[1][1] + B[1][1] * B[2][2] + B[2][2]*B[0][0]
		      -B[0][1] * B[1][0] - B[1][2] * B[2][1] - B[2][0]*B[0][2];
		
		for (int I=0; I<3; ++I) {
		  for (int J=0; J<3; ++J) {
		    Cauchy[I][J] = mus*(2.*B[I][J] - invB[I][J] )/3.
				  -mus*(2.*I1_B - I2_B )/3.*SolVAR[2*dim]*Id2th[I][J];
		    ;
		  }
		}
		
	      }	  		
	    }
	    //END build Chauchy Stress in moving domain
	    
	    //BEGIN v=0 + Momentum (Solid)      
	    {
	      for (unsigned i=0; i<nve; i++) {

		//BEGIN redidual v=0 in fixed domain 
		for(int idim=0; idim<dim; idim++) {
		  aRhs[indexVAR[idim]][i] += (-phi[i]*(-SolVAR[dim+idim]))*Weight_hat;
		}
		for (unsigned j=0; j<nve; j++) { //mass matrix preconditioner in the ALE equation
		   Bmass[i*nve+j] += mass_penalty*phi[i]*phi[j]*Weight_hat;
		}
                //END redidual v=0 in fixed domain 
                
                //BEGIN redidual Solid Momentum in moving domain 
		adept::adouble CauchyDIR[3]={0.,0.,0.};
		for(int idim=0.; idim<dim; idim++) {
		  for(int jdim=0.; jdim<dim; jdim++) {
		    CauchyDIR[idim]+= gradvoid IncompressibleFSIAssemblyAD(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assemble_matrix) {
       
    clock_t AssemblyTime=0;
    clock_t start_time, end_time;
  
    static adept::Stack s; 
    
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
    vector<adept::adouble> SolVAR(2*dim+1);
    vector<vector<adept::adouble> > GradSolVAR(2*dim);
    for(int i=0;i<2*dim;i++){
      GradSolVAR[i].resize(dim);
    }
    vector<vector<adept::adouble> > GradSolhatVAR(dim);
    for(int i=0;i<dim;i++){
      GradSolhatVAR[i].resize(dim);
    }
  
    vector <bool> solidmark;
    vector <double > phi_hat;
    vector <double > phi;
    vector <adept::adouble> gradphi;
    vector <adept::adouble> nablaphi;
    vector <double> gradphi_hat;
    vector <double> nablaphi_hat;
    
    solidmark.reserve(max_size);
    phi_hat.reserve(max_size);
    phi.reserve(max_size);
    gradphi.reserve(max_size*dim);
    nablaphi.reserve(max_size*3*(dim-1));
    gradphi_hat.reserve(max_size*dim);
    nablaphi_hat.reserve(max_size*3*(dim-1));
    
    const double *phi1;
    
    adept::adouble Weight=0.;
    double Weight_nojac=0.;
    double Weight_hat=0.;
  
    vector <vector < adept::adouble> > vx(dim);
    vector <vector < adept::adouble> > vx_face(dim);
    vector <vector < double> > vx_hat(dim);
  
    for(int i=0;i<dim;i++){
      vx[i].reserve(max_size);
      vx_face[i].resize(9);
      vx_hat[i].reserve(max_size);
    }
   
    vector< vector< adept::adouble > > Soli(2*dim+1);
    vector< vector< int > > dofsVAR(2*dim+1); 
    for(int i=0;i<2*dim+1;i++){
      Soli[i].reserve(max_size);
      dofsVAR[i].reserve(max_size);
    }
    
    vector< vector< double > > Rhs(2*dim+1);
    vector< vector< adept::adouble > > aRhs(2*dim+1);
    for(int i=0;i<2*dim+1;i++){
      aRhs[i].reserve(max_size);
      Rhs[i].reserve(max_size);
    }     
    
    vector < int > dofsVel;
    vector < int > dofsAll;
    dofsAll.reserve(max_size*(2+dim+1));
    dofsVel.reserve(max_size*dim);
    
    vector < double > B11;
    vector < double > B22;
    vector < double > B12;
    vector < double > Bmass;
    Bmass.reserve(max_size*max_size);
    B11.reserve(max_size*max_size);
    B22.reserve(max_size*max_size);
    B12.reserve(max_size*max_size);
    vector < double > Bm;
    Bm.reserve(max_size*dim*max_size*(2*dim+1));
        
    vector < double > Jac(max_size*dim*max_size*(2*dim+1));
    
    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof	 	= ml_prob.parameters.get<Fluid>("Fluid").get_density();             
    double mu_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_shear_modulus(); 
    double lambda_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_lambda();        
    double mus		= mu_lame/rhof;
    double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();   
    double lambda	= lambda_lame / rhof;
    double betans	= 1.;
    int    solid_model	= ml_prob.parameters.get<Solid>("Solid").get_physical_model();     
       
    bool incompressible=( 0.5 == ml_prob.parameters.get<Solid>("Solid").get_poisson_coeff() )?1:0;
    const bool penalty = ml_prob.parameters.get<Solid>("Solid").get_if_penalty();
    const bool mass_penalty = ml_prob.parameters.get<Solid>("Solid").get_if_mass_penalty();
        
    // gravity
    double _gravity[3]={0.,0.,0.};
     
    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));  
    unsigned end_ind2   = mymsh->GetEndIndex(SolType2);

    unsigned SolType1 = ml_sol->GetSolutionType(ml_sol->GetIndex("P"));  
    unsigned end_ind1   = mymsh->GetEndIndex(SolType1);

    // mesh and procs
    unsigned nel    = mymsh->GetElementNumber();
    unsigned igrid  = mymsh->GetLevel();
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
    
    // *** element loop ***
    for(int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

      unsigned kel        = mymsh->IS_Mts2Gmt_elem[iel]; 
      short unsigned kelt = myel->GetElementType(kel);
      unsigned nve        = myel->GetElementDofNumber(kel,end_ind2);
      unsigned nve1       = myel->GetElementDofNumber(kel,end_ind1);
      int flag_mat        = myel->GetElementMaterial(kel);

      // *******************************************************************************************************
    
      //initialization of everything is in common fluid and solid
    
      //Rhs
      for(int i=0; i<2*dim; i++) {
	dofsVAR[i].resize(nve);
	Soli[indexVAR[i]].resize(nve);
	aRhs[indexVAR[i]].resize(nve);
	Rhs[indexVAR[i]].resize(nve);
      }
      dofsVAR[2*dim].resize(nve1);
      Soli[indexVAR[2*dim]].resize(nve1);
      aRhs[indexVAR[2*dim]].resize(nve1);
      Rhs[indexVAR[2*dim]].resize(nve1);
      
      dofsAll.resize(0);
      dofsVel.resize(0);
      
      B11.resize(nve1*nve1);
      B22.resize(nve*nve);
      B12.resize(nve1*nve);
      Bm.resize(nve*dim*(2*dim*nve+nve1));
    
      Bmass.resize(nve*nve);
      memset(&Bmass[0],0,nve*nve*sizeof(double));
      
      // ----------------------------------------------------------------------------------------
      // coordinates, solutions, displacement, velocity dofs
    
      solidmark.resize(nve);
      phi_hat.resize(nve);
      phi.resize(nve);
      gradphi.resize(nve*dim);
      nablaphi.resize(nve*3*(dim-1));
      gradphi_hat.resize(nve*dim);
      nablaphi_hat.resize(nve*3*(dim-1));
      
      
      for(int i=0;i<dim;i++){
	vx[i].resize(nve);
	vx_hat[i].resize(nve);
      }
    
      for (unsigned i=0;i<nve;i++) {
	unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;
	unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
	// flag to know if the node "inode" lays on the fluid-solid interface
	solidmark[i]=myel->GetNodeRegion(inode); // to check
	for(int j=0; j<dim; j++) {
	  Soli[indexVAR[j]][i]     =  (*mysolution->_Sol[indVAR[j]])(inode_Metis);
	  Soli[indexVAR[j+dim]][i] =  (*mysolution->_Sol[indVAR[j+dim]])(inode_Metis);
	    	    
	  aRhs[indexVAR[j]][i]     = 0.;
	  aRhs[indexVAR[j+dim]][i] = 0.;
	   
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
	unsigned inode=(SolType1<3)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
	unsigned inode_Metis =mymsh->GetMetisDof(inode,SolType[2*dim]);
	dofsVAR[2*dim][i]=myLinEqSolver->GetKKDof(indVAR[2*dim],indexVAR[2*dim],inode);
	Soli[indexVAR[2*dim]][i] = (*mysolution->_Sol[indVAR[2*dim]])(inode_Metis);
	aRhs[indexVAR[2*dim]][i] = 0.;
      }
      
      // build dof ccomposition             
      for(int idim=0;idim<dim;idim++){
	dofsAll.insert( dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end() );
	dofsVel.insert( dofsVel.end(), dofsVAR[idim+dim].begin(), dofsVAR[idim+dim].end() );
      }
      dofsAll.insert( dofsAll.end(), dofsVel.begin(), dofsVel.end() );
      dofsAll.insert( dofsAll.end(), dofsVAR[2*dim].begin(), dofsVAR[2*dim].end() );
 
      if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {  
	
	s.new_recording();
	
	for (unsigned idim=0;idim<dim;idim++) {
	  for(int j=0; j<nve; j++) {
	    vx[idim][j]= vx_hat[idim][j] + Soli[indexVAR[idim]][j];
	  }
	}
	
	// Boundary integral
	{
	  double tau=0.;
	  vector<adept::adouble> normal(dim,0);
	       
	  // loop on faces
	  for(unsigned jface=0; jface<myel->GetElementFaceNumber(kel); jface++) {
		
	    // look for boundary faces
	    if(myel->GetFaceElementIndex(kel,jface)<0) {
	      unsigned int face = -(mymsh->el->GetFaceElementIndex(kel,jface)+1);	      
	      if( !ml_sol->_SetBoundaryConditionFunction(0.,0.,0.,"U",tau,face,0.) && tau!=0.){
		unsigned nve = mymsh->el->GetElementFaceDofNumber(kel,jface,SolType2);
		const unsigned felt = mymsh->el->GetElementFaceType(kel, jface); 
		  		  		  
		for(unsigned i=0; i<nve; i++) {
		  unsigned inode=mymsh->el->GetFaceVertexIndex(kel,jface,i)-1u;
		  unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
		  unsigned int ilocal = mymsh->el->GetLocalFaceVertexIndex(kel, jface, i);
		  for(unsigned idim=0; idim<dim; idim++) {
		    vx_face[idim][i]=(*mymsh->_coordinate->_Sol[idim])(inode_Metis) + Soli[indexVAR[idim]][ilocal];
		  }
		}
		for(unsigned igs=0; igs < ml_prob._ml_msh->_finiteElement[felt][SolType2]->GetGaussPointNumber(); igs++) {
		  (ml_prob._ml_msh->_finiteElement[felt][SolType2]->*(ml_prob._ml_msh->_finiteElement[felt][SolType2])->Jacobian_sur_AD_ptr)(vx_face,igs,Weight,gradphi,normal);
		  phi1 =ml_prob._ml_msh->_finiteElement[felt][SolType2]->GetPhi(igs);
		  // *** phi_i loop ***
		  for(unsigned i=0; i<nve; i++) {
		    adept::adouble value = - phi1[i]*tau/rhof*Weight;
		    unsigned int ilocal = mymsh->el->GetLocalFaceVertexIndex(kel, jface, i);
		    
		    for(unsigned idim=0; idim<dim; idim++) {
		      aRhs[indexVAR[dim+idim]][ilocal]   += value*normal[idim];
		    }
		    
		    
		    
		  }
		}
	      }
	    }
	  }    
	}
	  	  
	// *** Gauss point loop ***
	double area=1.;
	for (unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][SolType2]->GetGaussPointNumber(); ig++) {
	  // *** get Jacobian and test function and test function derivatives in the moving frame***
	  ml_prob._ml_msh->_finiteElement[kelt][SolType2]->Jacobian_AD(vx,ig,Weight,phi,gradphi,nablaphi);
	  ml_prob._ml_msh->_finiteElement[kelt][SolType2]->Jacobian(vx_hat,ig,Weight_hat,phi_hat,gradphi_hat,nablaphi_hat);
	  //phi =ml_prob._ml_msh->_finiteElement[kelt][SolType2]->GetPhi(ig);
	  phi1=ml_prob._ml_msh->_finiteElement[kelt][SolType1]->GetPhi(ig);
	  
	  if (flag_mat==2 || iel == mymsh->IS_Mts2Gmt_elem_offset[iproc]) {
	    if(ig==0){
	      double GaussWeight = ml_prob._ml_msh->_finiteElement[kelt][SolType2]->GetGaussWeight(ig);
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
	      SolVAR[i]+=phi[inode]*Soli[indexVAR[i]][inode];
	      for(int j=0; j<dim; j++) {
		GradSolVAR[i][j]+=gradphi[inode*dim+j]*Soli[indexVAR[i]][inode];
		if(i<dim){ 
		  GradSolhatVAR[i][j] += gradphi_hat[inode*dim+j]*Soli[indexVAR[i]][inode];
		}
	      }	      
	    }
	  } 
	  // pressure
	  SolVAR[2*dim]=0.;
	  for (unsigned inode=0; inode<nve1; inode++) {
	    adept::adouble soli = Soli[indexVAR[2*dim]][inode];
	    SolVAR[2*dim]+=phi1[inode]*soli;
	  }
	  // ---------------------------------------------------------------------------
	  //BEGIN FLUID ASSEMBLY ============
	  if(flag_mat==2){
	    //BEGIN ALE + Momentum (Navier-Stokes)
	    { 
	      for (unsigned i=0; i<nve; i++){
		
		//BEGIN redidual Laplacian ALE map in the reference domain
		adept::adouble LapmapVAR[3] = {0., 0., 0.};
		for(int idim=0; idim<dim; idim++) {
		  for(int jdim=0; jdim<dim; jdim++) {
		    LapmapVAR[idim] += (GradSolhatVAR[idim][jdim]*gradphi_hat[i*dim+jdim]) ;
		  }
		}
		for(int idim=0; idim<dim; idim++) {
		  aRhs[indexVAR[idim]][i]+=(!solidmark[i])*(-LapmapVAR[idim]*Weight_nojac);
		}
		//END redidual Laplacian ALE map in reference domain

		//BEGIN redidual Navier-Stokes in moving domain   
		adept::adouble LapvelVAR[3]={0.,0.,0.};
		adept::adouble AdvaleVAR[3]={0.,0.,0.};
		for(int idim=0.; idim<dim; idim++) {
		  for(int jdim=0.; jdim<dim; jdim++) {
		    LapvelVAR[idim]+=GradSolVAR[dim+idim][jdim]*gradphi[i*dim+jdim];
		    AdvaleVAR[idim]+=SolVAR[dim+jdim]*GradSolVAR[dim+idim][jdim]*phi[i];
		  }
		}
		for(int idim=0; idim<dim; idim++) {
		  aRhs[indexVAR[dim+idim]][i]+= (-AdvaleVAR[idim]      	      	      // advection term	
						 -IRe*LapvelVAR[idim]	 	      // viscous dissipation
						 +SolVAR[2*dim]*gradphi[i*dim+idim]   // pressure gradient
						 )*Weight;
		  //END redidual Navier-Stokes in moving domain    
		}
	      } 
	    } 
	    //END ALE + Momentum (Navier-Stokes)      
	    
	    //BEGIN continuity block 
	    {  	    
	      adept::adouble div_vel=0.;
	      for(int i=0; i<dim; i++) {
		div_vel+=GradSolVAR[dim+i][i];
	      }
	      for (unsigned i=0; i<nve1; i++) {
		aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*div_vel)*Weight;
	      }
	    }
	    //END continuity block ===========================
	  }   
	  //END FLUID ASSEMBLY ============
	    
	  //*******************************************************************************************************
	    
	  //BEGIN SOLID ASSEMBLY ============
	  else{	      
	    //BEGIN build Chauchy Stress in moving domain
	    //physical quantity
	    adept::adouble J_hat;
	    adept::adouble I_e;
	    adept::adouble Cauchy[3][3];
	    adept::adouble Id2th[3][3]= {{ 1.,0.,0.}, { 0.,1.,0.}, { 0.,0.,1.}};
	    
	    adept::adouble I1_B=0.;
	    adept::adouble I2_B=0.;
	    
	    if (solid_model==0) { // Saint-Venant
	      adept::adouble e[3][3];
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
		  Cauchy[i][j] = 2*mus*e[i][j]-2*mus*I_e*SolVAR[2*dim]*Id2th[i][j];
		  //+(penalty)*lambda*I_e*Id2th[i][j];
		}
	      }
	    }
	  
	    else { // hyperelastic non linear material
	      adept::adouble F[3][3]={{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
	      adept::adouble B[3][3];      
	    
	      for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
		  F[i][j]+=GradSolhatVAR[i][j];
		}
	      }
		
	      J_hat =   F[0][0]*F[1][1]*F[2][2] + F[0][1]*F[1][2]*F[2][0] + F[0][2]*F[1][0]*F[2][1]
		      - F[2][0]*F[1][1]*F[0][2] - F[2][1]*F[1][2]*F[0][0] - F[2][2]*F[1][0]*F[0][1];
	      for (int I=0; I<3; ++I) {
		for (int J=0; J<3; ++J) {
		  B[I][J]=0.;
		  for (int K=0; K<3; ++K) {
		    //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
		    B[I][J] += F[I][K]*F[J][K];
		  }
		}
	      }	  
	      if( solid_model <=4 ){ // Neo-Hookean
		I1_B = B[0][0] + B[1][1] + B[2][2];
		
	    	for (int I=0; I<3; ++I) {
		  for (int J=0; J<3; ++J) {
		    if	    ( 1 == solid_model ) Cauchy[I][J] = mus*B[I][J] 
							       -mus*I1_B*SolVAR[2*dim]*Id2th[I][J]; 	//Wood-Bonet J_hat  =1;
		    else if ( 2 == solid_model ) Cauchy[I][J] = mus/J_hat*B[I][J] 
							       -mus/J_hat*SolVAR[2*dim]*Id2th[I][J];    //Wood-Bonet J_hat !=1;
		    else if ( 3 == solid_model ) Cauchy[I][J] = mus*(B[I][J] - Id2th[I][J])/J_hat
							      + lambda/J_hat*log(J_hat)*Id2th[I][J]; 	//Wood-Bonet penalty
		    else if ( 4 == solid_model ) Cauchy[I][J] = mus*(B[I][J] - I1_B*Id2th[I][J]/3.)/pow(J_hat,5./3.)
						              + lambda*(J_hat-1.)*Id2th[I][J];  	  //Allan-Bower
		    
		  }
		}
	      }
	      else if ( 5 == solid_model ){ //Mooney-Rivlin
		adept::adouble detB=   B[0][0]* ( B[1][1]*B[2][2] - B[2][1]*B[1][2] ) 
				     - B[0][1]* ( B[2][2]*B[1][0] - B[1][2]*B[2][0] ) 
				     + B[0][2]* ( B[1][0]*B[2][1] - B[2][0]*B[1][1] );
		adept::adouble invdetB=1./detB;		      
		adept::adouble invB[3][3];
		
		invB[0][0] =  (B[1][1]*B[2][2]-B[2][1]*B[1][2])*invdetB;
		invB[1][0] = -(B[0][1]*B[2][2]-B[0][2]*B[2][1])*invdetB;
		invB[2][0] =  (B[0][1]*B[1][2]-B[0][2]*B[1][1])*invdetB;
		invB[0][1] = -(B[1][0]*B[2][2]-B[1][2]*B[2][0])*invdetB;
		invB[1][1] =  (B[0][0]*B[2][2]-B[0][2]*B[2][0])*invdetB;
		invB[2][1] = -(B[0][0]*B[1][2]-B[1][0]*B[0][2])*invdetB;
		invB[0][2] =  (B[1][0]*B[2][1]-B[2][0]*B[1][1])*invdetB;
		invB[1][2] = -(B[0][0]*B[2][1]-B[2][0]*B[0][1])*invdetB;
		invB[2][2] =  (B[0][0]*B[1][1]-B[1][0]*B[0][1])*invdetB;
		
		I1_B = B[0][0] + B[1][1] + B[2][2];
		I2_B = B[0][0] * B[1][1] + B[1][1] * B[2][2] + B[2][2]*B[0][0]
		      -B[0][1] * B[1][0] - B[1][2] * B[2][1] - B[2][0]*B[0][2];
		
		for (int I=0; I<3; ++I) {
		  for (int J=0; J<3; ++J) {
		    Cauchy[I][J] = mus*(2.*B[I][J] - invB[I][J] )/3.
				  -mus*(2.*I1_B - I2_B )/3.*SolVAR[2*dim]*Id2th[I][J];
		    ;
		  }
		}
		
	      }	  		
	    }
	    //END build Chauchy Stress in moving domain
	    
	    //BEGIN v=0 + Momentum (Solid)      
	    {
	      for (unsigned i=0; i<nve; i++) {

		//BEGIN redidual v=0 in fixed domain 
		for(int idim=0; idim<dim; idim++) {
		  aRhs[indexVAR[idim]][i] += (-phi[i]*(-SolVAR[dim+idim]))*Weight_hat;
		}
		for (unsigned j=0; j<nve; j++) { //mass matrix preconditioner in the ALE equation
		   Bmass[i*nve+j] += mass_penalty*phi[i]*phi[j]*Weight_hat;
		}
                //END redidual v=0 in fixed domain 
                
                //BEGIN redidual Solid Momentum in moving domain 
		adept::adouble CauchyDIR[3]={0.,0.,0.};
		for(int idim=0.; idim<dim; idim++) {
		  for(int jdim=0.; jdim<dim; jdim++) {
		    CauchyDIR[idim]+= gradphi[i*dim+jdim]*Cauchy[idim][jdim];
		  }
		}
		for(int idim=0; idim<dim; idim++) {
		  aRhs[indexVAR[dim+idim]][i] += (phi[i]*_gravity[idim] - CauchyDIR[idim])*Weight;			
		} 
		//END redidual Solid Momentum in moving domain 
	      }
	    }
	    //END v=0 + Momentum (Solid)   
	    
	    //BEGIN continuity block
	    {       
	        for (unsigned i=0; i<nve1; i++) {
		if(!penalty){		    
		  if (0 == solid_model) {
		    aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*(I_e + (!incompressible)/lambda*SolVAR[2*dim] ) )*Weight_hat;
		  }
		  else if (1 == solid_model || 5 == solid_model) {
		    aRhs[indexVAR[2*dim]][i] +=phi1[i]*(J_hat-1. + (!incompressible)/lambda*SolVAR[2*dim])*Weight_hat;
		  }
		  else if (2 == solid_model){
		    aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*( log(J_hat)/J_hat + (!incompressible)/lambda*SolVAR[2*dim] ) )*Weight_hat;
		  }
		}
		else if (3 == solid_model || 4 == solid_model){ // pressure = 0 in the solid
		  aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*( SolVAR[2*dim] ) )*Weight_hat;
		}
	      }
	    }
	    //END continuity block
	  }  
	  //END SOLID ASSEMBLY ============
	}
      }
	
      //BEGIN local to global assembly 	
      //copy adouble aRhs into double Rhs
      for (unsigned i=0;i<nve;i++) {
	for(int j=0; j<dim; j++) {
	  Rhs[indexVAR[j]][i]     = aRhs[indexVAR[j]][i].value();
	  Rhs[indexVAR[j+dim]][i] = aRhs[indexVAR[j+dim]][i].value();
	}
      }
      for (unsigned i=0;i<nve1;i++) {
	Rhs[indexVAR[2*dim]][i] = aRhs[indexVAR[2*dim]][i].value();
      }	
//       // store ALE mapping
//       for(int i=0; i<dim; i++) {
// 	myRES->add_vector_blocked(Rhs[indexVAR[i]],dofsVAR[i]);
// 	s.dependent(&aRhs[indexVAR[i]][0], nve);   
// 	s.independent(&Soli[indexVAR[i]][0], nve); 
// 	if(flag_mat!=2){ 
// 	  s.independent(&Soli[indexVAR[i+dim]][0], nve); 
// 	}
// 	s.jacobian(Jac);
// 	for (int inode=0;inode<nve;inode++){
// 	  for (int jnode=0;jnode<nve;jnode++){
// 	    B22[inode*nve+jnode]=-Jac[jnode*nve+inode];
// 	  }
// 	} 
// 	myKK ->add_matrix_blocked(B22,dofsVAR[i],dofsVAR[i]);    
// 	if(flag_mat!=2){ 
// 	  for (int inode=0;inode<nve;inode++){
// 	    for (int jnode=0;jnode<nve;jnode++){
// 	      B22[inode*nve+jnode]=-Jac[nve*nve+jnode*nve+inode];
// 	    }
// 	  }
// 	  myKK->add_matrix_blocked(B22,dofsVAR[i],dofsVAR[i+dim]);
// 	}
// 	s.clear_independents();
// 	s.clear_dependents();
//       }
            
      // store ALE mapping
      s.dependent(&aRhs[indexVAR[0]][0], nve);   
      s.independent(&Soli[indexVAR[0]][0], nve); 
      if(flag_mat!=2){ 
	s.independent(&Soli[indexVAR[dim]][0], nve); 
      }
      s.jacobian(&Jac[0]);
      
      for(int i=0; i<dim; i++) {
	myRES->add_vector_blocked(Rhs[indexVAR[i]],dofsVAR[i]);
	for (int inode=0;inode<nve;inode++){
	  for (int jnode=0;jnode<nve;jnode++){
	    B22[inode*nve+jnode]=-Jac[jnode*nve+inode]+Bmass[inode*nve+jnode];
	  }
	} 
	myKK ->add_matrix_blocked(B22,dofsVAR[i],dofsVAR[i]);    
	if(flag_mat!=2){ 
	  for (int inode=0;inode<nve;inode++){
	    for (int jnode=0;jnode<nve;jnode++){
	      B22[inode*nve+jnode]=-Jac[nve*nve+jnode*nve+inode];
	    }
	  }
	  myKK->add_matrix_blocked(B22,dofsVAR[i],dofsVAR[i+dim]);
	}
      }
      s.clear_independents();
      s.clear_dependents();      
      
      //Store Momentum equation
      for(int i=0; i<dim; i++) {
	myRES->add_vector_blocked(Rhs[indexVAR[dim+i]],dofsVAR[dim+i]);
	s.dependent(&aRhs[indexVAR[dim+i]][0], nve);   	  	  
      }
      for(int j=0; j<2*dim; j++) {  
	s.independent(&Soli[indexVAR[j]][0], nve); 
      }
      s.independent(&Soli[indexVAR[2*dim]][0], nve1);   
      s.jacobian(&Jac[0]);	
      unsigned nveAll=(2*dim*nve+nve1);
      unsigned nveVel=dim*nve;
      for (int inode=0;inode<nveVel;inode++){
	for (int jnode=0;jnode<nveAll;jnode++){
	   Bm[inode*nveAll+jnode]=-Jac[jnode*nveVel+inode];
	}
      }
      myKK->add_matrix_blocked(Bm,dofsVel,dofsAll);
      s.clear_independents();
      s.clear_dependents();
       
      //Store continuity equation
      myRES->add_vector_blocked(Rhs[indexVAR[2*dim]],dofsVAR[2*dim]);      
      s.dependent(&aRhs[indexVAR[2*dim]][0], nve1);   	
      for(int i=0; i<dim; i++) {
	if(flag_mat==2){ //Fluid only
	  s.independent(&Soli[indexVAR[dim+i]][0], nve); 
	}
	else{
	   s.independent(&Soli[indexVAR[i]][0], nve); 
	}
      }
      s.independent(&Soli[indexVAR[2*dim]][0], nve1);
      s.jacobian(&Jac[0]);
      for(int i=0; i<dim; i++) {
	if(flag_mat==2){ //Fluid only
	  for (int inode=0;inode<nve1;inode++){
	    for (int jnode=0;jnode<nve;jnode++){
	      B12[inode*nve+jnode]=-Jac[nve*nve1*i+jnode*nve1+inode];
	    }
	  }
	  myKK->add_matrix_blocked(B12,dofsVAR[2*dim],dofsVAR[dim+i]);
	}
	else{ //Solid only
	  for (int inode=0;inode<nve1;inode++){
	    for (int jnode=0;jnode<nve;jnode++){
	      B12[inode*nve+jnode]=-Jac[nve*nve1*i+jnode*nve1+inode];
	    }
	  }
	  myKK->add_matrix_blocked(B12,dofsVAR[2*dim],dofsVAR[i]);  
	} 
      }
      for (int inode=0;inode<nve1;inode++){
	for (int jnode=0;jnode<nve1;jnode++){
	  B11[inode*nve1+jnode]=-Jac[nve*nve1*dim+jnode*nve1+inode];
	}
      }
      myKK->add_matrix_blocked(B11,dofsVAR[2*dim],dofsVAR[2*dim]);
      
      s.clear_independents();	
      s.clear_dependents();      
      
      //END local to global assembly
   
    } //end list of elements loop

    myKK->close();
    myRES->close();
  
    delete area_elem_first; 
	  
    // *************************************
    end_time=clock();
    AssemblyTime+=(end_time-start_time);
    // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

  }
  
  
  
  
  
  
  
  phi[i*dim+jdim]*Cauchy[idim][jdim];
		  }
		}
		for(int idim=0; idim<dim; idim++) {
		  aRhs[indexVAR[dim+idim]][i] += (phi[i]*_gravity[idim] - CauchyDIR[idim])*Weight;			
		} 
		//END redidual Solid Momentum in moving domain 
	      }
	    }
	    //END v=0 + Momentum (Solid)   
	    
	    //BEGIN continuity block
	    {       
	        for (unsigned i=0; i<nve1; i++) {
		if(!penalty){		    
		  if (0 == solid_model) {
		    aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*(I_e + (!incompressible)/lambda*SolVAR[2*dim] ) )*Weight_hat;
		  }
		  else if (1 == solid_model || 5 == solid_model) {
		    aRhs[indexVAR[2*dim]][i] +=phi1[i]*(J_hat-1. + (!incompressible)/lambda*SolVAR[2*dim])*Weight_hat;
		  }
		  else if (2 == solid_model){
		    aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*( log(J_hat)/J_hat + (!incompressible)/lambda*SolVAR[2*dim] ) )*Weight_hat;
		  }
		}
		else if (3 == solid_model || 4 == solid_model){ // pressure = 0 in the solid
		  aRhs[indexVAR[2*dim]][i] += -(-phi1[i]*( SolVAR[2*dim] ) )*Weight_hat;
		}
	      }
	    }
	    //END continuity block
	  }  
	  //END SOLID ASSEMBLY ============
	}
      }
	
      //BEGIN local to global assembly 	
      //copy adouble aRhs into double Rhs
      for (unsigned i=0;i<nve;i++) {
	for(int j=0; j<dim; j++) {
	  Rhs[indexVAR[j]][i]     = aRhs[indexVAR[j]][i].value();
	  Rhs[indexVAR[j+dim]][i] = aRhs[indexVAR[j+dim]][i].value();
	}
      }
      for (unsigned i=0;i<nve1;i++) {
	Rhs[indexVAR[2*dim]][i] = aRhs[indexVAR[2*dim]][i].value();
      }	
//       // store ALE mapping
//       for(int i=0; i<dim; i++) {
// 	myRES->add_vector_blocked(Rhs[indexVAR[i]],dofsVAR[i]);
// 	s.dependent(&aRhs[indexVAR[i]][0], nve);   
// 	s.independent(&Soli[indexVAR[i]][0], nve); 
// 	if(flag_mat!=2){ 
// 	  s.independent(&Soli[indexVAR[i+dim]][0], nve); 
// 	}
// 	s.jacobian(Jac);
// 	for (int inode=0;inode<nve;inode++){
// 	  for (int jnode=0;jnode<nve;jnode++){
// 	    B22[inode*nve+jnode]=-Jac[jnode*nve+inode];
// 	  }
// 	} 
// 	myKK ->add_matrix_blocked(B22,dofsVAR[i],dofsVAR[i]);    
// 	if(flag_mat!=2){ 
// 	  for (int inode=0;inode<nve;inode++){
// 	    for (int jnode=0;jnode<nve;jnode++){
// 	      B22[inode*nve+jnode]=-Jac[nve*nve+jnode*nve+inode];
// 	    }
// 	  }
// 	  myKK->add_matrix_blocked(B22,dofsVAR[i],dofsVAR[i+dim]);
// 	}
// 	s.clear_independents();
// 	s.clear_dependents();
//       }
            
      // store ALE mapping
      s.dependent(&aRhs[indexVAR[0]][0], nve);   
      s.independent(&Soli[indexVAR[0]][0], nve); 
      if(flag_mat!=2){ 
	s.independent(&Soli[indexVAR[dim]][0], nve); 
      }
      s.jacobian(&Jac[0]);
      
      for(int i=0; i<dim; i++) {
	myRES->add_vector_blocked(Rhs[indexVAR[i]],dofsVAR[i]);
	for (int inode=0;inode<nve;inode++){
	  for (int jnode=0;jnode<nve;jnode++){
	    B22[inode*nve+jnode]=-Jac[jnode*nve+inode]+Bmass[inode*nve+jnode];
	  }
	} 
	myKK ->add_matrix_blocked(B22,dofsVAR[i],dofsVAR[i]);    
	if(flag_mat!=2){ 
	  for (int inode=0;inode<nve;inode++){
	    for (int jnode=0;jnode<nve;jnode++){
	      B22[inode*nve+jnode]=-Jac[nve*nve+jnode*nve+inode];
	    }
	  }
	  myKK->add_matrix_blocked(B22,dofsVAR[i],dofsVAR[i+dim]);
	}
      }
      s.clear_independents();
      s.clear_dependents();      
      
      //Store Momentum equation
      for(int i=0; i<dim; i++) {
	myRES->add_vector_blocked(Rhs[indexVAR[dim+i]],dofsVAR[dim+i]);
	s.dependent(&aRhs[indexVAR[dim+i]][0], nve);   	  	  
      }
      for(int j=0; j<2*dim; j++) {  
	s.independent(&Soli[indexVAR[j]][0], nve); 
      }
      s.independent(&Soli[indexVAR[2*dim]][0], nve1);   
      s.jacobian(&Jac[0]);	
      unsigned nveAll=(2*dim*nve+nve1);
      unsigned nveVel=dim*nve;
      for (int inode=0;inode<nveVel;inode++){
	for (int jnode=0;jnode<nveAll;jnode++){
	   Bm[inode*nveAll+jnode]=-Jac[jnode*nveVel+inode];
	}
      }
      myKK->add_matrix_blocked(Bm,dofsVel,dofsAll);
      s.clear_independents();
      s.clear_dependents();
       
      //Store continuity equation
      myRES->add_vector_blocked(Rhs[indexVAR[2*dim]],dofsVAR[2*dim]);      
      s.dependent(&aRhs[indexVAR[2*dim]][0], nve1);   	
      for(int i=0; i<dim; i++) {
	if(flag_mat==2){ //Fluid only
	  s.independent(&Soli[indexVAR[dim+i]][0], nve); 
	}
	else{
	   s.independent(&Soli[indexVAR[i]][0], nve); 
	}
      }
      s.independent(&Soli[indexVAR[2*dim]][0], nve1);
      s.jacobian(&Jac[0]);
      for(int i=0; i<dim; i++) {
	if(flag_mat==2){ //Fluid only
	  for (int inode=0;inode<nve1;inode++){
	    for (int jnode=0;jnode<nve;jnode++){
	      B12[inode*nve+jnode]=-Jac[nve*nve1*i+jnode*nve1+inode];
	    }
	  }
	  myKK->add_matrix_blocked(B12,dofsVAR[2*dim],dofsVAR[dim+i]);
	}
	else{ //Solid only
	  for (int inode=0;inode<nve1;inode++){
	    for (int jnode=0;jnode<nve;jnode++){
	      B12[inode*nve+jnode]=-Jac[nve*nve1*i+jnode*nve1+inode];
	    }
	  }
	  myKK->add_matrix_blocked(B12,dofsVAR[2*dim],dofsVAR[i]);  
	} 
      }
      for (int inode=0;inode<nve1;inode++){
	for (int jnode=0;jnode<nve1;jnode++){
	  B11[inode*nve1+jnode]=-Jac[nve*nve1*dim+jnode*nve1+inode];
	}
      }
      myKK->add_matrix_blocked(B11,dofsVAR[2*dim],dofsVAR[2*dim]);
      
      s.clear_independents();	
      s.clear_dependents();      
      
      //END local to global assembly
   
    } //end list of elements loop

    myKK->close();
    myRES->close();
  
    delete area_elem_first; 
	  
    // *************************************
    end_time=clock();
    AssemblyTime+=(end_time-start_time);
    // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

  }
  
  
  
  
  
  
  
  
  
  
   
   void IncompressibleFSIAssembly(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assemble_matrix) {
    
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
    
    vector <double> nablaphi;
    vector <double> nablaphi_hat;
  
    metis_node1.reserve(max_size);
    metis_node2.reserve(max_size);
    solidmark.reserve(max_size);
    phi.reserve(max_size);
    phi_hat.reserve(max_size);
    gradphi.reserve(max_size*dim);
    gradphi_hat.reserve(max_size*dim);
    
    nablaphi.reserve(max_size*3*(dim-1));
    nablaphi_hat.reserve(max_size*3*(dim-1));
  
    const double *phi1;
    
    double Weight=0.;
    double Weight_nojac=0.;
    double Weight_hat=0.;
  
    vector <vector < double> > vx(dim);
    vector <vector < double> > vx_face(dim);
    vector <vector < double> > vx_hat(dim);
  
    for(int i=0;i<dim;i++){
      vx[i].reserve(max_size);
      vx_face[i].resize(9);
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
    //double rhos 		= ml_prob.parameters.get<Solid>("Solid").get_density();            
    double mu_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_shear_modulus(); 
    double lambda_lame 	= ml_prob.parameters.get<Solid>("Solid").get_lame_lambda();        
    double mus		= mu_lame/rhof;
    double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();   
    double lambda	= lambda_lame / rhof;
    double betans	= 1.;
    int    solid_model	= ml_prob.parameters.get<Solid>("Solid").get_physical_model();     
       
    bool incompressible=( 0.5 == ml_prob.parameters.get<Solid>("Solid").get_poisson_coeff() )?1:0;
    const bool penalty = ml_prob.parameters.get<Solid>("Solid").get_if_penalty();
        
    //physical quantity
    double Jnp1_hat;
    double Jn_hat;
    double I_Bl;
    double I_e;
    double Cauchy[3][3];
    double Cauchy_old[3][3]; 
    double tg_stiff_matrix[3][3];
    //initialization C tensor: Saint-Venaint Kirchoff model : solid_model==0;
    const double Id2th[3][3]= {{ 1.,0.,0.}, { 0.,1.,0.}, { 0.,0.,1.}};
    double C_mat[3][3][3][3];
    for (int I=0; I<3; ++I) {
      for (int J=0; J<3; ++J) {
	for (int K=0; K<3; ++K) {
	  for (int L=0; L<3; ++L) {
	    C_mat[I][J][K][L] = 2.*mus*Id2th[I][K]*Id2th[J][L]+(penalty)*lambda*Id2th[I][J]*Id2th[K][L];
	  }
	}
      }
    }

    // ale map
    double _mu_ale[3] = {1.,1.,1.};
 
    // gravity
    double _gravity[3]={0.,0.,0.};
  
    // newton algorithm
    bool nwtn_alg = true;

    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));  
    unsigned end_ind2   = mymsh->GetEndIndex(SolType2);

    unsigned SolType1 = ml_sol->GetSolutionType(ml_sol->GetIndex("P"));  
    unsigned end_ind1   = mymsh->GetEndIndex(SolType1);

    // mesh and procs
    unsigned nel    = mymsh->GetElementNumber();
    unsigned igrid  = mymsh->GetLevel();
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
      nablaphi.resize(nve*3*(dim-1));
      nablaphi_hat.resize(nve*3*(dim-1));
      
        
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
	unsigned inode=(SolType1<3)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
	metis_node1[i]=mymsh->GetMetisDof(inode,SolType[2*dim]);
	dofsVAR[2*dim][i]=myLinEqSolver->GetKKDof(indVAR[2*dim],indexVAR[2*dim],inode);
      }
      // ----------------------------------------------------------------------------------------
       
      if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
	  
	  
	  
	  
	{
	  double tau=0.;
	  vector<double> normal(3.0);
	       
	  // loop on faces
	  for(unsigned jface=0; jface<myel->GetElementFaceNumber(kel); jface++) {
		
	    // look for boundary faces
	    if(myel->GetFaceElementIndex(kel,jface)<0) {
	      unsigned int face = -(mymsh->el->GetFaceElementIndex(kel,jface)+1);	      
	      if( !ml_sol->_SetBoundaryConditionFunction(0.,0.,0.,"U",tau,face,0.) && tau!=0.){
		unsigned nve = mymsh->el->GetElementFaceDofNumber(kel,jface,SolType2);
		const unsigned felt = mymsh->el->GetElementFaceType(kel, jface); 
		  		  		  
		for(unsigned i=0; i<nve; i++) {
		  unsigned inode=mymsh->el->GetFaceVertexIndex(kel,jface,i)-1u;
		  unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
		  for(unsigned idim=0; idim<dim; idim++) {
		    vx_face[idim][i]=(*mymsh->_coordinate->_Sol[idim])(inode_Metis)+(*mysolution->_Sol[indVAR[idim]])(inode_Metis);;
		  }
		}
		for(unsigned igs=0; igs < ml_prob._ml_msh->_finiteElement[felt][SolType2]->GetGaussPointNumber(); igs++) {
		  (ml_prob._ml_msh->_finiteElement[felt][SolType2]->*(ml_prob._ml_msh->_finiteElement[felt][SolType2])->Jacobian_sur_ptr)(vx_face,igs,Weight,phi,gradphi,normal);
		  // *** phi_i loop ***
		  for(unsigned i=0; i<nve; i++) {
		    double value = - phi[i]*tau/rhof*Weight;
		    unsigned int ilocal = mymsh->el->GetLocalFaceVertexIndex(kel, jface, i);
		      		      
		    Rhs[indexVAR[dim]][ilocal]   += value*normal[0];
		    Rhs[indexVAR[dim+1]][ilocal] += value*normal[1];
		    Rhs[indexVAR[dim+2]][ilocal] += value*normal[2];
		  }
		}
	      }
	    }
	  }    
	}
	  
	  
	//  *** Gauss point loop ***
	double area=1.;
	for (unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][SolType2]->GetGaussPointNumber(); ig++) {

	  // *** get Jacobian and test function and test function derivatives in the moving frame***
	  ml_prob._ml_msh->_finiteElement[kelt][SolType2]->Jacobian(vx,ig,Weight,phi,gradphi,nablaphi);
	  ml_prob._ml_msh->_finiteElement[kelt][SolType2]->Jacobian(vx_hat,ig,Weight_hat,phi_hat,gradphi_hat,nablaphi_hat);
	  phi1=ml_prob._ml_msh->_finiteElement[kelt][SolType1]->GetPhi(ig);
	  if (flag_mat==2) {
	    if(ig==0){
	      double GaussWeight = ml_prob._ml_msh->_finiteElement[kelt][SolType2]->GetGaussWeight(ig);
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
		  GradSolhatVAR[i][j] += gradphi_hat[inode*dim+j]*soli;
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
		    //Advection term II
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
	      // *** phi_i loop ***
	      for (unsigned i=0; i<nve; i++,gradfi+=dim,fi++) {
		const double *fj=phi1;
		// *** phi_j loop ***
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
	    if (solid_model==0) { // Saint-Venant
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
		  Cauchy[i][j] = 2*mus*e[i][j]+(penalty)*lambda*I_e*Id2th[i][j];
		}
	      }
	    }
	  
	    else { // Neo-Hookean material
	      double F[3][3]={{1.,0.,0.},{0.,1.,0.},{0.,0.,1.}};
	      double Bl[3][3];      
	    
	      for(int i=0;i<dim;i++){
		for(int j=0;j<dim;j++){
		  F[i][j]+=GradSolhatVAR[i][j];
		}
	      }
		
	      Jnp1_hat =  F[0][0]*F[1][1]*F[2][2] + F[0][1]*F[1][2]*F[2][0] + F[0][2]*F[1][0]*F[2][1]
		- F[2][0]*F[1][1]*F[0][2] - F[2][1]*F[1][2]*F[0][0] - F[2][2]*F[1][0]*F[0][1];
			  
	      if(1 == solid_model){ // eugenio's incompressible formulation
		// computation of the the three deformation tensor b
		for (int I=0; I<3; ++I) {
		  for (int J=0; J<3; ++J) {
		    Bl[I][J]=0.;
		    for (int K=0; K<3; ++K) {
		      //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
		      Bl[I][J] += F[I][K]*F[J][K];
		    }
		    Cauchy[I][J] = mus*(Bl[I][J] - Id2th[I][J]);
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
	      else if(2 == solid_model){ // Bonet and Wood nearly incompressible formulation   
		// computation of the the three deformation tensor b
		for (int I=0; I<3; ++I) {
		  for (int J=0; J<3; ++J) {
		    Bl[I][J]=0.;
		    for (int K=0; K<3; ++K) {
		      //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
		      Bl[I][J] += F[I][K]*F[J][K];
		    }
		    Cauchy[I][J] = (mus/Jnp1_hat)*(Bl[I][J] - Id2th[I][J]);
		  }
		}
		I_Bl = Bl[0][0] + Bl[1][1] + Bl[2][2];
		//for the incompressible(nearly incompressible) case
		for (int ii=0; ii<3; ++ii) {
		  for (int jj=0; jj<3; ++jj) {
		    for (int kk=0; kk<3; ++kk) {
		      for (int ll=0; ll<3; ++ll) {
			C_mat[ii][jj][kk][ll] = 2.*mus*pow(Jnp1_hat,-1.6666666666666)*(
										       0.333333333333*I_Bl*Id2th[ii][kk]*Id2th[jj][ll]	//1/3*I_c*i
										       // 										  +0.111111111111*I_Bl*Id2th[ii][jj]*Id2th[kk][ll] 	//1/9*I_b*IxI
										       // 										  -0.333333333333*Bl[ii][jj]*Id2th[kk][ll]		//-1/3*b*I
										       // 										  -0.333333333333*Id2th[ii][jj]*Bl[kk][ll]		//-1/3*b*I
										       )
			  -SolVAR[2*dim]*(Id2th[ii][jj]*Id2th[kk][ll]-2.*Id2th[ii][kk]*Id2th[jj][ll] );  	// -p(IxI-2i)
		      }
		    }
		  }
		}
	      }
	      else if (3 == solid_model){ //Bonet and Wood penalty formulation 
		//computation of the the three deformation tensor b
		for (int I=0; I<3; ++I) {
		  for (int J=0; J<3; ++J) {
		    Bl[I][J]=0.;
		    for (int K=0; K<3; ++K) {
		      //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
		      Bl[I][J] += F[I][K]*F[J][K];
		    }
		    Cauchy[I][J] = (mus/Jnp1_hat)*(Bl[I][J] - Id2th[I][J])+lambda/Jnp1_hat*log(Jnp1_hat)*Id2th[I][J];
		  }
		}
	    
		I_Bl = Bl[0][0] + Bl[1][1] + Bl[2][2];
	    		  
		double mup = (mus-lambda*log(Jnp1_hat))/Jnp1_hat;
		double lambdap = (lambda/Jnp1_hat);
		for (int ii=0; ii<3; ++ii) {
		  for (int jj=0; jj<3; ++jj) {
		    for (int kk=0; kk<3; ++kk) {
		      for (int ll=0; ll<3; ++ll) {
			C_mat[ii][jj][kk][ll] = 2.*mup*Id2th[ii][kk]*Id2th[jj][ll]
			  +lambdap*Id2th[ii][jj]*Id2th[kk][ll];						  
		      }
		    }
		  }
		}
	      }
	      else if (4 == solid_model){ //Allan Bower penalty formulation
		//computation of the the three deformation tensor b
		for (int I=0; I<3; ++I) {
		  for (int J=0; J<3; ++J) {
		    Bl[I][J]=0.;
		    for (int K=0; K<3; ++K) {
		      //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
		      Bl[I][J] += F[I][K]*F[J][K];
		    }
		  }
		}
	    
		I_Bl = Bl[0][0] + Bl[1][1] + Bl[2][2];
	    
		for (int I=0; I<3; ++I) {
		  for (int J=0; J<3; ++J) {
		    Cauchy[I][J] = mus /pow(Jnp1_hat,2./3.) *(Bl[I][J] - I_Bl*Id2th[I][J]/3.)
		      + lambda*Jnp1_hat*(Jnp1_hat-1.)*Id2th[I][J];
		  }
		}
		  
		//for the incompressible(nearly incompressible) case
		for (int ii=0; ii<3; ++ii) {
		  for (int jj=0; jj<3; ++jj) {
		    for (int kk=0; kk<3; ++kk) {
		      for (int ll=0; ll<3; ++ll) {								  
			C_mat[ii][jj][kk][ll] =  mus/pow(Jnp1_hat,2./3.)*( 
									  Id2th[ii][kk]*Bl[jj][ll]+Bl[ii][ll]*Id2th[jj][kk]
									  -(2./3.)*(Bl[ii][jj]*Id2th[kk][ll]+Id2th[ii][jj]*Bl[kk][ll])
									  +(2./3.)*I_Bl*Id2th[ii][jj]*Id2th[kk][ll]/3.) 
			  + lambda*(2.*Jnp1_hat-1.)*Jnp1_hat*Id2th[ii][jj]*Id2th[kk][ll];
			  
		      }
		    }
		  }
		}
	      }	  		
	    }
	    //----------------------------------------------------------------------------------------------------------------------------

	      

	    // Mass + Stiffness operator
	    {
	      const double *gradfi=&gradphi[0];
	      const double *fi=&phi[0];

	      //  *** phi_i loop ***
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
		  Rhs[indexVAR[dim+idim]][i] += (	phi[i]*_gravity[idim]*Weight
							-CauchyDIR[idim]*Weight
							+(!penalty)*SolVAR[2*dim]*gradphi[i*dim+idim]*Weight
							);
		}
              
		//---------------------------------------------------------------------------------------------------------------------------------

		//END RESIDUALS A + Bt block ===========================

		const double *gradfj=&gradphi[0];
		const double *gradfj_hat=&gradphi_hat[0];
		const double *fj=&phi[0];
		// *** phi_j loop ***
		for (unsigned j=0; j<nve; j++,gradfj+=dim,gradfj_hat+=dim,fj++) {
			
		  //  Kinematic equation v = du/dt --> In the steady state we write \deltau^n+1 - \deltav^n+1 = v - 0
		  for(int idim=0; idim<dim; idim++) {
		    // -(v_n+1,eta)
		    B[indexVAR[idim]][indexVAR[dim+idim]][i*nve+j] -= (*(fi))*(*(fj))*Weight_hat;
		    //  (u_n+1,eta)
		    B[indexVAR[idim]][indexVAR[idim]][i*nve+j] 	 += (*(fi))*(*(fj))*Weight_hat;
		  }
			
		  if(1 == solid_model){
		    //  Stiffness operator -- Elasticity equation (Linear or not)
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
		  else{
		    // tangent stiffness matrix
		    for (int idim=0; idim<dim; ++idim) {
		      for (int jdim=0; jdim<dim; ++jdim) {
			tg_stiff_matrix[idim][jdim] = 0.;
			for (int kdim=0; kdim<dim; ++kdim) {
			  for (int ldim=0; ldim<dim; ++ldim) {
			    tg_stiff_matrix[idim][jdim] += (*(gradfi+kdim))*0.25*(
										  C_mat[idim][kdim][jdim][ldim]+C_mat[idim][kdim][ldim][jdim]
										  +C_mat[kdim][idim][jdim][ldim]+C_mat[kdim][idim][ldim][jdim]
										  )*(*(gradfj+ldim));
			  }	
			}
		      }
		    }
                
		    // geometric tangent stiffness matrix
		    double geom_tg_stiff_matrix = 0.;
		    for(int idim=0; idim<dim; ++idim) {
		      for(int jdim=0; jdim<dim; ++jdim) {
			geom_tg_stiff_matrix += (*(gradfi+idim))*Cauchy[idim][jdim]*(*(gradfj+jdim));
		      }
		    }

		    // Stiffness operator -- Elasticity equation (Linear or not)
		    for(int idim=0; idim<dim; idim++) {
		      B[indexVAR[dim+idim]][indexVAR[idim]][i*nve+j] += geom_tg_stiff_matrix*Weight;
		      for(int jdim=0; jdim<dim; jdim++) {
			B[indexVAR[dim+idim]][indexVAR[jdim]][i*nve+j] += tg_stiff_matrix[idim][jdim]*Weight;
		      }
		    }
		  }	
		}
	      }
	    }
	    // ********************************************
	    if(!penalty){
	      { //Gradient of Pressure
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
	    }
		
	    {       
	      const double *fi=phi1;
	      // *** phi_i loop ***
	      for (unsigned i=0; i<nve1; i++,fi++) {
	    
		if(!penalty){
		  //BEGIN RESIDUALS B block ===========================
		    
		  if (0 == solid_model) {
		    Rhs[indexVAR[2*dim]][i] += -(-((*fi))*(I_e + (!incompressible)/lambda*SolVAR[2*dim] ) )*Weight_hat;
		  }
		  else if (1 == solid_model) {
		    Rhs[indexVAR[2*dim]][i] +=(*fi)*(Jnp1_hat-1. + (!incompressible)/lambda*SolVAR[2*dim] )*Weight_hat;
		  }
		  else if (2 == solid_model){
		    Rhs[indexVAR[2*dim]][i] += -(-((*fi))*( log(Jnp1_hat)/Jnp1_hat + (!incompressible)/lambda*SolVAR[2*dim] ) )*Weight_hat;
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
		else{ // no pressure in the solid
		  Rhs[indexVAR[2*dim]][i] += -(-((*fi))*( SolVAR[2*dim] ) )*Weight_hat;
		}
	      }
	    }
	    //  
	    {  // Pressure Mass term
	      const double *fi=phi1;
	      // *** phi_i loop ***
	      for (unsigned i=0; i<nve1; i++,fi++) {
		const double *fj=phi1;
		// *** phi_j loop ***
		for (unsigned j=0; j<nve1; j++,fj++) {
		  if(!penalty){
		    B[indexVAR[2*dim]][indexVAR[2*dim]][i*nve1+j] -= (!incompressible)/lambda*((*fi)*(*fj))*Weight_hat;
		  }
		  else{ // no pressure in the solid
		    B[indexVAR[2*dim]][indexVAR[2*dim]][i*nve1+j] -= ((*fi)*(*fj))*Weight_hat;
		  }
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

    //exit(0);
  }  
  
   
   
   
   
  */ 
   
   
}

#endif
