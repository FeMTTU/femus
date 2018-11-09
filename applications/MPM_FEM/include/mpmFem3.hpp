
using namespace femus;

double beta = 0.25;
double Gamma = 0.5;
double gravity[3] = {0., -9.81, 0.};
double scalingFactor1 = 1.e-2;
double scalingFactor2 = 1.e-6;
Line* linea;

void AssembleMPMSys(MultiLevelProblem& ml_prob) {
  
  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled
  
  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;
  
  //pointers and references
  
  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FEM");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = ml_sol->GetSolutionLevel(level);     // pointer to the solution (level) object
  LinearEquationSolver* myLinEqSolver = my_nnlin_impl_sys._LinSolver[level];  // pointer to the equation (level) object
  
  Mesh* mymsh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* myel = mymsh->el;   // pointer to the elem object in msh (level)
  SparseMatrix* myKK = myLinEqSolver->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* myRES =  myLinEqSolver->_RES;  // pointer to the global residual vector object in pdeSys (level)
    
  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;
   
  const unsigned dim = mymsh->GetDimension();
  
  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  
  // data
  unsigned iproc  = mymsh->processor_id();
  
  vector < double > phi;
  vector < double > phiHat;
  vector < adept::adouble> gradphi;
  vector < double > gradphiHat;
  
  phi.reserve(max_size);
  phiHat.reserve(max_size);
  
  gradphi.reserve(max_size * dim);
  gradphiHat.reserve(max_size * dim);
  
  vector <vector < adept::adouble> > vx(dim); //vx is coordX in assembly of ex30
  vector <vector < double> > vxHat(dim);
  
  vector< vector< adept::adouble > > solD(dim);  // local solution (displacement)
  vector< vector< adept::adouble > > solV(dim);  // local solution (displacement)
  
  vector< vector< double > > solDOld(dim);       // local solution (displacement)
  vector< vector< double > > solVOld(dim);
  vector< vector< double > > solAOld(dim);
   
  vector< vector< adept::adouble > > aRhsD(dim);     // local redidual vector
  vector< vector< adept::adouble > > aRhsV(dim);     // local redidual vector
    
  vector < int > dofsAll;
  
  vector < double > Jac;
  
  adept::adouble weight;
  double weightHat = 0.;
  
  //reading parameters for MPM body
  double density_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_density();
  double E_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_young_module();
  double mu_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_shear_modulus();
  double nu_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_poisson_coeff();
  double lambda_MPM = ml_prob.parameters.get<Solid> ("SolidMPM").get_lame_lambda();
  double K_MPM = E_MPM / (3.*(1. - 2. * nu_MPM)); //bulk modulus
  
  //reading parameters for FEM body
  double density_FEM = ml_prob.parameters.get<Solid> ("SolidFEM").get_density();
  double E_FEM = ml_prob.parameters.get<Solid> ("SolidFEM").get_young_module();
  double mu_FEM = ml_prob.parameters.get<Solid> ("SolidFEM").get_lame_shear_modulus();
  double nu_FEM = ml_prob.parameters.get<Solid> ("SolidFEM").get_poisson_coeff();
  double lambda_FEM = ml_prob.parameters.get<Solid> ("SolidFEM").get_lame_lambda();
  double K_FEM = E_FEM / (3.*(1. - 2. * nu_FEM)); //bulk modulus
  
  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  
  //vector < adept::adouble >* nullAdoublePointer = NULL;
  // nullDoublePointer = NULL;
  
  //variable-name handling
  const char varname[10][5] = {"DX", "DY", "DZ", "VX", "VY", "VZ", "AX", "AY", "AZ", "Mat"};
  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexSolA(dim);
  vector <unsigned> indexPdeD(dim);
  vector <unsigned> indexPdeV(dim);
  unsigned solType = ml_sol->GetSolutionType(&varname[0][0]);
  
  unsigned indexSolM = ml_sol->GetIndex("M");
  vector < bool > solidFlag;  
  
  for(unsigned k = 0; k < dim; k++) {
    indexSolD[k] = ml_sol->GetIndex(&varname[k][0]);
    indexSolV[k] = ml_sol->GetIndex(&varname[k + 3][0]);  
    if(ml_sol->GetIfFSI()){
      indexSolA[k] = ml_sol->GetIndex(&varname[k + 6][0]);
    }
    indexPdeD[k] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[k][0]);
    indexPdeV[k] = my_nnlin_impl_sys.GetSolPdeIndex(&varname[k + 3][0]);
  }
  
  unsigned indexSolNF =  ml_sol->GetIndex("NF");
  
  unsigned indexSolMat = ml_sol->GetIndex(&varname[9][0]);
  unsigned solTypeMat = ml_sol->GetSolutionType(&varname[9][0]);
  
  start_time = clock();
  
  myKK->zero();
  myRES->zero();
    
  //BEGIN loop on elements (to initialize the "soft" stiffness matrix)
  for(int iel = mymsh->_elementOffset[iproc]; iel < mymsh->_elementOffset[iproc + 1]; iel++) {
    
    short unsigned ielt = mymsh->GetElementType(iel);
    
    int material = mymsh->GetElementMaterial(iel);
    
    unsigned nDofs = mymsh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    
    // resize local arrays
    dofsAll.resize(2 * dim * nDofs);
    
    solidFlag.resize(nDofs);
    
    for(unsigned  k = 0; k < dim; k++) {
      solD[k].resize(nDofs);
      solDOld[k].resize(nDofs);
      
      solV[k].resize(nDofs);
      solVOld[k].resize(nDofs);
      
      vx[k].resize(nDofs);
      vxHat[k].resize(nDofs);
      if( material == 4 ) {
        solAOld[k].resize(nDofs);
      }
    }
       
    for(unsigned  k = 0; k < dim; k++) {
      aRhsD[k].assign(nDofs,0.);  
      aRhsV[k].assign(nDofs,0.);   
    }
    
    // local storage of local to global mapping and solution
    for(unsigned i = 0; i < nDofs; i++) {
      unsigned idof = mymsh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
      unsigned idofX = mymsh->GetSolutionDof(i, iel, 2);    // global to global mapping between solution node and solution dof
      
      solidFlag[i] = ( (*mysolution->_Sol[indexSolM])(idof) > 0.5 || mymsh->GetSolidMark(idof) ) ? true:false;
      
      for(unsigned  k = 0; k < dim; k++) {
        solD[k][i] = (*mysolution->_Sol[indexSolD[k]])(idof);      // global extraction and local storage for the solution
        solDOld[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idof);      // global extraction and local storage for the solution
        
        solV[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
        solVOld[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof); 
        
        dofsAll[i + k * nDofs] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel);    // global to global mapping between
        dofsAll[i + (dim + k) * nDofs] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);    // global to global mapping between
                
        vxHat[k][i] = (*mymsh->_topology->_Sol[k])(idofX);
               
        if( material == 4 ) {
          solAOld[k][i] = (*mysolution->_Sol[indexSolA[k]])(idof);
        }
      }
    }
    
    // start a new recording of all the operations involving adept::adouble variables
    s.new_recording();
   
    for(unsigned i = 0; i < nDofs; i++) { 
      for(unsigned k = 0; k < dim; k++) {
        vx[k][i] = vxHat[k][i] + solD[k][i];
      }
    }
    
    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < mymsh->_finiteElement[ielt][solType]->GetGaussPointNumber(); ig++) {
            
      mymsh->_finiteElement[ielt][solType]->Jacobian(vxHat, ig, weightHat, phiHat, gradphiHat);
           
      vector < vector < adept::adouble > > GradSolDgHat(dim);
      vector < vector < adept::adouble > > GradSolVgHat(dim);
      for(unsigned  k = 0; k < dim; k++) {
        GradSolDgHat[k].assign(dim,0.);
        GradSolVgHat[k].assign(dim,0.);
      }
      for(unsigned i = 0; i < nDofs; i++) {
        for(unsigned j = 0; j < dim; j++) {
          for(unsigned  k = 0; k < dim; k++) {
            GradSolDgHat[k][j] += gradphiHat[i * dim + j] * solD[k][i];
            GradSolVgHat[k][j] += gradphiHat[i * dim + j] * solV[k][i];
          }
        }
      }
      
      if(material == 2) {
        unsigned idofMat = mymsh->GetSolutionDof(0, iel, solTypeMat);
        double  MPMmaterial = (*mysolution->_Sol[indexSolMat])(idofMat);
        double scalingFactor = 0;
        if(MPMmaterial < 5) scalingFactor = scalingFactor1;
        else if(MPMmaterial < 9) scalingFactor = scalingFactor2;
        
        double mu = mu_MPM; 
        
        for(unsigned i = 0; i < nDofs; i++) {
          vector < adept::adouble > DStiffness(dim, 0.);
          vector < adept::adouble > VStiffness(dim, 0.);
          
          for(unsigned j = 0; j < dim; j++) {
            for(unsigned  k = 0; k < dim; k++) {
              DStiffness[k]   +=  gradphiHat[i * dim + j] * (GradSolDgHat[k][j] + GradSolDgHat[j][k]);
              VStiffness[k]   +=  gradphiHat[i * dim + j] * (GradSolVgHat[k][j] + GradSolVgHat[j][k]);
            }
          }
          for(unsigned  k = 0; k < dim; k++) { //Soft stiffness matrix
            if(MPMmaterial >= 2){
              aRhsD[k][i] += - DStiffness[k] * weightHat * mu * scalingFactor;
            }
            else if( !solidFlag[i] ){ //Mass Matrix
              aRhsD[k][i] += DStiffness[k] * weightHat;
              aRhsV[k][i] += VStiffness[k] * weightHat;
              
//               std::cout<<"A";
            }
          }
        }
      }
      else {
        
        mymsh->_finiteElement[ielt][solType]->Jacobian(vx, ig, weight, phi, gradphi);
        
        vector < adept::adouble > SolDg(dim, 0);
        vector < double > SolDgOld(dim, 0);
        vector < adept::adouble > SolVg(dim, 0);
        vector < double > SolVgOld(dim, 0);
        vector < double > SolAgOld(dim, 0);
        for(unsigned i = 0; i < nDofs; i++) {
          for(unsigned  k = 0; k < dim; k++) {
            SolDg[k] += phi[i] * solD[k][i];
            SolDgOld[k] += phi[i] * solDOld[k][i];
            SolVg[k] += phi[i] * solV[k][i];
            SolVgOld[k] += phi[i] * solVOld[k][i];
            SolAgOld[k] += phi[i] * solAOld[k][i];
          }
        }
        
        adept::adouble F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        adept::adouble B[3][3];
        adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
        adept::adouble Cauchy[3][3];
        
        for(int i = 0; i < dim; i++) {
          for(int j = 0; j < dim; j++) {
            F[i][j] += GradSolDgHat[i][j];
          }
        }
        
        adept::adouble J_hat =  F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
        - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];
        
        for(int i = 0; i < 3; i++) {
          for(int j = 0; j < 3; j++) {
            B[i][j] = 0.;
            
            for(int k = 0; k < 3; k++) {//left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
              B[i][j] += F[i][k] * F[j][k];
            }
          }
        }
        
        adept::adouble I1_B = B[0][0] + B[1][1] + B[2][2];
        
        for(int i = 0; i < 3; i++) {
          for(int j = 0; j < 3; j++) {
            Cauchy[i][j] = lambda_FEM * log(J_hat) / J_hat * Id2th[i][j] + mu_FEM / J_hat * (B[i][j] - Id2th[i][j]); 
          }
        }
        
        for(unsigned i = 0; i < nDofs; i++) {
          adept::adouble CauchyDIR[3] = {0., 0., 0.};
          
          for(int j = 0.; j < dim; j++) {
            for(int k = 0.; k < dim; k++) {
              CauchyDIR[j] += gradphi[i * dim + k] * Cauchy[j][k];
            }
          }
          
          for(int k = 0; k < dim; k++) {
            aRhsD[k][i] += (phi[i] * density_FEM / J_hat * gravity[k] * 0. - CauchyDIR[k] 
            - phi[i] * ( SolVg[k] - SolVgOld[k] ) / dt 
            ) * weight;
            
            aRhsV[k][i] += phi[i] * ( 0.5 * (solV[k][i] + solVOld[k][i]) - (solD[k][i] - solDOld[k][i])/dt ) * weight;
          }
        }
      }
    } // end gauss point loop
    
    std::vector<double> Rhs(2 * dim * nDofs);  //resize
    
    for(int i = 0; i < nDofs; i++) {
      for(unsigned  k = 0; k < dim; k++) {
        Rhs[ i +  k * nDofs ] = -aRhsD[k][i].value();
        Rhs[ i +  (k + dim) * nDofs ] = -aRhsV[k][i].value();
      }
    }
    
    myRES->add_vector_blocked(Rhs, dofsAll);
    
   
    Jac.resize(4 * dim * dim * nDofs * nDofs);
    // define the dependent variables
      
    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&aRhsD[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.dependent(&aRhsV[k][0], nDofs);
    }
      
    // define the independent variables
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solD[k][0], nDofs);
    }
    for(unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofs);
    }
      
    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0] , true);
    
//     for(unsigned i = 0; i< 2 * dim * nDofs; i++){
//       for(unsigned j = 0; j< 2 * dim* nDofs; j++){
//         std::cout << Jac[i* ( 2 * dim * nDofs) + j] << " ";
//       }
//       std::cout<<std::endl;
//     }
//     std::cout<<std::endl;
    
    myKK->add_matrix_blocked(Jac, dofsAll, dofsAll);
      
    s.clear_independents();
    s.clear_dependents();
  }
  //END building "soft" stiffness matrix
  
  
  //BEGIN loop on particles (used as Gauss points)
  std::vector < std::vector < adept::adouble > > solD1(dim);
  vector <vector < adept::adouble> > vx1(dim); //vx is coordX in assembly of ex30
  
  std::vector < std::vector < adept::adouble > > solD2(dim);
  vector <vector < adept::adouble> > vx2(dim); //vx is coordX in assembly of ex30
  
  //line instances
  std::vector<unsigned> markerOffset = linea->GetMarkerOffset();
  std::vector<Marker*> particles = linea->GetParticles();
  
  
  vector < double > phi1;
  vector < adept::adouble> gradphi1;
  adept::adouble weight1;
  
  vector < double > phi2;
  vector < adept::adouble> gradphi2;
  adept::adouble weight2;
  
  //initialization of iel
  unsigned ielOld = UINT_MAX;
  for(unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) { //element of particle iMarker
      
    unsigned iel = particles[iMarker]->GetMarkerElement();
    
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofs;
      
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {
        
        ielt = mymsh->GetElementType(iel);
        nDofs = mymsh->GetElementDofNumber(iel, solType);    // number of solution element dofs
    
        // resize local arrays
        dofsAll.resize(2 * dim * nDofs);
    
        for(unsigned  k = 0; k < dim; k++) {
          solD[k].resize(nDofs);
          solDOld[k].resize(nDofs);
          
          solV[k].resize(nDofs);
          solVOld[k].resize(nDofs);
          
          vx1[k].resize(nDofs);
          vx2[k].resize(nDofs);
          vxHat[k].resize(nDofs);
        }
       
        for(unsigned  k = 0; k < dim; k++) {
          aRhsD[k].assign(nDofs,0.);    //resize
          aRhsV[k].assign(nDofs,0.);    //resize
        }
    
        // local storage of global mapping and solution
        for(unsigned i = 0; i < nDofs; i++) {
          unsigned idof = mymsh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
          unsigned idofX = mymsh->GetSolutionDof(i, iel, 2);    // global to global mapping between solution node and solution dof
            
          for(unsigned  k = 0; k < dim; k++) {
            solDOld[k][i] = (*mysolution->_SolOld[indexSolD[k]])(idof);
            solD[k][i] = (*mysolution->_Sol[indexSolD[k]])(idof) + solDOld[k][i] ; 
            
            solV[k][i] = (*mysolution->_Sol[indexSolV[k]])(idof);
            solVOld[k][i] = (*mysolution->_SolOld[indexSolV[k]])(idof); 
                        
            dofsAll[i + k * nDofs] = myLinEqSolver->GetSystemDof(indexSolD[k], indexPdeD[k], i, iel);    // global to global mapping between
            dofsAll[i + (dim + k) * nDofs] = myLinEqSolver->GetSystemDof(indexSolV[k], indexPdeV[k], i, iel);    // global to global mapping between
            
            vxHat[k][i] = (*mymsh->_topology->_Sol[k])(idofX) + solDOld[k][i] ;
               
          }
        }
               
        // start a new recording of all the operations involving adept::adouble variables
        s.new_recording();
        
        for(int k = 0; k < dim; k++) {
          solD1[k].resize(nDofs);
          solD2[k].resize(nDofs);
          for(unsigned i = 0; i < nDofs; i++){
            solD1[k][i] = solD[k][i];
            solD2[k][i] = solD[k][i];
          }
        }
        
        unsigned ii[9][3][3] = { 
          { {0,3,7}, {1,2,5}, {4,6,8}},
          { {1,0,4}, {2,3,6}, {5,7,8}},
          { {2,1,5}, {3,0,7}, {6,4,8}},
          { {3,2,6}, {0,1,4}, {7,5,8}} };
                
        double neumannFactor =  1. * (*mysolution->_Sol[indexSolNF])(iel);        
        if( neumannFactor > 1.0e-10 ){
          
          for(unsigned iface = 0; iface < 4; iface++ ){
            int faceIndex = myel->GetBoundaryIndex(iel, iface);
            bool solidMark = false;
            if( ml_sol->GetIfFSI() ){
              unsigned idof = mymsh->GetSolutionDof(ii[iface][2][0], iel, solType);
              solidMark = mymsh->GetSolidMark(idof);
            }
            if( faceIndex == 1) {
              //std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAA";
              for(unsigned i = 0;i < 3; i++){
                for(int k = 0; k < dim; k++) {
                  unsigned i0 = ii[iface][i][0];
                  unsigned i1 = ii[iface][i][1];
                  unsigned i2 = ii[iface][i][2];  
                  solD1[k][i0] = neumannFactor * (- 1./3. * solD[k][i1] + 4./3 * solD[k][ i2 ] ) + (1. - neumannFactor) * solD[k][i0];
                }
              }
              break;
            }
            if( solidMark ){
              
              //std::cout << "BBBBBBBBBBBBBBBBBBBBBBBBBBB\n";
              
              for(unsigned i = 0;i < 3; i++){
                for(int k = 0; k < dim; k++) {
                  unsigned i0 = ii[iface][i][0];
                  unsigned i1 = ii[iface][i][1];
                  unsigned i2 = ii[iface][i][2];  
                  solD1[k][i0] = neumannFactor * (- 1./3. * solD[k][i1] + 4./3 * solD[k][ i2 ] ) + (1. - neumannFactor) * solD[k][i0];
                }
              }
              
              for(unsigned i = 0;i < 3; i++){
                for(int k = 0; k < dim; k++) {
                  unsigned i0 = ii[iface][i][0];
                  unsigned i1 = ii[iface][i][1];
                  unsigned i2 = ii[iface][i][2];  
                  solD2[k][i1] = (1. - neumannFactor) * solD[k][i1] + neumannFactor * solD[k][i0];
                  solD2[k][i2] = (1. - neumannFactor) * solD[k][i2] + neumannFactor * solD[k][i0];
                }
              }
              break;
            }
            
          }
          
        }
       
        for(int j = 0; j < dim; j++) {
          for(unsigned i = 0; i < nDofs; i++){
            vx1[j][i] = vxHat[j][i] + solD1[j][i];
            vx2[j][i] = vxHat[j][i] + solD2[j][i];
          }
        }
      
      }
      
      // the local coordinates of the particles are the Gauss points in this context
      std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();
      mymsh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradphiHat);
      
      std::vector <double> SolVpOld(dim);
      particles[iMarker]->GetMarkerVelocity(SolVpOld);
      
      std::vector <double> SolApOld(dim);
      particles[iMarker]->GetMarkerAcceleration(SolApOld);
      
      double mass = particles[iMarker]->GetMarkerMass();
      
      mymsh->_finiteElement[ielt][solType]->Jacobian(vx1, xi, weight1, phi1, gradphi1); //function to evaluate at the particles
      mymsh->_finiteElement[ielt][solType]->Jacobian(vx2, xi, weight2, phi2, gradphi2); //function to evaluate at the particles  
        
      // displacement and velocity
      //BEGIN evaluates SolDp1 at the particle iMarker
      vector<adept::adouble> SolDp1(dim,0.);
      vector<vector < adept::adouble > > GradSolDpHat1(dim);
      vector<adept::adouble> SolDp2(dim,0.);
      vector<vector < adept::adouble > > GradSolDpHat2(dim);
      
      vector<adept::adouble> SolVp(dim,0.);
      
      for(int k = 0; k < dim; k++) {
        GradSolDpHat1[k].assign(dim,0.);
        GradSolDpHat2[k].assign(dim,0.);
      }
        
      for(int k = 0; k < dim; k++) {
        for(unsigned i = 0; i < nDofs; i++) {
          SolDp1[k] += phi1[i] * solD[k][i];
          SolDp2[k] += phi2[i] * solD[k][i];
          SolVp[k] += phi1[i] * solV[k][i];
          for(int j = 0; j < dim; j++) {
            GradSolDpHat1[k][j] +=  gradphiHat[i * dim + j] * solD1[k][i];
            GradSolDpHat2[k][j] +=  gradphiHat[i * dim + j] * solD2[k][i];
          }
        }
      }
      //END evaluates SolDp at the particle iMarker
        
        
        
      //BEGIN computation of the Cauchy Stress
      std::vector < std::vector < double > > FpOld;
      FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient
      
      adept::adouble Id2th[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};
      
      adept::adouble FpNew1[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      adept::adouble F1[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble B1[3][3];
      adept::adouble Cauchy1[3][3];
      
      adept::adouble FpNew2[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
      adept::adouble F2[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
      adept::adouble B2[3][3];
      adept::adouble Cauchy2[3][3];
        
      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          FpNew1[i][j] += GradSolDpHat1[i][j];
          FpNew2[i][j] += GradSolDpHat2[i][j];
        }
      }
        
      for(int i = 0; i < dim; i++) {
        for(int j = 0; j < dim; j++) {
          for(int k = 0; k < dim; k++) {
            F1[i][j] += FpNew1[i][k] * FpOld[k][j];
            F2[i][j] += FpNew2[i][k] * FpOld[k][j];
          }
        }
      }
        
      if(dim == 2) F1[2][2] = 1.;
      if(dim == 2) F2[2][2] = 1.;
      
      adept::adouble J_hat1 =  F1[0][0] * F1[1][1] * F1[2][2] + 
                               F1[0][1] * F1[1][2] * F1[2][0] + 
                               F1[0][2] * F1[1][0] * F1[2][1] - 
                               F1[2][0] * F1[1][1] * F1[0][2] - 
                               F1[2][1] * F1[1][2] * F1[0][0] - 
                               F1[2][2] * F1[1][0] * F1[0][1];
                               
      adept::adouble J_hat2 =  F2[0][0] * F2[1][1] * F2[2][2] + 
                               F2[0][1] * F2[1][2] * F2[2][0] + 
                               F2[0][2] * F2[1][0] * F2[2][1] - 
                               F2[2][0] * F2[1][1] * F2[0][2] - 
                               F2[2][1] * F2[1][2] * F2[0][0] - 
                               F2[2][2] * F2[1][0] * F2[0][1];
        
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          B1[i][j] = 0.;
          B2[i][j] = 0.;
          
          for(int k = 0; k < 3; k++) {
            //left Cauchy-Green deformation tensor or Finger tensor (B = F*F^T)
            B1[i][j] += F1[i][k] * F1[j][k];
            B2[i][j] += F2[i][k] * F2[j][k];
          }
        }
      }
        
      adept::adouble I1_B1 = B1[0][0] + B1[1][1] + B1[2][2];
      adept::adouble I1_B2 = B2[0][0] + B2[1][1] + B2[2][2];
        
      double mu = mu_MPM;
      double lambda = lambda_MPM; 
      
      for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
          Cauchy1[i][j] = lambda * log(J_hat1) / J_hat1 * Id2th[i][j] + mu / J_hat1 * (B1[i][j] - Id2th[i][j]); //alternative formulation
          Cauchy2[i][j] = lambda * log(J_hat2) / J_hat2 * Id2th[i][j] + mu / J_hat2 * (B2[i][j] - Id2th[i][j]); //alternative formulation
        }
      }
      //END computation of the Cauchy Stress
        
      //BEGIN redidual Solid Momentum in moving domain
      for(unsigned i = 0; i < nDofs; i++) {
        adept::adouble CauchyDIR1[3] = {0., 0., 0.};
        adept::adouble CauchyDIR2[3] = {0., 0., 0.};
        
        for(int j = 0.; j < dim; j++) {
          for(int k = 0.; k < dim; k++) {
            CauchyDIR1[j] += gradphi1[i * dim + k] * Cauchy1[j][k];
            CauchyDIR2[j] += gradphi2[i * dim + k] * Cauchy2[j][k];
          }
        }
        
        unsigned idof = mymsh->GetSolutionDof(i, iel, solType);
        bool solidMark = mymsh->GetSolidMark(idof);
        for(int k = 0; k < dim; k++) {
          if(!solidMark){
            aRhsD[k][i] += (phi1[i] * gravity[k] - J_hat1 * CauchyDIR1[k] / density_MPM 
              - phi1[i] * (SolVp[k] - SolVpOld[k])/dt
              ) * mass;
          }
          else{
            aRhsD[k][i] += (phi2[i] * gravity[k] - J_hat2 * CauchyDIR2[k] / density_MPM 
               - phi2[i] * (SolVp[k] - SolVpOld[k])/dt
               ) * mass;
          }
          aRhsV[k][i] += phiHat[i] * ( 0.5 * (solV[k][i] + solVOld[k][i]) - solD[k][i] / dt) * weightHat;
        }
      }
      //END redidual Solid Momentum in moving domain
     
        
      if(iMarker == markerOffset[iproc + 1] - 1 || iel != particles[iMarker + 1]->GetMarkerElement()) {
          
        
        std::vector<double> Rhs(2 * dim * nDofs);  //resize
        
        for(int i = 0; i < nDofs; i++) {
          for(unsigned  k = 0; k < dim; k++) {
            Rhs[ i +  k * nDofs ] = -aRhsD[k][i].value();
            Rhs[ i +  (k + dim) * nDofs ] = -aRhsV[k][i].value();
          }
        }
        
        myRES->add_vector_blocked(Rhs, dofsAll);
        
        
        Jac.resize(4 * dim * dim * nDofs * nDofs);
        // define the dependent variables
        
        for(unsigned  k = 0; k < dim; k++) {
          s.dependent(&aRhsD[k][0], nDofs);
        }
        for(unsigned  k = 0; k < dim; k++) {
          s.dependent(&aRhsV[k][0], nDofs);
        }
        
        // define the independent variables
        for(unsigned  k = 0; k < dim; k++) {
          s.independent(&solD[k][0], nDofs);
        }
        for(unsigned  k = 0; k < dim; k++) {
          s.independent(&solV[k][0], nDofs);
        }
        
        // get the and store jacobian matrix (row-major)
        s.jacobian(&Jac[0] , true);
        myKK->add_matrix_blocked(Jac, dofsAll, dofsAll);
        
//         for(unsigned i = 0; i< 2 * dim * nDofs; i++){
//           for(unsigned j = 0; j< 2 * dim* nDofs; j++){
//             std::cout << Jac[i* ( 2 * dim * nDofs) + j] << " ";
//           }
//           std::cout<<std::endl;
//         }
//         std::cout<<std::endl;
        
        s.clear_independents();
        s.clear_dependents();
    
      }
      //END local to global assembly
        
      ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on particles
  
  myRES->close();
  mysolution->_Sol[indexSolMat]->close();
  
  myKK->close();
  
  // *************************************
  end_time = clock();
  AssemblyTime += (end_time - start_time);
  // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************
  
}


void SetNeumannFactor(MultiLevelSolution* ml_sol) {
  
  //pointers and references
    
  unsigned level = ml_sol->_mlMesh->GetNumberOfLevels() - 1;
  
  Solution  *sol = ml_sol->GetSolutionLevel(level);
  Mesh      *msh = ml_sol->_mlMesh->GetLevel(level);
  
  elem* el = msh->el;   // pointer to the elem object in msh (level)
  
  const unsigned dim = msh->GetDimension();
  
  //variable-name handling
  unsigned solNFType = ml_sol->GetSolutionType("NF");
  unsigned solNFIndex = ml_sol->GetIndex("NF");
  
  unsigned solType = ml_sol->GetSolutionType("DX");
  unsigned SolVIndex;
  if(ml_sol->GetIfFSI()){
    unsigned SolVIndex = ml_sol->GetIndex("VY");
  }
    
  unsigned iproc  = msh->processor_id();
    
  //line instances
  std::vector<unsigned> markerOffset = linea->GetMarkerOffset();
  unsigned markerOffset1 = markerOffset[iproc];
  unsigned markerOffset2 = markerOffset[iproc + 1];
  std::vector<Marker*> particles = linea->GetParticles();
  
  //initialization of iel
  unsigned ielOld = UINT_MAX;
  
  //BEGIN loop on particles (used as Gauss points)
  std::vector < bool > solidMark;
  std::vector < double > velMeshOld;
    
  sol->_Sol[solNFIndex]->zero();
  
  double massNF = 0;
  double massElement = 1;
  
  for(unsigned iMarker = markerOffset1; iMarker < markerOffset2; iMarker++) {
    //element of particle iMarker
    unsigned iel = particles[iMarker]->GetMarkerElement();
    if(iel != UINT_MAX) {
      short unsigned ielt;
      unsigned nDofsV;
      
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {
        
        massNF = 0.;
        massElement = 0.;
        
        ielt = msh->GetElementType(iel);
        nDofsV = msh->GetElementDofNumber(iel, solType);
                        
        if( ml_sol->GetIfFSI() ){
          solidMark.resize(nDofsV);
          velMeshOld.resize(nDofsV);
          for(unsigned i = 0; i < nDofsV; i++) {
            unsigned idof = msh->GetSolutionDof(i, iel, solType);
            solidMark[i] = msh->GetSolidMark(idof);
            velMeshOld[i] = (*sol->_Sol[SolVIndex])(idof);
          }
        }
        else {
          solidMark.assign(nDofsV,false);
          velMeshOld.assign(nDofsV,0.);
        }
        //END
      }
          
      std::vector <double> SolVpOld(dim);
      particles[iMarker]->GetMarkerVelocity(SolVpOld);
           
      double massParticle = particles[iMarker]->GetMarkerMass();
                
      massElement += massParticle;  
      
      unsigned ii[9][3][3] = { 
        { {0,3,7}, {1,2,5}, {4,6,8}},
        { {1,0,4}, {2,3,6}, {5,7,8}},
        { {2,1,5}, {3,0,7}, {6,4,8}},
        { {3,2,6}, {0,1,4}, {7,5,8}} };
      
      for(unsigned iface = 0; iface < 4; iface++){
        int faceIndex = el->GetBoundaryIndex(iel, iface);
        unsigned im = ii[iface][2][0];
        if( ( faceIndex == 1  &&  SolVpOld[1] > 0 ) 
           || ( solidMark[im] && (SolVpOld[1] - velMeshOld[im] > 0 ) ) ){
          massNF += massParticle;
          break;
        }
      }
      
      sol->_Sol[solNFIndex]->set(iel, massNF/massElement);
      
      ielOld = iel;
    }
    else {
      break;
    }       
  }
  sol->_Sol[solNFIndex]->close();
}


void GridToParticlesProjection(MultiLevelProblem & ml_prob, Line & linea) {
  
  // ml_prob is the global object from/to where get/set all the data
  // level is the level of the PDE system to be assembled
  // levelMax is the Maximum level of the MultiLevelProblem
  // assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled
  
  clock_t AssemblyTime = 0;
  clock_t start_time, end_time;
  
  //pointers and references
  
  TransientNonlinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<TransientNonlinearImplicitSystem> ("MPM_FEM");
  //NonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem> ("MPM_FEM");
  const unsigned  level = my_nnlin_impl_sys.GetLevelToAssemble();
  MultiLevelSolution* ml_sol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* mysolution = ml_sol->GetSolutionLevel(level);     // pointer to the solution (level) object
  
  Mesh* mymsh = ml_prob._ml_msh->GetLevel(level);     // pointer to the mesh (level) object
  elem* myel = mymsh->el;   // pointer to the elem object in msh (level)
  
  double dt =  my_nnlin_impl_sys.GetIntervalTime();
  
  const unsigned dim = mymsh->GetDimension();
  
  // data
  unsigned iproc  = mymsh->processor_id();
  
  // local objects
  vector< vector < double > > solD(dim);
  vector< vector < double > > solDOld(dim);
  
  
  vector < double > phiHat;
  vector < double > gradphiHat;
  vector < double > nablaphiHat;
  
  vector <vector < double> > vxHat(dim); //vx is coordX in assembly of ex30
  
  double weightHat;
  
  //variable-name handling
  const char varname[9][3] = {"DX", "DY", "DZ", "VX", "VY", "VW", "AX", "AY", "AW"};
  vector <unsigned> indexSolD(dim);
  vector <unsigned> indexSolV(dim);
  vector <unsigned> indexSolA(dim);
  unsigned solType = ml_sol->GetSolutionType(&varname[0][0]);
  
  for(unsigned k = 0; k < dim; k++) {
    indexSolD[k] = ml_sol->GetIndex(&varname[k][0]);
    if(ml_sol->GetIfFSI()){
      indexSolV[k] = ml_sol->GetIndex(&varname[k + 3][0]);
      indexSolA[k] = ml_sol->GetIndex(&varname[k + 6][0]);
    }
  }
  
  unsigned indexSolNF =  ml_sol->GetIndex("NF");
  
  //line instances
  std::vector<unsigned> markerOffset = linea.GetMarkerOffset();
  std::vector<Marker*> particles = linea.GetParticles();
    
  std::vector < std::vector < double > > solD1(dim);
  unsigned ielOld = UINT_MAX;
  
  //BEGIN loop on particles
  for(unsigned iMarker = markerOffset[iproc]; iMarker < markerOffset[iproc + 1]; iMarker++) {
    
    //element of particle iMarker
    unsigned iel = particles[iMarker]->GetMarkerElement();
        
    if(iel != UINT_MAX) {
      
      short unsigned ielt;
      unsigned nDofs;
      
      //update element related quantities only if we are in a different element
      if(iel != ielOld) {
        
        ielt = mymsh->GetElementType(iel);
        nDofs = mymsh->GetElementDofNumber(iel, solType);
        
        for(int i = 0; i < dim; i++) {
          solD[i].resize(nDofs);
          solDOld[i].resize(nDofs);
          vxHat[i].resize(nDofs);
        }
        
        //BEGIN copy of the value of Sol at the dofs idof of the element iel
        for(unsigned inode = 0; inode < nDofs; inode++) {
          unsigned idof = mymsh->GetSolutionDof(inode, iel, solType); //local 2 global solution
          unsigned idofX = mymsh->GetSolutionDof(inode, iel, 2); //local 2 global solution
          
          for(int i = 0; i < dim; i++) {
            solDOld[i][inode] = (*mysolution->_SolOld[indexSolD[i]])(idof);
            solD[i][inode] = (*mysolution->_Sol[indexSolD[i]])(idof) - solDOld[i][inode];
            
            //moving domain
            vxHat[i][inode] = (*mymsh->_topology->_Sol[i])(idofX) + solDOld[i][inode];
          }
        }
        //END
        
        
        
        for(int k = 0; k < dim; k++) {
          solD1[k].resize(nDofs);
          for(unsigned inode = 0; inode < nDofs; inode++){
            solD1[k][inode] = solD[k][inode];
          }
        }
        
        unsigned ii[9][3][3] = { 
          { {0,3,7}, {1,2,5}, {4,6,8}},
          { {1,0,4}, {2,3,6}, {5,7,8}},
          { {2,1,5}, {3,0,7}, {6,4,8}},
          { {3,2,6}, {0,1,4}, {7,5,8}} };
          
          double neumannFactor =  0. * (*mysolution->_Sol[indexSolNF])(iel);        
          if( neumannFactor > 1.0e-10 ){
                      
            for(unsigned iface = 0; iface < 4; iface++ ){
              int faceIndex = myel->GetBoundaryIndex(iel, iface);
              bool solidMark = false;
              if( ml_sol->GetIfFSI() ){
                unsigned idof = mymsh->GetSolutionDof(ii[iface][2][0], iel, solType);
                solidMark = mymsh->GetSolidMark(idof);
              }
              if( faceIndex == 1 || solidMark){
                for(unsigned inode = 0;inode < 3; inode++){
                  for(int k = 0; k < dim; k++) {
                    unsigned i0 = ii[iface][inode][0];
                    unsigned i1 = ii[iface][inode][1];
                    unsigned i2 = ii[iface][inode][2];  
                    solD1[k][i0] = neumannFactor * (- 1./3. * solD[k][i1] + 4./3 * solD[k][ i2 ] ) + (1. - neumannFactor) * solD[k][i0];
                  }
                }
                break;
              }
            }
          }
        }
        std::vector <double> xi = particles[iMarker]->GetMarkerLocalCoordinates();
      
        mymsh->_finiteElement[ielt][solType]->Jacobian(vxHat, xi, weightHat, phiHat, gradphiHat, nablaphiHat); //function to evaluate at the particles
      
        std::vector <double> solVpOld(dim);
        particles[iMarker]->GetMarkerVelocity(solVpOld);
      
        std::vector <double> solApOld(dim);
        particles[iMarker]->GetMarkerAcceleration(solApOld);
      

        
        std::vector <double> solDp(dim, 0.);
        //update displacement and acceleration
        for(int i = 0; i < dim; i++) {
          for(unsigned inode = 0; inode < nDofs; inode++) {
            solDp[i] += phiHat[inode] * solD[i][inode];
          }
        }
        
        particles[iMarker]->SetMarkerDisplacement(solDp);
        particles[iMarker]->UpdateParticleCoordinates();
        
        std::vector <double> solAp(dim);
        std::vector <double> solVp(dim);
        for(unsigned i = 0; i < dim; i++) {
          solAp[i] = 1. / (beta * dt * dt) * solDp[i] - 1. / (beta * dt) * solVpOld[i] - (1. - 2.* beta) / (2. * beta) * solApOld[i];
          solVp[i] = solVpOld[i] + dt * ((1. - Gamma) * solApOld[i] + Gamma * solAp[i]);
        }
        
        particles[iMarker]->SetMarkerVelocity(solVp);
        particles[iMarker]->SetMarkerAcceleration(solAp);
        
        //   update the deformation gradient
        
        vector< vector < double > > GradSolDpHat(dim);
  
        for(int i = 0; i < dim; i++) {
          GradSolDpHat[i].resize(dim);
        }
        
        for(int i = 0; i < dim; i++) {
          for(int j = 0; j < dim; j++) {
            GradSolDpHat[i][j] = 0.;
            for(unsigned inode = 0; inode < nDofs; inode++) {
              GradSolDpHat[i][j] +=  gradphiHat[inode * dim + j] * solD1[i][inode];
            }
          }
        }
        
        std::vector < std::vector < double > > FpOld;
        FpOld = particles[iMarker]->GetDeformationGradient(); //extraction of the deformation gradient
        
        double FpNew[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
        std::vector < std::vector < double > > Fp(dim);
        
        for(unsigned i = 0; i < dim; i++) {
          for(unsigned j = 0; j < dim; j++) {
            FpNew[i][j] += GradSolDpHat[i][j];
          }
        }
        
        for(unsigned i = 0; i < dim; i++) {
          Fp[i].resize(dim);
          for(unsigned j = 0; j < dim; j++) {
            Fp[i][j] = 0.;
            for(unsigned k = 0; k < dim; k++) {
              Fp[i][j] += FpNew[i][k] * FpOld[k][j];
            }
          }
        }
        
        particles[iMarker]->SetDeformationGradient(Fp);
        
        ielOld = iel;
    }
    else {
      break;
    }
  }
  //END loop on particles
  
  
  //BEGIN loop on elements to update grid velocity and acceleration
  for (unsigned idof = mymsh->_dofOffset[solType][iproc]; idof < mymsh->_dofOffset[solType][iproc + 1]; idof++) {
    
    if( mymsh->GetSolidMark(idof) ){
      
      for(int i = 0; i < dim; i++) {
        double disp = (*mysolution->_Sol[indexSolD[i]])(idof);
        double dispOld = (*mysolution->_SolOld[indexSolD[i]])(idof);
        double velOld = (*mysolution->_Sol[indexSolV[i]])(idof);
        double accOld = (*mysolution->_Sol[indexSolA[i]])(idof);
        double accNew = 1. / (beta * dt * dt) * (disp - dispOld) - 1. / (beta * dt) * velOld - (1. - 2.* beta) / (2. * beta) * accOld;
        double velNew = velOld + dt * ((1. - Gamma) * accOld + Gamma * accNew);
        mysolution->_Sol[indexSolV[i]]->set(idof, velNew);
        mysolution->_Sol[indexSolA[i]]->set(idof, accNew);
      }
    }
    else{
      for(int i = 0; i < dim; i++) {
        mysolution->_Sol[indexSolD[i]]->set(idof, 0.);
      }
    }
  }
  
  for(int i = 0; i < dim; i++) {
    mysolution->_Sol[indexSolD[i]]->close();
    if( ml_sol->GetIfFSI() ){    
      mysolution->_Sol[indexSolV[i]]->close();
      mysolution->_Sol[indexSolA[i]]->close();
    }
  }
  //END loop on elements to update grid velocity and acceleration
  
  linea.UpdateLineMPM();
  
  linea.GetParticlesToGridMaterial();
}


unsigned getNumberOfLayers(const double &a, const double &fac){
  double da = 1./fac; 
  double b =  da;
  unsigned n = 1;
  
  while(b < a){
    da /= fac;
    b += da;
    n++;
    if(n >= 100){
      std::cout<<"Error: number of layer is unbounded, try with a smaller factor\n";
      abort();
    }
  }
  return n;
}

