// /** started from file Ex6.cpp

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NumericVector.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"


using namespace femus;

  double force[3] = {1.,0.,0.}; 
  double forceMag[3] = {1.,2.,7.}; 

  std::string concatenate(const std::string str1, const unsigned num1) {
   
      std::ostringstream mystream("");
      
      mystream << str1 << "_" << num1;
      
      return mystream.str();
      
  }
  

bool SetBoundaryConditionBox(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  //1: bottom  //2: right  //3: top  //4: left
  
  bool dirichlet = false;
   value = 0.;
  
// LEFT ==========================  
      if (facename == 4) {
       if (!strcmp(SolName, "U_0"))    { dirichlet = false; }
  else if (!strcmp(SolName, "U_1"))    { dirichlet = true;      value = 0.; } 
      }
      
// RIGHT ==========================  
     if (facename == 2) {
       if (!strcmp(SolName, "U_0"))    {  dirichlet = false; }
  else if (!strcmp(SolName, "U_1"))    {  dirichlet = true;  value = 0.;  } 
  
  
      }
      
// BOTTOM ==========================  
      if (facename == 1) {
       if (!strcmp(SolName, "U_0"))    { dirichlet = true; value = 0.; }
  else if (!strcmp(SolName, "U_1"))    { dirichlet = true; value = 0.; } 
     
      }
      
// TOP ==========================  
      if (facename == 3) {
       if (!strcmp(SolName, "U_0"))    { dirichlet = true; value = 0.;  }
  else if (!strcmp(SolName, "U_1"))    { dirichlet = true; value = 0.;  } 
     
      }
      
      
//if the boundary integral is implemented below, you need to use these conditions on pressure to retrieve the Poiseuille flow (and set the volume force to zero)      
//       if (!strcmp(SolName, "P"))  { 
// 	 dirichlet = false;
//            if (facename == 4)  value = 1.; 
//            if (facename == 2)  value = 0.;
//    
//       }
      
  return dirichlet;
}




void AssembleNS_AD(MultiLevelProblem& ml_prob);    //, unsigned level, const unsigned &levelMax, const bool &assembleMatrix );


int main(int argc, char** args) {



  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

    // ======= Files =========================
    Files files;
    files.CheckIODirectories();
    
    const bool redirect_cout = false;
    std::string output_path;
    
    if (redirect_cout) {
        files.RedirectCout();
        output_path = files.GetOutputPath();
    }
    else {
       output_path = DEFAULT_OUTPUTDIR; 
    }

    // ======= Quad Rule ========================
    std::string fe_quad_rule("seventh");
    
  // define multilevel mesh
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  mlMsh.GenerateCoarseBoxMesh(2, 2, 0, 0., 1., 0., 1., 0., 0., QUAD9, fe_quad_rule.c_str());
//   mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
//   //mlMsh.ReadCoarseMesh ( "./input/square_quad.neu", "seventh", scalingFactor );
//   /* "seventh" is the order of accuracy that is used in the gauss integration scheme
//      probably in the furure it is not going to be an argument of this function   */
  unsigned dim = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  for (unsigned int d = 0; d < dim; d++) {
      const std::string  unknown_name = concatenate("U", d);
  mlSol.AddSolution(unknown_name.c_str(), LAGRANGE, SECOND);
  }
  
  mlSol.AddSolution("P", LAGRANGE, FIRST);

    const unsigned aux_mag_length = dim + 1;

  for (unsigned int d = 0; d < aux_mag_length; d++) {
      const std::string  unknown_name = concatenate("B", d);
       if (d <  aux_mag_length - 1) mlSol.AddSolution(unknown_name.c_str(), LAGRANGE, SECOND);
  else if (d == aux_mag_length - 1) mlSol.AddSolution(unknown_name.c_str(), LAGRANGE, FIRST);
  }

  
  mlSol.Initialize("All");

  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryConditionBox);
  mlSol.GenerateBdc("All");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);

  mlProb.SetQuadratureRuleAllGeomElems(fe_quad_rule);
  mlProb.SetFilesHandler(&files);

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("NS");

  // add solution "u" to system
  for (unsigned int d = 0; d < dim; d++) {
      const std::string  unknown_name =  concatenate("U", d);
  system.AddSolutionToSystemPDE(unknown_name.c_str());
  }  
  
  system.AddSolutionToSystemPDE("P");

  
  for (unsigned int d = 0; d < aux_mag_length; d++) {
      const std::string  unknown_name =  concatenate("B", d);
  system.AddSolutionToSystemPDE(unknown_name.c_str());
  }  
  
  
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleNS_AD);

  // initialize and solve the system
  system.init();
  system.SetOuterSolver(PREONLY);
  system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  mlSol.SetWriter(VTK);
  mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write(output_path, "biquadratic", variablesToBePrinted);

  return 0;
}




void AssembleNS_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  mlSol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*    sol        = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
  SparseMatrix*    KK         = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*   RES          = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  
    
  const unsigned aux_ns_length = dim + 1;
  const unsigned aux_mag_length = dim + 1;
  
  
  //solution variable
  vector < unsigned > solVIndex(dim);
   for (unsigned int d = 0; d < dim; d++) {
      const std::string  unknown_name =  concatenate("U", d);
  solVIndex[d] = mlSol->GetIndex(unknown_name.c_str());  // get the position of "U" in the ml_sol object
  }

  vector < unsigned > solMagIndex(aux_mag_length);
   for (unsigned int d = 0; d < solMagIndex.size(); d++) {
      const std::string  unknown_name =  concatenate("B", d);
  solMagIndex[d] = mlSol->GetIndex(unknown_name.c_str());  // get the position of "U" in the ml_sol object
  }
  

   //all unknowns ==========================
  vector < unsigned > solIndex_all(aux_ns_length + aux_mag_length);
  
  solIndex_all[0] = mlSol->GetIndex("U_0");
  solIndex_all[1] = mlSol->GetIndex("U_1");
  solIndex_all[2] = mlSol->GetIndex("P");
  solIndex_all[3] = mlSol->GetIndex("B_0");
  solIndex_all[4] = mlSol->GetIndex("B_1");
  solIndex_all[5] = mlSol->GetIndex("B_2");

  if (dim == 3) { std::cout << "Change variable inputs above"; abort(); } 
  
  vector < unsigned > solType_all(aux_ns_length + aux_mag_length);
  
   for (unsigned int d = 0; d < solType_all.size(); d++) {
       
       solType_all[d] = mlSol->GetSolutionType(solIndex_all[d]);
       
   }
   //all unknowns ==========================
  
  
  unsigned solVType = mlSol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  unsigned solPIndex = mlSol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = mlSol->GetSolutionType(solPIndex);    // get the finite element type for "u"

  std::vector < unsigned > solMagType_vec(aux_mag_length); 
  
    for (unsigned int d = 0; d < solMagType_vec.size(); d++) {
        solMagType_vec[d] = mlSol->GetSolutionType(solMagIndex[d]);
    }

  
  
  vector < unsigned > solVPdeIndex(dim);
  for (unsigned int d = 0; d < dim; d++) {
      const std::string  unknown_name =  concatenate("U", d);
  solVPdeIndex[d] = mlPdeSys->GetSolPdeIndex(unknown_name.c_str());
  }

  unsigned solPPdeIndex;
  solPPdeIndex = mlPdeSys->GetSolPdeIndex("P");    // get the position of "P" in the pdeSys object

   vector < unsigned > solMagPdeIndex(aux_mag_length);
   for (unsigned int d = 0; d < solMagPdeIndex.size(); d++) {
      const std::string  unknown_name =  concatenate("B", d);
      solMagPdeIndex[d] = mlPdeSys->GetSolPdeIndex(unknown_name.c_str());
   } 
  
  
  
  vector < vector < adept::adouble > >  solV(dim);
  vector < adept::adouble >  solP;

  vector < vector < adept::adouble > >  solMag(aux_mag_length);

  
  vector< vector < adept::adouble > > aResV(dim);    
  vector< adept::adouble > aResP;

  
  vector< vector < adept::adouble > > aResMag(aux_mag_length);    // local redidual vector
  
  
  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  for (unsigned  k = 0; k < dim; k++) {
    solV[k].reserve(maxSize);
    aResV[k].reserve(maxSize);
    coordX[k].reserve(maxSize);
  }

    for (unsigned  k = 0; k < solMag.size(); k++) {
      solMag[k].reserve(maxSize);
    }
    for (unsigned  k = 0; k < aResMag.size(); k++) {
      aResMag[k].reserve(maxSize);
    }
  
  solP.reserve(maxSize);
  aResP.reserve(maxSize);


  vector <double> phiV;  // local test function
  vector <double> phiV_x; // local test function first order partial derivatives
  vector <double> phiV_xx; // local test function second order partial derivatives

  phiV.reserve(maxSize);
  phiV_x.reserve(maxSize * dim);
  phiV_xx.reserve(maxSize * dim2);

  double* phiP;
  
  
  vector< vector <double> > phiMag(aux_mag_length);  // local test function
  vector< vector <double> > phiMag_x(aux_mag_length); // local test function first order partial derivatives
    for (unsigned  k = 0; k < aux_mag_length; k++) {
       phiMag[k].reserve(maxSize);
       phiMag_x[k].reserve(maxSize * dim);
    }
  
  
  
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve(2 * aux_mag_length * maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve(2 * aux_mag_length  * maxSize);

  vector < double > Jac;
  Jac.reserve(2 * aux_mag_length * maxSize * 2 * aux_mag_length * maxSize);

  KK->zero(); // Set to zero all the entries of the Global Matrix


  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    short unsigned ielGeom = msh->GetElementType(iel);

    unsigned nDofsV = msh->GetElementDofNumber(iel, solVType);    // number of solution element dofs
    unsigned nDofsP = msh->GetElementDofNumber(iel, solPType);    // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs
    

    std::vector < unsigned int > nDofs_all(aux_ns_length + aux_mag_length);
       for (unsigned int d = 0; d < nDofs_all.size(); d++) {
          nDofs_all[d] = msh->GetElementDofNumber(iel, solType_all[d]);
       }
       
    std::vector < unsigned int > nDofsMag_vec(aux_mag_length);
       for (unsigned int d = 0; d < nDofsMag_vec.size(); d++) {
          nDofsMag_vec[d] = msh->GetElementDofNumber(iel, solMagType_vec[d]);
       }
       
       
      std::vector < unsigned int > nDofsMag_vec_offset(aux_mag_length + 1);
     
      nDofsMag_vec_offset[0] = 0;
       for (unsigned int d = 1; d < nDofsMag_vec_offset.size(); d++)  {
           nDofsMag_vec_offset[d] = nDofsMag_vec_offset[d-1] + nDofsMag_vec[d-1];
       }
      
    
    unsigned nDofsVP = dim * nDofsV + nDofsP;
    
    unsigned nDofsMag = 0;
       for (unsigned int d = 0; d < nDofsMag_vec.size(); d++)  nDofsMag += nDofsMag_vec[d];

    unsigned nDofsTot = nDofsVP + nDofsMag;
       
       
    // resize local arrays
    sysDof.resize(nDofsTot);

    for (unsigned  k = 0; k < dim; k++) {
      solV[k].resize(nDofsV);
      coordX[k].resize(nDofsX);
    }

    solP.resize(nDofsP);

        for (unsigned  k = 0; k < solMag.size(); k++) {
       solMag[k].resize(nDofsMag_vec[k]);
       }
    
    for (unsigned  k = 0; k < dim; k++) {
      aResV[k].resize(nDofsV);    //resize
      std::fill(aResV[k].begin(), aResV[k].end(), 0.);    //set aRes to zero
    }

    aResP.resize(nDofsP);    //resize
    std::fill(aResP.begin(), aResP.end(), 0.);    //set aRes to zero
    
    for (unsigned  k = 0; k < aResMag.size(); k++) {
      aResMag[k].resize( nDofsMag_vec[k] );    //resize
      std::fill(aResMag[k].begin(), aResMag[k].end(), 0.);    //set aRes to zero
    }
    
    
    
    // local storage of global mapping and solution *****************************
    for (unsigned i = 0; i < nDofsV; i++) {
      unsigned solVDof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < dim; k++) {
        solV[k][i] = (*sol->_Sol[solVIndex[k]])(solVDof);      // global extraction and local storage for the solution
        sysDof[i + k * nDofsV] = pdeSys->GetSystemDof(solVIndex[k], solVPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }

    const unsigned press_offset = dim * nDofsV;
    
    for (unsigned i = 0; i < nDofsP; i++) {
      unsigned solPDof = msh->GetSolutionDof(i, iel, solPType);    // global to global mapping between solution node and solution dof
      solP[i] = (*sol->_Sol[solPIndex])(solPDof);      // global extraction and local storage for the solution
      sysDof[i + press_offset] = pdeSys->GetSystemDof(solPIndex, solPPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dof
    }
    
    
    const unsigned mag_offset = dim * nDofsV + nDofsP;
    
      for (unsigned  k = 0; k < solMag.size(); k++) {

          for (unsigned i = 0; i < nDofsMag_vec[k]; i++) {
      unsigned solMagDof = msh->GetSolutionDof(i, iel, solMagType_vec[k]);    // global to global mapping between solution node and solution dof

        solMag[k][i] = (*sol->_Sol[solMagIndex[k]])(solMagDof);      // global extraction and local storage for the solution
        sysDof[mag_offset + i + nDofsMag_vec_offset[k] ] = pdeSys->GetSystemDof(solMagIndex[k], solMagPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }   
    // local storage of global mapping and solution *****************************
    
    
    

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
        
    for (unsigned  k = 0; k < aux_mag_length; k++) {
      msh->_finiteElement[ielGeom][solMagType_vec[k]]->Jacobian(coordX, ig, weight, phiMag[k], phiMag_x[k], boost::none);
    }
    
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);
      
    
      vector < adept::adouble > solV_gss(dim, 0.);
      vector < vector < adept::adouble > > gradSolV_gss(dim);

      for (unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].resize(dim);
        std::fill(gradSolV_gss[k].begin(), gradSolV_gss[k].end(), 0.);
      }

      for (unsigned i = 0; i < nDofsV; i++) {
        for (unsigned  k = 0; k < dim; k++) {
          solV_gss[k] += phiV[i] * solV[k][i];
        }

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            gradSolV_gss[k][j] += phiV_x[i * dim + j] * solV[k][i];
          }
        }
      }

      adept::adouble solP_gss = 0;

      for (unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }

      
        vector < adept::adouble > solMag_qp(aux_mag_length, 0.);
      for (unsigned  k = 0; k < aux_mag_length; k++) {
          solMag_qp[k] = 0.;
        for (unsigned i = 0; i < nDofsMag_vec[k]; i++) {
          solMag_qp[k] += phiMag[k][i] * solMag[k][i];
        }
      
      } //this was missing!!
      
      double nu = 1.;

      // *** phiV_i loop ***
      for (unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > NSV(dim, 0.);

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            NSV[k]   +=  nu * phiV_x[i * dim + j] * (gradSolV_gss[k][j] + gradSolV_gss[j][k]);  //diffusion
            NSV[k]   +=  phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]);                                //advection
          }
        }

        for (unsigned  k = 0; k < dim; k++) {
          NSV[k] += - solP_gss * phiV_x[i * dim + k];   //pressure gradient
        }

        for (unsigned  k = 0; k < dim; k++) {
          aResV[k][i] += ( force[k] * phiV[i] - NSV[k] ) * weight;
        }
        
        
   
      } // end phiV_i loop

      // *** phiP_i loop ***
      for (unsigned i = 0; i < nDofsP; i++) {
        for (int k = 0; k < dim; k++) {
          aResP[i] += - (gradSolV_gss[k][k]) * phiP[i]  * weight;
        }
      } // end phiP_i loop
      
      for (unsigned  k = 0; k < aux_mag_length; k++) {
         for (unsigned i = 0; i < nDofsMag_vec[k]; i++) {
          aResMag[k][i] += (  forceMag[k] * phiMag[k][i] - phiMag[k][i] * solMag_qp[k] ) * weight;
        }   
       }


      

    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsTot);    //resize

    for (int i = 0; i < nDofsV; i++) {
      for (unsigned  k = 0; k < dim; k++) {
        Res[ i +  k * nDofsV ] = - aResV[k][i].value();
      }
    }

    for (int i = 0; i < nDofsP; i++) {
      Res[ i + press_offset ] = -aResP[i].value();

    }
    
    for (unsigned  k = 0; k < aux_mag_length; k++) {
      for (int i = 0; i < aResMag[k].size(); i++) {
        Res[ mag_offset + i +  nDofsMag_vec_offset[k] ] = - aResMag[k][i].value();
      }
    }


    RES->add_vector_blocked(Res, sysDof);

    //Extarct and store the Jacobian

    Jac.resize(nDofsTot * nDofsTot);
    // define the dependent variables

    // how to fill element Jacobian with AD - dependent variables ****************
    for (unsigned  k = 0; k < dim; k++) {
      s.dependent(&aResV[k][0], nDofsV);
    }

    s.dependent(&aResP[0], nDofsP);

    for (unsigned  k = 0; k < aResMag.size(); k++) {
      s.dependent(&aResMag[k][0], nDofsMag_vec[k]);
    }

    // how to fill element Jacobian with AD - dependent variables - end ****************
    
    // how to fill element Jacobian with AD - independent variables ****************
    for (unsigned  k = 0; k < dim; k++) {
      s.independent(&solV[k][0], nDofsV);
    }

    s.independent(&solP[0], nDofsP);

    for (unsigned  k = 0; k < aux_mag_length; k++) {
      s.independent(&solMag[k][0], nDofsMag_vec[k]);
    }
    // how to fill element Jacobian with AD - independent variables - end ****************

    // get the and store jacobian matrix (row-major)
    s.jacobian(&Jac[0] , true);
    
    KK->add_matrix_blocked(Jac, sysDof, sysDof);

    
    s.clear_independents();
    s.clear_dependents();

//          assemble_jacobian<double,double>::print_element_residual(iel, Res, nDofs_all, 10, 5);
//          assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, nDofs_all, 10, 5);
   
    
    
  } //end element loop for each process

  RES->close();

  KK->close();

  // ***************** END ASSEMBLY *******************
}

