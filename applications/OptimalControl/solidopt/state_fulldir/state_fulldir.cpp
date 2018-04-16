// /** started from file Ex6.cpp

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "LinearImplicitSystem.hpp"
#include "Parameter.hpp"
#include "Solid.hpp"
#include "Files.hpp"


// #define PRESS 1


using namespace femus;

  double force[3] = {0.,0.,0.}; 
  
  
double InitialValueDX(const std::vector < double >& x) {
  return 0.;
}

double InitialValueDY(const std::vector < double >& x) {
  return 0.;
}

double InitialValueDZ(const std::vector < double >& x) {
  return 0.;
}

double InitialValueP(const std::vector < double >& x) {
  return 0.;
}  
  
  
  
  

bool SetBoundaryConditionBox(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  //1: bottom  //2: right  //3: top  //4: left
  
  bool dirichlet = true;
   value = 0.;
  
// TOP ==========================  
      if (facename == 3) {
       if (!strcmp(SolName, "DX"))    { value = 70.; } //lid - driven
  else if (!strcmp(SolName, "DY"))    { value = 0.; } 
  	
      }
      
  return dirichlet;
}




void AssembleSolidMech_AD(MultiLevelProblem& ml_prob);

void AssembleSolidMech_nonAD(MultiLevelProblem& ml_prob);  //to double check in the future



int main(int argc, char** args) {



  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  Files files;
        files.CheckIODirectories();
        files.RedirectCout();


  MultiLevelMesh mlMsh;  // define multilevel mesh
  double scalingFactor = 1.;  // read coarse level mesh and generate finers level meshes
  
    //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;
 // *** apparently needed by non-AD assemble only **********************
  // add fluid material
  Parameter par(Lref,Uref);
  
   // Generate Solid Object
  double E = 1500000;
  double ni = 0.5;
  double rhos = 1000;
  Solid solid;
  solid = Solid(par,E,ni,rhos,"Mooney-Rivlin");

  std::cout << "Solid properties: " << std::endl;
  std::cout << solid << std::endl;

  
  mlMsh.GenerateCoarseBoxMesh(2,2,0,0.,1.,0.,1.,0.,0.,QUAD9,"seventh");
//   mlMsh.ReadCoarseMesh("./input/cube_hex.neu", "seventh", scalingFactor);
//   //mlMsh.ReadCoarseMesh ( "./input/square_quad.neu", "seventh", scalingFactor );
//   /* "seventh" is the order of accuracy that is used in the gauss integration scheme
//      probably in the furure it is not going to be an argument of this function   */
  unsigned dimension = mlMsh.GetDimension();

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  MultiLevelSolution ml_sol(&mlMsh);

  // add variables to ml_sol
  ml_sol.AddSolution("DX",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("DY",LAGRANGE,SECOND,1);
  if ( dimension == 3 ) ml_sol.AddSolution("DZ",LAGRANGE,SECOND,1);
  ml_sol.AddSolution("P",DISCONTINOUS_POLYNOMIAL,FIRST,1);

  ml_sol.Initialize("All");
  ml_sol.Initialize("DX", InitialValueDX);
  ml_sol.Initialize("DY", InitialValueDY);
  ml_sol.Initialize("DZ", InitialValueDZ);
  ml_sol.Initialize("P", InitialValueP);

  // attach the boundary condition function and generate boundary data
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionBox);
  ml_sol.GenerateBdc("All");

  // define the multilevel problem attach the ml_sol object to it
  MultiLevelProblem mlProb(&ml_sol);

  mlProb.parameters.set<Solid>("Solid") = solid;

  // add system Poisson in mlProb as a Linear Implicit System
  NonLinearImplicitSystem& system = mlProb.add_system < NonLinearImplicitSystem > ("SolidMech");

  // add solution "u" to system
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");
  if (dimension == 3) system.AddSolutionToSystemPDE("DZ");
  system.AddSolutionToSystemPDE("P");


  // attach the assembling function to system
  system.SetAssembleFunction(AssembleSolidMech_AD);

  // initilaize and solve the system
  system.init();
  
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  system.MLsolve();
//   system.MGsolve();

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);
  ml_sol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/,"biquadratic", variablesToBePrinted);
 
  //Destroy all the new systems
  mlProb.clear();
  
  return 0;
}




void AssembleSolidMech_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;

  //  extract pointers to the several objects that we are going to use
  NonLinearImplicitSystem& mlPdeSys   = ml_prob.get_system<NonLinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys.GetLevelToAssemble();
  bool assembleMatrix = mlPdeSys.GetAssembleMatrix(); 

  Mesh*          msh          = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*          el         = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*  ml_sol        = ml_prob._ml_sol;  // pointer to the multilevel solution object
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
  solVIndex[0] = ml_sol->GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = ml_sol->GetIndex("V");    // get the position of "V" in the ml_sol object

  if (dim == 3) solVIndex[2] = ml_sol->GetIndex("W");      // get the position of "V" in the ml_sol object

  unsigned solVType = ml_sol->GetSolutionType(solVIndex[0]);    // get the finite element type for "u"


  vector < unsigned > solVPdeIndex(dim);
  solVPdeIndex[0] = mlPdeSys.GetSolPdeIndex("U");    // get the position of "U" in the pdeSys object
  solVPdeIndex[1] = mlPdeSys.GetSolPdeIndex("V");    // get the position of "V" in the pdeSys object

  if (dim == 3) solVPdeIndex[2] = mlPdeSys.GetSolPdeIndex("W");

// #if PRESS == 1
  unsigned solPIndex;
  solPIndex = ml_sol->GetIndex("P");    // get the position of "P" in the ml_sol object
  unsigned solPType = ml_sol->GetSolutionType(solPIndex);    // get the finite element type for "u"
  
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
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
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
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solVType]->GetGaussPointNumber(); ig++) {
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);
// #if PRESS == 1
      phiP = msh->_finiteElement[ielGeom][solPType]->GetPhi(ig);
// #endif
      
      vector < adept::adouble > solV_gss(dim, 0);
      vector < vector < adept::adouble > > gradSolV_gss(dim);

      for (unsigned  k = 0; k < dim; k++) {
        gradSolV_gss[k].resize(dim);
        std::fill(gradSolV_gss[k].begin(), gradSolV_gss[k].end(), 0);
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

// #if PRESS == 1
      adept::adouble solP_gss = 0;
      for (unsigned i = 0; i < nDofsP; i++) {
        solP_gss += phiP[i] * solP[i];
      }
// #endif

      // *** phiV_i loop ***
      for (unsigned i = 0; i < nDofsV; i++) {
        vector < adept::adouble > NSV(dim, 0.);

        for (unsigned j = 0; j < dim; j++) {
          for (unsigned  k = 0; k < dim; k++) {
            NSV[k]   += /* nu*/  phiV_x[i * dim + j] * (gradSolV_gss[k][j] /*+ gradSolV_gss[j][k]*/);
//             NSV[k]   +=  phiV[i] * (solV_gss[j] * gradSolV_gss[k][j]);
          }
        }

// #if PRESS == 1
        for (unsigned  k = 0; k < dim; k++) {
          NSV[k] += -solP_gss * phiV_x[i * dim + k];
        }
// #endif

        for (unsigned  k = 0; k < dim; k++) {
          aResV[k][i] += ( force[k] * phiV[i] - NSV[k] ) * weight;
        }
      } // end phiV_i loop

// #if PRESS == 1
      // *** phiP_i loop ***
      for (unsigned i = 0; i < nDofsP; i++) {
        for (int k = 0; k < dim; k++) {
          aResP[i] +=  (gradSolV_gss[k][k]) * phiP[i]  * weight;
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
        Res[ i +  k * nDofsV ] = -aResV[k][i].value();
      }
    }

// #if PRESS == 1
    for (int i = 0; i < nDofsP; i++) {
      Res[ i + dim * nDofsV ] = -aResP[i].value();
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
//   std::ostringstream mat_out; mat_out << "matrix_ad" << mlPdeSys._nonliniteration  << ".txt";
//   KK->print_matlab(mat_out.str(),"ascii");
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


