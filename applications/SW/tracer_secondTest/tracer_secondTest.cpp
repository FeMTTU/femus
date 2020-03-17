/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution varables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "PetscMatrix.hpp"

#include "TransientSystem.hpp"
#include "LinearImplicitSystem.hpp"

#include "slepceps.h"
#include <slepcmfn.h>
#include <slepcsvd.h>

using namespace femus;

double dt = 1.;

double k_v = (2.5) * (0.00001);

double pi = acos (-1.);
//double k_h = 1. / ( 10 * pi );
double k_h = 0.0001 ;

const unsigned NumberOfLayers = 4;
unsigned RK_order = 4;

unsigned counter = 0;
unsigned counter2 = 0;

clock_t start_time = clock();

bool wave = false;
bool twostage = false;
bool splitting = true;
bool assembly = true; //assembly must be left always true

const double hRest[4] = {2.5, 2.5, 2.5, 2.5};

//const double hRest[20] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

// const double hRest[40] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
//                           0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25
//                          };

// const double hRest[80] = {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
//                           0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
//                           0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
//                           0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125
//                          };

double InitalValueVi (const std::vector < double >& x , const unsigned &i)
{
  double psi1 = 1. - (x[0] - 5.) * (x[0] - 5.) * (x[0] - 5.) * (x[0] - 5.) / (5.*5.*5.*5.);
  //double psi1 = 1. - ( pow(( x[0] - 50. ), 40) / pow(50., 40) );
  //double psi1 = 1. - ( pow(( x[0] - 10. ), 16) / pow(10., 16) );
  double z = -10. + hRest[0] / 2. + hRest[0] * (NumberOfLayers - i);
  double d_psi2 = (- (2.*z + 10.)) / (5 * 5);
  return (psi1 * d_psi2);
}

//  double InitalValueV0 ( const std::vector < double >& x ) {
//   //double psi1 = ( 10. - x[0] ) * x[0] / ( 5.*5. );
//   double psi1 = 1. - ( x[0] - 5. ) * ( x[0] - 5. ) * ( x[0] - 5. ) * ( x[0] - 5. ) / ( 5.*5.*5.*5. );
//   double z = -10 + hRest[0] / 2 + hRest[0] * ( NumberOfLayers - 1 );
//   double d_psi2 = ( - ( 2.*z + 10. ) ) / ( 5 * 5 );
//   //return ( -psi1 * d_psi2 );
//   return ( psi1 * d_psi2 );
// }



double InitalValueV0 (const std::vector < double >& x)
{
  return InitalValueVi (x, 1);
}
double InitalValueV1 (const std::vector < double >& x)
{
  return InitalValueVi (x, 2);
}
double InitalValueV2 (const std::vector < double >& x)
{
  return InitalValueVi (x, 3);
}
double InitalValueV3 (const std::vector < double >& x)
{
  return InitalValueVi (x, 4);
}
double InitalValueV4 (const std::vector < double >& x)
{
  return InitalValueVi (x, 5);
}
double InitalValueV5 (const std::vector < double >& x)
{
  return InitalValueVi (x, 6);
}
double InitalValueV6 (const std::vector < double >& x)
{
  return InitalValueVi (x, 7);
}
double InitalValueV7 (const std::vector < double >& x)
{
  return InitalValueVi (x, 8);
}
double InitalValueV8 (const std::vector < double >& x)
{
  return InitalValueVi (x, 9);
}
double InitalValueV9 (const std::vector < double >& x)
{
  return InitalValueVi (x, 10);
}
double InitalValueV10 (const std::vector < double >& x)
{
  return InitalValueVi (x, 11);
}
double InitalValueV11 (const std::vector < double >& x)
{
  return InitalValueVi (x, 12);
}
double InitalValueV12 (const std::vector < double >& x)
{
  return InitalValueVi (x, 13);
}
double InitalValueV13 (const std::vector < double >& x)
{
  return InitalValueVi (x, 14);
}
double InitalValueV14 (const std::vector < double >& x)
{
  return InitalValueVi (x, 15);
}
double InitalValueV15 (const std::vector < double >& x)
{
  return InitalValueVi (x, 16);
}
double InitalValueV16 (const std::vector < double >& x)
{
  return InitalValueVi (x, 17);
}
double InitalValueV17 (const std::vector < double >& x)
{
  return InitalValueVi (x, 18);
}
double InitalValueV18 (const std::vector < double >& x)
{
  return InitalValueVi (x, 19);
}
double InitalValueV19 (const std::vector < double >& x)
{
  return InitalValueVi (x, 20);
}
double InitalValueV20 (const std::vector < double >& x)
{
  return InitalValueVi (x, 21);
}
double InitalValueV21 (const std::vector < double >& x)
{
  return InitalValueVi (x, 22);
}
double InitalValueV22 (const std::vector < double >& x)
{
  return InitalValueVi (x, 23);
}
double InitalValueV23 (const std::vector < double >& x)
{
  return InitalValueVi (x, 24);
}
double InitalValueV24 (const std::vector < double >& x)
{
  return InitalValueVi (x, 25);
}
double InitalValueV25 (const std::vector < double >& x)
{
  return InitalValueVi (x, 26);
}
double InitalValueV26 (const std::vector < double >& x)
{
  return InitalValueVi (x, 27);
}
double InitalValueV27 (const std::vector < double >& x)
{
  return InitalValueVi (x, 28);
}
double InitalValueV28 (const std::vector < double >& x)
{
  return InitalValueVi (x, 29);
}
double InitalValueV29 (const std::vector < double >& x)
{
  return InitalValueVi (x, 30);
}
double InitalValueV30 (const std::vector < double >& x)
{
  return InitalValueVi (x, 31);
}
double InitalValueV31 (const std::vector < double >& x)
{
  return InitalValueVi (x, 32);
}
double InitalValueV32 (const std::vector < double >& x)
{
  return InitalValueVi (x, 33);
}
double InitalValueV33 (const std::vector < double >& x)
{
  return InitalValueVi (x, 34);
}
double InitalValueV34 (const std::vector < double >& x)
{
  return InitalValueVi (x, 35);
}
double InitalValueV35 (const std::vector < double >& x)
{
  return InitalValueVi (x, 36);
}
double InitalValueV36 (const std::vector < double >& x)
{
  return InitalValueVi (x, 37);
}
double InitalValueV37 (const std::vector < double >& x)
{
  return InitalValueVi (x, 38);
}
double InitalValueV38 (const std::vector < double >& x)
{
  return InitalValueVi (x, 39);
}
double InitalValueV39 (const std::vector < double >& x)
{
  return InitalValueVi (x, 40);
}
double InitalValueV40 (const std::vector < double >& x)
{
  return InitalValueVi (x, 41);
}
double InitalValueV41 (const std::vector < double >& x)
{
  return InitalValueVi (x, 42);
}
double InitalValueV42 (const std::vector < double >& x)
{
  return InitalValueVi (x, 43);
}
double InitalValueV43 (const std::vector < double >& x)
{
  return InitalValueVi (x, 44);
}
double InitalValueV44 (const std::vector < double >& x)
{
  return InitalValueVi (x, 45);
}
double InitalValueV45 (const std::vector < double >& x)
{
  return InitalValueVi (x, 46);
}
double InitalValueV46 (const std::vector < double >& x)
{
  return InitalValueVi (x, 47);
}
double InitalValueV47 (const std::vector < double >& x)
{
  return InitalValueVi (x, 48);
}
double InitalValueV48 (const std::vector < double >& x)
{
  return InitalValueVi (x, 49);
}
double InitalValueV49 (const std::vector < double >& x)
{
  return InitalValueVi (x, 50);
}
double InitalValueV50 (const std::vector < double >& x)
{
  return InitalValueVi (x, 51);
}
double InitalValueV51 (const std::vector < double >& x)
{
  return InitalValueVi (x, 52);
}
double InitalValueV52 (const std::vector < double >& x)
{
  return InitalValueVi (x, 53);
}
double InitalValueV53 (const std::vector < double >& x)
{
  return InitalValueVi (x, 54);
}
double InitalValueV54 (const std::vector < double >& x)
{
  return InitalValueVi (x, 55);
}
double InitalValueV55 (const std::vector < double >& x)
{
  return InitalValueVi (x, 56);
}
double InitalValueV56 (const std::vector < double >& x)
{
  return InitalValueVi (x, 57);
}
double InitalValueV57 (const std::vector < double >& x)
{
  return InitalValueVi (x, 58);
}
double InitalValueV58 (const std::vector < double >& x)
{
  return InitalValueVi (x, 59);
}
double InitalValueV59 (const std::vector < double >& x)
{
  return InitalValueVi (x, 60);
}
double InitalValueV60 (const std::vector < double >& x)
{
  return InitalValueVi (x, 61);
}
double InitalValueV61 (const std::vector < double >& x)
{
  return InitalValueVi (x, 62);
}
double InitalValueV62 (const std::vector < double >& x)
{
  return InitalValueVi (x, 63);
}
double InitalValueV63 (const std::vector < double >& x)
{
  return InitalValueVi (x, 64);
}
double InitalValueV64 (const std::vector < double >& x)
{
  return InitalValueVi (x, 65);
}
double InitalValueV65 (const std::vector < double >& x)
{
  return InitalValueVi (x, 66);
}
double InitalValueV66 (const std::vector < double >& x)
{
  return InitalValueVi (x, 67);
}
double InitalValueV67 (const std::vector < double >& x)
{
  return InitalValueVi (x, 68);
}
double InitalValueV68 (const std::vector < double >& x)
{
  return InitalValueVi (x, 69);
}
double InitalValueV69 (const std::vector < double >& x)
{
  return InitalValueVi (x, 70);
}
double InitalValueV70 (const std::vector < double >& x)
{
  return InitalValueVi (x, 71);
}
double InitalValueV71 (const std::vector < double >& x)
{
  return InitalValueVi (x, 72);
}
double InitalValueV72 (const std::vector < double >& x)
{
  return InitalValueVi (x, 73);
}
double InitalValueV73 (const std::vector < double >& x)
{
  return InitalValueVi (x, 74);
}
double InitalValueV74 (const std::vector < double >& x)
{
  return InitalValueVi (x, 75);
}
double InitalValueV75 (const std::vector < double >& x)
{
  return InitalValueVi (x, 76);
}
double InitalValueV76 (const std::vector < double >& x)
{
  return InitalValueVi (x, 77);
}
double InitalValueV77 (const std::vector < double >& x)
{
  return InitalValueVi (x, 78);
}
double InitalValueV78 (const std::vector < double >& x)
{
  return InitalValueVi (x, 79);
}
double InitalValueV79 (const std::vector < double >& x)
{
  return InitalValueVi (x, 80);
}



double InitalValueH (const std::vector < double >& x)
{
  return hRest[0];
}

double InitalValueT (const std::vector < double >& x)
{
  double pi = acos (-1.);

  if (x[0] < 5) {
    return 5.;

  }

  else {
    return 30.;
  }

//   return (sin(pi*x[0]));
//  return 10.;
}


double InitalValueB (const std::vector < double >& x)
{
  return 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );
}

double InitialValueK_RK (const std::vector < double >& x)
{
  return 0.;
}

bool SetBoundaryCondition (const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time)
{
  bool dirichlet = false;

  if (!strcmp (SolName, "HT")) {
    if (facename == 1 || facename == 2) {
      dirichlet = true;
      value = 0.;
    }
  }

  if (!strcmp (SolName, "T")) {
    if (facename == 1 || facename == 2) {
      dirichlet = true;
      value = 0.;
    }
  }

  return dirichlet;
}


void ETD (MultiLevelProblem& ml_prob, const unsigned & numberOfTimeSteps);

void RK4 (MultiLevelProblem& ml_prob, const bool & implicitEuler, const unsigned & numberOfTimeSteps);

void RK (MultiLevelProblem& ml_prob, const unsigned & numberOfTimeSteps);

void RKe (MultiLevelProblem& ml_prob, const unsigned & numberOfTimeSteps);

int main (int argc, char** args)
{

  SlepcInitialize (&argc, &args, PETSC_NULL, PETSC_NULL);
 
  // init Petsc-MPI communicator
  FemusInit mpinit (argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;

  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;

  unsigned nx = static_cast<unsigned> (floor (pow (2.,/*8*/ 4) + 0.5));       //Grid cell size = 3.90625 m
  nx += 3;
//     nx*=4;
//     std::cout <<" nx = " << nx << std::endl;

  double length = 10.; //2 * 1465700.;

  //mlMsh.GenerateCoarseBoxMesh ( nx, 0, 0, -length / 2, length / 2, 0., 0., 0., 0., EDGE3, "seventh" );
  mlMsh.GenerateCoarseBoxMesh (nx, 0, 0, 0, length, 0., 0., 0., 0., EDGE3, "seventh");

  //mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol (&mlMsh);

  for (unsigned i = 0; i < NumberOfLayers; i++) {
    char name[10];
    sprintf (name, "h%d", i);
    mlSol.AddSolution (name, DISCONTINUOUS_POLYNOMIAL, ZERO, 2);
    sprintf (name, "v%d", i);
    mlSol.AddSolution (name, LAGRANGE, FIRST, 2);
    sprintf (name, "T%d", i);
    mlSol.AddSolution (name, DISCONTINUOUS_POLYNOMIAL, ZERO, 2);
    sprintf (name, "HT%d", i);
    mlSol.AddSolution (name, DISCONTINUOUS_POLYNOMIAL, ZERO, 2);
    
    for(unsigned j = 0; j<RK_order;j++){
      sprintf (name, "K%d_%d",j+1, i);
      mlSol.AddSolution (name, DISCONTINUOUS_POLYNOMIAL, ZERO);
    }
  }

  mlSol.AddSolution ("b", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false);

  mlSol.AddSolution ("eta", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, false);

  mlSol.Initialize ("All");

  mlSol.Initialize ("v0", InitalValueV0);
  mlSol.Initialize ("v1", InitalValueV1);
  mlSol.Initialize ("v2", InitalValueV2);
  mlSol.Initialize ("v3", InitalValueV3);
//     mlSol.Initialize ( "v4", InitalValueV4 );
//     mlSol.Initialize ( "v5", InitalValueV5 );
//     mlSol.Initialize ( "v6", InitalValueV6 );
//     mlSol.Initialize ( "v7", InitalValueV7 );
//     mlSol.Initialize ( "v8", InitalValueV8 );
//     mlSol.Initialize ( "v9", InitalValueV9 );
//     mlSol.Initialize ( "v10", InitalValueV10 );
//     mlSol.Initialize ( "v11", InitalValueV11 );
//     mlSol.Initialize ( "v12", InitalValueV12 );
//     mlSol.Initialize ( "v13", InitalValueV13 );
//     mlSol.Initialize ( "v14", InitalValueV14 );
//     mlSol.Initialize ( "v15", InitalValueV15 );
//     mlSol.Initialize ( "v16", InitalValueV16 );
//     mlSol.Initialize ( "v17", InitalValueV17 );
//     mlSol.Initialize ( "v18", InitalValueV18 );
//     mlSol.Initialize ( "v19", InitalValueV19 );
//     if ( NumberOfLayers>39 ) {
//         mlSol.Initialize ( "v20", InitalValueV20 );
//         mlSol.Initialize ( "v21", InitalValueV21 );
//         mlSol.Initialize ( "v22", InitalValueV22 );
//         mlSol.Initialize ( "v23", InitalValueV23 );
//         mlSol.Initialize ( "v24", InitalValueV24 );
//         mlSol.Initialize ( "v25", InitalValueV25 );
//         mlSol.Initialize ( "v26", InitalValueV26 );
//         mlSol.Initialize ( "v27", InitalValueV27 );
//         mlSol.Initialize ( "v28", InitalValueV28 );
//         mlSol.Initialize ( "v29", InitalValueV29 );
//         mlSol.Initialize ( "v30", InitalValueV30 );
//         mlSol.Initialize ( "v31", InitalValueV31 );
//         mlSol.Initialize ( "v32", InitalValueV32 );
//         mlSol.Initialize ( "v33", InitalValueV33 );
//         mlSol.Initialize ( "v34", InitalValueV34 );
//         mlSol.Initialize ( "v35", InitalValueV35 );
//         mlSol.Initialize ( "v36", InitalValueV36 );
//         mlSol.Initialize ( "v37", InitalValueV37 );
//         mlSol.Initialize ( "v38", InitalValueV38 );
//         mlSol.Initialize ( "v39", InitalValueV39 );
//         if ( NumberOfLayers>79 ) {
//             mlSol.Initialize ( "v40", InitalValueV40 );
//             mlSol.Initialize ( "v41", InitalValueV41 );
//             mlSol.Initialize ( "v42", InitalValueV42 );
//             mlSol.Initialize ( "v43", InitalValueV43 );
//             mlSol.Initialize ( "v44", InitalValueV44 );
//             mlSol.Initialize ( "v45", InitalValueV45 );
//             mlSol.Initialize ( "v46", InitalValueV46 );
//             mlSol.Initialize ( "v47", InitalValueV47 );
//             mlSol.Initialize ( "v48", InitalValueV48 );
//             mlSol.Initialize ( "v49", InitalValueV49 );
//             mlSol.Initialize ( "v50", InitalValueV50 );
//             mlSol.Initialize ( "v51", InitalValueV51 );
//             mlSol.Initialize ( "v52", InitalValueV52 );
//             mlSol.Initialize ( "v53", InitalValueV53 );
//             mlSol.Initialize ( "v54", InitalValueV54 );
//             mlSol.Initialize ( "v55", InitalValueV55 );
//             mlSol.Initialize ( "v56", InitalValueV56 );
//             mlSol.Initialize ( "v57", InitalValueV57 );
//             mlSol.Initialize ( "v58", InitalValueV58 );
//             mlSol.Initialize ( "v59", InitalValueV59 );
//             mlSol.Initialize ( "v60", InitalValueV60 );
//             mlSol.Initialize ( "v61", InitalValueV61 );
//             mlSol.Initialize ( "v62", InitalValueV62 );
//             mlSol.Initialize ( "v63", InitalValueV63 );
//             mlSol.Initialize ( "v64", InitalValueV64 );
//             mlSol.Initialize ( "v65", InitalValueV65 );
//             mlSol.Initialize ( "v66", InitalValueV66 );
//             mlSol.Initialize ( "v67", InitalValueV67 );
//             mlSol.Initialize ( "v68", InitalValueV68 );
//             mlSol.Initialize ( "v69", InitalValueV69 );
//             mlSol.Initialize ( "v70", InitalValueV70 );
//             mlSol.Initialize ( "v71", InitalValueV71 );
//             mlSol.Initialize ( "v72", InitalValueV72 );
//             mlSol.Initialize ( "v73", InitalValueV73 );
//             mlSol.Initialize ( "v74", InitalValueV74 );
//             mlSol.Initialize ( "v75", InitalValueV75 );
//             mlSol.Initialize ( "v76", InitalValueV76 );
//             mlSol.Initialize ( "v77", InitalValueV77 );
//             mlSol.Initialize ( "v78", InitalValueV78 );
//             mlSol.Initialize ( "v79", InitalValueV79 );
//         }
//     }


  for (unsigned i = 0; i < NumberOfLayers; i++) {
    char name[10];
    sprintf (name, "h%d", i);
    mlSol.Initialize (name, InitalValueH);
  }

  for (unsigned i = 0; i < NumberOfLayers; i++) {
    char name[10];
    sprintf (name, "T%d", i);
    mlSol.Initialize (name, InitalValueT);
  }

  mlSol.Initialize ("b", InitalValueB);

  mlSol.AttachSetBoundaryConditionFunction (SetBoundaryCondition);
  mlSol.GenerateBdc ("All");

  MultiLevelProblem ml_prob (&mlSol);

  // ******* Add FEM system to the MultiLevel problem *******
  TransientLinearImplicitSystem& system = ml_prob.add_system < TransientLinearImplicitSystem > ("SWt");

  for (unsigned i = 0; i < NumberOfLayers; i++) {
    char name[10];
    sprintf (name, "HT%d", i);
    system.AddSolutionToSystemPDE (name);
  }

  system.init();

  mlSol.SetWriter (VTK);
  std::vector<std::string> print_vars;
  print_vars.push_back ("All");
  //mlSol.GetWriter()->SetDebugOutput(true);
  mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "linear", print_vars, 0);

  unsigned numberOfTimeSteps = 10000; //RK4: dt=0.5, numberOfTimeSteps = 16001
  dt = 0.5;
  bool implicitEuler = true;

  for (unsigned i = 0; i < numberOfTimeSteps; i++) {
    if (wave == true) {
      assembly = (i == 0) ? true : false;
    }

    system.CopySolutionToOldSolution();
    counter = i;
    //ETD ( ml_prob, numberOfTimeSteps );
    //RK4 ( ml_prob, implicitEuler, numberOfTimeSteps );
    RKe (ml_prob, numberOfTimeSteps);
    mlSol.GetWriter()->Write (DEFAULT_OUTPUTDIR, "linear", print_vars, (i + 1) / 1);
    
    counter2++;
  }

  std::cout << " TOTAL TIME:\t" << \
            static_cast<double> (clock() - start_time) / CLOCKS_PER_SEC << std::endl;
  return 0;
}


void ETD (MultiLevelProblem& ml_prob, const unsigned & numberOfTimeSteps)
{

  const unsigned& NLayers = NumberOfLayers;

  adept::Stack& s = FemusInit::_adeptStack;

  if (assembly) {
    s.continue_recording();

  }

  else {
    s.pause_recording();
  }

  TransientLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<TransientLinearImplicitSystem> ("SWt");   // pointer to the linear implicit system named "Poisson"

  unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1u;

  Mesh* msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  SparseMatrix* KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector* RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  NumericVector* EPS = pdeSys->_EPS; // pointer to the global residual vector object in pdeSys (level)

  NumericVector* RES2;
  RES2 = NumericVector::build().release();
  RES2->init (*RES);

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solIndexh (NLayers);

  std::vector < unsigned > solIndexv (NLayers);

  std::vector < unsigned > solIndexHT (NLayers);
  std::vector < unsigned > solPdeIndexHT (NLayers);

  std::vector < unsigned > solIndexT (NLayers);

  vector< int > l2GMapRow; // local to global mapping
  vector< int > l2GMapColumn; // local to global mapping

  for (unsigned i = 0; i < NLayers; i++) {
    char name[10];
    sprintf (name, "h%d", i);
    solIndexh[i] = mlSol->GetIndex (name);   // get the position of "hi" in the sol object

    sprintf (name, "v%d", i);
    solIndexv[i] = mlSol->GetIndex (name);   // get the position of "vi" in the sol object

    sprintf (name, "HT%d", i);
    solIndexHT[i] = mlSol->GetIndex (name);   // get the position of "Ti" in the sol object
    solPdeIndexHT[i] = mlPdeSys->GetSolPdeIndex (name);   // get the position of "Ti" in the pdeSys object

    sprintf (name, "T%d", i);
    solIndexT[i] = mlSol->GetIndex (name);   // get the position of "Ti" in the sol object

  }

  unsigned solTypeh = mlSol->GetSolutionType (solIndexh[0]);   // get the finite element type for "hi"
  unsigned solTypev = mlSol->GetSolutionType (solIndexv[0]);   // get the finite element type for "vi"
  unsigned solTypeHT = mlSol->GetSolutionType (solIndexHT[0]);   // get the finite element type for "Ti"

  if (assembly) {
    KK->zero();
  }

  RES->zero();

  MatSetOption ( (static_cast<PetscMatrix*> (KK))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  for (unsigned k = 0; k < NumberOfLayers; k++) {
    for (unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++) {
//             double valueT = ( *sol->_SolOld[solIndexT[k]] ) ( i );
      double valueT = (*sol->_Sol[solIndexT[k]]) (i);
      double valueH = (*sol->_Sol[solIndexh[k]]) (i);

      double valueHT = valueT * valueH;

      sol->_Sol[solIndexHT[k]]->set (i, valueHT);
    }

    sol->_Sol[solIndexHT[k]]->close();
  }

  std::vector < double > maxW (NLayers, -1.e6);
  maxW[0] = 0.;

  unsigned start = msh->_dofOffset[solTypeHT][iproc];
  unsigned end = msh->_dofOffset[solTypeHT][iproc + 1];

  for (unsigned i =  start; i <  end; i++) {
    //EPS->zero();

    vector < double > solhm (NLayers);
    vector < double > solh (NLayers);   // local coordinates
    vector < double > solhp (NLayers);
    vector < double > solvm (NLayers);   // local coordinates
    vector < double > solvp (NLayers);   // local coordinates
    vector < adept::adouble > solHTm (NLayers);   // local coordinates
    vector < adept::adouble > solHT (NLayers);   // local coordinates
    vector < adept::adouble > solHTp (NLayers);   // local coordinates

    vector < adept::adouble > solHTmm (NLayers);   // local coordinates
    vector < adept::adouble > solHTpp (NLayers);   // local coordinates

    vector< adept::adouble > aResHT (NLayers);

    unsigned bc1 = (i == start) ? 0 : 1;
    unsigned bc2 = (i == end - 1) ? 0 : 1;

    unsigned bc3 = (i > start + 1) ? 1 : 0;
    unsigned bc4 = (i < end - 2) ? 1 : 0;

    l2GMapRow.resize (NLayers);
    l2GMapColumn.resize ( (1 + bc1 + bc2) * NLayers);

    std::fill (aResHT.begin(), aResHT.end(), 0);   //set aRes to zero

    for (unsigned j = 0; j < NLayers; j++) {

      solh[j] = (*sol->_Sol[solIndexh[j]]) (i);
      solHT[j] = (*sol->_Sol[solIndexHT[j]]) (i);
      l2GMapRow[/*NLayers +*/ j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i);

      l2GMapColumn[/*NLayers +*/ j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i);

      solvm[j] = (*sol->_Sol[solIndexv[j]]) (i);
      solvp[j] = (*sol->_Sol[solIndexv[j]]) (i + 1);

      if (i > start) {
        solhm[j] = (*sol->_Sol[solIndexh[j]]) (i - 1);
        solHTm[j] = (*sol->_Sol[solIndexHT[j]]) (i - 1);

        l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i - 1);

      }

      if (i < end - 1) {
        solhp[j] = (*sol->_Sol[solIndexh[j]]) (i + 1);
        solHTp[j] = (*sol->_Sol[solIndexHT[j]]) (i + 1);

        l2GMapColumn[ (1 + bc1) * NLayers + j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i + 1);
      }

//       if ( i > start + 1 ) {
//         solHTmm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 2 );
//         if (i == end - 1) l2GMapColumn[( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//         else l2GMapColumn[( (1 + bc1) + bc3 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//       }
//
//       if ( i < end - 2 ) {
//         solHTpp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 2 );
//         l2GMapColumn[( (1 + bc1) + bc3 + bc4 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 2 );
//       }

    }

    if (assembly) {
      s.new_recording();
    }

    vector < double > x (2);   // local coordinates

    for (unsigned j = 0; j < 2; j++) {
      unsigned xDof  = msh->GetSolutionDof (j, i, 2);   // global to global mapping between coordinates node and coordinate dof
      x[j] = (*msh->_topology->_Sol[0]) (xDof);     // global extraction and local storage for the element coordinates
    }

    double dx = x[1] - x[0];

    double b = 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );

    std::vector < double > w (NLayers + 1, 0.);
    //w[0] = 1.;

//     std::vector < double > zMid ( NLayers );
//     for ( unsigned k = 0; k < NLayers; k++ ) {
//       zMid[k] = -b + solh[k] / 2.;
//       for ( unsigned i = k + 1; i < NLayers; i++ ) {
//         zMid[k] += solh[i];
//       }
//     }

    std::vector < double > zTop (NLayers);
    zTop[0] = 0;

    for (unsigned k = 1; k < NLayers; k++) {
      zTop[k] = zTop[k - 1] - solh[k];
    }

    std::vector < double > psi2 (NLayers);

    for (unsigned k = 0; k < NLayers; k++) {
      //psi2[k] = ( - ( zMid[k] + 10 ) * zMid[k] ) / 25;
      psi2[k] = 1. - (zTop[k] + 5.) * (zTop[k] + 5.) / (25.);
    }

    double xmid = 0.5 * (x[1] + x[0]);

    for (unsigned k = NLayers; k > 1; k--) {
      w[k - 1] = - (- 4. / 625.* (xmid - 5) * (xmid - 5) * (xmid - 5)) * psi2[k - 1];
      //w[k - 1] = - ( - 40./(pow(50.,40)) * pow((xmid - 50.), 39) ) * psi2[k - 1];
      //w[k - 1] = - ( - 16./(pow(10.,16)) * pow((xmid - 10.), 15) ) * psi2[k - 1];
      //std::cout << w[k-1] << " ";

      //w[k - 1] = ( ( 10. - 2. * xmid ) / 25. ) * psi2[k - 1];
      if (maxW[k - 1] < w[k - 1]) {
        maxW[k - 1] = w[k - 1];
      }
    }

    //std::cout << std::endl;

    for (unsigned k = 0; k < NLayers; k++) {

      //BEGIN FIRST ORDER
      if (splitting) {
        if (i > start) {
          if (solvm[k] > 0) {
            aResHT[k] += solHTm[k].value() * solvm[k] / dx;

          }

          else {
            aResHT[k] += solHT[k].value() * solvm[k] / dx;
          }
        }

        if (i < end - 1) {
          if (solvp[k] > 0) {
            aResHT[k] -= solHT[k].value() * solvp[k] / dx; //first order upwind

          }

          else {
            aResHT[k] -= solHTp[k].value() * solvp[k] / dx; //first order upwind
          }
        }

      }

      else {
        if (i > start) {
          if (solvm[k] > 0) {
            aResHT[k] += solHTm[k] * solvm[k] / dx;

          }

          else {
            aResHT[k] += solHT[k] * solvm[k] / dx;
          }
        }

        if (i < end - 1) {
          if (solvp[k] > 0) {
            aResHT[k] -= solHT[k] * solvp[k] / dx; //first order upwind

          }

          else {
            aResHT[k] -= solHTp[k] * solvp[k] / dx; //first order upwind
          }
        }
      }

      //END

      //BEGIN THIRD ORDER
//       if ( i > start ) {
//         aResHT[k] += 0.5 * ( solHTm[k].value() + solHT[k].value() ) * solvm[k] / dx;
//         if ( solvm[k] > 0 ) {
//           if ( i > start + 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHT[k].value() - 2.*solHTm[k].value() + solHTmm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//       }
//       if ( i < end - 1 ) {
//         aResHT[k] -= 0.5 * ( solHTp[k].value() + solHT[k].value() ) * solvp[k] / dx;
//         if ( solvp[k] > 0 ) {
//           if (i > start) {
//             aResHT[k] -= - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvp[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 2 ) {
//             aResHT[k] -= - 1. / 6. * ( solHTpp[k].value() - 2.*solHTp[k].value() + solHT[k].value() ) * solvp[k]  / dx;
//           }
//         }
//       }
      //END

      if (k < NLayers - 1) {  //bottom
        //aResHT[k] += w[k + 1] * 0.5 * ( solHT[k] / solh[k] + solHT[k + 1] / solh[k + 1] );
        if (w[k + 1] > 0) {
          aResHT[k] += w[k + 1] * (solHT[k + 1] / solh[k + 1]);

        }

        else {
          aResHT[k] += w[k + 1] * (solHT[k] / solh[k]);
        }
      }

      if (k > 0) {  //top
        //aResHT[k] -= w[k] * 0.5 * ( solHT[k - 1] / solh[k - 1] + solHT[k] / solh[k] );
        if (w[k] > 0) {
          aResHT[k] -= w[k] * (solHT[k] / solh[k]);

        }

        else {
          aResHT[k] -= w[k] * (solHT[k - 1] / solh[k - 1]);
        }
      }

      //BEGIN MAGHEGGIONE
//       if ( k == 0 ) {
//         aResHT[k] -= w[k] * ( solHT[k] / solh[k] );
//       }
//       if ( k == NLayers - 1 ) {
//         aResHT[k] -= w[k] * ( solHT[k - 2] / solh[k - 2] );
//       }
      //END


      adept::adouble deltaZt = 0.;
      adept::adouble deltaZb = 0.;
      adept::adouble ht = 0.;
      adept::adouble hb = 0.;

      if (k > 0) {
//         ht = ( solhm[k - 1] + solhm[k] + solhp[k - 1] + solhp[k] ) / 4.;
        ht = (solh[k - 1] + solh[k]) / 2.;
        deltaZt = (solHT[k - 1] / solh[k - 1] - solHT[k] / solh[k]) / ht;

      }

      else {
//         ht = 0.5 * ( solhm[k] + solhp[k] );
        ht = solh[k];
        deltaZt = 0.* (0. - solHT[k]) / ht;
      }

      if (k < NLayers - 1) {
//         hb = ( solhm[k] + solhm[k + 1] + solhp[k] + solhp[k + 1] ) / 4.;
        hb = (solh[k] + solh[k + 1]) / 2.;
        deltaZb = (solHT[k] / solh[k] - solHT[k + 1] / solh[k + 1]) / hb;

      }

      else {
//         hb = 0.5 * ( solhm[k] + solhp[k] );
        hb = solh[k];
        deltaZb = 0.* (solHT[k] - 0.) / hb;
      }

      //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<deltaZt - deltaZb<<std::endl;

      aResHT[k] += solh[k] * k_v * (deltaZt - deltaZb) / ( (ht + hb) / 2.);      // vertical diffusion

      //aResHT[k] += ((solhp[k] - solhm[k]) * k_h * (solHTp[k] - solHTm[k])) / (dx*dx); // horizontal diffusion
//       if (splitting){
//         if(i > start){
//           aResHT[k] += k_h * (0.5 * (solhm[k] + solh[k])) * (solHTm[k].value() / solhm[k] - solHT[k].value() / solh[k])/(dx*dx); // horizontal diffusion
//         }
//         if(i < end-1){
//           aResHT[k] += k_h * (0.5 * (solhp[k] + solh[k])) * (solHTp[k].value() / solhp[k] - solHT[k].value() / solh[k])/(dx*dx); // horizontal diffusion
//         }
//       }
//       else{
//         if(i > start){
//           aResHT[k] += k_h * (0.5 * (solhm[k] + solh[k])) * (solHTm[k] / solhm[k] - solHT[k] / solh[k])/(dx*dx); // horizontal diffusion
//         }
//         if(i < end-1){
//           aResHT[k] += k_h * (0.5 * (solhp[k] + solh[k])) * (solHTp[k] / solhp[k] - solHT[k] / solh[k])/(dx*dx); // horizontal diffusion
//         }
//       }

    }

    vector< double > Res (NLayers);   // local redidual vector
    vector< double > solht (NLayers);   // local redidual vector

    for (unsigned k = 0; k < NLayers; k++) {
      Res[k] =  aResHT[k].value();
      solht[k] = solHT[k].value();
    }

    RES->add_vector_blocked (Res, l2GMapRow);

    if (assembly) {
      s.dependent (&aResHT[0], NLayers);

      // define the independent variables
      s.independent (&solHT[0], NLayers);

      if (!splitting) {
        if (i > start) {
          s.independent (&solHTm[0], NLayers);
        }

        if (i < end - 1) {
          s.independent (&solHTp[0], NLayers);
        }

        /*    if ( i > start + 1) {
        s.independent ( &solHTmm[0], NLayers );
        }
        if ( i < end - 2 ) {
          s.independent ( &solHTpp[0], NLayers );
        } */
        // get the jacobian matrix (ordered by row major )
        vector < double > Jac (NLayers * NLayers * (1 + bc1 + bc2));
        s.jacobian (&Jac[0], true);

        //store K in the global matrix KK
        KK->add_matrix_blocked (Jac, l2GMapRow, l2GMapColumn);

      }

      else {
        vector < double > Jac (NLayers * NLayers);
        s.jacobian (&Jac[0], true);

        //store K in the global matrix KK
        KK->add_matrix_blocked (Jac, l2GMapRow, l2GMapRow);
      }

//         MatCreate(MPI_COMM_SELF,&A);
//         MatSetSizes(A,NLayers,NLayers,NLayers,NLayers);
//         MatSetType(A,MATSEQDENSE);
//         MatSeqDenseSetPreallocation(A,Jac);
//         MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
//         MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

      s.clear_independents();
      s.clear_dependents();
    }
  }

  RES->close();

  if (assembly) {
    KK->close();
  }

  // printing of the max value of w in every layer
//   for ( unsigned k = 0; k < NLayers; k++ ) {
//     std::cout << "layer " << k << " " << maxW[k] << std::endl;
//   }

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//   PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//   PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(),viewer);
//   double a;
//   std::cin>>a;

//  abort();
  //std::cout << "dt = " << dt << std::endl;
  //std::cout << "first stage " << std::endl;

  MFN mfn;
  Mat A = (static_cast<PetscMatrix*> (KK))->mat();
  FN f;

  Vec v = (static_cast< PetscVector* > (RES))->vec();
  Vec y = (static_cast< PetscVector* > (EPS))->vec();

  //BEGIN condition number of A
//     SVD svd;
//     PetscReal sigma_1,sigma_n;
//     PetscInt nconv1,nconv2;
//
//     SVDCreate(PETSC_COMM_WORLD,&svd);
//     SVDSetOperator(svd,A);
//     SVDSetFromOptions(svd);
//     SVDSetDimensions(svd,1,PETSC_DEFAULT,PETSC_DEFAULT);
//     //First request a singular value from one end of the spectrum
//     SVDSetWhichSingularTriplets(svd,SVD_LARGEST);
//     SVDSolve(svd);
//     //Get number of converged singular values
//     SVDGetConverged(svd,&nconv1);
//     //Get converged singular values: largest singular value is stored in sigma_1.
//     if (nconv1 > 0) {
//       SVDGetSingularTriplet(svd,0,&sigma_1,NULL,NULL);
//     }
//     else {
//       PetscPrintf(PETSC_COMM_WORLD," Unable to compute large singular value!\n\n");
//     }
//     //Request a singular value from the other end of the spectrum
//     SVDSetWhichSingularTriplets(svd,SVD_SMALLEST);
//     SVDSolve(svd);
//     //Get number of converged singular triplets
//     SVDGetConverged(svd,&nconv2);
//     //Get converged singular values: smallest singular value is stored in sigma_n.
//     if (nconv2 > 0) {
//       SVDGetSingularTriplet(svd,0,&sigma_n,NULL,NULL);
//     }
//     else {
//       PetscPrintf(PETSC_COMM_WORLD," Unable to compute small singular value!\n\n");
//     }
//     //Display solution and clean up
//     if (nconv1 > 0 && nconv2 > 0) {
//       PetscPrintf(PETSC_COMM_WORLD," Computed singular values: sigma_1=%.4f, sigma_n=%.4f\n",(double)sigma_1,(double)sigma_n);
//       PetscPrintf(PETSC_COMM_WORLD," Estimated condition number: sigma_1/sigma_n=%.4f\n\n",(double)(sigma_1/sigma_n));
//     }
//     //Free work space
//     SVDDestroy(&svd);
  //END

  PetscReal tol;
  PetscInt ncv, maxit, its;

  MFNCreate (PETSC_COMM_WORLD, &mfn);

  MFNSetOperator (mfn, A);
  MFNGetFN (mfn, &f);

  FNPhiSetIndex (f, 1);
  FNSetType (f, FNPHI);
// FNView(f,PETSC_VIEWER_STDOUT_WORLD);

  FNSetScale (f, dt, dt);
  MFNSetDimensions (mfn, 15);
  tol = 1e-6;
  MFNSetTolerances (mfn, tol, PETSC_DEFAULT);
  MFNSetFromOptions (mfn);

  clock_t etd_time;

  if (counter2 == 0) {
    etd_time = clock();
//         PetscInt vSize;
//         VecGetSize(v, &vSize);
//         std::cout<< "--------------------- v --------------------- " <<std::endl;
//         for ( PetscInt kk=0; kk<vSize; kk++ ) {
//             PetscScalar valueHT = 0.;
//             VecGetValues ( v, 1, &kk, &valueHT );
//             std::cout<< valueHT << std::endl;
//         }
  }

  MFNSolve (mfn, v, y);

  if (counter2 == 0) {
    std::cout << " ETD TIME :\t" << static_cast<double> (clock() - etd_time) / CLOCKS_PER_SEC << std::endl;
//         PetscInt ySize;
//         VecGetSize(y, &ySize);
//         std::cout<< "--------------------- y --------------------- " <<std::endl;
//         for ( PetscInt kk=0; kk<ySize; kk++ ) {
//             PetscScalar valueHT = 0.;
//             VecGetValues ( y, 1, &kk, &valueHT );
//             std::cout<< valueHT << std::endl;
//         }
//         std::cout<< "--------------------- END --------------------- " <<std::endl;
  }

  //BEGIN Get some information from the solver and display it
//     MFNGetIterationNumber(mfn,&its);
//     PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);
//     MFNGetDimensions(mfn,&ncv);
//     PetscPrintf(PETSC_COMM_WORLD," Subspace dimension: %D\n",ncv);
//     MFNGetTolerances(mfn,&tol,&maxit);
//     PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);
  //END

  MFNDestroy (&mfn);

  //VecView(y,PETSC_VIEWER_STDOUT_WORLD);

  //std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";

  sol->UpdateSol (mlPdeSys->GetSolPdeIndex(), EPS, pdeSys->KKoffset);

  if (twostage == true) {

//         std::cout << "second stage " << std::endl;

    RES2->zero();

    for (unsigned i =  start; i <  end; i++) {

      vector < double > solhm (NLayers);
      vector < double > solh (NLayers);   // local coordinates
      vector < double > solhp (NLayers);
      vector < double > solvm (NLayers);   // local coordinates
      vector < double > solvp (NLayers);   // local coordinates
      vector < double > solHTm (NLayers);   // local coordinates
      vector < double > solHT (NLayers);   // local coordinates
      vector < double > solHTp (NLayers);   // local coordinates

      vector < double > solHTmm (NLayers);   // local coordinates
      vector < double > solHTpp (NLayers);   // local coordinates

      vector< double > aResHT (NLayers, 0.);

      unsigned bc1 = (i == start) ? 0 : 1;
      unsigned bc2 = (i == end - 1) ? 0 : 1;

      unsigned bc3 = (i > start + 1) ? 1 : 0;
      unsigned bc4 = (i < end - 2) ? 1 : 0;

      l2GMapRow.resize (NLayers);
      l2GMapColumn.resize ( (1 + bc1 + bc2) * NLayers);

      //std::fill ( aResHT.begin(), aResHT.end(), 0 ); //set aRes to zero

      for (unsigned j = 0; j < NLayers; j++) {

        solh[j] = (*sol->_Sol[solIndexh[j]]) (i);
        solHT[j] = (*sol->_Sol[solIndexHT[j]]) (i);
        l2GMapRow[/*NLayers +*/ j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i);

        l2GMapColumn[/*NLayers +*/ j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i);

        solvm[j] = (*sol->_Sol[solIndexv[j]]) (i);
        solvp[j] = (*sol->_Sol[solIndexv[j]]) (i + 1);


        if (i > start) {
          solhm[j] = (*sol->_Sol[solIndexh[j]]) (i - 1);
          solHTm[j] = (*sol->_Sol[solIndexHT[j]]) (i - 1);

          l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i - 1);

        }

        if (i < end - 1) {
          solhp[j] = (*sol->_Sol[solIndexh[j]]) (i + 1);
          solHTp[j] = (*sol->_Sol[solIndexHT[j]]) (i + 1);

          l2GMapColumn[ (1 + bc1) * NLayers + j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i + 1);
        }

//       if ( i > start + 1 ) {
//         solHTmm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 2 );
//         if (i == end - 1) l2GMapColumn[( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//         else l2GMapColumn[( (1 + bc1) + bc3 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//       }
//
//       if ( i < end - 2 ) {
//         solHTpp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 2 );
//         l2GMapColumn[( (1 + bc1) + bc3 + bc4 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 2 );
//       }

      }

      // s.new_recording();

      vector < double > x (2);   // local coordinates

      for (unsigned j = 0; j < 2; j++) {
        unsigned xDof  = msh->GetSolutionDof (j, i, 2);   // global to global mapping between coordinates node and coordinate dof
        x[j] = (*msh->_topology->_Sol[0]) (xDof);     // global extraction and local storage for the element coordinates
      }

      double dx = x[1] - x[0];

      double b = 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );

      std::vector < double > w (NLayers + 1, 0.);
      //w[0] = 1.;

      std::vector < double > zTop (NLayers);
      zTop[0] = 0;

      for (unsigned k = 1; k < NLayers; k++) {
        zTop[k] = zTop[k - 1] - solh[k];
      }

      std::vector < double > psi2 (NLayers);

      for (unsigned k = 0; k < NLayers; k++) {
        //psi2[k] = ( - ( zMid[k] + 10 ) * zMid[k] ) / 25;
        psi2[k] = 1. - (zTop[k] + 5.) * (zTop[k] + 5.) / (25.);
      }

      double xmid = 0.5 * (x[1] + x[0]);

      for (unsigned k = NLayers; k > 1; k--) {
        w[k - 1] = - (-4. / 625.* (xmid - 5.) * (xmid - 5.) * (xmid - 5.)) * psi2[k - 1];

        //w[k - 1] = - ( - 40./(pow(50.,40)) * pow((xmid - 50.), 39) ) * psi2[k - 1];
        //w[k - 1] = - ( - 16./(pow(10.,16)) * pow((xmid - 10.), 15) ) * psi2[k - 1];
        //w[k - 1] = ( ( 10. - 2. * xmid ) / 25. ) * psi2[k - 1];
        if (maxW[k - 1] < w[k - 1]) {
          maxW[k - 1] = w[k - 1];
        }
      }


      for (unsigned k = 0; k < NLayers; k++) {

        //BEGIN FIRST ORDER
        if (i > start) {
          //aResHT[k] += 0.5 * (solHTm[k] + solHT[k]) * solvm[k]  / dx; //second order
          if (solvm[k] > 0) {
            aResHT[k] += solHTm[k] * solvm[k] / dx;

          }

          else {
            aResHT[k] += solHT[k] * solvm[k] / dx;
          }
        }

        if (i < end - 1) {
          //aResHT[k] -= 0.5 * (solHT[k] + solHTp[k]) * solvp[k]  / dx; //second order
          if (solvp[k] > 0) {
            aResHT[k] -= solHT[k] * solvp[k] / dx; //first order upwind

          }

          else {
            aResHT[k] -= solHTp[k] * solvp[k] / dx; //first order upwind
          }
        }

//       else{
//         aResHT[k] -= solHT[k] /*.value()*/  / dx; //first order upwind
//       }
        //END

        //BEGIN THIRD ORDER
//       if ( i > start ) {
//         aResHT[k] += 0.5 * ( solHTm[k].value() + solHT[k].value() ) * solvm[k] / dx;
//         if ( solvm[k] > 0 ) {
//           if ( i > start + 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHT[k].value() - 2.*solHTm[k].value() + solHTmm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//       }
//       if ( i < end - 1 ) {
//         aResHT[k] -= 0.5 * ( solHTp[k].value() + solHT[k].value() ) * solvp[k] / dx;
//         if ( solvp[k] > 0 ) {
//           if (i > start) {
//             aResHT[k] -= - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvp[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 2 ) {
//             aResHT[k] -= - 1. / 6. * ( solHTpp[k].value() - 2.*solHTp[k].value() + solHT[k].value() ) * solvp[k]  / dx;
//           }
//         }
//       }
        //END

        if (k < NLayers - 1) {  //bottom
          //aResHT[k] += w[k + 1] * 0.5 * ( solHT[k] / solh[k] + solHT[k + 1] / solh[k + 1] );
          if (w[k + 1] > 0) {
            aResHT[k] += w[k + 1] * (solHT[k + 1] / solh[k + 1]);

          }

          else {
            aResHT[k] += w[k + 1] * (solHT[k] / solh[k]);
          }
        }

        if (k > 0) {  //top
          //aResHT[k] -= w[k] * 0.5 * ( solHT[k - 1] / solh[k - 1] + solHT[k] / solh[k] );
          if (w[k] > 0) {
            aResHT[k] -= w[k] * (solHT[k] / solh[k]);

          }

          else {
            aResHT[k] -= w[k] * (solHT[k - 1] / solh[k - 1]);
          }
        }

        //BEGIN MAGHEGGIONE
//         if ( k == 0 ) {
//           aResHT[k] -= w[k] * ( solHT[k] / solh[k] );
//         }
//         if ( k == NLayers - 1 ) {
//           aResHT[k] -= w[k] * ( solHT[k - 2] / solh[k - 2] );
//         }
        //END

        double deltaZt = 0.;
        double deltaZb = 0.;
        double ht = 0.;
        double hb = 0.;

        if (k > 0) {
//        ht = ( solhm[k - 1] + solhm[k] + solhp[k - 1] + solhp[k] ) / 4.;
          ht = (solh[k - 1] + solh[k]) / 2.;
          deltaZt = (solHT[k - 1] / solh[k - 1] - solHT[k] / solh[k]) / ht;

        }

        else {
//        ht = 0.5 * ( solhm[k] + solhp[k] );
          ht = solh[k];
          deltaZt = 0.* (0. - solHT[k]) / ht;
        }

        if (k < NLayers - 1) {
//        hb = ( solhm[k] + solhm[k + 1] + solhp[k] + solhp[k + 1] ) / 4.;
          hb = (solh[k] + solh[k + 1]) / 2.;
          deltaZb = (solHT[k] / solh[k] - solHT[k + 1] / solh[k + 1]) / hb;

        }

        else {
//        hb = 0.5 * ( solhm[k] + solhp[k] );
          hb = solh[k];
          deltaZb = 0.* (solHT[k] - 0.) / hb;
        }

        //std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<deltaZt - deltaZb<<std::endl;

        aResHT[k] += solh[k] * k_v * (deltaZt - deltaZb) / ( (ht + hb) / 2.);      // vertical diffusion

//       //aResHT[k] += ((solhp[k] - solhm[k]) * k_h * (solHTp[k] - solHTm[k])) / (dx*dx); // horizontal diffusion
//         if(i > start){
//           aResHT[k] += k_h * (0.5 * (solhm[k] + solh[k])) * (solHTm[k] / solhm[k] - solHT[k] / solh[k])/(dx*dx); // horizontal diffusion
//         }
//         if(i < end-1){
//           aResHT[k] += k_h * (0.5 * (solhp[k] + solh[k])) * (solHTp[k] / solhp[k] - solHT[k] / solh[k])/(dx*dx); // horizontal diffusion
//         }

      }

      RES2->add_vector_blocked (aResHT, l2GMapRow);

    }

    RES2->close();


    //BEGIN ASSEMBLY: R2 = RES2 - RES - KK * EPS = RESnew - RESold - KK * (Vnew - Vold) = (ResNew - KK * Vnew) - (ResOld - KK * Vold) = 0 - 0 ;
    RES2->scale (-1.);
    RES2->add (*RES);
    RES2->add_vector (*EPS, *KK);
    RES2->scale (-1.);
    //END ASSEMBLY R2

    EPS->zero();

    //SLEPC
    MFN mfn;
    Mat A2 = (static_cast<PetscMatrix*> (KK))->mat();
    FN f2;

    //std::cout << "dt = " << dt << std::endl;

    Vec v2 = (static_cast< PetscVector* > (RES2))->vec();

    //VecView(v2,PETSC_VIEWER_STDOUT_WORLD);

    Vec y2 = (static_cast< PetscVector* > (EPS))->vec();

    MFNCreate (PETSC_COMM_WORLD, &mfn);

    MFNSetOperator (mfn, A2);
    MFNGetFN (mfn, &f2);

    FNPhiSetIndex (f2, 2);
    FNSetType (f2, FNPHI);
    //FNView(f,PETSC_VIEWER_STDOUT_WORLD);

    FNSetScale (f2, dt, dt);
    MFNSetDimensions (mfn, 15);
    MFNSetTolerances (mfn, 1e-6, PETSC_DEFAULT);
    MFNSetFromOptions (mfn);

    MFNSolve (mfn, v2, y2);

    MFNDestroy (&mfn);

    //VecView(y2,PETSC_VIEWER_STDOUT_WORLD);

    sol->UpdateSol (mlPdeSys->GetSolPdeIndex(), EPS, pdeSys->KKoffset);

  }

  //PARAVIEW
  unsigned solIndexeta = mlSol->GetIndex ("eta");
  unsigned solIndexb = mlSol->GetIndex ("b");
  sol->_Sol[solIndexeta]->zero();

  for (unsigned k = 0; k < NumberOfLayers; k++) {
    sol->_Sol[solIndexeta]->add (*sol->_Sol[solIndexh[k]]);
  }

  sol->_Sol[solIndexeta]->add (-1, *sol->_Sol[solIndexb]);

  for (unsigned k = 0; k < NumberOfLayers; k++) {

    //if(counter == numberOfTimeSteps-2) std::cout << "T" << k << "  ----------------------------------------------------" << std::endl;
    for (unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++) {
      double valueHT = (*sol->_Sol[solIndexHT[k]]) (i);
      double valueH = (*sol->_Sol[solIndexh[k]]) (i);

      double valueT = valueHT / valueH;

      if (counter == numberOfTimeSteps - 1) {
        std::cout.precision (14);
        std::cout << valueT << std::endl;
      }

      //if ( i == 10 ) std::cout << "temperature " << valueT << std::endl;
      //if (i == 0) valueT = 15.;
      //if (i == msh->_dofOffset[solTypeHT][iproc + 1] - 1 ) valueT = 15.;
      sol->_Sol[solIndexT[k]]->set (i, valueT);
    }

    sol->_Sol[solIndexT[k]]->close();

  }

  delete RES2;
}


void RK4 (MultiLevelProblem& ml_prob, const bool & implicitEuler, const unsigned & numberOfTimeSteps)
{

  const unsigned& NLayers = NumberOfLayers;

  adept::Stack& s = FemusInit::_adeptStack;

  TransientLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<TransientLinearImplicitSystem> ("SWt");   // pointer to the linear implicit system named "Poisson"

  unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1u;

  Mesh* msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solIndexh (NLayers);

  std::vector < unsigned > solIndexv (NLayers);

  std::vector < unsigned > solIndexHT (NLayers);
  std::vector < unsigned > solPdeIndexHT (NLayers);

  std::vector < unsigned > solIndexT (NLayers);

  vector< int > l2GMapRow; // local to global mapping
  vector< int > l2GMapColumn; // local to global mapping

  for (unsigned i = 0; i < NLayers; i++) {
    char name[10];
    sprintf (name, "h%d", i);
    solIndexh[i] = mlSol->GetIndex (name);   // get the position of "hi" in the sol object

    sprintf (name, "v%d", i);
    solIndexv[i] = mlSol->GetIndex (name);   // get the position of "vi" in the sol object

    sprintf (name, "HT%d", i);
    solIndexHT[i] = mlSol->GetIndex (name);   // get the position of "Ti" in the sol object
    solPdeIndexHT[i] = mlPdeSys->GetSolPdeIndex (name);   // get the position of "Ti" in the pdeSys object

    sprintf (name, "T%d", i);
    solIndexT[i] = mlSol->GetIndex (name);   // get the position of "Ti" in the sol object

  }

  unsigned solTypeh = mlSol->GetSolutionType (solIndexh[0]);   // get the finite element type for "hi"
  unsigned solTypev = mlSol->GetSolutionType (solIndexv[0]);   // get the finite element type for "vi"
  unsigned solTypeHT = mlSol->GetSolutionType (solIndexHT[0]);   // get the finite element type for "Ti"

  for (unsigned k = 0; k < NumberOfLayers; k++) {
    for (unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++) {
      double valueT = (*sol->_SolOld[solIndexT[k]]) (i);
      //double valueT = ( *sol->_Sol[solIndexT[k]] ) ( i );
      double valueH = (*sol->_Sol[solIndexh[k]]) (i);

      double valueHT = valueT * valueH;

      sol->_Sol[solIndexHT[k]]->set (i, valueHT);
    }

    sol->_Sol[solIndexHT[k]]->close();
  }

  std::vector < double > maxW (NLayers, -1.e6);
  maxW[0] = 0.;

  unsigned start = msh->_dofOffset[solTypeHT][iproc];
  unsigned end = msh->_dofOffset[solTypeHT][iproc + 1];

  for (unsigned i =  start; i <  end; i++) {

    vector < double > solhm (NLayers);
    vector < double > solh (NLayers);   // local coordinates
    vector < double > solhp (NLayers);
    vector < double > solvm (NLayers);   // local coordinates
    vector < double > solvp (NLayers);   // local coordinates
    vector < double > solHTm (NLayers);   // local coordinates
    vector < double > solHT (NLayers);   // local coordinates
    vector < double > solHTp (NLayers);   // local coordinates

    vector < double > solHTmm (NLayers);   // local coordinates
    vector < double > solHTpp (NLayers);   // local coordinates

    //vector< adept::adouble > aResHT ( NLayers );

    unsigned bc1 = (i == start) ? 0 : 1;
    unsigned bc2 = (i == end - 1) ? 0 : 1;

//     unsigned bc3 = ( i > start + 1 ) ? 1 : 0;
//     unsigned bc4 = ( i < end - 2 ) ? 1 : 0;

    l2GMapRow.resize (NLayers);
    l2GMapColumn.resize ( (1 + bc1 + bc2) * NLayers);

    //std::fill ( aResHT.begin(), aResHT.end(), 0 ); //set aRes to zero

    for (unsigned j = 0; j < NLayers; j++) {

      solh[j] = (*sol->_Sol[solIndexh[j]]) (i);
      solHT[j] = (*sol->_Sol[solIndexHT[j]]) (i);
      l2GMapRow[/*NLayers +*/ j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i);

      l2GMapColumn[/*NLayers +*/ j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i);

      solvm[j] = (*sol->_Sol[solIndexv[j]]) (i);
      solvp[j] = (*sol->_Sol[solIndexv[j]]) (i + 1);

      if (i > start) {
        solhm[j] = (*sol->_Sol[solIndexh[j]]) (i - 1);
        solHTm[j] = (*sol->_Sol[solIndexHT[j]]) (i - 1);

        l2GMapColumn[NLayers + j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i - 1);

      }

      if (i < end - 1) {
        solhp[j] = (*sol->_Sol[solIndexh[j]]) (i + 1);
        solHTp[j] = (*sol->_Sol[solIndexHT[j]]) (i + 1);

        l2GMapColumn[ (1 + bc1) * NLayers + j] = pdeSys->GetSystemDof (solIndexHT[j], solPdeIndexHT[j], 0, i + 1);
      }

//       if ( i > start + 1 ) {
//         solHTmm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 2 );
//         if (i == end - 1) l2GMapColumn[( 1 + bc1 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//         else l2GMapColumn[( (1 + bc1) + bc3 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i - 2 );
//       }
//
//       if ( i < end - 2 ) {
//         solHTpp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 2 );
//         l2GMapColumn[( (1 + bc1) + bc3 + bc4 ) * NLayers + j] = pdeSys->GetSystemDof ( solIndexHT[j], solPdeIndexHT[j], 0, i + 2 );
//       }

    }

//   s.new_recording();

    vector < double > x (2);   // local coordinates

    for (unsigned j = 0; j < 2; j++) {
      unsigned xDof  = msh->GetSolutionDof (j, i, 2);   // global to global mapping between coordinates node and coordinate dof
      x[j] = (*msh->_topology->_Sol[0]) (xDof);     // global extraction and local storage for the element coordinates
    }

    double dx = x[1] - x[0];

    double b = 10.; //( H_shelf + H_0 / 2 * (1 + tanh(hh / phi)) );

    std::vector < double > w (NLayers + 1, 0.);
    //w[0] = 1.;

//     std::vector < double > zMid ( NLayers );
//     for ( unsigned k = 0; k < NLayers; k++ ) {
//       zMid[k] = -b + solh[k] / 2.;
//       for ( unsigned i = k + 1; i < NLayers; i++ ) {
//         zMid[k] += solh[i];
//       }
//     }

    std::vector < double > zTop (NLayers);
    zTop[0] = 0;

    for (unsigned k = 1; k < NLayers; k++) {
      zTop[k] = zTop[k - 1] - solh[k];
    }

    std::vector < double > psi2 (NLayers);

    for (unsigned k = 0; k < NLayers; k++) {
      //psi2[k] = ( - ( zMid[k] + 10 ) * zMid[k] ) / 25;
      psi2[k] = 1. - (zTop[k] + 5.) * (zTop[k] + 5.) / (25.);
    }

    double xmid = 0.5 * (x[1] + x[0]);

    for (unsigned k = NLayers; k > 1; k--) {
      w[k - 1] = - (-4. / 625.* (xmid - 5) * (xmid - 5) * (xmid - 5)) * psi2[k - 1];

      //w[k - 1] = - ( - 40./(pow(50.,40)) * pow((xmid - 50.), 39) ) * psi2[k - 1];
      //w[k - 1] = - ( - 16./(pow(10.,16)) * pow((xmid - 10.), 15) ) * psi2[k - 1];
      //w[k - 1] = ( ( 10. - 2. * xmid ) / 25. ) * psi2[k - 1];
      if (maxW[k - 1] < w[k - 1]) {
        maxW[k - 1] = w[k - 1];
      }
    }


    std::vector < double > k1_RK (NLayers, 0.);
    std::vector < double > k2_RK (NLayers, 0.);
    std::vector < double > k3_RK (NLayers, 0.);
    std::vector < double > k4_RK (NLayers, 0.);

    for (unsigned RK_step = 0; RK_step < 4; RK_step++) {
      for (unsigned k = 0; k < NLayers; k++) {
        double RHS = 0.;
        double addition = 0.;

        if (RK_step == 1) {
          addition = k1_RK[k] * 0.5;

        }

        else if (RK_step == 2) {
          addition = k2_RK[k] * 0.5;

        }

        else if (RK_step == 3) {
          addition = k3_RK[k];
        }

        //BEGIN FIRST ORDER
        if (i > start) {
          if (solvm[k] > 0) {
            RHS += (solHTm[k] + addition) * solvm[k] / dx;

          }

          else {
            RHS += (solHT[k] + addition) * solvm[k] / dx;
          }
        }

        if (i < end - 1) {
          if (solvp[k] > 0) {
            RHS -= (solHT[k] + addition) * solvp[k] / dx;   //first order upwind

          }

          else {
            RHS -= (solHTp[k] + addition) * solvp[k] / dx;   //first order upwind
          }
        }

        //END

        //BEGIN THIRD ORDER
//       if ( i > start ) {
//         aResHT[k] += 0.5 * ( solHTm[k].value() + solHT[k].value() ) * solvm[k] / dx;
//         if ( solvm[k] > 0 ) {
//           if ( i > start + 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHT[k].value() - 2.*solHTm[k].value() + solHTmm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 1 ) {
//             aResHT[k] += - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvm[k]  / dx;
//           }
//         }
//       }
//       if ( i < end - 1 ) {
//         aResHT[k] -= 0.5 * ( solHTp[k].value() + solHT[k].value() ) * solvp[k] / dx;
//         if ( solvp[k] > 0 ) {
//           if (i > start) {
//             aResHT[k] -= - 1. / 6. * ( solHTp[k].value() - 2.*solHT[k].value() + solHTm[k].value() ) * solvp[k]  / dx;
//           }
//         }
//         else {
//           if ( i < end - 2 ) {
//             aResHT[k] -= - 1. / 6. * ( solHTpp[k].value() - 2.*solHTp[k].value() + solHT[k].value() ) * solvp[k]  / dx;
//           }
//         }
//       }
        //END

        if (k < NLayers - 1) {
//           RHS += w[k + 1] * 0.5 * ( ( solHT[k] + addition ) / solh[k] + ( solHT[k + 1] + addition ) / solh[k + 1] );
          if (w[k + 1] > 0) {
            RHS += w[k + 1] * ( (solHT[k + 1] / solh[k + 1] + addition));

          }

          else {
            RHS += w[k + 1] * ( (solHT[k] / solh[k] + addition));
          }
        }

        if (k > 0) {
//           RHS -= w[k] * 0.5 * ( ( solHT[k - 1] + addition ) / solh[k - 1] + ( solHT[k] + addition ) / solh[k] );
          if (w[k] > 0) {
            RHS -= w[k] * ( (solHT[k] / solh[k] + addition));

          }

          else {
            RHS -= w[k] * ( (solHT[k - 1] / solh[k - 1] + addition));
          }
        }

//         if(i > start){
//           RHS += k_h * (0.5 * (solhm[k] + solh[k])) * ((solHTm[k] / solhm[k] + addition)  - (solHT[k] / solh[k] + addition))/(dx*dx); // horizontal diffusion
//         }
//         if(i < end-1){
//           RHS += k_h * (0.5 * (solhp[k] + solh[k])) * ((solHTp[k] / solhp[k] + addition)  - (solHT[k] / solh[k] + addition))/(dx*dx); // horizontal diffusion
//         }


        if (RK_step == 0) {
          k1_RK[k] = RHS * dt;

        }

        else if (RK_step == 1) {
          k2_RK[k] = RHS * dt;

        }

        else if (RK_step == 2) {
          k3_RK[k] = RHS * dt;

        }

        else {
          k4_RK[k] = RHS * dt;
        }

      }

    }

    for (unsigned k = 0; k < NumberOfLayers; k++) {

//       std::cout << k1_RK[k] << " " << k2_RK[k] << " " << k3_RK[k] << " " << k4_RK[k] << std::endl;

      double valueHT = solHT[k] + 1. / 6. * (k1_RK[k] + 2.*k2_RK[k] + 2.*k3_RK[k] + k4_RK[k]);
      sol->_Sol[solIndexHT[k]]->set (i, valueHT);
      sol->_Sol[solIndexHT[k]]->close();
    }


    if (implicitEuler == false) {
      //BEGIN forward Euler for vertical diffusion
      std::vector < double > vert_diff (NLayers, 0.);

      for (unsigned k = 0; k < NumberOfLayers; k++) {
        double deltaZt = 0.;
        double deltaZb = 0.;
        double ht = 0.;
        double hb = 0.;

        if (k > 0) {
//        ht = ( solhm[k - 1] + solhm[k] + solhp[k - 1] + solhp[k] ) / 4.;
          ht = (solh[k - 1] + solh[k]) / 2.;
          deltaZt = (solHT[k - 1] / solh[k - 1] - solHT[k] / solh[k]) / ht;

        }

        else {
//        ht = 0.5 * ( solhm[k] + solhp[k] );
          ht = solh[k];
          deltaZt = 0.* (0. - solHT[k]) / ht;
        }

        if (k < NLayers - 1) {
//        hb = ( solhm[k] + solhm[k + 1] + solhp[k] + solhp[k + 1] ) / 4.;
          hb = (solh[k] + solh[k + 1]) / 2.;
          deltaZb = (solHT[k] / solh[k] - solHT[k + 1] / solh[k + 1]) / hb;

        }

        else {
//        hb = 0.5 * ( solhm[k] + solhp[k] );
          hb = solh[k];
          deltaZb = 0.* (solHT[k] - 0.) / hb;
        }

        //std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAA" << deltaZt - deltaZb << std::endl;
        vert_diff[k] = solh[k] * k_v * (deltaZt - deltaZb) / ( (ht + hb) / 2.);      // vertical diffusion
      }

      for (unsigned k = 0; k < NumberOfLayers; k++) {
        double valueHT = (*sol->_Sol[solIndexHT[k]]) (i);
        double valueH = (*sol->_Sol[solIndexh[k]]) (i);

        double valueT = valueHT / valueH;
        valueT = valueT + dt * vert_diff[k];

        sol->_Sol[solIndexT[k]]->set (i, valueT);
        sol->_Sol[solIndexT[k]]->close();
      }

      //END
    }

    else if (implicitEuler == true) {

      std::vector <double> Trhs (NLayers, 0.);
      std::vector < std::vector < double > > sysMatrix (NLayers);

      for (unsigned k = 0; k < NLayers; k++) {
        sysMatrix[k].assign (NLayers, 0.);

        double A = 0.;
        double C = 0.;
        double ht = 0.;
        double hb = 0.;

        if (k > 0) {
          //ht = ( solhm[k - 1] + solhm[k] + solhp[k - 1] + solhp[k] ) / 4.;
          ht = (solh[k - 1] + solh[k]) / 2.;
          A = solh[k] * k_v / ht ;

        }

        else {
          //ht = 0.5 * ( solhm[k] + solhp[k] );
          ht = solh[k];
        }

        if (k < NLayers - 1) {
          //hb = ( solhm[k] + solhm[k + 1] + solhp[k] + solhp[k + 1] ) / 4.;
          hb = (solh[k] + solh[k + 1]) / 2.;
          C = solh[k] * k_v / hb;
          C /= (ht + hb) * 0.5 ;

          if (k > 0) {
            A /= (ht + hb) * 0.5 ;
          }

        }

        else {
          //hb = 0.5 * ( solhm[k] + solhp[k] );
          hb = solh[k];
          A /= (ht + hb) * 0.5 ;
        }

        double B = 1. - A - C;

        sysMatrix[k][k] = B;

        if (k > 0) {
          sysMatrix[k][k - 1] = A;
        }

        if (k < NLayers - 1) {
          sysMatrix[k][k + 1] = C;
        }

        Trhs[k] = (*sol->_Sol[solIndexHT[k]]) (i) / (*sol->_Sol[solIndexh[k]]) (i);

      }

      //risolvere il sistema Nlayer X Nlayer
      KSP                solver;
      Mat                triDiagA;
      Vec                b, x;
      PetscInt           k, j, nlayers;
      PetscErrorCode     ierr;

      nlayers = static_cast<PetscInt> (NLayers);
      ierr = VecCreate (PETSC_COMM_WORLD, &x);
      ierr = VecSetSizes (x, PETSC_DECIDE, nlayers);
      ierr = VecSetFromOptions (x);
      ierr = VecDuplicate (x, &b);
      ierr = MatCreate (PETSC_COMM_WORLD, &triDiagA);
      ierr = MatSetSizes (triDiagA, PETSC_DECIDE, PETSC_DECIDE, nlayers, nlayers);
      ierr = MatSetFromOptions (triDiagA);
      ierr = MatSetUp (triDiagA);

      for (k = 0; k < nlayers; k++) {
        MatSetValues (triDiagA, 1, &k, 1, &k, &sysMatrix[k][k], INSERT_VALUES);
        VecSetValues (b, 1, &k, &Trhs[k], INSERT_VALUES);

        if (k > 0) {
          j = k - 1;
          MatSetValues (triDiagA, 1, &k, 1, &j, & sysMatrix[k - 1][k], INSERT_VALUES);
        }

        if (k < nlayers - 1) {
          j = k + 1;
          MatSetValues (triDiagA, 1, &k, 1, &j, &sysMatrix[k][k + 1], INSERT_VALUES);
        }
      }

      ierr = MatAssemblyBegin (triDiagA, MAT_FINAL_ASSEMBLY);
      ierr = MatAssemblyEnd (triDiagA, MAT_FINAL_ASSEMBLY);


//       PetscViewer    viewer;
//       PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//       PetscObjectSetName((PetscObject)viewer,"implicit Euler matrix");
//       PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//       MatView(triDiagA,viewer);
//       double a;
//       std::cin>>a;


      ierr = KSPCreate (PETSC_COMM_WORLD, &solver);
      ierr = KSPSetOperators (solver, triDiagA, triDiagA);
      ierr = KSPSetType (solver, KSPRICHARDSON);
      ierr = KSPSolve (solver, b, x);

      //1. aggiornare solT con x
      for (k = 0; k < nlayers; k++) {
        PetscScalar valueT = 0.;
        ierr = VecGetValues (x, 1, &k, &valueT);

        if (counter == numberOfTimeSteps - 1) {
          std::cout.precision (14);
          std::cout << valueT << std::endl;
        }

        sol->_Sol[solIndexT[k]]->set (i, valueT);
        sol->_Sol[solIndexT[k]]->close();

        //2. aggiornare solHT
        double valueH = (*sol->_Sol[solIndexh[k]]) (i);
        double valueHT = valueH * valueT;
        sol->_Sol[solIndexHT[k]]->set (i, valueHT);
        sol->_Sol[solIndexHT[k]]->close();
      }

      MatDestroy (&triDiagA);
      VecDestroy (&b);
      VecDestroy (&x);
      KSPDestroy (&solver);

    }

  }

  //BEGIN no vertical diffusion
//   for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
//     for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
//       double valueHT = ( *sol->_Sol[solIndexHT[k]] ) ( i );
//       double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );
//
//       double valueT = valueHT / valueH;
//       std::cout.precision(14);
//       if(counter==numberOfTimeSteps-2) std::cout << valueT << std::endl;
//
//       sol->_Sol[solIndexT[k]]->set ( i, valueT );
//     }
//
//     sol->_Sol[solIndexT[k]]->close();
//
//   }
  //END

//   PetscViewer    viewer;
//   PetscViewerDrawOpen(PETSC_COMM_WORLD,NULL,NULL,0,0,900,900,&viewer);
//   PetscObjectSetName((PetscObject)viewer,"FSI matrix");
//   PetscViewerPushFormat(viewer,PETSC_VIEWER_DRAW_LG);
//   MatView((static_cast<PetscMatrix*>(KK))->mat(),viewer);
//   double a;
//   std::cin>>a;

//  abort();

//   unsigned solIndexeta = mlSol->GetIndex ( "eta" );
//   unsigned solIndexb = mlSol->GetIndex ( "b" );
//   sol->_Sol[solIndexeta]->zero();
//   for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
//     sol->_Sol[solIndexeta]->add ( *sol->_Sol[solIndexh[k]] );
//   }
//   sol->_Sol[solIndexeta]->add ( -1, *sol->_Sol[solIndexb] );

}


void RK (MultiLevelProblem& ml_prob, const unsigned & numberOfTimeSteps)
{

  const unsigned& NLayers = NumberOfLayers;

  adept::Stack& s = FemusInit::_adeptStack;

  TransientLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<TransientLinearImplicitSystem> ("SWt");   // pointer to the linear implicit system named "Poisson"

  unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1u;

  Mesh* msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solIndexh (NLayers);

  std::vector < unsigned > solIndexv (NLayers);

  std::vector < unsigned > solIndexHT (NLayers);
  //std::vector < unsigned > solPdeIndexHT ( NLayers );

  std::vector < unsigned > solIndexT (NLayers);

  std::vector < unsigned > solIndexK1 (NLayers);
  std::vector < unsigned > solIndexK2 (NLayers);
  std::vector < unsigned > solIndexK3 (NLayers);
  std::vector < unsigned > solIndexK4 (NLayers);


  for (unsigned i = 0; i < NLayers; i++) {

    char name[10];
    sprintf (name, "h%d", i);
    solIndexh[i] = mlSol->GetIndex (name);   // get the position of "hi" in the sol object

    sprintf (name, "v%d", i);
    solIndexv[i] = mlSol->GetIndex (name);   // get the position of "vi" in the sol object

    sprintf (name, "HT%d", i);
    solIndexHT[i] = mlSol->GetIndex (name);   // get the position of "HTi" in the sol object
    //solPdeIndexHT[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "Ti" in the pdeSys object

    sprintf (name, "T%d", i);
    solIndexT[i] = mlSol->GetIndex (name);   // get the position of "Ti" in the sol object

    sprintf (name, "K1_%d", i);
    solIndexK1[i] = mlSol->GetIndex (name);   // get the position of "K1i" in the sol object

    if (RK_order > 1) {
      sprintf (name, "K2_%d", i);
      solIndexK2[i] = mlSol->GetIndex (name);   // get the position of "K2i" in the sol object

      if (RK_order > 2) {
        sprintf (name, "K3_%d", i);
        solIndexK3[i] = mlSol->GetIndex (name);   // get the position of "K3i" in the sol object

        if (RK_order > 3) {
          sprintf (name, "K4_%d", i);
          solIndexK4[i] = mlSol->GetIndex (name);   // get the position of "K4i" in the sol object
        }
      }
    }

  }

  unsigned solTypeh = mlSol->GetSolutionType (solIndexh[0]);   // get the finite element type for "hi"
  unsigned solTypev = mlSol->GetSolutionType (solIndexv[0]);   // get the finite element type for "vi"
  unsigned solTypeHT = mlSol->GetSolutionType (solIndexHT[0]);   // get the finite element type for "Ti"

//     for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
//         for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
//             //double valueT = ( *sol->_SolOld[solIndexT[k]] ) ( i );
//             double valueT = ( *sol->_Sol[solIndexT[k]] ) ( i );
//             double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );
// 
//             double valueHT = valueT * valueH;
// 
//             sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
//         }
//         sol->_Sol[solIndexHT[k]]->close();
//     }

  std::vector < double > maxW (NLayers, -1.e6);
  maxW[0] = 0.;

  //BEGIN RK loop
  for (unsigned RK_step = 0; RK_step < RK_order; RK_step++) {
    unsigned start = msh->_dofOffset[solTypeHT][iproc];
    unsigned end = msh->_dofOffset[solTypeHT][iproc + 1];

    //BEGIN i loop
    for (unsigned i =  start; i <  end; i++) {

      vector < double > solhm (NLayers);
      vector < double > solh (NLayers);   // local coordinates
      vector < double > solhp (NLayers);
      vector < double > solvm (NLayers);   // local coordinates
      vector < double > solvp (NLayers);   // local coordinates
      vector < double > solTm (NLayers);   // local coordinates
      vector < double > solT (NLayers);   // local coordinates
      vector < double > solTp (NLayers);   // local coordinates

      vector < double > solTmm (NLayers);   // local coordinates
      vector < double > solTpp (NLayers);   // local coordinates
      vector < double > solhmm (NLayers);   // local coordinates
      vector < double > solhpp (NLayers);   // local coordinates

      vector < double > solK1 (NLayers);   // local coordinates
      vector < double > solK1p (NLayers);   // local coordinates
      vector < double > solK1m (NLayers);   // local coordinates
      vector < double > solK1pp (NLayers);   // local coordinates
      vector < double > solK1mm (NLayers);   // local coordinates

      vector < double > solK2 (NLayers);   // local coordinates
      vector < double > solK2p (NLayers);   // local coordinates
      vector < double > solK2m (NLayers);   // local coordinates
      vector < double > solK2pp (NLayers);   // local coordinates
      vector < double > solK2mm (NLayers);   // local coordinates

      vector < double > solK3 (NLayers);   // local coordinates
      vector < double > solK3p (NLayers);   // local coordinates
      vector < double > solK3m (NLayers);   // local coordinates
      vector < double > solK3pp (NLayers);   // local coordinates
      vector < double > solK3mm (NLayers);   // local coordinates

      vector < double > solK4 (NLayers);   // local coordinates
      vector < double > solK4p (NLayers);   // local coordinates
      vector < double > solK4m (NLayers);   // local coordinates
      vector < double > solK4pp (NLayers);   // local coordinates
      vector < double > solK4mm (NLayers);   // local coordinates

//        vector < double > solHTmm ( NLayers ); // local coordinates
//        vector < double > solHTpp ( NLayers ); // local coordinates
//        vector < double > solHT ( NLayers ); // local coordinates
//        vector < double > solHTp ( NLayers ); // local coordinates
//        vector < double > solHTm ( NLayers ); // local coordinates

      unsigned bc1 = (i == start) ? 0 : 1;
      unsigned bc2 = (i == end - 1) ? 0 : 1;
      unsigned bc3 = (i > start + 1) ? 1 : 0;
      unsigned bc4 = (i < end - 2) ? 1 : 0;

      for (unsigned j = 0; j < NLayers; j++) {

        solh[j] = (*sol->_Sol[solIndexh[j]]) (i);
        solT[j] = (*sol->_Sol[solIndexT[j]]) (i);
//         solHT[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i );

        solK1[j] = (*sol->_Sol[solIndexK1[j]]) (i);

        if (RK_order > 1) {
          solK2[j] = (*sol->_Sol[solIndexK2[j]]) (i);

          if (RK_order > 2) {
            solK3[j] = (*sol->_Sol[solIndexK3[j]]) (i);

            if (RK_order > 3) {
              solK4[j] = (*sol->_Sol[solIndexK4[j]]) (i);
            }
          }
        }

        solvm[j] = (*sol->_Sol[solIndexv[j]]) (i);
        solvp[j] = (*sol->_Sol[solIndexv[j]]) (i + 1);

        if (i > start) {
          solhm[j] = (*sol->_Sol[solIndexh[j]]) (i - 1);
          solTm[j] = (*sol->_Sol[solIndexT[j]]) (i - 1);
//           solHTm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 1 );

          solK1m[j] = (*sol->_Sol[solIndexK1[j]]) (i - 1);

          if (RK_order > 1) {
            solK2m[j] = (*sol->_Sol[solIndexK2[j]]) (i - 1);

            if (RK_order > 2) {
              solK3m[j] = (*sol->_Sol[solIndexK3[j]]) (i - 1);

              if (RK_order > 3) {
                solK4m[j] = (*sol->_Sol[solIndexK4[j]]) (i - 1);
              }
            }
          }
        }

        if (i < end - 1) {
          solhp[j] = (*sol->_Sol[solIndexh[j]]) (i + 1);
          solTp[j] = (*sol->_Sol[solIndexT[j]]) (i + 1);
//           solHTp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 1 );

          solK1p[j] = (*sol->_Sol[solIndexK1[j]]) (i + 1);

          if (RK_order > 1) {
            solK2p[j] = (*sol->_Sol[solIndexK2[j]]) (i + 1);

            if (RK_order > 2) {
              solK3p[j] = (*sol->_Sol[solIndexK3[j]]) (i + 1);

              if (RK_order > 3) {
                solK4p[j] = (*sol->_Sol[solIndexK4[j]]) (i + 1);
              }
            }
          }
        }

//         if (i > start + 1) {
//           solTmm[j] = (*sol->_Sol[solIndexT[j]]) (i - 2);
//           solhmm[j] = (*sol->_Sol[solIndexh[j]]) (i - 2);
//           //solHTmm[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i - 2 );
// 
//           solK1mm[j] = (*sol->_Sol[solIndexK1[j]]) (i - 2);
// 
//           if (RK_order > 1) {
//             solK2mm[j] = (*sol->_Sol[solIndexK2[j]]) (i - 2);
// 
//             if (RK_order > 2) {
//               solK3mm[j] = (*sol->_Sol[solIndexK3[j]]) (i - 2);
// 
//               if (RK_order > 3) {
//                 solK4mm[j] = (*sol->_Sol[solIndexK4[j]]) (i - 2);
//               }
//             }
//           }
//         }

//         if (i < end - 2) {
//           solTpp[j] = (*sol->_Sol[solIndexT[j]]) (i + 2);
//           solhpp[j] = (*sol->_Sol[solIndexh[j]]) (i + 2);
//           //solHTpp[j] = ( *sol->_Sol[solIndexHT[j]] ) ( i + 2 );
// 
//           solK1pp[j] = (*sol->_Sol[solIndexK1[j]]) (i + 2);
// 
//           if (RK_order > 1) {
//             solK2pp[j] = (*sol->_Sol[solIndexK2[j]]) (i + 2);
// 
//             if (RK_order > 2) {
//               solK3pp[j] = (*sol->_Sol[solIndexK3[j]]) (i + 2);
// 
//               if (RK_order > 3) {
//                 solK4pp[j] = (*sol->_Sol[solIndexK4[j]]) (i + 2);
//               }
//             }
//           }
//         }

      }

      vector < double > x (2);   // local coordinates

      for (unsigned j = 0; j < 2; j++) {
        unsigned xDof  = msh->GetSolutionDof (j, i, 2);   // global to global mapping between coordinates node and coordinate dof
        x[j] = (*msh->_topology->_Sol[0]) (xDof);     // global extraction and local storage for the element coordinates
      }

      double dx = x[1] - x[0];

      double b = 10.;

      std::vector < double > w (NLayers + 1, 0.);

      std::vector < double > zTop (NLayers);
      zTop[0] = 0;

      for (unsigned k = 1; k < NLayers; k++) {
        zTop[k] = zTop[k - 1] - solh[k];
      }

      std::vector < double > psi2 (NLayers);

      for (unsigned k = 0; k < NLayers; k++) {
        //psi2[k] = ( - ( zMid[k] + 10 ) * zMid[k] ) / 25;
        psi2[k] = 1. - (zTop[k] + 5.) * (zTop[k] + 5.) / (25.);
      }

      double xmid = 0.5 * (x[1] + x[0]);

      for (unsigned k = NLayers; k > 1; k--) {
        w[k - 1] = - (-4. / 625.* (xmid - 5) * (xmid - 5) * (xmid - 5)) * psi2[k - 1];

        //w[k - 1] = - ( - 40./(pow(50.,40)) * pow((xmid - 50.), 39) ) * psi2[k - 1];
        //w[k - 1] = - ( - 16./(pow(10.,16)) * pow((xmid - 10.), 15) ) * psi2[k - 1];
        //w[k - 1] = ( ( 10. - 2. * xmid ) / 25. ) * psi2[k - 1];
        if (maxW[k - 1] < w[k - 1]) {
          maxW[k - 1] = w[k - 1];
        }
      }


      for (unsigned k = 0; k < NLayers; k++) {
        double RHS = 0.;
        double addition = 0.;
        double addition_m = 0.;
        double addition_p = 0.;
//         double addition_mm = 0.;
//         double addition_pp = 0.;
        double addition_t = 0.;
        double addition_b = 0.;

        if (RK_step == 1) {
          if (RK_order == 2) {
            addition = solK1[k];
            addition_m = solK1m[k];
            addition_p = solK1p[k];
//             addition_mm = solK1mm[k];
//             addition_pp = solK1pp[k];
            if (k > 0) addition_t = solK1[k-1]; 
            if (k < NLayers - 1) addition_b = solK1[k+1]; 
          }

          else if (RK_order == 3 || RK_order == 4) {
            addition = solK1[k] * 0.5;
            addition_m = solK1m[k] * 0.5;
            addition_p = solK1p[k] * 0.5;
//             addition_mm = solK1mm[k] * 0.5;
//             addition_pp = solK1pp[k] * 0.5;
            if (k > 0) addition_t = solK1[k-1] * 0.5; 
            if (k < NLayers - 1) addition_b = solK1[k+1] * 0.5;
          }

        }

        else if (RK_step == 2) {
          if (RK_order == 4) {
            addition = solK2[k] * 0.5;
            addition_m = solK2m[k] * 0.5;
            addition_p = solK2p[k] * 0.5;
//             addition_mm = solK2mm[k] * 0.5;
//             addition_pp = solK2pp[k] * 0.5;
            if (k > 0) addition_t = solK2[k-1] * 0.5; 
            if (k < NLayers - 1) addition_b = solK2[k+1] * 0.5;
          }

          else if (RK_order == 3) {
            addition = - solK1[k] + solK2[k] * 2.;
            addition_m = - solK1m[k] + solK2m[k] * 2.;
            addition_p = - solK1p[k] + solK2p[k] * 2.;
/*            addition_mm = - solK1mm[k] + solK2mm[k] * 2.;
            addition_pp = - solK1pp[k] + solK2pp[k] * 2.; */           
            if (k > 0) addition_t = solK2[k-1] * 2.; 
            if (k < NLayers - 1) addition_b = solK2[k+1] * 2.;
          }

        }

        else if (RK_step == 3) {
          addition = solK3[k];
          addition_m = solK3m[k];
          addition_p = solK3p[k];
//           addition_mm = solK3mm[k];
//           addition_pp = solK3pp[k];
          if (k > 0) addition_t = solK3[k-1]; 
          if (k < NLayers - 1) addition_b = solK3[k+1];
        }


        //BEGIN HORIZONTAL ADVECTION


        //BEGIN CENTRAL DIFF h and T separate
//                 if ( i > start ) {
//                      RHS +=  ( 0.5 * ( solhm[k] + solh[k] ) ) * solvm[k] * ( 0.5 * (solTm[k] + addition_m + solT[k] + addition ) ) / dx;
//                 }
//                 if ( i < end - 1 ) {
//                       RHS -= ( 0.5 * ( solh[k] + solhp[k] ) ) * solvp[k] * ( 0.5 * (solT[k] + addition + solTp[k] + addition_p ) ) / dx;
//                 }
        //END


        //BEGIN CENTRAL DIFF hT multiplied
//                         if ( i > start ) {
//                             RHS += 0.5 * ( solHTm[k] + addition_m + solHT[k] + addition ) * solvm[k]  / dx;
//                         }
//                         if ( i < end - 1 ) {
//                             RHS -= 0.5 * (solHT[k] + addition + solHTp[k] + addition_p ) * solvp[k]  / dx;
//                         }
        //END


        //BEGIN FIRST ORDER UPWINDING h and T separate
        if (i > start) {
          if (solvm[k] > 0) {
            RHS += ( (solTm[k] + addition_m) * solhm[k]) * solvm[k] / dx;
          }
          else {
            RHS += ( (solT[k] + addition) * solh[k]) * solvm[k] / dx;
          }
        }

        if (i < end - 1) {
          if (solvp[k] > 0) {
            RHS -= ( (solT[k] + addition) * solh[k]) * solvp[k] / dx;
          }
          else {
            RHS -= ( (solTp[k] + addition_p) * solhp[k]) * solvp[k] / dx;
          }
        }
        //END


        //BEGIN FIRST ORDER UPWINDING hT multiplied
//                 if ( i > start ) {
//                     if ( solvm[k] > 0 ) {
//                         RHS += ( solHTm[k] + addition_m ) * solvm[k] / dx;
//                     } 
//                     else {
//                         RHS += ( solHT[k] + addition ) * solvm[k] / dx;
//                     }
//                 }
//                 if ( i < end - 1 ) {
//                     if ( solvp[k] > 0 ) {
//                         RHS -= ( solHT[k] + addition ) * solvp[k] / dx;
//                     } 
//                     else {
//                         RHS -= ( solHTp[k] + addition_p ) * solvp[k] / dx;
//                     }
//                 }
        //END


        //BEGIN THIRD ORDER UPWINDING h and T separate

//                 if ( i > start ) {
//                   RHS += 0.5 * ( (solTm[k] + addition_m) * solhm[k] + (solT[k] + addition) * solh[k] ) * solvm[k] / dx;
//                     if ( solvm[k] > 0 ) {
//                       if ( i > start + 1 ) {
//                         RHS += - 1. / 6. * ( (solT[k] + addition) * solh[k] - 2.* (solTm[k] + addition_m) * solhm[k] + (solTmm[k] + addition_mm) * solhmm[k] ) *
//                                solvm[k]  / dx;
//                         }
//                     } 
//                     else {
//                       if ( i < end - 1 ) {
//                           RHS += - 1. / 6. * ( (solTp[k] + addition_p) * solhp[k]  - 2.* (solT[k] + addition) * solh[k] + (solTm[k] + addition_m) * solhm[k] ) * 
//                                  solvm[k]  / dx;
//                       }
//                     }
//                 }
//                 if ( i < end - 1 ) {
//                   RHS -= 0.5 * ( (solTp[k] + addition_p) * solhp[k] + (solT[k] + addition) * solh[k] ) * solvp[k] / dx;
//                     if ( solvp[k] > 0 ) {
//                       if ( i > start ) {
//                         RHS -= - 1. / 6. * ( (solTp[k] + addition_p) * solhp[k] - 2.* (solT[k] + addition) * solh[k] + (solTm[k] + addition_m) * solhm[k] ) * solvp[k]
//                                / dx;
//                         }
//                     } 
//                     else {
//                       if ( i < end - 2 ) {
//                         RHS -= - 1. / 6. * ( (solTpp[k] + addition_pp) * solhpp[k] - 2.* (solTp[k] + addition_p) * solhp[k] + (solT[k] + addition) * solh[k] ) * 
//                                solvp[k]  / dx;
//                         }
//                     }
//                 }

        //END


        //END HORIZONTAL ADVECTION



        //BEGIN VERTICAL ADVECTION


        //BEGIN CENTRAL DIFF h and T separate
//                 if ( k < NLayers - 1 ) {
//                    RHS += w[k + 1] * 0.5 * ( solT[k] + addition + solT[k + 1] + addition_b );
//                 }
//                 if ( k > 0 ) {
//                   RHS -= w[k] * 0.5 * ( solT[k - 1] + addition_t + solT[k] + addition ) ;
//                 }
        //END

        //BEGIN CENTRAL DIFF hT multiplied
//                 if ( k < NLayers - 1 ) {
//                   RHS += w[k + 1] * 0.5 * ( ( solHT[k] + addition ) / solh[k] + ( solHT[k + 1] + addition_b ) / solh[k + 1] );
//                 }
//                 if ( k > 0 ) {
//                   RHS -= w[k] * 0.5 * ( ( solHT[k - 1] + addition_t ) / solh[k - 1] + ( solHT[k] + addition ) / solh[k] );
//                 }
        //END

        //BEGIN FIRST ORDER UPWIND h and T separate
        if (k < NLayers - 1) {
          if (w[k + 1] > 0) {
            RHS += w[k + 1] * (solT[k + 1] + addition_b) ;
          }
          else {
            RHS += w[k + 1] * (solT[k] + addition) ;
          }
        }

        if (k > 0) {
          if (w[k] > 0) {
            RHS -= w[k] * (solT[k] + addition);
          }
          else {
            RHS -= w[k] * (solT[k - 1] + addition_t);
          }
        }
        //END

        //BEGIN FIRST ORDER UPWIND hT multiplied
//                 if ( k < NLayers - 1 ) {
//                   if ( w[k + 1] > 0 ) {
//                     RHS += w[k + 1] * ( ( solHT[k + 1] + addition_b ) / solh[k + 1]  );
//                     } 
//                     else {
//                       RHS += w[k + 1] * ( ( solHT[k] + addition ) / solh[k] );
//                     }
//                 }
//                 if ( k > 0 ) {
//                   if ( w[k] > 0 ) {
//                     RHS -= w[k] * ( ( solHT[k] + addition ) / solh[k] );
//                     } 
//                     else {
//                       RHS -= w[k] * ( ( solHT[k - 1] + addition_t ) / solh[k - 1] );
//                     }
//                 }
        //END



        //END VERTICAL ADVECTION

        // horizontal diffussion
        if (i > start) {
          RHS -= k_h * ( (0.5 * (solhm[k] + solh[k])) * ( (solT[k] + addition)  - (solTm[k] + addition_m)) / dx) / dx;
//           RHS += k_h * (0.5 * (solhm[k] + solh[k])) * ((solHTm[k] / solhm[k] + addition_m)  - (solHT[k] / solh[k] + addition))/(dx*dx); 
        }

        if (i < end - 1) {
          RHS += k_h * ( (0.5 * (solh[k] + solhp[k] ))  * ( (solTp[k] + addition_p)  - (solT[k] + addition)) / dx) / dx; 
//           RHS += k_h * (0.5 * (solhp[k] + solh[k])) * ((solHTp[k] / solhp[k] + addition_p)  - (solHT[k] / solh[k] + addition))/(dx*dx);
        }


        if (RK_step == 0) {
          double valueK1 = RHS * dt;
          sol->_Sol[solIndexK1[k]]->set (i, valueK1);
          sol->_Sol[solIndexK1[k]]->close();
        }
        else if (RK_step == 1) {
          double valueK2 = RHS * dt;
          sol->_Sol[solIndexK2[k]]->set (i, valueK2);
          sol->_Sol[solIndexK2[k]]->close();
        }
        else if (RK_step == 2) {
          double valueK3 = RHS * dt;
          sol->_Sol[solIndexK3[k]]->set (i, valueK3);
          sol->_Sol[solIndexK3[k]]->close();
        }
        else {
          double valueK4 = RHS * dt;
          sol->_Sol[solIndexK4[k]]->set (i, valueK4);
          sol->_Sol[solIndexK4[k]]->close();
        }

      }
      
      //END layer loop

    }

    //END i loop


  }

  //END RK step loop


  //BEGIN UPDATE SOL
  for (unsigned k = 0; k < NumberOfLayers; k++) {
    for (unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++) {

//       double valueHT = 0.; 
      double valueT = 0.;

      if (RK_order == 1) {
        double h = (*sol->_Sol[solIndexh[k]]) (i);
        double K1 = (*sol->_Sol[solIndexK1[k]]) (i);

        double Told = (*sol->_Sol[solIndexT[k]]) (i);
//         double HTold = ( *sol->_Sol[solIndexHT[k]] ) ( i );

        valueT = Told + K1 / h;
//         valueHT = HTold + K1 ;
//         valueT = valueHT / h;
      }

      else if (RK_order == 2) {
        double h = (*sol->_Sol[solIndexh[k]]) (i);
        double K1 = (*sol->_Sol[solIndexK1[k]]) (i);
        double K2 = (*sol->_Sol[solIndexK2[k]]) (i);

        double Told = (*sol->_Sol[solIndexT[k]]) (i);
//         double HTold = ( *sol->_Sol[solIndexHT[k]] ) ( i );

        valueT = Told  + (0.5 * (K1 + K2)) / h;
//         valueHT = HTold + 0.5 * ( K1 + K2 );
//         valueT = valueHT / h;
      }

      else if (RK_order == 3) {
        double h = (*sol->_Sol[solIndexh[k]]) (i);
        double K1 = (*sol->_Sol[solIndexK1[k]]) (i);
        double K2 = (*sol->_Sol[solIndexK2[k]]) (i);
        double K3 = (*sol->_Sol[solIndexK3[k]]) (i);

        double Told = (*sol->_Sol[solIndexT[k]]) (i);
//         double HTold = ( *sol->_Sol[solIndexHT[k]] ) ( i );

        valueT = Told + (1. / 6. * (K1 + 4. * K2 + K3)) / h;
//         valueHT = HTold + ( 1. / 6. * ( K1 + 4. * K2 + K3 ) ) ;
//         valueT = valueHT / h;
      }

      else if (RK_order == 4) {
        double h = (*sol->_Sol[solIndexh[k]]) (i);
        double K1 = (*sol->_Sol[solIndexK1[k]]) (i);
        double K2 = (*sol->_Sol[solIndexK2[k]]) (i);
        double K3 = (*sol->_Sol[solIndexK3[k]]) (i);
        double K4 = (*sol->_Sol[solIndexK4[k]]) (i);

        double Told = (*sol->_Sol[solIndexT[k]]) (i);
//         double HTold = ( *sol->_Sol[solIndexHT[k]] ) ( i );

        valueT = Told + (1. / 6. * (K1 + 2.*K2 + 2.*K3 + K4)) /  h;
//         valueHT = HTold + 1. / 6. * ( K1 + 2.*K2 + 2.*K3 + K4 );
//         valueT = valueHT / h;
      }

//       sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
      sol->_Sol[solIndexT[k]]->set (i, valueT);

//       sol->_Sol[solIndexHT[k]]->close();
      sol->_Sol[solIndexT[k]]->close();
    }
  }

  //END UPDATE SOL

  //BEGIN no vertical diffusion
  for (unsigned k = 0; k < NumberOfLayers; k++) {
    for (unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++) {
      double valueT = (*sol->_Sol[solIndexT[k]]) (i);

      std::cout.precision (14);

      if (counter == numberOfTimeSteps - 1) {
        std::cout << valueT << std::endl;
      }
    }
  }

  //END


}

void RKe (MultiLevelProblem& ml_prob, const unsigned & numberOfTimeSteps)
{

  const unsigned& NLayers = NumberOfLayers;

  adept::Stack& s = FemusInit::_adeptStack;

  TransientLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<TransientLinearImplicitSystem> ("SWt");   // pointer to the linear implicit system named "Poisson"

  unsigned level = ml_prob._ml_msh->GetNumberOfLevels() - 1u;

  Mesh* msh = ml_prob._ml_msh->GetLevel (level);   // pointer to the mesh (level) object
  elem* el = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution* mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution* sol = ml_prob._ml_sol->GetSolutionLevel (level);   // pointer to the solution (level) object

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem

  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
  unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


  //solution variable
  std::vector < unsigned > solIndexh (NLayers);
  std::vector < unsigned > solIndexv (NLayers);
  std::vector < unsigned > solIndexHT (NLayers);
  //std::vector < unsigned > solPdeIndexHT ( NLayers );

  std::vector < unsigned > solIndexT (NLayers);

  std::vector < std::vector < unsigned > > solIndexK(RK_order);
  for(unsigned i = 0; i < RK_order; i++){
    solIndexK[i].resize(NLayers);  
  }
  
  for (unsigned i = 0; i < NLayers; i++) {

    char name[10];
    sprintf (name, "h%d", i);
    solIndexh[i] = mlSol->GetIndex (name);   // get the position of "hi" in the sol object

    sprintf (name, "v%d", i);
    solIndexv[i] = mlSol->GetIndex (name);   // get the position of "vi" in the sol object

    sprintf (name, "HT%d", i);
    solIndexHT[i] = mlSol->GetIndex (name);   // get the position of "HTi" in the sol object
    //solPdeIndexHT[i] = mlPdeSys->GetSolPdeIndex ( name ); // get the position of "Ti" in the pdeSys object

    sprintf (name, "T%d", i);
    solIndexT[i] = mlSol->GetIndex (name);   // get the position of "Ti" in the sol object

    for(unsigned j = 0; j < RK_order; j++){
      sprintf (name, "K%d_%d",j+1 , i);
      solIndexK[j][i] = mlSol->GetIndex (name);   // get the position of "K1i" in the sol object
    }

  }

  unsigned solTypeh = mlSol->GetSolutionType (solIndexh[0]);   // get the finite element type for "hi"
  unsigned solTypev = mlSol->GetSolutionType (solIndexv[0]);   // get the finite element type for "vi"
  unsigned solTypeHT = mlSol->GetSolutionType (solIndexHT[0]);   // get the finite element type for "Ti"

//     for ( unsigned k = 0; k < NumberOfLayers; k++ ) {
//         for ( unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++ ) {
//             //double valueT = ( *sol->_SolOld[solIndexT[k]] ) ( i );
//             double valueT = ( *sol->_Sol[solIndexT[k]] ) ( i );
//             double valueH = ( *sol->_Sol[solIndexh[k]] ) ( i );
// 
//             double valueHT = valueT * valueH;
// 
//             sol->_Sol[solIndexHT[k]]->set ( i, valueHT );
//         }
//         sol->_Sol[solIndexHT[k]]->close();
//     }

  std::vector < double > maxW (NLayers, -1.e6);
  maxW[0] = 0.;

  //BEGIN RK loop
  for (unsigned RK_step = 0; RK_step < RK_order; RK_step++) {
    unsigned start = msh->_dofOffset[solTypeHT][iproc];
    unsigned end = msh->_dofOffset[solTypeHT][iproc + 1];

    //BEGIN i loop
    for (unsigned i =  start; i <  end; i++) {

      vector < double > solhm (NLayers);
      vector < double > solh (NLayers);   // local coordinates
      vector < double > solhp (NLayers);
      vector < double > solvm (NLayers);   // local coordinates
      vector < double > solvp (NLayers);   // local coordinates
      vector < double > solTm (NLayers);   // local coordinates
      vector < double > solT (NLayers);   // local coordinates
      vector < double > solTp (NLayers);   // local coordinates

      vector < double > solTmm (NLayers);   // local coordinates
      vector < double > solTpp (NLayers);   // local coordinates
      vector < double > solhmm (NLayers);   // local coordinates
      vector < double > solhpp (NLayers);   // local coordinates

      vector < vector < double > > solK(RK_step);   // local coordinates
      vector < vector < double > > solKp(RK_step);   // local coordinates
      vector < vector < double > > solKm(RK_step);   // local coordinates
      vector < vector < double > > solKpp(RK_step);   // local coordinates
      vector < vector < double > > solKmm(RK_step);   // local coordinates
      
      for(unsigned j = 0; j < RK_step; j++){
        solK[j].resize(NLayers);  
        solKm[j].resize(NLayers);  
        solKp[j].resize(NLayers);  
        solKpp[j].resize(NLayers);  
        solKmm[j].resize(NLayers);  
      }
    
      unsigned bc1 = (i == start) ? 0 : 1;
      unsigned bc2 = (i == end - 1) ? 0 : 1;
      unsigned bc3 = (i > start + 1) ? 1 : 0;
      unsigned bc4 = (i < end - 2) ? 1 : 0;

      for (unsigned k = 0; k < NLayers; k++) {

        solh[k] = (*sol->_Sol[solIndexh[k]]) (i);
        solT[k] = (*sol->_Sol[solIndexT[k]]) (i);

        for(unsigned j = 0; j < RK_step; j++){
          solK[j][k] = (*sol->_Sol[solIndexK[j][k]]) (i);
        }

        solvm[k] = (*sol->_Sol[solIndexv[k]]) (i);
        solvp[k] = (*sol->_Sol[solIndexv[k]]) (i + 1);

        if (i > start) {
          solhm[k] = (*sol->_Sol[solIndexh[k]]) (i - 1);
          solTm[k] = (*sol->_Sol[solIndexT[k]]) (i - 1);

          for(unsigned j = 0; j < RK_step; j++){
            solKm[j][k] = (*sol->_Sol[solIndexK[j][k]]) (i - 1);
          }
        }

        if (i < end - 1) {
          solhp[k] = (*sol->_Sol[solIndexh[k]]) (i + 1);
          solTp[k] = (*sol->_Sol[solIndexT[k]]) (i + 1);
          
          for(unsigned j = 0; j < RK_step; j++){
            solKp[j][k] = (*sol->_Sol[solIndexK[j][k]]) (i + 1);
          }
        }
      }

      vector < double > x (2);   // local coordinates

      for (unsigned j = 0; j < 2; j++) {
        unsigned xDof  = msh->GetSolutionDof (j, i, 2);   // global to global mapping between coordinates node and coordinate dof
        x[j] = (*msh->_topology->_Sol[0]) (xDof);     // global extraction and local storage for the element coordinates
      }

      double dx = x[1] - x[0];

      double b = 10.;

      std::vector < double > w (NLayers + 1, 0.);

      std::vector < double > zTop (NLayers);
      zTop[0] = 0;

      for (unsigned k = 1; k < NLayers; k++) {
        zTop[k] = zTop[k - 1] - solh[k];
      }

      std::vector < double > psi2 (NLayers);

      for (unsigned k = 0; k < NLayers; k++) {
        //psi2[k] = ( - ( zMid[k] + 10 ) * zMid[k] ) / 25;
        psi2[k] = 1. - (zTop[k] + 5.) * (zTop[k] + 5.) / (25.);
      }

      double xmid = 0.5 * (x[1] + x[0]);

      for (unsigned k = NLayers; k > 1; k--) {
        w[k - 1] = - (-4. / 625.* (xmid - 5) * (xmid - 5) * (xmid - 5)) * psi2[k - 1];

        //w[k - 1] = - ( - 40./(pow(50.,40)) * pow((xmid - 50.), 39) ) * psi2[k - 1];
        //w[k - 1] = - ( - 16./(pow(10.,16)) * pow((xmid - 10.), 15) ) * psi2[k - 1];
        //w[k - 1] = ( ( 10. - 2. * xmid ) / 25. ) * psi2[k - 1];
        if (maxW[k - 1] < w[k - 1]) {
          maxW[k - 1] = w[k - 1];
        }
      }


      for (unsigned k = 0; k < NLayers; k++) {
        double RHS = 0.;
        double addition = 0.;
        double addition_m = 0.;
        double addition_p = 0.;
        double addition_t = 0.;
        double addition_b = 0.;
        
        double a[4][3][3]={
            {{}},
            {{1.}},
            {{0.5},{-1.,2.}},
            {{0.5},{0.,0.5},{0.,0.,1.}}
        };
        
        for(unsigned j = 0; j < RK_step; j++){
          addition += a[RK_order - 1][RK_step - 1][j] * solK[j][k]/solh[k];
          if(i > start)        addition_m += a[RK_order - 1][RK_step - 1][j] * solKm[j][k] / solhm[k];
          if(i < end - 1)      addition_p += a[RK_order - 1][RK_step - 1][j] * solKp[j][k] / solhp[k];
          if (k > 0)           addition_t += a[RK_order - 1][RK_step - 1][j] * solK[j][k-1] / solh[k-1];
          if (k < NLayers - 1) addition_b += a[RK_order - 1][RK_step - 1][j] * solK[j][k+1] / solh[k+1];
        }
              
        //BEGIN HORIZONTAL ADVECTION
      
        //BEGIN FIRST ORDER UPWINDING h and T separate
        if (i > start) {
          if (solvm[k] > 0) {
            RHS += ( (solTm[k] + addition_m) * solhm[k]) * solvm[k] / dx;
          }
          else {
            RHS += ( (solT[k] + addition) * solh[k]) * solvm[k] / dx;
          }
        }

        if (i < end - 1) {
          if (solvp[k] > 0) {
            RHS -= ( (solT[k] + addition) * solh[k]) * solvp[k] / dx;
          }
          else {
            RHS -= ( (solTp[k] + addition_p) * solhp[k]) * solvp[k] / dx;
          }
        }
        //END

        //END HORIZONTAL ADVECTION
        

        //BEGIN VERTICAL ADVECTION
        
        //BEGIN FIRST ORDER UPWIND h and T separate
        if (k < NLayers - 1) {
          if (w[k + 1] > 0) {
            RHS += w[k + 1] * (solT[k + 1] + addition_b) ;
          }
          else {
            RHS += w[k + 1] * (solT[k] + addition) ;
          }
        }

        if (k > 0) {
          if (w[k] > 0) {
            RHS -= w[k] * (solT[k] + addition);
          }
          else {
            RHS -= w[k] * (solT[k - 1] + addition_t);
          }
        }
        //END

        //END VERTICAL ADVECTION

        // horizontal diffussion
        if (i > start) {
          RHS -= k_h * ( (0.5 * (solhm[k] + solh[k])) * ( (solT[k] + addition)  - (solTm[k] + addition_m)) / dx) / dx;
        }

        if (i < end - 1) {
          RHS += k_h * ( (0.5 * (solh[k] + solhp[k] ))  * ( (solTp[k] + addition_p)  - (solT[k] + addition)) / dx) / dx; 
        }

        sol->_Sol[solIndexK[RK_step][k]]->set (i, RHS * dt);
      }
      //END layer loop
    } 
    //END i loop
    for (unsigned k = 0; k < NLayers; k++) {
      sol->_Sol[solIndexK[RK_step][k]]->close();
    }
  }

  //END RK step loop


  //BEGIN UPDATE SOL
  
  double b[4][4] = {{1.},{0.5, 0.5}, {1./6., 2./3., 1./6.},{1./6., 1./3., 1./3., 1./6.}};
  for (unsigned k = 0; k < NumberOfLayers; k++) {
    for (unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++) {

      double valueT = 0.;
      double h = (*sol->_Sol[solIndexh[k]]) (i); 
      for(unsigned j = 0 ; j < RK_order; j++){
        valueT += b[RK_order - 1][j] * (*sol->_Sol[solIndexK[j][k]]) (i) / h;
      }
      sol->_Sol[solIndexT[k]]->add (i, valueT);
    }
    sol->_Sol[solIndexT[k]]->close();
  }

  //END UPDATE SOL

  //BEGIN no vertical diffusion
  for (unsigned k = 0; k < NumberOfLayers; k++) {
    for (unsigned i =  msh->_dofOffset[solTypeHT][iproc]; i <  msh->_dofOffset[solTypeHT][iproc + 1]; i++) {
      if (counter == numberOfTimeSteps - 1) {
        double valueT = (*sol->_Sol[solIndexT[k]]) (i);
        std::cout.precision (14);  
        std::cout << valueT << std::endl<<std::flush;
      }
    }
  }

  //END


}
