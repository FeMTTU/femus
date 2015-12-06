/** tutorial/Ex1
 * Project 1
 * In this project Navier-Stokes equation will be solved
 * The problem is a pressurized hydraulic fracture 
 * Pressure distribution du to hydauli load will be solved as a function 
 **/

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"

using namespace femus;

double InitalValueU(const std::vector < double >& x) {
  return x[0] + x[1];
}


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh ml_msh;
  double scalingFactor = 1.;
  // read coarse level mesh and generate finers level meshes
    ml_msh.ReadCoarseMesh("./input/hydr_frac.med", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
    ml_msh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
    ml_msh.PrintInfo();

  ml_msh.SetWriter(VTK);
  ml_msh.GetWriter()->write(DEFAULT_OUTPUTDIR,"biquadratic");

  
  
//   // define the multilevel solution and attach the mlMsh object to it
//   MultiLevelSolution mlSol(&mlMsh);
// 
//   // add variables to mlSol
//   mlSol.AddSolution("u", LAGRANGE, FIRST);
//  
//   mlSol.Initialize("All");    // initialize all varaibles to zero
// 
//   mlSol.Initialize("u", InitalValueU);
//  
//   // print solutions
//   std::vector < std::string > variablesToBePrinted;
//   variablesToBePrinted.push_back("u");
// 
//   VTKWriter vtkIO(&mlSol);
//   vtkIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);


  return 0;
}




