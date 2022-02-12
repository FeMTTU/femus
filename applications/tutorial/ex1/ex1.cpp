/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution variables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"

using namespace femus;

double InitalValueU(const std::vector < double >& x) {
  return x[0] + x[1];
}

double InitalValueP(const std::vector < double >& x) {
  return x[0];
}

double InitalValueT(const std::vector < double >& x) {
  return x[1];
}

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh; // Consider the "mlMsh" object of "MultiLevel" class.
  double scalingFactor = 1.;
  // read coarse level mesh and generate finers level meshes 
  mlMsh.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor); // Let mlMsh read the coarse mesh.
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 3; // We apply uniform refinement.
  unsigned numberOfSelectiveLevels = 0; // we may want to see solutions on different level of meshes. 
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL); // Add those refined meshed in mlMsh object. 
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh); // Define "mlSol" object with argument mlMsh of "MultiLevelSolution" class. 
  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, FIRST);
  mlSol.AddSolution("V", LAGRANGE, SERENDIPITY);
  mlSol.AddSolution("W", LAGRANGE, SECOND);
  mlSol.AddSolution("P", DISCONTINUOUS_POLYNOMIAL, ZERO);
  mlSol.AddSolution("T", DISCONTINUOUS_POLYNOMIAL, FIRST);

  mlSol.Initialize("All");    // initialize all varaibles to zero

  mlSol.Initialize("U", InitalValueU);
  mlSol.Initialize("P", InitalValueP);
  mlSol.Initialize("T", InitalValueT);    // note that this initialization is the same as piecewise constant element

  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("U");
  variablesToBePrinted.push_back("P");
  variablesToBePrinted.push_back("T");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  GMVWriter gmvIO(&mlSol);
  variablesToBePrinted.push_back("all");
  gmvIO.SetDebugOutput(false);
  gmvIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  return 0;
}




