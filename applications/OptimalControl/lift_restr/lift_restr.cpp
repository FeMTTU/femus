#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "LinearImplicitSystem.hpp"

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

bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
  bool dirichlet = true; //dirichlet
  value = 0;

  if (faceName == 2)
    dirichlet = false;

  return dirichlet;
}



int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  // read coarse level mesh and generate finers level meshes
  mlMsh.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 3;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("Thom", LAGRANGE, SECOND);
  mlSol.AddSolution("ThomAdj", LAGRANGE, SECOND);
  mlSol.AddSolution("Tcont", LAGRANGE, SECOND);

  mlSol.Initialize("All");    // initialize all varaibles to zero

//   mlSol.Initialize("Thom", InitalValueU);
//   mlSol.Initialize("ThomAdj", InitalValueP);
//   mlSol.Initialize("Tcont", InitalValueT);    // note that this initialization is the same as piecewise constant element
 
  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("Thom");
  mlSol.GenerateBdc("ThomAdj");
  mlSol.GenerateBdc("Tcont");

  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  
 // add system  in mlProb as a Linear Implicit System
  LinearImplicitSystem& system = mlProb.add_system < LinearImplicitSystem > ("Optsys");
 
  system.AddSolutionToSystemPDE("Thom");  
  system.AddSolutionToSystemPDE("ThomAdj");  
  system.AddSolutionToSystemPDE("Tcont");  
  
  
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("Thom");
  variablesToBePrinted.push_back("ThomAdj");
  variablesToBePrinted.push_back("Tcont");

  VTKWriter vtkIO(&mlSol);
  vtkIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  GMVWriter gmvIO(&mlSol);
  variablesToBePrinted.push_back("all");
  gmvIO.SetDebugOutput(false);
  gmvIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  return 0;
}