 #include "MultiLevelProblem.hpp"
// #include "TransientSystem.hpp"
// #include "NumericVector.hpp"
// #include "Fluid.hpp"
// #include "Parameter.hpp"

// #include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
// #include "NonLinearImplicitSystem.hpp"
// #include "FElemTypeEnum.hpp"
// #include "Files.hpp"

// using std::cout;
// using std::endl;


#include "FemusInit.hpp"
using namespace femus;

double InitVariableU(const double &x, const double &y, const double &z){
  return x+y;
}

double InitVariableP(const double &x, const double &y, const double &z){
  return x;
}

double InitVariableT(const double &x, const double &y, const double &z){
  return y;
}

int main(int argc,char **args) {
  
  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);
 
  // read coarse level mesh and generate finers level meshes
  MultiLevelMesh mlMsh;
  double scalingFactor=1.; 
  mlMsh.ReadCoarseMesh("./input/square.neu","seventh",scalingFactor); 
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfLevels=2;
  unsigned numberOfSelectiveLevels=0;
  mlMsh.RefineMesh(numberOfLevels , numberOfLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();
  
  // define and initialize variables
  MultiLevelSolution mlSol(&mlMsh);
  
  mlSol.AddSolution("U",LAGRANGE, FIRST);
  mlSol.AddSolution("V",LAGRANGE, SERENDIPITY);
  mlSol.AddSolution("W",LAGRANGE, SECOND);
  mlSol.AddSolution("P",DISCONTINOUS_POLYNOMIAL, ZERO);
  mlSol.AddSolution("T",DISCONTINOUS_POLYNOMIAL, FIRST);
  
  mlSol.Initialize("All");  
  mlSol.Initialize("U",InitVariableU);
  mlSol.Initialize("P",InitVariableP);
  mlSol.Initialize("T",InitVariableT); //note that this initialization is the same as piecewise constant element
  
  // Print solutions
  std::vector<std::string> variablesToBePrinted;
  variablesToBePrinted.push_back("U");
  variablesToBePrinted.push_back("P");
  variablesToBePrinted.push_back("T");

  VTKWriter vtkIO(mlSol);
  vtkIO.write_system_solutions(DEFAULT_OUTPUTDIR,"biquadratic",variablesToBePrinted);

  GMVWriter gmvIO(mlSol);
  gmvIO.write_system_solutions(DEFAULT_OUTPUTDIR,"biquadratic",variablesToBePrinted);
   
  return 1;
}

