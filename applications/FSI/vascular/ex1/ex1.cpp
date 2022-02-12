#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"

using namespace femus;

bool SetBoundaryCondition(const std::vector < double >& x, const char solName[], double& value, const int faceName, const double time) {
  bool dirichlet = false; //dirichlet
  value = 0;

  if(!strcmp(solName,"U")){
    if (faceName == 1){
      dirichlet = true;
    }
  }
  else if(!strcmp(solName,"V")){
    if (faceName == 2){
      dirichlet = true;
    }
  }
  else if(!strcmp(solName,"W")){
    if (faceName == 3){
      dirichlet = true;
    }
  }
  else if(!strcmp(solName,"X")){
    if (faceName == 4){
      dirichlet = true;
    }
  }
  else if(!strcmp(solName,"Y")){
    if (faceName == 5){
      dirichlet = true;
    }
  }
  else if(!strcmp(solName,"Z")){
    if (faceName == 6){
      dirichlet = true;
    }
  }
  else if(!strcmp(solName,"A")){
    if (faceName == 7){
      dirichlet = true;
    }
  }
  else if(!strcmp(solName,"B")){
    if (faceName == 8){
      dirichlet = true;
    }
  }
  else if(!strcmp(solName,"C")){
    if (faceName == 9){
      dirichlet = true;
    }
  }
  else if(!strcmp(solName,"D")){
    if (faceName == 10){
      dirichlet = true;
    }
  }
  else if(!strcmp(solName,"E")){
    if (faceName == 11){
      dirichlet = true;
    }
  }

  return dirichlet;
}

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  // read coarse level mesh and generate finers level meshes
  //mlMsh.ReadCoarseMesh("./input/aneurysm.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/junction.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/aneurysm_one_ball.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/aneurisma_aorta_9Oscaled.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("./input/Turek_stents.neu", "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the furure it is not going to be an argument of this function   */
  unsigned numberOfUniformLevels = 5;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  mlSol.AddSolution("W", LAGRANGE, SECOND);
  mlSol.AddSolution("X", LAGRANGE, SECOND);
  mlSol.AddSolution("Y", LAGRANGE, SECOND);
  mlSol.AddSolution("Z", LAGRANGE, SECOND);
  mlSol.AddSolution("A", LAGRANGE, SECOND);
  mlSol.AddSolution("B", LAGRANGE, SECOND);
  mlSol.AddSolution("C", LAGRANGE, SECOND);
  mlSol.AddSolution("D", LAGRANGE, SECOND);
  mlSol.AddSolution("E", LAGRANGE, SECOND);

  mlSol.Initialize("All");    // initialize all varaibles to zero

  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("U");
  mlSol.GenerateBdc("V");
  mlSol.GenerateBdc("W");
  mlSol.GenerateBdc("X");
  mlSol.GenerateBdc("Y");
  mlSol.GenerateBdc("Z");
  mlSol.GenerateBdc("A");
  mlSol.GenerateBdc("B");
  mlSol.GenerateBdc("C");
  mlSol.GenerateBdc("D");
  mlSol.GenerateBdc("E");


  // print solutions
  std::vector < std::string > variablesToBePrinted;
  /*variablesToBePrinted.push_back("U");
  variablesToBePrinted.push_back("V");
  variablesToBePrinted.push_back("W");
  variablesToBePrinted.push_back("X");
  variablesToBePrinted.push_back("Y");
  variablesToBePrinted.push_back("Z");
  */
  variablesToBePrinted.push_back("All");


  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

    return 0;
}
