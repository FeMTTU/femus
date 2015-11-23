
#include <sstream>
#include "FemusDefault.hpp"
#include "FemusInit.hpp"
#include "MultiLevelMesh.hpp"
#include "WriterEnum.hpp"

using namespace femus;

// Test for SalomeIO reading


int main(int argc,char **args) {

  FemusInit init(argc,args,MPI_COMM_WORLD);
  
  std::string med_file = "RectFracDom.med";
  std::ostringstream mystream; mystream << "./" << DEFAULT_INPUTDIR << "/" << med_file;
  const std::string infile = mystream.str();
 
  //Adimensional
  double Lref = 1.;
  
  MultiLevelMesh ml_msh;
  ml_msh.ReadCoarseMesh(infile.c_str(),"fifth",Lref);
  
  ml_msh.SetWriter(XDMF);
  ml_msh.GetWriter()->write(DEFAULT_OUTPUTDIR,"biquadratic");
  ml_msh.SetWriter(GMV);
  ml_msh.GetWriter()->write(DEFAULT_OUTPUTDIR,"biquadratic");
  ml_msh.SetWriter(VTK);
  ml_msh.GetWriter()->write(DEFAULT_OUTPUTDIR,"biquadratic");

  return 0;
}

















































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

// // // #include "FemusInit.hpp"
// // // #include "MultiLevelProblem.hpp"
// // // #include "VTKWriter.hpp"
// // // #include "GMVWriter.hpp"
// // // 
// // // using namespace femus;
// // // 
// // // double InitalValueU(const std::vector < double >& x) {
// // //   return x[0] + x[1];
// // // }
// // // 
// // // double InitalValueP(const std::vector < double >& x) {
// // //   return x[0];
// // // }
// // // 
// // // double InitalValueT(const std::vector < double >& x) {
// // //   return x[1];
// // // }
// // // 
// // // int main(int argc, char** args) {
// // // 
// // //   // init Petsc-MPI communicator
// // //   FemusInit mpinit(argc, args, MPI_COMM_WORLD);
// // // 
// // //   // define multilevel mesh
// // //   MultiLevelMesh mlMsh;
// // //   double scalingFactor = 1.;
// // //   // read coarse level mesh and generate finers level meshes
// // //   mlMsh.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor);
// // //   /* "seventh" is the order of accuracy that is used in the gauss integration scheme
// // //       probably in the furure it is not going to be an argument of this function   */
// // //   unsigned numberOfUniformLevels = 3;
// // //   unsigned numberOfSelectiveLevels = 0;
// // //   mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
// // //   mlMsh.PrintInfo();
// // // 
// // //   // define the multilevel solution and attach the mlMsh object to it
// // //   MultiLevelSolution mlSol(&mlMsh);
// // // 
// // //   // add variables to mlSol
// // //   mlSol.AddSolution("U", LAGRANGE, FIRST);
// // //   mlSol.AddSolution("V", LAGRANGE, SERENDIPITY);
// // //   mlSol.AddSolution("W", LAGRANGE, SECOND);
// // //   mlSol.AddSolution("P", DISCONTINOUS_POLYNOMIAL, ZERO);
// // //   mlSol.AddSolution("T", DISCONTINOUS_POLYNOMIAL, FIRST);
// // // 
// // //   mlSol.Initialize("All");    // initialize all varaibles to zero
// // // 
// // //   mlSol.Initialize("U", InitalValueU);
// // //   mlSol.Initialize("P", InitalValueP);
// // //   mlSol.Initialize("T", InitalValueT);    // note that this initialization is the same as piecewise constant element
// // // 
// // //   // print solutions
// // //   std::vector < std::string > variablesToBePrinted;
// // //   variablesToBePrinted.push_back("U");
// // //   variablesToBePrinted.push_back("P");
// // //   variablesToBePrinted.push_back("T");
// // // 
// // //   VTKWriter vtkIO(&mlSol);
// // //   vtkIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);
// // // 
// // //   GMVWriter gmvIO(&mlSol);
// // //   variablesToBePrinted.push_back("all");
// // //   gmvIO.SetDebugOutput(false);
// // //   gmvIO.write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);
// // // 
// // //   return 0;
// // // }




