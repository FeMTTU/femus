// tests the inclusion algorithm for all 2D and 3D elements


#include "FemusInit.hpp"
#include "Marker.hpp"
#include "MultiLevelMesh.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"


using namespace femus;

double InitalValueU(const std::vector < double >& x) {
  return 0.5;
}

double InitalValueV(const std::vector < double >& x) {
  return 0.5;
}

double InitalValueW(const std::vector < double >& x) {
  return 0.5;
}

bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level) {

  bool refine = 0;

  if(elemgroupnumber == 6 && level < 4) refine = 1;
  if(elemgroupnumber == 7 && level < 5) refine = 1;
  if(elemgroupnumber == 8 && level < 6) refine = 1;

  return refine;

}


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  std::vector < double > x(3, 0); // marker
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  std::vector < std::string > variablesToBePrinted;

  /* element types
  0 = HEX
  1 = TET
  2 = WEDGE
  3 = QUAD
  4 = TRI
   */

  unsigned solType = 2;

  std::cout << " --------------------------------------------------     TEST     --------------------------------------------------" << std::endl;

  mlMsh.ReadCoarseMesh("./input/prism3D.neu", "seventh", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

  unsigned dim = mlMsh.GetDimension();


 
  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND);
  mlSol.AddSolution("V", LAGRANGE, SECOND);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND);
  mlSol.Initialize("U" , InitalValueU);
  mlSol.Initialize("V" , InitalValueV);
  if(dim == 3) mlSol.Initialize("W", InitalValueW);

//   //Test 1 (QUAD):
//   
//     NOTE tests ran with 4 procs
//   x[0] = -0.46875; //the marker is in element 191 (proc 3 of 4)
//   x[1] = -0.5;
//   x[2] = 0.;
  

//Test 1 (TET):  element 20
    //NOTE Tests ran with 2 procs
      x[0] = -0.5;
      x[1] = 0.;
      x[2] = 0.;
  
  std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
  Marker a1Quad(x, VOLUME, mlMsh.GetLevel(0), solType, true);
  //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
  std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
  std::cout << " The marker type is " <<  a1Quad.GetMarkerType() << std::endl;
  
  double T = 1.0;
  unsigned n  = 20;
  a1Quad.Advection(mlSol.GetLevel(0), n, T);

  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);





  return 0;
}




