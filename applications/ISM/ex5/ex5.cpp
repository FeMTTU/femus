// tests the inclusion algorithm for all 2D and 3D elements


#include "FemusInit.hpp"
#include "Marker.hpp"
#include "Line.hpp"
#include "MultiLevelMesh.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "MyVector.hpp"

using namespace femus;


// 2D CASE with vorticity
double pi = acos(-1.);

double InitalValueU(const std::vector < double >& x) {
  double time = (x.size() == 4) ? x[3] : 0.;
  return 2. * sin(pi * (x[0] + 0.5)) * sin(pi * (x[0] + 0.5)) * sin(pi * (x[1] + 0.5)) * cos(pi * (x[1] + 0.5)) * cos(time);
}

double InitalValueV(const std::vector < double >& x) {
  double time = (x.size() == 4) ? x[3] : 0.;
  return -2. * sin(pi * (x[1] + 0.5)) * sin(pi * (x[1] + 0.5)) * sin(pi * (x[0] + 0.5)) * cos(pi * (x[0] + 0.5)) * cos(time);
}

double InitalValueW(const std::vector < double >& x) {
  double time = (x.size() == 4) ? x[3] : 0.;
  return 0.;
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

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  //unsigned numberOfUniformLevels = 3; //for refinement in 3D
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

  //mlMsh.ReadCoarseMesh("./input/prism3D.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/tri2.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cubeMixed.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/test3Dbis.neu", "seventh", scalingFactor);
  mlMsh.ReadCoarseMesh("./input/test2Dbis.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/test2D.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cubeTet.neu", "seventh", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

  unsigned dim = mlMsh.GetDimension();


  unsigned size = 100;
  std::vector < std::vector < double > > x; // marker
  std::vector < MarkerType > markerType;
  
  x.resize(size);
  markerType.resize(size);
  
  for(unsigned j = 0; j < size; j++) {
    x[j].assign(dim, 0);
    markerType[j] = VOLUME;
  }

  double pi = acos(-1.);
  for(unsigned j = 0; j < size; j++) {
    x[j][0] = 0. + 0.125 * cos(2.*pi / size * j);
    x[j][1] = .25 + 0.125 * sin(2.*pi / size * j);
    x[j][2] = 0.;
  }
  
  Line linea(x, markerType, mlMsh.GetLevel(numberOfUniformLevels - 1), solType, size, true);


  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("V", LAGRANGE, SECOND, 2);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND, 2);
  mlSol.Initialize("U" , InitalValueU);
  mlSol.Initialize("V" , InitalValueV);
  if(dim == 3) mlSol.Initialize("W", InitalValueW);

  std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
// Marker a1Quad(x, VOLUME, mlMsh.GetLevel(0), solType, true);
  //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
  //std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
  //std::cout << " The marker type is " <<  a1Quad.GetMarkerType() << std::endl;

  double T = 2 * acos(-1.);
  unsigned n  = 1;


  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  PrintLine(DEFAULT_OUTPUTDIR, linea.GetLine(), false, 0);

  return 0;
}



