// tests the inclusion algorithm for all 2D and 3D elements


#include "FemusInit.hpp"
#include "Marker.hpp"
#include "MultiLevelMesh.hpp"
#include "MultiLevelSolution.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "MyVector.hpp"

using namespace femus;


// 2D CASE rigid rotation
// double InitalValueU(const std::vector < double >& x) {
//   return -x[1];
// }
//
// double InitalValueV(const std::vector < double >& x) {
//   return x[0];
// }
//
// double InitalValueW(const std::vector < double >& x) {
//   return 0.;
// }


//3D CASE  rotation
// double InitalValueU(const std::vector < double >& x) {
//   return (-x[1]+x[2])/sqrt(3);
// }
//
// double InitalValueV(const std::vector < double >& x) {
//   return (x[0]-x[2])/sqrt(3);
// }
//
// double InitalValueW(const std::vector < double >& x) {
//   return (x[1]-x[0])/sqrt(3);
// }


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



// 3D CASE with vorticity
// double pi = acos(-1.);
//
// double InitalValueU(const std::vector < double >& x) {
//   double time = (x.size() == 4) ? x[3] : 0.;
//   return
//     2.*(sin(pi * (x[0] + 0.5)) * sin(pi * (x[0] + 0.5)) *
//     ( sin(pi * (x[1] + 0.5)) * cos(pi * (x[1] + 0.5)) - sin(pi * (x[2] + 0.5)) * cos(pi * (x[2] + 0.5)) )
//     )* cos(time);
// }
//
// double InitalValueV(const std::vector < double >& x) {
//   double time = (x.size() == 4) ? x[3] : 0.;
//   return
//     2.*(sin(pi * (x[1] + 0.5)) * sin(pi * (x[1] + 0.5)) *
//     ( sin(pi * (x[2] + 0.5)) * cos(pi * (x[2] + 0.5)) - sin(pi * (x[0] + 0.5)) * cos(pi * (x[0] + 0.5)) )
//     )* cos(time);
// }
//
// double InitalValueW(const std::vector < double >& x) {
//   double time = (x.size() == 4) ? x[3] : 0.;
//   return
//     2.*( sin(pi * (x[2] + 0.5)) * sin(pi * (x[2] + 0.5)) *
//     ( sin(pi * (x[0] + 0.5)) * cos(pi * (x[0] + 0.5)) - sin(pi * (x[1] + 0.5)) * cos(pi * (x[1] + 0.5)) )
//     )* cos(time);
//
//   return 0.;
// }



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



  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("V", LAGRANGE, SECOND, 2);
  if(dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND, 2);
  mlSol.Initialize("U" , InitalValueU);
  mlSol.Initialize("V" , InitalValueV);
  if(dim == 3) mlSol.Initialize("W", InitalValueW);

//   //Test 1 (QUAD):

//   x[0] = 0.125;
//   x[1] = 0.125;
//   x[2] = -0.25;


  std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
// Marker a1Quad(x, VOLUME, mlMsh.GetLevel(0), solType, true);
  //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
  //std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
  //std::cout << " The marker type is " <<  a1Quad.GetMarkerType() << std::endl;

  double T = 2 * acos(-1.);
  unsigned n  = 1;


//   std::vector< std::vector < std::vector < double > > > xn(1);
//   xn[0].resize(n + 1);
//   for(unsigned k = 0; k < n; k++) {
//     a1Quad.GetMarkerCoordinates(xn[0][k]);
//     a1Quad.Advection(mlSol.GetLevel(0), 2, T / n);
//   }
//   a1Quad.GetMarkerCoordinates(xn[0][n]);
//   for(unsigned i = 0;  i < xn[0].size(); i++) {
//     for(unsigned d = 0; d < xn[0][i].size(); d++) {
//       //   std::cout << xn[0][i][d] << " ";
//     }
//     // std::cout << std::endl;
//   }
//
//
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

// PrintLine(DEFAULT_OUTPUTDIR, xn);


  clock_t start_time = clock();
  clock_t init_time = clock();


  unsigned pSize = 100;
  std::vector < Marker*> particle(pSize);


//uncomment to do solid body rotation and vortex test
  double pi = acos(-1.);
  for(unsigned j = 0; j < pSize; j++) {
    x[0] = 0. + 0.125 * cos(2.*pi / pSize * j);
    x[1] = .25 + 0.125 * sin(2.*pi / pSize * j);
    x[2] = 0.;
    particle[j] = new Marker(x, 0., VOLUME, mlSol.GetLevel(numberOfUniformLevels - 1), solType, true);
  }


  std::vector < std::vector < std::vector < double > > > line(1);
  line[0].resize(pSize + 1);

  for(unsigned j = 0; j < pSize; j++) {
    particle[j]->GetMarkerCoordinates(line[0][j]);
  }
  particle[0]->GetMarkerCoordinates(line[0][pSize]);

  std::vector < std::vector < std::vector < double > > > line0 = line; // saves the initial position
  PrintLine(DEFAULT_OUTPUTDIR, line, false, 0);

  n = 30;

  std::cout << std::endl << " init in  " << std::setw(11) << std::setprecision(6) << std::fixed
            << static_cast<double>((clock() - init_time)) / CLOCKS_PER_SEC << " s" << std::endl;


  clock_t advection_time = clock();


  for(unsigned k = 1; k <= n; k++) {
    //uncomment for  vortex test
    mlSol.CopySolutionToOldSolution();
    mlSol.UpdateSolution("U" , InitalValueU, pi * k / n);
    mlSol.UpdateSolution("V" , InitalValueV, pi * k / n);
    if(dim == 3) mlSol.UpdateSolution("W" , InitalValueW, pi * k / n);


    //uncomment for vortex test and rigid rotation
    for(unsigned j = 0; j < pSize; j++) {
//       std::cout << j <<" " << k << std::endl<<std::flush;
      particle[j]->Advection( n, T / n,  mlSol.GetLevel(numberOfUniformLevels - 1));
      particle[j]->GetMarkerCoordinates(line[0][j]);
    }
    particle[0]->GetMarkerCoordinates(line[0][pSize]);
    PrintLine(DEFAULT_OUTPUTDIR, line, false, k);
  }


  std::cout << std::endl << " advection in: " << std::setw(11) << std::setprecision(6) << std::fixed
            << static_cast<double>((clock() - advection_time)) / CLOCKS_PER_SEC << " s" << std::endl;

  std::cout << std::endl << " RANNA in: " << std::setw(11) << std::setprecision(6) << std::fixed
            << static_cast<double>((clock() - start_time)) / CLOCKS_PER_SEC << " s" << std::endl;


  //computing the geometric error
  double error = 0.;
  for(unsigned j = 0; j < pSize + 1; j++) {
    double tempError = 0.;
    for(unsigned i = 0; i < dim; i++) {
      tempError += (line0[0][j][i] - line[0][j][i]) * (line0[0][j][i] - line[0][j][i]);
    }
    error += sqrt(tempError);
  }

  error = error / pSize;

   std::cout << " ERROR = " << std::setprecision(15) << error << std::endl;

 

  for(unsigned j = 0; j < pSize; j++) {
    delete particle[j];
  }



  return 0;
}



