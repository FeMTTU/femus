// tests the inclusion algorithm for all 2D and 3D elements


#include "FemusInit.hpp"
#include "Marker.hpp"
#include "Line.hpp"
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
// double InitalValueU(const std::vector < double >& x)
// {
//   return (-x[1] + x[2]) / sqrt(3);
// }
// 
// double InitalValueV(const std::vector < double >& x)
// {
//   return (x[0] - x[2]) / sqrt(3);
// }
// 
// double InitalValueW(const std::vector < double >& x)
// {
//   return (x[1] - x[0]) / sqrt(3);
// }


// //2D CASE with vorticity
// double pi = acos(-1.);
// 
// double InitalValueU(const std::vector < double >& x) {
//   double time = (x.size() == 4) ? x[3] : 0.;
//   return 2. * sin(pi * (x[0] + 0.5)) * sin(pi * (x[0] + 0.5)) * sin(pi * (x[1] + 0.5)) * cos(pi * (x[1] + 0.5)) * cos(time);
// }
// 
// double InitalValueV(const std::vector < double >& x) {
//   double time = (x.size() == 4) ? x[3] : 0.;
//   return -2. * sin(pi * (x[1] + 0.5)) * sin(pi * (x[1] + 0.5)) * sin(pi * (x[0] + 0.5)) * cos(pi * (x[0] + 0.5)) * cos(time);
// }
// 
// double InitalValueW(const std::vector < double >& x) {
//   double time = (x.size() == 4) ? x[3] : 0.;
//   return 0.;
// }



// 3D CASE with vorticity
double pi = acos(-1.);

double InitalValueU(const std::vector < double >& x) {
  double time = (x.size() == 4) ? x[3] : 0.;
  return
    2.*(sin(pi * (x[0] + 0.5)) * sin(pi * (x[0] + 0.5)) *
    ( sin(pi * (x[1] + 0.5)) * cos(pi * (x[1] + 0.5)) - sin(pi * (x[2] + 0.5)) * cos(pi * (x[2] + 0.5)) )
    )* cos(time);
}

double InitalValueV(const std::vector < double >& x) {
  double time = (x.size() == 4) ? x[3] : 0.;
  return
    2.*(sin(pi * (x[1] + 0.5)) * sin(pi * (x[1] + 0.5)) *
    ( sin(pi * (x[2] + 0.5)) * cos(pi * (x[2] + 0.5)) - sin(pi * (x[0] + 0.5)) * cos(pi * (x[0] + 0.5)) )
    )* cos(time);
}

double InitalValueW(const std::vector < double >& x) {
  double time = (x.size() == 4) ? x[3] : 0.;
  return
    2.*( sin(pi * (x[2] + 0.5)) * sin(pi * (x[2] + 0.5)) *
    ( sin(pi * (x[0] + 0.5)) * cos(pi * (x[0] + 0.5)) - sin(pi * (x[1] + 0.5)) * cos(pi * (x[1] + 0.5)) )
    )* cos(time);

  return 0.;
}


bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level)
{

  bool refine = 0;

  if (elemgroupnumber == 6 && level < 4) refine = 1;
  if (elemgroupnumber == 7 && level < 5) refine = 1;
  if (elemgroupnumber == 8 && level < 6) refine = 1;

  return refine;

}


int main(int argc, char** args)
{

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  unsigned numberOfUniformLevels = 3; //for refinement in 3D
  //unsigned numberOfUniformLevels = 1;
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
  mlMsh.ReadCoarseMesh("./input/test3Dbis.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/test2Dbis.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/test2D.neu", "seventh", scalingFactor);
  //mlMsh.ReadCoarseMesh("./input/cubeTet.neu", "seventh", scalingFactor);
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

  unsigned dim = mlMsh.GetDimension();

//   unsigned size = 100;
//   std::vector < std::vector < double > > x; // marker
//   std::vector < MarkerType > markerType;
//
//   x.resize(size);
//   markerType.resize(size);
//
//   for(unsigned j = 0; j < size; j++) {
//     x[j].assign(dim, 0);
//     markerType[j] = VOLUME;
//   }
//
//   double pi = acos(-1.);
//   for(unsigned j = 0; j < size; j++) {
//     x[j][0] = 0. + 0.125 * cos(2.*pi / size * j);
//     x[j][1] = .25 + 0.125 * sin(2.*pi / size * j);
//    // x[j][2] = 0.;
//   }
//
//   Line linea(x, markerType, mlMsh.GetLevel(numberOfUniformLevels - 1), solType);


  //BEGIN TEST FOR UPDATE LINE: to use it, make _particles public and initialize as the identity printList in UpdateLine
//   std::vector <Marker* > particles(size);
//
//   for(unsigned j = 0; j < size; j++) {
//     std::vector <double> y(dim);
//     y[0] = -0.15 + 0.2 * cos(2.*pi / size * j);
//     y[1] = 0. + 0.2 * sin(2.*pi / size * j);
//     particles[j] = new Marker(y, VOLUME, mlMsh.GetLevel(numberOfUniformLevels - 1), solType, true);
//   }
//
//   for(unsigned i=0; i<size; i++){
//     linea._particles[i] = particles[i];
//   }
//
//   linea.UpdateLine();

  //END

//   PrintLine(DEFAULT_OUTPUTDIR, linea.GetLine(), false, 0);

  MultiLevelSolution mlSol(&mlMsh);
  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND, 2);
  mlSol.AddSolution("V", LAGRANGE, SECOND, 2);
  if (dim == 3) mlSol.AddSolution("W", LAGRANGE, SECOND, 2);
  mlSol.Initialize("U" , InitalValueU);
  mlSol.Initialize("V" , InitalValueV);
  if (dim == 3) mlSol.Initialize("W", InitalValueW);

  std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
// Marker a1Quad(x, VOLUME, mlMsh.GetLevel(0), solType, true);
  //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
  //std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
  //std::cout << " The marker type is " <<  a1Quad.GetMarkerType() << std::endl;

  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted);

  clock_t start_time = clock();
  clock_t init_time = clock();


  unsigned size = 10000;

  std::vector < std::vector < double > > x; // marker
  std::vector < MarkerType > markerType;

  x.resize(size);
  markerType.resize(size);

  std::vector < std::vector < std::vector < double > > > line(1);
  std::vector < std::vector < std::vector < double > > > line0(1);

  for (unsigned j = 0; j < size; j++) {
    x[j].assign(dim, 0.);
    markerType[j] = VOLUME;
  }

   srand(2); //TODO 3D rotation n=10, problem at iteration 6 with seed srand(1);    FIXED
             //TODO 3D vortex n=16, problem at iteration 11 with seed srand(2);     FIXED
             //TODO 3D vortex srand(2) gives different errors
  double pi = acos(-1.);
  for (unsigned j = 0; j < size; j++) {

    //BEGIN random initialization
//     double r_rad = static_cast <double> (rand()) / RAND_MAX;
//     r_rad = 0.4 * (1. - r_rad * r_rad * r_rad);
//     double r_theta = static_cast <double> (rand()) / RAND_MAX * 2 * pi;
//     double r_phi = static_cast <double> (rand()) / RAND_MAX * pi;
// 
//     
//     
//     x[j][0] = r_rad * sin(r_phi) * cos(r_theta);
//     x[j][1] = r_rad * sin(r_phi) * sin(r_theta);
//     if (dim == 3) {
//       x[j][2] = r_rad * cos(r_phi);
//     }
    //END
 
 //  if(j==1012) std::cout<< std::setprecision(14)<<x[j][0] <<" "<< x[j][1]<<" "<<x[j][2]<<std::endl;

//BEGIN ordered initialization
    x[j][0] = 0. + 0.125 * cos(2.*pi / size * j);
    x[j][1] = .25 + 0.125 * sin(2.*pi / size * j);
    if (dim == 3) {
      x[j][2] = 0.;
    }
//END 
    
    
  } //end for on initialization

     
 //exit(0);

  Line linea(x, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);

  linea.GetLine(line0[0]);
  PrintLine(DEFAULT_OUTPUTDIR, line0, false, 0);

  double T = 2 * acos(-1.);

  unsigned n = 4;

  std::cout << std::endl << " init in  " << std::setw(11) << std::setprecision(6) << std::fixed
            << static_cast<double>((clock() - init_time)) / CLOCKS_PER_SEC << " s" << std::endl;


  clock_t advection_time = clock();

  for (unsigned k = 1; k <= n; k++) {
    std::cout << "Iteration = " << k << std::endl;
    //uncomment for  vortex test
    mlSol.CopySolutionToOldSolution();
    mlSol.UpdateSolution("U" , InitalValueU, pi * k / n);
    mlSol.UpdateSolution("V" , InitalValueV, pi * k / n);
    if (dim == 3) mlSol.UpdateSolution("W" , InitalValueW, pi * k / n);
    linea.AdvectionParallel(40, T / n, 4);
    linea.GetLine(line[0]);
    PrintLine(DEFAULT_OUTPUTDIR, line, false, k);
  }


  std::cout << std::endl << " advection in: " << std::setw(11) << std::setprecision(6) << std::fixed
            << static_cast<double>((clock() - advection_time)) / CLOCKS_PER_SEC << " s" << std::endl;

  std::cout << std::endl << " RANNA in: " << std::setw(11) << std::setprecision(6) << std::fixed
            << static_cast<double>((clock() - start_time)) / CLOCKS_PER_SEC << " s" << std::endl;

  //computing the geometric error
  double error = 0.;
  for (unsigned j = 0; j < size + 1; j++) {
    double tempError = 0.;
    for (unsigned i = 0; i < dim; i++) {
      tempError += (line0[0][j][i] - line[0][j][i]) * (line0[0][j][i] - line[0][j][i]);
    }
    error += sqrt(tempError);
  }

  error = error / size;

  std::cout << " ERROR = " << std::setprecision(15) << error << std::endl;

//   for(unsigned j = 0; j < size; j++) {
//     std::vector <double> trial(dim);
//     trial = linea._particles[linea._printList[j]]->GetIprocMarkerCoordinates();
//     for(unsigned i=0; i<dim; i++){
//       std::cout << " x[" << j << "][" << i << "]=" << trial[i] << std::endl;
//     }
//   }



  return 0;
}



