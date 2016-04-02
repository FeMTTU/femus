
#include "FemusInit.hpp"
#include "Marker.hpp"
#include "MultiLevelMesh.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"


using namespace femus;

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  mlMsh.ReadCoarseMesh( "./input/square.neu", "seventh", scalingFactor );
  

  std::vector < double > x(3);
  
//   x[0]=1.61638e-05;  //should be on element 49 of square2.neu
//   x[1]=7.81317e-05;
//   x[2]=0.;

  x[0]=0.249999999;  //should be on element 27 of square.neu
  x[1]=0.21875;
  x[2]=0.;
  
  Marker a( x, VOLUME, mlMsh.GetLevel(0) );
  
  
  std::vector < double > y = a.GetMarkerCoordinates();
  
  std::cout<< " The coordinates of the marker are " << y[0] << " ," << y[1] << " ," <<y[2]<<std::endl;
  std::cout << " The marker type is " <<  a.GetMarkerType() <<std::endl;
  
  
  // print mesh
  
  MultiLevelSolution mlSol( &mlMsh );
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back( "All" );

  VTKWriter vtkIO( &mlSol );
  vtkIO.SetDebugOutput( true );
  vtkIO.Write( DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted );
  
  return 0;
}




