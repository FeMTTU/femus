
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
  mlMsh.ReadCoarseMesh( "./input/square2.neu", "seventh", scalingFactor );
  

  std::vector < double > x(3);
  
  x[0]=-5e-06; 
  x[1]=9.05928e-05;
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




