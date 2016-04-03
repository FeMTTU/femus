
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

/*  TESTS FOR THE MESH square.neu */  
  
//Test 1: the marker is just a little bit BELOW an horizontal edge (the lower edge of element 1 in mesh square.neu)  
//  x[0]=0.33375;  //the marker is outside of the mesh (of square.neu)
//  x[1]=-0.062500000001; //if the y was -0.625 then it is on the lower edge of element 1
//  x[2]=0.;  //  
//WARNING With more than 7 zeros the marker is considered to be ON the edge and so it would belong to element 1 anyways.
// NOTE the code does what it is supposed to do, the error is 1e-12.

////////////////////////////////////////////////////////////////////////////////////////////

//Test 2: the marker is just a little bit ABOVE an horizontal edge (the lower edge of element 1 in mesh square.neu)  
//   x[0]=0.33375;  //the marker is inside the mesh (of square.neu), it is in element 1
//   x[1]=-0.062499999999; //if the y was -0.625 then it is on the lower edge of element 1, 
//   x[2]=0.;   //
//WARNING With more than 9 nines the code thinks the marker is actually ON the edge. It gives the right result but for the wrong reason.
//NOTE the code does what it is supposed to do , the error is (1e-12). 

//////////////////////////////////////////////////////////////////////////////////////////

//Test 3: the marker is just a little bit to the LEFT of a vertical edge (the left edge of element 1 in mesh square.neu) 
// x[0]=0.312499999999;  //the marker is on element 0
// x[1]=-0.05555555; // if the x is 0.3125 then it is on the vertical edge of element 1
// x[2]=0.; //
// WARNING With more than 8 nines the code considers the marker to be ON the edge and so it is considered to be on element 0
// only because that edge appears first in element 0 than in element 1.
// In this case it is ok because the element numbers increase from left to right.
// NOTE With 8 nines the code does what it is supposed to do (the error is 1e-12).

/////////////////////////////////////////////////////////////////////////////////////////

//Test 4: the marker is just a little bit to the RIGHT of a vertical edge (the left edge of element 1 in mesh square.neu) 
//x[0]=0.3125000000001;  //the marker is on element 1
//x[1]=-0.05555555; // if the x is 0.3125 then it is on the vertical edge of element 1
//x[2]=0.; 
//WARNING With more than 7 zeros the marker is considered to be ON the vertical edge, which belongs first to element 0
// this is why it is considered to be on element 0 instead of on element 1.
// NOTE With 7 zeros, the code does what it is supposed to do, the error is (1e-12).

/////////////////////////////////////////////////////////////////////////////////////////

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




