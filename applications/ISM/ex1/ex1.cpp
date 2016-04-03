
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
  
// /* TESTS ON THE MESH square2.neu*/  
  
//   // Test 1: the maker is VERY close to a vertex of element 49 (point 37)
//   x[0]=1.61638e-05;  //it is on element 49 
//   x[1]=7.81317e-05;
//   x[2]=0.;
//   // WARNING these are actually NOT the exact coordinates of the vertex (the first 5 decimal digits are equal).
//   // NOTE I put a tollerance of 1e-04 for equality and it seems to work nice.
  
/////////////////////////////////////////////////////////////////////////////////////////////  

//   // Test 2: the maker is VERY close to a midpoint of an edge of element 60 (point 144)
//   x[0]=5.77398e-05;  //it is on element 60 
//   x[1]=7.23015e-05;
//   x[2]=0.;
//   // WARNING these are actually NOT the exact coordinates of the vertex (the first 5 decimal digits are equal).
//   // NOTE I put a tollerance of 1e-04 for equality and it seems to work nice.
  
//////////////////////////////////////////////////////////////////////////////////////////////

// Test 3: the maker has the "same" x coordinate of a midpoint of an edge of element 60 (point 144)
//  x[0]=5.77298e-05;  //it is on element 60 
//  x[1]=7.23015e-05;
//  x[2]=0.;
// WARNING these are actually NOT the exact coordinates of the vertex (the first 5 decimal digits are equal).
// NOTE I put a tollerance of 1e-04 for equality and it seems to work nice.
  
//Test 4: the maker is VERY close to a vertex of element 15 (point 19)
 x[0]=-2.72848e-05;  //it is on element 15 
 x[1]=2.62822e-05;
 x[2]=0.;
//WARNING these are actually NOT the exact coordinates of the vertex (the first 5 decimal digits are equal).
//NOTE I put a tollerance of 1e-04 for equality and it seems to work nice.
  
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */

/*  TESTS FOR THE MESH square.neu */  
  
//Test 1: the marker is just a little bit BELOW an horizontal edge (the lower edge of element 1 in mesh square.neu)  
//  x[0]=0.33375;  //the marker is outside of the mesh (of square.neu)
//  x[1]=-0.062500000001; //if the y was -0.625 then it is on the lower edge of element 1
//  x[2]=0.;  //  
// WARNING With more than 7 zeros the marker is considered to be ON the edge and so it would belong to element 1 anyways.
// NOTE the code does what it is supposed to do, the error is 1e-12.

////////////////////////////////////////////////////////////////////////////////////////////

//Test 2: the marker is just a little bit ABOVE an horizontal edge (the lower edge of element 1 in mesh square.neu)  
//   x[0]=0.33375;  //the marker is inside the mesh (of square.neu), it is in element 1
//   x[1]=-0.062499999999; //if the y was -0.625 then it is on the lower edge of element 1, 
//   x[2]=0.;   //
// WARNING With more than 9 nines the code thinks the marker is actually ON the edge. It gives the right result but for the wrong reason.
// NOTE the code does what it is supposed to do , the error is (1e-12). 

////////////////////////////////////////////////////////////////////////////////////////////

//Test 3: the marker is just a little bit to the LEFT of a vertical edge (the left edge of element 1 in mesh square.neu) 
// x[0]=0.312499999999;  //the marker is on element 0
// x[1]=-0.05555555; // if the x is 0.3125 then it is on the vertical edge of element 1
// x[2]=0.; //
// WARNING With more than 8 nines the code considers the marker to be ON the edge and so it is considered to be on element 0
// only because that edge appears first in element 0 than in element 1.
// In this case it is ok because the element numbers increase from left to right.
// NOTE With 8 nines the code does what it is supposed to do (the error is 1e-12).

///////////////////////////////////////////////////////////////////////////////////////////

//Test 4: the marker is just a little bit to the RIGHT of a vertical edge (the left edge of element 1 in mesh square.neu) 
//x[0]=0.3125000000001;  //the marker is on element 1
//x[1]=-0.05555555; // if the x is 0.3125 then it is on the vertical edge of element 1
//x[2]=0.; 
// WARNING With more than 7 zeros the marker is considered to be ON the vertical edge, which belongs first to element 0
// this is why it is considered to be on element 0 instead of on element 1.
// NOTE With 7 zeros, the code does what it is supposed to do, the error is (1e-12).

//Test 5: the marker is on the lower half of the RIGHT edge of element 0 
//x[0]=0.3125;  //
//x[1]=-0.05555555; // 
//x[2]=0.; 
// NOTE The code does what it is supposed to do.

////////////////////////////////////////////////////////////////////////////////////////

//Test 6: the marker is the north-east vertex of element 35
//x[0]=0.25;  //
//x[1]=0.3125; // 
//x[2]=0.; 
// NOTE The code does what it is supposed to do.

////////////////////////////////////////////////////////////////////////////////////////

//Test 7: the marker is on element 35 but the y coordinate is the same as the MIDPOINT of the RIGHT edge of element 35 
//x[0]=0.23;  //
//x[1]=0.28125; // 
//x[2]=0.; 
// NOTE The code does what it is supposed to do.


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

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




