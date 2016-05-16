
#include "FemusInit.hpp"
#include "Marker.hpp"
#include "MultiLevelMesh.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"


using namespace femus;



bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level) {

  bool refine = 0;

  if (elemgroupnumber == 6 && level < 4) refine = 1;
  if (elemgroupnumber == 7 && level < 5) refine = 1;
  if (elemgroupnumber == 8 && level < 6) refine = 1;

  return refine;

}


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
 mlMsh.ReadCoarseMesh( "./input/square2.neu", "seventh", scalingFactor );
  
//   mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
  unsigned numberOfUniformLevels = 1;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);
  
  std::vector < double > x(3);

/////////////////////////////////////////////////////////////////////////////////////////////

// WARNING changes to the tollerances have been made after these Tests have been performed so the comments may not reflect
// exactly what is really happening
// /* TESTS ON THE MESH square2.neu*/  
  
//   // Test 1: the maker is VERY close to a vertex of element 49 (point 37)
//   x[0]=1.61638e-05;  //it is on element 49 
//   x[1]=7.81317e-05;
//   x[2]=0.;
//   // WARNING these are actually NOT the exact coordinates of the vertex (the first 5 decimal digits are equal).
//   // NOTE I put a tollerance of 1e-04 for equality and it seems to work nice.
    // I changed the tollerance to 1-e10 because with 1e-04 it was not working with square.neu however now the code, in this limit case, gives the
    // right answer but because it considers the edge PASSING THROUGH the marker and not anymore that the vertex is ON the marker
  
/////////////////////////////////////////////////////////////////////////////////////////////  

//   // Test 2: the maker is VERY close to a midpoint of an edge of element 60 (point 144)
//   x[0]=5.77398e-05;  //it is on element 60 
//   x[1]=7.23015e-05;
//   x[2]=0.;
//   // WARNING these are actually NOT the exact coordinates of the vertex (the first 5 decimal digits are equal).
//   // NOTE I put a tollerance of 1e-04 for equality and it seems to work nice.
     // I changed the tollerance to 1-e10 because with 1e-04 it was not working with square.neu however now the code, in this limit case, gives the
     // right answer but because it considers the edge PASSING THROUGH the marker and not anymore that the vertex is ON the marker
  
//////////////////////////////////////////////////////////////////////////////////////////////

// Test 3: the maker has the "same" y coordinate of a midpoint of an edge of element 60 (point 144)
//  x[0]=5.77298e-05;  //it is on element 60 
//  x[1]=7.23015e-05;
//  x[2]=0.;
  
//   x[0]=-4.9e-05;  //it is on element 153 proc=2
//   x[1]=2e-05;
//   x[2]=0.;
  
// WARNING these are actually NOT the exact coordinates of the vertex (the first 5 decimal digits are equal).
// NOTE I put a tollerance of 1e-04 for equality and it seems to work nice.
// I changed the tollerance to 1-e10 because with 1e-04 it was not working with square.neu however now the code, in this limit case, gives the
// right answer but it only considers the edge to intersect the x axis, not at the midpoint.
  
//////////////////////////////////////////////////////////////////////////////////////////////  
  
//Test 4: the maker is VERY close to a vertex of element 15 (point 19)
//  x[0]=-2.72848e-05;  //it is on element 15 
//  x[1]=2.62822e-05;
//  x[2]=0.;
//WARNING these are actually NOT the exact coordinates of the vertex (the first 5 decimal digits are equal).
//NOTE I put a tollerance of 1e-04 for equality and it seems to work nice.
 // I changed the tollerance to 1-e10 because with 1e-04 it was not working with square.neu however now the code, in this limit case, gives the
 // right answer but because it considers the edge PASSING THROUGH the marker and not anymore that the vertex is ON the marker
  
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  */
// WARNING changes to the tollerances have been made after these Tests have been performed so the comments may not reflect
// exactly what is really happening

/*  TESTS FOR THE MESH square.neu */  
  
//Test 1: the marker is just a little bit BELOW an horizontal edge (the lower edge of element 1 in mesh square.neu)  
/* x[0]=0.33375;  //the marker is outside of the mesh (of square.neu)
 x[1]=-0.062500000001; //if the y was -0.625 then it is on the lower edge of element 1
 x[2]=0.;  // */ 
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
// x[0]=0.3124999999999;  //the marker is on element 0
// x[1]=-0.05555555; // if the x is 0.3125 then it is on the vertical edge of element 1
// x[2]=0.; //
// WARNING With more than 8 nines the code considers the marker to be ON the edge and so it is considered to be on element 0
// only because that edge appears first in element 0 than in element 1.
// In this case it is ok because the element numbers increase from left to right.
// NOTE With 8 nines the code does what it is supposed to do (the error is 1e-12).

///////////////////////////////////////////////////////////////////////////////////////////

//Test 4: the marker is just a little bit to the RIGHT of a vertical edge (the left edge of element 1 in mesh square.neu) 
// x[0]=0.312500000001;  //the marker is on element 1
// x[1]=-0.05555555; // if the x is 0.3125 then it is on the vertical edge of element 1
// x[2]=0.; 
// WARNING With more than 7 zeros the marker is considered to be ON the vertical edge, which belongs first to element 0
// this is why it is considered to be on element 0 instead of on element 1.
// NOTE With 7 zeros, the code does what it is supposed to do, the error is (1e-12).

//Test 5: the marker is on the lower half of the RIGHT edge of element 0 
// x[0]=0.3125;  //
// x[1]=-0.05555555; // 
// x[2]=0.; 
// NOTE The code does what it is supposed to do.

////////////////////////////////////////////////////////////////////////////////////////

//Test 6: the marker is the north-east vertex of element 43
// x[0]=0.25;  //
// x[1]=0.3125; // 
// x[2]=0.; 
// NOTE The code does what it is supposed to do.

////////////////////////////////////////////////////////////////////////////////////////

// //Test 7: the marker is on element 35 but the y coordinate is the same as the MIDPOINT of the RIGHT edge of element 35 
// x[0]=0.23;  //
// x[1]=0.28125; // 
// x[2]=0.; 
// // NOTE The code does what it is supposed to do.


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

// TESTS ON THE MESH cubeHex.neu (tests run with 2 procs)

//Test 1
// x[0]=-0.3;  //  point 1409 shared by element 97,112 (proc 2) and 45 and 55 (proc 1) It says the marker is on element 45, yes!
// x[1]=0.1; // 
// x[2]=0.5; 

////////////////////////////////////////////////////////////////////////////////////////

// //Test 2: the marker is on the midpoint of an edge of element 60 (the trick was that the intersection point was actually ON a midpoint of an edge) 
// x[0]=-0.5;  //
// x[1]=0.5; // 
// x[2]=0.; 
 
 ////////////////////////////////////////////////////////////////////////////////////////
 
//  //Test 3: the marker is to the right of the midpoint of an edge of element 60 
//  x[0]=-0.5;  //
//  x[1]=0.5; // 
//  x[2]=0.005; 

  ////////////////////////////////////////////////////////////////////////////////////////
 
  // //Test 4: the marker is a vertex of element 60 
//  x[0]=-0.5;  //
//  x[1]=0.5; // 
//  x[2]=-0.1; //
 
  ////////////////////////////////////////////////////////////////////////////////////////
 
   // //Test 5: the marker is close to a vertex of element 60 but it belongs to element 62
/* x[0]=-0.5;  //
 x[1]=0.5; // 
 x[2]=-0.10000001; */// If you add one more zero then it will be in element 60.
 
  ////////////////////////////////////////////////////////////////////////////////////////
  
     // //Test 6: the marker is a little below a vertex of element 60 (it goes outside the domain)
//  x[0]=-0.5000001;  //
//  x[1]=0.5; // 
//  x[2]=-0.05; //

  ////////////////////////////////////////////////////////////////////////////////////////
    
    // //Test 7: the marker is over the lower edge of element 57 (it is not on the domain because the x is too small)
//  x[0]=-0.6;  //
//  x[1]=0.1; // 
//  x[2]=0.2; //


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

// TESTS ON THE MESH prism3D.neu (tests run with 2 procs)

//Test 1 :  element 20   OK
// x[0] = -0.5;
// x[1] = 0.;
// x[2] = 0.;

/////////////////////////////////////////////////////////////////////////////////////////

//Test 2 :  22
// x[0] = 0.;
// x[1] = 0.;
// x[2] = 0.001;
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */



  Marker a( x, VOLUME, mlMsh.GetLevel(0), true );
  //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
  
  std::vector <double> xc(3);
  xc[0]=0; // x coordinate of vertex 1
  xc[1]=0; // y coordinate of vertex 1 (in the reference frame)
  xc[2]=0; // z coordinate of vertex 1 (in the reference frame)
  a.InverseMappingTEST(xc);  

  
  
  
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




