// tests the inclusion algorithm for all 2D and 3D elements


#include "FemusInit.hpp"
#include "Marker.hpp"
#include "MultiLevelMesh.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"


using namespace femus;



bool SetRefinementFlag(const std::vector < double >& x, const int& elemgroupnumber, const int& level) {

  bool refine = 0;

  if(elemgroupnumber == 6 && level < 4) refine = 1;
  if(elemgroupnumber == 7 && level < 5) refine = 1;
  if(elemgroupnumber == 8 && level < 6) refine = 1;

  return refine;

}


int main(int argc, char** args) {

//   // init Petsc-MPI communicator
//   FemusInit mpinit(argc, args, MPI_COMM_WORLD);
// 
//   std::vector < double > x(3, 0); // marker
//   MultiLevelMesh mlMsh;
//   double scalingFactor = 1.;
//   unsigned numberOfUniformLevels = 1;
//   unsigned numberOfSelectiveLevels = 0;
//   std::vector < std::string > variablesToBePrinted;
// 
//   /* element types
//   0 = HEX
//   1 = TET
//   2 = WEDGE
//   3 = QUAD
//   4 = TRI
//    */
// 
//   int elementType = 2; // this decides what elements to test
// 
//   unsigned solType = 2; 
// 
//   switch(elementType) {
//     case 3: { // QUAD
// 
//       std::cout << " --------------------------------------------------     QUAD     --------------------------------------------------" << std::endl;
// 
//       mlMsh.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor);
// //   mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
//       mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);
// 
//       //NOTE tests ran with 4 procs
// //Test 1 (QUAD):
//       x[0] = 0.33375; //the marker is in element 117 (proc 1)
//       x[1] = -0.0627;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a1QUAD(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a1QUAD.GetMarkerType() << std::endl;
// 
// //Test 2 (QUAD):
//       x[0] = -0.0625; //the marker is in element 245 / 246 (proc 3)
//       x[1] = -0.09375;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a2QUAD(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a2QUAD.GetMarkerType() << std::endl;
// 
// //Test 3 (QUAD):
//       x[0] = 0; //the marker is shared by 4 elements, each one in a different proc
//       x[1] = 0.0625;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a3QUAD(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a3QUAD.GetMarkerType() << std::endl;
// 
// //Test 4 (QUAD):
//       x[0] = 0.4377; //the marker is on element 63, it is on the boundary of the domain
//       x[1] = 0.5;
//       x[2] = 0.;
// 
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a4QUAD(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a4QUAD.GetMarkerType() << std::endl;
// 
// //Test 5 (QUAD): the marker is on the lower half of the RIGHT edge of element 0
//       x[0] = 0.3125;
//       x[1] = -0.05555555;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a5QUAD(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a5QUAD.GetMarkerType() << std::endl;
// 
// //Test 6 (QUAD): the marker is the north-east vertex of element 43
//       x[0] = 0.25;
//       x[1] = 0.3125;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a6QUAD(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a6QUAD.GetMarkerType() << std::endl;
// 
// // Test 7 (QUAD): the marker is on element 35 and the y coordinate is the same as the MIDPOINT of the RIGHT edge of element 35
//       x[0] = 0.23;
//       x[1] = 0.28125;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a7QUAD(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a7QUAD.GetMarkerType() << std::endl;
// 
// // Test 8 (QUAD): the marker the central face point of element 27
//       x[0] = 0.21875;
//       x[1] = 0.21875;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a8QUAD(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a8QUAD.GetMarkerType() << std::endl;
// 
// // Test 9 (QUAD): the marker is outside the domain
// 
//         x[0]=1.;
//         x[1]=0.;
//         x[2]=0.;
// 
// 
//         std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//         Marker a9QUAD( x, VOLUME, mlMsh.GetLevel(0), solType, true );
//         //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//         std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
//         std::cout << " The marker type is " <<  a9QUAD.GetMarkerType() <<std::endl;
// 
//         //print mesh
//         MultiLevelSolution mlSol( &mlMsh );
//         variablesToBePrinted.push_back( "All" );
// 
//         VTKWriter vtkIO( &mlSol );
//         vtkIO.SetDebugOutput( true );
//         vtkIO.Write( Files::_application_output_directory, "biquadratic", variablesToBePrinted );
// 
// 
// 
//     }
//     break;
// 
// 
//     case 4: { // TRI
//       std::cout << " --------------------------------------------------     TRI      --------------------------------------------------" << std::endl;
// 
//       mlMsh.ReadCoarseMesh("./input/square2.neu", "seventh", scalingFactor);
// //   mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
//       mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);
// 
// //Test 1 (TRI):
//       x[0] = 5.e-05; // it is on element 12
//       x[1] = 0.00015;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a1TRI(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a1TRI.GetMarkerType() << std::endl;
//       
// //Test 2 (TRI):
//       x[0] = 8.29839e-05; //it is on element 115
//       x[1] = -2.38334e-05;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a2TRI(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a2TRI.GetMarkerType() << std::endl;
// 
// // Test 3 (TRI):
//       x[0] = -5.e-05; //it is on element 153
//       x[1] = 2.9e-05;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a3TRI(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a3TRI.GetMarkerType() << std::endl;
// 
// // Test 4 (TRI):
//       x[0] = -5.5e-05; //it is outside the domain
//       x[1] = 2.9e-05;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a4TRI(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a4TRI.GetMarkerType() << std::endl;
// 
// 
//       //print mesh
//       MultiLevelSolution mlSol(&mlMsh);
//       variablesToBePrinted.push_back("All");
// 
//       VTKWriter vtkIO(&mlSol);
//       vtkIO.SetDebugOutput(true);
//       vtkIO.Write(Files::_application_output_directory, "biquadratic", variablesToBePrinted);
// 
//     }
//     break;
// 
// 
//     case 0: // HEX
// 
//     {
// 
//       std::cout << " --------------------------------------------------     HEX      --------------------------------------------------" << std::endl;
// 
//      mlMsh.ReadCoarseMesh("./input/cubeHex.neu", "seventh", scalingFactor);
// //   mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
//       
//       mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);
// 
// //NOTE Tests ran with 2 procs
// //Test 1 (HEX):
//       x[0] = -0.3; //  point 98 shared by element 97,112 (proc 1) and 45 and 55 (proc 0). Proc 0 says the marker is in 45, proc 1 says it is is 112.
//       x[1] = 0.1; //
//       x[2] = 0.5;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a1HEX(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a1HEX.GetMarkerType() << std::endl;
// 
// //Test 2 (HEX): the marker is on the midpoint of an edge of element 60
//       x[0] = -0.5;
//       x[1] = 0.5;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a2HEX(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a2HEX.GetMarkerType() << std::endl;
// 
// //Test 3 (HEX): the marker is to the right of the midpoint of an edge of element 60
//       x[0] = -0.5;
//       x[1] = 0.5;
//       x[2] = 0.005;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a3HEX(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a3HEX.GetMarkerType() << std::endl;
// 
// //Test 4 (HEX): the marker is a vertex of element 60
//       x[0] = -0.5;
//       x[1] = 0.5;
//       x[2] = -0.1;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a4HEX(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a4HEX.GetMarkerType() << std::endl;
// 
//       //Test 5 (HEX): the marker is close to a vertex of element 60 but it belongs to element 62
//       x[0] = -0.5;
//       x[1] = 0.5;
//       x[2] = -0.10000001; // If you add one more zero then it will be in element 60.
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a5HEX(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a5HEX.GetMarkerType() << std::endl;
// 
// //Test 6 (HEX): the marker is a little below a vertex of element 60 (it goes outside the domain)
//       x[0] = -0.5000001;
//       x[1] = 0.5;
//       x[2] = -0.05;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a6HEX(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a6HEX.GetMarkerType() << std::endl;
// 
// //Test 7 (HEX): the marker is over the lower edge of element 57 (it is not on the domain because the x is too small)
//       x[0] = -0.6;
//       x[1] = 0.1;
//       x[2] = 0.2;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a7HEX(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a7HEX.GetMarkerType() << std::endl;
// 
// 
//       //print mesh
//       MultiLevelSolution mlSol(&mlMsh);
//       variablesToBePrinted.push_back("All");
// 
//       VTKWriter vtkIO(&mlSol);
//       vtkIO.SetDebugOutput(true);
//       vtkIO.Write(Files::_application_output_directory, "biquadratic", variablesToBePrinted);
// 
//     }
//     break;
// 
//     case 1: { // TET
//       std::cout << " --------------------------------------------------     TET      --------------------------------------------------" << std::endl;
// 
//       mlMsh.ReadCoarseMesh("./input/prism3D.neu", "seventh", scalingFactor);
// //   mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
//       mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);
// 
// //NOTE Tests ran with 2 procs
// //Test 1 (TET):  element 20
//       x[0] = -0.5;
//       x[1] = 0.;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a1TET(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a1TET.GetMarkerType() << std::endl;
// 
// //Test 2 (TET):  element 22
//       x[0] = 0.;
//       x[1] = 0.;
//       x[2] = 0.001;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a2TET(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a2TET.GetMarkerType() << std::endl;
// 
// 
//       //print mesh
//       MultiLevelSolution mlSol(&mlMsh);
//       variablesToBePrinted.push_back("All");
// 
//       VTKWriter vtkIO(&mlSol);
//       vtkIO.SetDebugOutput(true);
//       vtkIO.Write(Files::_application_output_directory, "biquadratic", variablesToBePrinted);
// 
//     }
//     break;
// 
//     case 2: { //
// 
//       std::cout << " --------------------------------------------------    WEDGE     --------------------------------------------------" << std::endl;
// 
//       mlMsh.ReadCoarseMesh("./input/wedge.neu", "seventh", scalingFactor);
// //   mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
//       mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);
// 
// //NOTE Tests ran with 2 procs
// //Test 1 (WEDGE):
//       x[0] = 0.0815329; //element 10
//       x[1] = 0.23;
//       x[2] = 0.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a1WEDGE(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a1WEDGE.GetMarkerType() << std::endl;
// 
// //Test 2 (WEDGE):
//       x[0] = 0.36; // if z =5 is in element 59
//       x[1] = 1.4;
//       x[2] = 5.;
// 
//       std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;
//       Marker a2WEDGE(x, VOLUME, mlMsh.GetLevel(0), solType, true);
//       //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
//       std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
//       std::cout << " The marker type is " <<  a2WEDGE.GetMarkerType() << std::endl;
// 
// 
//       //print mesh
//       MultiLevelSolution mlSol(&mlMsh);
//       variablesToBePrinted.push_back("All");
// 
//       VTKWriter vtkIO(&mlSol);
//       vtkIO.SetDebugOutput(true);
//       vtkIO.Write(Files::_application_output_directory, "biquadratic", variablesToBePrinted);
//     }
//     break;
//   }

  return 0;
}




