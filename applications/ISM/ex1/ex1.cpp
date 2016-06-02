// tests the inclusion algorithm for all 2D and 3D elements


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

    std::vector < double > x(3,0); // marker
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
    
    int elementType = 1; // this decides what elements to test

    switch(elementType) {
    case 3: // QUAD
    {

        std::cout << " --------------------------------------------------     QUAD     --------------------------------------------------" << std::endl;

        mlMsh.ReadCoarseMesh( "./input/square.neu", "seventh", scalingFactor );
//   mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
        mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

//Test 1 (QUAD): the marker is just a little bit BELOW an horizontal edge
        x[0]=0.33375;  //the marker is outside the mesh
        x[1]=-0.062500000001; //if the y was -0.625 then it is on the lower edge of element 1
        x[2]=0.;
//WARNING With more than 7 zeros the marker is considered to be ON the edge and so it would belong to element 1 anyways.

        Marker a1QUAD( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a1QUAD.GetMarkerType() <<std::endl;

//Test 2 (QUAD): the marker is just a little bit ABOVE an horizontal edge
        x[0]=0.33375;  //the marker is inside the mesh (of square.neu), it is in element 1
        x[1]=-0.062499999999; //if the y was -0.625 then it is on the lower edge of element 1,
        x[2]=0.;
// WARNING With more than 9 nines the code thinks the marker is actually ON the edge.

        Marker a2QUAD( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a2QUAD.GetMarkerType() <<std::endl;

//Test 3 (QUAD): the marker is just a little bit to the LEFT of a vertical edge
        x[0]=0.3124999999999;  //the marker is on element 0
        x[1]=-0.05555555; // if the x is 0.3125 then it is on the vertical edge of element 1
        x[2]=0.; //
// WARNING With more than 8 nines the code considers the marker to be ON the edge and so it is considered to be on element 0
// only because that edge appears first in element 0 than in element 1.
// In this case it is ok because the element numbers increase from left to right.

        Marker a3QUAD( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a3QUAD.GetMarkerType() <<std::endl;

//Test 4 (QUAD): the marker is just a little bit to the RIGHT of a vertical edge
        x[0]=0.312500000001;  //the marker is on element 1
        x[1]=-0.05555555; // if the x is 0.3125 then it is on the vertical edge of element 1
        x[2]=0.;
// WARNING With more than 7 zeros the marker is considered to be ON the vertical edge, which belongs first to element 0
// this is why it is considered to be on element 0 instead of on element 1.

        Marker a4QUAD( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a4QUAD.GetMarkerType() <<std::endl;

//Test 5 (QUAD): the marker is on the lower half of the RIGHT edge of element 0
        x[0]=0.3125;
        x[1]=-0.05555555;
        x[2]=0.;

        Marker a5QUAD( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a5QUAD.GetMarkerType() <<std::endl;

//Test 6 (QUAD): the marker is the north-east vertex of element 43
        x[0]=0.25;
        x[1]=0.3125;
        x[2]=0.;

        Marker a6QUAD( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a6QUAD.GetMarkerType() <<std::endl;

// Test 7 (QUAD): the marker is on element 35 but the y coordinate is the same as the MIDPOINT of the RIGHT edge of element 35
        x[0]=0.23;
        x[1]=0.28125;
        x[2]=0.;

        Marker a7QUAD( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a7QUAD.GetMarkerType() <<std::endl;


        //print mesh
        MultiLevelSolution mlSol( &mlMsh );
        variablesToBePrinted.push_back( "All" );

        VTKWriter vtkIO( &mlSol );
        vtkIO.SetDebugOutput( true );
        vtkIO.Write( DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted );
	
	

    }
    break;

    case 4: // TRI
    {
        std::cout << " --------------------------------------------------     TRI      --------------------------------------------------" << std::endl;

        mlMsh.ReadCoarseMesh( "./input/square2.neu", "seventh", scalingFactor );
//   mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
        mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

//Test 1 (TRI): the maker is VERY close to a vertex of element 49 (point 37)
        x[0]=1.61638e-05;  // it is on element 49
        x[1]=7.81317e-05;
        x[2]=0.;

        Marker a1TRI( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a1TRI.GetMarkerType() <<std::endl;

//Test 2 (TRI): the maker is VERY close to a midpoint of an edge of element 60 (point 144)
        x[0]=5.77398e-05;  //it is on element 60
        x[1]=7.23015e-05;
        x[2]=0.;

        Marker a2TRI( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a2TRI.GetMarkerType() <<std::endl;

// Test 3 (TRI): the maker has the "same" y coordinate of a midpoint of an edge of element 60 (point 144)
        x[0]=5.77298e-05;  //it is on element 60
        x[1]=7.23015e-05;
        x[2]=0.;

        Marker a3TRI( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a3TRI.GetMarkerType() <<std::endl;

// Test 4 (TRI): the maker is VERY close to a vertex of element 15 (point 19)
        x[0]=-2.72848e-05;  //it is on element 15
        x[1]=2.62822e-05;
        x[2]=0.;

        Marker a4TRI( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a4TRI.GetMarkerType() <<std::endl;


        //print mesh
        MultiLevelSolution mlSol( &mlMsh );
        variablesToBePrinted.push_back( "All" );

        VTKWriter vtkIO( &mlSol );
        vtkIO.SetDebugOutput( true );
        vtkIO.Write( DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted );

    }
 break;
    
    
    case 0: // HEX

    {

        std::cout << " --------------------------------------------------     HEX      --------------------------------------------------" << std::endl;

        mlMsh.ReadCoarseMesh( "./input/cubeHex.neu", "seventh", scalingFactor );
//   mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
        mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

//NOTE Tests ran with 2 procs
//Test 1 (HEX):
        x[0]=-0.3;  //  point 1409 shared by element 97,112 (proc 2) and 45 and 55 (proc 1) It says the marker is on element 45, yes!
        x[1]=0.1; //
        x[2]=0.5;

        Marker a1HEX( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a1HEX.GetMarkerType() <<std::endl;

//Test 2 (HEX): the marker is on the midpoint of an edge of element 60 (the trick was that the intersection point was actually ON a midpoint of an edge)
        x[0]=-0.5;
        x[1]=0.5;
        x[2]=0.;

        Marker a2HEX( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a2HEX.GetMarkerType() <<std::endl;

//Test 3 (HEX): the marker is to the right of the midpoint of an edge of element 60
        x[0]=-0.5;
        x[1]=0.5;
        x[2]=0.005;

        Marker a3HEX( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a3HEX.GetMarkerType() <<std::endl;

//Test 4 (HEX): the marker is a vertex of element 60
        x[0]=-0.5;
        x[1]=0.5;
        x[2]=-0.1;

        Marker a4HEX( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a4HEX.GetMarkerType() <<std::endl;

        //Test 5 (HEX): the marker is close to a vertex of element 60 but it belongs to element 62
        x[0]=-0.5;
        x[1]=0.5;
        x[2]=-0.10000001;  // If you add one more zero then it will be in element 60.

        Marker a5HEX( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a5HEX.GetMarkerType() <<std::endl;

//Test 6 (HEX): the marker is a little below a vertex of element 60 (it goes outside the domain)
        x[0]=-0.5000001;
        x[1]=0.5;
        x[2]=-0.05;

        Marker a6HEX( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a6HEX.GetMarkerType() <<std::endl;

//Test 7 (HEX): the marker is over the lower edge of element 57 (it is not on the domain because the x is too small)
        x[0]=-0.6;
        x[1]=0.1;
        x[2]=0.2;

        Marker a7HEX( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a7HEX.GetMarkerType() <<std::endl;


        //print mesh
        MultiLevelSolution mlSol( &mlMsh );
        variablesToBePrinted.push_back( "All" );

        VTKWriter vtkIO( &mlSol );
        vtkIO.SetDebugOutput( true );
        vtkIO.Write( DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted );

    }
break;

    case 1: // TET
    {
        std::cout << " --------------------------------------------------     TET      --------------------------------------------------" << std::endl;

        mlMsh.ReadCoarseMesh( "./input/prism3D.neu", "seventh", scalingFactor );
//   mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
        mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

//NOTE Tests ran with 2 procs
//Test 1 (TET):  element 20
        x[0] = -0.5;
        x[1] = 0.;
        x[2] = 0.;

        Marker a1TET( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a1TET.GetMarkerType() <<std::endl;

//Test 2 (TET):  element 22
        x[0] = 0.;
        x[1] = 0.;
        x[2] = 0.001;

        Marker a2TET( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a2TET.GetMarkerType() <<std::endl;


        //print mesh
        MultiLevelSolution mlSol( &mlMsh );
        variablesToBePrinted.push_back( "All" );

        VTKWriter vtkIO( &mlSol );
        vtkIO.SetDebugOutput( true );
        vtkIO.Write( DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted );

    }
break;
    
    case 2: //WEDGE
    {

        std::cout << " --------------------------------------------------    WEDGE     --------------------------------------------------" << std::endl;

        mlMsh.ReadCoarseMesh( "./input/wedge.neu", "seventh", scalingFactor );
//   mlMsh.ReadCoarseMesh("./input/quadAMR.neu", "seventh", scalingFactor);
        mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

//NOTE Tests ran with 2 procs
//Test 1 (WEDGE):
        x[0] = -0.5;
        x[1] = 0.;
        x[2] = 0.;

        Marker a1WEDGE( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a1WEDGE.GetMarkerType() <<std::endl;

//Test 2 (WEDGE):
        x[0] = 0.;
        x[1] = 0.;
        x[2] = 0.001;

        Marker a2WEDGE( x, VOLUME, mlMsh.GetLevel(0), true );
        //Marker a( x, VOLUME, mlMsh.GetLevel(numberOfUniformLevels + numberOfSelectiveLevels -1) );
        std::cout<< " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," <<x[2]<<std::endl;
        std::cout << " The marker type is " <<  a2WEDGE.GetMarkerType() <<std::endl;


        //print mesh
        MultiLevelSolution mlSol( &mlMsh );
        variablesToBePrinted.push_back( "All" );

        VTKWriter vtkIO( &mlSol );
        vtkIO.SetDebugOutput( true );
        vtkIO.Write( DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted );
    }
    break;
    }

    return 0;
}




