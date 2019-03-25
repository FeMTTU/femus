// tests the inclusion algorithm for all 2D and 3D elements


#include "FemusInit.hpp"
#include "Marker.hpp"
#include "MultiLevelMesh.hpp"
#include "MultiLevelSolution.hpp"
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

    // init Petsc-MPI communicator
    FemusInit mpinit(argc, args, MPI_COMM_WORLD);

    std::vector < double > x(3, 0); // marker
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
    unsigned solType = 2;

    for(int elementType = 0; elementType < 1; elementType++) {
        if(elementType == 0) {
            std::cout << " --------------------------------------------------     HEX      --------------------------------------------------" << std::endl;
//             mlMsh.ReadCoarseMesh("./input/cubeHex.neu", "seventh", scalingFactor);
// 	    x[0] = -0.3; //  point 98 shared by element 97,112 (proc 1) and 45 and 55 (proc 0). Proc 1 says the marker is in 45, proc 1 says it is is 112.
//             x[1] = 0.1; //
//             x[2] = 0.5;
	    
	    mlMsh.ReadCoarseMesh("./input/aneurysm_Sara_5.neu", "seventh", scalingFactor);
            x[0] = 0.; 
            x[1] = 0.; 
            x[2] = 0.;
        }
        else if(elementType == 1) {
            std::cout << " --------------------------------------------------     TET      --------------------------------------------------" << std::endl;
            //mlMsh.ReadCoarseMesh("./input/prism3D.neu", "seventh", scalingFactor);
            mlMsh.ReadCoarseMesh("./input/torus.neu", "seventh", scalingFactor); //this is a torus
	    x[0] = -0.5;
            x[1] = 0.;
            x[2] = 0.;
        }
        else if(elementType == 2) {
            std::cout << " --------------------------------------------------    WEDGE     --------------------------------------------------" << std::endl;
            mlMsh.ReadCoarseMesh("./input/wedge.neu", "seventh", scalingFactor);
            x[0] = 0.0815329; //element 10
            x[1] = 0.23;
            x[2] = 0.;
        }
        else if(elementType == 3) {
            std::cout << " --------------------------------------------------     QUAD     --------------------------------------------------" << std::endl;
//             mlMsh.ReadCoarseMesh("./input/square.neu", "seventh", scalingFactor);
//             x[0] = 0.33375; //the marker is in element 117 (proc 1)
//             x[1] = -0.0627;
//             x[2] = 0.;
	    
// 	    mlMsh.ReadCoarseMesh("./input/distortedQuad.neu", "seventh", scalingFactor);
//             x[0] = 0.; 
//             x[1] = -0.;
//             x[2] = 0.;
// 	    
	    mlMsh.ReadCoarseMesh("./input/rotatedQuad.neu", "seventh", scalingFactor);
            x[0] = 0.; 
            x[1] = -0.;
            x[2] = 0.;
	    
        }
        else if(elementType == 4) {
            std::cout << " --------------------------------------------------     TRI      --------------------------------------------------" << std::endl;
            mlMsh.ReadCoarseMesh("./input/square2.neu", "seventh", scalingFactor);
            x[0] = 5.e-05; // it is on element 12
            x[1] = 0.00015;
            x[2] = 0.;
        }
        mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , SetRefinementFlag);

        


        //print mesh
        MultiLevelSolution mlSol(&mlMsh);
	
	
	
	Marker mrk(x, 0.,  VOLUME, mlSol.GetLevel(0), solType, true);
        std::cout << " The coordinates of the marker are " << x[0] << " ," << x[1] << " ," << x[2] << std::endl;
        std::cout << " The marker type is " <<  mrk.GetMarkerType() << std::endl;


       // mrk.InverseMappingTEST(x, mlSol.GetLevel(0));
	
	
	
	
        variablesToBePrinted.push_back("All");

        VTKWriter vtkIO(&mlSol);
        vtkIO.SetDebugOutput(true);
        vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted,elementType);

    } // this decides what elements to test

    return 0;
}




