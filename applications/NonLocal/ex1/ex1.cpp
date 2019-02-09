
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

#include "slepceps.h"

#include "../include/nonlocal_assembly.hpp"


//FIRST NONLOCAL EX IN FEMUS: nonlocal diffusion for a body with different material properties

using namespace femus;


bool SetBoundaryCondition ( const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time )
{

    bool dirichlet = true;
//     value = 0.;
    value = x[0];

    if ( facename == 2 ) {
        bool dirichlet = false; //Neumann at the interface boundaries
    }

    return dirichlet;
}

unsigned numberOfUniformLevels = 2;

int main ( int argc, char** argv )
{


   // init Petsc-MPI communicator
    FemusInit mpinit ( argc, argv, MPI_COMM_WORLD );

    MultiLevelMesh mlMsh;
    double scalingFactor = 1.;
    unsigned numberOfSelectiveLevels = 0;
//     mlMsh.ReadCoarseMesh ( "../input/nonlocal_boundary_test.neu", "second", scalingFactor );
//     mlMsh.ReadCoarseMesh ( "../input/interface.neu", "second", scalingFactor );
    mlMsh.ReadCoarseMesh ( "../input/maxTest1.neu", "second", scalingFactor );
//     mlMsh.ReadCoarseMesh ( "../input/martaTest1.neu", "second", scalingFactor );
//     mlMsh.ReadCoarseMesh ( "../input/martaTest2.neu", "second", scalingFactor );
//         mlMsh.ReadCoarseMesh ( "../input/martaTest3.neu", "second", scalingFactor );
//     mlMsh.ReadCoarseMesh ( "../input/trial1.neu", "second", scalingFactor );
//     mlMsh.ReadCoarseMesh ( "../input/trial2.neu", "second", scalingFactor );
    mlMsh.RefineMesh ( numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL );
    
//         mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);
//     numberOfUniformLevels = 1;

    unsigned dim = mlMsh.GetDimension();

    MultiLevelSolution mlSol ( &mlMsh );

    // add variables to mlSol
    mlSol.AddSolution ( "u", LAGRANGE, FIRST, 2 );

    mlSol.Initialize ( "All" );

    mlSol.AttachSetBoundaryConditionFunction ( SetBoundaryCondition );

    // ******* Set boundary conditions *******
    mlSol.GenerateBdc ( "All" );

    MultiLevelProblem ml_prob ( &mlSol );

    // ******* Add FEM system to the MultiLevel problem *******
    LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ( "NonLocal" );
    system.AddSolutionToSystemPDE ( "u" );

    // ******* System FEM Assembly *******
    system.SetAssembleFunction ( AssembleNonLocalSys );
    system.SetMaxNumberOfLinearIterations ( 1 );
    // ******* set MG-Solver *******
    system.SetMgType ( V_CYCLE );

    system.SetAbsoluteLinearConvergenceTolerance ( 1.e-50 );
    //   system.SetNonLinearConvergenceTolerance(1.e-9);
//   system.SetMaxNumberOfNonLinearIterations(20);

    system.SetNumberPreSmoothingStep ( 1 );
    system.SetNumberPostSmoothingStep ( 1 );

    // ******* Set Preconditioner *******
    system.SetMgSmoother ( GMRES_SMOOTHER );

    system.SetSparsityPatternMultiplyingFactor ( 1000u ); //TODO tune 10u

    system.init();

    // ******* Set Smoother *******
    system.SetSolverFineGrids ( GMRES );

    system.SetPreconditionerFineGrids ( ILU_PRECOND );

    system.SetTolerances ( 1.e-20, 1.e-20, 1.e+50, 100 );

    system.MGsolve();

    // ******* Print solution *******
    mlSol.SetWriter ( VTK );
    std::vector<std::string> print_vars;
    print_vars.push_back ( "All" );
    mlSol.GetWriter()->SetDebugOutput ( true );
    mlSol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0 );

    return 0;

} //end main










