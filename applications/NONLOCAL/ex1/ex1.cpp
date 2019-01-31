
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
    bool dirichlet = true; //Neumann
    value = 0.;
    return dirichlet;
}

double L = 0.1;

void GetEigenPair ( MultiLevelProblem& ml_prob, const int& numberOfEigPairs, std::vector < std::pair<double, double> >& eigenvalues );

void IsElementInBall ( unsigned &check, const std::vector<double> &ballCenter, const double &ballRadius, const std::vector <double> &elementCoordinates );

unsigned numberOfUniformLevels = 2;

int main ( int argc, char** argv )
{


    //BEGIN eigenvalue problem instances
    PetscErrorCode ierr;
    ierr = SlepcInitialize ( &argc, &argv, PETSC_NULL, PETSC_NULL );

//   numberOfEigPairs = 5; //number of eigenpairs desired

    eigenvalues.resize ( numberOfEigPairs ); //this is where we store the eigenvalues

    //END


    //BEGIN deterministic FEM instances

    // init Petsc-MPI communicator
    FemusInit mpinit ( argc, argv, MPI_COMM_WORLD );

    MultiLevelMesh mlMsh;
    double scalingFactor = 1.;
    unsigned numberOfSelectiveLevels = 0;
//   mlMsh.ReadCoarseMesh("../input/square.neu", "fifth", scalingFactor);
    mlMsh.ReadCoarseMesh ( "../input/nonlocal_boundary_test.neu", "fifth", scalingFactor );
    mlMsh.RefineMesh ( numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL );

    unsigned dim = mlMsh.GetDimension();

    MultiLevelSolution mlSol ( &mlMsh );

    // add variables to mlSol
    mlSol.AddSolution ( "u", LAGRANGE, FIRST, 2 );

    for ( unsigned i = 0; i < numberOfEigPairs; i++ ) {
        char name[10];
        sprintf ( name, "egnf%d", i );
        mlSol.AddSolution ( name, LAGRANGE, FIRST, 0, false );
    }

    mlSol.Initialize ( "All" );

    mlSol.AttachSetBoundaryConditionFunction ( SetBoundaryCondition );

    // ******* Set boundary conditions *******
    mlSol.GenerateBdc ( "All" );

    MultiLevelProblem ml_prob ( &mlSol );

    // ******* Add FEM system to the MultiLevel problem *******
    LinearImplicitSystem& system = ml_prob.add_system < LinearImplicitSystem > ( "UQ" );
    system.AddSolutionToSystemPDE ( "u" );

    // ******* System FEM Assembly *******
    system.SetAssembleFunction ( AssembleNonlocalSys );
    system.SetMaxNumberOfLinearIterations ( 1 );
    //system.SetAssembleFunction(AssembleFEM);
    // ******* set MG-Solver *******
    system.SetMgType ( V_CYCLE );

    system.SetAbsoluteLinearConvergenceTolerance ( 1.e-50 );
    //   system.SetNonLinearConvergenceTolerance(1.e-9);
//   system.SetMaxNumberOfNonLinearIterations(20);

    system.SetNumberPreSmoothingStep ( 1 );
    system.SetNumberPostSmoothingStep ( 1 );

    // ******* Set Preconditioner *******
    system.SetMgSmoother ( GMRES_SMOOTHER );

    system.init();

    // ******* Set Smoother *******
    system.SetSolverFineGrids ( GMRES );

    system.SetPreconditionerFineGrids ( ILU_PRECOND );

    system.SetTolerances ( 1.e-20, 1.e-20, 1.e+50, 100 );
    //END

//     GetEigenPair ( ml_prob, numberOfEigPairs, eigenvalues ); //solve the generalized eigenvalue problem and compute the eigenpairs

    for ( int i = 0; i < numberOfEigPairs; i++ ) {
        std::cout << eigenvalues[i].first << " " << eigenvalues[i].second << std::endl;
    }

    system.MGsolve();


    // ******* Print solution *******
    mlSol.SetWriter ( VTK );
    std::vector<std::string> print_vars;
    print_vars.push_back ( "All" );
    mlSol.GetWriter()->SetDebugOutput ( true );
    mlSol.GetWriter()->Write ( DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0 );

    //ierr = SlepcFinalize();
    //CHKERRQ(ierr);

    return 0;

} //end main

void GetEigenPair ( MultiLevelProblem& ml_prob, const int& numberOfEigPairs, std::vector < std::pair<double, double> >& eigenvalues )
{

    LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ( "UQ" ); // pointer to the linear implicit system named "Poisson"

    unsigned level = numberOfUniformLevels - 1;

    double varianceInput = stdDeviationInput * stdDeviationInput;

    Mesh*                    msh = ml_prob._ml_msh->GetLevel ( level ); // pointer to the mesh (level) object
    elem*                     el = msh->el;  // pointer to the elem object in msh (level)

    MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
    Solution*                sol = ml_prob._ml_sol->GetSolutionLevel ( level ); // pointer to the solution (level) object

    LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
    SparseMatrix*             MM = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)

    const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
    const unsigned maxSize = static_cast< unsigned > ( ceil ( pow ( 3, dim ) ) ); // conservative: based on line3, quad9, hex27

    unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
    unsigned    nprocs = msh->n_processors(); // get the process_id (for parallel computation)


    //solution variable
    unsigned soluIndex;
    soluIndex = mlSol->GetIndex ( "u" ); // get the position of "u" in the ml_sol object
    unsigned solType = mlSol->GetSolutionType ( soluIndex ); // get the finite element type for "u"

    unsigned soluPdeIndex;
    soluPdeIndex = mlPdeSys->GetSolPdeIndex ( "u" ); // get the position of "u" in the pdeSys object

    unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

    vector < vector < double > > x1 ( dim ); // local coordinates
    vector < vector < double > > x2 ( dim ); // local coordinates

    for ( unsigned k = 0; k < dim; k++ ) {
        x1[k].reserve ( maxSize );
        x2[k].reserve ( maxSize );
    }

    vector <double> phi_x; // local test function first order partial derivatives

    phi_x.reserve ( maxSize * dim );

    vector< int > l2GMap1; // local to global mapping
    vector< int > l2GMap2; // local to global mapping
    l2GMap1.reserve ( maxSize );
    l2GMap2.reserve ( maxSize );

    vector < double > MMlocal;
    MMlocal.reserve ( maxSize * maxSize );

    vector < double > CClocal;
    CClocal.reserve ( maxSize * maxSize );

    MM->zero(); // Set to zero all the entries of the Global Matrix

    int MM_size = msh->_dofOffset[solType][nprocs];
    int MM_local_size = msh->_dofOffset[solType][iproc + 1] - msh->_dofOffset[solType][iproc];

    SparseMatrix* CC;
    CC = SparseMatrix::build().release();
    CC->init ( MM_size, MM_size, MM_local_size, MM_local_size, MM_local_size, MM_size - MM_local_size );
    CC->zero();

    for ( int kproc = 0; kproc < nprocs; kproc++ ) {
        for ( int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++ ) {

            short unsigned ielGeom2;
            unsigned nDof2;
            unsigned nDofx2;

            if ( iproc == kproc ) {
                ielGeom2 = msh->GetElementType ( jel );
                nDof2  = msh->GetElementDofNumber ( jel, solType ); // number of solution element dofs
                nDofx2 = msh->GetElementDofNumber ( jel, xType ); // number of coordinate element dofs
            }

            MPI_Bcast ( &ielGeom2, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD );
            MPI_Bcast ( &nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD );
            MPI_Bcast ( &nDofx2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD );

            // resize local arrays
            l2GMap2.resize ( nDof2 );

            for ( int k = 0; k < dim; k++ ) {
                x2[k].resize ( nDofx2 );
            }

            // local storage of global mapping and solution
            if ( iproc == kproc ) {
                for ( unsigned j = 0; j < nDof2; j++ ) {
                    l2GMap2[j] = pdeSys->GetSystemDof ( soluIndex, soluPdeIndex, j, jel ); // global to global mapping between solution node and pdeSys dof
                }
            }

            MPI_Bcast ( &l2GMap2[0], nDof2, MPI_UNSIGNED, kproc, MPI_COMM_WORLD );

            // local storage of coordinates
            if ( iproc == kproc ) {
                for ( unsigned j = 0; j < nDofx2; j++ ) {
                    unsigned xDof  = msh->GetSolutionDof ( j, jel, xType ); // global to global mapping between coordinates node and coordinate dof

                    for ( unsigned k = 0; k < dim; k++ ) {
                        x2[k][j] = ( *msh->_topology->_Sol[k] ) ( xDof ); // global extraction and local storage for the element coordinates
                    }
                }
            }

            for ( unsigned k = 0; k < dim; k++ ) {
                MPI_Bcast ( & x2[k][0], nDofx2, MPI_DOUBLE, kproc, MPI_COMM_WORLD );
            }

            unsigned jgNumber = msh->_finiteElement[ielGeom2][solType]->GetGaussPointNumber();
            vector < vector < double > > xg2 ( jgNumber );
            vector <double> weight2 ( jgNumber );
            vector < vector <double> > phi2 ( jgNumber ); // local test function

            for ( unsigned jg = 0; jg < jgNumber; jg++ ) {
                msh->_finiteElement[ielGeom2][solType]->Jacobian ( x2, jg, weight2[jg], phi2[jg], phi_x );

                xg2[jg].assign ( dim, 0. );

                for ( unsigned j = 0; j < nDof2; j++ ) {
                    for ( unsigned k = 0; k < dim; k++ ) {
                        xg2[jg][k] += x2[k][j] * phi2[jg][j];
                    }
                }
            }

            // element loop: each process loops only on the elements that owns
            for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

                short unsigned ielGeom1 = msh->GetElementType ( iel );
                unsigned nDof1  = msh->GetElementDofNumber ( iel, solType ); // number of solution element dofs
                unsigned nDofx1 = msh->GetElementDofNumber ( iel, xType ); // number of coordinate element dofs

                // resize local arrays
                l2GMap1.resize ( nDof1 );
                //std::vector<bool>bdcDirichlet(nDof1);

                for ( int k = 0; k < dim; k++ ) {
                    x1[k].resize ( nDofx1 );
                }

                // local storage of global mapping and solution
                for ( unsigned i = 0; i < nDof1; i++ ) {
                    l2GMap1[i] = pdeSys->GetSystemDof ( soluIndex, soluPdeIndex, i, iel ); // global to global mapping between solution node and pdeSys dof
                    //unsigned solDof = msh->GetSolutionDof(i, iel, solType);    // global to global mapping between solution node and solution dof
                    //bdcDirichlet[i] = ( (*sol->_Bdc[soluIndex])(solDof) < 1.5)? false:false;
                }

                // local storage of coordinates
                for ( unsigned i = 0; i < nDofx1; i++ ) {
                    unsigned xDof  = msh->GetSolutionDof ( i, iel, xType ); // global to global mapping between coordinates node and coordinate dof

                    for ( unsigned k = 0; k < dim; k++ ) {
                        x1[k][i] = ( *msh->_topology->_Sol[k] ) ( xDof ); // global extraction and local storage for the element coordinates
                    }
                }

                if ( iel == jel ) MMlocal.assign ( nDof1 * nDof1, 0. ); //resize

                CClocal.assign ( nDof1 * nDof2, 0. ); //resize

                // *** Gauss point loop ***
                unsigned igNumber = msh->_finiteElement[ielGeom1][solType]->GetGaussPointNumber();
                double weight1;
                vector <double> phi1;  // local test function

                for ( unsigned ig = 0; ig < igNumber; ig++ ) {

                    msh->_finiteElement[ielGeom1][solType]->Jacobian ( x1, ig, weight1, phi1, phi_x );

                    // evaluate the solution, the solution derivatives and the coordinates in the gauss point
                    vector < double > xg1 ( dim, 0. );

                    for ( unsigned i = 0; i < nDof1; i++ ) {
                        for ( unsigned k = 0; k < dim; k++ ) {
                            xg1[k] += x1[k][i] * phi1[i];
                        }
                    }

                    if ( iel == jel ) {
                        for ( unsigned i = 0; i < nDof1; i++ ) {
                            for ( unsigned i1 = 0; i1 < nDof1; i1++ ) {
                                MMlocal[ i * nDof1 + i1 ] += phi1[i] * phi1[i1] * weight1;
                            }
                        }
                    }

                    for ( unsigned jg = 0; jg < jgNumber; jg++ ) {
                        double dist = 0.;

                        for ( unsigned k = 0; k < dim; k++ ) {
                            dist += fabs ( xg1[k] - xg2[jg][k] );
                        }

                        double C = varianceInput * exp ( - dist / L );

                        for ( unsigned i = 0; i < nDof1; i++ ) {
                            for ( unsigned j = 0; j < nDof2; j++ ) {
                                CClocal[i * nDof2 + j] += weight1 * phi1[i] * C * phi2[jg][j] * weight2[jg];
                            }//endl j loop
                        } //endl i loop
                    } //endl jg loop
                } //endl ig loop

                if ( iel == jel ) MM->add_matrix_blocked ( MMlocal, l2GMap1, l2GMap1 );

                CC->add_matrix_blocked ( CClocal, l2GMap1, l2GMap2 );
            } // end iel loop
        } //end jel loop
    } //end kproc loop

    MM->close();
    CC->close();

    //BEGIN solve the eigenvalue problem

    int ierr;
    EPS eps;
    PetscInt convergedSolns, numberOfIterations;

    ierr = EPSCreate ( PETSC_COMM_WORLD, &eps );
    CHKERRABORT ( MPI_COMM_WORLD, ierr );
    ierr = EPSSetOperators ( eps, ( static_cast<PetscMatrix*> ( CC ) )->mat(), ( static_cast<PetscMatrix*> ( MM ) )->mat() );
    CHKERRABORT ( MPI_COMM_WORLD, ierr );
    ierr = EPSSetFromOptions ( eps );
    CHKERRABORT ( MPI_COMM_WORLD, ierr );

    //ierr = EPSSetDimensions(eps, numberOfEigPairs, 8 * numberOfEigPairs, 600);
    ierr = EPSSetDimensions ( eps, numberOfEigPairs, PETSC_DEFAULT, PETSC_DEFAULT );
    CHKERRABORT ( MPI_COMM_WORLD, ierr );
    ierr = EPSSetWhichEigenpairs ( eps, EPS_LARGEST_MAGNITUDE );
    CHKERRABORT ( MPI_COMM_WORLD, ierr );

    //ierr = EPSSetTolerances(eps,1.0e-10,1000);
    CHKERRABORT ( MPI_COMM_WORLD, ierr );

    ierr = EPSSolve ( eps );
    CHKERRABORT ( MPI_COMM_WORLD, ierr );

    //ierr = EPSView(eps, PETSC_VIEWER_STDOUT_SELF);

    std::cout << " -----------------------------------------------------------------" << std::endl;

    ierr = EPSGetConverged ( eps, &convergedSolns );
    CHKERRABORT ( MPI_COMM_WORLD, ierr );
    ierr = PetscPrintf ( PETSC_COMM_WORLD, " Number of converged eigenpairs: %D\n\n", convergedSolns );
    CHKERRABORT ( MPI_COMM_WORLD, ierr );

    if ( convergedSolns > 0 ) {

        for ( unsigned i = 0; i < numberOfEigPairs; i++ ) {

            char name[10];
            sprintf ( name, "egnf%d", i );
            soluIndex = mlSol->GetIndex ( name ); // get the position of "u" in the ml_sol object

            // Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and ki (imaginary part)

            ierr = EPSGetEigenpair ( eps, i, &eigenvalues[i].first, &eigenvalues[i].second, ( static_cast<PetscVector*> ( sol->_Sol[soluIndex] ) )->vec(), NULL );
            CHKERRABORT ( MPI_COMM_WORLD, ierr );

        }
    }

    ierr = EPSDestroy ( &eps );
    CHKERRABORT ( MPI_COMM_WORLD, ierr );

    delete CC;


    //BEGIN GRAM SCHMIDT ORTHONORMALIZATION

    std::vector <unsigned> eigfIndex ( numberOfEigPairs );
    char name[10];

    for ( unsigned i = 0; i < numberOfEigPairs; i++ ) {
        sprintf ( name, "egnf%d", i );
        eigfIndex[i] = mlSol->GetIndex ( name ); // get the position of "u" in the ml_sol object
    }

    vector < double >  eigenFunction ( numberOfEigPairs ); // local solution
    vector < double >  eigenFunctionOld ( numberOfEigPairs ); // local solution

    std::vector < std::vector < double > > coeffsGS_local ( numberOfEigPairs );
    std::vector < std::vector < double > > coeffsGS_global ( numberOfEigPairs );

    for ( unsigned i = 0; i < numberOfEigPairs; i++ ) {
        coeffsGS_local[i].assign ( numberOfEigPairs, 0. );
        coeffsGS_global[i].assign ( numberOfEigPairs, 0. );
    }

    for ( unsigned iGS = 0; iGS < numberOfEigPairs; iGS++ ) {

        if ( iGS > 0 ) {

            for ( unsigned jGS = 0; jGS < iGS; jGS++ ) {

                //BEGIN COMPUTE coeffsGS LOCAL

                for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

                    short unsigned ielGeom = msh->GetElementType ( iel );
                    unsigned nDofu  = msh->GetElementDofNumber ( iel, solType ); // number of solution element dofs
                    unsigned nDofx = msh->GetElementDofNumber ( iel, xType ); // number of coordinate element dofs

                    eigenFunction.resize ( nDofu );
                    eigenFunctionOld.resize ( nDofu );

                    for ( int i = 0; i < dim; i++ ) {
                        x1[i].resize ( nDofx );
                    }

                    // local storage of global mapping and solution
                    for ( unsigned i = 0; i < nDofu; i++ ) {
                        unsigned solDof = msh->GetSolutionDof ( i, iel, solType ); // global to global mapping between solution node and solution dof
                        eigenFunction[i] = ( *sol->_Sol[eigfIndex[iGS]] ) ( solDof );
                        eigenFunctionOld[i] = ( *sol->_Sol[eigfIndex[jGS]] ) ( solDof );
                    }

                    // local storage of coordinates
                    for ( unsigned i = 0; i < nDofx; i++ ) {
                        unsigned xDof  = msh->GetSolutionDof ( i, iel, xType ); // global to global mapping between coordinates node and coordinate dof

                        for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                            x1[jdim][i] = ( *msh->_topology->_Sol[jdim] ) ( xDof ); // global extraction and local storage for the element coordinates
                        }
                    }

                    double weight;
                    vector <double> phi;  // local test function

                    // *** Gauss point loop ***
                    for ( unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++ ) {
                        // *** get gauss point weight, test function and test function partial derivatives ***
                        msh->_finiteElement[ielGeom][solType]->Jacobian ( x1, ig, weight, phi, phi_x );
                        double eigenFunction_gss = 0.;
                        double eigenFunction_gss_old = 0.;

                        for ( unsigned i = 0; i < nDofu; i++ ) {
                            eigenFunction_gss += phi[i] * eigenFunction[i];
                            eigenFunction_gss_old += phi[i] * eigenFunctionOld[i];
                        }

                        coeffsGS_local[iGS][jGS] += eigenFunction_gss * eigenFunction_gss_old * weight;
                    }
                }

                //END COMPUTE coeffsGS LOCAL

                MPI_Allreduce ( &coeffsGS_local[iGS][jGS], &coeffsGS_global[iGS][jGS], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

            }

            for ( unsigned idof = msh->_dofOffset[solType][iproc]; idof < msh->_dofOffset[solType][iproc + 1]; idof++ ) {
                double sum = 0.;

                for ( unsigned jGS = 0; jGS < iGS; jGS++ ) {
                    sum += coeffsGS_global[iGS][jGS] * ( *sol->_Sol[eigfIndex[jGS]] ) ( idof );
                }

                double valueToSet = ( *sol->_Sol[eigfIndex[iGS]] ) ( idof ) - sum;
                sol->_Sol[eigfIndex[iGS]]->set ( idof, valueToSet );
            }

        }

        sol->_Sol[eigfIndex[iGS]]->close();

        double local_norm2 = 0.;

        for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

            short unsigned ielGeom = msh->GetElementType ( iel );
            unsigned nDofu  = msh->GetElementDofNumber ( iel, solType ); // number of solution element dofs
            unsigned nDofx = msh->GetElementDofNumber ( iel, xType ); // number of coordinate element dofs

            eigenFunction.resize ( nDofu );

            for ( int i = 0; i < dim; i++ ) {
                x1[i].resize ( nDofx );
            }

            // local storage of global mapping and solution
            for ( unsigned i = 0; i < nDofu; i++ ) {
                unsigned solDof = msh->GetSolutionDof ( i, iel, solType ); // global to global mapping between solution node and solution dof
                eigenFunction[i] = ( *sol->_Sol[eigfIndex[iGS]] ) ( solDof );
            }

            // local storage of coordinates
            for ( unsigned i = 0; i < nDofx; i++ ) {
                unsigned xDof  = msh->GetSolutionDof ( i, iel, xType ); // global to global mapping between coordinates node and coordinate dof

                for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                    x1[jdim][i] = ( *msh->_topology->_Sol[jdim] ) ( xDof ); // global extraction and local storage for the element coordinates
                }
            }

            double weight;
            vector <double> phi;  // local test function

            // *** Gauss point loop ***
            for ( unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++ ) {
                // *** get gauss point weight, test function and test function partial derivatives ***
                msh->_finiteElement[ielGeom][solType]->Jacobian ( x1, ig, weight, phi, phi_x );
                double eigenFunction_gss = 0.;

                for ( unsigned i = 0; i < nDofu; i++ ) {
                    eigenFunction_gss += phi[i] * eigenFunction[i];
                }

                local_norm2 += eigenFunction_gss * eigenFunction_gss * weight;
            }
        }

        double norm2 = 0.;
        MPI_Allreduce ( &local_norm2, &norm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
        double norm = sqrt ( norm2 );
        std::cout << "norm = " << norm << std::endl;
        sol->_Sol[eigfIndex[iGS]]->scale ( 1. / norm );

        sol->_Sol[eigfIndex[iGS]]->close();

    }

    //END GRAM SCHMIDT ORTHONORMALIZATION

    //BEGIN GRAM SCHMIDT CHECK
    vector < double >  eigenFunctionCheck ( numberOfEigPairs ); // local solution
    vector < double >  eigenFunctionOldCheck ( numberOfEigPairs ); // local solution


    for ( unsigned i1 = 0; i1 < numberOfEigPairs; i1++ ) {
        for ( unsigned j1 = 0; j1 < numberOfEigPairs; j1++ ) {

            double integral = 0.;

            for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

                short unsigned ielGeom = msh->GetElementType ( iel );
                unsigned nDofu  = msh->GetElementDofNumber ( iel, solType ); // number of solution element dofs
                unsigned nDofx = msh->GetElementDofNumber ( iel, xType ); // number of coordinate element dofs

                eigenFunctionCheck.resize ( nDofu );
                eigenFunctionOldCheck.resize ( nDofu );

                for ( int i = 0; i < dim; i++ ) {
                    x1[i].resize ( nDofx );
                }

                // local storage of global mapping and solution
                for ( unsigned i = 0; i < nDofu; i++ ) {
                    unsigned solDof = msh->GetSolutionDof ( i, iel, solType ); // global to global mapping between solution node and solution dof
                    eigenFunctionCheck[i] = ( *sol->_Sol[eigfIndex[i1]] ) ( solDof );
                    eigenFunctionOldCheck[i] = ( *sol->_Sol[eigfIndex[j1]] ) ( solDof );
                }

                // local storage of coordinates
                for ( unsigned i = 0; i < nDofx; i++ ) {
                    unsigned xDof  = msh->GetSolutionDof ( i, iel, xType ); // global to global mapping between coordinates node and coordinate dof

                    for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                        x1[jdim][i] = ( *msh->_topology->_Sol[jdim] ) ( xDof ); // global extraction and local storage for the element coordinates
                    }
                }

                double weight;
                vector <double> phi;  // local test function

                // *** Gauss point loop ***
                for ( unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solType]->GetGaussPointNumber(); ig++ ) {
                    // *** get gauss point weight, test function and test function partial derivatives ***
                    msh->_finiteElement[ielGeom][solType]->Jacobian ( x1, ig, weight, phi, phi_x );
                    double eigenFunction_gss = 0.;
                    double eigenFunction_gss_old = 0.;

                    for ( unsigned i = 0; i < nDofu; i++ ) {
                        eigenFunction_gss += phi[i] * eigenFunctionCheck[i];
                        eigenFunction_gss_old += phi[i] * eigenFunctionOldCheck[i];
                    }

                    integral += eigenFunction_gss * eigenFunction_gss_old * weight;
                }
            }

            double globalIntegral = 0.;
            MPI_Allreduce ( &integral, &globalIntegral, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
            std::cout << "i = " << i1 << " , " << "j = " << j1 << " , " << "integral = " << globalIntegral << std::endl;
        }
    }

    //END GRAM SCHMIDT CHECK

    // ***************** END ASSEMBLY *******************
}


void IsElementInBall ( unsigned &check, const std::vector<double> &ballCenter, const double &ballRadius, const std::vector <double> &elementCoordinates )
{

    unsigned dim = ballCenter.size();

    std::vector<double> rescaledCoordinates ( dim, 0. );

    for ( unsigned n = 0; n < dim; n++ ) {
        rescaledCoordinates[n] = elementCoordinates[n] - ballCenter[n]; //rescaling so that the center is the origin
    }




}





