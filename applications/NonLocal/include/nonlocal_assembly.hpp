#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "MultiLevelSolution.hpp"


//THIS IS THE 1D ASSEMBLY FOR THE NONLOCAL INTERFACE PROBLEM

using namespace femus;

// double GetExactSolutionValue(const std::vector < double >& x) {
//   double pi = acos(-1.);
//   return cos(pi * x[0]) * cos(pi * x[1]);
// };
//
// void GetExactSolutionGradient(const std::vector < double >& x, vector < double >& solGrad) {
//   double pi = acos(-1.);
//   solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
//   solGrad[1] = -pi * cos(pi * x[0]) * sin(pi * x[1]);
// };
//
// double GetExactSolutionLaplace ( const std::vector < double >& x )
// {
//     double pi = acos ( -1. );
//     return -pi * pi * cos ( pi * x[0] ) * cos ( pi * x[1] ) - pi * pi * cos ( pi * x[0] ) * cos ( pi * x[1] );
// };

bool nonLocalAssembly = true;
//DELTA sizes: martaTest1: 0.4, martaTest2: 0.01, martaTest3: 0.53, martaTest4: 0.2, maxTest1: both 0.4, maxTest2: both 0.01, maxTest3: both 0.53, maxTest4: both 0.2, maxTest5: both 0.1, maxTest6: both 0.8,  maxTest7: both 0.05, maxTest8: both 0.025, maxTest9: both 0.0125, maxTest10: both 0.00625
double delta1 = 0.1; //DELTA SIZES (w 2 refinements): interface: delta1 = 0.4, delta2 = 0.2, nonlocal_boundary_test.neu: 0.0625 * 4
double delta2 = 0.1;
double epsilon = ( delta1 > delta2 ) ? delta1 : delta2;
double leftBound = - 1.2;
double rightBound = 1.2;

void GetBoundaryFunctionValue ( double &value, const std::vector < double >& x )
{
//     value = 0.;
//     value = x[0];
//     value = x[0] * x[0];
    value = x[0] * x[0] * x[0];
//     value = x[0] * x[0] * x[0] * x[0] + 0.1 * x[0] * x[0];
//     value = x[0] * x[0] * x[0] * x[0];
//     value =  2 * x[0] + x[0] * x[0] * x[0] * x[0] * x[0]; //this is 2x + x^5


}

void ReorderElement ( std::vector < int > &dofs, std::vector < double > & sol, std::vector < std::vector < double > > & x );

void RectangleAndBallRelation ( bool &theyIntersect, const std::vector<double> &ballCenter, const double &ballRadius, const std::vector < std::vector < double> > &elementCoordinates,  std::vector < std::vector < double> > &newCoordinates );

void RectangleAndBallRelation2 ( bool & theyIntersect, const std::vector<double> &ballCenter, const double & ballRadius, const std::vector < std::vector < double> > &elementCoordinates, std::vector < std::vector < double> > &newCoordinates );

const elem_type *fem = new const elem_type_2D ( "quad", "linear", "second" ); //to use a different quadrature rule in the inner integral

void AssembleNonLocalSys ( MultiLevelProblem& ml_prob )
{
    adept::Stack& s = FemusInit::_adeptStack;

    LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ( "NonLocal" );
    const unsigned level = mlPdeSys->GetLevelToAssemble();

    Mesh*                    msh = ml_prob._ml_msh->GetLevel ( level );
    elem*                     el = msh->el;

    MultiLevelSolution*    mlSol = ml_prob._ml_sol;
    Solution*                sol = ml_prob._ml_sol->GetSolutionLevel ( level );

    LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
    SparseMatrix*             KK = pdeSys->_KK;
    NumericVector*           RES = pdeSys->_RES;

    const unsigned  dim = msh->GetDimension();
    unsigned dim2 = ( 3 * ( dim - 1 ) + ! ( dim - 1 ) );
    const unsigned maxSize = static_cast< unsigned > ( ceil ( pow ( 3, dim ) ) ); // conservative: based on line3, quad9, hex27

    unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
    unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

    unsigned soluIndex;
    soluIndex = mlSol->GetIndex ( "u" ); // get the position of "u" in the ml_sol object
    unsigned soluType = mlSol->GetSolutionType ( soluIndex ); // get the finite element type for "u"

    unsigned soluPdeIndex;
    soluPdeIndex = mlPdeSys->GetSolPdeIndex ( "u" ); // get the position of "u" in the pdeSys object

    vector < adept::adouble >  solu; // local solution for the local assembly (it uses adept)
    solu.reserve ( maxSize );

    vector < double >  solu1; // local solution for the nonlocal assembly
    vector < double >  solu2; // local solution for the nonlocal assembly
    solu1.reserve ( maxSize );
    solu2.reserve ( maxSize );

    unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

    vector < vector < double > > x1 ( dim );
    vector < vector < double > > x2 ( dim );

    for ( unsigned k = 0; k < dim; k++ ) {
        x1[k].reserve ( maxSize );
        x2[k].reserve ( maxSize );
    }

    vector <double> phi;  // local test function
    vector <double> phi_x; // local test function first order partial derivatives
    vector <double> phi_xx; // local test function second order partial derivatives
    double weight; // gauss point weight

    phi.reserve ( maxSize );
    phi_x.reserve ( maxSize * dim );
    phi_xx.reserve ( maxSize * dim2 );

    vector< adept::adouble > aRes; // local redidual vector
    aRes.reserve ( maxSize );

    vector< int > l2GMap1; // local to global mapping
    vector< int > l2GMap2; // local to global mapping
    l2GMap1.reserve ( maxSize );
    l2GMap2.reserve ( maxSize );

    vector< double > Res1; // local redidual vector
    Res1.reserve ( maxSize );
    vector< double > Res2; // local redidual vector
    Res1.reserve ( maxSize );

    vector < double > Jac11;
    Jac11.reserve ( maxSize * maxSize );
    vector < double > Jac12;
    Jac12.reserve ( maxSize * maxSize );

    vector < double > Jac21;
    Jac21.reserve ( maxSize * maxSize );
    vector < double > Jac22;
    Jac22.reserve ( maxSize * maxSize );


    KK->zero(); // Set to zero all the entries of the Global Matrix

    //BEGIN nonlocal assembly

    //create element groups

    //BEGIN
    unsigned numberOfElements = msh->GetNumberOfElements();

    std::vector < unsigned > elementGroups ( numberOfElements ) ;

    for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

        short unsigned ielGeom = msh->GetElementType ( iel );
        unsigned nDof1  = msh->GetElementDofNumber ( iel, soluType );

        for ( int k = 0; k < dim; k++ ) {
            x1[k].resize ( nDof1 );
        }

        for ( unsigned i = 0; i < nDof1; i++ ) {
            unsigned xDof  = msh->GetSolutionDof ( i, iel, xType );

            for ( unsigned k = 0; k < dim; k++ ) {
                x1[k][i] = ( *msh->_topology->_Sol[k] ) ( xDof );
            }
        }

        double xMax = x1[0][1];
        double xMin = x1[0][0];

        if ( xMax <= leftBound && xMin >= leftBound - delta1 )  elementGroups[iel] = 5;

        else if ( xMax <= 0. && xMin >= leftBound )  elementGroups[iel] = 7;

        else if ( xMax <= rightBound && xMin >= 0. )  elementGroups[iel] = 8;

        else if ( xMax <= rightBound + delta2 && xMin >= rightBound )  elementGroups[iel] = 6;




    }

    for ( unsigned iel = 0; iel < elementGroups.size(); iel++ ) {
        unsigned group = 0;
        MPI_Allreduce ( &elementGroups[iel], &group, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );
        elementGroups[iel] = group;
    }

//         for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {
//               for ( int kproc = 0; kproc < nprocs; kproc++ ) {
//             MPI_Bcast ( &elementGroups[iel], 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD );
//         }
//     }

    for ( unsigned iel = 0; iel < elementGroups.size(); iel++ ) {
        std::cout << " iel = " << iel << " , " << "elementGroups[" << iel << "] = " << elementGroups[iel] << std::endl;
    }

    //END

    //loop to change _Bdc in the boundary elements and assign the BoundaryFunctionValue to their nodes
    //BEGIN
    for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

        short unsigned ielGroup = elementGroups[iel];

        if ( ielGroup == 5 || ielGroup == 6 ) { //5 and 6 are the boundary surfaces

            unsigned nDofu  = msh->GetElementDofNumber ( iel, soluType );
            std::vector <double> dofCoordinates ( dim );

            for ( unsigned i = 0; i < nDofu; i++ ) {
                unsigned solDof = msh->GetSolutionDof ( i, iel, soluType );
                unsigned xDof = msh->GetSolutionDof ( i, iel, xType );
                sol->_Bdc[soluIndex]->set ( solDof, 0. );

                for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                    dofCoordinates[jdim] = ( *msh->_topology->_Sol[jdim] ) ( xDof );
                }

                double bdFunctionValue;
                GetBoundaryFunctionValue ( bdFunctionValue, dofCoordinates );
                sol->_Sol[soluIndex]->set ( solDof, bdFunctionValue );

            }

        }

    }

    sol->_Bdc[soluIndex]->close();
    sol->_Sol[soluIndex]->close();
    //END

    for ( int kproc = 0; kproc < nprocs; kproc++ ) {
        for ( int jel = msh->_elementOffset[kproc]; jel < msh->_elementOffset[kproc + 1]; jel++ ) {

            short unsigned jelGeom;
            short unsigned jelGroup = elementGroups[jel];
            unsigned nDof2;
            //unsigned nDofx2;

            if ( iproc == kproc ) {
                jelGeom = msh->GetElementType ( jel );
                nDof2  = msh->GetElementDofNumber ( jel, soluType );
            }

            MPI_Bcast ( &jelGeom, 1, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD );
            MPI_Bcast ( &nDof2, 1, MPI_UNSIGNED, kproc, MPI_COMM_WORLD );

            l2GMap2.resize ( nDof2 );
            solu2.resize ( nDof2 );

            for ( int k = 0; k < dim; k++ ) {
                x2[k].resize ( nDof2 );
            }

            if ( iproc == kproc ) {
                for ( unsigned j = 0; j < nDof2; j++ ) {
                    l2GMap2[j] = pdeSys->GetSystemDof ( soluIndex, soluPdeIndex, j, jel );
                    unsigned solDof = msh->GetSolutionDof ( j, jel, soluType );
                    solu2[j] = ( *sol->_Sol[soluIndex] ) ( solDof );
                    unsigned xDof  = msh->GetSolutionDof ( j, jel, xType );

                    for ( unsigned k = 0; k < dim; k++ ) {
                        x2[k][j] = ( *msh->_topology->_Sol[k] ) ( xDof );
                    }
                }

//                 ReorderElement ( l2GMap2, solu2, x2 );
            }

            MPI_Bcast ( &l2GMap2[0], nDof2, MPI_UNSIGNED, kproc, MPI_COMM_WORLD );
            MPI_Bcast ( &solu2[0], nDof2, MPI_DOUBLE, kproc, MPI_COMM_WORLD );

            for ( unsigned k = 0; k < dim; k++ ) {
                MPI_Bcast ( & x2[k][0], nDof2, MPI_DOUBLE, kproc, MPI_COMM_WORLD );
            }

            for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

                short unsigned ielGeom = msh->GetElementType ( iel );
                short unsigned ielGroup = elementGroups[iel];
                unsigned nDof1  = msh->GetElementDofNumber ( iel, soluType );

                l2GMap1.resize ( nDof1 );
                solu1.resize ( nDof1 );

                Jac11.assign ( nDof1 * nDof1, 0. );
                Jac12.assign ( nDof1 * nDof2, 0. );
                Jac21.assign ( nDof2 * nDof1, 0. );
                Jac22.assign ( nDof2 * nDof2, 0. );
                Res1.assign ( nDof1, 0. );
                Res2.assign ( nDof2, 0. );

                for ( int k = 0; k < dim; k++ ) {
                    x1[k].resize ( nDof1 );
                }

                for ( unsigned i = 0; i < nDof1; i++ ) {
                    l2GMap1[i] = pdeSys->GetSystemDof ( soluIndex, soluPdeIndex, i, iel );
                    unsigned solDof = msh->GetSolutionDof ( i, iel, soluType );
                    solu1[i] = ( *sol->_Sol[soluIndex] ) ( solDof );
                    unsigned xDof  = msh->GetSolutionDof ( i, iel, xType );

                    for ( unsigned k = 0; k < dim; k++ ) {
                        x1[k][i] = ( *msh->_topology->_Sol[k] ) ( xDof );
                    }
                }

//                 ReorderElement ( l2GMap1, solu1, x1 );

                unsigned igNumber = msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber();
                vector < vector < double > > xg1 ( igNumber );
                vector <double> weight1 ( igNumber );
                vector < vector <double> > phi1x ( igNumber );

                for ( unsigned ig = 0; ig < igNumber; ig++ ) {
                    msh->_finiteElement[ielGeom][soluType]->Jacobian ( x1, ig, weight1[ig], phi1x[ig], phi_x );

                    xg1[ig].assign ( dim, 0. );

                    for ( unsigned i = 0; i < nDof1; i++ ) {
                        for ( unsigned k = 0; k < dim; k++ ) {
                            xg1[ig][k] += x1[k][i] * phi1x[ig][i];
                        }
                    }
                }

                double radius;

                if ( ( ielGroup == 5 || ielGroup == 7 ) && ( jelGroup == 5 || jelGroup == 7 ) ) radius = delta1; //both x and y are in Omega_1

                else if ( ( ielGroup == 5 || ielGroup == 7 ) && ( jelGroup == 6 || jelGroup == 8 ) ) radius = epsilon; // x is in Omega_1 and y is in Omega_2

                else if ( ( ielGroup == 6 || ielGroup == 8 ) && ( jelGroup == 5 || jelGroup == 7 ) ) radius = epsilon; // x is in Omega_2 and y is in Omega_1

                else if ( ( ielGroup == 6 || ielGroup == 8 ) && ( jelGroup == 6 || jelGroup == 8 ) ) radius = delta2; // both x and y are in Omega_2

                bool ifAnyIntersection = false;

                for ( unsigned ig = 0; ig < igNumber; ig++ ) {

                    if ( iel == jel ) {
                        for ( unsigned i = 0; i < nDof1; i++ ) {
//                                 Res1[i] -= 0. * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = 0)
//                             Res1[i] -=  - 2. * weight1[ig]  * phi1x[ig][i]; //Ax - f (so f = - 2)
                            Res1[i] -=  - 6. * xg1[ig][0] * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 6 x )
//                                 Res1[i] -= ( - 12. * xg1[ig][0] * xg1[ig][0] - 6. / 5. * radius * radius - 2. * radius ) * weight1[ig] * phi1x[ig][i];  //Ax - f (so f = - 12x^2 - 6/5 * delta^2 - 2 delta)
//                                      Res1[i] -=  - 20. * ( xg1[ig][0] * xg1[ig][0] * xg1[ig][0] ) * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 20 x^3 )
//                                 Res1[i] -=  - 12. * ( xg1[ig][0] * xg1[ig][0] ) * weight1[ig] * phi1x[ig][i]; //Ax - f (so f = - 12 x^2 )
                        }
                    }

                    std::vector< std::vector < double > > x2New;
                    bool theyIntersect;
                    RectangleAndBallRelation2 ( theyIntersect, xg1[ig], radius, x2, x2New );

                    if ( theyIntersect ) {

                        ifAnyIntersection = true;

                        unsigned jgNumber = msh->_finiteElement[jelGeom][soluType]->GetGaussPointNumber();

//                             unsigned jgNumber = fem->GetGaussPointNumber();

                        for ( unsigned jg = 0; jg < jgNumber; jg++ ) {

                            vector <double>  phi2y;
                            double weight2;

                            msh->_finiteElement[jelGeom][soluType]->Jacobian ( x2New, jg, weight2, phi2y, phi_x );
//                                 fem->Jacobian ( x2New, jg, weight2, phi2y, phi_x );

                            std::vector< double > xg2 ( dim, 0. );

                            for ( unsigned j = 0; j < nDof2; j++ ) {
                                for ( unsigned k = 0; k < dim; k++ ) {
                                    xg2[k] += x2New[k][j] * phi2y[j];
                                }
                            }

                            std::vector <double> xg2Local ( dim );

                            for ( unsigned k = 0; k < dim; k++ ) {
                                xg2Local[k] = - 1. + 2. * ( xg2[k] - x2[k][k] ) / ( x2[k][k + 1] - x2[k][k] );
                            }

                            double weightTemp;
                            msh->_finiteElement[jelGeom][soluType]->Jacobian ( x2, xg2Local, weightTemp, phi2y, phi_x );
//                                 fem->Jacobian ( x2, xg2Local, weightTemp, phi2y, phi_x );

                            for ( unsigned i = 0; i < nDof1; i++ ) {
                                for ( unsigned j = 0; j < nDof1; j++ ) {
                                    double jacValue11 = weight1[ig] * weight2 * 3. / 2. * ( 1. / pow ( radius, 3. ) ) * ( phi1x[ig][i] ) * phi1x[ig][j];
                                    Jac11[i * nDof1 + j] -= jacValue11;
                                    Res1[i] +=  jacValue11 * solu1[j];
                                }

                                for ( unsigned j = 0; j < nDof2; j++ ) {
                                    double jacValue12 = - weight1[ig] * weight2 * 3. / 2. * ( 1. / pow ( radius, 3. ) ) * ( phi1x[ig][i] ) * phi2y[j];
                                    Jac12[i * nDof2 + j] -= jacValue12;
                                    Res1[i] +=  jacValue12 * solu2[j];
                                }//endl j loop
                            }

                            for ( unsigned i = 0; i < nDof2; i++ ) {
                                for ( unsigned j = 0; j < nDof1; j++ ) {
                                    double jacValue21 = weight1[ig] * weight2 * 3. / 2. * ( 1. / pow ( radius, 3. ) ) * ( - phi2y[i] ) * phi1x[ig][j];
                                    Jac21[i * nDof1 + j] -= jacValue21;
                                    Res2[i] +=  jacValue21 * solu1[j];
                                }

                                for ( unsigned j = 0; j < nDof2; j++ ) {
                                    double jacValue22 = - weight1[ig] * weight2 * 3. / 2. * ( 1. / pow ( radius, 3. ) ) * ( - phi2y[i] ) * phi2y[j];
                                    Jac22[i * nDof2 + j] -= jacValue22;
                                    Res2[i] +=  jacValue22 * solu2[j];
                                }//endl j loop
                            } //endl i loop
                        }//end jg loop
                    }
                }//end ig loop

                if ( ifAnyIntersection ) {
                    KK->add_matrix_blocked ( Jac11, l2GMap1, l2GMap1 );
                    KK->add_matrix_blocked ( Jac12, l2GMap1, l2GMap2 );
                    RES->add_vector_blocked ( Res1, l2GMap1 );

                    KK->add_matrix_blocked ( Jac21, l2GMap2, l2GMap1 );
                    KK->add_matrix_blocked ( Jac22, l2GMap2, l2GMap2 );
                    RES->add_vector_blocked ( Res2, l2GMap2 );
                }
            } //end iel loop
        } // end jel loop
    } //end kproc loop

    RES->close();

    KK->close();

//     Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
//     MatAssemblyBegin ( A, MAT_FINAL_ASSEMBLY );
//     MatAssemblyEnd ( A, MAT_FINAL_ASSEMBLY );
//     PetscViewer viewer;
//     MatView ( A, viewer );

//     Vec v = ( static_cast< PetscVector* > ( RES ) )->vec();
//     VecView(v,PETSC_VIEWER_STDOUT_WORLD);

    // ***************** END ASSEMBLY *******************
}




void AssembleLocalSys ( MultiLevelProblem& ml_prob )
{
    adept::Stack& s = FemusInit::_adeptStack;

    LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ( "Local" );
    const unsigned level = mlPdeSys->GetLevelToAssemble();

    Mesh*                    msh = ml_prob._ml_msh->GetLevel ( level );
    elem*                     el = msh->el;

    MultiLevelSolution*    mlSol = ml_prob._ml_sol;
    Solution*                sol = ml_prob._ml_sol->GetSolutionLevel ( level );

    LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
    SparseMatrix*             KK = pdeSys->_KK;
    NumericVector*           RES = pdeSys->_RES;

    const unsigned  dim = msh->GetDimension();
    unsigned dim2 = ( 3 * ( dim - 1 ) + ! ( dim - 1 ) );
    const unsigned maxSize = static_cast< unsigned > ( ceil ( pow ( 3, dim ) ) ); // conservative: based on line3, quad9, hex27

    unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)
    unsigned    nprocs = msh->n_processors(); // get the noumber of processes (for parallel computation)

    unsigned soluIndex;
    soluIndex = mlSol->GetIndex ( "u_local" ); // get the position of "u" in the ml_sol object
    unsigned soluType = mlSol->GetSolutionType ( soluIndex ); // get the finite element type for "u"

    unsigned soluPdeIndex;
    soluPdeIndex = mlPdeSys->GetSolPdeIndex ( "u_local" ); // get the position of "u" in the pdeSys object

    vector < adept::adouble >  solu; // local solution for the local assembly (it uses adept)
    solu.reserve ( maxSize );

    unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

    vector < vector < double > > x1 ( dim );

    for ( unsigned k = 0; k < dim; k++ ) {
        x1[k].reserve ( maxSize );
    }

    vector <double> phi;  // local test function
    vector <double> phi_x; // local test function first order partial derivatives
    vector <double> phi_xx; // local test function second order partial derivatives
    double weight; // gauss point weight

    phi.reserve ( maxSize );
    phi_x.reserve ( maxSize * dim );
    phi_xx.reserve ( maxSize * dim2 );

    vector< adept::adouble > aRes; // local redidual vector
    aRes.reserve ( maxSize );

    vector< int > l2GMap1; // local to global mapping
    l2GMap1.reserve ( maxSize );

    vector< double > Res1; // local redidual vector
    Res1.reserve ( maxSize );

    vector < double > Jac11;
    Jac11.reserve ( maxSize * maxSize );

    KK->zero(); // Set to zero all the entries of the Global Matrix

    //BEGIN local assembly

    //BEGIN
    unsigned numberOfElements = msh->GetNumberOfElements();

    std::vector < unsigned > elementGroups ( numberOfElements ) ;

    for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

        short unsigned ielGeom = msh->GetElementType ( iel );
        unsigned nDof1  = msh->GetElementDofNumber ( iel, soluType );

        for ( int k = 0; k < dim; k++ ) {
            x1[k].resize ( nDof1 );
        }

        for ( unsigned i = 0; i < nDof1; i++ ) {
            unsigned xDof  = msh->GetSolutionDof ( i, iel, xType );

            for ( unsigned k = 0; k < dim; k++ ) {
                x1[k][i] = ( *msh->_topology->_Sol[k] ) ( xDof );
            }
        }

        double xMax = x1[0][1];
        double xMin = x1[0][0];

        if ( xMax <= leftBound && xMin >= leftBound - delta1 )  elementGroups[iel] = 5;

        else if ( xMax <= 0. && xMin > leftBound )  elementGroups[iel] = 6;

        else if ( xMax < rightBound && xMin >= 0. )  elementGroups[iel] = 7;

        else if ( xMax <= rightBound + delta2 && xMin >= rightBound )  elementGroups[iel] = 8;

    }

    for ( int kproc = 0; kproc < nprocs; kproc++ ) {
        unsigned elementsInIproc = msh->_elementOffset[iproc + 1] - msh->_elementOffset[iproc];
        MPI_Bcast ( &elementGroups[iproc], elementsInIproc, MPI_UNSIGNED_SHORT, kproc, MPI_COMM_WORLD );
    }

    //END

    for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

        short unsigned ielGroup = elementGroups[iel];

        if ( ielGroup == 5 || ielGroup == 6 ) { //5 and 6 are the boundary surfaces

            unsigned nDofu  = msh->GetElementDofNumber ( iel, soluType );
            std::vector <double> dofCoordinates ( dim );

            for ( unsigned i = 0; i < nDofu; i++ ) {
                unsigned solDof = msh->GetSolutionDof ( i, iel, soluType );
                unsigned xDof = msh->GetSolutionDof ( i, iel, xType );
                sol->_Bdc[soluIndex]->set ( solDof, 0. );

                for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                    dofCoordinates[jdim] = ( *msh->_topology->_Sol[jdim] ) ( xDof );
                }

                double bdFunctionValue;
                GetBoundaryFunctionValue ( bdFunctionValue, dofCoordinates );
                sol->_Sol[soluIndex]->set ( solDof, bdFunctionValue );

            }

        }

    }

    sol->_Bdc[soluIndex]->close();
    sol->_Sol[soluIndex]->close();


    // element loop: each process loops only on the elements that owns
    for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

        short unsigned ielGeom = msh->GetElementType ( iel );
        unsigned nDofu  = msh->GetElementDofNumber ( iel, soluType ); // number of solution element dofs
        unsigned nDofx = msh->GetElementDofNumber ( iel, xType ); // number of coordinate element dofs

        // resize local arrays
        l2GMap1.resize ( nDofu );
        solu.resize ( nDofu );

        for ( int i = 0; i < dim; i++ ) {
            x1[i].resize ( nDofx );
        }

        aRes.resize ( nDofu ); //resize
        std::fill ( aRes.begin(), aRes.end(), 0 ); //set aRes to zero

        // local storage of global mapping and solution
        for ( unsigned i = 0; i < nDofu; i++ ) {
            unsigned solDof = msh->GetSolutionDof ( i, iel, soluType ); // global to global mapping between solution node and solution dof
            solu[i] = ( *sol->_Sol[soluIndex] ) ( solDof ); // global extraction and local storage for the solution
            l2GMap1[i] = pdeSys->GetSystemDof ( soluIndex, soluPdeIndex, i, iel ); // global to global mapping between solution node and pdeSys dof
        }

        // local storage of coordinates
        for ( unsigned i = 0; i < nDofx; i++ ) {
            unsigned xDof  = msh->GetSolutionDof ( i, iel, xType ); // global to global mapping between coordinates node and coordinate dof

            for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                x1[jdim][i] = ( *msh->_topology->_Sol[jdim] ) ( xDof ); // global extraction and local storage for the element coordinates
            }
        }

        // start a new recording of all the operations involving adept::adouble variables
        s.new_recording();

        // *** Gauss point loop ***
        for ( unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++ ) {
            // *** get gauss point weight, test function and test function partial derivatives ***
            msh->_finiteElement[ielGeom][soluType]->Jacobian ( x1, ig, weight, phi, phi_x, boost::none );

            // evaluate the solution, the solution derivatives and the coordinates in the gauss point

            vector < adept::adouble > gradSolu_gss ( dim, 0. );
            vector < double > x_gss ( dim, 0. );

            for ( unsigned i = 0; i < nDofu; i++ ) {
                for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                    gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
                    x_gss[jdim] += x1[jdim][i] * phi[i];
                }
            }


            double aCoeff = 1.;

            // *** phi_i loop ***
            for ( unsigned i = 0; i < nDofu; i++ ) {

                adept::adouble laplace = 0.;

                for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                    laplace   +=  aCoeff * phi_x[i * dim + jdim] * gradSolu_gss[jdim];
                }

                double srcTerm =  6. * x_gss[0] ; // so f = - 6 x
//                 double srcTerm =  12. * x_gss[0] * x_gss[0] ; // so f = - 12 x^2
//                 double srcTerm =  2. ; // so f = - 2
//                 double srcTerm =  - 1. ; // so f = 1
//                 double srcTerm =  0./*- GetExactSolutionLaplace(x_gss)*/ ;
                aRes[i] += ( srcTerm * phi[i] + laplace ) * weight;

            } // end phi_i loop
        } // end gauss point loop

        //--------------------------------------------------------------------------------------------------------
        // Add the local Matrix/Vector into the global Matrix/Vector

        //copy the value of the adept::adoube aRes in double Res and store
        Res1.resize ( nDofu ); //resize

        for ( int i = 0; i < nDofu; i++ ) {
            Res1[i] = - aRes[i].value();
        }

        RES->add_vector_blocked ( Res1, l2GMap1 );

        // define the dependent variables
        s.dependent ( &aRes[0], nDofu );

        // define the independent variables
        s.independent ( &solu[0], nDofu );

        // get the jacobian matrix (ordered by row major )
        Jac11.resize ( nDofu * nDofu ); //resize
        s.jacobian ( &Jac11[0], true );

        //store K in the global matrix KK
        KK->add_matrix_blocked ( Jac11, l2GMap1, l2GMap1 );

        s.clear_independents();
        s.clear_dependents();

    } //end element loop for each process

    //END local assembly

    RES->close();

    KK->close();

//     Mat A = ( static_cast<PetscMatrix*> ( KK ) )->mat();
//     MatAssemblyBegin ( A, MAT_FINAL_ASSEMBLY );
//     MatAssemblyEnd ( A, MAT_FINAL_ASSEMBLY );
//     PetscViewer viewer;
//     MatView ( A, viewer );

//     Vec v = ( static_cast< PetscVector* > ( RES ) )->vec();
//     VecView(v,PETSC_VIEWER_STDOUT_WORLD);

    // ***************** END ASSEMBLY *******************
}








void RectangleAndBallRelation ( bool & theyIntersect, const std::vector<double> &ballCenter, const double & ballRadius, const std::vector < std::vector < double> > &elementCoordinates, std::vector < std::vector < double> > &newCoordinates )
{

    //elementCoordinates are the coordinates of the vertices of the element

    theyIntersect = false; //by default we assume the two sets are disjoint

    unsigned dim = 1;
    unsigned nDofs = elementCoordinates[0].size();

    std::vector< std::vector < double > > ballVerticesCoordinates ( dim );
    newCoordinates.resize ( dim );


    for ( unsigned n = 0; n < dim; n++ ) {
        newCoordinates[n].resize ( nDofs );
        ballVerticesCoordinates[n].resize ( 2 );

        for ( unsigned i = 0; i < nDofs; i++ ) {
            newCoordinates[n][i] = elementCoordinates[n][i]; //this is just an initalization, it will be overwritten
        }
    }

    double xMinElem = elementCoordinates[0][0];
    double xMaxElem = elementCoordinates[0][1];


    for ( unsigned i = 0; i < 2; i++ ) {
        if ( elementCoordinates[0][i] < xMinElem ) xMinElem = elementCoordinates[0][i];

        if ( elementCoordinates[0][i] > xMaxElem ) xMaxElem = elementCoordinates[0][i];
    }

    ballVerticesCoordinates[0][0] =  ballCenter[0] - ballRadius;

    ballVerticesCoordinates[0][1] = ballCenter[0] + ballRadius;

    newCoordinates[0][0] = ( ballVerticesCoordinates[0][0] >= xMinElem ) ? ballVerticesCoordinates[0][0] : xMinElem;

    newCoordinates[0][1] = ( ballVerticesCoordinates[0][1] >= xMaxElem ) ? xMaxElem : ballVerticesCoordinates[0][1];

    if ( newCoordinates[0][0] < newCoordinates[0][1] ) { //ball and rectangle intersect

        theyIntersect = true;

        if ( nDofs == 3 ) newCoordinates[0][2] = 0.5 * ( newCoordinates[0][0] + newCoordinates[0][1] );

    }

}

const unsigned swap[4][9] = {
    {0, 1, 2, 3, 4, 5, 6, 7, 8},
    {3, 0, 1, 2, 7, 4, 5, 6, 8},
    {2, 3, 0, 1, 6, 7, 4, 5, 8},
    {1, 2, 3, 0, 5, 6, 7, 4, 8}
};

void ReorderElement ( std::vector < int > &dofs, std::vector < double > & sol, std::vector < std::vector < double > > & x )
{

    unsigned type = 0;

    if ( fabs ( x[0][0] - x[0][1] ) > 1.e-10 ) {
        if ( x[0][0] - x[0][1] > 0 ) {
            type = 2;
        }
    }

    else {
        type = 1;

        if ( x[1][0] - x[1][1] > 0 ) {
            type = 3;
        }
    }

    if ( type != 0 ) {
        std::vector < int > dofsCopy = dofs;
        std::vector < double > solCopy = sol;
        std::vector < std::vector < double > > xCopy = x;

        for ( unsigned i = 0; i < dofs.size(); i++ ) {
            dofs[i] = dofsCopy[swap[type][i]];
            sol[i] = solCopy[swap[type][i]];

            for ( unsigned k = 0; k < x.size(); k++ ) {
                x[k][i] = xCopy[k][swap[type][i]];
            }
        }
    }
}


void RectangleAndBallRelation2 ( bool & theyIntersect, const std::vector<double> &ballCenter, const double & ballRadius, const std::vector < std::vector < double> > &elementCoordinates, std::vector < std::vector < double> > &newCoordinates )
{

    theyIntersect = false; //by default we assume the two sets are disjoint

    unsigned dim = 1;
    unsigned nDofs = elementCoordinates[0].size();

    newCoordinates.resize ( dim );

    for ( unsigned i = 0; i < dim; i++ ) {
        newCoordinates[i].resize ( nDofs );
    }

    double xMin = elementCoordinates[0][0];
    double xMax = elementCoordinates[0][1];

    if ( xMin > xMax ) {
        std::cout << "error, the nodes are not ordered in the right order" << std::endl;

        for ( unsigned i = 0; i < nDofs; i++ ) {
            std::cout <<  elementCoordinates[0][i]  << std::endl;
        }

        exit ( 0 );
    }


    double xMinBall = ballCenter[0] - ballRadius;
    double xMaxBall = ballCenter[0] + ballRadius;


    xMin = ( xMin > xMinBall ) ? xMin : xMinBall;
    xMax = ( xMax < xMaxBall ) ? xMax : xMaxBall;

    if ( xMin < xMax ) { //ball and rectangle intersect

        theyIntersect = true;

        newCoordinates[0][0] = xMin;
        newCoordinates[0][1] = xMax;

        if ( nDofs == 3 )  newCoordinates[0][2] = 0.5 * ( xMin + xMax );

    }

}

