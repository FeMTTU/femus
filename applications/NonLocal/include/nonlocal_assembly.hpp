#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

//THIS IS THE ASSEMBLY FOR THE NONLOCAL INTERACTION PROBLEM

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

//BEGIN to remove
unsigned quadratureType = 0;
int numberOfEigPairs = 2; //dimension of the stochastic variable
std::vector < std::pair<double, double> > eigenvalues ( numberOfEigPairs );
double stdDeviationInput = 0.8;  //standard deviation of the normal distribution (it is the same as the standard deviation of the covariance function in GetEigenPair)
double meanInput = 0.;
//END to remove

void GetBoundaryFunctionValue ( double &value, const std::vector < double >& x )
{
    value = 1.;

}

void IsRectangleInBall ( unsigned &check, const std::vector<double> &ballCenter, const double &ballRadius, const std::vector < std::vector < double> > &elementCoordinates,  std::vector < std::vector < double> > &newCoordinates );

void AssembleNonlocalSys ( MultiLevelProblem& ml_prob )
{
    adept::Stack& s = FemusInit::_adeptStack;

    //  extract pointers to the several objects that we are going to use

    LinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<LinearImplicitSystem> ( "UQ" ); // pointer to the linear implicit system named "Poisson"
    const unsigned level = mlPdeSys->GetLevelToAssemble();

    Mesh*                    msh = ml_prob._ml_msh->GetLevel ( level ); // pointer to the mesh (level) object
    elem*                     el = msh->el;  // pointer to the elem object in msh (level)

    MultiLevelSolution*    mlSol = ml_prob._ml_sol;  // pointer to the multilevel solution object
    Solution*                sol = ml_prob._ml_sol->GetSolutionLevel ( level ); // pointer to the solution (level) object

    LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level]; // pointer to the equation (level) object
    SparseMatrix*             KK = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
    NumericVector*           RES = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

    const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
    unsigned dim2 = ( 3 * ( dim - 1 ) + ! ( dim - 1 ) ); // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
    const unsigned maxSize = static_cast< unsigned > ( ceil ( pow ( 3, dim ) ) ); // conservative: based on line3, quad9, hex27

    unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

    //solution variable
    unsigned soluIndex;
    soluIndex = mlSol->GetIndex ( "u" ); // get the position of "u" in the ml_sol object
    unsigned soluType = mlSol->GetSolutionType ( soluIndex ); // get the finite element type for "u"

    unsigned soluPdeIndex;
    soluPdeIndex = mlPdeSys->GetSolPdeIndex ( "u" ); // get the position of "u" in the pdeSys object

    vector < adept::adouble >  solu; // local solution
    solu.reserve ( maxSize );


    vector < vector < double > > x ( dim ); // local coordinates
    unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

    for ( unsigned i = 0; i < dim; i++ ) {
        x[i].reserve ( maxSize );
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

    vector< int > l2GMap; // local to global mapping
    l2GMap.reserve ( maxSize );
    vector< double > Res; // local redidual vector
    Res.reserve ( maxSize );
    vector < double > Jac;
    Jac.reserve ( maxSize * maxSize );

    KK->zero(); // Set to zero all the entries of the Global Matrix

    //loop to change _bdc in the boundary elements and assign the BoundaryFunctionValue to their nodes
    //BEGIN
    for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

        short unsigned ielGroup = msh->GetElementGroup ( iel );

        if ( ielGroup == 5 ) { //5 is the boundary surface in nonlocal_boundary_test.neu

            unsigned nDofu  = msh->GetElementDofNumber ( iel, soluType );
            std::vector <double> dofCoordinates ( dim );

            for ( unsigned i = 0; i < nDofu; i++ ) {
                unsigned solDof = msh->GetSolutionDof ( i, iel, soluType );
                sol->_Bdc[0]->set ( solDof, 0. ); //TODO not sure about _Bdc[0] solution but it seems to work

                for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                    dofCoordinates[jdim] = ( *msh->_topology->_Sol[jdim] ) ( solDof );
                }

                double bdFunctionValue;
                GetBoundaryFunctionValue ( bdFunctionValue, dofCoordinates );
                sol->_Sol[soluIndex]->set ( solDof, bdFunctionValue );

            }

        }

    }

    //END

    // element loop: each process loops only on the elements that owns
    for ( int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++ ) {

        short unsigned ielGeom = msh->GetElementType ( iel );
        unsigned nDofu  = msh->GetElementDofNumber ( iel, soluType ); // number of solution element dofs
        unsigned nDofx = msh->GetElementDofNumber ( iel, xType ); // number of coordinate element dofs

        // resize local arrays
        l2GMap.resize ( nDofu );
        solu.resize ( nDofu );

        for ( int i = 0; i < dim; i++ ) {
            x[i].resize ( nDofx );
        }

        aRes.resize ( nDofu ); //resize
        std::fill ( aRes.begin(), aRes.end(), 0 ); //set aRes to zero

        // local storage of global mapping and solution
        for ( unsigned i = 0; i < nDofu; i++ ) {
            unsigned solDof = msh->GetSolutionDof ( i, iel, soluType ); // global to global mapping between solution node and solution dof
            solu[i] = ( *sol->_Sol[soluIndex] ) ( solDof ); // global extraction and local storage for the solution

            l2GMap[i] = pdeSys->GetSystemDof ( soluIndex, soluPdeIndex, i, iel ); // global to global mapping between solution node and pdeSys dof
        }

        // local storage of coordinates
        for ( unsigned i = 0; i < nDofx; i++ ) {
            unsigned xDof  = msh->GetSolutionDof ( i, iel, xType ); // global to global mapping between coordinates node and coordinate dof

            for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                x[jdim][i] = ( *msh->_topology->_Sol[jdim] ) ( xDof ); // global extraction and local storage for the element coordinates
            }
        }


        // start a new recording of all the operations involving adept::adouble variables
        s.new_recording();

        // *** Gauss point loop ***
        for ( unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++ ) {
            // *** get gauss point weight, test function and test function partial derivatives ***
            msh->_finiteElement[ielGeom][soluType]->Jacobian ( x, ig, weight, phi, phi_x, boost::none );

            // evaluate the solution, the solution derivatives and the coordinates in the gauss point

            vector < adept::adouble > gradSolu_gss ( dim, 0. );
            vector < double > x_gss ( dim, 0. );

            for ( unsigned i = 0; i < nDofu; i++ ) {
                for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                    gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
                    x_gss[jdim] += x[jdim][i] * phi[i];
                }
            }


            //BEGIN testing IsRectangleInBall
            if ( iel == 73 ) {

                std::cout << "--------------------------------------------------------------------" << std::endl;

                std::vector< double > centerCoordinates ( dim );


                for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                    centerCoordinates[jdim] = x[jdim][2] ; //node 98, elem 73 is the center, mesh size is 0.0625 with numberOfUnifRef = 2
                }

                for ( unsigned jdim = 0; jdim < dim; jdim++ ) {

                    std::cout << centerCoordinates[jdim] << " "; // elem 73, node 98 is the center, mesh size is 0.0625 with numberOfUnifRef = 2
                }

                std::cout << std::endl;


//                 double radius = 0.09375; // (1.5 of the mesh size)
//                 double radius = 0.0625; // (the mesh size)
                double radius = 0.03125; // (0.5 of the mesh size)


                for ( int jel = msh->_elementOffset[iproc]; jel < msh->_elementOffset[iproc + 1]; jel++ ) {

                    std::vector< std::vector < double > > newCoordinates;
                    unsigned check;

                    unsigned nDofx2 = msh->GetElementDofNumber ( jel, xType ); // number of coordinate element dofs

                    std::vector< std::vector < double >> x2 ( dim );

                    for ( int ii = 0; ii < dim; ii++ ) {
                        x2[ii].resize ( nDofx2 );
                    }

                    // local storage of coordinates
                    for ( unsigned ii = 0; ii < nDofx2; ii++ ) {
                        unsigned xDof2  = msh->GetSolutionDof ( ii, jel, xType ); // global to global mapping between coordinates node and coordinate dof

                        for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                            x2[jdim][ii] = ( *msh->_topology->_Sol[jdim] ) ( xDof2 ); // global extraction and local storage for the element coordinates
                        }
                    }

                    IsRectangleInBall ( check, centerCoordinates, radius, x2, newCoordinates );

                    if ( check == 0 || check == 1 ) {
                        std::cout << "elem = " << jel << " , " << "check = " << check << std::endl;

                        for ( unsigned ivertex = 0; ivertex < newCoordinates[0].size(); ivertex++ ) {
                            for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                                std::cout << "xOld[" << jdim << "][" << ivertex << "]= " << x2[jdim][ivertex];
                            }

                            std::cout << std::endl;
                        }

                        for ( unsigned ivertex = 0; ivertex < newCoordinates[0].size(); ivertex++ ) {
                            for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                                std::cout << "xNew[" << jdim << "][" << ivertex << "]= " << newCoordinates[jdim][ivertex];
                            }

                            std::cout << std::endl;
                        }

                        std::cout << std::endl;
                    }

                }

                std::cout << "--------------------------------------------------------------------" << std::endl;

            }

            //END testing IsRectangleInBall

            double aCoeff = 1.;

            // *** phi_i loop ***
            for ( unsigned i = 0; i < nDofu; i++ ) {

                adept::adouble laplace = 0.;

                for ( unsigned jdim = 0; jdim < dim; jdim++ ) {
                    laplace   +=  aCoeff * phi_x[i * dim + jdim] * gradSolu_gss[jdim];
                }

                double srcTerm = 1./*- GetExactSolutionLaplace(x_gss)*/ ;
                aRes[i] += ( srcTerm * phi[i] + laplace ) * weight;

            } // end phi_i loop
        } // end gauss point loop

        //--------------------------------------------------------------------------------------------------------
        // Add the local Matrix/Vector into the global Matrix/Vector

        //copy the value of the adept::adoube aRes in double Res and store
        Res.resize ( nDofu ); //resize

        for ( int i = 0; i < nDofu; i++ ) {
            Res[i] = - aRes[i].value();
        }

        RES->add_vector_blocked ( Res, l2GMap );

        // define the dependent variables
        s.dependent ( &aRes[0], nDofu );

        // define the independent variables
        s.independent ( &solu[0], nDofu );

        // get the jacobian matrix (ordered by row major )
        Jac.resize ( nDofu * nDofu ); //resize
        s.jacobian ( &Jac[0], true );

        //store K in the global matrix KK
        KK->add_matrix_blocked ( Jac, l2GMap, l2GMap );

        s.clear_independents();
        s.clear_dependents();

    } //end element loop for each process

    RES->close();

    KK->close();

    // ***************** END ASSEMBLY *******************
}


void IsRectangleInBall ( unsigned &check, const std::vector<double> &ballCenter, const double &ballRadius, const std::vector < std::vector < double> > &elementCoordinates, std::vector < std::vector < double> > &newCoordinates )
{

    //check = 0 : element entirely contained in ball
    //check = 1 : element and ball intersect
    //check = 2 : element and ball are disjoint

    check = 2; //by default we assume the two sets are disjoint

    unsigned dim = 2;

    std::vector< std::vector < double > > rescaledCoordinates ( dim );
    newCoordinates.resize ( dim );


    for ( unsigned n = 0; n < dim; n++ ) {
        rescaledCoordinates[n].resize ( 4 );
        newCoordinates[n].resize ( 4 );

        for ( unsigned i = 0; i < 4; i++ ) {
            rescaledCoordinates[n][i] = elementCoordinates[n][i] - ballCenter[n]; //rescaling so that the center is the origin
            newCoordinates[n][i] = 0.;
        }
    }

    std::vector<unsigned> vertexCheck ( 4, 0. );

    for ( unsigned i = 0; i < 4; i++ ) {
        double lInfinityNorm = 0.;

        for ( unsigned n = 0; n < dim; n++ ) {
            if ( fabs ( rescaledCoordinates[n][i] ) >= lInfinityNorm ) lInfinityNorm = fabs ( rescaledCoordinates[n][i] );
        }

        if ( lInfinityNorm <= ballRadius ) vertexCheck[i] = 1;
    }

    unsigned insideCheck = 0;

    for ( unsigned i = 0; i < 4; i++ ) {
        insideCheck += vertexCheck[i];
    }

    if ( insideCheck == 4 ) check = 0; //the element is entirely contained in the ball

    else { //either the element is outside or they intersect

        //bottom left corner of ball (south west)
        double xBallSW =  - ballRadius;
        double yBallSW =  - ballRadius;

        //top right corner of ball (north east)
        double xBallNE = ballRadius;
        double yBallNE = ballRadius;

        //bottom left corner of rectangle (south west)
        double xRecSW = rescaledCoordinates[0][0]; //ordering according to Elem.hpp, 0 is SW node
        double yRecSW = rescaledCoordinates[1][0];

        //top right corner of rectangle (north east)
        double xRecNE = rescaledCoordinates[0][2]; //ordering according to Elem.hpp, 2 is NE node
        double yRecNE = rescaledCoordinates[1][2];

        newCoordinates[0][0] = ( xBallSW >= xRecSW ) ? xBallSW : xRecSW;
        newCoordinates[1][0] = ( yBallSW >= yRecSW ) ? yBallSW : yRecSW;

        newCoordinates[0][2] = ( xBallNE >= xRecNE ) ? xRecNE : xBallNE;
        newCoordinates[1][2] = ( yBallNE >= yRecNE ) ? yRecNE : yBallNE;

        if ( newCoordinates[0][0] < newCoordinates[0][2] && newCoordinates[1][0] < newCoordinates[1][2] ) { //ball and rectangle intersect

            check = 1;

            newCoordinates[0][1] = newCoordinates[0][2];
            newCoordinates[1][1] = newCoordinates[1][0];

            newCoordinates[0][3] = newCoordinates[0][0];
            newCoordinates[1][3] = newCoordinates[1][2];

            //translate back
            for ( unsigned n = 0; n < dim; n++ ) {
                for ( unsigned i = 0; i < 4; i++ ) {
                    newCoordinates[n][i] += ballCenter[n] ;
                }
            }

        }

    }

}


