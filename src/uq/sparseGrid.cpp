#include <math.h>
#include "sparseGrid.hpp"
#include <iostream>
#include <cstdlib>
#include <numeric>
#include <vector>

namespace femus
{

    sparseGrid::sparseGrid ( std::vector < std::vector < double > >  &samples, const bool & output )
    {
        _output = output;

        _N = samples[0].size();
        _M = samples.size();
        _L = 2 /*static_cast<unsigned> ( log10 ( _M )  )*/; //NOTE this might change if we see fit

        _intervals.resize ( _N );
        _hs.resize ( _N );
        _nodes.resize ( _N );
        _hierarchicalDofs.resize ( _N );

        for ( unsigned n = 0; n < _N; n++ ) {
            _intervals[n].resize ( 2 );
            _hs[n].resize ( _L );
            _nodes[n].resize ( _L );
            _hierarchicalDofs[n].resize ( _L );

            for ( unsigned l = 0; l < _L; l++ ) {
                unsigned dofsFullGrid = static_cast<unsigned> ( pow ( 2, l + 1 ) - 1 );
                _nodes[n][l].resize ( dofsFullGrid );

                unsigned   dofsHierarchical = static_cast<unsigned> ( ( ( pow ( 2, l + 1 ) - 1 ) + 1. ) * 0.5 );

                _hierarchicalDofs[n][l].resize ( dofsHierarchical );
//                 _hierarchicalDofs[n][l].resize ( dofsFullGrid ); //this is to use the method with standard FEM

                if ( _output ) std::cout << "dofsFullGrid = " << dofsFullGrid << " , " << "dofsHierarchical = " << dofsHierarchical << std::endl;

            }
        }

        double an;
        double bn;

        //BEGIN figuring out the one dimensional intervals, mesh sizes and node coordinates
        for ( unsigned n = 0; n < _N; n++ ) {
            an = - 1.; //temporary left extremum of the n-th interval [an,bn]
            bn = 1.;   //temporary right extremum of the n-th interval [an,bn]

            for ( unsigned m = 0; m < _M; m++ ) {
                if ( samples[m][n] <= an ) {
                    an = samples[m][n];
                }

                else if ( samples[m][n] >= bn ) {
                    bn = samples[m][n];
                }
            }

            _intervals[n][0] = an - 0.5; //the 0.5 is just to add some tolerance
            _intervals[n][1] = bn + 0.5;

            if ( _output )  std::cout << "a = " << _intervals[n][0] << " , " << "b = " << _intervals[n][1] << std::endl;

            for ( unsigned l = 0; l < _L; l++ ) {
                _hs[n][l] = fabs ( _intervals[n][1] - _intervals[n][0] ) / pow ( 2, l + 1 ); // mesh size of the n-th interval

                if ( _output )   std::cout << "h[" << n << "][" << l << "]= " << _hs[n][l] << std::endl;

                unsigned hierarchicalCounter = 0;

                for ( unsigned i = 0; i < _nodes[n][l].size(); i++ ) {
                    _nodes[n][l][i] = _intervals[n][0] + ( i + 1 ) * _hs[n][l]; //coordinates of the nodes of the n-th interval

                    if ( _output )    std::cout << "node[" << n << "][" << l << "][" << i << "]= " << _nodes[n][l][i] << std::endl;

                    if ( i % 2 == 0 ) {
                        _hierarchicalDofs[n][l][hierarchicalCounter] = i;
                        hierarchicalCounter++;
                    }
                    
//                     _hierarchicalDofs[n][l][i] = i; //this is to use the method with standard FEM
                }
            }

        }

        if ( _output ) {
            for ( unsigned n = 0; n < _N; n++ ) {
                for ( unsigned l = 0; l < _L; l++ ) {
                    for ( unsigned i = 0; i < _hierarchicalDofs[n][l].size(); i++ ) {
                        std::cout << "_hierarchicalDofs[" << n << "][" << l << "][" << i << "] = " << _hierarchicalDofs[n][l][i] << " " ;
                    }

                    std::cout << std::endl;
                }
            }
        }

        //END


        //BEGIN computation of the sparse grid index set

        //this is done by first computing the tensor product space of the complete space and then eliminating those indices that don't satisfy the requirement: | |_1 <= _L + N - 1 (this condition holds assuming all indices start from 1)


        //here we compute the tensor produc index set for the full grid
        std::vector < std::vector < unsigned > > Tp;

        ComputeTensorProductSet ( Tp, _L, _N );

        unsigned tensorProductDim = Tp.size();

        //here we only consider the indices of Tp that satisfy the requirement
        unsigned indexCounter = 0;

        for ( unsigned i = 0; i < tensorProductDim; i++ ) {

            unsigned sum = 0;

            for ( unsigned j = 0; j < _N; j++ ) {
                sum += ( Tp[i][j] + 1 ); //this is because our indices start from 0 but the requirement assumes they start from 1
            }

            if ( sum <= _L + _N - 1 ) {  
//                if ( sum == _L * _N ) { //this is to use the method with standard FEM

                _indexSetW.resize ( indexCounter + 1 );
                _indexSetW[indexCounter].resize ( _N );

                for ( unsigned j = 0; j < _N; j++ ) {
                    _indexSetW[indexCounter][j] = Tp[i][j];
                }

                indexCounter++;
            }
        }

        if ( _output ) {
            for ( unsigned i = 0; i < _indexSetW.size(); i++ ) {
                for ( unsigned j = 0; j < _N; j++ ) {
                    std::cout << "_indexSetW[" << i << "][" << j << "]= " << _indexSetW[i][j];
                }

                std::cout << std::endl;
            }
        }

        _numberOfWs = _indexSetW.size();


        //END


        //BEGIN construction of dofIdentifier and storing memory for _nodalValuesPDF

        _dofIdentifier.resize ( _numberOfWs );
        _nodalValuesPDF.resize ( _numberOfWs );

        if ( _output )  std::cout << "-------------------------- Number of W sets = " << _numberOfWs <<  "-------------------------- " << std::endl;

        for ( unsigned w = 0; w < _numberOfWs; w++ ) {
            unsigned identifiersOfW = 1;

            for ( unsigned n = 0; n < _N; n++ ) {

                identifiersOfW *= _hierarchicalDofs[n][_indexSetW[w][n]].size();

            }

            _dofIdentifier[w].resize ( identifiersOfW );
            _nodalValuesPDF[w].resize ( identifiersOfW );

            if ( _output ) std::cout << "identifiersOfW = " << identifiersOfW << std::endl;
        }

        //Here we create the dofs for each W
        std::vector<std::vector< std::vector < unsigned > > > dofsW;
        dofsW.resize ( _numberOfWs );

        std::vector < unsigned> counterW ( _numberOfWs, 0 );

        for ( unsigned w = 0; w < _numberOfWs; w++ ) {

            dofsW[w].resize ( _dofIdentifier[w].size() );

            std::vector< std::vector <int> > inputCartesian ( _N );

            for ( unsigned n = 0; n < _N; n++ ) {
                inputCartesian[n].resize ( _hierarchicalDofs[n][_indexSetW[w][n]].size() );

                for ( unsigned i1 = 0; i1 < _hierarchicalDofs[n][_indexSetW[w][n]].size(); i1++ ) {
                    inputCartesian[n][i1] = _hierarchicalDofs[n][_indexSetW[w][n]][i1];
                }
            }

            CartesianProduct ( inputCartesian, dofsW[w] );


            if ( _output ) {
                for ( unsigned i = 0; i < dofsW[w].size(); i++ ) {
                    for ( unsigned n = 0; n < _N; n++ ) {
                        std::cout << "dofsW[" << w << "][" << i << "][" << n << "] = " << dofsW[w][i][n] << " ";
                    }

                    std::cout << std::endl;
                }
            }

        }


        //now we need to store these dofs in _dofIdentifier
        for ( unsigned w = 0; w < _numberOfWs; w++ ) {
            for ( unsigned i = 0; i < _dofIdentifier[w].size(); i++ ) {
                _dofIdentifier[w][i].resize ( _N );

                for ( unsigned n = 0; n < _N; n++ ) {
                    _dofIdentifier[w][i][n].resize ( 3 );
                    _dofIdentifier[w][i][n][0] = n;
                    _dofIdentifier[w][i][n][1] = _indexSetW[w][n];
                    _dofIdentifier[w][i][n][2] = dofsW[w][i][n];
                }
            }
        }

        if ( _output ) {
            for ( unsigned w = 0; w < _numberOfWs; w++ ) {
                std::cout << " ------------------------- w = " << w << " ------------------------- " << std::endl;

                for ( unsigned i = 0; i < _dofIdentifier[w].size(); i++ ) {
                    std::cout << " ------------------------- i = " << i << " ------------------------- " << std::endl;

                    for ( unsigned n = 0; n < _N; n++ ) {
                        std::cout << " ------------------------- dim = " << n << " ------------------------- " << std::endl;

                        for ( unsigned j = 0; j < 3; j++ ) {
                            std::cout << "_dofIdentifier[" << w << "][" << i << "][" << n << "][" << j << "] = " << _dofIdentifier[w][i][n][j] << " " ;
                        }

                        std::cout << std::endl;
                    }
                }
            }
        }

        //END


        //BEGIN construction of _hierarchicalDofsCoordinates


        _hierarchicalDofsCoordinates.resize ( _numberOfWs );

        for ( unsigned w = 0; w < _numberOfWs; w++ ) {

            _hierarchicalDofsCoordinates[w].resize ( _dofIdentifier[w].size() );

            for ( unsigned i = 0; i < _dofIdentifier[w].size(); i++ ) {

                _hierarchicalDofsCoordinates[w][i].resize ( _N );

                for ( unsigned n = 0; n < _N; n++ ) {

                    _hierarchicalDofsCoordinates[w][i][n] = _nodes[_dofIdentifier[w][i][n][0]][_dofIdentifier[w][i][n][1]][_dofIdentifier[w][i][n][2]];

                }
            }
        }

        //END



    }

    void sparseGrid::EvaluateOneDimensionalPhi ( double &phi, const double &x, const unsigned &n, const unsigned &l, const unsigned &i, const bool &scale )
    {

        //i tells you that you are going to compute the phi associated with the i-th node

        double leftBoundOfSupport = _nodes[n][l][i] - _hs[n][l]; // this is xi - h
        double rightBoundOfSupport = _nodes[n][l][i] + _hs[n][l]; // this is xi + h

        if ( x >= leftBoundOfSupport && x <= rightBoundOfSupport ) {

            phi = 1 - fabs ( ( x - _nodes[n][l][i] ) / _hs[n][l] );

            if ( scale == true ) { //this is to be used when building the nodal values of the PDF approximation
                phi /= _hs[n][l];
            }

        }

        else phi = 0.;
    }

    void sparseGrid::EvaluatePhi ( double &phi, const std::vector <double> &x, std::vector < std::vector < unsigned > > identifier, const bool &scale )
    {

        //identifier tells us what one dimensional phis to multiply in order to get the desired phi
        //identifier[n][0] = n (dim), identifier[n][1] = l (level), identifier[n][2] = i (node)

        std::vector < double > oneDimPhi ( _N, 0. );

        phi = 1.;

        for ( unsigned n = 0; n < _N; n++ ) {
            EvaluateOneDimensionalPhi ( oneDimPhi[n], x[n], identifier[n][0], identifier[n][1], identifier[n][2], scale );
            phi *= oneDimPhi[n];
        }

    }

    void sparseGrid::EvaluateNodalValuesPDF ( std::vector < std::vector < double > >  &samples )
    {
        
        for ( unsigned w = 0; w < _numberOfWs; w++ ) {
            for ( unsigned i = 0; i < _nodalValuesPDF[w].size(); i++ ) {
                _nodalValuesPDF[w][i] = 0.;
            }
        }

        for ( unsigned m = 0; m < _M; m++ ) {

            double sumDenom = 1.;
//             double sumDenom = 0.;
//             for ( unsigned w1 = 0; w1 < _numberOfWs; w1++ ) {
//                 for ( unsigned i1 = 0; i1 < _nodalValuesPDF[w1].size(); i1++ ) {
//                     double valuePhiDenom;
//                     EvaluatePhi ( valuePhiDenom, samples[m], _dofIdentifier[w1][i1], false );
//                     sumDenom += valuePhiDenom;
//                 }
//             }
            
            
            for ( unsigned w = 0; w < _numberOfWs; w++ ) {
                for ( unsigned i = 0; i < _nodalValuesPDF[w].size(); i++ ) {
                    double valuePhi;
                    EvaluatePhi ( valuePhi, samples[m], _dofIdentifier[w][i], false );

                    _nodalValuesPDF[w][i] += ( valuePhi / sumDenom ) / _M;
                }
            }


        }

//         //to remove, this is just a check
//         double sumOfNodalValues = 0.;
// 
//         for ( unsigned w = 0; w < _numberOfWs; w++ ) {
//             for ( unsigned i = 0; i < _nodalValuesPDF[w].size(); i++ ) {
//                 sumOfNodalValues += _nodalValuesPDF[w][i];
//             }
//         }
// 
//         std::cout << " sumOfNodalValues =" << sumOfNodalValues << std::endl;


    }

    void sparseGrid::PrintNodalValuesPDF()
    {

        for ( unsigned w = 0; w < _numberOfWs; w++ ) {
            for ( unsigned i = 0; i < _nodalValuesPDF[w].size(); i++ ) {

                for ( unsigned n = 0; n < _N; n++ ) {

                    std::cout << _hierarchicalDofsCoordinates[w][i][n] << " " ;

                }

                std::cout << _nodalValuesPDF[w][i] << std::endl;
            }
        }

    }

    void sparseGrid::EvaluatePDF ( double &pdfValue, std::vector < double >  &x )
    {
//         pdfValue = 0.;

        std::vector<std::vector<double>> phiFncts(_numberOfWs);
        for ( unsigned w = 0; w < _numberOfWs; w++ ) {
            phiFncts[w].resize(_nodalValuesPDF[w].size());
            for ( unsigned i = 0; i < _nodalValuesPDF[w].size(); i++ ) {
                double valuePhi;
                EvaluatePhi ( valuePhi, x, _dofIdentifier[w][i], true );
//                 pdfValue += _nodalValuesPDF[w][i] * valuePhi;
              phiFncts[w][i]  = _nodalValuesPDF[w][i] * valuePhi;
              std::cout<< phiFncts[w][i] << "," ;
            }
        }        
        std::cout<<std::endl;
    }

    void sparseGrid::ComputeTensorProductSet ( std::vector< std::vector <unsigned>> &Tp, const unsigned &T1, const unsigned &T2 )
    {

        unsigned tensorProductDim = pow ( T1, T2 );

        Tp.resize ( tensorProductDim );

        for ( unsigned i = 0; i < tensorProductDim; i++ ) {
            Tp[i].resize ( T2 );
        }

        unsigned index = 0;
        unsigned counters[T2 + 1];
        memset ( counters, 0, sizeof ( counters ) );

        while ( !counters[T2] ) {

            for ( unsigned j = 0; j < T2; j++ ) {
                Tp[index][j] = counters[T2 - 1 - j];

                if ( _output )  std::cout << " Tp[" << index << "][" << j << "]= " << Tp[index][j] ;
            }

            if ( _output )  std::cout << std::endl;

            index++;

            unsigned i;

            for ( i = 0; counters[i] == T1 - 1; i++ ) { // inner loops that are at maxval restart at zero
                counters[i] = 0;
            }

            ++counters[i];  // the innermost loop that isn't yet at maxval, advances by 1
        }

    }

    void sparseGrid::CartesianProduct ( std::vector<std::vector<int> >& inputCartesian, std::vector<std::vector<unsigned> >& dofsWi )
    {
        unsigned counterDofs = 0;

        auto product = [] ( long long a, std::vector<int>& b ) {
            return a * b.size();
        };
        const long long N = accumulate ( inputCartesian.begin(), inputCartesian.end(), 1LL, product );
        std::vector<int> u ( inputCartesian.size() );

        for ( long long n = 0 ; n < N ; ++n ) {
            lldiv_t q { n, 0 };

            for ( long long i = inputCartesian.size() - 1 ; 0 <= i ; --i ) {
                q = div ( q.quot, inputCartesian[i].size() );
                u[i] = inputCartesian[i][q.rem];
            }

            dofsWi.resize ( counterDofs + 1 );
            dofsWi[counterDofs].resize ( u.size() );

            for ( int j = 0; j < u.size(); j++ ) {

                dofsWi[counterDofs][j] = u[j];

            }

            counterDofs++;

        }
    }


}













