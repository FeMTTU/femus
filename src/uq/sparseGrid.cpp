#include <math.h>
#include "sparseGrid.hpp"
#include <iostream>

namespace femus
{

    sparseGrid::sparseGrid (const unsigned &N, const unsigned &M, std::vector < std::vector < double > >  &samples)
    {

        _N = N;
        _M = M;
        _L = static_cast<unsigned> (log10 (M) + 1);

        _intervals.resize (_N);
        _hs.resize (_N);
        _nodes.resize (_N);
        _hierarchicalDofs.resize (_N);

        for (unsigned n = 0; n < _N; n++) {
            _intervals[n].resize (2);
            _hs[n].resize (_L);
            _nodes[n].resize (_L);
            _hierarchicalDofs[n].resize (_L);
            for (unsigned l = 0; l < _L; l++) {
                unsigned dofsFullGrid = static_cast<unsigned> (pow (2, l + 1) - 1);
                _nodes[n][l].resize (dofsFullGrid);
                unsigned dofsHierarchical;

                if (dofsFullGrid % 2 != 0) { //odd number
                    dofsHierarchical = static_cast<unsigned> ( ( (pow (2, l + 1) - 1) + 1.) * 0.5);
                }
                else { //even number
                    dofsHierarchical = static_cast<unsigned> ( ( (pow (2, l + 1) - 1)) * 0.5);
                }

                _hierarchicalDofs[n][l].resize (dofsHierarchical);

                std::cout << "dofsFullGrid = " << dofsFullGrid << " , " << "dofsHierarchical = " << dofsHierarchical << std::endl;

            }
        }

        double an;
        double bn;

        //BEGIN figuring out the one dimensional intervals, mesh sizes and node coordinates
        for (unsigned n = 0; n < _N; n++) {
            an = - 1.; //temporary left extremum of the n-th interval [an,bn]
            bn = 1.;   //temporary right extremum of the n-th interval [an,bn]
            for (unsigned m = 0; m < _M; m++) {
                if (samples[m][n] <= an) {
                    an = samples[m][n];
                }
                else if (samples[m][n] >= bn) {
                    bn = samples[m][n];
                }
            }

            _intervals[n][0] = an - 0.5; //the 0.5 is just to add some tolerance
            _intervals[n][1] = bn + 0.5;

            std::cout << "a = " << _intervals[n][0] << " , " << "b = " << _intervals[n][1] << std::endl;

            for (unsigned l = 0; l < _L; l++) {
                _hs[n][l] = fabs (_intervals[n][1] - _intervals[n][0]) / pow (2, l + 1); // mesh size of the n-th interval
                std::cout << "h[" << n << "][" << l << "]= " << _hs[n][l] << std::endl;
                unsigned hierarchicalCounter = 0;
                for (unsigned i = 0; i < _nodes[n][l].size(); i++) {
                    _nodes[n][l][i] = _intervals[n][0] + (i + 1) * _hs[n][l]; //coordinates of the nodes of the n-th interval
                    std::cout << "node[" << n << "][" << l << "][" << i << "]= " << _nodes[n][l][i] << std::endl;
                    if (i % 2 == 0) {
                        _hierarchicalDofs[n][l][hierarchicalCounter] = i;
                        hierarchicalCounter++;
                    }
                }
            }

        }

        for (unsigned n = 0; n < _N; n++) {
            for (unsigned l = 0; l < _L; l++) {
                for (unsigned i = 0; i < _hierarchicalDofs[n][l].size(); i++) {
                    std::cout << "_hierarchicalDofs[" << n << "][" << l << "][" << i << "] = " << _hierarchicalDofs[n][l][i] << " " ;
                }
                std::cout << std::endl;
            }
        }

        //END


        //BEGIN computation of the sparse grid index set

        //this is done by first computing the tensor product space of the complete space and then eliminating those indices that don't satisfy the requirement: | |_1 <= _L + N - 1 (this condition holds assuming all indices start from 1)


        //here we compute the tensor produc index set for the full grid
        std::vector < std::vector < unsigned > > Tp;

        ComputeTensorProductSet (Tp, _L, _N);

        unsigned tensorProductDim = Tp.size();

        //here we only consider the indices of Tp that satisfy the requirement
        unsigned indexCounter = 0;

        for (unsigned i = 0; i < tensorProductDim; i++) {

            unsigned sum = 0;

            for (unsigned j = 0; j < _N; j++) {
                sum += (Tp[i][j] + 1); //this is because our indices start from 0 but the requirement assumes they start from 1
            }

            if (sum <= _L + _N - 1) {

                _indexSet.resize (indexCounter + 1);
                _indexSet[indexCounter].resize (_N);

                for (unsigned j = 0; j < _N; j++) {
                    _indexSet[indexCounter][j] = Tp[i][j];
                }
                indexCounter++;
            }
        }

        for (unsigned i = 0; i < _indexSet.size(); i++) {
            for (unsigned j = 0; j < _N; j++) {
                std::cout << "_indexSet[" << i << "][" << j << "]= " << _indexSet[i][j];
            }
            std::cout << std::endl;
        }

        _numberOfWs = _indexSet.size();


        //END


        //BEGIN construction of dofIdentifier

        _dofIdentifier.resize (_numberOfWs);

        std::vector< unsigned > maxDofs (_numberOfWs);

        for (unsigned w = 0; w < _numberOfWs; w++) {
            unsigned identifiersOfW = 1;
            for (unsigned i = 0; i < _N; i++) {

                identifiersOfW *= _hierarchicalDofs[i][_indexSet[w][i]].size();

                if (_nodes[i][_indexSet[w][i]].size() >= maxDofs[w]) maxDofs[w] = _nodes[i][_indexSet[w][i]].size();

            }

            _dofIdentifier[w].resize (identifiersOfW);
        }

        for (unsigned w = 0; w < _numberOfWs; w++) {
            if (maxDofs[w] < 2) maxDofs[w] = 2; //otherwise ComputeTensorProductSet doesn't work
            std::cout << "maxDofs[" << w << "] = " << maxDofs[w] << std::endl;
        }


        //TODO continue from here...
        for (unsigned w = 0; w < _numberOfWs; w++) {

            std::vector <std::vector <unsigned> > Vdofs;
            ComputeTensorProductSet (Vdofs, maxDofs[w], _N);

            std::cout << " ---------------------------------------------------- " << std::endl;

            std::vector <std::vector <unsigned> >  nodesToKeep (_N);

            unsigned counterKeep = 0;
            for (unsigned i = 0; i < Vdofs.size(); i++) {
                bool skip = false;
                for (unsigned n = 0; n < _N; n++) {
                    for (unsigned j = 0; j < _hierarchicalDofs[n][_indexSet[w][n]].size(); j++) {

                        if (Vdofs[i][n] != _hierarchicalDofs[n][_indexSet[w][n]][j]) {
                            skip = true;
                        }

                        if (skip) break;

                    }

                    if (skip) break;
                }
                
                if (!skip){
                    
                    nodesToKeep.resize(counterKeep + 1);
                    for (unsigned n = 0; n < _N; n++) {
                        nodesToKeep[counterKeep][n] = Vdofs[i][n];
                    }
                    counterKeep++;
                    
                }
                
            }
            
            //to erase
            for (unsigned ii = 0; ii < nodesToKeep.size(); ii++) {
                for (unsigned n = 0; n < _N; n++) {
                    std::cout << "nodesToKeep[" << ii << "][" << n << "] = " << nodesToKeep[ii][n];
                }
                std::cout << std::endl;
            }
            

        }



        //END

    }

    void sparseGrid::EvaluateOneDimensionalPhi (double &phi, const double &x, const unsigned &n, const unsigned &l, const unsigned &i, const bool &scale)
    {

        //i tells you that you are going to compute the phi associated with the i-th node

        double leftBoundOfSupport = _nodes[n][l][i] - _hs[n][l]; // this is xi - h
        double rightBoundOfSupport = _nodes[n][l][i] + _hs[n][l]; // this is xi + h

        if (x >= leftBoundOfSupport && x <= rightBoundOfSupport) {

            phi = 1 - fabs ( (x - _nodes[n][l][i]) / _hs[n][l]);

            if (scale == true) {  //this is to be used when building the nodal values of the PDF approximation
                phi /= _hs[n][l];
            }

        }

        else phi = 0.;
    }

    void sparseGrid::EvaluatePhi (double &phi, const std::vector <double> &x, std::vector < std::vector < unsigned > > identifier)
    {

        //identifier tells us what one dimensional phis to multiply in order to get the desired phi
        //identifier[n][0] = n (dim), identifier[n][1] = l (level), identifier[n][2] = i (node)

        std::vector < double > oneDimPhi (_N, 0.);

        phi = 1.;

        for (unsigned n = 0; n < _N; n++) {
            EvaluateOneDimensionalPhi (oneDimPhi[n], x[n], identifier[n][0], identifier[n][1], identifier[n][2], false);
            phi *= oneDimPhi[n];
        }

    }

    void sparseGrid::ComputeTensorProductSet (std::vector< std::vector <unsigned>> &Tp, const unsigned &T1, const unsigned &T2)
    {

        unsigned tensorProductDim = pow (T1, T2);

        Tp.resize (tensorProductDim);
        for (unsigned i = 0; i < tensorProductDim; i++) {
            Tp[i].resize (T2);
        }

        unsigned index = 0;
        unsigned counters[T2 + 1];
        memset (counters, 0, sizeof (counters));

        while (!counters[T2]) {

            for (unsigned j = 0; j < T2; j++) {
                Tp[index][j] = counters[T2 - 1 - j];
                std::cout << " Tp[" << index << "][" << j << "]= " << Tp[index][j] ;
            }
            std::cout << std::endl;
            index++;

            unsigned i;
            for (i = 0; counters[i] == T1 - 1; i++) {   // inner loops that are at maxval restart at zero
                counters[i] = 0;
            }
            ++counters[i];  // the innermost loop that isn't yet at maxval, advances by 1
        }

    }


}











