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

        for (unsigned n = 0; n < _N; n++) {
            _intervals[n].resize (2);
            _hs[n].resize (_L);
            _nodes[n].resize (_L);
            for (unsigned l = 0; l < _L; l++) {
                unsigned size = static_cast<unsigned> (pow (2, l + 1) - 1);
                _nodes[n][l].resize (size);
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
                for (unsigned i = 0; i < _nodes[n][l].size(); i++) {
                    _nodes[n][l][i] = _intervals[n][0] + (i + 1) * _hs[n][l]; //coordinates of the nodes of the n-th interval
                    std::cout << "node[" << n << "][" << l << "][" << i << "]= " << _nodes[n][l][i] << std::endl;
                }
            }

        }

        //END

    }
    
    void sparseGrid::EvaluateOneDimensionalPhi (double &phi, const double &x, const unsigned &n, const unsigned &l, const unsigned &i, const bool &scale){
        
        //i tells you that you are going to compute the phi associated with the i-th node
        
        double leftBoundOfSupport = _nodes[n][l][i] - _hs[n][l]; // this is xi - h
        double rightBoundOfSupport = _nodes[n][l][i] + _hs[n][l]; // this is xi + h
        
        if( x >= leftBoundOfSupport && x <= rightBoundOfSupport ) {
            
            phi = 1 - fabs( (x - _nodes[n][l][i]) / _hs[n][l] );
            
            if(scale == true) {   //this is to be used when building the nodal values of the PDF approximation
                phi /= _hs[n][l];
            }
            
        }
        
        else phi = 0.;
    }


}










