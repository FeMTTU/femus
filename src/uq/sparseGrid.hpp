
#ifndef __sparseGrid_hpp__
#define __sparseGrid_hpp__

#include <vector>
#include <map>

#include <boost/random.hpp>

namespace femus
{

    class sparseGrid
    {

    public:

        sparseGrid ( std::vector < std::vector < double > >  &samples, const bool &output );

        void EvaluateOneDimensionalPhi ( double &phi, const double &x, const unsigned &n, const unsigned &l, const unsigned &i, const bool &scale );

        void EvaluatePhi ( double &phi, const std::vector <double> &x, std::vector < std::vector < unsigned > > identifier, const bool &scale );

        void EvaluateNodalValuesPDF ( std::vector < std::vector < double > >  &samples );
        
        void EvaluatePDF ( std::vector < double >  &x );

        void ComputeTensorProductSet ( std::vector< std::vector <unsigned>> &Tp, const unsigned &T1, const unsigned &T2 );
        
        void PrintNodalValuesPDF();

    private:
        //defining parameters
        unsigned _N; //number of dimensions of the parameter space
        unsigned _M; //number of samples for each dimension, has to be expressed as M=pow(10,alpha), with unsigned alpha
        unsigned _L; //number of levels 

        //monodimensional quantities
        std::vector < std::vector < double > > _intervals;
        std::vector < std::vector < double > > _hs;
        std::vector < std::vector < std::vector < double > > > _nodes;
        std::vector < std::vector < std::vector < double > > > _hierarchicalDofs;

        //sparse grid quantities
        unsigned _numberOfWs; // V = direct sum_{i=1,...,_numberOfWs} W_i, see (3.61) in Bungartz and Griebel, Acta Numerica 2004
        std::vector < std::vector < unsigned > > _indexSetW; //these are the indices of the W subspaces
        std::vector < std::vector < std::vector < std::vector < unsigned > > > > _dofIdentifier; //contains all identifiers
        //_dofIdentifier[*][][][]:ranges of W subspaces, _dofIdentifier[][*][][]: ranges over the dofs of the W subspace,
        //_dofIdentifier[][][*][]:ranges over the identifier of W, _dofIdentifier[][][][*]: ranges from 0 to 2, 0=dim, 1=level, 2=dof
        std::vector < std::vector < double > > _nodalValuesPDF; //nodal values of the approximate PDF

        //output parameter
        bool _output;

    };

}

#endif
