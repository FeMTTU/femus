
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

        sparseGrid ( std::vector < std::vector < double > >  &samples, const double &xmin, const double &xmax, const bool &output );

        void EvaluateOneDimensionalPhi ( double &phi, const double &x, const unsigned &n, const unsigned &l, const unsigned &i, const bool &scale );

        void EvaluatePhi ( double &phi, const std::vector <double> &x, std::vector < std::vector < unsigned > > identifier, const bool &scale );
        
        void InSupportOneDimensional ( unsigned &maybeThere, const double &x, const unsigned &n, const unsigned &l, const unsigned &i );
        
        void InSupport ( unsigned &isIt, const std::vector <double> &x, std::vector < std::vector < unsigned > > identifier );
        
        void PiecewiseConstPhi( double &phi, const std::vector <double> &x, std::vector < std::vector < unsigned > > identifier );

        void EvaluateNodalValuesPDF ( std::vector < std::vector < double > >  &samples );
        
        void EvaluatePDF (double &pdfValue, std::vector < double >  &x, const bool &print);
        
        void ComputeAvgL2Error( double &aL2E, std::vector < std::vector < double > >  &samples, const unsigned &analyticPdfType);
        
        void EvaluatePDFIntegral (double &integral);
        
        void SupportIsContained ( bool &check, const unsigned &w, const unsigned &i, const unsigned &wLower, const unsigned &iLower);

        void ComputeTensorProductSet ( std::vector< std::vector <unsigned>> &Tp, const unsigned &T1, const unsigned &T2 );
        
        void PrintNodalValuesPDF();
        
        void CartesianProduct ( std::vector<std::vector<int> >& inputCartesian, std::vector<std::vector<unsigned> >& dofsWi );

    private:
        //defining parameters
        unsigned _N; //number of dimensions of the parameter space
        unsigned _M; //number of samples for each dimension, has to be expressed as M=pow(10,alpha), with unsigned alpha
        unsigned _L; //number of levels 

        //monodimensional quantities
        std::vector < std::vector < double > > _intervals;
        std::vector < std::vector < double > > _hs;
        std::vector < std::vector < std::vector < double > > > _nodes; //number of internal nodes of the grid
        std::vector < std::vector < std::vector < double > > > _hierarchicalDofs;
        std::vector < std::vector < std::vector < double > > > _hierarchicalDofsCoordinates;

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
