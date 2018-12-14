
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

    sparseGrid(const unsigned &_N, const unsigned &_M, std::vector < std::vector < double > >  &samples);
    
    //NOTE: this function evaluates any nodal phi, but for the Ws, we only need those associated with an even index i
    void EvaluateOneDimensionalPhi (double &phi, const double &x, const unsigned &n, const unsigned &l, const unsigned &i, const bool &scale);
    
    void EvaluatePhi(double &phi, const std::vector <double> &x, std::vector < std::vector < unsigned > > identifier);
    
    void ComputeTensorProductSet (std::vector< std::vector <unsigned>> &Tp, const unsigned &T1, const unsigned &T2);
    
private:
    //defining parameters
    unsigned _N; //number of dimensions of the parameter space
    unsigned _M; //number of samples for each dimension, has to be expressed as M=pow(10,alpha), with unsigned alpha
    unsigned _L; //number of levels (L=log10(M)+1) I have decided to relate M and L in this way
    
    //monodimensional quantities
    std::vector < std::vector < double > > _intervals; 
    std::vector < std::vector < double > > _hs;
    std::vector < std::vector < std::vector < double > > > _nodes;
    std::vector < std::vector < std::vector < double > > > _hierarchicalDofs;
    
    //sparse grid quantities
    unsigned _numberOfWs; // V = direct sum_{i=1,...,_numberOfWs} W_i, see (3.61) in Bungartz and Griebel, Acta Numerica 2004
    std::vector < std::vector < unsigned > > _indexSetW; //these are the indices of the W sets
    std::vector < std::vector < std::vector < std::vector < unsigned > > > > _dofIdentifier; //contains all identifiers 
    
};

}

#endif
