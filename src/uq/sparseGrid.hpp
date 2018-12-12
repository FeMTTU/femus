



#ifndef __sparseGrid_hpp__
#define __sparseGrid_hpp__

#include <vector>
#include <map>


namespace femus
{

class sparseGrid
{

public:

    sparseGrid(const unsigned &_N, const unsigned &_M, std::vector < std::vector < double > >  &samples);
    
    void EvaluateOneDimensionalPhi (double &phi, const double &x, const unsigned &n, const unsigned &l, const unsigned &i, const bool &scale);
    
private:
    //defining parameters
    unsigned _N; //number of dimensions of the parameter space
    unsigned _M; //number of samples for each dimension, has to be expressed as M=pow(10,alpha), with unsigned alpha
    unsigned _L; //number of levels (L=log10(M)+1) I have decided to relate M and L in this way
    
    //monodimensional quantities
    std::vector < std::vector < double > > _intervals; 
    std::vector < std::vector < double > > _hs;
    std::vector < std::vector < std::vector < double > > > _nodes;
    
};

}

#endif
