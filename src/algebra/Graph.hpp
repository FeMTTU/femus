#ifndef __femus_algebra_Graph_hpp__
#define __femus_algebra_Graph_hpp__

#include "Typedefs.hpp"

#include <vector>


namespace femus {



typedef std::vector<unsigned int> Row;

class Graph : public std::vector<Row> {

  public:

  uint _m;
  uint _n; 
  uint _ml;
  uint _nl; 
  uint _ndz;
  uint _noz; 
  uint _ml_start;
  
  Graph() : _m(0),_n(0),_ml(0),_nl(0),_ndz(0),_noz(0),_ml_start(0){}

  void print();
  
};


} //end namespace femus



#endif