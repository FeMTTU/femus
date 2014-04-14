
#include <iostream>
#include "Graph.hpp"


namespace femus {






void Graph::print() {

  for (uint i=0; i< (*this).size(); i++) {
    std::cout << "  "<< i << "  " << (*this)[i].size() << ": ";
    for (uint j=0; j< (*this)[i].size(); j++) {
      std::cout <<  (*this)[i][j] << "  ";
    }
    std::cout << " \n";
  }



}


} //end namespace femus


