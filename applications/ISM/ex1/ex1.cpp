
#include "FemusInit.hpp"
#include "Marker.hpp"

using namespace femus;

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  std::vector < double > x(3);
  x[0]=1;
  x[1]=3;
  x[2]=4;
  Marker a(x,INTERSECTION);
  
  
  std::vector < double > y = a.GetMarkerCoordinates();
  
  std::cout<<y[0] << " " << y[1] << " " <<y[2]<<std::endl;
  std::cout << a.GetMarkerType() <<std::endl;
  
  return 0;
}




