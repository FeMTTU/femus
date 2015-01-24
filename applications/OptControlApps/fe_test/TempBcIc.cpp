//C++
#include <cmath>

//library headers
#include "Box.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "MultiLevelMeshTwo.hpp" 
#include "NormTangEnum.hpp"
#include "TimeLoop.hpp"

//application headers
#include "TempQuantities.hpp"
#include "EqnT.hpp"


// In these files, just like in the equation file, you need to have an EXPLICIT KNOWLEDGE 
// of the ORDER of the UNKNOWNS

//If for a given equation you decide to change the NUMBER of UNKNOWNS then it is a little bit of a problem
//in the sense that you cannot do it totally automatically, you need to revisit the EQUATION and the BC/IC.


void EqnT::ic_read(const double xp[], double u_value[], const double el_xm[]) const {

  Box* box = static_cast<Box*>(_mesh.GetDomain());
  
  double* x_rotshift = new double[_mesh.get_dim()];
  _mesh._domain->TransformPointToRef(xp,x_rotshift); 

    u_value[0] = 3.; 
    u_value[1] = 4.;  
    u_value[2] = 5.;

  delete[] x_rotshift;
  
  return;
  
}