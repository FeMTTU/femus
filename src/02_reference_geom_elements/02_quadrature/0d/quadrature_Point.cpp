

#include "quadrature_Point.hpp"


namespace femus {


// ************** POINT - BEGIN ***************
  const unsigned point_gauss::GaussPoints[5] = {1, 1, 1, 1, 1};
  const double * point_gauss::Gauss[5] = { Gauss0[0], Gauss1[0], Gauss2[0], Gauss3[0], Gauss4[0]};


//first row-weights, second row: x-coordinates
  const double point_gauss::Gauss0[2][1] = {{1}, {0}};

  const double point_gauss::Gauss1[2][1] = {{1}, {0}};

  const double point_gauss::Gauss2[2][1] = {{1}, {0}};

  const double point_gauss::Gauss3[2][1] = {{1}, {0}};

  const double point_gauss::Gauss4[2][1] = {{1}, {0}};
// ************** POINT - END ***************

}
