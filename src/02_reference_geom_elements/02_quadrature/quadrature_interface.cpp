/*=========================================================================

  Program: FEMUS
  Module: Gauss
  Authors: Eugenio Aulisa, Giorgio Bornia

  Copyright (c) FEMTTU
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/

#include "quadrature_interface.hpp"

#include "quadrature_Point.hpp"
#include "quadrature_Line.hpp"
#include "quadrature_Triangle.hpp"
#include "quadrature_Quadrangle.hpp"
#include "quadrature_Hexahedron.hpp"
#include "quadrature_Tetrahedron.hpp"
#include "quadrature_Wedge.hpp"



#include <iostream>
#include <cstdlib>
#include <cstring>

namespace femus
{

  // ************** CONSTRUCTOR - BEGIN ***************
  Gauss::Gauss(const char *geom_elem, const char *order_gauss) : _order(order_gauss)
  {
    if (!strcmp(order_gauss, "zero")  || !strcmp(order_gauss, "first")) {
      gauss_order = 0;
    }
    else if (!strcmp(order_gauss, "second") || !strcmp(order_gauss, "third") ) {
      gauss_order = 1;
    }
    else if (!strcmp(order_gauss, "fourth") || !strcmp(order_gauss, "fifth") ) {
      gauss_order = 2;
    }
    else if (!strcmp(order_gauss, "sixth")  || !strcmp(order_gauss, "seventh") ) {
      gauss_order = 3;
    }
    else if (!strcmp(order_gauss, "eighth") || !strcmp(order_gauss, "ninth") ) {
      gauss_order = 4;
    }
    else if (!strcmp(order_gauss, "tenth") || !strcmp(order_gauss, "eleventh") ) {
      gauss_order = 5;
    }
    else {
      std::cout << order_gauss << "is not a valid option for the Gauss points of" << geom_elem << std::endl;
      abort();
    }

    if (!strcmp(geom_elem, "hex")) {
      GaussWeight = hex_gauss::Gauss[gauss_order];
      GaussPoints = hex_gauss::GaussPoints[gauss_order];
    }
    else if (!strcmp(geom_elem, "wedge")) {
      GaussWeight = wedge_gauss::Gauss[gauss_order];
      GaussPoints = wedge_gauss::GaussPoints[gauss_order];
    }
    else if (!strcmp(geom_elem, "tet")) {
      GaussWeight = tet_gauss::Gauss[gauss_order];
      GaussPoints = tet_gauss::GaussPoints[gauss_order];
    }
    else if (!strcmp(geom_elem, "quad")) {
      GaussWeight = quad_gauss::Gauss[gauss_order];
      GaussPoints = quad_gauss::GaussPoints[gauss_order];
    }
    else if (!strcmp(geom_elem, "tri")) {
      GaussWeight = tri_gauss::Gauss[gauss_order];
      GaussPoints = tri_gauss::GaussPoints[gauss_order];
    }
    else if (!strcmp(geom_elem, "line")) {
      GaussWeight = line_gauss::Gauss[gauss_order];
      GaussPoints = line_gauss::GaussPoints[gauss_order];
    }
    else if (!strcmp(geom_elem, "point")) {
      GaussWeight = point_gauss::Gauss[gauss_order];
      GaussPoints = point_gauss::GaussPoints[gauss_order];
    }
    else {
      std::cout << geom_elem << " is not a valid option" << std::endl;
      abort();
    }

  }
  // ************** CONSTRUCTOR - END ***************










} //end namespace femus
