#ifndef __mgqrule_h__
#define __mgqrule_h__


#include <string>

#include "Typedefs_conf.hpp"
#include "VBTypeEnum.hpp"




//the choice of the quadrature rule depends on the 
//chosen FE and also on external functions you use... 
//it depends on the OPERATORS appearing in MATRIX and RHS

// The choice of the quadrature rule depends on the ORDER of the OPERATORS in the equation.
// The user makes this choice.
// For simplicity we choose to have ONLY ONE QUADRATURE RULE for EACH EQUATION.
// Clearly, the QRule will depend on the GEOMETRIC ELEMENT
// Also, the FE element will depend on the GEOMETRIC ELEMENT
// So, when you attach a quadrature rule to a FE, you have to make sure that they are built on top of the SAME GeomElement
// Or, you can first build the FE with the GEomEl, and then you EXTRACT THAT GEOMEL and with that you build the Quadrature rule.
// Then with that new Qrule, you keep filling the phimap and dphidxezmap

// clearly, the quadrature rule does not know what FE will be there
// it actually doesn't necessarily need to know the GEOM ELEMENT. It can get it LATER with the ASSOCIATED FE
// So, the thing is that I should put the POINTS as well... but the points depend on the GEOMETRY...
// Since what we have are the EVALUATIONS at the points, we don't necessarily need info about the points themselves...
// but we need to know THE GEOMETRY of the element (reference cube, reference triangle, etc.)
// because that changes the number of points, point coordinates and weights as well!!!

// I guess all of this happens in the function ATTACH_QUADRATURE_RULE, which takes the GeomEl from FE

class GeomEl;


class QRule  {

  public:

//GeomEl =========
   GeomEl * _geomel;

//Quadrature =========
   std::string  _qrule_type = "Gauss5th";
    uint         _NoGaussVB[VB];
    double*       _weightVB[VB];
    
    
     QRule(GeomEl* geomel_in);
    
    ~QRule();
    
    
    
};

#endif