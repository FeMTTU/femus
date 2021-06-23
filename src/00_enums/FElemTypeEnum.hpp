#ifndef __femus_enums_FElemTypeEnum_hpp__
#define __femus_enums_FElemTypeEnum_hpp__

enum  FEFamily {
    LAGRANGE = 0,
    DISCONTINUOUS_POLYNOMIAL
};

enum FEOrder {
    ZERO = 0,
    FIRST,
    SERENDIPITY,
    SECOND
};

#endif
