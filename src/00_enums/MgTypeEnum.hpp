#ifndef __femus_enums_MgTypeEnum_hpp__
#define __femus_enums_MgTypeEnum_hpp__

enum MgType {
    F_CYCLE=0,
    V_CYCLE,
    M_CYCLE
};

enum MgSmootherType {
    FULL=0,
    MULTIPLICATIVE,
    ADDITIVE,
    KASKADE
};

#endif
