#ifndef __femus_enums_MgSmootherEnum_hpp__
#define __femus_enums_MgSmootherEnum_hpp__

enum MgSmoother {
    FEMuS_DEFAULT = 0,
    FEMuS_ASM,
    FEMuS_FIELDSPLIT,
};

enum CoarseLevelInclude{
  INCLUDE_COARSE_LEVEL_FALSE = 0,
  INCLUDE_COARSE_LEVEL_TRUE 
};

#endif
