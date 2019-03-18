#ifndef __femus_enums_MgSmootherEnum_hpp__
#define __femus_enums_MgSmootherEnum_hpp__

enum MgSmoother {
    FEMuS_DEFAULT_SMOOTHER = 0,
    ASM_SMOOTHER,
    FIELDSPLIT_SMOOTHER,
};

enum CoarseLevelInclude{
  INCLUDE_COARSE_LEVEL_FALSE = 0,
  INCLUDE_COARSE_LEVEL_TRUE 
};

#endif
