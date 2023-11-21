#ifndef __femus_enums_FSIEnum_hpp__
#define __femus_enums_FSIEnum_hpp__


enum  FSIMaterialIndex_SolidPorousFluid { 
  SOLID_FLAG_INDEX = 0,
  POROUS_FLAG_INDEX,
  FLUID_FLAG_INDEX
};

enum  FSIMaterialIndex_SolidFluid { 
  SOLID_FLAG_INDEX_SolidFluid = 0,
  FLUID_FLAG_INDEX_SolidFluid = 1
};


///These values follow the Gambit material conventions
enum  FSIMaterialValue { 
  SOLID_FLAG_VALUE  = 4,
  POROUS_FLAG_VALUE = 3,
  FLUID_FLAG_VALUE  = 2
};

  constexpr unsigned int FSIMaterialIndex_SolidPorousFluid_count = FLUID_FLAG_INDEX - SOLID_FLAG_INDEX + 1;
     
  constexpr unsigned int FSIMaterialIndex_SolidFluid_count = FLUID_FLAG_INDEX_SolidFluid - SOLID_FLAG_INDEX_SolidFluid + 1;

#endif
