#ifndef __femus_mesh_MED_IO_group_info_hpp__
#define __femus_mesh_MED_IO_group_info_hpp__


#include "GeomElemBase.hpp"

#include <string>


namespace femus
{

 #define TYPE_FOR_INT_DATASET  int  //do not use "unsigned": in fact, in MED 
                                  // these numbers are NEGATIVE for elements,
                                  //    while they are POSITIVE for nodes!
 #define TYPE_FOR_REAL_DATASET  double

    
// Auxiliary struct to store group information
  struct GroupInfo {
      std::string        _user_defined_group_string;
      std::string        _node_or_cell_group;
      TYPE_FOR_INT_DATASET _med_flag;
      int                _user_defined_flag;
      unsigned int       _user_defined_property;
      GeomElemBase*      _geom_el;
      int _size;
  };
  
}




#endif
