#ifndef __femus_app_specifics_hpp__
#define __femus_app_specifics_hpp__


#include <vector>
#include <string>


class app_specifics {
 
  public:
      
      app_specifics() {  _mesh_files.resize(2); }
  
   std::vector< std::string >   _mesh_files;   
      
    
};



#endif
