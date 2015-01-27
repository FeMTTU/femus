#include "FEElemBase.hpp"

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "Files.hpp"

#include "GeomElTypeEnum.hpp"
#include "FETypeEnum.hpp"
#include "FEQuad1.hpp"
#include "FEQuad4.hpp"
#include "FEQuad9.hpp"
#include "FEHex1.hpp"
#include "FEHex8.hpp"
#include "FEHex27.hpp"
#include "FETri1.hpp"
#include "FETri3.hpp"
#include "FETri6.hpp"
#include "FETet1.hpp"
#include "FETet4.hpp"
#include "FETet10.hpp"
#include "FEEdge1.hpp"
#include "FEEdge2.hpp"
#include "FEEdge3.hpp"



namespace femus {



FEElemBase::FEElemBase() { }

FEElemBase::~FEElemBase() { }


//this build class allows me to return a pointer to a child of this class
//even if i am a father
//inside this class the children must be INSTANTIATED, and then they are returned.
//These instantiations are never destroyed until you explicitly delete them
//the build() function returns a POINTER

FEElemBase* FEElemBase::build(const std::string geomel_id_in, const uint fe_family_in) {

  if ( fe_family_in != QQ && fe_family_in != LL && fe_family_in != KK ) {
    std::cout << "FE::FE: FE family " << fe_family_in << " not supported" << std::endl;
    abort();
  }
  
  
     if (!strcmp(geomel_id_in.c_str(),"hex")) {  
       
      switch(fe_family_in) {
      case(QQ):
        return new  FEHex27()  ;
      case(LL):
        return new  FEHex8()  ;
      case(KK):
        return new  FEHex1()  ;
      }
      
      }
      else if (!strcmp(geomel_id_in.c_str(),"wedge")) {
           std::cout << "Not supported yet" << std::endl; abort();
      }
      else if (!strcmp(geomel_id_in.c_str(),"tet")) {
	
      switch(fe_family_in) {
      case(QQ):
        return new  FETet10()  ;
      case(LL):
        return new  FETet4()  ;
      case(KK):
        return new  FETet1()  ;
      }

      }
      else if (!strcmp(geomel_id_in.c_str(),"quad")) {
	
      switch(fe_family_in) {
      case(QQ):
        return new  FEQuad9()  ;
      case(LL):
        return new  FEQuad4()  ;
      case(KK):
        return new  FEQuad1()  ;
      }
      
      }
      else if (!strcmp(geomel_id_in.c_str(),"tri")) {
	
      switch(fe_family_in) {
      case(QQ):
        return new  FETri6()  ;
      case(LL):
        return new  FETri3()  ;
      case(KK):
        return new  FETri1()  ;
      }
	
      }
      
      else if (!strcmp(geomel_id_in.c_str(),"line")) { 
	
      switch(fe_family_in) {
      case(QQ):
        return new  FEEdge3()  ;
      case(LL):
        return new  FEEdge2()  ;
      case(KK):
        return new  FEEdge1()  ; 
      }
      
      }
      else {std::cout << "Geometric element not recognized" << std::endl; abort();}
  
  
}


} //end namespace femus


