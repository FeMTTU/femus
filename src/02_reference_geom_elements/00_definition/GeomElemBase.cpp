
#include "GeomElemBase.hpp"



#include "GeomElTypeEnum.hpp"
#include "GeomElemEdge2.hpp"
#include "GeomElemEdge3.hpp"
#include "GeomElemQuad4.hpp"
#include "GeomElemQuad9.hpp"
#include "GeomElemTri3.hpp"
#include "GeomElemTri6.hpp"
#include "GeomElemHex8.hpp"
#include "GeomElemHex27.hpp"
#include "GeomElemTet4.hpp"
#include "GeomElemTet10.hpp"


#include "FElemTypeEnum_list.hpp"

#include <iostream>
#include <cstdlib>
#include <cstring>


namespace femus {



  GeomElemBase::GeomElemBase() { }

  GeomElemBase::~GeomElemBase() { }



//this build class allows me to return a pointer to a child of this class
//even if i am a father
//inside this class the children must be INSTANTIATED, and then they are returned.
//These instantiations are never destroyed until you explicitly delete them
//the build() function returns a POINTER


  GeomElemBase* GeomElemBase::build (const std::string geomel_id_in, const uint fe_family_in) {

    if (fe_family_in >= NFE_FAMS_C_ZERO_LAGRANGE) {
      std::cout << "FE::FE: FE family " << fe_family_in << " not supported" << std::endl;
      abort();
    }


    if (!strcmp (geomel_id_in.c_str(), "hex")) {

      switch (fe_family_in) {
        case (0) :
          return new  GeomElemHex8()  ;
        case (1) :
          abort();
        case (2) :
          return new  GeomElemHex27()  ;
      }

    }

    else if (!strcmp (geomel_id_in.c_str(), "tet")) {

      switch (fe_family_in) {
        case (0) :
          return new  GeomElemTet4()  ;
        case (1) :
          abort();
        case (2) :
          return new  GeomElemTet10() ;
      }


    }

    else if (!strcmp (geomel_id_in.c_str(), "wedge")) {
      std::cout << "Not supported yet" << std::endl;
      abort();
    }

    else if (!strcmp (geomel_id_in.c_str(), "quad")) {

      switch (fe_family_in) {
        case (0) :
          return new  GeomElemQuad4()  ;
        case (1) :
          abort();
        case (2) :
          return new  GeomElemQuad9() ;
      }
    }

    else if (!strcmp (geomel_id_in.c_str(), "tri")) {

      switch (fe_family_in) {
        case (0) :
          return new GeomElemTri3();
        case (1) :
          abort();
        case (2) :
          return new  GeomElemTri6();
      }

    }

    else if (!strcmp (geomel_id_in.c_str(), "line")) {

      switch (fe_family_in) {
        case (0) :
          return new GeomElemEdge2();
        case (1) :
          abort();
        case (2) :
          return new GeomElemEdge3();
      }

    }

    else {
      std::cout << "Geometric element not recognized" << std::endl;
      abort();

    }

  }


} //end namespace femus


