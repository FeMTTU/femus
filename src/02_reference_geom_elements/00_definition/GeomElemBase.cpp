
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

// definition needed for static constexpr - BEGIN

      constexpr unsigned GeomElemBase::_n_face_types_max;
      constexpr unsigned GeomElemBase::_index_for_quadrilateral_faces;
      constexpr unsigned GeomElemBase::_index_for_all_faces          ;

      constexpr unsigned GeomElemBase::_max_space_dimension;

// definition needed for static constexpr - END


//this build class allows me to return a pointer to a child of this class
//even if i am a father
//inside this class the children must be INSTANTIATED, and then they are returned.
//These instantiations are never destroyed until you explicitly delete them
//the build() function returns a POINTER


  std::unique_ptr<GeomElemBase> GeomElemBase::build (const std::string geomel_id_in, const uint fe_family_in) {

    if (fe_family_in >= NFE_FAMS_C_ZERO_LAGRANGE) {
      std::cout << "FE::FE: FE family " << fe_family_in << " not supported" << std::endl;
      abort();
    }


    if (!strcmp (geomel_id_in.c_str(), "hex")) {

      switch (fe_family_in) {
        case (0) : 
        {
          std::unique_ptr<GeomElemBase>  el_ptr( new  GeomElemHex8() )  ;
          return el_ptr;
        }
        case (1) :
        {
          abort();
        }
        case (2) :
        {
          std::unique_ptr<GeomElemBase>  el_ptr( new  GeomElemHex27() ) ;
          return el_ptr;
        }
      }

    }

    else if (!strcmp (geomel_id_in.c_str(), "tet")) {

      switch (fe_family_in) {
        case (0) :
        {
          std::unique_ptr<GeomElemBase>  el_ptr( new  GeomElemTet4()  );
          return el_ptr;
        }
        case (1) :
        {
          abort();
        }
        case (2) :
        {
          std::unique_ptr<GeomElemBase>  el_ptr( new  GeomElemTet10() );
          return el_ptr;
        }
        
      }


    }

    else if (!strcmp (geomel_id_in.c_str(), "wedge")) {
      std::cout << "Not supported yet" << std::endl;
      abort();
    }

    else if (!strcmp (geomel_id_in.c_str(), "quad")) {

      switch (fe_family_in) {
        case (0) :
        {
          std::unique_ptr<GeomElemBase>  el_ptr( new  GeomElemQuad4()  );
          return el_ptr;
        }
        case (1) :
        {
          abort();
        }
        case (2) :
        {
          std::unique_ptr<GeomElemBase>  el_ptr( new  GeomElemQuad9() );
          return el_ptr;
        }
      }
    }

    else if (!strcmp (geomel_id_in.c_str(), "tri")) {

      switch (fe_family_in) {
        case (0) :
        {
          std::unique_ptr<GeomElemBase>  el_ptr( new GeomElemTri3() );
          return el_ptr;
        }
        case (1) :
        {
          abort();
        }
        case (2) :
        {
          std::unique_ptr<GeomElemBase>  el_ptr( new  GeomElemTri6() );
          return el_ptr;
        }
      }

    }

    else if (!strcmp (geomel_id_in.c_str(), "line")) {

      switch (fe_family_in) {
        case (0) :
        {
          std::unique_ptr<GeomElemBase>  el_ptr( new GeomElemEdge2() );
          return el_ptr;
        }
        case (1) :
        {
          abort();
        }
        case (2) :
        {
          std::unique_ptr<GeomElemBase>  el_ptr( new GeomElemEdge3() );
          return el_ptr;
        }
      }

    }

    else {
      std::cout << "Geometric element not recognized" << std::endl;
      abort();

    }

  }


} //end namespace femus


