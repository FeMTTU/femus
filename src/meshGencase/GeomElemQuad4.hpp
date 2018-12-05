#ifndef __femus_meshGencase_GeomElemQuad4_hpp__
#define __femus_meshGencase_GeomElemQuad4_hpp__




#include "GeomElemBase.hpp"


namespace femus {



class GeomElemQuad4 : public GeomElemBase  {

public:
  
     GeomElemQuad4();
     
    ~GeomElemQuad4();
  
    
    unsigned int  get_dimension() const { return 2; };
    unsigned int n_nodes()        const { return 4; };
    std::string   get_name_med()  const { return "QU4"; };
    std::string   get_name_xdmf() const { return "Quadrilateral"; };

    float get_embedding_matrix(const uint,const uint,const uint);

               double get_prol(const uint j) {return _Prol[j];};
     static const double _Prol[/*NNDS*/9*4/*NNDSL*/];
     
     //Shapes at quadrature points
     //The quadrature is known at runtime... so, we should SWITCH in each element with all the possible
     //quadrature rules... beh, intanto possiamo fare la porcata e fare una routine di lettura con tutti gli switch nella classe padre

private:
    
      static const float _embedding_matrix[4][4][4];   // (volume)
    
};


} //end namespace femus



#endif
