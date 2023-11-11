#ifndef __femus_meshGencase_GeomElemQuad4_hpp__
#define __femus_meshGencase_GeomElemQuad4_hpp__




#include "GeomElemBase.hpp"


namespace femus {



class GeomElemQuad4 : public GeomElemBase  {

public:
  
     GeomElemQuad4();
     
    ~GeomElemQuad4();
  
    
    unsigned int  get_dimension() const { return 2; };
    unsigned int n_nodes_linear() const { return 4; };
    
    unsigned int n_nodes()        const { return 4; };
    
    std::string   get_name_med()  const { return "QU4"; };
    std::string   get_name_xdmf() const { return "Quadrilateral"; };

// Refinement - BEGIN ===
public:
      
    float get_embedding_matrix(const uint,const uint,const uint);

               double get_prol(const uint j) {return _Prol[j];};
     
private:
      
      static const double _Prol[/*NNDS*/9*4/*NNDSL*/];
   
      static const float _embedding_matrix[4][4][4];   // (volume)
// Refinement - END ===
    
};


} //end namespace femus



#endif
