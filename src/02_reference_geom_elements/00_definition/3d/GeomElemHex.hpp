#ifndef __femus_mesh_GeomElemHex_hpp__
#define __femus_mesh_GeomElemHex_hpp__


#include "GeomElemBase.hpp"


namespace femus {



class GeomElemHex : public GeomElemBase  {

public:
    
    GeomElemHex() : GeomElemBase() { };

    
    unsigned int  get_dimension() const { return _dim; };
    unsigned int n_nodes_linear() const { return _n_vertices; };


   static constexpr unsigned int _dim        =  3;
   static constexpr unsigned int _n_vertices =  8;
    
   virtual const unsigned int num_non_triangular_faces()  const {  return _number_of_non_triangular_faces; }
   virtual const unsigned int num_triangular_faces()  const {  return _number_of_triangular_faces; }
   
private:
   static constexpr unsigned int _number_of_non_triangular_faces = 6;
   
   static constexpr unsigned int _number_of_triangular_faces = 0;

};


} //end namespace femus



#endif
