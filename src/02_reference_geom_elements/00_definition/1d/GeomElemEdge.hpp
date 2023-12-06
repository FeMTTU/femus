#ifndef __femus_mesh_GeomElemEdge_hpp__
#define __femus_mesh_GeomElemEdge_hpp__


#include "GeomElemBase.hpp"


namespace femus {



class GeomElemEdge : public GeomElemBase  {

public:
    
    GeomElemEdge() : GeomElemBase() {
       set_faceNumber_offsets();
   }
         
    unsigned int  get_dimension() const { return _dim; };
    unsigned int n_nodes_linear() const { return _n_vertices; };

   static constexpr unsigned int _dim        =  1;
   static constexpr unsigned int _n_vertices =  2;
   
   virtual const unsigned int num_quadrilateral_faces()  const {  return _number_of_quadrilateral_faces; }
   virtual const unsigned int num_non_quadrilateral_faces()  const {  return _number_of_non_quadrilateral_faces; }
   
private:
   
   static constexpr unsigned int _number_of_quadrilateral_faces = 0;
   
   static constexpr unsigned int _number_of_non_quadrilateral_faces = 2;

};


} //end namespace femus



#endif
