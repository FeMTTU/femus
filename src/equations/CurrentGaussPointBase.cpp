#include "CurrentGaussPointBase.hpp"

#include "MultiLevelProblemTwo.hpp"
#include "MeshTwo.hpp"
#include "GeomEl.hpp"

#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"
#include "ElemType.hpp"


namespace femus {





//I need to hold the equations map pointer, because i also need a mesh pointer
// for the geometric element
//maybe later on i'd just pass the GeomElement(GeomEl) and the MathElement(FE)
//by the way, with the MultiLevelProblemTwo I reach the Utils, the Mesh, and so the GeomEl, and so on...
CurrentGaussPointBase::CurrentGaussPointBase(const CurrentElem & curr_el_in, MultiLevelProblemTwo& e_map_in ):
 _current_elem(curr_el_in),
       _eqnmap(e_map_in),
    _elem_type(e_map_in._elem_type[curr_el_in.GetDim() - 1]),
        _qrule(e_map_in.    _qrule[curr_el_in.GetDim() - 1]) {
  
  _IntDim[VV] = _eqnmap._mesh.get_dim();
  _IntDim[BB] = _eqnmap._mesh.get_dim() - 1; 
  
  //TODO probabilmente anche qui si puo' fare del TEMPLATING!!!
  //BISOGNA STARE ATTENTI CHE SE FAI DEL TEMPLATING con le ALLOCAZIONI STATICHE allora ti diverti poco con i DOPPI o TRIPLI ARRAY

     for (int fe = 0; fe < QL; fe++) {
          _phi_ndsQLVB_g[fe] =  new double[                       _elem_type[fe]->GetNDofs() ];
   _dphidxyz_ndsQLVB_g3D[fe] =  new double[ 3                   * _elem_type[fe]->GetNDofs() ];
     _dphidxyz_ndsQLVB_g[fe] =  new double[ curr_el_in.GetDim() * _elem_type[fe]->GetNDofs() ];   
  _dphidxezeta_ndsQLVB_g[fe] =  new double[ curr_el_in.GetDim() * _elem_type[fe]->GetNDofs() ];     
   }
  
  //Jacobian matrices, normals, tangents
  
  _normal_g  = new double[_IntDim[VV]];
  _tangent_g = new double*[_IntDim[BB]];
  _InvJac_g  = new double*[_IntDim[VV]];
    for (int i = 0; i < _IntDim[VV]; i++) {  _InvJac_g[i] = new double[_IntDim[VV]]; }
    for (int i = 0; i < _IntDim[BB]; i++) { _tangent_g[i] = new double[_IntDim[VV]]; }
  
  
  
}


CurrentGaussPointBase::~CurrentGaussPointBase() {
  

      for (int j=0;j< QL;j++) {
	delete []  _dphidxyz_ndsQLVB_g[j];
	delete []  _dphidxyz_ndsQLVB_g3D[j];
	delete [] _dphidxezeta_ndsQLVB_g[j];
      }

          for (int j=0;j< QL;j++) delete [] _phi_ndsQLVB_g[j];
    
       for (int i = 0; i < _IntDim[BB]; i++) { delete [] _tangent_g[i];}
       for (int i = 0; i < _IntDim[VV]; i++) { delete [] _InvJac_g[i];}
   delete [] _tangent_g;
   delete [] _InvJac_g;
   delete [] _normal_g;
  
}



//this is what allows RUNTIME selection of the templates!!!
   CurrentGaussPointBase& CurrentGaussPointBase::build(const CurrentElem & elem_in, MultiLevelProblemTwo& eqmap_in, const uint dim_in) {
      
      
      switch(dim_in) {
	
	case(2):  return *(new  CurrentGaussPoint<2>(elem_in,eqmap_in));
	
	case(3):  return *(new  CurrentGaussPoint<3>(elem_in,eqmap_in)); 
	
	default: {std::cout << "CurrentGaussPointBase: Only 2D and 3D" << std::endl; abort();}
	  
      } //dim_in
      
      
    }


} //end namespace femus


