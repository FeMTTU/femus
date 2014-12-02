#include "CurrGaussPointBase.hpp"

#include "QRule.hpp"
#include "EquationsMap.hpp"
#include "FEElemBase.hpp"
#include "MeshTwo.hpp"
#include "GeomEl.hpp"

#include "CurrGaussPoint.hpp"


namespace femus {





//I need to hold the equations map pointer, because i also need a mesh pointer
// for the geometric element
//maybe later on i'd just pass the GeomElement(GeomEl) and the MathElement(FE)
//by the way, with the EquationsMap I reach the Utils, the Mesh, and so the GeomEl, and so on...
CurrGaussPointBase::CurrGaussPointBase( EquationsMap& e_map_in ):
    _eqnmap(e_map_in),_qrule(e_map_in._qrule) {
  
  _IntDim[VV] = _eqnmap._mesh.get_dim();
  _IntDim[BB] = _eqnmap._mesh.get_dim() - 1; 
  
  //TODO probabilmente anche qui si puo' fare del TEMPLATING!!!
  //BISOGNA STARE ATTENTI CHE SE FAI DEL TEMPLATING con le ALLOCAZIONI STATICHE allora ti diverti poco con i DOPPI o TRIPLI ARRAY

  for (int fe = 0; fe < QL; fe++)  _AbsFEVect[fe] = e_map_in._AbstractFE[VV][fe];   

  for (int vb = 0; vb < VB; vb++) {
     for (int fe = 0; fe < QL; fe++) {
   _dphidxyz_ndsQLVB_g3D[vb][fe] =  new double[ 3 * _AbsFEVect[fe]->_ndof[vb]]; //both VV and BB are 3 in general (vector product, or ONE?!?)
     _dphidxyz_ndsQLVB_g[vb][fe] =  new double[ _IntDim[vb] * _AbsFEVect[fe]->_ndof[vb]];   
  _dphidxezeta_ndsQLVB_g[vb][fe] =  new double[ _IntDim[vb] * _AbsFEVect[fe]->_ndof[vb]];     
          _phi_ndsQLVB_g[vb][fe] =  new double[ _AbsFEVect[fe]->_ndof[vb]];     
   }
 }  
  
  //Jacobian matrices, normals, tangents
  
  _normal_g  = new double[_IntDim[VV]];
  _tangent_g = new double*[_IntDim[BB]];
  _InvJac_g  = new double*[_IntDim[VV]];
    for (int i = 0; i < _IntDim[VV]; i++) {  _InvJac_g[i] = new double[_IntDim[VV]]; }
    for (int i = 0; i < _IntDim[BB]; i++) { _tangent_g[i] = new double[_IntDim[VV]];}
  
  
  
}


CurrGaussPointBase::~CurrGaussPointBase() {
  
    for (int i=0;i< VB;i++){
      for (int j=0;j< QL;j++) {
	delete [] _phi_ndsQLVB_g[i][j];
	delete []  _dphidxyz_ndsQLVB_g[i][j];
	delete []  _dphidxyz_ndsQLVB_g3D[i][j];
	delete [] _dphidxezeta_ndsQLVB_g[i][j];
      }
    }
   
       for (int i = 0; i < _IntDim[BB]; i++) { delete [] _tangent_g[i];}
       for (int i = 0; i < _IntDim[VV]; i++) { delete [] _InvJac_g[i];}
   delete [] _tangent_g;
   delete [] _InvJac_g;
   delete [] _normal_g;
  
}



//this is what allows RUNTIME selection of the templates!!!
   CurrGaussPointBase& CurrGaussPointBase::build(EquationsMap& eqmap_in, const uint dim_in) {
      
      
      switch(dim_in) {
	
	case(2):  return *(new  CurrGaussPoint<2>(eqmap_in));
	
	case(3):  return *(new CurrGaussPoint<3>(eqmap_in)); 
	
	default: {std::cout << "CurrGaussPointBase: Only 2D and 3D" << std::endl; abort();}
	  
      } //dim_in
      
      
    }


} //end namespace femus


