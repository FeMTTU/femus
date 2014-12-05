#include "FEElemBase.hpp"

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "Files.hpp"

#include "GeomElTypeEnum.hpp"
#include "GeomEl.hpp"

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



FEElemBase::FEElemBase(const GeomEl & geomel_in )   : _geomel(geomel_in) {

  _myelems.resize(VB);  //cannot go in build function, because it is static

}

FEElemBase::~FEElemBase() {

  //TODO done in init(), do appropriately!!!!!!!!!!!!!!!
//     for (uint vb=0; vb < VB; vb++)  {
//
//     delete []       _phi_mapVB[vb];
//     delete []  _dphidxez_mapVB[vb];
//
//   }

}


//static function: it cannot act on the data of each instantiation ...

//this build class allows me to return a pointer to a child of this class
//even if i am a father
//inside this class the children must be INSTANTIATED, and then they are returned.
//These instantiations are never destroyed until you explicitly delete them
//the build() function returns a POINTER

FEElemBase* FEElemBase::build(const GeomEl & geomel_in, const uint fe_family_in) {

  if ( fe_family_in != QQ && fe_family_in != LL && fe_family_in != KK ) {
    std::cout << "FE::FE: FE family " << fe_family_in << " not supported" << std::endl;
    abort();
  }
  
  
     if (!strcmp(geomel_in._geomel_id.c_str(),"hex")) {  
       
      switch(fe_family_in) {
      case(QQ):
        return new  FEHex27(geomel_in)  ;
      case(LL):
        return new  FEHex8(geomel_in)  ;
      case(KK):
        return new  FEHex1(geomel_in)  ;
      }
      
      }
      else if (!strcmp(geomel_in._geomel_id.c_str(),"wedge")) {
           std::cout << "Not supported yet" << std::endl; abort();
      }
      else if (!strcmp(geomel_in._geomel_id.c_str(),"tet")) {
	
      switch(fe_family_in) {
      case(QQ):
        return new  FETet10(geomel_in)  ;
      case(LL):
        return new  FETet4(geomel_in)  ;
      case(KK):
        return new  FETet1(geomel_in)  ;
      }

      }
      else if (!strcmp(geomel_in._geomel_id.c_str(),"quad")) {
	
      switch(fe_family_in) {
      case(QQ):
        return new  FEQuad9(geomel_in)  ;  //FELagrange2D order2 on quadr
      case(LL):
        return new  FEQuad4(geomel_in)  ;  //FELagrange2D order1 on quadr
      case(KK):
        return new  FEQuad1(geomel_in)  ;  //FELagrange2D order0 on quadr
      }
      
      }
      else if (!strcmp(geomel_in._geomel_id.c_str(),"tri")) {
	
      switch(fe_family_in) {
      case(QQ):
        return new  FETri6(geomel_in)  ;
      case(LL):
        return new  FETri3(geomel_in)  ;
      case(KK):
        return new  FETri1(geomel_in)  ;
      }
	
      }
      
      else if (!strcmp(geomel_in._geomel_id.c_str(),"line")) { 
	
      switch(fe_family_in) {
      case(QQ):
        return new  FEEdge3(geomel_in)  ;
      case(LL):
        return new  FEEdge2(geomel_in)  ;
      case(KK):
        return new  FEEdge1(geomel_in)  ; 
      }
      
      }
      else {std::cout << "Geometric element not recognized" << std::endl; abort();}
  
  
}



//======================
// =======================
// on VB it is a loop
// QUADR o TRIANG is not a loop
// DIM is not a loop
// VB is a loop
// the order is not a loop

// The QRULE has been assigned already by the user based on the  MAXIMUM FE Order
//but of course you dont want to template on the NUMBER OF GAUSS POINTS...
// TODO here the thing is that HERE THERE is a MIXING of GEOMETRY and FE which should be removed!
// The FE should be geometry-agnostic!!!

//Il problema e' che ci sono molte incongruenze:
//GeomEl e' solo in 2D e 3D pero' ha delle cose basate su VB
//FEElemBase e' in 1-2-3D
// la embedding_matrix dovrebbe essere una cosa di GeomEl anziche' di FE

//Inoltre, FEelem = FEElem(GeomEl )
//          Qrule = Qrule(GeomEl )

// In realta', siccome la Qrule e' asservita a FEElem, quando si associa la Qrule a FEElem
// si puo' COSTRUIRE la Qrule che altrimenti rimane "nel limbo.."

// Allora noi nell'assemblaggio avremo il CurrentGaussPoint che prendera' dal FE le i valori delle phi e dphi nei nodi.

// The thing is that first you set the QRule and then you assign it to the FE class
// so in this init i have to say "if the quadrature is the Fifth, and the geomel is this (Quad9, Hex27, Tet10),
// then these are the evaluations of the quadrature points

// See, here is the thing: we are doing cases based on the ELEMENT  TYPE. There is no point in doing the cases for the element type TWICE,
// first in the quadrature rule and then in the FEElem...
// well, the determination of the QUADRATURE POINTS and the respective WEIGHT only depends on the GEOM ELEMENT and the given QUADRATURE ORDER.

//In Libmesh the FE gives EVERYTHING: shape at points, weights of points, Jacobian of the transformation, and so on (of course it did that through the Qrule Previously)
// Difference with libmesh: we are NOT ADAPTIVE.
// Also, libmesh handles DIFFERENT GEOMETRIC ELEMENTS, so every time he has to RESIZE the vector of LOCAL SHAPE FUNCTION and LOCAL GAUSS POINTS
// LIBMESH: all geom elements
// DEAL II: only quadrilateral and hexahedrals
// NETGEN: only triangular
// LIFE V:
// SUNDANCE
// FENICS

// QUI c'e' un loop NGAUSS[vb,dim[vb],quadtri] x NDOFS[vb,dim[vb], feorder,quadtri]


void FEElemBase::evaluate_shape_at_qp(const uint fe_family_in, const char* gauss_order) {

  uint space_dim = _geomel._dim;

// ================================================================================
// ============================ begin switch fe order ===============================
// ================================================================================
  switch(fe_family_in) {  //what leads is the FE ORDER, the quadrature rule is usually mathematically chosen with respect to that

// ================================================================================
  case(QQ):  {  /*"biquadratic"*/
    
     if (!strcmp(_geomel._geomel_id.c_str(),"hex")) {  
        _myelems[VV] = new elem_type_3D("hex","biquadratic", gauss_order);
        _myelems[BB] = new elem_type_2D("quad","biquadratic", gauss_order);
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"wedge")) {
           std::cout << "Not supported yet" << std::endl; abort();
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"tet")) {
        _myelems[VV] = new elem_type_3D("tet","biquadratic", gauss_order);
        _myelems[BB] = new elem_type_2D("tri","biquadratic", gauss_order);
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"quad")) {
        _myelems[VV] = new elem_type_2D("quad","biquadratic", gauss_order); //TODO valgrind
        _myelems[BB] = new elem_type_1D("line","biquadratic", gauss_order); //TODO valgrind
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"tri")) {
        _myelems[VV] = new elem_type_2D("tri","biquadratic", gauss_order); //TODO valgrind
        _myelems[BB] = new elem_type_1D("line","biquadratic", gauss_order); //TODO valgrind
      }
      
      else if (!strcmp(_geomel._geomel_id.c_str(),"line")) { 
        _myelems[VV] = new elem_type_1D("line","biquadratic", gauss_order); //TODO valgrind
      }
      else {std::cout << "Geometric element not recognized" << std::endl; abort();}
    
    
    
    break;
  } //end feorder QQ

// ================================================================================

// ================================================================================
  case(LL): {  /*"linear"*/
    
    
     if (!strcmp(_geomel._geomel_id.c_str(),"hex")) {  
        _myelems[VV] = new elem_type_3D("hex","linear", gauss_order);
	_myelems[BB] = new elem_type_2D("quad","linear", gauss_order);
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"wedge")) {
           std::cout << "Not supported yet" << std::endl; abort();
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"tet")) {
        _myelems[VV] = new elem_type_3D("tet","linear", gauss_order);
	_myelems[BB] = new elem_type_2D("tri","linear", gauss_order);
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"quad")) {
        _myelems[VV] = new elem_type_2D("quad","linear", gauss_order);
	_myelems[BB] = new elem_type_1D("line","linear", gauss_order);
	
      
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"tri")) {
        _myelems[VV] = new elem_type_2D("tri","linear",  gauss_order);
	_myelems[BB] = new elem_type_1D("line","linear", gauss_order);
      }
      
      else if (!strcmp(_geomel._geomel_id.c_str(),"line")) { 
	_myelems[VV] = new elem_type_1D("line","linear", gauss_order);
      }
      else {std::cout << "Geometric element not recognized" << std::endl; abort();}
    

    
    
    
    
    break;
  } //end feorder LL
// ================================================================================


// ================================================================================
  case(KK): {  /*"constant"*/
    

     if (!strcmp(_geomel._geomel_id.c_str(),"hex")) {  
        _myelems[VV] = new elem_type_3D("hex","constant", gauss_order);
	_myelems[BB] = new elem_type_2D("quad","constant", gauss_order);
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"wedge")) {
           std::cout << "Not supported yet" << std::endl; abort();
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"tet")) {
        _myelems[VV] = new elem_type_3D("tet","constant", gauss_order);
	_myelems[BB] = new elem_type_2D("tri","constant", gauss_order);
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"quad")) {
        _myelems[VV] = new elem_type_2D("quad","constant", gauss_order);
	_myelems[BB] = new elem_type_1D("line","constant", gauss_order);
	
      
      }
      else if (!strcmp(_geomel._geomel_id.c_str(),"tri")) {
        _myelems[VV] = new elem_type_2D("tri","constant",  gauss_order);
	_myelems[BB] = new elem_type_1D("line","constant", gauss_order);
      }
      
      else if (!strcmp(_geomel._geomel_id.c_str(),"line")) { 
	_myelems[VV] = new elem_type_1D("line","constant", gauss_order);
      }
      else {std::cout << "Geometric element not recognized" << std::endl; abort();}
    
    
    
    break;
  }  //end feorder KK

// ================================================================================
  default: {
    std::cout << "FE family not implemented" << std::endl;
    abort();
    break;
  }

  }  //end switch fe order
// ================================================================================
// ============================ end switch fe order ===============================
// ================================================================================


// ============== allocate canonical shape ==================================================================
for (int vb=0; vb<VB; vb++) {

    uint dim = space_dim - vb;

    _phi_mapVBGD[vb] = new double*[_myelems[vb]->GetGaussRule().GetGaussPointsNumber()];// TODO valgrind, remember to DEALLOCATE THESE
    _dphidxez_mapVBGD[vb] = new double*[_myelems[vb]->GetGaussRule().GetGaussPointsNumber()];

    for (int g = 0; g < _myelems[vb]->GetGaussRule().GetGaussPointsNumber(); g++) {
      _phi_mapVBGD[vb][g] = new double[_myelems[vb]->GetNDofs()];
      _dphidxez_mapVBGD[vb][g] = new double[_myelems[vb]->GetNDofs()*dim];
    }

  } //end VB
// ============== allocate canonical shape ==================================================================


  // loop ===========================
  // loop ===========================
  for (int vb=0; vb<VB; vb++) {

 
// HEX 27 CASE ========================================== 
// HEX 27 CASE ========================================== 
// HEX 27 CASE ========================================== 
// TRICK: for the HEX27 let me put the MAP for converting to my connectivity    
// from eu connectivity to my (=libmesh) connectivity
const unsigned map_hex27[27] = {0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15,24,20,21,22,23,25,26};

if ( vb == VV && fe_family_in == QQ && space_dim == 3  && (!strcmp(_geomel._geomel_id.c_str(),"hex")) ) {
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << "REMEMBER THAT ONLY HEX27 HAS A DIFFERENT CONNECTIVITY MAP"  << std::endl;

  
      for (int ig = 0; ig < _myelems[vb]->GetGaussRule().GetGaussPointsNumber(); ig++) {

      for (int idof=0; idof < _myelems[vb]->GetNDofs(); idof++) {
//                 std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << vb << " " << ig << " " << idof << std::endl;
        _phi_mapVBGD[vb][ig][idof] = _myelems[vb]->GetPhi(ig)[ map_hex27[idof] ];
// 	std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << vb << " " << ig << " " << dof << " phi " << _phi_mapVBGD[vb][ig][dof] << std::endl;

// derivatives in canonical element
        uint dim = space_dim - vb;
        for (uint idim = 0; idim < dim; idim++) {
// 		 double* temp =  ( _myelems[vb]->*(_myelems[vb]->Dphiptr[vb][idim]) )(ig);  //how to access a pointer to member function
          double* tempTwo =  ( _myelems[vb]->*(_myelems[vb]->_DPhiXiEtaZetaPtr[idim]) )(ig);  //how to access a pointer to member function
          _dphidxez_mapVBGD[vb][ig][ idof + idim*_myelems[vb]->GetNDofs()] =  tempTwo[ map_hex27[idof] ];
          std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << vb << " " << ig << " " << idof << " " << idim << " dphi         " << _dphidxez_mapVBGD[vb][ig][ idof + idim*_myelems[vb]->GetNDofs()]  << "                                      "  << std::endl;

        }

      }
      
    }  // end gauss
  
  
  
}
// HEX 27 CASE ========================================== 
// HEX 27 CASE ========================================== 
// HEX 27 CASE ========================================== 


   else { 

    for (int ig = 0; ig < _myelems[vb]->GetGaussRule().GetGaussPointsNumber(); ig++) {

      for (int idof=0; idof < _myelems[vb]->GetNDofs(); idof++) {
//                 std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << vb << " " << ig << " " << idof << std::endl;
        _phi_mapVBGD[vb][ig][idof] = _myelems[vb]->GetPhi(ig)[idof];
// 	std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << vb << " " << ig << " " << dof << " phi " << _phi_mapVBGD[vb][ig][dof] << std::endl;

// derivatives in canonical element
        uint dim = space_dim - vb;
        for (uint idim = 0; idim < dim; idim++) {
// 		 double* temp =  ( _myelems[vb]->*(_myelems[vb]->Dphiptr[vb][idim]) )(ig);  //how to access a pointer to member function
          double* tempTwo =  ( _myelems[vb]->*(_myelems[vb]->_DPhiXiEtaZetaPtr[idim]) )(ig);  //how to access a pointer to member function
          _dphidxez_mapVBGD[vb][ig][ idof + idim*_myelems[vb]->GetNDofs()] =  tempTwo[idof];
          std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << vb << " " << ig << " " << idof << " " << idim << " dphi         " << _dphidxez_mapVBGD[vb][ig][ idof + idim*_myelems[vb]->GetNDofs()]  << "                                      "  << std::endl;

        }

      }
      
    }  // end gauss
    
    
  } //else trick 
    
  } //end VB for
// loop ===========================
// loop  ===========================

  return;
}



} //end namespace femus


