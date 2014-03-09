#include "FEElemBase.hpp"

#include <iostream>
#include <cstdlib>
#include <cstring>

#include "Utils.hpp"
#include "Files.hpp"
#include "QRule.hpp"

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

#include "ElemType.hpp"

FEElemBase::FEElemBase(GeomEl* geomel_in ) {

  _geomel = geomel_in;
  _n_children = _geomel->n_se[VV];

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


FEElemBase* FEElemBase::build(GeomEl* geomel_in, const uint order) {


  switch(geomel_in->_dim) {

  case(2): {

    switch(geomel_in->_geomel_type) {
    case(QUADR): {
      switch(order) {
      case(QQ):
        return new  FEQuad9(geomel_in)  ;  //FELagrange2D order2 on quadr
      case(LL):
        return new  FEQuad4(geomel_in)  ;  //FELagrange2D order1 on quadr
      case(KK):
        return new  FEQuad1(geomel_in)  ;  //FELagrange2D order0 on quadr
      }
    }
    case(TRIANG): {
      switch(order) {
      case(QQ):
        return new  FETri6(geomel_in)  ;
      case(LL):
        return new  FETri3(geomel_in)  ;
      case(KK):
        return new  FETri1(geomel_in)  ;
      }
    }
    }

  } //dim 2

  case(3): {

    switch(geomel_in->_geomel_type) {
    case(QUADR): {
      switch(order) {
      case(QQ):
        return new  FEHex27(geomel_in)  ;
      case(LL):
        return new  FEHex8(geomel_in)  ;
      case(KK):
        return new  FEHex1(geomel_in)  ;
      }
    }
    case(TRIANG): {
      switch(order) {
      case(QQ):
        return new  FETet10(geomel_in)  ;
      case(LL):
        return new  FETet4(geomel_in)  ;
      case(KK):
        return new  FETet1(geomel_in)  ;
      }
    }
    }

  } //dim 3

  default: {
    std::cout << "FEElemBase: Mistaken build" << std::endl;
    abort();
  }

  } //dim


}



//this build class allows me to return a pointer to a child of this class
//even if i am a father
//inside this class the children must be INSTANTIATED, and then they are returned.
//These instantiations are never destroyed until you explicitly delete them
//the build() function returns a POINTER



void FEElemBase::AssociateQRule(QRule* qrule_in)  {
  _qrule = qrule_in;
}
void FEElemBase::SetOrder(uint fe)  {
  _order = fe;
}
void FEElemBase::SetUtils(Utils * utils_in)  {
  _utils = utils_in;
}



void FEElemBase::init() {

  if ( _geomel->_geomel_type != QUADR && _geomel->_geomel_type != TRIANG  ) {
    std::cout << "FE::FE: GeomEl type " << _geomel->_geomel_type << " not supported" << std::endl;
    abort();
  }

  if ( _order != QQ && _order != LL && _order != KK ) {
    std::cout << "FE::FE: FE family " << _order << " not supported" << std::endl;
    abort();
  }


  std::cout << " Reading FE:: " << std::endl;

  std::string femus_dir = getenv("FEMUS_DIR");  //TODO remove this
  std::string    ext_in = DEFAULT_EXT_IN;
  std::ostringstream file;

  std::string    f_shape = "shape";
  std::string  geomel[2];
  geomel[QUADR]  =  "quadr_";
  geomel[TRIANG] = "triang_";
  uint space_dim = _geomel->_dim;

  if ( _order == QQ || _order == LL || _order == KK) {

    for (int vb=0; vb<VB; vb++) {

      _phi_mapVB[vb] = new double[_qrule->_NoGaussVB[vb]*_ndof[vb]];
      _dphidxez_mapVB[vb] = new double[_qrule->_NoGaussVB[vb]*_ndof[vb]*(space_dim-vb)];

      uint dim = space_dim - vb;
      file.str("");
      file << femus_dir << "/" << DEFAULT_CONTRIBDIR << "/" << DEFAULT_FEMDIR << "/" << geomel[_geomel->_geomel_type]
           << f_shape << dim << "D_" << _qrule->_NoGaussVB[vb] << "." << _ndof[vb] << ext_in;
      std::ifstream infile(file.str().c_str());
      readVB(vb,dim,infile);               //TODO    AAAAAAAAAAAAAAAAAAAAAAAAAA  notice that one dimension is the SPACE and another one is the MANIFOLD!!!

      
    }

  }

  else {
    std::cout << "FE::init: FE order " << _order << " not supported" << std::endl;
    abort();
  }

  return;
}




//TODO OK, in this routine we have THE SAME dimension of the ELEMENT
//as the SAME NUMBER OF DERIVATIVES
//This is because  WE ARE NOT ON A MANIFOLD
//If we were on a manifold we could consider CURVED SURFACE ELEMENTS and maybe
//we would need to distinguish between the dimension of the manifold and the dimension
//of the embedding space
//Since we are in a single abstract element, this routine will be called twice, one for VV and one for BB

//we have to fill the VV and BB arrays
//so we have to pick the 3 and 2d, or the 2 and 1d

void FEElemBase::readVB(const uint vb, const uint dim_in, std::istream& infile) {

  if (!infile) {
    std::cout << "Read: Gauss Input file "<< infile << " not opened."
              << std::endl;
    exit(3);
  }
  const int  bufLen = 256;
  char  buf[bufLen+1];
  sprintf(buf,"0");
  uint ng;
  double dummy;
  uint ngss = _qrule->_NoGaussVB[vb];
  uint nsh  = _ndof[vb];
  std::cout << " Reading elem with vb = " << vb << " , " << ngss << " gaussian points and "
            << nsh << " shape functions \n" << std::endl;
  while (strncmp(buf,"gpoints",7) != 0)    infile >> buf;
  infile >> ng;
  if ( ngss != ng) {
    std::cout << " We are reading from the wrong file" << std::endl;
    abort();
  }

  // Reading 2d shape function at gauss points
  for (uint k=0; k<ngss; k++) {
    infile >> buf;
    infile >> dummy;
//     _qrule->_weightVB[vb][k] = dummy;  //TODO AAA if you read this again it doesnt work!!! do not read this again
    for (uint s=0; s< nsh; s++) {
      infile >> dummy;
      _phi_mapVB[vb][s*ngss+k]=dummy;
      for (uint idim=0; idim< dim_in; idim++) {
        infile >> dummy;
        _dphidxez_mapVB[vb][ idim*nsh*ngss +s*ngss +k ] = dummy;
      }
    }
  }

  return;
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


void FEElemBase::init_switch() {

  if ( _geomel->_geomel_type != QUADR && _geomel->_geomel_type != TRIANG  ) {
    std::cout << "FE::FE: GeomEl type " << _geomel->_geomel_type << " not supported" << std::endl;
    abort();
  }

  if ( _order != QQ && _order != LL && _order != KK ) {
    std::cout << "FE::FE: FE family " << _order << " not supported" << std::endl;
    abort();
  }


  std::string  geomel[2];
  geomel[QUADR]  =  "quadr_";
  geomel[TRIANG] = "triang_";
  uint space_dim = _geomel->_dim;


  for (int vb=0; vb<VB; vb++) {

    uint dim = space_dim - vb;

// TODO remember to DEALLOCATE THESE
    _phi_mapVBGD[vb] = new double*[_qrule->_NoGaussVB[vb]];
    _dphidxez_mapVBGD[vb] = new double*[_qrule->_NoGaussVB[vb]];

    for (int g=0; g < _qrule->_NoGaussVB[vb]; g++) {
      _phi_mapVBGD[vb][g] = new double[_ndof[vb]];
      _dphidxez_mapVBGD[vb][g] = new double[_ndof[vb]*dim];
    }

  } //end VB



  if ( _qrule->_qrule_type == "Gauss5th") {     //only quadrature rule implemented so far
    std::string gauss_ord = "fifth";
    
    switch(_order) {  //what leads is the FE ORDER, the quadrature rule is usually mathematically chosen with respect to that

// ================================================================================
    case(QQ):  {

      switch(space_dim) {

      case(2): {

        switch(_geomel->_geomel_type)  {

        case(QUADR): {  //QUADR-2D-QQ  ========
          elem_type* myelems[VB];
          myelems[VV] = new elem_type("quad","biquadratic",gauss_ord.c_str());
          myelems[BB] = new elem_type("line","biquadratic",gauss_ord.c_str());

	  
	//here i can construct beforehand an array of function pointers for VV and BB, so that i'll loop on it later
	//the declaration of function pointers doesn't seem to me to be well understood by the compiler...  
	//he basically accepts any declaration... but then he checks when you do the ASSIGNMENT...  
	// i have to understand if i should make these functions STATIC or not  
	// it seems like you don't need STATIC
	//STATIC function DOES NOT MEAN a CONST function (a function that modifies the data)
	//   
	//i want to make an array of arrays of function pointers, EXTERNALLY STATIC [VB] and INTERNALLY DYNAMIC[ dim or (dim-1)]
	  
//        *FunctionPointer is a function.
// 	  hence, FunctionPointer is a pointer to function
// 	  you want an array of pointers to pointers to function,
// 	  so your array will be FunctionPointer*
// 	  you have to think that one star is embedded in the "Pointer" word
// 	  typedef double* (elem_type::*FunctionPointer)(const unsigned & ig) const; //declaring the FunctionPointer type
	  // I guess I have to set this definition at the class level
	  
	  //ok, allora, che cosa voglio fare:
	  //ho due istanziazioni di una classe
	  //da ciascuna istanziazione voglio prendere diversi puntatori a funzione
	  
	  
// // // 	  elem_type::FunctionPointer* Dphiptr[VB];   // OR ALSO double* (elem_type::**Dphiptr[VB])(const unsigned & ig) const;  //static two-array of function pointers 
// // // 	  Dphiptr[VV] = new elem_type::FunctionPointer[2];
// // // 	  Dphiptr[BB] = new elem_type::FunctionPointer[1];
// // // // 	  THIS ONE DOESN'T WORK double* ((elem_type::**prova))(const unsigned & ig) const  myfptr = new  double* (elem_type::*prova)(const unsigned & ig) const [2];
// // // 	  Dphiptr[VV][0] = &elem_type::GetDPhiDXi;
// // // 	  Dphiptr[VV][1] = &elem_type::GetDPhiDEta;
// // // 	  Dphiptr[BB][0] = &elem_type::GetDPhiDXi;
// // // 	  
// // // 	  
// // // 	  double* (elem_type::*DphiptrVV[2])(const unsigned & ig) const;  //static array of pointers
// // // 	  double* (elem_type::*DphiptrBB[1])(const unsigned & ig) const;

	  for (int vb=0; vb<VB; vb++) {

            if ( myelems[vb]->GetGaussPointNumber() != _qrule->_NoGaussVB[vb]) {
              std::cout << "Wrong gauss points" << std::endl;
              abort();
            }

            for (int ig = 0; ig < _qrule->_NoGaussVB[vb]; ig++) {

              for (int idof=0; idof < _ndof[vb]; idof++) {
//                 std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << vb << " " << ig << " " << idof << std::endl;
                _phi_mapVBGD[vb][ig][idof] = myelems[vb]->GetPhi(ig)[idof];
// 	std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << vb << " " << ig << " " << dof << " phi " << _phi_mapVBGD[vb][ig][dof] << std::endl;

// derivatives in canonical element
		uint dim = space_dim - vb;
		for (uint idim = 0; idim < dim; idim++) {
		  
		_dphidxez_mapVBGD[vb][ig][ idof + idim*_ndof[vb]]/* =  myelems[vb]->Dphiptr[vb][idim](ig)[idof]*/;

		}
 		
	      }
            }
          } //end VB
          

          
          break; //quadr
        }

        case(TRIANG): {




          break;
        }  //end TRIANG-2D-QQ  ======


        } //end geomel_type

        break;
      }  //end 2D
      case(3): {

        switch(_geomel->_geomel_type)  {

        case(QUADR): {  //QUADR-3D-QQ
          elem_type myelem("hex","biquadratic",gauss_ord.c_str());
          elem_type myelem_b("quad","biquadratic",gauss_ord.c_str());

          //boundary element: Quad9
          break;
        } //end //QUADR-3D-QQ

        case(TRIANG): {  //TRIANG-3D-QQ
          elem_type myelem("tet","biquadratic",gauss_ord.c_str());
          elem_type myelem_b("tri","biquadratic",gauss_ord.c_str());

          break;
        }  //end TRIANG-3D-QQ

        } //end geomel_type

        break;
      }  //end 2D
      default: {
        std::cout << "Space_dim ONE not implemented" << std::endl;
        abort();
        break;
      } //end 1D

      }  //end switch dimension

      break;
    } //end feorder QQ

// ================================================================================

// ================================================================================
    case(LL): {
      switch(space_dim) {
      case(2): {
        switch(_geomel->_geomel_type)  {
        case(QUADR): {  //QUADR-2D-LL
          elem_type myelem("quad","linear",gauss_ord.c_str());
          elem_type myelem_b("line","linear",gauss_ord.c_str());

          break;
        } //end //QUADR-2D-LL

        case(TRIANG): { //TRIANG-2D-LL
          elem_type myelem("tri","linear",gauss_ord.c_str());
          elem_type myelem_b("line","linear",gauss_ord.c_str());


          break;
        }  //end TRIANG-2D-LL

        } //end switch geomel_type
        break;
      }  //end 2D
      case(3): {
        switch(_geomel->_geomel_type)  {
        case(QUADR): { //QUADR-3D-LL
          elem_type myelem("hex","linear",gauss_ord.c_str());
          elem_type myelem_b("quad","linear",gauss_ord.c_str());

          break;
        } //end //QUADR-3D-LL

        case(TRIANG): { //TRIANG-3D-LL
          elem_type myelem("tet","linear",gauss_ord.c_str());
          elem_type myelem_b("tri","linear",gauss_ord.c_str());

          break;
        }  //end //TRIANG-3D-LL

        } //end geomel_type

        break;
      }  //end 3D

      default: {
        std::cout << "Space_dim ONE not implemented" << std::endl;
        abort();
        break;
      }

      }  //end switch dimension

      break;
    } //end feorder LL
// ================================================================================


// ================================================================================
    case(KK): {
      switch(space_dim) {
      case(2): {
        switch(_geomel->_geomel_type)  {

        case(QUADR): {  //QUADR-2D-KK
          elem_type myelem("quad","constant",gauss_ord.c_str());
          elem_type myelem_b("line","constant",gauss_ord.c_str());

        } //end //QUADR-2D-KK

        case(TRIANG): { //TRIANG-2D-KK
          std::cout << "Not implemented yet" << std::endl;
          abort();

        }  //end //TRIANG-2D-KK

        } //end geomel_type
      }  //end 2D
      case(3): {
        switch(_geomel->_geomel_type)  {
        case(QUADR): { //QUADR-3D-KK
          elem_type myelem("hex","constant",gauss_ord.c_str());
          elem_type myelem_b("quad","constant",gauss_ord.c_str());

          break;
        } //end  //QUADR-3D-KK

        case(TRIANG): {  //TRIANG-3D-KK
          std::cout << "Not implemented yet" << std::endl;
          abort();

          break;
        }  //end //TRIANG-3D-KK

        } //end geomel_type

        break;
      }  //end 2D
      default: {
        std::cout << "Space_dim ONE not implemented" << std::endl;
        abort();
        break;
      }


      }  //end switch dimension

      break;
    }  //end feorder KK
// ================================================================================

    }  //end switch fe order


  } //end if Gauss5th

  else {
    std::cout << "Quadrature rule not implemented" << std::endl;
    abort();
  }


  return;
}