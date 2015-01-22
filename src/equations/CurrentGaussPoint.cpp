#include "CurrentGaussPoint.hpp"

#include "MultiLevelProblemTwo.hpp"
#include "Math.hpp"
#include "CurrentQuantity.hpp"
#include "MeshTwo.hpp"
#include "GeomEl.hpp"

#include <cmath>


namespace femus {



//I need to hold the equations map pointer, because i also need a mesh pointer
// for the geometric element
//maybe later on i'd just pass the GeomElement(GeomEl) and the MathElement(FE)
//by the way, with the MultiLevelProblemTwo I reach the Utils, the Mesh, and so the GeomEl, and so on...
template <unsigned int FM_DIM>
CurrentGaussPoint<FM_DIM>::CurrentGaussPoint(const CurrentElem & curr_el_in, MultiLevelProblemTwo& e_map_in ): 
        CurrentGaussPointBase(curr_el_in,e_map_in) {
  
 
}


template <unsigned int FM_DIM>
CurrentGaussPoint<FM_DIM>::~CurrentGaussPoint() {
  
//     for (int i=0;i< VB;i++){
//       for (int j=0;j< QL;j++) {
// 	delete [] _phi_ndsQLVB_g[i][j];
// 	delete []  _dphidxyz_ndsQLVB_g[i][j];
// 	delete []  _dphidxyz_ndsQLVB_g3D[i][j];
// 	delete [] _dphidxezeta_ndsQLVB_g[i][j];
//       }
//     }
   
  
}



template <unsigned int FM_DIM>
void CurrentGaussPoint<FM_DIM>::SetPhiElDofsFEVB_g(const uint qlflag, const uint qp) {
  
    const uint el_nnodes =  _elem_type[qlflag]->GetNDofs();
    const uint el_ngauss =  _qrule.GetGaussPointsNumber();
   
   
    for (uint eln=0; eln<el_nnodes; eln++)    { 
              uint lqp=eln*el_ngauss+qp;
              _phi_ndsQLVB_g[qlflag][eln] = _elem_type[qlflag]->GetPhi(qp,eln);
         }

return;

}


//TODO this is DIFFERENT on the boundary!
//The single derivatives of the domain transformation are NOT in an INVJac
//because the transformation is not identified by a square matrix,
//but a matrix dimension x bdry_dimension
//well, actually we can do an _InvJac_g[vb]
// _InvJac_g[VV] is dimensionxdimension
// _InvJac_g[BB] is dimensionxbdry_dimension
//when computing the absolute derivative, starting from the canonical derivatives,
//you must loop over the CANONICAL variables of the element, 
//so over 'dimension' variables or 'bdry_dimension' variables
//then the rest should be alright
 template <unsigned int FM_DIM>
void  CurrentGaussPoint<FM_DIM>::SetDPhiDxyzElDofsFEVB_g(const uint qlflag, const uint qp) {
    
    
  const uint ndim      = _current_elem.GetDim();
  const uint el_nnodes = _elem_type[qlflag]->GetNDofs();
  const uint el_ngauss =     _qrule.GetGaussPointsNumber();
  const uint goffset   = el_nnodes*el_ngauss;

  std::vector<double> dphidxi_g(ndim);  //the dimension of this should be _IntDim[vbflag], but we make it static like this, it doesnt hurt

         for (uint eln=0; eln<el_nnodes; eln++)    { 
              uint lqp=eln*el_ngauss+qp;
  
           for (uint idim=0; idim<ndim; idim++)  dphidxi_g[idim] = _elem_type[qlflag]->GetDPhiDxez(qp, eln + idim*el_nnodes);    /*->_dphidxez_mapVB[vbflag][lqp+idim*goffset];*/
	    
	    for (uint idim=0; idim<ndim; idim++) {
	    double sum = 0.;
	     for (uint jdim=0; jdim<ndim; jdim++) sum += _InvJac_g[jdim][idim]*dphidxi_g[jdim];
	        _dphidxyz_ndsQLVB_g[qlflag][eln+idim*el_nnodes] = sum;
	  }
  
  
	 }

//CHECK in NS
//   int eln=0;
//   std::cout << "%%%%%%%%%%%% dphidcsi_g for node " << eln << " for qp " << qp << std::endl;
//   for (uint idim=0; idim<ndim; idim++)  std::cout <<  _dphidxezeta_ndsQLVB_g[vb][FE_VEL][eln+idim*el_ndof_q] << " ";
//   std::cout << "%%%%%%%%%%%% " << std::endl;
  

	 
  return;
 }





//=====================
//I dont need qp here
//also, i only do this function AFTER FILLING _dphidxyz_ndsQLVB_g!
template <unsigned int FM_DIM>
void  CurrentGaussPoint<FM_DIM>::ExtendDphiDxyzElDofsFEVB_g(const uint qlflag/*, const uint qp*/) {

  //AAA: valid from dimension to 3

  const uint ndim     = _current_elem.GetDim();
  const uint el_ndofs = _elem_type[qlflag]->GetNDofs();
  
//set to zero  
   for (uint eln=0; eln<el_ndofs; eln++)  {
      for (uint i=0; i<3; i++) {
             _dphidxyz_ndsQLVB_g3D[qlflag][eln+i*el_ndofs]=0.; 
            }
       }
//extend
   for (uint eln=0; eln<el_ndofs; eln++)    { 
       for (uint idim=0; idim<ndim; idim++) {
                 _dphidxyz_ndsQLVB_g3D[qlflag][eln+idim*el_ndofs] = _dphidxyz_ndsQLVB_g[qlflag][eln+idim*el_ndofs];
	         }
	      }
  
  
  
 return; 
}
 

//canonical derivatives
template <unsigned int FM_DIM>
void  CurrentGaussPoint<FM_DIM>::SetDPhiDxezetaElDofsFEVB_g(const uint qlflag, const uint qp) {
    
           const uint ndim = _current_elem.GetDim();
         const uint elndof = _elem_type[qlflag]->GetNDofs();
    const uint   el_ngauss =     _qrule.GetGaussPointsNumber();
      const uint   goffset = elndof*el_ngauss;


    for (uint idim=0; idim < ndim; idim++) {

         for (uint eln=0; eln<elndof; eln++)    { 
              uint lqp=eln*el_ngauss+qp;
  
	         _dphidxezeta_ndsQLVB_g[qlflag][eln+idim*elndof] = _elem_type[qlflag]->GetDPhiDxez(qp,eln + idim*elndof);  
	   }
	 }

//here you set _is_ready_for_Jac
//_dphidxez_mapVB is long 27x3x27:    for every dimension(direction derivative)
//                                    and every shape function 
//                                    and every gauss point
//the first problem that happens to us is that
//the structure of the InvJac is not the same as the element orientation, 
//while it must be like that... it seems like the columns are shifted by one... no ...
//so, for every given qp, we must pick 27x3 values:
// all the shape derivatives in the xi direction, in the canonical connectivity order
//then all the shape derivatives in the eta direction, in the canonical connectivity order
//finally all the shape derivatives in the zeta direction, ,in the canonical connectivity order
// and the dphimap in FE how is it filled? 
//It is dphimap [eln*el_ngauss + qp + idim*goffset] = dphimap [eln*el_ngauss + idim*elndof*el_ngauss + qp ]
//    = dphimap [(eln+idim*elndof)*el_ngauss + qp ] = dphimap [(eln+idim*elndof)*el_ngauss + qp ]
// so, for every gauss point, you actually jump a lot...
//in fact you run over eln and idim
//So, in a SINGLE CALL of qp, I must fill 27x3=81 values, running over BOTH "eln" and "idim"
//Changing qp only means shifting by ONE all the runs
//But the problem is: did we correctly READ the VALUES from the FILES?
//so let us go back to the files.
//In the files the values for the 27 gauss points are written.
//first, the gauss points are separated.
//then, does it matter in what order the gauss points are read?
//I dont think so: THE GAUSS ORDER is NOT IMPORTANT, just FIX ONE and do everything with that one FIXED.


  return;
 }


//this routine computes the surface jacobian and the normal and tangent vectors at the gauss point qp

//the boundary mesh in twoD runs in counterclockwise sense, therefore the tangent runs in counterclockwise sense      
//the normal tangent basis we construct here is RIGHT HANDED,
//made of (tg,norm) or (tg1,tg2,norm)
//in all cases the normal is OUTGOING from the domain
//this happens only THANKS to the way in which the BOUNDARY NODES 
//are ORDERED in the boundary mesh.
//They are ordered in counterclockwise manner,
// watching from OUTSIDE of the domain.


// *****************************************
// SPECIALIZATION of SUR JAC for DIM = 2
// *****************************************
template <>
void CurrentGaussPoint<2>::ComputeJacBB() {
  
  const uint dim_in = 2;
  
  
  _normal_g[0] =   _dxyzdxieta_g[0][1];// dydxi ;
  _normal_g[1] = - _dxyzdxieta_g[0][0];// -dxdxi ; //check these signs

  _tangent_g[0][0] = _dxyzdxieta_g[0][0]; //dxdxi
  _tangent_g[0][1] = _dxyzdxieta_g[0][1]; //dydxi

double fact = sqrt(Math::dot(_tangent_g[0],_tangent_g[0],dim_in)); //the vector length is dimension
Math::normalize(_tangent_g[0],fact,dim_in);
  
  return;
}


// *****************************************
// SPECIALIZATION of SUR JAC for DIM = 3
// *****************************************
template <>
void CurrentGaussPoint<3>::ComputeJacBB() {
  
  const uint dim_in = 3;
  
   _tangent_g[0][0] = _dxyzdxieta_g[0][0]; //dxdxi
   _tangent_g[0][1] = _dxyzdxieta_g[0][1]; //dydxi
   _tangent_g[0][2] = _dxyzdxieta_g[0][2]; //dzdxi

   _tangent_g[1][0] = _dxyzdxieta_g[1][0]; //dxdeta
   _tangent_g[1][1] = _dxyzdxieta_g[1][1]; //dydeta
   _tangent_g[1][2] = _dxyzdxieta_g[1][2]; //dzdeta

   double fact0 = sqrt(Math::dot(_tangent_g[0],_tangent_g[0],dim_in)); //the vector length is dimension
   double fact1 = sqrt(Math::dot(_tangent_g[1],_tangent_g[1],dim_in)); //the vector length is dimension
   Math::normalize(_tangent_g[0],fact0,dim_in);
   Math::normalize(_tangent_g[1],fact1,dim_in);

//if i normalize the tangent first then i can use it to compute the normal without normalizing the normal again
//TODO WRONG!!! if the tangents are normalized then the normal need not be normalized!!!
//this happens only if the two tangents are ORTHOGONAL to ONE ANOTHER!
//So i need to normalize the normal!!! (in fact i did that)
//so, beware that if you have a non straight boundary element in 3D
//then the two tangent vectors need not be ORTHOGONAL TO EACH OTHER!!!

   Math::cross(_dxyzdxieta_g[0],_dxyzdxieta_g[1],_normal_g);

///print the tangent vectors 
//clearly for the computations in the matrix i'll need its components WRT the PHYSICAL FRAME (x,y,z)
// 	for (uint j=0; j<bdry_dim; j++){
// std::cout << "=====Tangent vector " << j << " ";
// 	  for (uint i=0; i<dimension; i++){
// std::cout << _tangent_g[j][i] << " ";
// 	}
// }
// std::cout << std::endl;
  
  return;
  
}




template <unsigned int FM_DIM>
double CurrentGaussPoint<FM_DIM>::JacVectBB_g(CurrentQuantity& xyz )/* const*/ {
  
  //here you check assert(_is_ready_for_Jac);

//     const uint spacedim = _eqnmap._mesh._dim;

    const uint    Order = xyz._FEord;  //order of the coordinate transformation 
    const uint     xoff = xyz._ndof;

    //     double dxyzdxieta_g[FM_DIM - 1][FM_DIM]; 
    //TODO TODO TODO NOW I'LL PUT IT  CLASS MEMBER, (IT SHOULD STAY TOGETHER WITH _normal, _tangent and _jacobian...) 
    // I don't know if it wins to have the STATIC ALLOCATION outside or not...
    // because now that this is a CLASS VARIABLE, you have to imagine that there is like a THIS pointer in front of every occurrence,
    // which is like a FURTHER SUM every time!!!
    
    
    for (uint i=0; i< FM_DIM - 1; i++){
         for (uint j=0; j< FM_DIM; j++){ _dxyzdxieta_g[i][j]=0.; } }

      for (uint i=0; i< FM_DIM; i++){
	for (uint j=0; j<FM_DIM - 1; j++){
          for (uint s=0;s<xoff;s++) {
	  _dxyzdxieta_g[j][i] += xyz._val_dofs[s +i*xoff]*_dphidxezeta_ndsQLVB_g[Order][s+j*xoff];
	  }
	}
      }

      
    double JacSur=0.;

    ComputeJacBB();

    
        JacSur = sqrt(Math::dot(_normal_g,_normal_g,FM_DIM));
	Math::normalize(_normal_g,JacSur,FM_DIM);
   
   return JacSur;

}

// *****************************************
// SPECIALIZATION of JAC for DIM = 2
// *****************************************
template <>
double CurrentGaussPoint<2>::ComputeJacVV() {

double  det=0.;

  //here you could do an extended 3D matrix and use the general formula
  //now, we do like that to make things faster
  //we should have done a 3D matrix with 1 in the south east position and zeros 
  
  det = _dxyzdxezeta_g[0][0]*_dxyzdxezeta_g[1][1] - _dxyzdxezeta_g[0][1]*_dxyzdxezeta_g[1][0]; 
    double invdet=1./det;

  _InvJac_g[0][0] =   invdet*_dxyzdxezeta_g[1][1];
  _InvJac_g[1][1] =   invdet*_dxyzdxezeta_g[0][0];
  _InvJac_g[0][1] = - invdet*_dxyzdxezeta_g[0][1];
  _InvJac_g[1][0] = - invdet*_dxyzdxezeta_g[1][0];
  //much faster than the whole formula

//     _InvJac[0]=_InvJac_g[0][0];    // dxi dx
//   _InvJac[1]=_InvJac_g[1][0];    //deta dx
//   _InvJac[2]= _InvJac_g[0][1];   // dxi dy
//   _InvJac[3]= _InvJac_g[1][1];     //deta dy
 

return det;  
  
}

// *****************************************
// SPECIALIZATION of JAC for DIM = 3
// *****************************************
template <>
double CurrentGaussPoint<3>::ComputeJacVV() {

  const uint dim_in = 3;
  
double  det = 0.;

  int row=0;
  int rowp1=(row+1)%3;
  int rowp2=(row+2)%3;
  
  for (uint j=0; j< dim_in; j++) {
    int jdimp1=(j+1)%3;
    int jdimp2=(j+2)%3;
    det += _dxyzdxezeta_g[row][j]*(  _dxyzdxezeta_g[rowp1][jdimp1]*_dxyzdxezeta_g[rowp2][jdimp2] 
                                   - _dxyzdxezeta_g[rowp1][jdimp2]*_dxyzdxezeta_g[rowp2][jdimp1]  );
    
  }
  
    double invdet=1./det;

/////// TODO check where this is wrong    
//    for (uint j=0;j<dimension;j++) {
//       int jp1= (j+1)%3;
//       int jp2= (j+2)%3;
//     for (uint i=0;i<dimension;i++) {
//       int ip1= (i+1)%3;
//       int ip2= (i+2)%3;
//       int coeff=1-2*((i+j+2)%2);  //there must be an error in this one (-1)^(i+j)
//         _InvJac_g[j][i]  =   dxyzdxezeta_g[ip1][jp1]*dxyzdxezeta_g[ip2][jp2] 
//                            - dxyzdxezeta_g[ip1][jp2]*dxyzdxezeta_g[ip2][jp1];
//         _InvJac_g[j][i] *= coeff*invdet;
// //see if sthg can be optimised
//       }
//     }
/////// TODO check where this is wrong    

//  
 //the old index spans the matrix _InvJac column-wise
//  _InvJac[0] = _InvJac_g[0][0];                //  (y_eta *z_zeta -y_zeta*z_eta) *invdet; //was ok
//  _InvJac[3]=  _InvJac_g[0][1];                // -(x_eta *z_zeta -x_zeta*z_eta) *invdet; //was wrong
//  _InvJac[6]=  _InvJac_g[0][2];                //  (y_zeta*x_eta  -x_zeta*y_eta) *invdet; //was ok
//  _InvJac[1]=  _InvJac_g[1][0];                // -(y_xi  *z_zeta -y_zeta*z_xi)  *invdet;  //was wrong
//  _InvJac[4]=  _InvJac_g[1][1];                //  (x_xi  *z_zeta -x_zeta*z_xi)  *invdet; //was ok
//  _InvJac[7]=  _InvJac_g[1][2];                // -(y_zeta*x_xi   -y_xi  *x_zeta)*invdet;//was wrong
//  _InvJac[2]=  _InvJac_g[2][0];                //  (z_eta *y_xi   -y_eta *z_xi)  *invdet; //was ok
// _InvJac[5]=  _InvJac_g[2][1];                //  -(x_xi  *z_eta  -x_eta *z_xi)  *invdet;//was wrong
//  _InvJac[8]=  _InvJac_g[2][2];                //  (x_xi  *y_eta  -y_xi  *x_eta) *invdet; //was ok

////////////// BY HAND first, and try to understand what is the wrong order
// i think that the dphidcsi order is CORRECT 
const uint XXX=0; const uint CSI=0;
const uint YYY=1; const uint ETA=1;
const uint ZZZ=2; const uint ZETA=2;

   _InvJac_g[0][0] = (  _dxyzdxezeta_g[YYY][ETA]*_dxyzdxezeta_g[ZZZ][ZETA] - _dxyzdxezeta_g[YYY][ZETA] * _dxyzdxezeta_g[ZZZ][ETA])*invdet;
   _InvJac_g[0][1] = -( _dxyzdxezeta_g[XXX][ETA]*_dxyzdxezeta_g[ZZZ][ZETA] - _dxyzdxezeta_g[XXX][ZETA] * _dxyzdxezeta_g[ZZZ][ETA])*invdet;
   _InvJac_g[0][2] = (  _dxyzdxezeta_g[YYY][ZETA]*_dxyzdxezeta_g[XXX][ETA] - _dxyzdxezeta_g[XXX][ZETA] * _dxyzdxezeta_g[YYY][ETA])*invdet;
   _InvJac_g[1][0] = -( _dxyzdxezeta_g[YYY][CSI]*_dxyzdxezeta_g[ZZZ][ZETA] - _dxyzdxezeta_g[YYY][ZETA] * _dxyzdxezeta_g[ZZZ][CSI])*invdet;
   _InvJac_g[1][1] = (  _dxyzdxezeta_g[XXX][CSI]*_dxyzdxezeta_g[ZZZ][ZETA] - _dxyzdxezeta_g[XXX][ZETA] * _dxyzdxezeta_g[ZZZ][CSI])*invdet;
   _InvJac_g[1][2] = -( _dxyzdxezeta_g[YYY][ZETA]*_dxyzdxezeta_g[XXX][CSI] - _dxyzdxezeta_g[YYY][CSI ] * _dxyzdxezeta_g[XXX][ZETA])*invdet;
   _InvJac_g[2][0] = (  _dxyzdxezeta_g[ZZZ][ETA]*_dxyzdxezeta_g[YYY][CSI ] - _dxyzdxezeta_g[YYY][ ETA] * _dxyzdxezeta_g[ZZZ][CSI])*invdet;
   _InvJac_g[2][1] = -( _dxyzdxezeta_g[XXX][CSI]*_dxyzdxezeta_g[ZZZ][ ETA] - _dxyzdxezeta_g[XXX][ ETA] * _dxyzdxezeta_g[ZZZ][CSI])*invdet; 
   _InvJac_g[2][2] = (  _dxyzdxezeta_g[XXX][CSI]*_dxyzdxezeta_g[YYY][ ETA] - _dxyzdxezeta_g[YYY][CSI ] * _dxyzdxezeta_g[XXX][ETA])*invdet;
///////////END BY HAND
   
///print of InvJac, very useful for debugging
// // //       std::cout << "******* InvJac for qp " << qp << std::endl;
// // //    for (uint j=0;j<dimension;j++) {
// // //         for (uint i=0;i<dimension;i++) {
// // // 	std::cout <<  _InvJac_g[j][i] << " ";
// // // 	  
// // // 	}
// // // 	 std::cout << std::endl;
// // //    }
// // //          std::cout << "******* END InvJac for qp " << qp << std::endl;


return det;  
  
  
}



//We should make no distinction and put everything in the same part
//this function must be called after the CANONICAL DERIVATIVES and the DOMAIN
//it fills _InvJac_g: NOTICE that THIS IS UNIQUE, there is no distinction between LINEAR and QUADRATIC.

template <unsigned int FM_DIM>
double CurrentGaussPoint<FM_DIM>::JacVectVV_g(CurrentQuantity& xyz )/* const*/ {

const uint Order    = xyz._FEord;  //order of the coordinate transformation
const uint xoff     = xyz._ndof;
  
  
      for (uint i=0; i< FM_DIM; i++){
         for (uint j=0; j< FM_DIM; j++){ _dxyzdxezeta_g[i][j]=0.; } }

      for (uint i=0; i< FM_DIM; i++){
	for (uint j=0; j< FM_DIM; j++){
          for (uint s=0;s<(uint)xoff;s++) {
	  _dxyzdxezeta_g[i][j] += xyz._val_dofs[s + i * xoff]*_dphidxezeta_ndsQLVB_g[Order][s+j*xoff];
	  }
	}
      }
      
 // _dphidxezeta_ndsQLVB_g[vb][Order][s+j*(elnshape)]
// ok, so I must be sure that there is a correspondence between the nodes and the dphidxi of the CORRESPONDING nodes at the GAUSS POINT qp
// first, CHECK the NODE CORRESPONDENCE
//second, check the GAUSS CORRESPONDENCE WITH WHAT YOU MULTIPLY in the following
//_dphidxezeta_ndsQLVB_g[vb][Order]
//let us see what is stored in this pointer. For a given Order and vb,
//we store a number of values given by NSHAPExdimension (= 27x3)
//so the function is f[Order][vb][ NSHAPE(Order,vb) x dimension ]
//All of this happens for a given gauss point qp, that is not even used any longer 
// in this routine, it is fixed from before.
//So now the point is: the _dphidxezeta_ndsQLVB_g vector (= 27x3)
// is filled starting from a  vector which is 27x3x27
//and every time,for a varying qp, I get 27x3 different points
//so let us see how i fill this thing

  
  return ComputeJacVV();
  
}

//***************************************
//explicit instantiations for 2D and 3D 
// (this way i restrict the set of available template values to the ones i ask for an explicit instantiation)
//****************************************
template class CurrentGaussPoint<2>;
template class CurrentGaussPoint<3>;




} //end namespace femus


