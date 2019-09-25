
#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"

#include "Math.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "Domain.hpp"
#include "MultiLevelProblem.hpp"
#include "FETypeEnum.hpp"
#include "GeomElTypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "Quantity.hpp"
#include "TimeLoop.hpp"
#include "CurrentElem.hpp"

#include "FemusDefault.hpp"

#include "OptLoop.hpp"

using namespace femus;


 void GenMatRhsNSAD(MultiLevelProblem &ml_prob) {

  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_NSAD");
  const unsigned Level = my_system.GetLevelToAssemble();
  
    // ========= parameters
  const double alphaVel = ml_prob.GetInputParser().get("alphaVel");
  const double IRe      = 1./ml_prob.GetInputParser().get("Re");

  //=========== Operators 
  double dphijdx_g[DIMENSION];
  double dphiidx_g[DIMENSION];
    double    curlxiXB_g3D[3];

// //======== GEOMETRICAL ELEMENT =======
  const uint space_dim = ml_prob._ml_msh->GetDimension();

        my_system._LinSolver[Level]->_KK->zero();
        my_system._LinSolver[Level]->_RESC->zero();
  
// ==========================================  
  Mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(Level);
  elem		*myel		=  mymsh->el;
  const unsigned myproc  = mymsh->processor_id();
	
// ==========================================  
// ==========================================  
   {//BEGIN VOLUME    
   
   const uint mesh_vb = VV;

    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level];
    
  for (uint iel=0; iel < (nel_e - nel_b); iel++) {
    
    CurrentElem<double>       currelem(iel,myproc,Level,VV,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType(),mymsh);
    //   CurrentGaussPointBase & currgp = //   CurrentGaussPointBase::build(currelem,ml_prob.GetQuadratureRule(currelem.GetDim()));
   
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
      CurrentQuantity VelAdjOldX(currelem);
    VelAdjOldX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_VelocityAdj0"); 
    VelAdjOldX.VectWithQtyFillBasic();
    VelAdjOldX.Allocate();
    
    CurrentQuantity VelAdjOldY(currelem);
    VelAdjOldY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_VelocityAdj1"); 
    VelAdjOldY.VectWithQtyFillBasic();
    VelAdjOldY.Allocate();
    
    CurrentQuantity VelAdjOldZ(currelem);
    VelAdjOldZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_VelocityAdj2"); 
    VelAdjOldZ.VectWithQtyFillBasic();
    VelAdjOldZ.Allocate();
    
    std::vector<CurrentQuantity*> VelAdjOld_vec;
    VelAdjOld_vec.push_back(&VelAdjOldX);
    VelAdjOld_vec.push_back(&VelAdjOldY);
    VelAdjOld_vec.push_back(&VelAdjOldZ);

    CurrentQuantity PressAdjOld(currelem);
    PressAdjOld._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_PressureAdj");
    PressAdjOld.VectWithQtyFillBasic();
    PressAdjOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================
  
//=========EXTERNAL QUANTITIES (couplings) =====

//========= DOMAIN MAPPING
    CurrentQuantity xyz(currelem);
    xyz._dim      = space_dim;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

    CurrentQuantity VelX(currelem);
    VelX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity0"); 
    VelX.VectWithQtyFillBasic();
    VelX.Allocate();
    
    CurrentQuantity VelY(currelem);
    VelY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity1"); 
    VelY.VectWithQtyFillBasic();
    VelY.Allocate();
    
    CurrentQuantity VelZ(currelem);
    VelZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity2"); 
    VelZ.VectWithQtyFillBasic();
    VelZ.Allocate();
    
    std::vector<CurrentQuantity*> Vel_vec;   
    Vel_vec.push_back(&VelX);
    Vel_vec.push_back(&VelY);
    Vel_vec.push_back(&VelZ);


    
    
  CurrentQuantity VelDesX(currelem);
    VelDesX._dim        = 1;
    VelDesX._FEord      = VelX._FEord;
    VelDesX._ndof       = VelX._ndof;
    VelDesX.Allocate();

  CurrentQuantity VelDesY(currelem);
    VelDesY._dim        = 1;
    VelDesY._FEord      = VelX._FEord;
    VelDesY._ndof       = VelX._ndof;
    VelDesY.Allocate();

  CurrentQuantity VelDesZ(currelem);
    VelDesZ._dim        = 1;
    VelDesZ._FEord      = VelX._FEord;
    VelDesZ._ndof       = VelX._ndof;
    VelDesZ.Allocate();

  std::vector<CurrentQuantity*> VelDes_vec;   
    VelDes_vec.push_back(&VelDesX);
    VelDes_vec.push_back(&VelDesY);
    VelDes_vec.push_back(&VelDesZ);
    
    
    
    CurrentQuantity BhomX(currelem);
    BhomX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom0"); 
    BhomX.VectWithQtyFillBasic();
    BhomX.Allocate();
    
    CurrentQuantity BhomY(currelem);
    BhomY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom1"); 
    BhomY.VectWithQtyFillBasic();
    BhomY.Allocate();
    
    CurrentQuantity BhomZ(currelem);
    BhomZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom2"); 
    BhomZ.VectWithQtyFillBasic();
    BhomZ.Allocate();
    
    std::vector<CurrentQuantity*> Bhom_vec;   
    Bhom_vec.push_back(&BhomX);
    Bhom_vec.push_back(&BhomY);
    Bhom_vec.push_back(&BhomZ);
 
    //Bhom only to retrieve the dofs

    CurrentQuantity BextX(currelem);
    BextX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt0"); 
    BextX.VectWithQtyFillBasic();
    BextX.Allocate();
    
    CurrentQuantity BextY(currelem);
    BextY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt1"); 
    BextY.VectWithQtyFillBasic();
    BextY.Allocate();
    
    CurrentQuantity BextZ(currelem);
    BextZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt2"); 
    BextZ.VectWithQtyFillBasic();
    BextZ.Allocate();
    
    std::vector<CurrentQuantity*> Bext_vec;   
    Bext_vec.push_back(&BextX);
    Bext_vec.push_back(&BextY);
    Bext_vec.push_back(&BextZ);
 
//========= auxiliary, must be AFTER Bhom!   //TODO this doesnt have any associated quantity!
  CurrentQuantity Bmag(currelem); //total
    Bmag._dim        = Bhom_vec.size();
    Bmag._FEord      = BhomX._FEord;
    Bmag._ndof       = BhomX._ndof;
    Bmag.Allocate();
    
//===============
    CurrentQuantity BhomAdjX(currelem);
    BhomAdjX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj0"); 
    BhomAdjX.VectWithQtyFillBasic();
    BhomAdjX.Allocate();
    
    CurrentQuantity BhomAdjY(currelem);
    BhomAdjY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj1"); 
    BhomAdjY.VectWithQtyFillBasic();
    BhomAdjY.Allocate();
    
    CurrentQuantity BhomAdjZ(currelem);
    BhomAdjZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj2"); 
    BhomAdjZ.VectWithQtyFillBasic();
    BhomAdjZ.Allocate();
    
    std::vector<CurrentQuantity*> BhomAdj_vec;
    BhomAdj_vec.push_back(&BhomAdjX);
    BhomAdj_vec.push_back(&BhomAdjY);
    BhomAdj_vec.push_back(&BhomAdjZ);
    
    CurrentQuantity BhomAdj_vecQuant(currelem);   //without quantity, nor equation
    BhomAdj_vecQuant._dim     = BhomAdj_vec.size();
    BhomAdj_vecQuant._FEord   = BhomAdj_vec[0]->_FEord;
    BhomAdj_vecQuant._ndof    = BhomAdj_vec[0]->_ndof;
    BhomAdj_vecQuant.Allocate();
//========= END EXTERNAL QUANTITIES =================


    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.SetDofobjConnCoords();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);

    currelem.SetElDofsBc();
    

    
 for (uint idim=0; idim < space_dim; idim++)    {
    VelAdjOld_vec[idim]->GetElemDofs();  
          Vel_vec[idim]->GetElemDofs();
      BhomAdj_vec[idim]->GetElemDofs();
         Bhom_vec[idim]->GetElemDofs();
         Bext_vec[idim]->GetElemDofs();
    VelDesired(&ml_prob,*VelDes_vec[idim],currelem,idim);
 }
    
    
    BhomAdj_vecQuant.GetElemDofs(BhomAdj_vec);  //this must be AFTER filling the scalar components
    PressAdjOld.GetElemDofs();

 
//======SUM Bhom and Bext  //from now on, you'll only use Bmag //Bmag,Bext and Bhom must have the same orders!
    Math::zeroN(&Bmag._val_dofs[0],Bmag._dim*Bmag._ndof);

    for (uint ivarq=0; ivarq < Bmag._dim; ivarq++)    { //ivarq is like idim
          for (uint d=0; d < Bmag._ndof; d++)    {
          const uint     indxq  =         d + ivarq*Bmag._ndof;
          Bmag._val_dofs[indxq] = Bext_vec[ivarq]->_val_dofs[d] + Bhom_vec[ivarq]->_val_dofs[d];
	  }
    }    
//=======    
    
///optimal control
    xyz.SetElemAverage();
  int el_flagdom = ElFlagControl(xyz._el_average,ml_prob._ml_msh);


   const uint el_ngauss = ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();

   for (uint qp = 0; qp < el_ngauss; qp++) {
//=======here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     { 
//   currgp.SetPhiElDofsFEVB_g (fe,qp);  
//   currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}  
	  
const double      det = 1.; // currgp.JacVectVV_g(xyz);
const double dtxJxW_g = det*ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { 
//   currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp);
//   currgp.ExtendDphiDxyzElDofsFEVB_g(fe); 
}
//=======end of the "COMMON SHAPE PART"==================

     BhomAdj_vecQuant.curl_g();
     Bmag.val_g();
 
  for (uint idim=0; idim < space_dim; idim++) {
    VelAdjOld_vec[idim]->val_g();
    VelDes_vec[idim]->val_g();
    Vel_vec[idim]->val_g();
    Vel_vec[idim]->grad_g();
  }
  
//vector product
        Math::extend(&Bmag._val_g[0],&Bmag._val_g3D[0],space_dim);
        Math::cross(&BhomAdj_vecQuant._curl_g3D[0],&Bmag._val_g3D[0],curlxiXB_g3D);

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
// // //        for (uint i = 0; i < VelAdjOldX._ndof; i++)     {
// // // 
// // // 	const double phii_g = currgp._phi_ndsQLVB_g[VelAdjOldX._FEord][i];
// // //         for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim]= currgp._dphidxyz_ndsQLVB_g[VelAdjOldX._FEord][i+idim*VelAdjOldX._ndof];
// // // 
// // //          for (uint idim=0; idim<space_dim; idim++) {
// // //             const uint irowq=i+idim*VelAdjOldX._ndof;  //quadratic rows index
// // //            currelem.Rhs()(irowq) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                            - curlxiXB_g3D[idim]*phii_g                    //this is due to the variation of velocity in the MAGNETIC ADVECTION, so it is due to a   NONLINEAR COUPLING "u times B", "MY_STATE times OTHER_STATE"
// // //                            - alphaVel*el_flagdom*(Vel_vec[idim]->_val_g[0] - VelDes_vec[idim]->_val_g[0])*phii_g    //this is the dependence that counts
// // //                            )
// // //                            + (1-currelem.GetBCDofFlag()[irowq])*detb*VelAdjOld_vec[idim]->_val_dofs[i]; //Dirichlet bc
// // // 	   }
// // // 
// // //   for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
// // //           const uint irowq = i+idim*VelAdjOldX._ndof;
// // //           currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
// // //         }// end filling diagonal for Dirichlet bc
// // // 	 
// // //         for (uint j=0; j<VelAdjOldX._ndof; j++) {// A element matrix
// // // 	  
// // //            double                                  phij_g       =      currgp._phi_ndsQLVB_g[VelAdjOldX._FEord][j];
// // //            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[VelAdjOldX._FEord][j+idim*VelAdjOldX._ndof];
// // // 
// // //           double     Lap_g = Math::dot(dphijdx_g,dphiidx_g,space_dim);
// // // 	  double Advphii_g = 0.;
// // // 	  for (uint idim=0; idim<space_dim; idim++) Advphii_g += Vel_vec[idim]->_val_g[0]*dphiidx_g[idim];   //TODO can put it outside
// // //           
// // //           for (uint idim=0; idim<space_dim; idim++) { //filled in as 1-2-3 // 4-5-6 // 7-8-9
// // //             int irowq = i+idim*VelAdjOldX._ndof;
// // //             // diagonal blocks [1-5-9]
// // //             currelem.Mat()(irowq,j+idim*VelAdjOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                    + IRe*(dphijdx_g[idim]*dphiidx_g[idim] + Lap_g)      //Adjoint of D is D, Adjoint of Laplacian is Laplacian
// // //                    + phij_g*phii_g*/*dveldx_g*/Vel_vec[idim]->_grad_g[0][idim]  //Adjoint of Advection 1 delta(u) DOT grad(u): adj of nonlinear stuff has 2 TERMS (well, not always)
// // //                    + phij_g*Advphii_g                                   //Adjoint of Advection 2 u DOT grad (delta(u)): adj of nonlinear stuff has 2 TERMS
// // //                );
// // //             // block +1 [2-6-7]
// // //             int idimp1=(idim+1)%space_dim;
// // //             currelem.Mat()(irowq,j+idimp1*VelAdjOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                    + IRe*(dphijdx_g[idim]*dphiidx_g[idimp1])
// // //                    + phij_g*phii_g*/*dveldx_g*/Vel_vec[idimp1]->_grad_g[0][idim]
// // //                );
// // // #if (DIMENSION==3)
// // //             // block +2 [3-4-8]
// // //             int idimp2=(idim+2)%space_dim;
// // //             currelem.Mat()(irowq,j+idimp2*VelAdjOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                   + IRe*(dphijdx_g[idim]*dphiidx_g[idimp2])
// // //                   + phij_g*phii_g*/*dveldx_g*/Vel_vec[idimp2]->_grad_g[0][idim]
// // //                );
// // // #endif
// // //           }
// // //  
// // //         } 
// // //                                       // end A element matrix
// // //       
// // //        for (uint j=0; j<PressAdjOld._ndof; j++) {// B^T element matrix ( p*div(v) ) (NS eq)
// // //           const double psij_g =  currgp._phi_ndsQLVB_g[PressAdjOld._FEord][j];
// // //           const int jclml= j + space_dim*VelAdjOldX._ndof;
// // //           for (uint idim=0; idim<space_dim; idim++) {
// // //             uint irowq = i+idim*VelAdjOldX._ndof;
// // //             currelem.Mat()(irowq,jclml) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(-psij_g*dphiidx_g[idim]);
// // //            }
// // //         }
// // //                                      // end B^T element matrix
// // // 
// // //           if (i<PressAdjOld._ndof) {//  pressure equation (KOMP dp/dt=rho*div) 
// // //           double psii_g = currgp._phi_ndsQLVB_g[PressAdjOld._FEord][i];
// // // 	  const uint irowl = i+space_dim*VelAdjOldX._ndof;  //vertical offset
// // //           currelem.Rhs()(irowl)=0.;  // rhs
// // // 
// // //           for (uint j=0; j<VelAdjOldX._ndof; j++) { // B element matrix q*div(u)
// // //             for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[VelAdjOldX._FEord][j+idim*VelAdjOldX._ndof];
// // //             for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*VelAdjOldX._ndof) += -dtxJxW_g*psii_g*dphijdx_g[idim]; 
// // //                 }
// // // 
// // //         }
// // //                          // end pressure eq (cont)
// // //       }
//===================================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) =====================
//===================================================================
      
    }
    // end element gaussian integration loop
    
    ///  Add element matrix and rhs to the global ones.
    my_system._LinSolver[Level]->_KK->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    my_system._LinSolver[Level]->_RESC->add_vector(currelem.Rhs(),currelem.GetDofIndices());
    
  } 
  // end of element loop

 
  }//END VOLUME
  
        my_system._LinSolver[Level]->_KK->close();
        my_system._LinSolver[Level]->_RESC->close();
  
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << my_system.name() << ": assembled  Level " << Level
           << " with " << my_system._LinSolver[Level]->_KK->m() << " dofs" << std::endl;
#endif    

return;
}

