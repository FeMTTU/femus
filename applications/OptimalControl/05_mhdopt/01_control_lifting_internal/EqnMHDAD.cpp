
#include "FemusDefault.hpp"

#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"

#include "Math.hpp"
#include "Quantity.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "MultiLevelProblem.hpp"
#include "FETypeEnum.hpp"
#include "GeomElTypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "TimeLoop.hpp"

#include "CurrentElem.hpp"

#include "OptLoop.hpp"

using namespace femus;

  void GenMatRhsMHDAD(MultiLevelProblem &ml_prob){

  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_MHDAD");
  const unsigned Level = my_system.GetLevelToAssemble();
  
  //========= parameters
   double IRem =  1./ml_prob.GetInputParser().get("Rem");
   double S    = ml_prob.GetInputParser().get("S");
  
//======== GEOMETRICAL ELEMENT =======
  const uint space_dim =       ml_prob._ml_msh->GetDimension();

  //=========== Operators 
  std::vector<double> VelAdj_vec_val_g(space_dim);
  std::vector<double> VelAdj_vec_val_g3D(3);
  double dphijdx_g[DIMENSION];
  double dphiidx_g[DIMENSION];
  double  curlBXlambda_g3D[3]; 
    
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
    CurrentQuantity BhomAdjOldX(currelem);
    BhomAdjOldX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj0"); 
    BhomAdjOldX.VectWithQtyFillBasic();
    BhomAdjOldX.Allocate();
    
    CurrentQuantity BhomAdjOldY(currelem);
    BhomAdjOldY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj1"); 
    BhomAdjOldY.VectWithQtyFillBasic();
    BhomAdjOldY.Allocate();
    
    CurrentQuantity BhomAdjOldZ(currelem);
    BhomAdjOldZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj2"); 
    BhomAdjOldZ.VectWithQtyFillBasic();
    BhomAdjOldZ.Allocate();
    
    std::vector<CurrentQuantity*> BhomAdjOld_vec;   
    BhomAdjOld_vec.push_back(&BhomAdjOldX);
    BhomAdjOld_vec.push_back(&BhomAdjOldY);
    BhomAdjOld_vec.push_back(&BhomAdjOldZ);
 
    CurrentQuantity BhomLagMultAdjOld(currelem);
    BhomLagMultAdjOld._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomLagMultAdj");
    BhomLagMultAdjOld.VectWithQtyFillBasic();
    BhomLagMultAdjOld.Allocate();    /// @todo probably this Allocate not needed here
//========= END INTERNAL QUANTITIES (unknowns of the equation) ================= 

//=========EXTERNAL QUANTITIES (couplings) =====

 //========= DOMAIN MAPPING
    CurrentQuantity xyz(currelem);
    xyz._dim      = DIMENSION;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();    

  //==========
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
 
  //==========
      CurrentQuantity VelAdjX(currelem);
    VelAdjX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_VelocityAdj0"); 
    VelAdjX.VectWithQtyFillBasic();
    VelAdjX.Allocate();
    
    CurrentQuantity VelAdjY(currelem);
    VelAdjY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_VelocityAdj1"); 
    VelAdjY.VectWithQtyFillBasic();
    VelAdjY.Allocate();
    
    CurrentQuantity VelAdjZ(currelem);
    VelAdjZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_VelocityAdj2"); 
    VelAdjZ.VectWithQtyFillBasic();
    VelAdjZ.Allocate();
    
    std::vector<CurrentQuantity*> VelAdj_vec;
    VelAdj_vec.push_back(&VelAdjX);
    VelAdj_vec.push_back(&VelAdjY);
    VelAdj_vec.push_back(&VelAdjZ);
 
    //==========    
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
 
//=========
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
 
//========= auxiliary, must be AFTER Bhom!
    CurrentQuantity Bmag(currelem); //total
    Bmag._dim        = Bhom_vec.size();
    Bmag._FEord      = BhomX._FEord;
    Bmag._ndof       = BhomX._ndof;
    Bmag.Allocate();    
    
//========= END EXTERNAL QUANTITIES =================

//=====================
//=====================    

    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.SetDofobjConnCoords();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);    

    currelem.SetElDofsBc();
    
    for (uint idim=0; idim < space_dim; idim++)    {
           BhomAdjOld_vec[idim]->GetElemDofs();
                  Vel_vec[idim]->GetElemDofs();
               VelAdj_vec[idim]->GetElemDofs();
                 Bhom_vec[idim]->GetElemDofs();
                 Bext_vec[idim]->GetElemDofs();
    }
    

//======SUM Bhom and Bext  //from now on, you'll only use Bmag //Bmag,Bext and Bhom must have the same orders!
    Math::zeroN(&Bmag._val_dofs[0],Bmag._dim*Bmag._ndof);

    for (uint ivarq=0; ivarq < Bmag._dim; ivarq++)    { //ivarq is like idim
          for (uint d=0; d < Bmag._ndof; d++)    {
          const uint     indxq  =         d + ivarq*Bmag._ndof;
          Bmag._val_dofs[indxq] = Bext_vec[ivarq]->_val_dofs[d] + Bhom_vec[ivarq]->_val_dofs[d];
	  }
    }    
//=======


    const uint el_ngauss = ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();
    
    for (uint qp = 0; qp < el_ngauss; qp++) {
//=======here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     {          
//   currgp.SetPhiElDofsFEVB_g (fe,qp); 
//   currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}  
	  
const double      det = 1.; //currgp.JacVectVV_g(xyz);
const double dtxJxW_g = det * ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { 
//   currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp);
//   currgp.ExtendDphiDxyzElDofsFEVB_g(fe);
}
//=======end of the "COMMON SHAPE PART"==================


      Bmag.curl_g();
      Bmag.val_g();
      
      for (uint idim=0; idim<space_dim; idim++)  {
          BhomAdjOld_vec[idim]->val_g(); 
                 Vel_vec[idim]->val_g(); 
              VelAdj_vec[idim]->val_g(); 
	VelAdj_vec_val_g[idim] = VelAdj_vec[idim]->_val_g[0];
      }

 //vector product
          Math::extend(&VelAdj_vec_val_g[0],&VelAdj_vec_val_g3D[0],space_dim);
          Math::cross(&Bmag._curl_g3D[0],&VelAdj_vec_val_g3D[0],curlBXlambda_g3D);

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//============================================================== 

// // //        for (uint i=0; i < BhomAdjOldX._ndof; i++)     {
// // // //======="COMMON tEST PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
// // //         const double phii_g = currgp._phi_ndsQLVB_g[BhomAdjOldX._FEord][i];
// // //         for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BhomAdjOldX._FEord][i+idim*BhomAdjOldX._ndof];
// // // //======= END "COMMON tEST PART for QTYZERO" ==========
// // // 
// // //    	  double BDdphii_g      = Math::dot(  &Bmag._val_g[0],dphiidx_g,space_dim);
// // // 	  double lambdaDdphii_g = Math::dot(&VelAdj_vec_val_g[0],dphiidx_g,space_dim);
// // // 	  
// // //          for (uint idim=0; idim<space_dim; idim++) {
// // //             const uint irowq = i+idim*BhomAdjOldX._ndof;
// // //            currelem.Rhs()(irowq) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                             - S*curlBXlambda_g3D[idim]*phii_g                             //from NS
// // //                             + S*(BDdphii_g*VelAdj_vec_val_g[idim] - lambdaDdphii_g*Bmag._val_g[idim])     //from NS
// // //                         )
// // //                            + (1-currelem.GetBCDofFlag()[irowq])*detb*BhomAdjOld_vec[idim]->_val_dofs[i]; //Dirichlet bc
// // // 	   }
// // // 
// // //         for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
// // //           const uint irowq = i+idim*BhomAdjOldX._ndof;
// // //           currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
// // //         }        // end filling diagonal for Dirichlet bc
// // // 	 
// // //         for (uint j=0; j<BhomAdjOldX._ndof; j++) {// A element matrix
// // // //======="COMMON SHAPE PART for QTYZERO": ==========
// // //            double                                  phij_g       =      currgp._phi_ndsQLVB_g[BhomAdjOldX._FEord][j];
// // //            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BhomAdjOldX._FEord][j+idim*BhomAdjOldX._ndof];
// // // //======= END "COMMON SHAPE PART for QTYZERO" ==========
// // // 
// // //            double Lap_g    = Math::dot(dphijdx_g,dphiidx_g,space_dim);
// // // 	  double Advphij_g = 0.;
// // // 	  for (uint idim=0; idim<space_dim; idim++) Advphij_g += Vel_vec[idim]->_val_g[0]*dphijdx_g[idim];
// // //           
// // //           for (uint idim=0; idim<space_dim; idim++) { //filled in as 1-2-3 // 4-5-6 // 7-8-9
// // //             int irowq = i+idim*BhomAdjOldX._ndof;
// // //             // diagonal blocks [1-5-9]
// // //             currelem.Mat()(irowq,j+idim*BhomAdjOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                  + LAP_MHD*IRem*(Lap_g)
// // //                  + (1-LAP_MHD)*IRem*(   Lap_g - dphijdx_g[idim]* dphiidx_g[idim] )
// // //                  - phii_g*(         Advphij_g - dphijdx_g[idim]*Vel_vec[idim]->_val_g[0] )
// // //                );
// // //             // block +1 [2-6-7]
// // //             int idimp1=(idim+1)%space_dim;
// // //             currelem.Mat()(irowq,j+idimp1*BhomAdjOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                  + (1-LAP_MHD)*IRem*(  -dphijdx_g[idim]* dphiidx_g[idimp1] )
// // //                  - phii_g*(            -dphijdx_g[idim]*Vel_vec[idimp1]->_val_g[0] )
// // //                );
// // // #if (DIMENSION==3)
// // //             // block +2 [3-4-8]
// // //             int idimp2=(idim+2)%space_dim;
// // //             currelem.Mat()(irowq,j+idimp2*BhomAdjOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                   + (1-LAP_MHD)*IRem*(-dphijdx_g[idim]* dphiidx_g[idimp2] )
// // //                   - phii_g*(          -dphijdx_g[idim]*Vel_vec[idimp2]->_val_g[0] )
// // //                  );
// // // #endif
// // //           }
// // //  
// // //         } 
// // //                                       // end A element matrix
// // //       
// // //        for (uint j=0; j < BhomLagMultAdjOld._ndof; j++) {// B^T element matrix ( p*div(v) )
// // //           const double psij_g =  currgp._phi_ndsQLVB_g[BhomLagMultAdjOld._FEord][j];
// // //           const int jclml= j + space_dim*BhomAdjOldX._ndof;
// // //           for (uint idim=0; idim<space_dim; idim++) {
// // //             uint irowq = i+idim*BhomAdjOldX._ndof;
// // //             currelem.Mat()(irowq,jclml) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(-psij_g*dphiidx_g[idim]);
// // //            }
// // //         }
// // //                                      // end B^T element matrix
// // // 
// // //           if (i<BhomLagMultAdjOld._ndof) {//  pressure equation (KOMP dp/dt=rho*div) 
// // //           double psii_g = currgp._phi_ndsQLVB_g[BhomLagMultAdjOld._FEord][i];
// // // 	  const uint irowl = i+space_dim*BhomAdjOldX._ndof;
// // //           currelem.Rhs()(irowl)=0.;  // rhs
// // //  //             currelem.Mat()(irowl,j+space_dim*el_ndof_q)  += dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt;
// // // 
// // //           for (uint j=0; j<BhomAdjOldX._ndof; j++) { // B element matrix q*div(u)
// // //             for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BhomAdjOldX._FEord][j+idim*BhomAdjOldX._ndof];
// // //             for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*BhomAdjOldX._ndof) += -dtxJxW_g*psii_g*dphijdx_g[idim]; 
// // //                 }
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
  
    // *****************************************************************
    // *****************************************************************
    // *****************************************************************

        my_system._LinSolver[Level]->_KK->close();
        my_system._LinSolver[Level]->_RESC->close();
  
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << my_system.name() << ": assembled  Level " << Level
           << " with " << my_system._LinSolver[Level]->_KK->m() << " dofs " << std::endl;
#endif     

  return;
}



