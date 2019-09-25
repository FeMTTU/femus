
//library includes
#include "FemusDefault.hpp"

#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"

#include "Math.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "MultiLevelProblem.hpp"
#include "FETypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "GeomElTypeEnum.hpp"
#include "Quantity.hpp"
#include "TimeLoop.hpp"
#include "CurrentElem.hpp"


// application
#include "OptLoop.hpp"

using namespace femus;

  
  void GenMatRhsMHDCONT(MultiLevelProblem &ml_prob){

  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_MHDCONT");
  const unsigned Level = my_system.GetLevelToAssemble();
  
//======== GEOMETRICAL ELEMENT =======
  const uint space_dim =       ml_prob._ml_msh->GetDimension();

  //====== Physics
  double       IRem = 1./ml_prob.GetInputParser().get("Rem");
  double          S = ml_prob.GetInputParser().get("S");
  const double betaL2   = ml_prob.GetInputParser().get("betaL2");
  const double gammaLap = ml_prob.GetInputParser().get("gammaLap");
//===========================

  //==== Operators @ gauss ========
  std::vector<double> Vel_vec_val_g(space_dim);
  std::vector<double> Vel_vec_val_g3D(3);
  std::vector<double> VelAdj_vec_val_g(space_dim);
  std::vector<double> VelAdj_vec_val_g3D(3);
  double dphijdx_g[DIMENSION];
  double dphiidx_g[DIMENSION];
  double      dphiidx_g3D[3];
  double curlxiXdphii_g3D[3];
  double   curlxiXvel_g3D[3];
  double curlbXlambda_g3D[3];
  double        Lapxi_g[DIMENSION];
//===========================

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
    CurrentQuantity BextOldX(currelem);
    BextOldX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt0"); 
    BextOldX.VectWithQtyFillBasic();
    BextOldX.Allocate();
    
    CurrentQuantity BextOldY(currelem);
    BextOldY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt1"); 
    BextOldY.VectWithQtyFillBasic();
    BextOldY.Allocate();
    
    CurrentQuantity BextOldZ(currelem);
    BextOldZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt2"); 
    BextOldZ.VectWithQtyFillBasic();
    BextOldZ.Allocate();
    
    std::vector<CurrentQuantity*> BextOld_vec;   
    BextOld_vec.push_back(&BextOldX);
    BextOld_vec.push_back(&BextOldY);
    BextOld_vec.push_back(&BextOldZ);

    CurrentQuantity LagMultOld(currelem);
    LagMultOld._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExtLagMult");
    LagMultOld.VectWithQtyFillBasic();
    LagMultOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currelem);
    xyz._dim      = DIMENSION;
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
 
    CurrentQuantity Bhom_vecQuant(currelem);   //without quantity, nor equation
    Bhom_vecQuant._dim     = Bhom_vec.size();
    Bhom_vecQuant._FEord   = Bhom_vec[0]->_FEord;
    Bhom_vecQuant._ndof    = Bhom_vec[0]->_ndof;
    Bhom_vecQuant.Allocate();
   
//===============
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
  //=========END EXTERNAL QUANTITIES (couplings) =====

 
    currelem.Mat().zero();
    currelem.Rhs().zero(); 
     
    currelem.SetDofobjConnCoords();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);

    currelem.SetElDofsBc();
    
 
    for (uint idim=0; idim < space_dim; idim++)    {
          BextOld_vec[idim]->GetElemDofs(); 
              Vel_vec[idim]->GetElemDofs();
           VelAdj_vec[idim]->GetElemDofs();
          BhomAdj_vec[idim]->GetElemDofs();
             Bhom_vec[idim]->GetElemDofs();
    }
    
    
    BhomAdj_vecQuant.GetElemDofs(BhomAdj_vec);
    Bhom_vecQuant.GetElemDofs(Bhom_vec);

//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
   const uint el_ngauss = ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();
    
    for (uint qp = 0; qp < el_ngauss; qp++) {

//======= here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     {          
//   currgp.SetPhiElDofsFEVB_g (fe,qp);
//   currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}

 const double      det = 1.;//currgp.JacVectVV_g(xyz);
 const double dtxJxW_g = det*ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussWeight(qp);
 const double     detb = det/el_ngauss;

for (uint fe = 0; fe < QL; fe++)     {    
//   currgp.SetDPhiDxyzElDofsFEVB_g (fe,qp);
//   currgp.ExtendDphiDxyzElDofsFEVB_g (fe); 
}
//======= end of the "COMMON SHAPE PART"==================

//========preparation for things that are independent of (i,j), dofs of test and shape =====================
       BhomAdj_vecQuant.curl_g();
       Bhom_vecQuant.curl_g();

   for (uint idim=0; idim<space_dim; idim++)   {
      BhomAdj_vec[idim]->grad_g();
         Bhom_vec[idim]->val_g();
          Vel_vec[idim]->val_g();
       VelAdj_vec[idim]->val_g();
          Vel_vec_val_g[idim] =  Vel_vec[idim]->_val_g[0];
       VelAdj_vec_val_g[idim] =  VelAdj_vec[idim]->_val_g[0];
   }
	
// vector product
      Math::extend(   &Vel_vec_val_g[0],   &Vel_vec_val_g3D[0],space_dim);
      Math::extend(&VelAdj_vec_val_g[0],&VelAdj_vec_val_g3D[0],space_dim);

      Math::cross(&BhomAdj_vecQuant._curl_g3D[0],   &Vel_vec_val_g3D[0],curlxiXvel_g3D ); 
      Math::cross(   &Bhom_vecQuant._curl_g3D[0],&VelAdj_vec_val_g3D[0],curlbXlambda_g3D );
//========end preparation for things that are independent of (i,j) dofs of test and shape =====================
  
//================================
//========= FILLING ELEMENT MAT/RHS
//=================================
//============ QTYZERO dofs (rows) ============
//=================================
// // // //actually this is used for both QTYZERO and QTYONE,
// // // //because the dofs of pressure are FEWER...
// // //        for (uint i = 0; i < BextOldX._ndof; i++)     {
// // // 
// // // //============ preparation for (i) ============
// // // //(i) is used only for the RHS or for both MATRIX and RHS
// // // //first, you put whatever is related to the Test Functions of THESE rows:
// // // //test function value, test function first derivative
// // //         const double                             phii_g       =      currgp._phi_ndsQLVB_g[BextOldX._FEord][i];
// // //         for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BextOldX._FEord][i+idim*BextOldX._ndof];
// // // 
// // // 	 Math::extend(dphiidx_g,dphiidx_g3D,space_dim);
// // // 	 Math::cross(&BhomAdj_vecQuant._curl_g3D[0],dphiidx_g3D,curlxiXdphii_g3D);
// // // 
// // // 	  double bDdphii_g = 0.;
// // // 	  for (uint idim=0; idim<space_dim; idim++)  bDdphii_g += Bhom_vec[idim]->_val_g[0]*dphiidx_g[idim];
// // // 	  
// // // 	  double lambdaDdphii_g = Math::dot(&VelAdj_vec_val_g[0],dphiidx_g,space_dim);
// // // 
// // // 	  for (uint idim=0; idim<space_dim; idim++) Lapxi_g[idim]=0.;
// // // 	  for (uint idim=0; idim<space_dim; idim++)  {
// // // 	      for (uint jdim=0; jdim<space_dim; jdim++) {
// // //                     Lapxi_g[idim] += BhomAdj_vec[idim]->_grad_g[0][jdim]*dphiidx_g[jdim];  //TODO CHECK THIS
// // //                           }
// // // 	  }
// // // //============end preparation for (i) ============
// // // 
// // // //============ QTYZERO rhs ============
// // //          for (uint idim=0; idim<space_dim; idim++) {
// // //             const uint irowq = i+idim*BextOldX._ndof;
// // //             
// // //             currelem.Rhs()(irowq) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                            -(1-LAP_MHD)*IRem*curlxiXdphii_g3D[idim]                 //from MHD  1
// // //                                -LAP_MHD*IRem*Lapxi_g[idim]                              //from MHD  2
// // //                            + curlxiXvel_g3D[idim]*phii_g                            //from MHD
// // //                            - S*curlbXlambda_g3D[idim]*phii_g                             //from NS
// // //                            + S*(bDdphii_g*VelAdj_vec[idim]->_val_g[0] - lambdaDdphii_g*Bhom_vec[idim]->_val_g[0])     //from NS
// // //                          )
// // //                          +(1-currelem.GetBCDofFlag()[irowq])*detb*BextOld_vec[idim]->_val_dofs[i]; //Dirichlet bc
// // // 	   }
// // // //============ END QTYZERO rhs ============
// // // 
// // //         for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
// // //           const uint irowq = i+idim*BextOldX._ndof;
// // //           currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
// // //         }
// // //                                            // end filling diagonal for Dirichlet bc
// // // 
// // // //============ QTYZERO x QTYZERO dofs matrix (A matrix) ============
// // //         for (uint j=0; j<BextOldX._ndof; j++) {
// // // //============preparation for (j) ============
// // // //here you put things that depend either on (j) or on (i,j)
// // //           double                                  phij_g       =      currgp._phi_ndsQLVB_g[BextOldX._FEord][j];
// // //           for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BextOldX._FEord][j+idim*BextOldX._ndof];
// // // 
// // // 	  double Lap_g = Math::dot(dphijdx_g,dphiidx_g,space_dim);
// // // 	  double lambdaDdphii_g = Math::dot(&VelAdj_vec_val_g[0],dphiidx_g,space_dim);  //TODO why did i recompute it?
// // // 	  double lambdaDdphij_g = Math::dot(&VelAdj_vec_val_g[0],dphijdx_g,space_dim);
// // // //============ end preparation for (j) ============
// // // 
// // //           for (uint idim=0; idim<space_dim; idim++) {
// // //             int irowq = i+idim*BextOldX._ndof;
// // //             // diagonal blocks [1-5-9]
// // //             currelem.Mat()(irowq,j+idim*BextOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                  + betaL2*phij_g*phii_g
// // //                  + gammaLap*Lap_g
// // //                  + S*phij_g*(lambdaDdphii_g  - VelAdj_vec[idim]->_val_g[0]*dphiidx_g[idim] )
// // //                  + S*phii_g*(lambdaDdphij_g  - VelAdj_vec[idim]->_val_g[0]*dphijdx_g[idim] )
// // //                );
// // //             // block +1 [2-6-7]
// // //             int idimp1=(idim+1)%space_dim;
// // //             currelem.Mat()(irowq,j+idimp1*BextOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                  + S*phij_g*( -VelAdj_vec[idim]->_val_g[0]*dphiidx_g[idimp1])
// // //                  + S*phii_g*( -VelAdj_vec[idimp1]->_val_g[0]*dphijdx_g[idim])
// // //                );
// // // #if (DIMENSION==3)
// // //             // block +2 [3-4-8]
// // //             int idimp2=(idim+2)%space_dim;
// // //             currelem.Mat()(irowq,j+idimp2*BextOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                   + S*phij_g*( -VelAdj_vec[idim]->_val_g[0]*dphiidx_g[idimp2]) 
// // //                   + S*phii_g*( -VelAdj_vec[idimp2]->_val_g[0]*dphijdx_g[idim])
// // //                );
// // // #endif
// // //           }
// // //                     
// // //         } 
// // // //============ END QTYZERO x QTYZERO dofs matrix (A matrix) ============
// // // 
// // // //============ QTYZERO x QTYONE dofs matrix (B^T matrix) // ( p*div(v) ) (NS eq) ============
// // //        for (uint j=0; j < LagMultOld._ndof; j++) {
// // //           const double psij_g = currgp._phi_ndsQLVB_g[LagMultOld._FEord][j];
// // //           const int jclml = j+space_dim*BextOldX._ndof;
// // //           for (uint idim=0; idim<space_dim; idim++) {
// // //             uint irowq=i+idim*BextOldX._ndof;
// // //             currelem.Mat()(irowq,jclml) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(-psij_g*dphiidx_g[idim]);
// // //            }
// // //         }
// // // //============ END QTYZERO x QTYONE dofs matrix (B^T matrix) ============
// // // 
// // // //===================================================================
// // // //============ END QTYZERO dofs (rows) ============
// // // //===================================================================
// // // 
// // // //============================================
// // // //============ QTYONE dofs (rows) ============
// // // //===============================================
// // //           if ( i < LagMultOld._ndof ) {
// // // //============ preparation for (i) ============
// // //           double psii_g = currgp._phi_ndsQLVB_g[LagMultOld._FEord][i];
// // // 	  const uint irowl = i+space_dim*BextOldX._ndof;
// // // 
// // // //============ QTYONE rhs ============
// // //           currelem.Rhs()(irowl)=0.;
// // // //============ END QTYONE rhs ============
// // // 
// // // //============ QTYONE x QTYONE ============
// // //  //             currelem.Mat()(irowl,j+space_dim*el_ndof_q)  += dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt; // (KOMP dp/dt=rho*div) 
// // // //============ END QTYONE x QTYONE ============
// // // 
// // // //============ QTYONE x QTYZERO (B matrix) q*div(u) ============
// // //           for (uint j = 0; j < BextOldX._ndof; j++) { // B element matrix 
// // // //============ preparation for (j) ============
// // //             for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim]= currgp._dphidxyz_ndsQLVB_g[BextOldX._FEord][j+idim*BextOldX._ndof];
// // //             for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*BextOldX._ndof) += -dtxJxW_g*psii_g*dphijdx_g[idim]; 
// // //                 }
// // // //============ END QTYONE x QTYZERO ============
// // // 
// // //         }
// // // //===============================================
// // // //============ END QTYONE dofs (rows) ============
// // // //================================================
// // //       }
//===================================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) =====================
//===================================================================

    }
//==============================================================
//================== END GAUSS LOOP (qp loop) ======================
//==============================================================
    
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
           << " with " << my_system._LinSolver[Level]->_KK->m() << " dofs " << std::endl;
#endif     

  return;
}


