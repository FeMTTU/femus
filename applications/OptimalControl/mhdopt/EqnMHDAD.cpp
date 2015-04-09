
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
#include "QTYnumEnum.hpp"
#include "TimeLoop.hpp"

#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"

#include "OptLoop.hpp"

using namespace femus;

  void GenMatRhsMHDAD(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assemble_matrix)  {

  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_MHDAD");
  
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
  
    CurrentElem       currelem(iel,myproc,Level,VV,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType(),mymsh);
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,ml_prob.GetQrule(currelem.GetDim()));
   
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity BhomAdjOldX(currgp);
    BhomAdjOldX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj0"); 
    BhomAdjOldX.VectWithQtyFillBasic();
    BhomAdjOldX.Allocate();
    
    CurrentQuantity BhomAdjOldY(currgp);
    BhomAdjOldY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj1"); 
    BhomAdjOldY.VectWithQtyFillBasic();
    BhomAdjOldY.Allocate();
    
    CurrentQuantity BhomAdjOldZ(currgp);
    BhomAdjOldZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj2"); 
    BhomAdjOldZ.VectWithQtyFillBasic();
    BhomAdjOldZ.Allocate();
    
    std::vector<CurrentQuantity*> BhomAdjOld_vec;   
    BhomAdjOld_vec.push_back(&BhomAdjOldX);
    BhomAdjOld_vec.push_back(&BhomAdjOldY);
    BhomAdjOld_vec.push_back(&BhomAdjOldZ);
 
    CurrentQuantity BhomLagMultAdjOld(currgp);
    BhomLagMultAdjOld._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomLagMultAdj");
    BhomLagMultAdjOld.VectWithQtyFillBasic();
    BhomLagMultAdjOld.Allocate();    
//========= END INTERNAL QUANTITIES (unknowns of the equation) ================= 

//=========EXTERNAL QUANTITIES (couplings) =====

 //========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();    

//========== Quadratic domain, auxiliary  
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = MESH_ORDER;
  xyz_refbox._ndof     = myel->GetElementDofNumber(ZERO_ELEM,BIQUADR_FE);
  xyz_refbox.Allocate();    

  //==========
      CurrentQuantity VelX(currgp);
    VelX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity0"); 
    VelX.VectWithQtyFillBasic();
    VelX.Allocate();
    
    CurrentQuantity VelY(currgp);
    VelY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity1"); 
    VelY.VectWithQtyFillBasic();
    VelY.Allocate();
    
    CurrentQuantity VelZ(currgp);
    VelZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity2"); 
    VelZ.VectWithQtyFillBasic();
    VelZ.Allocate();
    
    std::vector<CurrentQuantity*> Vel_vec;   
    Vel_vec.push_back(&VelX);
    Vel_vec.push_back(&VelY);
    Vel_vec.push_back(&VelZ);
 
  //==========
      CurrentQuantity VelAdjX(currgp);
    VelAdjX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_VelocityAdj0"); 
    VelAdjX.VectWithQtyFillBasic();
    VelAdjX.Allocate();
    
    CurrentQuantity VelAdjY(currgp);
    VelAdjY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_VelocityAdj1"); 
    VelAdjY.VectWithQtyFillBasic();
    VelAdjY.Allocate();
    
    CurrentQuantity VelAdjZ(currgp);
    VelAdjZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_VelocityAdj2"); 
    VelAdjZ.VectWithQtyFillBasic();
    VelAdjZ.Allocate();
    
    std::vector<CurrentQuantity*> VelAdj_vec;
    VelAdj_vec.push_back(&VelAdjX);
    VelAdj_vec.push_back(&VelAdjY);
    VelAdj_vec.push_back(&VelAdjZ);
 
    //==========    
    CurrentQuantity Bhom(currgp);
    Bhom._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom");
    Bhom.VectWithQtyFillBasic();
    Bhom.Allocate();    

//=========
    CurrentQuantity Bext(currgp);
    Bext._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt");
    Bext.VectWithQtyFillBasic();
    Bext.Allocate();    
  
//========= auxiliary, must be AFTER Bhom! //TODO this is an example of  Vect which is not associated to a Quantity
    CurrentQuantity Bmag(currgp); //total
    Bmag._dim        = Bhom._dim;               //same as Bhom
    Bmag._FEord      = Bhom._FEord;             //same as Bhom
    Bmag._ndof       = ml_prob.GetElemType()[currelem.GetDim()-1][Bmag._FEord]->GetNDofs();
    Bmag.Allocate();    
    
//========= END EXTERNAL QUANTITIES =================

//=====================
//=====================    

    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.SetDofobjConnCoords();
    currelem.SetMidpoint();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);    
    currelem.TransformElemNodesToRef(ml_prob._ml_msh->GetDomain(),&xyz_refbox._val_dofs[0]);    

    currelem.SetElDofsBc();
    
    BhomLagMultAdjOld.GetElemDofs();

    for (uint idim=0; idim < space_dim; idim++)    {
           BhomAdjOld_vec[idim]->GetElemDofs();
      
      if ( Vel_vec[idim]->_eqnptr != NULL )    Vel_vec[idim]->GetElemDofs();
    else                                       Vel_vec[idim]->_qtyptr->FunctionDof(*Vel_vec[idim],0.,&xyz_refbox._val_dofs[0]);
    if ( VelAdj_vec[idim]->_eqnptr != NULL) VelAdj_vec[idim]->GetElemDofs();
    else                                    VelAdj_vec[idim]->_qtyptr->FunctionDof(*VelAdj_vec[idim],0.,&xyz_refbox._val_dofs[0]);
    }
    
    if ( Bhom._eqnptr != NULL )     Bhom.GetElemDofs();
    else                            Bhom._qtyptr->FunctionDof(Bhom,0.,&xyz_refbox._val_dofs[0]);
    if ( Bext._eqnptr != NULL )     Bext.GetElemDofs();
    else                            Bext._qtyptr->FunctionDof(Bext,0.,&xyz_refbox._val_dofs[0]);

//======SUM Bhom and Bext  //from now on, you'll only use Bmag //Bmag,Bext and Bhom must have the same orders!
    Math::zeroN(&Bmag._val_dofs[0],Bmag._dim*Bmag._ndof);

    for (uint ivarq=0; ivarq < Bmag._dim; ivarq++)    { //ivarq is like idim
          for (uint d=0; d < Bmag._ndof; d++)    {
          const uint     indxq  =         d + ivarq*Bmag._ndof;
          Bmag._val_dofs[indxq] = Bext._val_dofs[indxq] + Bhom._val_dofs[indxq];
	  }
    }    
//=======


    const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();
    
    for (uint qp = 0; qp < el_ngauss; qp++) {
//=======here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     {          
  currgp.SetPhiElDofsFEVB_g (fe,qp); 
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}  
	  
const double      det = currgp.JacVectVV_g(xyz);   //InvJac: is the same for both QQ and LL!
const double dtxJxW_g = det * ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { 
  currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp);
  currgp.ExtendDphiDxyzElDofsFEVB_g(fe);
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

       for (uint i=0; i < BhomAdjOldX._ndof; i++)     {
//======="COMMON tEST PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
        const double phii_g = currgp._phi_ndsQLVB_g[BhomAdjOldX._FEord][i];
        for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BhomAdjOldX._FEord][i+idim*BhomAdjOldX._ndof];
//======= END "COMMON tEST PART for QTYZERO" ==========

   	  double BDdphii_g      = Math::dot(  &Bmag._val_g[0],dphiidx_g,space_dim);
	  double lambdaDdphii_g = Math::dot(&VelAdj_vec_val_g[0],dphiidx_g,space_dim);
	  
         for (uint idim=0; idim<space_dim; idim++) {
            const uint irowq = i+idim*BhomAdjOldX._ndof;
           currelem.Rhs()(irowq) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                            - S*curlBXlambda_g3D[idim]*phii_g                             //from NS
                            + S*(BDdphii_g*VelAdj_vec_val_g[idim] - lambdaDdphii_g*Bmag._val_g[idim])     //from NS
                        )
                           + (1-currelem.GetBCDofFlag()[irowq])*detb*BhomAdjOld_vec[idim]->_val_dofs[i]; //Dirichlet bc
	   }

        for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
          const uint irowq = i+idim*BhomAdjOldX._ndof;
          currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
        }        // end filling diagonal for Dirichlet bc
	 
        for (uint j=0; j<BhomAdjOldX._ndof; j++) {// A element matrix
//======="COMMON SHAPE PART for QTYZERO": ==========
           double                                  phij_g       =      currgp._phi_ndsQLVB_g[BhomAdjOldX._FEord][j];
           for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BhomAdjOldX._FEord][j+idim*BhomAdjOldX._ndof];
//======= END "COMMON SHAPE PART for QTYZERO" ==========

           double Lap_g    = Math::dot(dphijdx_g,dphiidx_g,space_dim);
	  double Advphij_g = 0.;
	  for (uint idim=0; idim<space_dim; idim++) Advphij_g += Vel_vec[idim]->_val_g[0]*dphijdx_g[idim];
          
          for (uint idim=0; idim<space_dim; idim++) { //filled in as 1-2-3 // 4-5-6 // 7-8-9
            int irowq = i+idim*BhomAdjOldX._ndof;
            // diagonal blocks [1-5-9]
            currelem.Mat()(irowq,j+idim*BhomAdjOldX._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                 + LAP_MHD*IRem*(Lap_g)
                 + (1-LAP_MHD)*IRem*(   Lap_g - dphijdx_g[idim]* dphiidx_g[idim] )
                 - phii_g*(         Advphij_g - dphijdx_g[idim]*Vel_vec[idim]->_val_g[0] )
               );
            // block +1 [2-6-7]
            int idimp1=(idim+1)%space_dim;
            currelem.Mat()(irowq,j+idimp1*BhomAdjOldX._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                 + (1-LAP_MHD)*IRem*(  -dphijdx_g[idim]* dphiidx_g[idimp1] )
                 - phii_g*(            -dphijdx_g[idim]*Vel_vec[idimp1]->_val_g[0] )
               );
#if (DIMENSION==3)
            // block +2 [3-4-8]
            int idimp2=(idim+2)%space_dim;
            currelem.Mat()(irowq,j+idimp2*BhomAdjOldX._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                  + (1-LAP_MHD)*IRem*(-dphijdx_g[idim]* dphiidx_g[idimp2] )
                  - phii_g*(          -dphijdx_g[idim]*Vel_vec[idimp2]->_val_g[0] )
                 );
#endif
          }
 
        } 
                                      // end A element matrix
      
       for (uint j=0; j < BhomLagMultAdjOld._ndof; j++) {// B^T element matrix ( p*div(v) )
          const double psij_g =  currgp._phi_ndsQLVB_g[BhomLagMultAdjOld._FEord][j];
          const int jclml= j + space_dim*BhomAdjOldX._ndof;
          for (uint idim=0; idim<space_dim; idim++) {
            uint irowq = i+idim*BhomAdjOldX._ndof;
            currelem.Mat()(irowq,jclml) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(-psij_g*dphiidx_g[idim]);
           }
        }
                                     // end B^T element matrix

          if (i<BhomLagMultAdjOld._ndof) {//  pressure equation (KOMP dp/dt=rho*div) 
          double psii_g = currgp._phi_ndsQLVB_g[BhomLagMultAdjOld._FEord][i];
	  const uint irowl = i+space_dim*BhomAdjOldX._ndof;
          currelem.Rhs()(irowl)=0.;  // rhs
 //             currelem.Mat()(irowl,j+space_dim*el_ndof_q)  += dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt;

          for (uint j=0; j<BhomAdjOldX._ndof; j++) { // B element matrix q*div(u)
            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BhomAdjOldX._FEord][j+idim*BhomAdjOldX._ndof];
            for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*BhomAdjOldX._ndof) += -dtxJxW_g*psii_g*dphijdx_g[idim]; 
                }
        }
                         // end pressure eq (cont)
      }
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

  {//BEGIN BOUNDARY  // *****************************************************************

    const uint mesh_vb = BB;
    
    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*ml_prob.GetMeshTwo()._iproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*ml_prob.GetMeshTwo()._iproc+Level];

 for (uint iel=0;iel < (nel_e - nel_b) ; iel++) {
  
    CurrentElem       currelem(iel,myproc,Level,BB,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType(),mymsh);
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,ml_prob.GetQrule(currelem.GetDim()));
   
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
      CurrentQuantity BhomAdjOldX(currgp);
    BhomAdjOldX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj0"); 
    BhomAdjOldX.VectWithQtyFillBasic();
    BhomAdjOldX.Allocate();
    
    CurrentQuantity BhomAdjOldY(currgp);
    BhomAdjOldY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj1"); 
    BhomAdjOldY.VectWithQtyFillBasic();
    BhomAdjOldY.Allocate();
    
    CurrentQuantity BhomAdjOldZ(currgp);
    BhomAdjOldZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj2"); 
    BhomAdjOldZ.VectWithQtyFillBasic();
    BhomAdjOldZ.Allocate();
    
    std::vector<CurrentQuantity*> BhomAdjOld_vec;   
    BhomAdjOld_vec.push_back(&BhomAdjOldX);
    BhomAdjOld_vec.push_back(&BhomAdjOldY);
    BhomAdjOld_vec.push_back(&BhomAdjOldZ);
 
    CurrentQuantity BhomLagMultAdjOld(currgp);
    BhomLagMultAdjOld._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomLagMultAdj");
    BhomLagMultAdjOld.VectWithQtyFillBasic();
    BhomLagMultAdjOld.Allocate();    
//========= END INTERNAL QUANTITIES (unknowns of the equation) ================= 

//=========EXTERNAL QUANTITIES (couplings) =====

 //========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();    

//========== Quadratic domain, auxiliary  
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = MESH_ORDER;
  xyz_refbox._ndof     = myel->GetElementFaceDofNumber(ZERO_ELEM,ZERO_FACE,BIQUADR_FE);
  xyz_refbox.Allocate();    
   
//========= END EXTERNAL QUANTITIES =================
  
//=====================
//=====================    


     currelem.Mat().zero();
     currelem.Rhs().zero();

     currelem.SetDofobjConnCoords();
     currelem.SetMidpoint();
     
     currelem.ConvertElemCoordsToMappingOrd(xyz);
     currelem.TransformElemNodesToRef(ml_prob._ml_msh->GetDomain(),&xyz_refbox._val_dofs[0]);    

     currelem.SetElDofsBc();
     
     BhomLagMultAdjOld.GetElemDofs();
   
//============ BC =======
       int press_fl = currelem.Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(BhomAdjOldX._ndof,space_dim,BhomLagMultAdjOld); 
//========END BC============
     
    const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();

    for (uint qp=0; qp< el_ngauss; qp++) {
//======= "COMMON SHAPE PART"============================
 for (uint fe = 0; fe < QL; fe++)     {
   currgp.SetPhiElDofsFEVB_g (fe,qp); 
   currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);
}

        const double det   = currgp.JacVectBB_g(xyz);
	const double dtxJxW_g = det * ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
//=======end "COMMON SHAPE PART"===================================   
      
      xyz_refbox.val_g();
      BhomLagMultAdjOld._qtyptr->Function_txyz(0.,&xyz_refbox._val_g[0]/*xyz._val_g*/,&BhomLagMultAdjOld._val_g[0]);  //i prefer using the function instead of the p_old vector
   
//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ===================
//==============================================================
           for (uint i=0; i < BhomAdjOldX._ndof; i++) {
 
	const double phii_g = currgp._phi_ndsQLVB_g[BhomAdjOldX._FEord][i];

        for (uint idim=0; idim< space_dim; idim++)    {
             uint irowq=i+idim*BhomAdjOldX._ndof;
            currelem.Rhs()(irowq)  += 
          currelem.GetBCDofFlag()[irowq]*           
           dtxJxW_g*(   -1.*press_fl*BhomLagMultAdjOld._val_g[0]*currgp.get_normal_ptr()[idim]*phii_g
        ); 

	   
 }
           //end of idim loop
	
	   } 
//==============================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) ===============
//==============================================================      
      
    }  //gauss
    
    my_system._LinSolver[Level]->_KK->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    my_system._LinSolver[Level]->_RESC->add_vector(currelem.Rhs(),currelem.GetDofIndices());

 }//elem loop

 
  }//END BOUNDARY ************************
  
        my_system._LinSolver[Level]->_KK->close();
        my_system._LinSolver[Level]->_RESC->close();
  
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << my_system.name() << ": assembled  Level " << Level
           << " with " << my_system._LinSolver[Level]->_KK->m() << " dofs " << std::endl;
#endif     

  return;
}



