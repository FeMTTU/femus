
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
#include "NormTangEnum.hpp"
#include "Quantity.hpp"
#include "QTYnumEnum.hpp"
#include "TimeLoop.hpp"
#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"

#include "FemusDefault.hpp"

#include "OptLoop.hpp"

using namespace femus;


 void GenMatRhsNSAD(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assemble_matrix) {

  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_NSAD");
  
    // ========= parameters
  const double alphaVel = ml_prob.GetInputParser().get("alphaVel");
  const double IRe      = 1./ml_prob.GetInputParser().get("Re");

  //=========== Operators 
  double dphijdx_g[DIMENSION];
  double dphiidx_g[DIMENSION];
    double    curlxiXB_g3D[3];

   const double time =  0.;  //ml_prob._timeloop._curr_time;
   
//======= TIME - STATIONARY OR NOT =======
const int NonStatNSAD = (int) ml_prob.GetInputParser().get("NonStatNSAD");
  const double   dt     = 1.;  //ml_prob._timeloop._timemap.get("dt");

// //======== GEOMETRICAL ELEMENT =======
  const uint space_dim = ml_prob.GetMeshTwo().get_dim();
  const uint  mesh_ord = (int) ml_prob.GetMeshTwo().GetRuntimeMap().get("mesh_ord");
  const uint    meshql = (int) ml_prob.GetMeshTwo().GetRuntimeMap().get("meshql");  //======== ELEMENT MAPPING =======

  
//========= BCHandling =========
  const double penalty_val = ml_prob.GetMeshTwo().GetRuntimeMap().get("penalty_val");    

        my_system._A[Level]->zero();
        my_system._b[Level]->zero();
  
   {//BEGIN VOLUME    
   
   const uint mesh_vb = VV;

    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*ml_prob.GetMeshTwo()._iproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*ml_prob.GetMeshTwo()._iproc+Level];
    
  for (uint iel=0; iel < (nel_e - nel_b); iel++) {
    
    CurrentElem       currelem(Level,VV,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType());
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,ml_prob.GetQrule(currelem.GetDim()));
   
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    //QTYZERO
    CurrentQuantity VelAdjOld(currgp);
    VelAdjOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYZERO];
    VelAdjOld.VectWithQtyFillBasic();
    VelAdjOld.Allocate();
  
    //QTYONE
    CurrentQuantity PressAdjOld(currgp);
    PressAdjOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYONE];
    PressAdjOld.VectWithQtyFillBasic();
    PressAdjOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================
  
//=========EXTERNAL QUANTITIES (couplings) =====

//========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

//========== Quadratic domain, auxiliary  
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = NVE[ ml_prob.GetMeshTwo()._geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];
  xyz_refbox.Allocate();
  
  CurrentQuantity Vel(currgp);
    Vel._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity");
    Vel.VectWithQtyFillBasic();
    Vel.Allocate();
  
  CurrentQuantity VelDes(currgp);
    VelDes._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_DesVelocity");
    VelDes.VectWithQtyFillBasic();
    VelDes.Allocate();

  CurrentQuantity Bhom(currgp);
    Bhom._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom");
    Bhom.VectWithQtyFillBasic();
    Bhom.Allocate();
 
  CurrentQuantity Bext(currgp);
    Bext._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt");
    Bext.VectWithQtyFillBasic();
    Bext.Allocate();

//========= auxiliary, must be AFTER Bhom!   //TODO this doesnt have any associated quantity!
  CurrentQuantity Bmag(currgp); //total
    Bmag._dim        = Bhom._dim;               //same as Bhom
    Bmag._FEord      = Bhom._FEord;             //same as Bhom
    Bmag._ndof       = ml_prob.GetElemType()[currelem.GetDim()-1][Bmag._FEord]->GetNDofs();
    Bmag.Allocate();
    
//===============
  CurrentQuantity BhomAdj(currgp); 
    BhomAdj._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj"); 
    BhomAdj.VectWithQtyFillBasic();
    BhomAdj.Allocate();
    
    //========= END EXTERNAL QUANTITIES =================


    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.SetDofobjConnCoords(ml_prob.GetMeshTwo()._iproc,iel);
    currelem.SetMidpoint();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);
    currelem.TransformElemNodesToRef(ml_prob.GetMeshTwo().GetDomain(),&xyz_refbox._val_dofs[0]);    

    currelem.SetElDofsBc();
    
    VelAdjOld.GetElemDofs();  
    PressAdjOld.GetElemDofs();

    if ( Vel._eqnptr != NULL )  Vel.GetElemDofs();
    else                         Vel._qtyptr->FunctionDof(Vel,time,&xyz_refbox._val_dofs[0]);    //give the Hartmann flow, if not solving NS
    if ( Bhom._eqnptr != NULL )  Bhom.GetElemDofs();
    else                         Bhom._qtyptr->FunctionDof(Bhom,time,&xyz_refbox._val_dofs[0]);
    if ( Bext._eqnptr != NULL )  Bext.GetElemDofs();
    else                         Bext._qtyptr->FunctionDof(Bext,time,&xyz_refbox._val_dofs[0]);
    if ( BhomAdj._eqnptr != NULL )  BhomAdj.GetElemDofs();
    else                            BhomAdj._qtyptr->FunctionDof(BhomAdj,time,&xyz_refbox._val_dofs[0]);    
    if ( VelDes._eqnptr != NULL )  VelDes.GetElemDofs();
    else                           VelDes._qtyptr->FunctionDof(VelDes,time,&xyz_refbox._val_dofs[0]);    

 
//======SUM Bhom and Bext  //from now on, you'll only use Bmag //Bmag,Bext and Bhom must have the same orders!
    Math::zeroN(&Bmag._val_dofs[0],Bmag._dim*Bmag._ndof);

    for (uint ivarq=0; ivarq < Bmag._dim; ivarq++)    { //ivarq is like idim
          for (uint d=0; d < Bmag._ndof; d++)    {
          const uint     indxq  =         d + ivarq*Bmag._ndof;
          Bmag._val_dofs[indxq] = Bext._val_dofs[indxq] + Bhom._val_dofs[indxq];
	  }
    }    
//=======    
    
//=======    
///optimal control
    xyz_refbox.SetElemAverage();
  int el_flagdom = ElFlagControl(xyz_refbox._el_average,&ml_prob.GetMeshTwo());
//=======    


   const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();

   for (uint qp = 0; qp < el_ngauss; qp++) {
//=======here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     { 
  currgp.SetPhiElDofsFEVB_g (fe,qp);  
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  }  
	  
const double      det = dt*currgp.JacVectVV_g(xyz);   //InvJac: is the same for both QQ and LL!
const double dtxJxW_g = det*ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { 
  currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp);
  currgp.ExtendDphiDxyzElDofsFEVB_g(fe); 
}
//=======end of the "COMMON SHAPE PART"==================

   VelAdjOld.val_g();
     BhomAdj.curl_g();
        Bmag.val_g();
         Vel.val_g();
      VelDes.val_g();
         Vel.grad_g();

//vector product
        Math::extend(&Bmag._val_g[0],&Bmag._val_g3D[0],space_dim);
        Math::cross(&BhomAdj._curl_g3D[0],&Bmag._val_g3D[0],curlxiXB_g3D);

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
       for (uint i = 0; i < VelAdjOld._ndof; i++)     {
//======="COMMON tEST PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
        const double phii_g = currgp._phi_ndsQLVB_g[VelAdjOld._FEord][i];
        for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim]= currgp._dphidxyz_ndsQLVB_g[VelAdjOld._FEord][i+idim*VelAdjOld._ndof];
//======= END "COMMON tEST PART for QTYZERO" ==========

         for (uint idim=0; idim<space_dim; idim++) {
            const uint irowq=i+idim*VelAdjOld._ndof;  //quadratic rows index
           currelem.Rhs()(irowq) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                            NonStatNSAD*VelAdjOld._val_g[idim]*phii_g/dt  //time
                           - curlxiXB_g3D[idim]*phii_g                    //this is due to the variation of velocity in the MAGNETIC ADVECTION, so it is due to a   NONLINEAR COUPLING "u times B", "MY_STATE times OTHER_STATE"
                           - alphaVel*el_flagdom*(Vel._val_g[idim] - VelDes._val_g[idim])*phii_g    //this is the dependence that counts
                           )
                           + (1-currelem.GetBCDofFlag()[irowq])*detb*VelAdjOld._val_dofs[irowq]; //Dirichlet bc
	   }

  for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
          const uint irowq = i+idim*VelAdjOld._ndof;
          currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
        }// end filling diagonal for Dirichlet bc
	 
        for (uint j=0; j<VelAdjOld._ndof; j++) {// A element matrix
//======="COMMON SHAPE PART for QTYZERO": ==========
           double                                  phij_g       =      currgp._phi_ndsQLVB_g[VelAdjOld._FEord][j];
           for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[VelAdjOld._FEord][j+idim*VelAdjOld._ndof];
//======= END "COMMON SHAPE PART for QTYZERO" ==========

          double     Lap_g = Math::dot(dphijdx_g,dphiidx_g,space_dim);
	  double Advphii_g = Math::dot(&Vel._val_g[0],dphiidx_g,space_dim);   //TODO can put it outside
          
          for (uint idim=0; idim<space_dim; idim++) { //filled in as 1-2-3 // 4-5-6 // 7-8-9
            int irowq = i+idim*VelAdjOld._ndof;
            // diagonal blocks [1-5-9]
            currelem.Mat()(irowq,j+idim*VelAdjOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                   NonStatNSAD* phij_g*phii_g/dt // time
                   + IRe*(dphijdx_g[idim]*dphiidx_g[idim] + Lap_g)      //Adjoint of D is D, Adjoint of Laplacian is Laplacian
                   + phij_g*phii_g*/*dveldx_g*/Vel._grad_g[idim][idim]  //Adjoint of Advection 1 delta(u) DOT grad(u): adj of nonlinear stuff has 2 TERMS (well, not always)
                   + phij_g*Advphii_g                                   //Adjoint of Advection 2 u DOT grad (delta(u)): adj of nonlinear stuff has 2 TERMS
               );
            // block +1 [2-6-7]
            int idimp1=(idim+1)%space_dim;
            currelem.Mat()(irowq,j+idimp1*VelAdjOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                   + IRe*(dphijdx_g[idim]*dphiidx_g[idimp1])
                   + phij_g*phii_g*/*dveldx_g*/Vel._grad_g[idimp1][idim]
               );
#if (DIMENSION==3)
            // block +2 [3-4-8]
            int idimp2=(idim+2)%space_dim;
            currelem.Mat()(irowq,j+idimp2*VelAdjOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                  + IRe*(dphijdx_g[idim]*dphiidx_g[idimp2])
                  + phij_g*phii_g*/*dveldx_g*/Vel._grad_g[idimp2][idim]
               );
#endif
          }
 
        } 
                                      // end A element matrix
      
       for (uint j=0; j<PressAdjOld._ndof; j++) {// B^T element matrix ( p*div(v) ) (NS eq)
          const double psij_g =  currgp._phi_ndsQLVB_g[PressAdjOld._FEord][j];
          const int jclml= j + space_dim*VelAdjOld._ndof;
          for (uint idim=0; idim<space_dim; idim++) {
            uint irowq = i+idim*VelAdjOld._ndof;
            currelem.Mat()(irowq,jclml) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(-psij_g*dphiidx_g[idim]);
           }
        }
                                     // end B^T element matrix

          if (i<PressAdjOld._ndof) {//  pressure equation (KOMP dp/dt=rho*div) 
          double psii_g = currgp._phi_ndsQLVB_g[PressAdjOld._FEord][i];
	  const uint irowl = i+space_dim*VelAdjOld._ndof;  //vertical offset
          currelem.Rhs()(irowl)=0.;  // rhs

          for (uint j=0; j<VelAdjOld._ndof; j++) { // B element matrix q*div(u)
            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[VelAdjOld._FEord][j+idim*VelAdjOld._ndof];
            for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*VelAdjOld._ndof) += -dtxJxW_g*psii_g*dphijdx_g[idim]; 
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
    my_system._A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    my_system._b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());
    
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
     
    CurrentElem       currelem(Level,BB,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType());

    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,ml_prob.GetQrule(currelem.GetDim()));
   
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    //QTYZERO
    CurrentQuantity VelAdjOld(currgp);
    VelAdjOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYZERO];
    VelAdjOld.VectWithQtyFillBasic();
    VelAdjOld.Allocate();
  
    //QTYONE
    CurrentQuantity PressAdjOld(currgp);
    PressAdjOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYONE];
    PressAdjOld.VectWithQtyFillBasic();
    PressAdjOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================
  
//=========EXTERNAL QUANTITIES (couplings) =====

//========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

//========== Quadratic domain, auxiliary  
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = NVE[ ml_prob.GetMeshTwo()._geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];
  xyz_refbox.Allocate();
  
//========= END EXTERNAL QUANTITIES =================


     currelem.Mat().zero();
     currelem.Rhs().zero();

     currelem.SetDofobjConnCoords(ml_prob.GetMeshTwo()._iproc,iel);
     currelem.SetMidpoint();
     
     currelem.ConvertElemCoordsToMappingOrd(xyz);
     currelem.TransformElemNodesToRef(ml_prob.GetMeshTwo().GetDomain(),&xyz_refbox._val_dofs[0]);    

     currelem.SetElDofsBc();

     VelAdjOld.GetElemDofs();
     PressAdjOld.GetElemDofs();

//============ BC =======
       int press_fl = currelem.Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(VelAdjOld,PressAdjOld); 
 //========END BC============
   
   //==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
   const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();
   
    for (uint qp=0; qp< el_ngauss; qp++) {
//======= "COMMON SHAPE PART"============================
 for (uint fe = 0; fe < QL; fe++)     { 
    currgp.SetPhiElDofsFEVB_g (fe,qp); 
    currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);
}

        const double det   = dt * currgp.JacVectBB_g(xyz);
	const double dtxJxW_g = det * ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
//=======end "COMMON SHAPE PART"===================================   
      
   xyz_refbox.val_g();
      PressAdjOld._qtyptr->Function_txyz(time,&xyz_refbox._val_g[0]/*xyz._val_g*/,&PressAdjOld._val_g[0]);  //i prefer using the function instead of the p_old vector
      

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
      for (uint i=0; i < VelAdjOld._ndof; i++) {

	const double phii_g = currgp._phi_ndsQLVB_g[VelAdjOld._FEord][i];

         for (uint idim=0; idim< space_dim; idim++)    {
             uint irowq=i+idim*VelAdjOld._ndof;
          currelem.Rhs()(irowq)  += 
            currelem.GetBCDofFlag()[irowq]*           
           dtxJxW_g*(   -1.*press_fl*PressAdjOld._val_g[0]*currgp.get_normal_ptr()[idim]*phii_g

	  );
	   
	 }
           //end of idim loop

     }
//==============================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================      
      
    }  //gauss
   
    my_system._A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    my_system._b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());

 }//elem loop
   
  }//END BOUNDARY ************************

        my_system._A[Level]->close();
        my_system._b[Level]->close();
  
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << my_system.name() << ": assembled  Level " << Level
           << " with " << my_system._A[Level]->m() << " dofs" << std::endl;
#endif    

return;
}

