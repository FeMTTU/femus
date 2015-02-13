
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
#include "Quantity.hpp"
#include "QTYnumEnum.hpp"
#include "TimeLoop.hpp"
#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"


// application
#include "OptLoop.hpp"

using namespace femus;

  
  void GenMatRhsMHDCONT(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assemble_matrix) {

  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_MHDCONT");
  
   const double time =  0.; //ml_prob._timeloop._curr_time;
   
//======= TIME - STATIONARY OR NOT =======
  const int NonStatMHDCONT = (int) ml_prob.GetInputParser().get("NonStatMHDCONT");
  const double        dt   = 1.; //ml_prob._timeloop._timemap.get("dt");

//======== GEOMETRICAL ELEMENT =======
  const uint space_dim =       ml_prob.GetMeshTwo().get_dim();
  const uint mesh_ord  = (int) ml_prob.GetMeshTwo().GetRuntimeMap().get("mesh_ord");
  const uint meshql    = (int) ml_prob.GetMeshTwo().GetRuntimeMap().get("meshql"); //======== ELEMENT MAPPING =======

  //========= BCHandling =========
  const double penalty_val = ml_prob.GetMeshTwo().GetRuntimeMap().get("penalty_val");    
    //========= reference values =========
  
  //====== Physics
  double       IRem = 1./ml_prob.GetInputParser().get("Rem");
  double          S = ml_prob.GetInputParser().get("S");
  const double betaL2   = ml_prob.GetInputParser().get("betaL2");
  const double gammaLap = ml_prob.GetInputParser().get("gammaLap");
//===========================

  //==== Operators @ gauss ========   //TODO USER
  double dphijdx_g[DIMENSION];
  double dphiidx_g[DIMENSION];
  double      dphiidx_g3D[3];
  double curlxiXdphii_g3D[3];
  double   curlxiXvel_g3D[3];
  double curlbXlambda_g3D[3];
  double        Lapxi_g[DIMENSION];
//===========================

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
    CurrentQuantity BeOld(currgp);
    BeOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYZERO];
    BeOld.VectWithQtyFillBasic();
    BeOld.Allocate();

    CurrentQuantity LagMultOld(currgp);
    LagMultOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYONE];
    LagMultOld.VectWithQtyFillBasic();
    LagMultOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

    //==================Quadratic domain, auxiliary
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = NVE[ ml_prob.GetMeshTwo()._geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];
  xyz_refbox.Allocate();
  
#if VELOCITY_QTY==1
    CurrentQuantity Vel(currgp);
    Vel._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity"); //an alternative cannot exist, because it is an Unknown of This Equation
    Vel.VectWithQtyFillBasic();
    Vel.Allocate();
#endif  

    //==================
    CurrentQuantity VelAdj(currgp);
    VelAdj._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_VelocityAdj");
    VelAdj.VectWithQtyFillBasic();
    VelAdj.Allocate();
  
    //==================
    CurrentQuantity Bhom(currgp); 
    Bhom._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom");
    Bhom.VectWithQtyFillBasic();
    Bhom.Allocate();
    
//===============
    CurrentQuantity BhomAdj(currgp); 
    BhomAdj._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomAdj"); 
    BhomAdj.VectWithQtyFillBasic();
    BhomAdj.Allocate();
  //=========END EXTERNAL QUANTITIES (couplings) =====

 
    currelem.Mat().zero();
    currelem.Rhs().zero(); 
     
    currelem.SetDofobjConnCoords(ml_prob.GetMeshTwo()._iproc,iel);
    currelem.SetMidpoint();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);
    currelem.TransformElemNodesToRef(ml_prob.GetMeshTwo().GetDomain(),&xyz_refbox._val_dofs[0]);    

    currelem.SetElDofsBc();
    
         BeOld.GetElemDofs();  
    LagMultOld.GetElemDofs();

 
    
    if ( Vel._eqnptr != NULL )        Vel.GetElemDofs();
    else                              Vel._qtyptr->FunctionDof(Vel,time,&xyz_refbox._val_dofs[0]);
    if ( VelAdj._eqnptr != NULL )  VelAdj.GetElemDofs();
    else                           VelAdj._qtyptr->FunctionDof(VelAdj,time,&xyz_refbox._val_dofs[0]);
    if ( Bhom._eqnptr != NULL )      Bhom.GetElemDofs();
    else                             Bhom._qtyptr->FunctionDof(Bhom,time,&xyz_refbox._val_dofs[0]);
    if ( BhomAdj._eqnptr != NULL ) BhomAdj.GetElemDofs();
    else                           BhomAdj._qtyptr->FunctionDof(BhomAdj,time,&xyz_refbox._val_dofs[0]);    
    

//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
   const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();
    
    for (uint qp = 0; qp < el_ngauss; qp++) {

//======= here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     {          
  currgp.SetPhiElDofsFEVB_g (fe,qp);
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}

 const double      det = dt*currgp.JacVectVV_g(xyz);   //InvJac: is unique!
 const double dtxJxW_g = det*ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
 const double     detb = det/el_ngauss;

for (uint fe = 0; fe < QL; fe++)     {    
  currgp.SetDPhiDxyzElDofsFEVB_g (fe,qp);
  currgp.ExtendDphiDxyzElDofsFEVB_g (fe); 
}
//======= end of the "COMMON SHAPE PART"==================

//========preparation for things that are independent of (i,j), dofs of test and shape =====================
       BhomAdj.curl_g();
          Bhom.curl_g();
       BeOld.val_g(); 
         Vel.val_g();
      VelAdj.val_g();
        Bhom.val_g();
     BhomAdj.grad_g();

// vector product
      Math::extend(   &Vel._val_g[0],   &Vel._val_g3D[0],space_dim);
      Math::extend(&VelAdj._val_g[0],&VelAdj._val_g3D[0],space_dim);

      Math::cross(&BhomAdj._curl_g3D[0],   &Vel._val_g3D[0],curlxiXvel_g3D ); 
      Math::cross(   &Bhom._curl_g3D[0],&VelAdj._val_g3D[0],curlbXlambda_g3D );
//========end preparation for things that are independent of (i,j) dofs of test and shape =====================
  
//================================
//========= FILLING ELEMENT MAT/RHS
//=================================
//============ QTYZERO dofs (rows) ============
//=================================
//actually this is used for both QTYZERO and QTYONE,
//because the dofs of pressure are FEWER...
       for (uint i = 0; i < BeOld._ndof; i++)     {

//============ preparation for (i) ============
//(i) is used only for the RHS or for both MATRIX and RHS
//first, you put whatever is related to the Test Functions of THESE rows:
//test function value, test function first derivative
        const double                             phii_g       =      currgp._phi_ndsQLVB_g[BeOld._FEord][i];
        for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BeOld._FEord][i+idim*BeOld._ndof];

	 Math::extend(dphiidx_g,dphiidx_g3D,space_dim);
	 Math::cross(&BhomAdj._curl_g3D[0],dphiidx_g3D,curlxiXdphii_g3D);

	  double bDdphii_g      = Math::dot(  &Bhom._val_g[0],dphiidx_g,space_dim);
	  double lambdaDdphii_g = Math::dot(&VelAdj._val_g[0],dphiidx_g,space_dim);

	  for (uint idim=0; idim<space_dim; idim++) Lapxi_g[idim]=0.;
	  for (uint idim=0; idim<space_dim; idim++)  {
	      for (uint jdim=0; jdim<space_dim; jdim++) {
                    Lapxi_g[idim] += BhomAdj._grad_g[idim][jdim]*dphiidx_g[jdim];  //TODO CHECK THIS
                          }
	  }
//============end preparation for (i) ============

//============ QTYZERO rhs ============
         for (uint idim=0; idim<space_dim; idim++) {
            const uint irowq = i+idim*BeOld._ndof;
            
            currelem.Rhs()(irowq) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                          NonStatMHDCONT*BeOld._val_g[idim]*phii_g/dt
                           -(1-LAP_MHD)*IRem*curlxiXdphii_g3D[idim]                 //from MHD  1
                               -LAP_MHD*IRem*Lapxi_g[idim]                              //from MHD  2
                           + curlxiXvel_g3D[idim]*phii_g                            //from MHD
                           - S*curlbXlambda_g3D[idim]*phii_g                             //from NS
                           + S*(bDdphii_g*VelAdj._val_g[idim] - lambdaDdphii_g*Bhom._val_g[idim])     //from NS
                         )
                         +(1-currelem.GetBCDofFlag()[irowq])*detb*BeOld._val_dofs[irowq]; //Dirichlet bc
	   }
//============ END QTYZERO rhs ============

        for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
          const uint irowq = i+idim*BeOld._ndof;
          currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
        }
                                           // end filling diagonal for Dirichlet bc

//============ QTYZERO x QTYZERO dofs matrix (A matrix) ============
        for (uint j=0; j<BeOld._ndof; j++) {
//============preparation for (j) ============
//here you put things that depend either on (j) or on (i,j)
          double                                  phij_g       =      currgp._phi_ndsQLVB_g[BeOld._FEord][j];
          for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BeOld._FEord][j+idim*BeOld._ndof];

	  double Lap_g = Math::dot(dphijdx_g,dphiidx_g,space_dim);
	  double lambdaDdphii_g = Math::dot(&VelAdj._val_g[0],dphiidx_g,space_dim);  //TODO why did i recompute it?
	  double lambdaDdphij_g = Math::dot(&VelAdj._val_g[0],dphijdx_g,space_dim);
//============ end preparation for (j) ============

          for (uint idim=0; idim<space_dim; idim++) {
            int irowq = i+idim*BeOld._ndof;
            // diagonal blocks [1-5-9]
            currelem.Mat()(irowq,j+idim*BeOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                 NonStatMHDCONT*phij_g*phii_g/dt// time
                 + betaL2*phij_g*phii_g
                 + gammaLap*Lap_g
                 + S*phij_g*(lambdaDdphii_g  - VelAdj._val_g[idim]*dphiidx_g[idim] )
                 + S*phii_g*(lambdaDdphij_g  - VelAdj._val_g[idim]*dphijdx_g[idim] )
               );
            // block +1 [2-6-7]
            int idimp1=(idim+1)%space_dim;
            currelem.Mat()(irowq,j+idimp1*BeOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                 + S*phij_g*( -VelAdj._val_g[idim]*dphiidx_g[idimp1])
                 + S*phii_g*( -VelAdj._val_g[idimp1]*dphijdx_g[idim])
               );
#if (DIMENSION==3)
            // block +2 [3-4-8]
            int idimp2=(idim+2)%space_dim;
            currelem.Mat()(irowq,j+idimp2*BeOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                  + S*phij_g*( -VelAdj._val_g[idim]*dphiidx_g[idimp2]) 
                  + S*phii_g*( -VelAdj._val_g[idimp2]*dphijdx_g[idim])
               );
#endif
          }
                    
        } 
//============ END QTYZERO x QTYZERO dofs matrix (A matrix) ============

//============ QTYZERO x QTYONE dofs matrix (B^T matrix) // ( p*div(v) ) (NS eq) ============
       for (uint j=0; j < LagMultOld._ndof; j++) {
          const double psij_g = currgp._phi_ndsQLVB_g[LagMultOld._FEord][j];
          const int jclml = j+space_dim*BeOld._ndof;
          for (uint idim=0; idim<space_dim; idim++) {
            uint irowq=i+idim*BeOld._ndof;
            currelem.Mat()(irowq,jclml) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(-psij_g*dphiidx_g[idim]);
           }
        }
//============ END QTYZERO x QTYONE dofs matrix (B^T matrix) ============

//===================================================================
//============ END QTYZERO dofs (rows) ============
//===================================================================

//============================================
//============ QTYONE dofs (rows) ============
//===============================================
          if ( i < LagMultOld._ndof ) {
//============ preparation for (i) ============
          double psii_g = currgp._phi_ndsQLVB_g[LagMultOld._FEord][i];
	  const uint irowl = i+space_dim*BeOld._ndof;

//============ QTYONE rhs ============
          currelem.Rhs()(irowl)=0.;
//============ END QTYONE rhs ============

//============ QTYONE x QTYONE ============
 //             currelem.Mat()(irowl,j+space_dim*el_ndof_q)  += dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt; // (KOMP dp/dt=rho*div) 
//============ END QTYONE x QTYONE ============

//============ QTYONE x QTYZERO (B matrix) q*div(u) ============
          for (uint j = 0; j < BeOld._ndof; j++) { // B element matrix 
//============ preparation for (j) ============
            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim]= currgp._dphidxyz_ndsQLVB_g[BeOld._FEord][j+idim*BeOld._ndof];
            for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*BeOld._ndof) += -dtxJxW_g*psii_g*dphijdx_g[idim]; 
                }
//============ END QTYONE x QTYZERO ============

        }
//===============================================
//============ END QTYONE dofs (rows) ============
//================================================
      }
//===================================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) =====================
//===================================================================

    }
//==============================================================
//================== END GAUSS LOOP (qp loop) ======================
//==============================================================
    
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
    CurrentQuantity BeOld(currgp);
    BeOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYZERO];
    BeOld.VectWithQtyFillBasic();
    BeOld.Allocate();

    CurrentQuantity LagMultOld(currgp);
    LagMultOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYONE];
    LagMultOld.VectWithQtyFillBasic();
    LagMultOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING7
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

    //==================Quadratic domain, auxiliary
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = NVE[ ml_prob.GetMeshTwo()._geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];
  xyz_refbox.Allocate();
  
  //=========END EXTERNAL QUANTITIES (couplings) =====


     currelem.Mat().zero();
     currelem.Rhs().zero();

     currelem.SetDofobjConnCoords(ml_prob.GetMeshTwo()._iproc,iel);
     currelem.SetMidpoint();
     
     currelem.ConvertElemCoordsToMappingOrd(xyz);
     currelem.TransformElemNodesToRef(ml_prob.GetMeshTwo().GetDomain(),&xyz_refbox._val_dofs[0]);    

     currelem.SetElDofsBc();
     
          BeOld.GetElemDofs();
     LagMultOld.GetElemDofs();

//============ BC =======
       int press_fl = currelem.Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(BeOld,LagMultOld); 
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

        const double det   = dt*currgp.JacVectBB_g(xyz);
	const double dtxJxW_g = det*ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
//=======end "COMMON SHAPE PART"===================================

   xyz_refbox.val_g();
      LagMultOld._qtyptr->Function_txyz(time,&xyz_refbox._val_g[0]/*xyz._val_g*/,&LagMultOld._val_g[0]);  //i prefer using the function instead of the p_old vector

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
//============ QTYZERO dofs (rows) ============
      for (uint i=0; i<BeOld._ndof; i++) {
//============ preparation for (i) ============
const double phii_g = currgp._phi_ndsQLVB_g[BeOld._FEord][i];

//============ QTYZERO rhs ============
               for (uint idim=0; idim< space_dim; idim++)    {
             uint irowq=i+idim*BeOld._ndof;
            currelem.Rhs()(irowq)  += 
          currelem.GetBCDofFlag()[irowq]*           
           dtxJxW_g*(   -1.*press_fl*LagMultOld._val_g[0]*currgp.get_normal_ptr()[idim]*phii_g );   
 
	 }
           //end of idim loop
      
        }
//==============================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================   
      
    } 
//==================================================================
//================== END GAUSS LOOP (qp loop) ======================
//==================================================================
   
    my_system._A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    my_system._b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());

  }
  //end bdry element loop
      
  
    }  
    
// END BOUNDARY  // **************************

        my_system._A[Level]->close();
        my_system._b[Level]->close();
  
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << my_system.name() << ": assembled  Level " << Level
           << " with " << my_system._A[Level]->m() << " dofs " << std::endl;
#endif     

  return;
}


