
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
#include "VBTypeEnum.hpp"
#include "GeomElTypeEnum.hpp"
#include "Domain.hpp"
#include "CurrentElem.hpp"
#include "Box.hpp"
  

//application
#include "OptLoop.hpp"
#include "TempQuantities.hpp"


//===================================================
/// This function assembles the matrix and the rhs:
 void GenMatRhsNS(MultiLevelProblem &ml_prob)  {

   SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_NS");
   const unsigned Level = my_system.GetLevelToAssemble();
    
   const uint   _AdvPic_fl = 1;
   const uint   _AdvNew_fl = 0;
   const uint   _Stab_fl = 0;
   const double _Komp_fac = 0.;

  //========== GEOMETRIC ELEMENT ========
  const uint           space_dim = ml_prob._ml_msh->GetDimension();

  //====== reference values ========================
//====== related to Quantities on which Operators act, and to the choice of the "LEADING" EQUATION Operator
  const double IRe = 1./ml_prob.GetInputParser().get("Re");
  const double IFr = 1./ml_prob.GetInputParser().get("Fr");
//================================================  

//================================================  
//=======Operators @ gauss =======================
  std::vector<double>    dphijdx_g(space_dim); // ShapeDer(): used for Laplacian,Divergence, ..., associated to an Unknown Quantity
  std::vector<double>    dphiidx_g(space_dim); // Test(): used for Laplacian,Advection,Divergence... associated to an Unknown Quantity
  std::vector<double>     AdvRhs_g(space_dim); //Operator: Adv(u,u,phi)
//================================================  

        my_system._LinSolver[Level]->_KK->zero();
        my_system._LinSolver[Level]->_RESC->zero();

// ==========================================  
  Mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(Level);
  elem		*myel		=  mymsh->el;
  const unsigned myproc  = mymsh->processor_id();
	
// ==========================================  
// ==========================================  

  {//BEGIN VOLUME
//========================
//========================
    const uint mesh_vb = VV;
  
    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level];

  for (int iel=0; iel < (nel_e - nel_b); iel++) {
    
    CurrentElem<double>       currelem(iel,myproc,Level,VV,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType(),mymsh);
    //   CurrentGaussPointBase & currgp = //   CurrentGaussPointBase::build(currelem,ml_prob.GetQuadratureRule(currelem.GetDim()));
 
  
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity VelOldX(currelem);
    VelOldX._qtyptr  = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity0");
    VelOldX._SolName = "Qty_Velocity0";
    VelOldX.VectWithQtyFillBasic();
    VelOldX.Allocate();

    CurrentQuantity VelOldY(currelem);
    VelOldY._qtyptr  = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity1");
    VelOldY._SolName = "Qty_Velocity1";
    VelOldY.VectWithQtyFillBasic();
    VelOldY.Allocate();

    std::vector<CurrentQuantity*> VelOld_vec;   
    VelOld_vec.push_back(&VelOldX);
    VelOld_vec.push_back(&VelOldY);
    
  
    const uint   qtyzero_ord  = VelOldX._FEord;
    const uint   qtyzero_ndof = VelOldX._ndof;  //same as Y 

//=========
    CurrentQuantity pressOld(currelem);
    pressOld._qtyptr  = ml_prob.GetQtyMap().GetQuantity("Qty_Pressure");
    pressOld._SolName = "Qty_Pressure";
    pressOld.VectWithQtyFillBasic();
    pressOld.Allocate();

   const uint qtyone_ord  = pressOld._FEord;
   const uint qtyone_ndof = pressOld._ndof; 

   //order
   const uint  qtyZeroToOne_DofOffset = VelOldX._ndof + VelOldY._ndof;//VelOld._ndof*VelOld._dim;
   
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currelem); //domain
    xyz._dim      = space_dim;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();
    
//other Physical constant Quantities
//=======gravity==================================
  CurrentQuantity gravity(currelem);
  gravity._dim = space_dim;
//   gravity.Allocate(); CANNOT DO THIS NOW BECAUSE NOT ALL THE DATA FOR THE ALLOCATION ARE FILLED
  gravity._val_g.resize(gravity._dim);
  gravity._val_g[0] = ml_prob.GetInputParser().get("dirgx");
  gravity._val_g[1] = ml_prob.GetInputParser().get("dirgy");
  if ( space_dim == 3 )   gravity._val_g[2] = ml_prob.GetInputParser().get("dirgz"); 

//========================
//========================
 
    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.SetDofobjConnCoords();
    currelem.SetElDofsBc();

    currelem.ConvertElemCoordsToMappingOrd(xyz);

//=======RETRIEVE the DOFS of the UNKNOWN QUANTITIES,i.e. MY EQUATION
     VelOldX.GetElemDofs();
     VelOldY.GetElemDofs();
    pressOld.GetElemDofs();

   const uint el_ngauss = ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();
   
    for (uint qp = 0; qp < el_ngauss; qp++) {  

//=======here starts the "COMMON SHAPE PART"==================
// the call of these things should be related to the Operator
//also, the choice of the QuadratureRule should be dependent of the Involved Operators
//and the FE orders on which these Operators act

//these  phi and dphi are used for the different stages:
//BEFORE i: by the interpolation functions
//INSIDE i: for the tEST functions and derivatives
//INSIDE j: for the SHAPE functions and derivatives
//again, it should be the Operator to decide what functions to be called,
//for what FE ORDER and what DERIVATIVE ORDER
//if we decide that the PREPARATION of the tEST and SHAPE 
//of a certain Unknown are COMMON TO ALL,
//Then we must only concentrate on preparing the OTHER involved quantities in that Operator
for (uint fe = 0; fe < QL; fe++)     {          
//   currgp.SetPhiElDofsFEVB_g (fe,qp);
//   currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
  }
	  
const double      det = 1.; //currgp.JacVectVV_g(xyz);
const double dtxJxW_g = det*ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { 
//   currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp);
//   currgp.ExtendDphiDxyzElDofsFEVB_g(fe);
  }
//=======end of the "COMMON SHAPE PART"==================

//now we want to fill the element matrix and rhs with ONLY values AT GAUSS POINTS
//we can divide the values at Gauss points into three parts:
//1-constant values
//2-values which depend on (i)
//3-values which depend on (i,j)
//1- can be used for filling both Ke and Fe
//2- can be used to fill Fe and also Ke
//3- can be used only for filling Ke

//Here, before entering the (i,j) loop, you compute quantities that DO NOT depend on (i,j),
//therefore the ELEMENT SUM, GAUSS SUM And tEST/SHAPE sums have ALL been performed
//but, these quantities depend on idim and jdim, because they are  involved in multiplications with tEST and SHAPE functions

//Internal Quantities
      VelOldX.val_g();
      VelOldX.grad_g();
      VelOldY.val_g();
      VelOldY.grad_g();
      
//Advection all VelOld
      for (uint idim=0; idim<space_dim; idim++) { AdvRhs_g[idim]=0.;}
      
       for (uint idim=0; idim<space_dim; idim++)  {
          for (uint b=0; b<space_dim; b++) {
	    AdvRhs_g[idim] += VelOld_vec[b]->_val_g[0]*VelOld_vec[idim]->_grad_g[0][b];
	  } 
      }   // grad is [ivar][idim], i.e. [u v w][x y z]	  
	    
//Divergence VelOld
	  double Div_g=0.;
          for (uint idim=0; idim<space_dim; idim++)    Div_g += VelOld_vec[idim]->_grad_g[0][idim];

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
// TODO according to the order we should switch DIM loop and DOF loop
//     for (uint i=0; i<qtyzero_ndof; i++)     {
// //======="COMMON tEST PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
//                                     const double phii_g       =      currgp._phi_ndsQLVB_g[qtyzero_ord][i];
//         for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[qtyzero_ord][i+idim*qtyzero_ndof];
// //======= END "COMMON tEST PART for QTYZERO" ==========
// 	
// 	  for (uint idim=0; idim<space_dim; idim++) {
//             const uint irowq=i+idim*qtyzero_ndof;  //  (i):       dof of the tEST function
//                                                   //(idim): component of the tEST function
//            currelem.Rhs()(irowq) += 
//          currelem.GetBCDofFlag()[irowq]*
//            dtxJxW_g*(       + _AdvNew_fl*          AdvRhs_g[idim]*phii_g     // NONLIN
//                             +            IFr*gravity._val_g[idim]*phii_g     // gravity                           
//                                )
//             + (1-currelem.GetBCDofFlag()[irowq])*detb*VelOld_vec[idim]->_val_dofs[i] //Dirichlet bc    
// 	;
//           }
//            // end filling element rhs u
// 
//   for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
//           const uint irowq = i+idim*qtyzero_ndof;
//           currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
//         }
//                                        // end filling diagonal for Dirichlet bc
// 
// //============ QTYZERO x QTYZERO dofs matrix (A matrix) ============
//         for (uint j=0; j< qtyzero_ndof; j++) {
// //======="COMMON SHAPE PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
// //   (j):       dof of the SHAPE function
//            double                                  phij_g       =      currgp._phi_ndsQLVB_g[qtyzero_ord][j];
//            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[qtyzero_ord][j+idim*qtyzero_ndof];
// //======= END "COMMON SHAPE PART for QTYZERO" ==========
//   
//           double Lap_g = Math::dot(&dphijdx_g[0],&dphiidx_g[0],space_dim);
// 	  double Adv_g=0.;
// 	  for (uint idim=0; idim<space_dim; idim++) Adv_g += VelOld_vec[idim]->_val_g[0] * dphijdx_g[idim]; // =Math::dot(&VelOld._val_g[0],dphijdx_g,space_dim);
//          
//           for (uint idim=0; idim<space_dim; idim++) { //filled in as 1-2-3 // 4-5-6 // 7-8-9
//             int irowq = i+idim*qtyzero_ndof;      //(i) is still the dof of the tEST functions
//                                                   //(idim): component of the tEST functions
// 
//             currelem.Mat()(irowq,j+idim*qtyzero_ndof)  // diagonal blocks [1-5-9] [idim(rows),idim(columns)]  //(idim): component of the SHAPE functions
//                += 
//             currelem.GetBCDofFlag()[irowq]*    
//             dtxJxW_g*(
//                  + _AdvPic_fl*                           Adv_g*phii_g                //TODO NONLIN
//                  + _AdvNew_fl*phij_g*VelOld_vec[idim]->_grad_g[0][idim]*phii_g        //TODO NONLIN
//                  + _AdvPic_fl*_Stab_fl*        0.5*Div_g*phij_g*phii_g                //TODO NONLIN
//                  +                         IRe*(      dphijdx_g[idim]*dphiidx_g[idim] + Lap_g)
//                );
// 
//             int idimp1=(idim+1)%space_dim;    // block +1 [2-6-7] [idim(rows),idim+1(columns)]  //(idimp1): component of the SHAPE functions
//             currelem.Mat()(irowq,j+idimp1*qtyzero_ndof)
//                +=
//             currelem.GetBCDofFlag()[irowq]*
//             dtxJxW_g*(
//                    _AdvNew_fl*phij_g*VelOld_vec[idim]->_grad_g[0][idimp1]*phii_g           //TODO NONLIN
//                               +            IRe*(     dphijdx_g[idim]*dphiidx_g[idimp1])
//                );
// 
//           }
//  
//         } 
//============ END QTYZERO x QTYZERO dofs matrix (A matrix) ============

//============ QTYZERO x QTYONE dofs matrix (B^T matrix) // ( p*div(v) ) (NS eq) ============
//        for (uint j=0; j<qtyone_ndof; j++) {
// //======="COMMON SHAPE PART for QTYONE" ==================
// 	 const double psij_g = currgp._phi_ndsQLVB_g[qtyone_ord][j];
// //======="COMMON SHAPE PART for QTYONE" - END ============
// 
//           const int jclml = j + qtyZeroToOne_DofOffset;
//           for (uint idim=0; idim<space_dim; idim++) {
//             uint irowq = i+idim*qtyzero_ndof;
//             currelem.Mat()(irowq,jclml) +=
//                currelem.GetBCDofFlag()[irowq]*                    
//                dtxJxW_g*(-psij_g*dphiidx_g[idim]);   /**   (-1.)*/
// 	    
//            } 
// 
//         }
// //============ END QTYZERO x QTYONE dofs matrix (B^T matrix) ============
// 
//           if (i < qtyone_ndof) {
// //======="COMMON tEST PART for QTYONE" ============
//           double psii_g = currgp._phi_ndsQLVB_g[qtyone_ord][i];
// //======= "COMMON tEST PART for QTYONE" - END ============
// 	  const uint irowl = i + qtyZeroToOne_DofOffset;
//           currelem.Rhs()(irowl)=0.;  // rhs
//  //             Mat()(irowl,j+space_dim*qtyzero_ndof)  += (1./dt)*dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt;  //no bc here (KOMP dp/dt=rho*div)
// 
//           for (uint j=0; j<qtyzero_ndof; j++) { // B element matrix q*div(u)
// //======="COMMON SHAPE PART for QTYZERO" ==================
//             for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[qtyzero_ord][j+idim*qtyzero_ndof];
// //======="COMMON SHAPE PART for QTYZERO" - END ============
// 	    
//             for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*qtyzero_ndof) += -/*(1./dt)**/dtxJxW_g*psii_g*dphijdx_g[idim]; 
//                 }
// 
//         }
//                          // end pressure eq (cont)
//       }
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
           << " with " << my_system._LinSolver[Level]->_KK->m() << " dofs" << std::endl;
#endif

    return;
}
