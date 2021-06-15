#ifndef __eqnt_hpp__
#define __eqnt_hpp__


#include "FemusDefault.hpp"

#include "NumericVector.hpp"
#include "DenseVector.hpp"
#include "SparseMatrix.hpp"
#include "DenseMatrix.hpp"
#include "LinearEquationSolver.hpp"

#include "Files.hpp"
#include "Math.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "MultiLevelProblem.hpp"
#include "SystemTwo.hpp"
#include "FETypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "VBTypeEnum.hpp"
#include "GeomElTypeEnum.hpp"
#include "Domain.hpp"
#include "TimeLoop.hpp"
#include "CurrentElem.hpp"
#include "OptLoop.hpp"
#include "paral.hpp"

// application
#include "TempQuantities.hpp"


// The question is: WHERE is the ORDER of the VARIABLES established?
// I mean, on one hand there is an ORDER REGARDLESS of the FE family;
//on the other hand there is an ORDER IN TERMS of FE.
// The ORDER is important for defining the BLOCKS in the VECTORS and MATRICES.
//In the initMeshToDof, i order QQ, LL and KK variables, BUT STILL YOU DO NOT KNOW
// to WHICH QUANTITIES those VARIABLES are ASSOCIATED!
// We should make a map that, for each SCALAR VARIABLE of EACH FE TYPE, 
// tells you TO WHICH QUANTITY it is ASSOCIATED.
// The point is then you have SCALAR variables associated to SCALAR QUANTITIES,
// or SCALAR VARIABLES associated to VECTOR QUANTITIES.
// So you should associate SCALAR VARIABLES to COMPONENTS of QUANTITIES.
// The first time when you start associating scalar variables to quantities
// is when you set the BC's


// ================================================================================
// Now, the idea in the construction of the ELEMENT MATRIX and ELEMENT RHS is this.
// You fill the element matrix, and each row (and column) of it will correspond
// to a certain row of the global matrix.
// HOW DO I KNOW THAT THIS CORRESPONDENCE is OK?
// Well, this is given by the el_dof_indices.
// every component of el_dof_indices tells me what is the GLOBAL ROW (and GLOBAL COLUMN)
// So, if the length of your el_dof_indices is 10, then you'll have 10x10 elements to add to the global matrix,
// given by all the possible couples "el_dof_indices[i] X el_dof_indices[j]"
// TODO: Now, how do you know that the phi of one line is the correct phi, so it is the phi of THAT DOF?
// I mean, how is the dependence on the REAL domain (NOT canonical) involved?
// For the phi, the values are the same both on the canonical and on the stretched domain.
// For the derivatives, I think the point is: you must pick the REAL dphidx in the SAME ORDER as you pick the CORRESPONDING DOFS.
// Now, my point is: on a given row, are you sure that the code picks the correct dphidx?


/// This function assembles the matrix and the rhs:
void  GenMatRhsT(MultiLevelProblem &ml_prob){

  //if we are just a function not inside a class, we have to retrieve ourselves...
  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_T");
  const unsigned Level = my_system.GetLevelToAssemble();
// ==========================================  
  Mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(Level);
  elem		*myel		=  mymsh->el;
  const unsigned myproc  = mymsh->processor_id();
	
  
  //======== ELEMENT MAPPING =======
  const uint space_dim =       ml_prob._ml_msh->GetDimension();

  //====== reference values ========================
  const double IRe = 1./ml_prob.GetInputParser().get("Re");
  const double IPr = 1./ml_prob.GetInputParser().get("Pr");
  const double alphaT  = ml_prob.GetInputParser().get("alphaT");
  const double alphaL2 = ml_prob.GetInputParser().get("alphaL2");
  const double alphaH1 = ml_prob.GetInputParser().get("alphaH1");
  
  //==== AUXILIARY ==============
    std::vector<double> dphijdx_g(space_dim);
    std::vector<double> dphiidx_g(space_dim);

        my_system._LinSolver[Level]->_KK->zero();
        my_system._LinSolver[Level]->_RESC->zero();

// ==========================================  
// ==========================================  
 {//BEGIN VOLUME
 
   const uint mesh_vb = VV;
   
   const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level];
   const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level+1];

//    const uint nel_beg = mymsh->_elementOffset[myproc];
//    const uint nel_end = mymsh->_elementOffset[myproc+1];
   
  for (uint iel = 0; iel < (nel_e - nel_b); iel++) {

//   for (uint iel_two = nel_beg; iel_two < nel_end; iel_two++) {
  
  CurrentElem<double>       currelem(iel,myproc,Level,VV,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType(),mymsh);    
  //   CurrentGaussPointBase & currgp = //   CurrentGaussPointBase::build(currelem,ml_prob.GetQuadratureRule(currelem.GetDim()));
  

//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Tempold(currelem);
    Tempold._SolName = "Qty_Temperature";
    Tempold._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_Temperature"); 
    Tempold.VectWithQtyFillBasic();
    Tempold.Allocate();

//====================================
    CurrentQuantity Tlift(currelem);
    Tlift._SolName = "Qty_TempLift";
    Tlift._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_TempLift");
    Tlift.VectWithQtyFillBasic();
    Tlift.Allocate();

//=====================================
    CurrentQuantity TAdj(currelem);
    TAdj._SolName = "Qty_TempAdj";
    TAdj._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_TempAdj"); 
    TAdj.VectWithQtyFillBasic();
    TAdj.Allocate();
   
//=========EXTERNAL QUANTITIES (couplings) =====
    //========= //DOMAIN MAPPING
  CurrentQuantity xyz(currelem);  //no quantity
    xyz._dim      = space_dim;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

  //==================
    CurrentQuantity velX(currelem);
    velX._SolName = "Qty_Velocity0";
    velX._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity0"); 
    velX.VectWithQtyFillBasic();
    velX.Allocate();
    
  //==================
    CurrentQuantity velY(currelem);
    velY._SolName = "Qty_Velocity1";
    velY._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity1"); 
    velY.VectWithQtyFillBasic();
    velY.Allocate();    
    
//===============Tdes=====================
    CurrentQuantity Tdes(currelem);
    Tdes._SolName = "Qty_TempDes";
    Tdes._dim      = Tempold._dim;
    Tdes._FEord    = Tempold._FEord;
    Tdes._ndof     = Tempold._ndof;
    Tdes.Allocate();

  
// ==========================================  
// ==========================================  

    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.SetDofobjConnCoords();
    currelem.SetElDofsBc();
    

  Tempold.GetElemDofs();
    Tlift.GetElemDofs();
     TAdj.GetElemDofs();
     
     velX.GetElemDofs();
     velY.GetElemDofs();
   
     TempDesired(Tdes,currelem);
     
    currelem.ConvertElemCoordsToMappingOrd(xyz);
    xyz.SetElemAverage();
    int domain_flag = ElFlagControl(xyz._el_average,ml_prob._ml_msh);
//====================    

   const uint el_ngauss = ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();
   
   for (uint qp=0; qp< el_ngauss; qp++) {

//======= "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)   { 
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
//======= end of the "COMMON SHAPE PART"==================

 	Tempold.val_g(); 
          Tlift.val_g(); 
           TAdj.val_g(); 
           velX.val_g(); 
           velY.val_g(); 
           Tdes.val_g();
 
// // //       for (uint i=0; i < Tempold._ndof; i++)     {
// // // 
// // //         const double phii_g = currgp._phi_ndsQLVB_g[Tempold._FEord][i];
// // // 
// // //         for (uint idim = 0; idim < space_dim; idim++) dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[Tempold._FEord][i+idim*Tempold._ndof];
// // // 
// // // //=========== FIRST ROW ===============
// // //         currelem.Rhs()(i) +=      
// // //            currelem.GetBCDofFlag()[i]*dtxJxW_g*( 0. )
// // // 	   + (1-currelem.GetBCDofFlag()[i])*detb*(Tempold._val_dofs[i]);
// // //         
// // //         currelem.Mat()(i,i) +=  (1-currelem.GetBCDofFlag()[i])*detb;
// // // 
// // // //========= SECOND ROW (CONTROL) =====================
// // // 	 int ip1 = i + /* 1* */Tempold._ndof;   //suppose that T' T_0 T_adj have the same order
// // // 	 currelem.Rhs()(ip1) +=      
// // //            currelem.GetBCDofFlag()[ip1]*dtxJxW_g*( 
// // //                      + alphaT*domain_flag*(Tdes._val_g[0])*phii_g // T_d delta T_0
// // // 	  )
// // // 	   + (1-currelem.GetBCDofFlag()[ip1])*detb*(Tlift._val_dofs[i]);
// // //         
// // //          currelem.Mat()(ip1,ip1) +=  (1-currelem.GetBCDofFlag()[ip1])*detb;
// // // 
// // // //======= THIRD ROW (ADJOINT) ===================================
// // // 	 int ip2 = i + 2 * Tempold._ndof;   //suppose that T' T_0 T_adj have the same order
// // //            currelem.Rhs()(ip2) +=      
// // //            currelem.GetBCDofFlag()[ip2]*dtxJxW_g*( 
// // //                 + alphaT*domain_flag*(Tdes._val_g[0])*phii_g // T_d delta T'
// // // 	     )
// // // 	   + (1-currelem.GetBCDofFlag()[ip2])*detb*(Tempold._val_dofs[i]);
// // //         
// // //         currelem.Mat()(ip2,ip2) +=  (1-currelem.GetBCDofFlag()[ip2])*detb;
// // // 
// // // 	 // Matrix Assemblying ---------------------------
// // //         for (uint j=0; j<Tempold._ndof; j++) {
// // //           double phij_g = currgp._phi_ndsQLVB_g[Tempold._FEord][j];
// // // 	  
// // //         for (uint idim = 0; idim < space_dim; idim++)  dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[Tempold._FEord][j+idim*Tempold._ndof]; 
// // //            
// // //    
// // //           double Lap_g   = Math::dot(&dphijdx_g[0],&dphiidx_g[0],space_dim);
// // //           double Advection = velX._val_g[0]*dphijdx_g[0] + velY._val_g[0]*dphijdx_g[1]; //Math::dot(&vel._val_g[0],&dphijdx_g[0],space_dim);
// // // 
// // // 	    int ip1 = i + Tempold._ndof;
// // // 	    int jp1 = j + Tempold._ndof;
// // // 	    int ip2 = i + 2*Tempold._ndof;
// // // 	    int jp2 = j + 2*Tempold._ndof;
// // // 
// // // // 	           T     T_0     T_adj
// // // 	    
// // // // 	    T      X      X       O
// // // 	     
// // // // 	    T_0   
// // // 	    
// // // // 	    T_adj
// // // 	    
// // // 	    
// // // //============ FIRST ROW state  delta T ===============
// // // //======= DIAGONAL =============================
// // // 	   currelem.Mat()(i,j) +=        
// // //             currelem.GetBCDofFlag()[i]*dtxJxW_g*( 
// // //             + Advection*phii_g
// // //             + IRe*IPr*Lap_g  
// // //             );
// // // 
// // // //===============================
// // //     //same operators for T and T_0
// // // 	    currelem.Mat()(i,jp1) +=        
// // //             currelem.GetBCDofFlag()[i]*dtxJxW_g*(    
// // //             + Advection*phii_g
// // //             + IRe*IPr*Lap_g
// // // 	    );
// // // 
// // // //====================================
// // // 	   currelem.Mat()(i,jp2) +=        
// // //             currelem.GetBCDofFlag()[i]*dtxJxW_g*( 
// // //                 0.
// // //             );
// // // 
// // // 	    
// // // //============= SECOND ROW (LIFTING) delta T_0 =============
// // // //===== DIAGONAL ===========================
// // //          currelem.Mat()(ip1,jp1) +=        
// // //             currelem.GetBCDofFlag()[ip1]*
// // //             dtxJxW_g*( 
// // //              + alphaL2*phij_g*phii_g  //L_2 control norm
// // //              + alphaH1*Lap_g          //H_1 control norm
// // //               + alphaT*domain_flag*(phij_g)*phii_g  //T_0 delta T_0  //ADDED///////////////
// // //             ); 
// // // //====================================
// // // 	   currelem.Mat()(ip1,j) +=        
// // //             currelem.GetBCDofFlag()[ip1]*
// // //             dtxJxW_g*( 
// // //                 + alphaT*domain_flag*(phij_g)*phii_g  //T' delta T_0     //ADDED///////////////
// // //             );
// // // //====================================
// // // 	   currelem.Mat()(ip1,jp2) +=        
// // //             currelem.GetBCDofFlag()[ip1]*
// // //              dtxJxW_g*( 
// // //                  -Advection*phii_g
// // //                 + IRe*IPr*Lap_g
// // //            );
// // // 
// // // //============= THIRD ROW (ADJOINT) =============
// // // //======= DIAGONAL ==================
// // //           currelem.Mat()(ip2,jp2) +=        
// // //             currelem.GetBCDofFlag()[ip2]*
// // //               dtxJxW_g*( 
// // //             - Advection*phii_g  //minus sign
// // //             + IRe*IPr*Lap_g
// // //                
// // //             ); 
// // // //====================================
// // // 	   currelem.Mat()(ip2,j) +=        
// // //             currelem.GetBCDofFlag()[ip2]*
// // //             dtxJxW_g*( 
// // //                + alphaT*domain_flag*(phij_g)*phii_g  //T' delta T'
// // //             );
// // // //====================================
// // // 	   currelem.Mat()(ip2,jp1) +=        
// // //             currelem.GetBCDofFlag()[ip2]*
// // //             dtxJxW_g*( 
// // //                + alphaT*domain_flag*(phij_g)*phii_g  //T_0 delta T'     ///ADDED///////
// // //             );
// // // 
// // //         }  //end j (col)
// // //       }   //end i (row)
    } // end of the quadrature point qp-loop

       my_system._LinSolver[Level]->_KK->add_matrix(currelem.Mat(),currelem.GetDofIndices());
       my_system._LinSolver[Level]->_RESC->add_vector(currelem.Rhs(),currelem.GetDofIndices());
  } // end of element loop
  // *****************************************************************

 
   }//END VOLUME
  
        my_system._LinSolver[Level]->_KK->close();
        my_system._LinSolver[Level]->_RESC->close();

 
#ifdef DEFAULT_PRINT_INFO
  std::cout << " Matrix and RHS assembled for equation " << my_system.name()
            << " Level "<< Level << " dofs " << my_system._LinSolver[Level]->_KK->n() << std::endl;
#endif

  return;
}




#endif
