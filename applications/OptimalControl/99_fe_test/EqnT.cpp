
#include "FemusDefault.hpp"

#include "NumericVector.hpp"
#include "DenseVector.hpp"
#include "SparseMatrix.hpp"
#include "DenseMatrix.hpp"

#include "Files.hpp"
#include "Math.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "SystemTwo.hpp"
#include "FETypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "VBTypeEnum.hpp"
#include "GeomElTypeEnum.hpp"
#include "Domain.hpp"
#include "TimeLoop.hpp"
#include "CurrentElem.hpp"
#include "paral.hpp"
#include "MultiLevelProblem.hpp"

// application
#include "TempQuantities.hpp"



 void GenMatRhsT(MultiLevelProblem &ml_prob) {
   
  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_T");
  const unsigned Level = my_system.GetLevelToAssemble();
  
  const double time =  0.;

  //======== ELEMENT MAPPING =======
  const uint space_dim =       ml_prob._ml_msh->GetDimension();
  
        my_system._LinSolver[Level]->_KK->zero();
        my_system._LinSolver[Level]->_RESC->zero();

// ==========================================  
  Mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(Level);
  elem		*myel		=  mymsh->el;
  const unsigned myproc  = mymsh->processor_id();
	
// ==========================================  
// ==========================================
  
  {//BEGIN VOLUME
  
  //==== AUXILIARY ==============
    double* dphijdx_g = new double[space_dim];
    double* dphiidx_g = new double[space_dim];

    
  const uint mesh_vb = VV;
  
   const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level+1];
   const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level];

  for (uint iel=0; iel < (nel_e - nel_b); iel++) {
    
  CurrentElem<double>       currelem(iel,myproc,Level,VV,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType(),mymsh);
  //   CurrentGaussPointBase & currgp = //   CurrentGaussPointBase::build(currelem,ml_prob.GetQuadratureRule(currelem.GetDim()));
  
//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Tempold(currelem);
    Tempold._qtyptr   = my_system.GetUnknownQuantitiesVector()[0]; 
    Tempold.VectWithQtyFillBasic();
    Tempold.Allocate();

    //=========EXTERNAL QUANTITIES (couplings) =====
    //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currelem);  //no quantity
    xyz._dim      = space_dim;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

   
    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.SetDofobjConnCoords();
    currelem.SetElDofsBc();

    currelem.ConvertElemCoordsToMappingOrd(xyz);

    Tempold.GetElemDofs();

//====================    
     const uint el_ngauss = ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();

    for (uint qp=0; qp< el_ngauss; qp++) {

//======= "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)   { 
//   currgp.SetPhiElDofsFEVB_g (fe,qp);
//   currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp); 
  }
	  
const double      det = 1.;//currgp.JacVectVV_g(xyz);
const double dtxJxW_g = det*ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     {
//   currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp); 
//   currgp.ExtendDphiDxyzElDofsFEVB_g(fe); 
  }
//======= end of the "COMMON SHAPE PART"==================

 	Tempold.val_g(); 
	   
	  
// // //       for (uint i=0; i < Tempold._ndof/*the maximum number is for biquadratic*/; i++)     {
// // // 
// // //         const double phii_g   = currgp._phi_ndsQLVB_g[Tempold._FEord][i];
// // // 
// // //         for (uint idim = 0; idim < space_dim; idim++) dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[Tempold._FEord][i+idim*Tempold._ndof];
// // // 
// // // //=========== FIRST ROW ===============
// // //         currelem.Rhs()(i) +=      
// // //            currelem.GetBCDofFlag()[i]*dtxJxW_g*( 
// // //                 7.*phii_g
// // // 	  )
// // // 	   + (1-currelem.GetBCDofFlag()[i])*detb*(Tempold._val_dofs[i]);
// // //         
// // //         currelem.Mat()(i,i) +=  (1-currelem.GetBCDofFlag()[i])*detb;
// // // 
// // // 	 // Matrix Assemblying ---------------------------
// // //         for (uint j=0; j<Tempold._ndof; j++) {
// // //           double phij_g   = currgp._phi_ndsQLVB_g[Tempold._FEord][j];
// // // 	  
// // //         for (uint idim = 0; idim < space_dim; idim++)   {
// // // 	  dphijdx_g  [idim] = currgp._dphidxyz_ndsQLVB_g[Tempold._FEord][j+idim*Tempold._ndof]; 
// // //           }
// // // 	  
// // // 	  
// // //           double Lap_g   = Math::dot(dphijdx_g,dphiidx_g,space_dim);
// // // 
// // // 	    int ip1 = i + Tempold._ndof;
// // // 	    int jp1 = j + Tempold._ndof;
// // // 
// // //  
// // // //============ FIRST ROW state  delta T ===============
// // // //======= DIAGONAL =============================
// // // 	   currelem.Mat()(i,j) +=        
// // //             currelem.GetBCDofFlag()[i]*dtxJxW_g*( 
// // //              Lap_g  
// // //             );
// // // 
// // //         }  //end j (col)
// // //       }   //end i (row)
      
      
    } // end of the quadrature point qp-loop

    currelem.Mat().print_scientific(std::cout);
    
       my_system._LinSolver[Level]->_KK->add_matrix(currelem.Mat(),currelem.GetDofIndices());
       my_system._LinSolver[Level]->_RESC->add_vector(currelem.Rhs(),currelem.GetDofIndices());
  } // end of element loop
  // *****************************************************************

  delete [] dphijdx_g;
  delete [] dphiidx_g;
  
  }//END VOLUME
  
        my_system._LinSolver[Level]->_KK->close();
        my_system._LinSolver[Level]->_RESC->close();

#ifdef DEFAULT_PRINT_INFO
  std::cout << " Matrix and RHS assembled for equation " << my_system.name()
            << " Level "<< Level << " dofs " << my_system._LinSolver[Level]->_KK->n() << std::endl;
#endif

  return;
}


