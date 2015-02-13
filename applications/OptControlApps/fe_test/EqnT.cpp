
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
#include "QTYnumEnum.hpp"
#include "Domain.hpp"
#include "TimeLoop.hpp"
#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"
#include "paral.hpp"
#include "MultiLevelProblem.hpp"

// application
#include "TempQuantities.hpp"



 void GenMatRhsT(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assemble_matrix) {
   
  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_T");

  const double time =  0.;

//========== PROCESSOR INDEX
  const uint myproc = ml_prob.GetMeshTwo()._iproc;

//========= BCHandling =========
  const double penalty_val = ml_prob.GetMeshTwo().GetRuntimeMap().get("penalty_val");    

  //======== ELEMENT MAPPING =======
  const uint space_dim =       ml_prob.GetMeshTwo().get_dim();
  const uint  meshql   = (int) ml_prob.GetMeshTwo().GetRuntimeMap().get("meshql");
  const uint  mesh_ord = (int) ml_prob.GetMeshTwo().GetRuntimeMap().get("mesh_ord");
  
        my_system._A[Level]->zero();
        my_system._b[Level]->zero();

  {//BEGIN VOLUME
  
  //==== AUXILIARY ==============
    double* dphijdx_g = new double[space_dim];
    double* dphiidx_g = new double[space_dim];
    double* dphijdx_gLL = new double[space_dim];
    double* dphiidx_gLL = new double[space_dim];
    double* dphijdx_gKK = new double[space_dim];
    double* dphiidx_gKK = new double[space_dim];

    
  const uint mesh_vb = VV;
  
   const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level+1];
   const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level];

  for (uint iel=0; iel < (nel_e - nel_b); iel++) {
    
  CurrentElem       currelem(Level,VV,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType());
  
  CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,ml_prob.GetQrule(currelem.GetDim()));
  
//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Tempold(currgp);
    Tempold._qtyptr   = my_system.GetUnknownQuantitiesVector()[0]; 
    Tempold.VectWithQtyFillBasic();
    Tempold.Allocate();

//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Temp2(currgp);
    Temp2._qtyptr   = my_system.GetUnknownQuantitiesVector()[1]; 
    Temp2.VectWithQtyFillBasic();
    Temp2.Allocate();

//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Temp3(currgp);
    Temp3._qtyptr   = my_system.GetUnknownQuantitiesVector()[2]; 
    Temp3.VectWithQtyFillBasic();
    Temp3.Allocate();
    
    //=========EXTERNAL QUANTITIES (couplings) =====
    //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);  //no quantity
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

    //==================Quadratic domain, auxiliary, must be QUADRATIC!!! ==========
  CurrentQuantity xyz_refbox(currgp);  //no quantity
    xyz_refbox._dim      = space_dim;
    xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
    xyz_refbox._ndof     = NVE[ ml_prob.GetMeshTwo()._geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];
    xyz_refbox.Allocate();

    
    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.SetDofobjConnCoords(myproc,iel);
    currelem.SetMidpoint();

    currelem.ConvertElemCoordsToMappingOrd(xyz);
    currelem.TransformElemNodesToRef(ml_prob.GetMeshTwo().GetDomain(),&xyz_refbox._val_dofs[0]);    

    
//MY EQUATION
//the elements are, for every level:
// 1)DOF INDICES
// 2)BC FLAGS
// 3)BC VALUES 
// 1) and 2) are taken in a single vector, 3) are considered separately
      
    currelem.SetElDofsBc();

  Tempold.GetElemDofs();
    Temp2.GetElemDofs();
    Temp3.GetElemDofs();
    
// ===============      
// Now the point is this: there are several functions of space
// which are expressed with respect to a reference frame
//ok, now that the dofs are filled for xyz_refbox, I can use the el_average
//Well, the alternative is to consider  the elem in the refbox as
    //either a Vect or a CurrentElem !
    //I could consider it as another element, but only with the geometrical part!

  xyz_refbox.SetElemAverage();
//====================    
    
//===== FILL the DOFS of the EXTERNAL QUANTITIES: you must assure that for every Vect the quantity is set correctly
   const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();

    for (uint qp=0; qp< el_ngauss; qp++) {

//======= "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)   { 
  currgp.SetPhiElDofsFEVB_g (fe,qp);
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp); 
  }
	  
const double      det = currgp.JacVectVV_g(xyz);
const double dtxJxW_g = det*ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     {
  currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp); 
  currgp.ExtendDphiDxyzElDofsFEVB_g(fe); 
  }
//======= end of the "COMMON SHAPE PART"==================

 	Tempold.val_g(); 
 	  Temp2.val_g(); 
 	  Temp3.val_g(); 
	   
	   // always remember to get the dofs for the variables you use!
           // The point is that you fill the dofs with different functions...
           // you should need a flag to check if the dofs have been correctly filled
 
	  
      for (uint i=0; i < Tempold._ndof/*the maximum number is for biquadratic*/; i++)     {

        const double phii_g   = currgp._phi_ndsQLVB_g[Tempold._FEord][i];
        const double phii_gLL = currgp._phi_ndsQLVB_g[Temp2._FEord][i];
        const double phii_gKK = currgp._phi_ndsQLVB_g[Temp3._FEord][i];

        for (uint idim = 0; idim < space_dim; idim++) dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[Tempold._FEord][i+idim*Tempold._ndof];

//=========== FIRST ROW ===============
        currelem.Rhs()(i) +=      
           currelem.GetBCDofFlag()[i]*dtxJxW_g*( 
                7.*phii_g
	  )
	   + (1-currelem.GetBCDofFlag()[i])*detb*(Tempold._val_dofs[i]);
        
        currelem.Mat()(i,i) +=  (1-currelem.GetBCDofFlag()[i])*detb;

//========= SECOND ROW =====================
	 int ip1 = i + Tempold._ndof; 
	 
	if (i < currelem.GetElemType(Temp2._FEord)->GetNDofs() ) { 
	 currelem.Rhs()(ip1) +=      
           currelem.GetBCDofFlag()[ip1]*dtxJxW_g*( 
                0.07*phii_gLL
	  )
	   + (1-currelem.GetBCDofFlag()[ip1])*detb*(Temp2._val_dofs[i]);
        
         currelem.Mat()(ip1,ip1) +=  (1-currelem.GetBCDofFlag()[ip1])*detb;
	}
	
//======= THIRD ROW ===================================
	 int ip2 = i + Tempold._ndof + Temp2._ndof;
	 
	if (i < currelem.GetElemType(Temp3._FEord)->GetNDofs() ) { 
           currelem.Rhs()(ip2) +=      
           currelem.GetBCDofFlag()[ip2]*dtxJxW_g*( 
                0.07*phii_gKK
	     )
	   + (1-currelem.GetBCDofFlag()[ip2])*detb*(Temp3._val_dofs[i]);
        
        currelem.Mat()(ip2,ip2) +=  (1-currelem.GetBCDofFlag()[ip2])*detb;
	}
	
	 // Matrix Assemblying ---------------------------
        for (uint j=0; j<Tempold._ndof; j++) {
          double phij_g   = currgp._phi_ndsQLVB_g[Tempold._FEord][j];
          double phij_gLL = currgp._phi_ndsQLVB_g[Temp2._FEord][j];
          double phij_gKK = currgp._phi_ndsQLVB_g[Temp3._FEord][j];
	  
        for (uint idim = 0; idim < space_dim; idim++)   {
	  dphijdx_g  [idim] = currgp._dphidxyz_ndsQLVB_g[Tempold._FEord][j+idim*Tempold._ndof]; 
	  dphijdx_gLL[idim] = currgp._dphidxyz_ndsQLVB_g[Temp2._FEord]  [j+idim*Temp2._ndof]; 
	  dphijdx_gKK[idim] = currgp._dphidxyz_ndsQLVB_g[Temp3._FEord]  [j+idim*Temp3._ndof]; 
          }
	  
	  
          double Lap_g   = Math::dot(dphijdx_g,dphiidx_g,space_dim);
          double Lap_gLL = Math::dot(dphijdx_gLL,dphiidx_gLL,space_dim);
          double Lap_gKK = Math::dot(dphijdx_gKK,dphiidx_gKK,space_dim);

	    int ip1 = i + Tempold._ndof;
	    int jp1 = j + Tempold._ndof;
	    int ip2 = i + Tempold._ndof + Temp2._ndof;
	    int jp2 = j + Tempold._ndof + Temp2._ndof;

 
//============ FIRST ROW state  delta T ===============
//======= DIAGONAL =============================
	   currelem.Mat()(i,j) +=        
            currelem.GetBCDofFlag()[i]*dtxJxW_g*( 
             Lap_g  
            );

//=========== SECOND ROW  =============
//===== DIAGONAL ===========================
 	if ( i < currelem.GetElemType(Temp2._FEord)->GetNDofs() )  { 
  	if ( j < currelem.GetElemType(Temp2._FEord)->GetNDofs() ) { 
       currelem.Mat()(ip1,jp1) +=        
            currelem.GetBCDofFlag()[ip1]*
            dtxJxW_g*(
// 	      phij_gLL*phii_gLL
           +  Lap_gLL
            ); 
	  }
	}
//============= THIRD ROW  =============
//======= DIAGONAL ==================
	if ( i < currelem.GetElemType(Temp3._FEord)->GetNDofs() )  { 
  	if ( j < currelem.GetElemType(Temp3._FEord)->GetNDofs() ) { 
          currelem.Mat()(ip2,jp2) +=        
            currelem.GetBCDofFlag()[ip2]*
              dtxJxW_g*( 
              phij_gKK*phii_gKK
            + Lap_gKK
            ); 
	  }
	}
	
	
	    
        }  //end j (col)
      }   //end i (row)
    } // end of the quadrature point qp-loop

    currelem.Mat().print_scientific(std::cout);
    
       my_system._A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
       my_system._b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());
  } // end of element loop
  // *****************************************************************

  delete [] dphijdx_g;
  delete [] dphiidx_g;
  delete [] dphijdx_gLL;
  delete [] dphiidx_gLL;
  delete [] dphijdx_gKK;
  delete [] dphiidx_gKK;
  
  }//END VOLUME
  
        my_system._A[Level]->close();
        my_system._b[Level]->close();

#ifdef DEFAULT_PRINT_INFO
  std::cout << " Matrix and RHS assembled for equation " << my_system.name()
            << " Level "<< Level << " dofs " << my_system._A[Level]->n() << std::endl;
#endif

  return;
}


