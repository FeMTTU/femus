#include "EqnT.hpp"

#include "FemusDefault.hpp"

#include "NumericVector.hpp"
#include "DenseVector.hpp"
#include "SparseMatrix.hpp"
#include "DenseMatrix.hpp"

#include "Files.hpp"
#include "Math.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "GeomEl.hpp"
#include "MultiLevelProblemTwo.hpp"
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
#include "MultiLevelProblemTwo.hpp"

// application
#include "TempQuantities.hpp"


namespace femus {

// ======================================================
EqnT::EqnT(  std::vector<Quantity*> int_map_in,
             MultiLevelProblemTwo& equations_map_in,
             std::string eqname_in,
             std::string varname_in):
    SystemTwo(int_map_in,equations_map_in,eqname_in,varname_in) {

//=======  _var_names: they are the names of the quantities which are unkwnowns to this equation  ===========
   for (uint i=0; i<int_map_in.size(); i++)  _var_names[i]=int_map_in[i]->_name;
  
}


//================ DESTRUCTOR    
      EqnT::~EqnT() {    }



 void  EqnT::GenMatRhs(const uint Level) {

  const double time =  0.;

//========== PROCESSOR INDEX
  const uint myproc = _iproc;

//========= BCHandling =========
  const double penalty_val = _mesh.GetRuntimeMap().get("penalty_val");    

  //======== ELEMENT MAPPING =======
  const uint space_dim =       _mesh.get_dim();
  const uint  meshql   = (int) _mesh.GetRuntimeMap().get("meshql");
  const uint  mesh_ord = (int) _mesh.GetRuntimeMap().get("mesh_ord");

  {//BEGIN VOLUME
  
  //==== AUXILIARY ==============
    double* dphijdx_g = new double[space_dim];
    double* dphiidx_g = new double[space_dim];
    double* dphijdx_gLL = new double[space_dim];
    double* dphiidx_gLL = new double[space_dim];
    double* dphijdx_gKK = new double[space_dim];
    double* dphiidx_gKK = new double[space_dim];

    
  const uint mesh_vb = VV;
  
  CurrentElem       currelem(VV,this,_mesh,_eqnmap._elem_type);
  CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
  
//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Tempold(currgp);
    Tempold._qtyptr   = _QtyInternalVector[0]; 
    Tempold.VectWithQtyFillBasic();
    Tempold.Allocate();

//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Temp2(currgp);
    Temp2._qtyptr   = _QtyInternalVector[1]; 
    Temp2.VectWithQtyFillBasic();
    Temp2.Allocate();

//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Temp3(currgp);
    Temp3._qtyptr   = _QtyInternalVector[2]; 
    Temp3.VectWithQtyFillBasic();
    Temp3.Allocate();
    
    //=========EXTERNAL QUANTITIES (couplings) =====
    //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);  //no quantity
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();

    //==================Quadratic domain, auxiliary, must be QUADRATIC!!! ==========
  CurrentQuantity xyz_refbox(currgp);  //no quantity
    xyz_refbox._dim      = space_dim;
    xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
    xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
    xyz_refbox.Allocate();

   const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level+1];
   const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level];

  for (uint iel=0; iel < (nel_e - nel_b); iel++) {
    
    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.set_el_nod_conn_lev_subd(Level,myproc,iel);
    currelem.SetMidpoint();

    currelem.ConvertElemCoordsToMappingOrd(xyz);
    _mesh.TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]);    

    
//MY EQUATION
//the elements are, for every level:
// 1)DOF INDICES
// 2)BC FLAGS
// 3)BC VALUES 
// 1) and 2) are taken in a single vector, 3) are considered separately
      
    currelem.SetElDofsBc(Level);

  Tempold.GetElemDofs(Level);
    Temp2.GetElemDofs(Level);
    Temp3.GetElemDofs(Level);
    
    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),Tempold._FEord,currelem.GetBCDofFlag()); //only the Qtyzero Part is modified!

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
   const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();

    for (uint qp=0; qp< el_ngauss; qp++) {

//======= "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)   { 
  currgp.SetPhiElDofsFEVB_g (fe,qp);
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp); 
  }
	  
const double      det = currgp.JacVectVV_g(xyz);
const double dtxJxW_g = det*_eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
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
	 
	if (i < _eqnmap._elem_type[currelem.GetDim()-1][ Temp2._FEord ]->GetNDofs() ) { 
	 currelem.Rhs()(ip1) +=      
           currelem.GetBCDofFlag()[ip1]*dtxJxW_g*( 
                0.07*phii_gLL
	  )
	   + (1-currelem.GetBCDofFlag()[ip1])*detb*(Temp2._val_dofs[i]);
        
         currelem.Mat()(ip1,ip1) +=  (1-currelem.GetBCDofFlag()[ip1])*detb;
	}
	
//======= THIRD ROW ===================================
	 int ip2 = i + Tempold._ndof + Temp2._ndof;
	 
	if (i < _eqnmap._elem_type[currelem.GetDim()-1][ Temp3._FEord ]->GetNDofs() ) { 
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
 	if ( i < _eqnmap._elem_type[currelem.GetDim()-1][ Temp2._FEord ]->GetNDofs() )  { 
  	if ( j < _eqnmap._elem_type[currelem.GetDim()-1][ Temp2._FEord ]->GetNDofs() ) { 
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
	if ( i < _eqnmap._elem_type[currelem.GetDim()-1][ Temp3._FEord ]->GetNDofs() )  { 
  	if ( j < _eqnmap._elem_type[currelem.GetDim()-1][ Temp3._FEord ]->GetNDofs() ) { 
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
    
       _A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
       _b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());
  } // end of element loop
  // *****************************************************************

  delete [] dphijdx_g;
  delete [] dphiidx_g;
  delete [] dphijdx_gLL;
  delete [] dphiidx_gLL;
  delete [] dphijdx_gKK;
  delete [] dphiidx_gKK;
  
  }//END VOLUME
  
  {//BEGIN BOUNDARY  // *****************************************************************
  
   //-----Nonhomogeneous Neumann------
 // Qflux = - k grad(T) by definition
//  QfluxDOTn>0: energy flows outside (cooling)  QfluxDOTn<0: energy flows inside (heating)
    double* Qflux_g = new double[space_dim];
    
  const uint mesh_vb = BB;
  
  CurrentElem       currelem(BB,this,_mesh,_eqnmap._elem_type);
  CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
  
//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Tempold(currgp);
    Tempold._qtyptr   = _QtyInternalVector[0]; 
    Tempold.VectWithQtyFillBasic();
    Tempold.Allocate();

//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Temp2(currgp);
    Temp2._qtyptr   = _QtyInternalVector[1]; 
    Temp2.VectWithQtyFillBasic();
    Temp2.Allocate();

//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Temp3(currgp);
    Temp3._qtyptr   = _QtyInternalVector[2]; 
    Temp3.VectWithQtyFillBasic();
    Temp3.Allocate();
    
    //=========EXTERNAL QUANTITIES (couplings) =====
    //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);  //no quantity
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();

    //==================Quadratic domain, auxiliary, must be QUADRATIC!!! ==========
  CurrentQuantity xyz_refbox(currgp);  //no quantity
    xyz_refbox._dim      = space_dim;
    xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
    xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
    xyz_refbox.Allocate();

   const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level+1];
   const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level];
   
     for (uint iel=0;iel < (nel_e - nel_b) ; iel++) {

      currelem.Mat().zero();
      currelem.Rhs().zero();

      currelem.set_el_nod_conn_lev_subd(Level,myproc,iel);
      currelem.SetMidpoint(); 

      currelem.ConvertElemCoordsToMappingOrd(xyz);
    _mesh.TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]);    
     
      currelem.SetElDofsBc(Level);
      
       Tempold.GetElemDofs(Level);
         Temp2.GetElemDofs(Level);
         Temp3.GetElemDofs(Level);

     if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),Tempold._FEord,currelem.GetBCDofFlag()); //only the Quadratic Part is modified!
  
 //============ FLAGS ================
     double el_penalty = 0.;
     int pen_sum=0;
     for (uint i=0; i< Tempold._ndof; i++)   pen_sum += currelem.GetBCDofFlag()[i];
//pen_sum == 0: all the nodes must be zero, so that ALL THE NODES of that face are set properly, even if with a penalty integral and not with the nodes!
//this is a "FALSE" NATURAL boundary condition, it is ESSENTIAL actually
     if (pen_sum == 0/*< el_n_dofs porcata*/) { el_penalty = penalty_val;   }  //strictly minor for Dirichlet penalty, i.e. AT LEAST ONE NODE to get THE WHOLE ELEMENT
                                                           //also for Neumann flag is the same
                                                           
//you should check that when you impose nodal bc flags you have to involve exactly one element for bc==0...
//but on the other hand you will not impose bc==1 on all nodes for the element
//So, it is recommended that for Dirichlet conditions you set ALL the NODES of that element to ZERO

//in order to do the flag here, since it is a "TRUE" NATURAL BOUNDARY CONDITION,
//it only suffices that SOME OF THE NODES ARE with bc=1, AT LEAST ONE
int el_Neum_flag=0;
     uint Neum_sum=0;
     for (uint i=0; i < Tempold._ndof; i++)   Neum_sum += currelem.GetBCDofFlag()[i];
     for (uint i=0; i < Tempold._ndof; i++)   Neum_sum += currelem.GetBCDofFlag()[i + Tempold._ndof];
            if ( Neum_sum == 2*Tempold._ndof )  { el_Neum_flag=1;  }

//====================================

   const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();

    for (uint qp=0; qp< el_ngauss; qp++) {

//======= "COMMON SHAPE PART"============================
  for (uint fe = 0; fe < QL; fe++)  {
    currgp.SetPhiElDofsFEVB_g (fe,qp);
    currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp); 
  }
        const double  det   = currgp.JacVectBB_g(xyz);
        const double dtxJxW_g = det * _eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
//=======end "COMMON SHAPE PART"===================================

       xyz.val_g();
       
       static_cast<Temperature*>(_eqnmap._qtymap.get_qty("Qty_Temperature"))->heatflux_txyz(time,&xyz._val_g[0],Qflux_g);
   
	Tempold.val_g(); //For the penalty Dirichlet //i need this for interpolating the old function at the gauss point

	   double QfluxDn_g=Math::dot( Qflux_g,currgp.get_normal_ptr(),space_dim );
	 
        for (uint i=0; i<Tempold._ndof; i++) {
   	const double phii_g =  currgp._phi_ndsQLVB_g[Tempold._FEord][i]; 
	
       currelem.Rhs()(i) +=
          0.*(currelem.GetBCDofFlag()[i]*
         el_Neum_flag*dtxJxW_g*(-QfluxDn_g)*phii_g    // beware of the sign  //this integral goes in the first equation
	 + el_penalty*dtxJxW_g*Tempold._val_g[0]*phii_g)  //clearly, if you continue using bc=0 for setting nodal Dirichlet, this must go outside
	 ; 
	 
         if (_Dir_pen_fl == 1) {
            for (uint j=0; j<Tempold._ndof; j++) {
               double phij_g = currgp._phi_ndsQLVB_g[Tempold._FEord][j];
	       currelem.Mat()(i,j) += 0.*el_penalty*dtxJxW_g*phij_g*phii_g;
	    } 
          }

	  //end of j loop
	}
          // end of i loop
    }
        // end BDRYelement gaussian integration loop
        
         _A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
         _b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());
   
  }
      // end of BDRYelement loop
    

  }//END BOUNDARY

  
#ifdef DEFAULT_PRINT_INFO
  std::cout << " Matrix and RHS assembled for equation " << _eqname
            << " Level "<< Level << " dofs " << _A[Level]->n() << std::endl;
#endif

  return;
}


} //end namespace femus
