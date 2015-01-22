#include "EqnT.hpp"

#include "FemusDefault.hpp"

#include "NumericVector.hpp"
#include "DenseVector.hpp"
#include "SparseMatrix.hpp"
#include "DenseMatrix.hpp"
#include "LinearEquationSolver.hpp"

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
#include "OptLoop.hpp"
#include "paral.hpp"

// application
#include "Temp_conf.hpp"
#include "TempQuantities.hpp"
#include "EqnNS.hpp"



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
// ======================================================
EqnT::EqnT(  const TimeLoop & time_loop_in,
	     std::vector<Quantity*> int_map_in,
             MultiLevelProblemTwo& equations_map_in,
             std::string eqname_in,
             std::string varname_in):
    SystemTwo(int_map_in,equations_map_in,eqname_in,varname_in),
    _my_timeloop(time_loop_in) {

//=======  _var_names[]  ===========
  _var_names[0]="T";
  _var_names[1]="Tlift";
  _var_names[2]="Tadj";

//========= MG solver ===================
  for (uint l=0;l<_NoLevels;l++)  _solver[l]->set_solver_type(SOLVERT);
  
//============= DIR PENALTY===============
   _Dir_pen_fl = TEMP_DIR_PENALTY;
  
}


//================ DESTRUCTOR    
      EqnT::~EqnT() {    }


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

//TODO  what happens for STANDARD OUTPUTS? can we do in such a way that EVERYTHING is printed TO FILE?
// we should REDIRECT TO THE *SAME* FILE ALL THE std output and std errors of ALL THE LIBRARIES!

//NOW, PAY ATTENTION: The "iel" written as "iel=0; iel < (nel_e - nel_b);" is used for PICKING the CONNECTIVITY from the ELEMENT CONNECTIVITY MAP!
// But, the iel as DofObject Index must be given in the correct form!
// So, I will distinguish iel into iel_mesh and iel_DofObj:
//In both cases we start from a "Geometrical Entity", the ELEMENT.
//In the mesh case, we only use the ELEMENT to pick its CONNECTIVITY.
//In the DofObj case,we use iel to pick the corresponding Element DOF from the node_dof map! 


/// This function assembles the matrix and the rhs:
void  EqnT::GenMatRhs(const uint Level) {

  const double time = 0.; // _eqnmap._timeloop._curr_time;

//========== PROCESSOR INDEX
  const uint myproc = _iproc;

//==========FLAG FOR STATIONARITY OR NOT
  const double    dt = 1.; //_eqnmap._timeloop._timemap.get("dt");
  const uint Nonstat = _phys.get("NonStatTEMP");
  
//========= BCHandling =========
  const double penalty_val = _mesh.GetRuntimeMap().get("penalty_val");    

  //======== ELEMENT MAPPING =======
  const uint space_dim =       _mesh.get_dim();
  const uint  meshql   = (int) _mesh.GetRuntimeMap().get("meshql");
  const uint  mesh_ord = (int) _mesh.GetRuntimeMap().get("mesh_ord");

  //====== reference values ========================
  const double IRe = 1./_phys.get("Re");
  const double IPr = 1./_phys.get("Pr");
  const double alphaT  = _phys.get("alphaT");
  const double alphaL2 = _phys.get("alphaL2");
  const double alphaH1 = _phys.get("alphaH1");
  
  //==== AUXILIARY ==============
    double* dphijdx_g = new double[space_dim];
    double* dphiidx_g = new double[space_dim];
   //-----Nonhomogeneous Neumann------
 // Qflux = - k grad(T) by definition
//  QfluxDOTn>0: energy flows outside (cooling)  QfluxDOTn<0: energy flows inside (heating)
    double* Qflux_g = new double[space_dim];

    int T4_ord = T4_ORD;
    
// ==========================================  
// ==========================================  
 {//BEGIN VOLUME
 const uint mesh_vb = VV;
  
  CurrentElem       currelem(VV,this,_mesh,_eqnmap._elem_type);    
  CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
  

//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Tempold(currgp);
    Tempold._qtyptr   = _QtyInternalVector[0]; 
    Tempold.VectWithQtyFillBasic();
    Tempold.Allocate();

//====================================
    CurrentQuantity Tlift(currgp);
    Tlift._qtyptr   = _eqnmap._qtymap.get_qty("Qty_TempLift");//_QtyInternalVector[1]; 
    Tlift.VectWithQtyFillBasic();
    Tlift.Allocate();

//=====================================
    CurrentQuantity TAdj(currgp);
    TAdj._qtyptr   = _eqnmap._qtymap.get_qty("Qty_TempAdj");//_QtyInternalVector[2]; 
    TAdj.VectWithQtyFillBasic();
    TAdj.Allocate();
    
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
  
  //==================
    CurrentQuantity vel(currgp);
    vel._qtyptr   = _eqnmap._qtymap.get_qty("Qty_Velocity"); 
    vel.VectWithQtyFillBasic();
    vel.Allocate();
    
//===============Tdes=====================
    CurrentQuantity Tdes(currgp);
    Tdes._qtyptr   = _eqnmap._qtymap.get_qty("Qty_TempDes"); 
    Tdes.VectWithQtyFillBasic();
    Tdes.Allocate();

  

   const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level+1];
   const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level];

// ==========================================  
// ==========================================  
   

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
    Tlift.GetElemDofs(Level);
     TAdj.GetElemDofs(Level);
     
    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),Tempold._FEord,currelem.GetBCDofFlag()); //only the Qtyzero Part is modified!

// ===============      
// Now the point is this: there are several functions of space
// which are expressed with respect to a reference frame
//ok, now that the dofs are filled for xyz_refbox, I can use the el_average
//Well, the alternative is to consider  the elem in the refbox as
    //either a Vect or a CurrentElem !
    //I could consider it as another element, but only with the geometrical part!

  xyz_refbox.SetElemAverage();
  
int domain_flag = ElFlagControl(xyz_refbox._el_average,&_mesh);
//====================    
    
//===== FILL the DOFS of the EXTERNAL QUANTITIES: you must assure that for every Vect the quantity is set correctly
// for every Vect it must be clear if it belongs to an equation or not, and which equation it belongs to;
// this is usually made clear by the related QUANTITY.
// Now the main thing to check is the difference between Vect WITH QUANTITY and Vect WITHOUT QUANTITY.
  // it is better to avoid using GetElDofs if the Vect is internal, only if external  
  //Do not use GetElDofs if you want to pick an intermediate dof...
      
   if ( vel._eqnptr != NULL )  vel.GetElemDofs(Level);
   else                        vel._qtyptr->FunctionDof(vel,time,&xyz_refbox._val_dofs[0]);

   if ( Tdes._eqnptr != NULL )  Tdes.GetElemDofs(Level);
   else                         Tdes._qtyptr->FunctionDof(Tdes,time,&xyz_refbox._val_dofs[0]);


   const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();
   
   for (uint qp=0; qp< el_ngauss; qp++) {

//======= "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)   { 
  currgp.SetPhiElDofsFEVB_g (fe,qp);
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp); 
}
	  
const double      det = dt*currgp.JacVectVV_g(xyz);
const double dtxJxW_g = det*_eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { 
  currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp);
  currgp.ExtendDphiDxyzElDofsFEVB_g(fe);
}
//======= end of the "COMMON SHAPE PART"==================

 	Tempold.val_g(); 
          Tlift.val_g(); 
           TAdj.val_g(); 
            vel.val_g(); 
           Tdes.val_g();
	   
	   // always remember to get the dofs for the variables you use!
           // The point is that you fill the dofs with different functions...
           // you should need a flag to check if the dofs have been correctly filled

      /// d) Local (element) assemblying energy equation
      for (uint i=0; i < Tempold._ndof; i++)     {

        const double phii_g = currgp._phi_ndsQLVB_g[Tempold._FEord][i];

        for (uint idim = 0; idim < space_dim; idim++) dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[Tempold._FEord][i+idim*Tempold._ndof];

//=========== FIRST ROW ===============
        currelem.Rhs()(i) +=      
           currelem.GetBCDofFlag()[i]*dtxJxW_g*( 
                Nonstat*Tempold._val_g[0]*phii_g/dt
	  )
	   + (1-currelem.GetBCDofFlag()[i])*detb*(Tempold._val_dofs[i]);
        
        currelem.Mat()(i,i) +=  (1-currelem.GetBCDofFlag()[i])*detb;

//========= SECOND ROW (CONTROL) =====================
	 int ip1 = i + /* 1* */Tempold._ndof;   //suppose that T' T_0 T_adj have the same order
	 currelem.Rhs()(ip1) +=      
           currelem.GetBCDofFlag()[ip1]*dtxJxW_g*( 
                Nonstat*Tempold._val_g[0]*phii_g/dt
                     + alphaT*domain_flag*(Tdes._val_g[0])*phii_g // T_d delta T_0    /////// ADDED /////
	  )
	   + (1-currelem.GetBCDofFlag()[ip1])*detb*(Tlift._val_dofs[i]);
        
         currelem.Mat()(ip1,ip1) +=  (1-currelem.GetBCDofFlag()[ip1])*detb;

//======= THIRD ROW (ADJOINT) ===================================
	 int ip2 = i + 2 * Tempold._ndof;   //suppose that T' T_0 T_adj have the same order
           currelem.Rhs()(ip2) +=      
           currelem.GetBCDofFlag()[ip2]*dtxJxW_g*( 
                Nonstat*Tempold._val_g[0]*phii_g/dt
                + alphaT*domain_flag*(Tdes._val_g[0])*phii_g // T_d delta T'
	     )
	   + (1-currelem.GetBCDofFlag()[ip2])*detb*(Tempold._val_dofs[i]);
        
        currelem.Mat()(ip2,ip2) +=  (1-currelem.GetBCDofFlag()[ip2])*detb;

#if FOURTH_ROW==1
	 int ip3 = i + 3 * Tempold._ndof;   //suppose that T' T_0 T_adj have the same order
	 
	 if (i < _eqnmap._elem_type[currelem.GetDim()-1][T4_ord]->GetNDofs()) { currelem.Rhs()(ip3) +=  currelem.GetBCDofFlag()[ip3]*dtxJxW_g*(currgp._phi_ndsQLVB_g[T4_ord][i]) + (1-currelem.GetBCDofFlag()[ip3])*detb*1300.;
	              currelem.Mat()(ip3,ip3)  += ( 1-currelem.GetBCDofFlag()[ip3] )*detb;  }
#endif
	 // Matrix Assemblying ---------------------------
        for (uint j=0; j<Tempold._ndof; j++) {
          double phij_g = currgp._phi_ndsQLVB_g[Tempold._FEord][j];
	  
        for (uint idim = 0; idim < space_dim; idim++)  dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[Tempold._FEord][j+idim*Tempold._ndof]; 
           
   
          double Lap_g   = Math::dot(dphijdx_g,dphiidx_g,space_dim);
          double Advection = Math::dot(&vel._val_g[0],dphijdx_g,space_dim);

	    int ip1 = i + Tempold._ndof;
	    int jp1 = j + Tempold._ndof;
	    int ip2 = i + 2*Tempold._ndof;
	    int jp2 = j + 2*Tempold._ndof;

// 	           T     T_0     T_adj
	    
// 	    T      X      X       O
	     
// 	    T_0   
	    
// 	    T_adj
	    
	    
//============ FIRST ROW state  delta T ===============
//======= DIAGONAL =============================
	   currelem.Mat()(i,j) +=        
            currelem.GetBCDofFlag()[i]*dtxJxW_g*( 
              Nonstat*phij_g*phii_g/dt 
            + Advection*phii_g
            + IRe*IPr*Lap_g  
            );

//===============================
    //same operators for T and T_0
	    currelem.Mat()(i,jp1) +=        
            currelem.GetBCDofFlag()[i]*dtxJxW_g*(    
              Nonstat*phij_g*phii_g/dt 
            + Advection*phii_g
            + IRe*IPr*Lap_g
	    );

//====================================
	   currelem.Mat()(i,jp2) +=        
            currelem.GetBCDofFlag()[i]*dtxJxW_g*( 
                0.
            );

	    
//============= SECOND ROW (LIFTING) delta T_0 =============
//===== DIAGONAL ===========================
         currelem.Mat()(ip1,jp1) +=        
            currelem.GetBCDofFlag()[ip1]*
            dtxJxW_g*( 
              Nonstat*phij_g*phii_g/dt
             + alphaL2*phij_g*phii_g  //L_2 control norm
             + alphaH1*Lap_g          //H_1 control norm
              + alphaT*domain_flag*(phij_g)*phii_g  //T_0 delta T_0  //ADDED///////////////
            ); 
//====================================
	   currelem.Mat()(ip1,j) +=        
            currelem.GetBCDofFlag()[ip1]*
            dtxJxW_g*( 
                + alphaT*domain_flag*(phij_g)*phii_g  //T' delta T_0     //ADDED///////////////
            );
//====================================
	   currelem.Mat()(ip1,jp2) +=        
            currelem.GetBCDofFlag()[ip1]*
             dtxJxW_g*( 
                 -Advection*phii_g
                + IRe*IPr*Lap_g
           );

//============= THIRD ROW (ADJOINT) =============
//======= DIAGONAL ==================
          currelem.Mat()(ip2,jp2) +=        
            currelem.GetBCDofFlag()[ip2]*
              dtxJxW_g*( 
              Nonstat*phij_g*phii_g/dt
            - Advection*phii_g  //minus sign
            + IRe*IPr*Lap_g
               
            ); 
//====================================
	   currelem.Mat()(ip2,j) +=        
            currelem.GetBCDofFlag()[ip2]*
            dtxJxW_g*( 
               + alphaT*domain_flag*(phij_g)*phii_g  //T' delta T'
            );
//====================================
	   currelem.Mat()(ip2,jp1) +=        
            currelem.GetBCDofFlag()[ip2]*
            dtxJxW_g*( 
               + alphaT*domain_flag*(phij_g)*phii_g  //T_0 delta T'     ///ADDED///////
            );

#if FOURTH_ROW==1
	 int ip3 = i + 3*Tempold._ndof;   //suppose that T' T_0 T_adj have the same order
// 	    int jp3 = j + 3*Tempold._ndof;
	   if (i < _eqnmap._elem_type[currelem.GetDim()-1][T4_ord]->GetNDofs() ) currelem.Mat()(ip3,ip3) += currelem.GetBCDofFlag()[ip3]*dtxJxW_g*(currgp._phi_ndsQLVB_g[ T4_ord ][/*j*/i]*currgp._phi_ndsQLVB_g[ T4_ord ][i]);   
#endif
	    
        }  //end j (col)
      }   //end i (row)
    } // end of the quadrature point qp-loop

       _A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
       _b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());
  } // end of element loop
  // *****************************************************************

 
   }//END VOLUME
  
   { //BEGIN BOUNDARY  // *****************************************************************
   
 const uint mesh_vb = BB;
  
  CurrentElem       currelem(BB,this,_mesh,_eqnmap._elem_type);    
  CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
  

//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    CurrentQuantity Tempold(currgp);
    Tempold._qtyptr   = _QtyInternalVector[0]; 
    Tempold.VectWithQtyFillBasic();
    Tempold.Allocate();

//====================================
    CurrentQuantity Tlift(currgp);
    Tlift._qtyptr   = _eqnmap._qtymap.get_qty("Qty_TempLift");//_QtyInternalVector[1]; 
    Tlift.VectWithQtyFillBasic();
    Tlift.Allocate();

//=====================================
    CurrentQuantity TAdj(currgp);
    TAdj._qtyptr   = _eqnmap._qtymap.get_qty("Qty_TempAdj");//_QtyInternalVector[2]; 
    TAdj.VectWithQtyFillBasic();
    TAdj.Allocate();
    
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
    
//===============Tdes=====================
    CurrentQuantity Tdes(currgp);
    Tdes._qtyptr   = _eqnmap._qtymap.get_qty("Qty_TempDes"); 
    Tdes.VectWithQtyFillBasic();
    Tdes.Allocate();

  

   const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level+1];
   const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level];

// ==========================================  
// ==========================================     
   
    
     for (uint iel=0;iel < (nel_e - nel_b) ; iel++) {

      currelem.Mat().zero();
      currelem.Rhs().zero();

      currelem.set_el_nod_conn_lev_subd(Level,myproc,iel);
      currelem.SetMidpoint(); 

      currelem.ConvertElemCoordsToMappingOrd(xyz);
    _mesh.TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]);    
     
      currelem.SetElDofsBc(Level);
      
       Tempold.GetElemDofs(Level);
         Tlift.GetElemDofs(Level);
          TAdj.GetElemDofs(Level);

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
        const double  det   = dt*currgp.JacVectBB_g(xyz);
        const double dtxJxW_g = det * _eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
//=======end "COMMON SHAPE PART"===================================

       xyz.val_g();
       
       static_cast<Temperature*>(_eqnmap._qtymap.get_qty("Qty_Temperature"))->heatflux_txyz(time,&xyz._val_g[0],Qflux_g);
   
	Tempold.val_g(); //For the penalty Dirichlet //i need this for interpolating the old function at the gauss point

	   double QfluxDn_g=Math::dot( Qflux_g,currgp.get_normal_ptr(),space_dim );
	 
        for (uint i=0; i<Tempold._ndof; i++) {
   	const double phii_g =  currgp._phi_ndsQLVB_g[Tempold._FEord][i]; 
	
       currelem.Rhs()(i) +=
          currelem.GetBCDofFlag()[i]*
         el_Neum_flag*dtxJxW_g*(-QfluxDn_g)*phii_g    // beware of the sign  //this integral goes in the first equation
	 + el_penalty*dtxJxW_g*Tempold._val_g[0]*phii_g  //clearly, if you continue using bc=0 for setting nodal Dirichlet, this must go outside
	 ; 
	 
         if (_Dir_pen_fl == 1) {
            for (uint j=0; j<Tempold._ndof; j++) {
               double phij_g = currgp._phi_ndsQLVB_g[Tempold._FEord][j];
	       currelem.Mat()(i,j) += el_penalty*dtxJxW_g*phij_g*phii_g;
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


//======================
      //TODO here I have to put the right offset back
     //This function must receive a BOUNDARY ELEMENT according to the common "FEMuS OFFSET NUMERATION LevxSubd",
     //and it must yield the VOLUME ELEMENT NUMBER again according to the common "FEMuS OFFSET NUMERATION LevxSubd",
     //so that it has embedded the information of "BELONGING TO SOME GIVEN SUBDOMAIN"
     //We start from a Boundary element belonging to some subdomain, and we want to get the CORRESPONDING VOLUME ELEMENT NUMBER,
     //which is the volume element number IN THE "FEMuS OFFSET NUMERATION LevxSubd" NUMERATION,
     //so that once you have it you already know to which processor you belong!
     // Now this vol_iel that we have is already ok, because it is in the ABSOLUTE FEMuS VOLUME ELEMENT ORDERING
     // Actually, in the mesh file generated by GenCase we never wrote so far the list of numbers of the elements,
     // we only needed to write the offsets,
     //then we reconstruct the numbering later,
     //then we reassociate that numbering to the connectivity.
     //It is like the connectivities were "orphan" of their respective numbers, but actually the connectivities 
     //were printed exactly in the new Femus ordering,
     //(as a matter of fact, the FEMuS mesh file forgets about the libmesh ordering but only considers the femus one),
     //so actually their order INTRINSICALLY GIVES the Femus ordering.
     
     //Now there is another story: the _el_bdry_to_vol returns an element number in the 
     // Femus ABSOLUTE Element numbering.
     // But, the ABSOLUTE element numbering is used to pick the CONNECTIVITIES.
     //To go from the ABSOLUTE to the RELATIVE Femus element numbering,
     //you have to do like this:
     // first go from the RELATIVE boundary to the ABSOLUTE boundary
     //then from the absolute boundary you get the ABSOLUTE VOLUME element number
     //finally, from the ABSOLUTE VOLUME YOU WANT TO GET THE DofObject
     //TODO the thing that does not seem very nice to me is that 
     //there is an _el_bdry_to_vol map for every LEVEL, 
     //but the elements are NOT numbered by level, but according to
     //the ABSOLUTE element numbering
     
      
      //this iel_femus is ABSOLUTE, but PER LEVEL
     // the following is ABSOLUTE, not even per level, COMPLETELY ABSOLUTE
     //maybe it would be better to have it absolute but per level
//======================
      
      
//======================
// Now the dof indices on the boundary must match with the 
//dofs on the volume
//If with constant elements you do not have any dof on the boundary,
//then you have to think like this:
//What is the meaning of ENFORCING a BOUNDARY CONDITION VALUE?
//The meaning is that you enforce the TRACE of the VOLUME SHAPE FUNCTIONS.
//Then, of course, enforcing the TRACE means enforcing the VOLUME DOF, for constant FE.
//With nodes on the boundary, you are enforcing the DOF ITSELF, because 
//among all the nodes there are some which are ONLY BOUNDARY NODES.
//Again, enforcing those boundary nodes can be seen as a form of ENFORCING the TRACE of the 
//VOLUME SHAPE FUNCTIONS associated to those BOUNDARY NODES.
//In fact, you enforce the DOF values so that you "stretch" your TRACE at the boundary.

//Now, the point here is that we have a VOLUME DOF, iel, associated to the CONSTANT FE,
//which may not belong to the current row actually
//Wait, so, do we need to have the VOLUME iel in the list of BOUNDARY DOF INDICES?
//"Boundary dof indices" are the 
// "dof_indices that are involved in considering the current BOUNDARY ELEMENT"
//DOF_INDICES means "ROWS and COLUMNS" of the matrix

//I dont care if i am in a BOUNDARY INTEGRAL and the dofs I involve
//are not "boundary dofs" strictly speaking...
//I am only interested in INVOLVING ALL THE DOFS THAT ARE INVOLVED .
//Now, I need a map that for every boundary element
//gives me the corresponding VOLUME element in the mesh.
//For every volume element,  i can give you the children, which can be one or two...
//for every boundary element, there is only one father,
//so let us do a simple map that goes from CHILDREN to FATHER.
//That is not gonna be in terms of DOF but in terms of MESH.
//So we go from BOUNDARY iel to VOLUME iel.
//We need to do this in the Gencase because it is there
//where we can use the Libmesh functions







	 //on the other hand consider that the penalty term does not have to appear whenever we are inside the domain
//the problem is that being inside the domain doesnt mean being Dirichlet, but here bc=0 may mean either being inside the domain
//or being a Dirichlet node. We want to put the penalty term where bc=0,but not inside the domain!!!
//so it takes a flag for setting dirichlet
//a flag for setting that we are inside the domain...no,but we already know that we are on the boundary, because
//the problem is that some nodes are BOTH VOLUME AND BOUNDARY. Or better, the SOME VOLUME NODES are also BOUNDARY NODES...
//no no


 //this should be removed I think: yes it should,
 //I mean its useless. With the PENALTY APPROACH you dont have a flag 0-1, but 0-\infty!
 //therefore, you only have to distinguish two types of elements:
 // - the ELEMENTS in which penalty=\infty, where the Dirichlet integral dominates over everything
 // - the ELEMENTS with penalty=0, for which we know that bc=1 so the bc_alldofs in front is not needed!
 //in the first ones the Dirichlet integral, AND ONLY IT,
 //in the second ones the Neumann integral, TOGETHER WITH THE REST OF THE EQUATION
 //instead in the NODAL DIRICHLET BCs, Neumann is still ELEMENTAL, so we must implement it 
 //with Neum_flag, NOT WITH bc_alldofs again!
 //in conclusion, bc_alldofs must be eliminated FOREVER, because it's NODAL, and no nodal thing must be used for the boundary INTEGRALS!
 //what remains is Neum_flag when you dont use penalty
 //Well,actually we HAVE to use bc_alldofs[i] in case of NODAL DIRICHLET. In fact, in that case
 //ALL THE ROWS OF THE LINEAR SYSTEM must be multiplied by bc_alldofs, because
 // the bc_qldof flag means "in this row we are putting or not a Dirichlet nodal value"
 // instead, el_Neum_flag means "For the ELEMENT we consider we must compute the integral"
//                                  (and this computation will involve more than one row, so we cannot isolate the rows)
//clearly, this situation creates a little mismatch, because it may happen
//that the nodal bc-s are 1-1-0,
//but the integral should be computed over ALL the element
//it seems like you have to interpolate with three functions but you stop the sum to 2 instead.
//The missing term should appear in equation "i" in some column "j", but equation "i" in the nodal case
//has all zero except 1 on the diagonal
//Removing that node corresponds to TRUNCATING the INTERPOLATION of the SHAPE FUNCTIONS (i index) at that gauss point.
//That doesnt give any problem because the "coefficient that doesn't multiply anything" would have multiplied 
// a zero, so no fear!
	 
	 
	 
	 
	 
	 
	 //PAY ATTENTION here because there may be a slight problem if some errors in the update of Tempold._val_g are there,
//like it happened with p_old for pressure in NS, which suggested us to use a FUNCTION instead of the previous value

//el_Neum_flag and el_penalty are clearly ALTERNATIVE: 
// 	 el_Neum_flag == 1 el_penalty == 0
// 	 el_Neum_flag == 0 el_penalty == 10e20
//This means that you shouldnt need multiplying by el_Neum_flag in the Neumann integral.
//Well actually I would like to keep some symmetry...
//In practice, the basis is writing the general equation... Then the issue is enforcing Dirichlet.
//In the NODAL case it is bc_alldofs that rules.
//In the PENALTY case it is el_penalty that rules
	//So, in both cases computing  el_Neum_flag is useless
	 //in the nodal case, el_Neum_flag is 1 when bc is 1
//You may want to use el_Neum_flag to avoid adding the small Neum integral to the 10e20 row, 
//if you're very fond of symmetry...
//yes, it is a minor thing but it is harmless
//well, actually it might be also good when the numbers of the heat flux are very high...
//So we'll leave it, it looks as stabilyzing...

// T' should be equal to -T_0
//The point is that T' is of course homogeneous at the boundary, by definition.
//The part that takes into account the boundary conditions is only T_0:
//now try to set the bc free for T': in this case you get exactly T' = -T_0;
// but this is not what you should get.
//Of course, if you have different boundary conditions for T' and T_0,
//they can never be equal...
//With the magnetic field we didnt have b' = - B_e
//now try to solve for different T_0 (change the equation, put the laplacian, get sthg else...)
//is it true that different T_0 of the SAME boundary conditions lead to the same solution T' + T_0 ?
//Given different lifting functions of the same boudnarty conditions,
//what is the condition sucj that we do not depend on the particular lifting function?
// For the magnetic field we had div B = 0, and we claimed that in that case 
//the final sum was the same
//if we had div B != 0, the final result in terms of u and p was different
//maybe that case was particular only for Hartmann, which is a linear, simplified, MHD problem?

//our previous claim was: the NS + MHD operators where in such a form that their action 
// on various B_e was leading to results that were only dependent on the boundary conditions.
// We should check this thing more THEORETICALLY: considering how the operators act on the 
//decomposition b+ B_e , where B_e  is or not divergence-free.

//in order to find the equivalence, i think that we have to check the OPERATOR we have.
//if the operator is not sensible to some class of functions, or better if the action of that 
//operator does not "filter out" the particular choice of the lifting function 
//                     and only retain the boundary condition values,
//then that is not a good operator for us.
// we need a sort of "independent operator", or "operator of boundary conditions"
//LAPLACIAN
//LAPLACIAN + TRNSPORT

// The point is: if the solution of the operator for T is unique,
// then whatever T_0 you set you get another T', but the sum (T' + T_0) is unique.
// And, of course, if the boundary conditions are the same, the sum (T' + T_0) is the same. 

//Now with the lifting of B you ALSO needed that the lifting function was DIVERGENCE-FREE:
// otherwise, you were getting different solutions for (u,p).
//But, in this case you dont have any particular constraint on T_0

//For the temperature, you may say that if the solution of the state equations is unique
//for a given boundary condition, then whatever T_0 you get, the sum T' + T_0 is unique.

// Numerically speaking, the only point would be to be sure that the matrix is diagonal dominant
// (coming from coercivity) for every lifting function: clearly we know that coercivity 
// is fulfilled only for certain lifting functions, for a certain epsilon value wrt the viscosity.

// PROBLEMA: Boundary conditions for the STATE, ADJOINT AND CONTROL variables
// Secondo la mia teoria, sia per T' sia per la sua aggiunta le condizioni al boundary devono essere FISSE.

// Come e' giusto per un controllo d boundary, l'effetto prevalente deve avvenire vicino al boundary. 
// Quindi il controllo di boundary agisce sul boundary del dominio, pertanto non puo' spingere molto all'interno.

//$$$$$$$$$$$$$$$$ Come posso aumentare la spinta all'interno?

//$$$$$$$$$$$$$ Il ruolo delle boundary conditions: slego T_0 al boundary, pero' le altre non sono slegate 
// al boundary. Non slegherei T', magari proverei a slegare l'aggiunta, 
// che da la spinta al controllo, che da la spinta allo stato

//$$$$$$$$$$$$$$$ Provare a cambiare con il target, o con la condizione iniziale di Stato/Controllo

//$$$$$$$$$$ E' chiaro che il controllo di boundary non e' cosi' efficace come il controllo distribuito

//CONDIZIONI di NEUMANN con il metodo del lifting: come si fanno a mettere?
//Ad esempio nella parte non controllata vorrei poter mettere dei pezzi adiabatici.
//Allora qual e' il punto: 
//    l'imposizione del controllo e l'imposizione della condizione di neumann
//    avverrebbero nello stesso modo,
//    cioe' non fissando il valore in quel nodo.
//    Allora come fa il sistema a capire se un pezzo di boundary della lifting function
//    e' un pezzo di CONTROLLO o un pezzo di NEUMANN OMOGENEO ?

// Aspetta: se tu nell'equazione della temperatura stai lasciando tutto libero e 
// non stai fissando l'integrale di boundary, allora stai fissando tutto uguale a zero...
//Poi pero' il pezzo dell'integrale di boundary per T_0 va a finire nell'equazione del controllo 
// quando fai la variazione delta T_0, quindi nell'equazione del controllo avresti 
//   l'integrale di boundary dell'aggiunta...
// e allora la domanda e': Qual e' la condizione di Neumann dell'aggiunta?



//Poi, come si fa a tenere la temperatura positiva? Un algoritmo serio dovrebbe evitare temperature negative...
// Se pero' l'equazione e' adimensionale, i valori adimensionali possono venire anche negativi, 
// dipende da come hai adimensionalizzato: l'importante e' che poi i valori finali siano positivi
// SE IO AGGIUNGO AI VINCOLI ANCHE UN VINCOLO DI DISUGUAGLIANZA DICENDO CHE LA TEMPERATURA SIA POSITIVA?



// Succede questo (quasi paradossale): se la regione di controllo e' vicino al boundary control,
// converge fino ad alpha massimo molto piccolo.
// Se la regione di controllo e' lontana dal boundary control, allora converge fino ad un alpha max molto piu' grande.

