#include "EqnT.hpp"

#include "FemusDefault.hpp"

#include "numeric_vectorM.hpp"
#include "dense_vectorM.hpp"
#include "sparse_matrixM.hpp"
#include "dense_matrixM.hpp"
#include "linear_solverM.hpp"

#include "Files.hpp"
#include "Utils.hpp"
#include "Physics.hpp"
#include "mesh.hpp"
#include "GeomEl.hpp"
#include "EquationsMap.hpp"
#include "EqnBase.hpp"
#include "FEType_enum.hpp"
#include "FEElemBase.hpp"
#include "QRule.hpp"
#include "NormTang_enum.hpp"
#include "VBType_enum.hpp"
#include "QTYnum_enum.hpp"
#include "Domain.hpp"
#include "TimeLoop.hpp"
#include "CurrGaussPoint.hpp"
#include "CurrElem.hpp"
#include "paral.hpp"

// application
#include "Temp_conf.hpp"
#include "TempQuantities.hpp"
#include "TempPhysics.hpp"
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
EqnT::EqnT(  std::vector<Quantity*> int_map_in,
             EquationsMap& equations_map_in,
             std::string eqname_in,
             std::string varname_in):
    EqnBase(int_map_in,equations_map_in,eqname_in,varname_in) {

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
void  EqnT::GenMatRhsVB(const uint vb, const double time,const uint Level) {

  CurrElem       currelem(*this,_eqnmap);
//   CurrGaussPoint   currgp(_eqnmap);    
  CurrGaussPointBase & currgp = CurrGaussPointBase::build(_eqnmap, _mesh._dim);
  
  TempPhysics* myphys; myphys = static_cast<TempPhysics*>(&_phys);

//========== PROCESSOR INDEX
  const uint myproc = _iproc;

//==========FLAG FOR STATIONARITY OR NOT
  const double    dt = _eqnmap._timeloop._timemap.get("dt");
  const uint Nonstat = myphys->_physrtmap.get("NonStatTEMP");
  
//========= BCHandling =========
  const double penalty_val = _utils._urtmap.get("penalty_val");    

  //======== ELEMENT MAPPING =======
  const uint space_dim = _mesh._dim;
  const uint  meshql   = (int) _utils._urtmap.get("meshql");
  const uint  mesh_ord = (int) _utils._urtmap.get("mesh_ord");

  
//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    QuantityLocal Tempold(currgp,currelem);
    Tempold._qtyptr   = _QtyInternalVector[0]; 
    Tempold.VectWithQtyFillBasic();
    Tempold._val_dofs = new double[Tempold._dim*Tempold._ndof[vb]];
    Tempold._val_g    = new double[Tempold._dim];

//====================================
    QuantityLocal Tlift(currgp,currelem);
    Tlift._qtyptr   = _eqnmap._qtymap.get_qty("Qty_TempLift");//_QtyInternalVector[1]; 
    Tlift.VectWithQtyFillBasic();
    Tlift._val_dofs = new double[Tlift._dim*Tlift._ndof[vb]];
    Tlift._val_g    = new double[Tlift._dim];
    Tlift._grad_g   = new double*[Tlift._dim];
    Tlift._grad_g[0]= new double[space_dim];

//=====================================
    QuantityLocal TAdj(currgp,currelem);
    TAdj._qtyptr   = _eqnmap._qtymap.get_qty("Qty_TempAdj");//_QtyInternalVector[2]; 
    TAdj.VectWithQtyFillBasic();
    TAdj._val_dofs = new double[TAdj._dim*TAdj._ndof[vb]];
    TAdj._val_g    = new double[TAdj._dim];
    TAdj._grad_g   = new double*[TAdj._dim];
    TAdj._grad_g[0]= new double[space_dim];    
    
//=========EXTERNAL QUANTITIES (couplings) =====
    //========= //DOMAIN MAPPING
    QuantityLocal xyz(currgp,currelem);  //no quantity
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof[VV] = _AbstractFE[xyz._FEord]->_ndof[VV];
    xyz._ndof[BB] = _AbstractFE[xyz._FEord]->_ndof[BB];
    xyz._val_dofs = new double[xyz._dim*xyz._ndof[vb]];
    xyz._val_g    = new double[xyz._dim];

    //==================Quadratic domain, auxiliary, must be QUADRATIC!!! ==========
  QuantityLocal xyz_refbox(currgp,currelem);  //no quantity
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof[VV] = _mesh._GeomEl._elnds[VV][xyz_refbox._FEord];
  xyz_refbox._ndof[BB] = _mesh._GeomEl._elnds[BB][xyz_refbox._FEord];
  xyz_refbox._val_dofs = new double[xyz_refbox._dim*xyz_refbox._ndof[vb]]; 
  xyz_refbox._val_g    = new double[xyz_refbox._dim];
  xyz_refbox._el_average.resize(VB);
  for (uint i=0; i<VB; i++)  xyz_refbox._el_average[i].resize(xyz_refbox._dim);
  //==================

    QuantityLocal vel(currgp,currelem);
    vel._qtyptr   = _eqnmap._qtymap.get_qty("Qty_Velocity"); 
    vel.VectWithQtyFillBasic();
    vel._val_dofs = new double[vel._dim*vel._ndof[vb]];
    vel._val_g    = new double[vel._dim];   
    
//===============Tdes=====================
    QuantityLocal Tdes(currgp,currelem);
    Tdes._qtyptr   = _eqnmap._qtymap.get_qty("Qty_TempDes"); 
    Tdes.VectWithQtyFillBasic();
    Tdes._val_dofs = new double[Tdes._dim*Tdes._ndof[vb]];
    Tdes._val_g    = new double[Tdes._dim];

  //==== AUXILIARY ==============
    double* dphijdx_g = new double[space_dim];
    double* dphiidx_g = new double[space_dim];
   //-----Nonhomogeneous Neumann------
 // Qflux = - k grad(T) by definition
//  QfluxDOTn>0: energy flows outside (cooling)  QfluxDOTn<0: energy flows inside (heating)
    double* Qflux_g = new double[space_dim];

  //====== reference values ========================
  const double IRe = 1./myphys->_Re;
  const double IPr = 1./myphys->_Pr;
  const double alphaT  = myphys->_physrtmap.get("alphaT");
  const double alphaL2 = myphys->_physrtmap.get("alphaL2");
  const double alphaH1 = myphys->_physrtmap.get("alphaH1");

    int T4_ord = T4_ORD;
  
  /// b) Element  Loop over the volume (n_elem)
   const uint el_ngauss = _eqnmap._qrule._NoGaussVB[vb];
   const uint nel_e = _mesh._off_el[vb][_NoLevels*myproc+Level+1];
   const uint nel_b = _mesh._off_el[vb][_NoLevels*myproc+Level];

  if (vb==VV)   {//BEGIN VOLUME
    
  for (uint iel=0; iel < (nel_e - nel_b); iel++) {
    
    currelem._KeM[vb].zero();
    currelem._FeM[vb].zero(); 

    currelem.get_el_nod_conn_lev_subd(vb,Level,myproc,iel);
    currelem.get_el_DofObj_lev_subd(vb,Level,myproc,iel);
    currelem.get_el_ctr(vb);

    currelem.ConvertElemCoordsToMappingOrd(vb,xyz);
    _mesh.TransformElemNodesToRef(vb,currelem._xx_nds[vb],xyz_refbox._val_dofs);    

    
//MY EQUATION
//the elements are, for every level:
// 1)DOF INDICES
// 2)BC FLAGS
// 3)BC VALUES 
// 1) and 2) are taken in a single vector, 3) are considered separately
      
    currelem.GetElDofsBc(vb,Level);

  Tempold.GetElDofsVect(vb,Level);
    Tlift.GetElDofsVect(vb,Level);
     TAdj.GetElDofsVect(vb,Level);
     
    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(vb,Tempold._FEord,currelem._bc_eldofs[vb]); //only the Qtyzero Part is modified!

// ===============      
// Now the point is this: there are several functions of space
// which are expressed with respect to a reference frame
//ok, now that the dofs are filled for xyz_refbox, I can use the el_average
//Well, the alternative is to consider  the elem in the refbox as
    //either a Vect or a CurrElem !
    //I could consider it as another element, but only with the geometrical part!

  xyz_refbox.SetElemAverage(vb);
  
int domain_flag = myphys->ElFlagControl(xyz_refbox._el_average[vb]);
//====================    
    
//===== FILL the DOFS of the EXTERNAL QUANTITIES: you must assure that for every Vect the quantity is set correctly
// for every Vect it must be clear if it belongs to an equation or not, and which equation it belongs to;
// this is usually made clear by the related QUANTITY.
// Now the main thing to check is the difference between Vect WITH QUANTITY and Vect WITHOUT QUANTITY.
  // it is better to avoid using GetElDofs if the Vect is internal, only if external  
  //Do not use GetElDofs if you want to pick an intermediate dof...
      
   if ( vel._eqnptr != NULL )  vel.GetElDofsVect(vb,Level);
   else                        vel._qtyptr->FunctionDof(vb,vel,time,xyz_refbox._val_dofs);

   if ( Tdes._eqnptr != NULL )  Tdes.GetElDofsVect(vb,Level);
   else                         Tdes._qtyptr->FunctionDof(vb,Tdes,time,xyz_refbox._val_dofs);


    for (uint qp=0; qp< el_ngauss; qp++) {

//======= "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)   { 
  currgp.SetPhiElDofsFEVB_g (vb,fe,qp);
  currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp); 
}
	  
const double      det = dt*currgp.JacVectVV_g(vb,xyz);
const double dtxJxW_g = det*_eqnmap._qrule._weightVB[vb][qp];
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { currgp.SetDPhiDxyzElDofsFEVB_g   (vb,fe,qp); }
for (uint fe = 0; fe < QL; fe++)     { currgp.ExtendDphiDxyzElDofsFEVB_g(vb,fe); }
//======= end of the "COMMON SHAPE PART"==================

 	Tempold.val_g(vb); 
          Tlift.val_g(vb); 
           TAdj.val_g(vb); 
            vel.val_g(vb); 
           Tdes.val_g(vb);
	   
	   // always remember to get the dofs for the variables you use!
           // The point is that you fill the dofs with different functions...
           // you should need a flag to check if the dofs have been correctly filled

      /// d) Local (element) assemblying energy equation
      for (uint i=0; i < Tempold._ndof[vb]; i++)     {

        const double phii_g = currgp._phi_ndsQLVB_g[vb][Tempold._FEord][i];

        for (uint idim = 0; idim < space_dim; idim++) dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[vb][Tempold._FEord][i+idim*Tempold._ndof[vb]];

//=========== FIRST ROW ===============
        currelem._FeM[vb](i) +=      
           currelem._bc_eldofs[vb][i]*dtxJxW_g*( 
                Nonstat*Tempold._val_g[0]*phii_g/dt
	  )
	   + (1-currelem._bc_eldofs[vb][i])*detb*(Tempold._val_dofs[i]);
        
        currelem._KeM[vb](i,i) +=  (1-currelem._bc_eldofs[vb][i])*detb;

//========= SECOND ROW (CONTROL) =====================
	 int ip1 = i + /* 1* */Tempold._ndof[vb];   //suppose that T' T_0 T_adj have the same order
	 currelem._FeM[vb](ip1) +=      
           currelem._bc_eldofs[vb][ip1]*dtxJxW_g*( 
                Nonstat*Tempold._val_g[0]*phii_g/dt
                     + alphaT*domain_flag*(Tdes._val_g[0])*phii_g // T_d delta T_0    /////// ADDED /////
	  )
	   + (1-currelem._bc_eldofs[vb][ip1])*detb*(Tlift._val_dofs[i]);
        
         currelem._KeM[vb](ip1,ip1) +=  (1-currelem._bc_eldofs[vb][ip1])*detb;

//======= THIRD ROW (ADJOINT) ===================================
	 int ip2 = i + 2 * Tempold._ndof[vb];   //suppose that T' T_0 T_adj have the same order
           currelem._FeM[vb](ip2) +=      
           currelem._bc_eldofs[vb][ip2]*dtxJxW_g*( 
                Nonstat*Tempold._val_g[0]*phii_g/dt
                + alphaT*domain_flag*(Tdes._val_g[0])*phii_g // T_d delta T'
	     )
	   + (1-currelem._bc_eldofs[vb][ip2])*detb*(Tempold._val_dofs[i]);
        
        currelem._KeM[vb](ip2,ip2) +=  (1-currelem._bc_eldofs[vb][ip2])*detb;

#if FOURTH_ROW==1
	 int ip3 = i + 3 * Tempold._ndof[vb];   //suppose that T' T_0 T_adj have the same order
	 
	 if (i < _AbstractFE[T4_ord]->_ndof[vb]) { currelem._FeM[vb](ip3) +=  currelem._bc_eldofs[vb][ip3]*dtxJxW_g*(currgp._phi_ndsQLVB_g[vb][T4_ord][i]) + (1-currelem._bc_eldofs[vb][ip3])*detb*1300.;
	              currelem._KeM[vb](ip3,ip3)  += ( 1-currelem._bc_eldofs[vb][ip3] )*detb;  }
#endif
	 // Matrix Assemblying ---------------------------
        for (uint j=0; j<Tempold._ndof[vb]; j++) {
          double phij_g = currgp._phi_ndsQLVB_g[vb][Tempold._FEord][j];
	  
        for (uint idim = 0; idim < space_dim; idim++)  dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[vb][Tempold._FEord][j+idim*Tempold._ndof[vb]]; 
           
   
          double Lap_g   = _utils.dot(dphijdx_g,dphiidx_g,space_dim);
          double Advection = _utils.dot(vel._val_g,dphijdx_g,space_dim);

	    int ip1 = i + Tempold._ndof[vb];
	    int jp1 = j + Tempold._ndof[vb];
	    int ip2 = i + 2*Tempold._ndof[vb];
	    int jp2 = j + 2*Tempold._ndof[vb];

// 	           T     T_0     T_adj
	    
// 	    T      X      X       O
	     
// 	    T_0   
	    
// 	    T_adj
	    
	    
//============ FIRST ROW state  delta T ===============
//======= DIAGONAL =============================
	   currelem._KeM[vb](i,j) +=        
            currelem._bc_eldofs[vb][i]*dtxJxW_g*( 
              Nonstat*phij_g*phii_g/dt 
            + Advection*phii_g
            + IRe*IPr*Lap_g  
            );

//===============================
    //same operators for T and T_0
	    currelem._KeM[vb](i,jp1) +=        
            currelem._bc_eldofs[vb][i]*dtxJxW_g*(    
              Nonstat*phij_g*phii_g/dt 
            + Advection*phii_g
            + IRe*IPr*Lap_g
	    );

//====================================
	   currelem._KeM[vb](i,jp2) +=        
            currelem._bc_eldofs[vb][i]*dtxJxW_g*( 
                0.
            );

	    
//============= SECOND ROW (LIFTING) delta T_0 =============
//===== DIAGONAL ===========================
         currelem._KeM[vb](ip1,jp1) +=        
            currelem._bc_eldofs[vb][ip1]*
            dtxJxW_g*( 
              Nonstat*phij_g*phii_g/dt
             + alphaL2*phij_g*phii_g  //L_2 control norm
             + alphaH1*Lap_g          //H_1 control norm
              + alphaT*domain_flag*(phij_g)*phii_g  //T_0 delta T_0  //ADDED///////////////
            ); 
//====================================
	   currelem._KeM[vb](ip1,j) +=        
            currelem._bc_eldofs[vb][ip1]*
            dtxJxW_g*( 
                + alphaT*domain_flag*(phij_g)*phii_g  //T' delta T_0     //ADDED///////////////
            );
//====================================
	   currelem._KeM[vb](ip1,jp2) +=        
            currelem._bc_eldofs[vb][ip1]*
             dtxJxW_g*( 
                 -Advection*phii_g
                + IRe*IPr*Lap_g
           );

//============= THIRD ROW (ADJOINT) =============
//======= DIAGONAL ==================
          currelem._KeM[vb](ip2,jp2) +=        
            currelem._bc_eldofs[vb][ip2]*
              dtxJxW_g*( 
              Nonstat*phij_g*phii_g/dt
            - Advection*phii_g  //minus sign
            + IRe*IPr*Lap_g
               
            ); 
//====================================
	   currelem._KeM[vb](ip2,j) +=        
            currelem._bc_eldofs[vb][ip2]*
            dtxJxW_g*( 
               + alphaT*domain_flag*(phij_g)*phii_g  //T' delta T'
            );
//====================================
	   currelem._KeM[vb](ip2,jp1) +=        
            currelem._bc_eldofs[vb][ip2]*
            dtxJxW_g*( 
               + alphaT*domain_flag*(phij_g)*phii_g  //T_0 delta T'     ///ADDED///////
            );

#if FOURTH_ROW==1
	 int ip3 = i + 3*Tempold._ndof[vb];   //suppose that T' T_0 T_adj have the same order
// 	    int jp3 = j + 3*Tempold._ndof[vb];
	   if (i < _AbstractFE[T4_ord]->_ndof[vb]) currelem._KeM[vb](ip3,ip3) += currelem._bc_eldofs[vb][ip3]*dtxJxW_g*(currgp._phi_ndsQLVB_g[vb][ T4_ord ][/*j*/i]*currgp._phi_ndsQLVB_g[vb][ T4_ord ][i]);   
#endif
	    
        }  //end j (col)
      }   //end i (row)
    } // end of the quadrature point qp-loop

       _A[Level]->add_matrix(currelem._KeM[vb],currelem._el_dof_indices[vb]);
       _b[Level]->add_vector(currelem._FeM[vb],currelem._el_dof_indices[vb]);
  } // end of element loop
  // *****************************************************************

  }//END VOLUME
  
  else if (vb==BB)  {//BEGIN BOUNDARY  // *****************************************************************
    
     for (uint iel=0;iel < (nel_e - nel_b) ; iel++) {

      currelem._KeM[vb].zero();
      currelem._FeM[vb].zero();

      currelem.get_el_nod_conn_lev_subd(vb,Level,myproc,iel);
      currelem.get_el_DofObj_lev_subd(vb,Level,myproc,iel);
      currelem.get_el_ctr(vb); 

      currelem.ConvertElemCoordsToMappingOrd(vb,xyz);
    _mesh.TransformElemNodesToRef(vb,currelem._xx_nds[vb],xyz_refbox._val_dofs);    
     
      currelem.GetElDofsBc(vb,Level);
      
       Tempold.GetElDofsVect(vb,Level);
         Tlift.GetElDofsVect(vb,Level);
          TAdj.GetElDofsVect(vb,Level);

     if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(vb,Tempold._FEord,currelem._bc_eldofs[vb]); //only the Quadratic Part is modified!
  
 //============ FLAGS ================
     double el_penalty = 0.;
     int pen_sum=0;
     for (uint i=0; i< Tempold._ndof[vb]; i++)   pen_sum += currelem._bc_eldofs[vb][i];
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
     for (uint i=0; i < Tempold._ndof[vb]; i++)   Neum_sum += currelem._bc_eldofs[vb][i];
     for (uint i=0; i < Tempold._ndof[vb]; i++)   Neum_sum += currelem._bc_eldofs[vb][i + Tempold._ndof[vb]];
            if ( Neum_sum == 2*Tempold._ndof[vb] )  { el_Neum_flag=1;  }

//====================================

    for (uint qp=0; qp< el_ngauss; qp++) {

//======= "COMMON SHAPE PART"============================
  for (uint fe = 0; fe < QL; fe++)  {
    currgp.SetPhiElDofsFEVB_g (vb,fe,qp);
    currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp); 
  }
        const double  det   = dt*currgp.JacVectBB_g(vb,xyz);
        const double dtxJxW_g = det * _eqnmap._qrule._weightVB[vb][qp];
//=======end "COMMON SHAPE PART"===================================

       xyz.val_g(vb);
       
       static_cast<Temperature*>(_eqnmap._qtymap.get_qty("Qty_Temperature"))->heatflux_txyz(time,xyz._val_g,Qflux_g);
   
	Tempold.val_g(vb); //For the penalty Dirichlet //i need this for interpolating the old function at the gauss point

	   double QfluxDn_g=_utils.dot( Qflux_g,currgp.get_normal_ptr(),space_dim );
	 
        for (uint i=0; i<Tempold._ndof[vb]; i++) {
   	const double phii_g =  currgp._phi_ndsQLVB_g[vb][Tempold._FEord][i]; 
	
       currelem._FeM[vb](i) +=
          currelem._bc_eldofs[vb][i]*
         el_Neum_flag*dtxJxW_g*(-QfluxDn_g)*phii_g    // beware of the sign  //this integral goes in the first equation
	 + el_penalty*dtxJxW_g*Tempold._val_g[0]*phii_g  //clearly, if you continue using bc=0 for setting nodal Dirichlet, this must go outside
	 ; 
	 
         if (_Dir_pen_fl == 1) {
            for (uint j=0; j<Tempold._ndof[vb]; j++) {
               double phij_g = currgp._phi_ndsQLVB_g[vb][Tempold._FEord][j];
	       currelem._KeM[vb](i,j) += el_penalty*dtxJxW_g*phij_g*phii_g;
	    } 
          }

	  //end of j loop
	}
          // end of i loop
    }
        // end BDRYelement gaussian integration loop
        
        _A[Level]->add_matrix(currelem._KeM[vb],currelem._el_dof_indices[vb]);
        _b[Level]->add_vector(currelem._FeM[vb],currelem._el_dof_indices[vb]);
   
  }
      // end of BDRYelement loop
    
  }//END BOUNDARY

else {   std::cout << " No line integrals yet... " << std::endl; abort();}

 
//=========== cleaning stage ==============
  delete []  Tdes._val_g;    delete []  Tdes._val_dofs;
  delete []  TAdj._val_g;    delete []  TAdj._val_dofs;
  delete []  Tlift._val_g;    delete []  Tlift._val_dofs;
  delete []  Tempold._val_g;    delete []  Tempold._val_dofs;
  delete []  xyz._val_g;        delete []  xyz._val_dofs;
  delete []  vel._val_g;        delete []  vel._val_dofs;
  
#ifdef DEFAULT_PRINT_INFO
  std::cout << " Matrix and RHS assembled for equation " << _eqname
            << " Level "<< Level << " dofs " << _A[Level]->n() << " vb = " << vb << std::endl;
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




//===================================
// The best place for this function is here,
//because we have all the structures for INTEGRATION;
//all in all, an equation means putting integrals
//in the entries of a matrix

//in this routine you dont need 
//DOF maps nor
//BC maps

// You only need DOF VALUES

// This function computes the integral only for the current processor

double EqnT::ComputeIntegral (const uint vb, const uint Level) {

    CurrElem       currelem(*this,_eqnmap);  //TODO in these functions you only need the GEOMETRIC PART, not the DOFS PART
    CurrGaussPointBase & currgp = CurrGaussPointBase::build(_eqnmap, _mesh._dim);

  //====== Physics cast
  TempPhysics *optphys; optphys = static_cast<TempPhysics*>(&_phys);

  //====== processor index
  const uint myproc = _iproc;
  const uint           space_dim = _mesh._dim;
  const uint mesh_ord = (int) _utils._urtmap.get("mesh_ord");  
  const uint meshql   = (int) _utils._urtmap.get("meshql");   //======== ELEMENT MAPPING =======
 
  //========== 
    QuantityLocal Tempold(currgp,currelem);
    Tempold._qtyptr   =  _eqnmap._qtymap.get_qty("Qty_Temperature"); 
    Tempold.VectWithQtyFillBasic();
    Tempold._val_dofs = new double[Tempold._dim*Tempold._ndof[vb]];
    Tempold._val_g    = new double[Tempold._dim];

  //========== 
    QuantityLocal Tlift(currgp,currelem);
    Tlift._qtyptr   =  _eqnmap._qtymap.get_qty("Qty_TempLift"); 
    Tlift.VectWithQtyFillBasic();
    Tlift._val_dofs = new double[Tlift._dim*Tlift._ndof[vb]];
    Tlift._val_g    = new double[Tlift._dim];
    
 //===========
    QuantityLocal Tdes(currgp,currelem);
         Tdes._qtyptr   = _eqnmap._qtymap.get_qty("Qty_TempDes"); 
    Tdes.VectWithQtyFillBasic();
    Tdes._val_dofs = new double[Tdes._dim*Tdes._ndof[vb]];
    Tdes._val_g    = new double[Tdes._dim];
  
//========= DOMAIN MAPPING
    QuantityLocal xyz(currgp,currelem);
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof[VV] = _AbstractFE[xyz._FEord]->_ndof[VV];
    xyz._ndof[BB] = _AbstractFE[xyz._FEord]->_ndof[BB];
    xyz._val_dofs = new double[xyz._dim*xyz._ndof[vb]];
    xyz._val_g    = new double[xyz._dim];

//========== Quadratic domain, auxiliary  
  QuantityLocal xyz_refbox(currgp,currelem);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof[VV] = _mesh._GeomEl._elnds[VV][xyz_refbox._FEord];
  xyz_refbox._ndof[BB] = _mesh._GeomEl._elnds[BB][xyz_refbox._FEord];
  xyz_refbox._val_dofs = new double[xyz_refbox._dim*xyz_refbox._ndof[vb]]; 
  xyz_refbox._val_g    = new double[xyz_refbox._dim];
  xyz_refbox._el_average.resize(VB);
  for (uint i=0; i<VB; i++)  xyz_refbox._el_average[i].resize(xyz_refbox._dim);
  
    
   double integral = 0.;

      const uint el_ngauss = _eqnmap._qrule._NoGaussVB[vb];

//parallel sum
    const uint nel_e = _mesh._off_el[vb][_NoLevels*myproc+Level+1];
    const uint nel_b = _mesh._off_el[vb][_NoLevels*myproc+Level];
  
    for (uint iel=0; iel < (nel_e - nel_b); iel++) {

      currelem.get_el_nod_conn_lev_subd(vb,Level,myproc,iel);
      currelem.get_el_DofObj_lev_subd(vb,Level,myproc,iel);
      currelem.get_el_ctr(vb);
      
      currelem.ConvertElemCoordsToMappingOrd(vb,xyz);
      _mesh.TransformElemNodesToRef(vb,currelem._xx_nds[vb],xyz_refbox._val_dofs);

// =============== 
      xyz_refbox.SetElemAverage(vb);
      int el_flagdom = optphys->ElFlagControl(xyz_refbox._el_average[vb]);
//====================     
 
    if ( Tempold._eqnptr != NULL )   Tempold.GetElDofsVect(vb,Level);
    else                             Tempold._qtyptr->FunctionDof(vb,Tempold,0.,xyz_refbox._val_dofs);
    if ( Tlift._eqnptr != NULL )       Tlift.GetElDofsVect(vb,Level);
    else                               Tlift._qtyptr->FunctionDof(vb,Tlift,0.,xyz_refbox._val_dofs);
    if ( Tdes._eqnptr != NULL )         Tdes.GetElDofsVect(vb,Level);
    else                                Tdes._qtyptr->FunctionDof(vb,Tdes,0.,xyz_refbox._val_dofs);    



    for (uint qp = 0; qp < el_ngauss; qp++) {

     for (uint fe = 0; fe < QL; fe++)     {  currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp);  }  
     
   const double  Jac_g = currgp.JacVectVV_g(vb,xyz);  //not xyz_refbox!      
   const double  wgt_g = _eqnmap._qrule._weightVB[vb][qp];

     for (uint fe = 0; fe < QL; fe++)     {          currgp.SetPhiElDofsFEVB_g (vb,fe,qp);  }

 Tempold.val_g(vb);
   Tlift.val_g(vb);
    Tdes.val_g(vb);

  double deltau_squarenorm_g = 0.;
   deltau_squarenorm_g += (Tempold._val_g[0] + Tlift._val_g[0] - Tdes._val_g[0])*
                          (Tempold._val_g[0] + Tlift._val_g[0] - Tdes._val_g[0]); 

  integral += el_flagdom*wgt_g*Jac_g*deltau_squarenorm_g;
   
    }//gauss loop
     
    }//element loop
    
////////////////////////////////////////        
       std::cout << "integral on processor 0: " << integral << std::endl;

   double J=0.;
#ifdef HAVE_MPI
      MPI_Allreduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  //THIS IS THE RIGHT ONE!!
#else
      J = integral;
#endif
   
    std::cout << "@@@@@@@@@@@@@@@@ functional value: " << J << std::endl;
//=====functional evaluation =======   

    
    ///////// let us also print the functional value in a unique file,
    /////////so that we explore the variation wrt alpha
    
    std::string intgr_fname = _eqnmap._utils._files.get_basepath() + "/" +
                              _eqnmap._utils._files.get_frtmap().get("OUTPUT_DIR") + "/" + "alpha";
  
	std::ofstream intgr_fstream;
    
    if (paral::get_rank() ==0 ){ 
      intgr_fstream.open(intgr_fname.c_str(),ios_base::app); 
      intgr_fstream << _eqnmap._utils._files.get_frtmap().get("OUTTIME_DIR") << " " << optphys->_physrtmap.get("alphaT") << " " << optphys->_physrtmap.get("injsuc")<< " "  << J << " " << std::endl ; 
      intgr_fstream.close();  //you have to close to disassociate the file from the stream
}
 
    
  return J;  
    
}


/////////////////////////

double EqnT::ComputeNormControl (const uint vb, const uint Level, const uint reg_ord ) {

  //reg_ord = 0: L2
  //reg_ord = 1: H1

    CurrElem       currelem(*this,_eqnmap);  //TODO in these functions you only need the GEOMETRIC PART, not the DOFS PART
    CurrGaussPointBase & currgp = CurrGaussPointBase::build(_eqnmap, _mesh._dim);
  
  // processor index
  const uint myproc = _iproc;
  const uint space_dim = _mesh._dim;
  const uint mesh_ord = (int) _utils._urtmap.get("mesh_ord");  
  const uint meshql   = (int) _utils._urtmap.get("meshql");    //======== ELEMENT MAPPING =======

  
//======Functions in the integrand ============

      //========== 
    QuantityLocal Tlift(currgp,currelem);
         Tlift._qtyptr   =  _eqnmap._qtymap.get_qty("Qty_TempLift"); 
    Tlift.VectWithQtyFillBasic();
    Tlift._val_dofs = new double[Tlift._dim*Tlift._ndof[vb]];
    Tlift._val_g    = new double[Tlift._dim];
    Tlift._grad_g       = new double*[1];
    Tlift._grad_g[0]    = new double[space_dim]; 

//========= DOMAIN MAPPING
    QuantityLocal xyz(currgp,currelem);
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof[VV] = _AbstractFE[xyz._FEord]->_ndof[VV];
    xyz._ndof[BB] = _AbstractFE[xyz._FEord]->_ndof[BB];
    xyz._val_dofs = new double[xyz._dim*xyz._ndof[vb]];
    xyz._val_g    = new double[xyz._dim];

//========== Quadratic domain, auxiliary  
  QuantityLocal xyz_refbox(currgp,currelem);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof[VV] = _mesh._GeomEl._elnds[VV][xyz_refbox._FEord];
  xyz_refbox._ndof[BB] = _mesh._GeomEl._elnds[BB][xyz_refbox._FEord];
  xyz_refbox._val_dofs = new double[xyz_refbox._dim*xyz_refbox._ndof[vb]]; 
  xyz_refbox._val_g    = new double[xyz_refbox._dim];
  
    
   double integral = 0.;

//loop over the geom el types
      const uint el_ngauss = _eqnmap._qrule._NoGaussVB[vb];   //elem gauss points

//parallel sum
    const uint nel_e = _mesh._off_el[vb][_NoLevels*myproc+Level+1];
    const uint nel_b = _mesh._off_el[vb][_NoLevels*myproc+Level];
  
    for (int iel=0; iel < (nel_e - nel_b); iel++) {

      currelem.get_el_nod_conn_lev_subd(vb,Level,myproc,iel);
      currelem.get_el_DofObj_lev_subd(vb,Level,myproc,iel);
      currelem.get_el_ctr(vb);

      currelem.ConvertElemCoordsToMappingOrd(vb,xyz);
      _mesh.TransformElemNodesToRef(vb,currelem._xx_nds[vb],xyz_refbox._val_dofs);
     
     Tlift.GetElDofsVect(vb,Level);


  for (uint qp = 0; qp < el_ngauss; qp++) {

     for (uint fe = 0; fe < QL; fe++)     {
       currgp.SetPhiElDofsFEVB_g (vb,fe,qp);
       currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp);  
    }  
     
      const double  Jac_g = currgp.JacVectVV_g(vb,xyz);  //not xyz_refbox!      
      const double  wgt_g = _eqnmap._qrule._weightVB[vb][qp];

  Tlift.val_g(vb);
  Tlift.grad_g(vb);

  double deltau_squarenorm_g = 0.;
   deltau_squarenorm_g += (1-reg_ord) * (Tlift._val_g[0])*(Tlift._val_g[0]); 

   for (uint idim = 0; idim < space_dim; idim++)  {   deltau_squarenorm_g  += reg_ord * (Tlift._grad_g[0][idim])*(Tlift._grad_g[0][idim])  ;   }

  integral += /*el_flagdom**/wgt_g*Jac_g*deltau_squarenorm_g;   //Do it over ALL THE DOMAIN!
   
    }//gauss loop
     
    }//element loop
    
       std::cout << "integral on processor 0: " << integral << std::endl;

   double J = 0.;
#ifdef HAVE_MPI
      MPI_Allreduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  //THIS IS THE RIGHT ONE!!
#else
      J = integral;
#endif
   
    std::cout << "@@@@@@@@@@@@@@@@ control norm: " << reg_ord << " " << J << std::endl;
    
  return J;  
  
}



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

