#include "EqnNS.hpp"

#include "FemusDefault.hpp"

#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "linear_solverM.hpp"

#include "Utils.hpp"
#include "Physics.hpp"
#include "mesh.hpp"
#include "GeomEl.hpp"
#include "EquationsMap.hpp"
#include "FETypeEnum.hpp"
#include "FEElemBase.hpp"
#include "QRule.hpp"
#include "NormTangEnum.hpp"
#include "VBType_enum.hpp"
#include "QTYnumEnum.hpp"
#include "Domain.hpp"
#include "TimeLoop.hpp"
#include "CurrGaussPoint.hpp"
#include "CurrElem.hpp"
  

//application
#include "Temp_conf.hpp"
#include "TempQuantities.hpp"
#include "TempPhysics.hpp"
#include "EqnT.hpp"



///=============== Constructor
  EqnNS::EqnNS(    std::vector<Quantity*> int_map_in,  //no reference!
	           EquationsMap& equations_map_in,
                   std::string eqname_in,
                   std::string varname_in):
           EqnBase(int_map_in,equations_map_in,eqname_in,varname_in),
     _AdvPic_fl(ADVPIC_NS),
     _AdvNew_fl(ADVNEW_NS),
     _Stab_fl(STAB_NS),
     _Komp_fac(KOMP_NS)   {


//=======  _var_names[]  ===========
                          _var_names[0] = "ux"; //variable names
                          _var_names[1] = "uy";
   if ( _mesh._dim == 3 ) _var_names[2] = "uz";
                 _var_names[_mesh._dim] = "up";
     
//=======  _refvalue[] ==============   
                          _refvalue[0] =  _QtyInternalVector[0]->_refvalue[0]; 
                          _refvalue[1] =  _QtyInternalVector[0]->_refvalue[1]; 
   if ( _mesh._dim == 3 ) _refvalue[2] =  _QtyInternalVector[0]->_refvalue[2]; 
                 _refvalue[_mesh._dim] =  _QtyInternalVector[1]->_refvalue[0];

//========= MG solver ===================
  for(uint l=0;l<_NoLevels;l++)  _solver[l]->set_solver_type(SOLVERNS);

//============= DIR PENALTY===============
   _Dir_pen_fl = NS_DIR_PENALTY;
   
    }
//====== END CONSTRUCTOR    
    
    
//================ DESTRUCTOR    
      EqnNS::~EqnNS() {    }
    
    
//===================================================
/// This function assembles the matrix and the rhs:
 void EqnNS::GenMatRhsVB(const uint vb,const double time,const uint Level)  {  ///< VB Assemblying.

    CurrElem       currelem(*this,_eqnmap);
//     CurrGaussPoint   currgp(_eqnmap);
//here, I have to put the BASE, because that is my RUN-TIME DIMENSION-INDEPENDENT GaussPoint interface!!!!    
    CurrGaussPointBase & currgp = CurrGaussPointBase::build(_eqnmap, _mesh._dim);


//========== PROCESSOR INDEX
  const uint myproc = _iproc;

//==========FLAG FOR STATIONARITY OR NOT
  const int NonStatNS = (int) _phys._physrtmap.get("NonStatNS");
  const double     dt = _eqnmap._timeloop._timemap.get("dt");

//========= BCHandling =========
  const double penalty_val = _utils._urtmap.get("penalty_val");    


  //========== GEOMETRIC ELEMENT ========
  const uint           space_dim = _mesh._dim /*_IntDim[vb]*/;
  //AAAAAAA! // is this the SPACE dimension, or the UNKNOWN dimension, or the INTEGRATION (MANIFOLD) dimension?
  //well, it depends. The fact is that some operators are defined only for certain dimensions of the vector
  //for instance, divergence acts on vectors whose dimension is like the domain dimension
  const uint mesh_ord = (int) _utils._urtmap.get("mesh_ord");
  const uint meshql   = (int) _utils._urtmap.get("meshql"); //======== ELEMENT MAPPING =======


  
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    QuantityLocal VelOld(currgp,currelem);
    VelOld._qtyptr   = _QtyInternalVector[QTYZERO]; //an alternative cannot exist, because it is an Unknown of This Equation
    VelOld.VectWithQtyFillBasic();
    VelOld._val_dofs = new double[VelOld._dim*VelOld._ndof[vb]];
    VelOld._val_g    = new double[VelOld._dim];
    VelOld._grad_g   = new double*[VelOld._dim];
  for (uint i=0; i< VelOld._dim;i++) {VelOld._grad_g[i] = new double[space_dim];}

   const uint   qtyzero_ord  = VelOld._FEord;
   const uint   qtyzero_ndof = VelOld._ndof[vb]; 
//     Velocity*  vel_castqtyptr = static_cast<Velocity*>(VelOld._qtyptr); //casting for quantity-specific functions

//=========
    QuantityLocal pressOld(currgp,currelem);
    pressOld._qtyptr   = _QtyInternalVector[QTYONE];
    pressOld.VectWithQtyFillBasic();
    pressOld._val_dofs = new double[pressOld._dim*pressOld._ndof[vb]];
    pressOld._val_g    = new double[pressOld._dim];

   const uint qtyone_ord  = pressOld._FEord;
   const uint qtyone_ndof = pressOld._ndof[vb]; 

   //order
   const uint  qtyZeroToOne_DofOffset = VelOld._ndof[vb]*VelOld._dim;
   
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    QuantityLocal xyz(currgp,currelem); //domain
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof[VV] = _AbstractFE[xyz._FEord]->_ndof[VV];
    xyz._ndof[BB] = _AbstractFE[xyz._FEord]->_ndof[BB];
    xyz._val_dofs = new double[xyz._dim*xyz._ndof[vb]];
    xyz._val_g    = new double[xyz._dim];
//=========
    //==================Quadratic domain, auxiliary, must be QUADRATIC!!! ==========
  QuantityLocal xyz_refbox(currgp,currelem);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof[VV] = _mesh._GeomEl._elnds[VV][xyz_refbox._FEord];
  xyz_refbox._ndof[BB] = _mesh._GeomEl._elnds[BB][xyz_refbox._FEord];
  xyz_refbox._val_dofs = new double[xyz_refbox._dim*xyz_refbox._ndof[vb]]; 
  xyz_refbox._val_g    = new double[xyz_refbox._dim];
  //==================
    
//======================== MAG WORLD ================================
//======================== MAG WORLD ================================

//===================TEMPERATURE WORLD=============================
#if TEMP_QTY==1
    QuantityLocal Temp(currgp,currelem);
    Temp._qtyptr   =  _eqnmap._qtymap.get_qty("Qty_Temperature");
    Temp.VectWithQtyFillBasic();
    Temp._val_dofs = new double[Temp._dim*Temp._ndof[vb]];
    Temp._val_g    = new double[Temp._dim];
#endif
//=================== TEMPERATURE WORLD============================

//other Physical constant Quantities

//=======gravity==================================
  QuantityLocal gravity(currgp,currelem);
  gravity._dim = space_dim;
  gravity._val_g    = new double[gravity._dim];
  gravity._val_g[0] = _phys._physrtmap.get("dirgx");
  gravity._val_g[1] = _phys._physrtmap.get("dirgy");
  if ( space_dim == 3 )   gravity._val_g[2] = _phys._physrtmap.get("dirgz"); 

  //=======density and viscosity===================
  const double rhof = _phys._physrtmap.get("rho0");
  const double  muf = _phys._physrtmap.get("mu0");
  const double rho1=0.; const double mu1=0.;
  //====== reference values ========================
//====== related to Quantities on which Operators act, and to the choice of the "LEADING" EQUATION Operator
  TempPhysics *optphys; optphys = static_cast<TempPhysics*>(&_phys);
  const double IRe = 1./optphys->_Re;
  const double IFr = 1./optphys->_Fr;
//================================================  

//================================================  
//=======Operators @ gauss =======================
  double*    dphijdx_g = new double[space_dim]; // ShapeDer(): used for Laplacian,Divergence, ..., associated to an Unknown Quantity
  double*    dphiidx_g = new double[space_dim]; // Test(): used for Laplacian,Advection,Divergence... associated to an Unknown Quantity
  double*     AdvRhs_g = new double[space_dim]; //Operator: Adv(u,u,phi)
//================================================  
  
 
  /// b) Element  Loop over the volume (n_elem = number of elements of each FEM type)

   const uint el_ngauss = _eqnmap._qrule._NoGaussVB[vb];        // element set up for that FEM type
    
    const uint nel_e = _mesh._off_el[vb][_NoLevels*myproc+Level+1];
    const uint nel_b = _mesh._off_el[vb][_NoLevels*myproc+Level];

  if (vb==VV)   {//BEGIN VOLUME

    for (int iel=0; iel < (nel_e - nel_b); iel++) {
 
    currelem._KeM[vb].zero();
    currelem._FeM[vb].zero(); 

    currelem.get_el_nod_conn_lev_subd(vb,Level,myproc,iel);
    currelem.get_el_DofObj_lev_subd(vb,Level,myproc,iel);
    currelem.get_el_ctr(vb);

    currelem.ConvertElemCoordsToMappingOrd(vb,xyz);
    _mesh.TransformElemNodesToRef(vb,currelem._xx_nds[vb],xyz_refbox._val_dofs);    

//=======RETRIEVE the DOFS of the UNKNOWN QUANTITIES,i.e. MY EQUATION
    currelem.GetElDofsBc(vb,Level);
      VelOld.GetElDofsVect(vb,Level);
    pressOld.GetElDofsVect(vb,Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(vb,qtyzero_ord,currelem._bc_eldofs[vb]); //only the Qtyzero Part is modified!

   
//=======RETRIEVE the DOFS of the COUPLED QUANTITIES    
#if (TEMP_QTY==1)
   if ( Temp._eqnptr != NULL )  Temp.GetElDofsVect(vb,Level);
     else                       Temp._qtyptr->FunctionDof(vb,Temp,time,xyz_refbox._val_dofs);
#endif

//======== TWO PHASE WORLD
    double phase=1.;
    double rho_nd =  phase+(1-phase)*rho1/rhof;
    double  mu_nd =  phase+(1-phase)*mu1/muf;
//======== TWO PHASE WORLD

    
//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
    for (uint qp = 0; qp < el_ngauss; qp++) {  

//=======here starts the "COMMON SHAPE PART"==================
// the call of these things should be related to the Operator
//also, the choice of the QuadratureRule should be dependent of the Involved Operators
//and the FE orders on which these Operators act
//after the COMMON SHAPE PART i dont want to have any dependence on qp anymore!
//these  phi and dphi are used for the different stages:
//BEFORE i: by the interpolation functions
//INSIDE i: for the tEST functions and derivatives
//INSIDE j: for the SHAPE functions and derivatives
//again, it should be the Operator to decide what functions to be called,
//for what FE ORDER and what DERIVATIVE ORDER
//if we decide that the PREPARATION of the tEST and SHAPE 
//of a certain Unknown are COMMON TO ALL,
//Then we must only concentrate on preparing the OTHER involved quantities in that Operator
for (uint fe = 0; fe < QL; fe++)     {          currgp.SetPhiElDofsFEVB_g (vb,fe,qp);  }
for (uint fe = 0; fe < QL; fe++)     {  currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp);  }  
	  
const double      det = dt*currgp.JacVectVV_g(vb,xyz);   //InvJac: is the same for both QQ and LL!
const double dtxJxW_g = det*_eqnmap._qrule._weightVB[vb][qp];
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { currgp.SetDPhiDxyzElDofsFEVB_g   (vb,fe,qp); }
for (uint fe = 0; fe < QL; fe++)     { currgp.ExtendDphiDxyzElDofsFEVB_g(vb,fe); }
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
      VelOld.val_g(vb);     //fills _val_g, needs _val_dofs
      VelOld.grad_g(vb);     //fills _grad_g, needs _val_dofs //TODO can we see the analogies between JacVectVV_g and grad_g?
//Advection all VelOld
      for (uint idim=0; idim<space_dim; idim++) { AdvRhs_g[idim]=0.;}
          for (uint idim=0; idim<space_dim; idim++)  {
	    for (uint b=0; b<space_dim; b++) {
	    AdvRhs_g[idim] += VelOld._val_g[b]*VelOld._grad_g[idim][b]; }   }   //TODO NONLIN // grad is [ivar][idim], i.e. [u v w][x y z]
//Divergence VelOld
	  double Div_g=0.;
          for (uint idim=0; idim<space_dim; idim++)    Div_g += VelOld._grad_g[idim][idim];     //TODO NONLIN

#if (TEMP_QTY==1)
     Temp.val_g(vb);
 #endif

       
#ifdef AXISYMX
      dtxJxW_g  *=yyg;  //what is the symmetry axis? I guess the x axis.
#endif

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
// TODO according to the order we should switch DIM loop and DOF loop
    for (uint i=0; i<qtyzero_ndof; i++)     {
//======="COMMON tEST PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
        //in this part we follow the Test Fncs of the FIRST QUANTITY of THIS Equation
        //since the Equation in this FE code is in WEAK FORM,
        //every Operator would have the test FUNCTIONS as ONE and ONLY ONE OF THEIR ARGUMENTS 
                                    const double phii_g       =      currgp._phi_ndsQLVB_g[vb][qtyzero_ord][i];
        for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[vb][qtyzero_ord][i+idim*qtyzero_ndof];
//======= END "COMMON tEST PART for QTYZERO" ==========
	
	  for (uint idim=0; idim<space_dim; idim++) {
            const uint irowq=i+idim*qtyzero_ndof;  //  (i):       dof of the tEST function
                                                  //(idim): component of the tEST function
           currelem._FeM[vb](irowq) += 
         currelem._bc_eldofs[vb][irowq]*
           dtxJxW_g*(          NonStatNS*rho_nd*     VelOld._val_g[idim]*phii_g/dt  //time
                            + _AdvNew_fl*rho_nd*          AdvRhs_g[idim]*phii_g     //TODO NONLIN
                            +            rho_nd*IFr*gravity._val_g[idim]*phii_g     // gravity                           
                               )
            + (1-currelem._bc_eldofs[vb][irowq])*detb*VelOld._val_dofs[irowq] //Dirichlet bc    
	;
          }
           // end filling element rhs u

if (_Dir_pen_fl == 0)  { //faster than multiplying by _Dir_pen_fl
//actually this is not needed because you add only zeros
//it is only needed to avoid doing the operations on Ke
// Bc_ConvertToDirichletPenalty puts all ONE on the quadratic dofs
  for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
          const uint irowq = i+idim*qtyzero_ndof;
          currelem._KeM[vb](irowq,irowq) += (1-currelem._bc_eldofs[vb][irowq])*detb;
        }
}                                         // end filling diagonal for Dirichlet bc

//============ QTYZERO x QTYZERO dofs matrix (A matrix) ============
        for (uint j=0; j< qtyzero_ndof; j++) {
//======="COMMON SHAPE PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
//   (j):       dof of the SHAPE function
           double                                  phij_g       =      currgp._phi_ndsQLVB_g[vb][qtyzero_ord][j];
           for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[vb][qtyzero_ord][j+idim*qtyzero_ndof];
//======= END "COMMON SHAPE PART for QTYZERO" ==========
  
          double Lap_g=_utils.dot(dphijdx_g,dphiidx_g,space_dim);
	  double Adv_g=_utils.dot(VelOld._val_g,dphijdx_g,space_dim);
          
          for (uint idim=0; idim<space_dim; idim++) { //filled in as 1-2-3 // 4-5-6 // 7-8-9
            int irowq = i+idim*qtyzero_ndof;      //(i) is still the dof of the tEST functions
                                                  //(idim): component of the tEST functions

            currelem._KeM[vb](irowq,j+idim*qtyzero_ndof)  // diagonal blocks [1-5-9] [idim(rows),idim(columns)]  //(idim): component of the SHAPE functions
               += 
            currelem._bc_eldofs[vb][irowq]*    
            dtxJxW_g*(
                   +NonStatNS*                           rho_nd*phij_g*phii_g/dt
                 + _AdvPic_fl*                           rho_nd* Adv_g*phii_g                //TODO NONLIN
                 + _AdvNew_fl*rho_nd*phij_g*VelOld._grad_g[idim][idim]*phii_g                //TODO NONLIN
                 + _AdvPic_fl*_Stab_fl*rho_nd*        0.5*Div_g*phij_g*phii_g                //TODO NONLIN
                 +                         mu_nd*IRe*(      dphijdx_g[idim]*dphiidx_g[idim] + Lap_g)
               );

            int idimp1=(idim+1)%space_dim;    // block +1 [2-6-7] [idim(rows),idim+1(columns)]  //(idimp1): component of the SHAPE functions
            currelem._KeM[vb](irowq,j+idimp1*qtyzero_ndof)
               +=
            currelem._bc_eldofs[vb][irowq]*
            dtxJxW_g*(
                   _AdvNew_fl*rho_nd*phij_g*VelOld._grad_g[idim][idimp1]*phii_g           //TODO NONLIN
                              +            mu_nd*IRe*(     dphijdx_g[idim]*dphiidx_g[idimp1])
               );
#if (DIMENSION==3)
	    
            int idimp2=(idim+2)%space_dim;   // block +2 [3-4-8] [idim(rows),idim+2(columns)]  //(idimp2): component of the SHAPE functions
            currelem._KeM[vb](irowq,j+idimp2*qtyzero_ndof)
               +=
           currelem._bc_eldofs[vb][irowq]*           
            dtxJxW_g*(
                   _AdvNew_fl*rho_nd*phij_g*VelOld._grad_g[idim][idimp2]*phii_g          //TODO NONLIN
                              +            mu_nd*IRe*(     dphijdx_g[idim]*dphiidx_g[idimp2])
               );
#endif
          }
 
        } 
//============ END QTYZERO x QTYZERO dofs matrix (A matrix) ============

//============ QTYZERO x QTYONE dofs matrix (B^T matrix) // ( p*div(v) ) (NS eq) ============
       for (uint j=0; j<qtyone_ndof; j++) {
//======="COMMON SHAPE PART for QTYONE" ==================
	 const double psij_g = currgp._phi_ndsQLVB_g[vb][qtyone_ord][j];
//======="COMMON SHAPE PART for QTYONE" - END ============

          const int jclml = j + qtyZeroToOne_DofOffset;
          for (uint idim=0; idim<space_dim; idim++) {
            uint irowq = i+idim*qtyzero_ndof;
            currelem._KeM[vb](irowq,jclml) +=
               currelem._bc_eldofs[vb][irowq]*                    
               dtxJxW_g*(-psij_g*dphiidx_g[idim]);   /**   (-1.)*/
	    
           } 

        }
//============ END QTYZERO x QTYONE dofs matrix (B^T matrix) ============

          if (i < qtyone_ndof) {
//======="COMMON tEST PART for QTYONE" ============
          double psii_g = currgp._phi_ndsQLVB_g[vb][qtyone_ord][i];
//======= "COMMON tEST PART for QTYONE" - END ============
	  const uint irowl = i + qtyZeroToOne_DofOffset;
          currelem._FeM[vb](irowl)=0.;  // rhs
 //             _KeM[vb](irowl,j+space_dim*qtyzero_ndof)  += (1./dt)*dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt;  //no bc here (KOMP dp/dt=rho*div)

          for (uint j=0; j<qtyzero_ndof; j++) { // B element matrix q*div(u)
//======="COMMON SHAPE PART for QTYZERO" ==================
            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[vb][qtyzero_ord][j+idim*qtyzero_ndof];
//======="COMMON SHAPE PART for QTYZERO" - END ============
	    
            for (uint idim=0; idim<space_dim; idim++) currelem._KeM[vb](irowl,j+idim*qtyzero_ndof) += -/*(1./dt)**/dtxJxW_g*psii_g*dphijdx_g[idim]; 
                }

        }
                         // end pressure eq (cont)
      }
//===================================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) =====================
//===================================================================
    }
//==============================================================
//================== END GAUSS LOOP (qp loop) ======================
//==============================================================
    
    ///  Add element matrix and rhs to the global ones.
                   _A[Level]->add_matrix(currelem._KeM[vb],currelem._el_dof_indices[vb]);//     std::cout << "KeM l1"<< vb << " " << _KeM[vb].l1_norm() << std::endl;
                   _b[Level]->add_vector(currelem._FeM[vb],currelem._el_dof_indices[vb]);//     std::cout << "FeM l2"<< vb << " " << _FeM[vb].l2_norm() << std::endl;

    
/////////HERE, AT THE END OF THE ELEMENT, YOU CAN SEARCH FOR MY ROW in THIS ELEMENT

    
  } 
  // end of element loop

  }//END VOLUME
  
    // *****************************************************************
    // *****************************************************************
    // *****************************************************************

    else if (vb==BB)  {//BEGIN BOUNDARY  // *****************************************************************

//=== auxiliary Operators at the boundary
//   double strainU_g['dimension']['dimension'];
//   double strainUtrDn_g['dimension'];


  for (uint iel=0; iel < (nel_e - nel_b) ; iel++) {

     currelem._KeM[vb].zero();  
     currelem._FeM[vb].zero();

     currelem.get_el_nod_conn_lev_subd(vb,Level,myproc,iel);
     currelem.get_el_DofObj_lev_subd(vb,Level,myproc,iel);
     currelem.get_el_ctr(vb);
     
     currelem.ConvertElemCoordsToMappingOrd(vb,xyz);
     _mesh.TransformElemNodesToRef(vb,currelem._xx_nds[vb],xyz_refbox._val_dofs);    

     currelem.GetElDofsBc(vb,Level);
     
     VelOld.GetElDofsVect(vb,Level);
     pressOld.GetElDofsVect(vb,Level);

     if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(vb,qtyzero_ord,currelem._bc_eldofs[vb]); //only the Quadratic Part is modified! /*OK DIR_PEN*/
       

//============ BC =======
       int     el_flag[NT] = {0,0}; //normal and tangential flag
#if DIMENSION==2
       double el_value[N1T1] = {0.,0.};  //1 normal and 1 tangential
#elif DIMENSION==3
       double el_value[N1T3] = {0.,0.,0.,0.}; //1 normal and 3 tangential
#endif
       double  dbl_pen[NT] = {0.,0.};  //normal and tangential penalty value
   
    
       Bc_GetElFlagValLevSubd(Level,myproc,iel,el_flag,el_value);

if (_Dir_pen_fl == 1)  { 
       if (el_flag[NN] == 1) {   dbl_pen[NN]=penalty_val; } //normal dirichlet
       if (el_flag[TT] == 1) {   dbl_pen[TT]=penalty_val; } //tangential dirichlet
   }

//checks
//TODO here i should check that the nodal bc dirichlet i put correspond to the element NT flags

       uint press_fl=0;
       Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(vb,currelem._bc_eldofs[vb],VelOld,pressOld,press_fl); //compute the PRESSURE FLAG with the PRESSURE nodal bc flags
 //only the LINEAR PART is USED!!
       
// // TODO  if ( (1-el_flag[NN]) != press_fl)  {std::cout << "Sthg wrong with press elflags" << std::endl;abort();}

//========END BC============

//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
    for (uint qp=0; qp< el_ngauss; qp++) {
            
//======= "COMMON SHAPE PART"============================
 for (uint fe = 0; fe < QL; fe++)     {        currgp.SetPhiElDofsFEVB_g (vb,fe,qp);  } //for velocity test functions AND for pressure shape functions

for (uint fe = 0; fe < QL; fe++)     {      currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp); }

        const double det   = dt*currgp.JacVectBB_g(vb,xyz);
	const double dtxJxW_g = det*_eqnmap._qrule._weightVB[vb][qp];
//=======end "COMMON SHAPE PART"===================================

//-------- pressure==============
    //"predof" VS "post gauss method"
    //pay attention to this fact: even if we are in the equation to which pressure is associated
    //(pressure is an unknown here!), we might prefer using the FUNCTION instead,
    //in order to get some values that DO NOT CHANGE with the solution
    //so i'm using the FUNCTION of an UNKNOWN, it is not an external quantity  
	  //clearly, this is ALTERNATIVE to the command
	  
// // //    val_g(vb,pressOld);
   
   xyz_refbox.val_g(vb);
      pressOld._qtyptr->Function_txyz(time,xyz_refbox._val_g/*xyz._val_g*/,pressOld._val_g);  //i prefer using the function instead of the p_old vector
       
//--- strain, derivative of velocity ============== 
      
// // // TODO    /*VelOld._qtyptr*/vel_castqtyptr->strain_txyz_box(time,xyz_refbox._val_g/*xyz._val_g*/,strainU_g);
// // //          for (uint idim=0; idim< space_dim; idim++)  strainUtrDn_g[idim]=0.;
// // // 	    for (uint idim=0; idim< space_dim; idim++)  {   //beware that this is space_dim, not IntDim
// // // 	      for (uint jdim=0; jdim< space_dim; jdim++)    {
// // //                 strainUtrDn_g[idim] += strainU_g[jdim][idim]*_normal_g[jdim];
// // // 	        }
// // // 	    }
//=================================================

      //-----old velocity=============
	  VelOld.val_g(vb);  //substitute el_value...


//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
      for (uint i=0; i<qtyzero_ndof; i++) { // rhs (NS eq) //interpolation of the velocity test functions on the boundary

	const double phii_g = currgp._phi_ndsQLVB_g[vb][qtyzero_ord][i];

         for (uint idim=0; idim< space_dim; idim++)    {
             uint irowq=i+idim*qtyzero_ndof;
            currelem._FeM[vb](irowq)  += 
         currelem._bc_eldofs[vb][irowq]*           
           dtxJxW_g*(   -1.*/*press_fl*/(1-el_flag[NN])*pressOld._val_g[0]*currgp.get_normal_ptr()[idim]*phii_g  //  //OLD VALUES //AAA multiplying int times uint!!!

// // //             TODO STRAIN AT THE BOUNDARY            + /*stress_fl*/el_flag[1]*IRe*strainUtrDn_g[idim]*phii_g 

	  )
                                //projection over the physical (x,y,z)
      + _Dir_pen_fl *dtxJxW_g*phii_g*(dbl_pen[NN]*el_value[0]*currgp.get_normal_ptr()[idim] 
                                    + dbl_pen[TT]*el_value[1]*currgp.get_tangent_ptr()[0][idim]  // VelOld._val_g[idim] instead of el_value...
               #if DIMENSION==3
		                    + dbl_pen[TT]*el_value[1]*currgp.get_tangent_ptr()[1][idim]    
               #endif   
          )
	   ;   
	   
//====================
if (_Dir_pen_fl == 1) {  //much faster than multiplying by _Dir_pen_fl=0 , and much better than removing the code with the #ifdef //  #if (NS_DIR_PENALTY==1)  
	   for (uint jdim=0; jdim< space_dim; jdim++)    {

	   for (uint j=0; j<qtyzero_ndof; j++) {
          const double phij_g = currgp._phi_ndsQLVB_g[vb][qtyzero_ord][j];

  currelem._KeM[vb](irowq,j+jdim*qtyzero_ndof) +=                //projection over the physical (x,y,z) 
      + /*_Dir_pen_fl**/dtxJxW_g*phii_g*phij_g*(dbl_pen[NN]*currgp.get_normal_ptr()[jdim]*currgp.get_normal_ptr()[idim]   //the PENALTY is BY ELEMENT, but the (n,t) is BY GAUSS because we cannot compute now a nodal normal
                                              + dbl_pen[TT]*currgp.get_tangent_ptr()[0][jdim]*currgp.get_tangent_ptr()[0][idim]
                 #if DIMENSION==3
                                              + dbl_pen[TT]*currgp.get_tangent_ptr()[1][jdim]*currgp.get_tangent_ptr()[1][idim]
                #endif   
                   );
	         } //end j
      
               } //end jdim
  
             }  //end penalty if
//====================

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
    
    _A[Level]->add_matrix(currelem._KeM[vb],currelem._el_dof_indices[vb]);//      std::cout << "KeM "<< vb << " " << _KeM[vb].l1_norm() << std::endl;
    _b[Level]->add_vector(currelem._FeM[vb],currelem._el_dof_indices[vb]);//      std::cout << "FeM "<< vb << " " << _FeM[vb].l2_norm() << std::endl;

    
  }
  // end of BDRYelement loop

  
    }
  // END BOUNDARY ******************************
  

  //DESTROY ALL THE Vect  
// ===============domain ========
    delete [] xyz_refbox._val_g; delete [] xyz_refbox._val_dofs;
    delete [] xyz._val_g; delete [] xyz._val_dofs;
//=================================  

//=========Internal Quantities: no ifdef for them =========
   for (uint i=0; i< VelOld._dim;i++) {delete [] VelOld._grad_g[i];}
      delete [] VelOld._grad_g;

   delete []   VelOld._val_g; delete [] VelOld._val_dofs;     
   delete [] pressOld._val_g; delete [] pressOld._val_dofs;
//==============================================
   
//======= inner Vect without Quantity ==========
//not an Unknown
//not from External Equation or Function
   delete [] gravity._val_g;
//======= inner Vect without Quantity ==========
   
#if TEMP_QTY==1  //================   
  delete [] Temp._val_g; delete [] Temp._val_dofs;
#endif         //================


#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << _eqname << ": assembled  Level " << Level
           << " with " << _A[Level]->m() << " dofs" << std::endl;
#endif

    return;
}

