#include "EqnMHDCONT.hpp"

//library includes
#include "FemusDefault.hpp"

#include "dense_matrixM.hpp"
#include "sparse_matrixM.hpp"
#include "DenseVector.hpp"
#include "numeric_vectorM.hpp"
#include "linear_solverM.hpp"

#include "Utils.hpp"
#include "Physics.hpp"
#include "mesh.hpp"
#include "GeomEl.hpp"
#include "EquationsMap.hpp"
#include "FEType_enum.hpp"
#include "FEElemBase.hpp"
#include "QRule.hpp"
#include "NormTang_enum.hpp"
#include "Quantity.hpp"
#include "QTYnum_enum.hpp"
#include "TimeLoop.hpp"
#include "CurrGaussPoint.hpp"
#include "CurrElem.hpp"


// application
#include "Opt_conf.hpp"
#include "OptPhysics.hpp"
#include "EqnNS.hpp"
#include "EqnMHD.hpp"
#include "EqnNSAD.hpp"
#include "EqnMHDAD.hpp"



//now, the idea is this:
//if you have specific data for this equation, then you may 
//define also specific member functions for it
  EqnMHDCONT::EqnMHDCONT( std::vector<Quantity*> int_map_in,
	           EquationsMap& equations_map_in,
                   std::string eqname_in,
                   std::string varname_in):
      EqnBase(int_map_in,equations_map_in,eqname_in,varname_in)
  {

//====== VARNAMES of the equation=================
   _var_names[0]="Becontx";
   _var_names[1]="Beconty";
#if (DIMENSION==3)
    _var_names[2]="Becontz";
#endif
    _var_names[DIMENSION]="Becontp";

    
//=========  REFVALUES of the Unknown quantities of the equation===============
   _refvalue[0] = _QtyInternalVector[0]->_refvalue[0];
   _refvalue[1] = _QtyInternalVector[0]->_refvalue[1];
#if (DIMENSION==3)
    _refvalue[2]= _QtyInternalVector[0]->_refvalue[2];
#endif
   _refvalue[DIMENSION]= _QtyInternalVector[1]->_refvalue[0];

//=========== Solver ================
     for(uint l=0;l<_NoLevels;l++)    _solver[l]->set_solver_type(SOLVERMHDCONT);

//====== Dirichlet penalty     ==================
    _Dir_pen_fl = MHDCONT_DIR_PENALTY;
    
    
  }

//=======================
//destructor
  EqnMHDCONT::~EqnMHDCONT()  {
    
    for (uint Level =0; Level<_NoLevels; Level++) {
          delete _x_oldopt[Level];
      }
    
      _x_oldopt.clear();

    
  }

void EqnMHDCONT::init_equation_data() {
  
//======equation-specific vectors     =====================
      _x_oldopt.resize(_NoLevels);
      
        for(uint  Level=0; Level<_NoLevels; Level++)  {
   uint n_glob=_Dim[Level]; //is it already filled? Now yes!!!!!!!!!
  _x_oldopt[Level] = NumericVectorM::build().release(); _x_oldopt[Level]->init(n_glob,false, SERIALM);
       }
 
  
 return; 
}
  
  
/// This function assembles the matrix and the rhs:

  void EqnMHDCONT::GenMatRhsVB(const uint vb,const double time,const uint Level)  {

    CurrElem       currelem(*this,_eqnmap);
    CurrGaussPointBase & currgp = CurrGaussPointBase::build(_eqnmap, _mesh._dim);
    
//======= TIME - STATIONARY OR NOT =======
  const int NonStatMHDCONT = (int) _phys._physrtmap.get("NonStatMHDCONT");
  const double        dt   = _eqnmap._timeloop._timemap.get("dt");

//======== GEOMETRICAL ELEMENT =======
  const uint space_dim = _mesh._dim;
  const uint mesh_ord = (int) _utils._urtmap.get("mesh_ord");
  const uint meshql = (int) _utils._urtmap.get("meshql"); //======== ELEMENT MAPPING =======

  //========= BCHandling =========
  const double penalty_val = _utils._urtmap.get("penalty_val");    

//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
  //QTYZERO
    QuantityLocal BeOld(currgp,currelem);
    BeOld._qtyptr   = _QtyInternalVector[QTYZERO];
    BeOld.VectWithQtyFillBasic();
    BeOld._val_dofs = new double[BeOld._dim*BeOld._ndof[vb]];
    BeOld._val_g    = new double[BeOld._dim];

  //QTYONE
    QuantityLocal LagMultOld(currgp,currelem);
    LagMultOld._qtyptr   = _QtyInternalVector[QTYONE];
    LagMultOld.VectWithQtyFillBasic();
    LagMultOld._val_dofs = new double[LagMultOld._dim*LagMultOld._ndof[vb]];
    LagMultOld._val_g    = new double[LagMultOld._dim];
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    QuantityLocal xyz(currgp,currelem);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof[VV] = _AbstractFE[xyz._FEord]->_ndof[VV];
    xyz._ndof[BB] = _AbstractFE[xyz._FEord]->_ndof[BB];
    xyz._val_dofs = new double[xyz._dim*xyz._ndof[vb]];
    xyz._val_g    = new double[xyz._dim];

    //==================Quadratic domain, auxiliary  //this is something that belongs to the user, so we must put it somewhere else,in the MeshUser or something
  QuantityLocal xyz_refbox(currgp,currelem);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof[VV] = _mesh._GeomEl._elnds[VV][xyz_refbox._FEord];
  xyz_refbox._ndof[BB] = _mesh._GeomEl._elnds[BB][xyz_refbox._FEord];
  xyz_refbox._val_dofs = new double[xyz_refbox._dim*xyz_refbox._ndof[vb]]; 
  xyz_refbox._val_g    = new double[xyz_refbox._dim];
  
#if VELOCITY_QTY==1
    QuantityLocal Vel(currgp,currelem);
    Vel._qtyptr      = _eqnmap._qtymap.get_qty("Qty_Velocity"); //an alternative cannot exist, because it is an Unknown of This Equation
    Vel.VectWithQtyFillBasic();
    Vel._val_dofs    = new double[Vel._dim*Vel._ndof[vb]]; 
    Vel._val_dofs3D  = new double[       3*Vel._ndof[vb]]; //needed for a vector product
    Vel._val_g       = new double[Vel._dim];
    Vel._val_g3D     = new double[3];   //needed for a vector product
#endif  

    //==================
    QuantityLocal VelAdj(currgp,currelem);
    VelAdj._qtyptr      = _eqnmap._qtymap.get_qty("Qty_VelocityAdj");
    VelAdj.VectWithQtyFillBasic();
    VelAdj._val_dofs    = new double[VelAdj._dim*VelAdj._ndof[vb]];
    VelAdj._val_dofs3D  = new double[       3*VelAdj._ndof[vb]]; //needed for a vector product
    VelAdj._val_g       = new double[VelAdj._dim];
    VelAdj._val_g3D     = new double[3];   //needed for a vector product 
  
    //==================
    QuantityLocal Bhom(currgp,currelem); 
    Bhom._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldHom");
    Bhom.VectWithQtyFillBasic();
    Bhom._val_dofs =   new double[Bhom._dim*Bhom._ndof[vb]];
    Bhom._val_dofs3D = new double[        3*Bhom._ndof[vb]];  //this is for the curl, not for val_dofs3D!!! 
    Bhom._val_g      = new double[Bhom._dim];
    Bhom._curl_g3D   = new double[3];

//===============
    QuantityLocal BhomAdj(currgp,currelem); 
    BhomAdj._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldHomAdj"); 
    BhomAdj.VectWithQtyFillBasic();
    BhomAdj._val_dofs =   new double[BhomAdj._dim*BhomAdj._ndof[vb]];
    BhomAdj._val_dofs3D = new double[           3*BhomAdj._ndof[vb]];
    BhomAdj._curl_g3D =  new double[3];
    BhomAdj._grad_g     = new double*[BhomAdj._dim];
    for (uint i=0; i< BhomAdj._dim;i++) {BhomAdj._grad_g[i] = new double[DIMENSION];}
  //=========END EXTERNAL QUANTITIES (couplings) =====

    //========= reference values =========
  //====== Physics
  OptPhysics *optphys; optphys = static_cast<OptPhysics*>(&_phys);
  double       IRem = 1./optphys->_Rem;
  double          S = optphys->_S;
  const double betaL2   = _phys._physrtmap.get("betaL2");
  const double gammaLap = _phys._physrtmap.get("gammaLap");
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

  
   const uint el_ngauss = _eqnmap._qrule._NoGaussVB[vb];   //elem gauss points
    
    const uint nel_e = _mesh._off_el[vb][_NoLevels*_iproc+Level+1];
    const uint nel_b = _mesh._off_el[vb][_NoLevels*_iproc+Level];

  if (vb==VV)   {//BEGIN VOLUME


  for (uint iel=0; iel < (nel_e - nel_b); iel++) {
 
    currelem._KeM[vb].zero();
    currelem._FeM[vb].zero(); 
     
    currelem.get_el_nod_conn_lev_subd(vb,Level,_iproc,iel);
    currelem.get_el_DofObj_lev_subd(vb,Level,_iproc,iel);
    currelem.get_el_ctr(vb);
    
    currelem.ConvertElemCoordsToMappingOrd(vb,xyz);
    _mesh.TransformElemNodesToRef(vb,currelem._xx_nds[vb],xyz_refbox._val_dofs); 

    currelem.GetElDofsBc(vb,Level);
    
         BeOld.GetElDofsVect(vb,Level);  
    LagMultOld.GetElDofsVect(vb,Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(vb,BeOld._FEord,currelem._bc_eldofs[vb]);  //only the Quadratic Part is modified!
 
    
    if ( Vel._eqnptr != NULL )        Vel.GetElDofsVect(vb,Level);
    else                              Vel._qtyptr->FunctionDof(vb,Vel,time,xyz_refbox._val_dofs);
    if ( VelAdj._eqnptr != NULL )  VelAdj.GetElDofsVect(vb,Level);
    else                           VelAdj._qtyptr->FunctionDof(vb,VelAdj,time,xyz_refbox._val_dofs);
    if ( Bhom._eqnptr != NULL )      Bhom.GetElDofsVect(vb,Level);
    else                             Bhom._qtyptr->FunctionDof(vb,Bhom,time,xyz_refbox._val_dofs);
    if ( BhomAdj._eqnptr != NULL ) BhomAdj.GetElDofsVect(vb,Level);
    else                           BhomAdj._qtyptr->FunctionDof(vb,BhomAdj,time,xyz_refbox._val_dofs);    
    

//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
    for (uint qp = 0; qp < el_ngauss; qp++) {

//======= here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     {          currgp.SetPhiElDofsFEVB_g (vb,fe,qp);  }
for (uint fe = 0; fe < QL; fe++)     {  currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp);  }

 const double      det = dt*currgp.JacVectVV_g(vb,xyz);   //InvJac: is unique!
 const double dtxJxW_g = det*_eqnmap._qrule._weightVB[vb][qp];
 const double     detb = det/el_ngauss;

for (uint fe = 0; fe < QL; fe++)     {    currgp.SetDPhiDxyzElDofsFEVB_g (vb,fe,qp); }
for (uint fe = 0; fe < QL; fe++)     { currgp.ExtendDphiDxyzElDofsFEVB_g (vb,fe); }
//======= end of the "COMMON SHAPE PART"==================

//========preparation for things that are independent of (i,j), dofs of test and shape =====================
       BhomAdj.curl_g(vb);
          Bhom.curl_g(vb);
       
                               //Bhom.curl(); but how can i automatically activate the "new" for curlVect? I can do them all 
                               //if i do a class I can do them ALL together in the constructor,
                               //and delete them automatically in the destructor.

       BeOld.val_g(vb); 
         Vel.val_g(vb);
      VelAdj.val_g(vb);
        Bhom.val_g(vb);
     BhomAdj.grad_g(vb);

// vector product
      _utils.extend(   Vel._val_g,   Vel._val_g3D,space_dim);
      _utils.extend(VelAdj._val_g,VelAdj._val_g3D,space_dim);

    _utils.cross(BhomAdj._curl_g3D,   Vel._val_g3D,curlxiXvel_g3D ); 
    _utils.cross(   Bhom._curl_g3D,VelAdj._val_g3D,curlbXlambda_g3D );
//========end preparation for things that are independent of (i,j) dofs of test and shape =====================
  
//================================
//========= FILLING ELEMENT MAT/RHS
//=================================
//============ QTYZERO dofs (rows) ============
//=================================
//actually this is used for both QTYZERO and QTYONE,
//because the dofs of pressure are FEWER...
       for (uint i = 0; i < BeOld._ndof[vb]; i++)     {

//============ preparation for (i) ============
//(i) is used only for the RHS or for both MATRIX and RHS
//first, you put whatever is related to the Test Functions of THESE rows:
//test function value, test function first derivative
        const double                             phii_g       =      currgp._phi_ndsQLVB_g[vb][BeOld._FEord][i];
        for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[vb][BeOld._FEord][i+idim*BeOld._ndof[vb]];

	 _utils.extend(dphiidx_g,dphiidx_g3D,space_dim);
	 _utils.cross(BhomAdj._curl_g3D,dphiidx_g3D,curlxiXdphii_g3D);

	  double bDdphii_g      = _utils.dot(  Bhom._val_g,dphiidx_g,space_dim);
	  double lambdaDdphii_g = _utils.dot(VelAdj._val_g,dphiidx_g,space_dim);

	  for (uint idim=0; idim<space_dim; idim++) Lapxi_g[idim]=0.;
	  for (uint idim=0; idim<space_dim; idim++)  {
	      for (uint jdim=0; jdim<space_dim; jdim++) {
                    Lapxi_g[idim] += BhomAdj._grad_g[idim][jdim]*dphiidx_g[jdim];  //TODO CHECK THIS
                          }
	  }
//============end preparation for (i) ============

//============ QTYZERO rhs ============
         for (uint idim=0; idim<space_dim; idim++) {
            const uint irowq = i+idim*BeOld._ndof[vb];
            
            currelem._FeM[vb](irowq) += currelem._bc_eldofs[vb][irowq]*dtxJxW_g*(
                          NonStatMHDCONT*BeOld._val_g[idim]*phii_g/dt
                           -(1-LAP_MHD)*IRem*curlxiXdphii_g3D[idim]                 //from MHD  1
                               -LAP_MHD*IRem*Lapxi_g[idim]                              //from MHD  2
                           + curlxiXvel_g3D[idim]*phii_g                            //from MHD
                           - S*curlbXlambda_g3D[idim]*phii_g                             //from NS
                           + S*(bDdphii_g*VelAdj._val_g[idim] - lambdaDdphii_g*Bhom._val_g[idim])     //from NS
                         )
                         +(1-currelem._bc_eldofs[vb][irowq])*detb*BeOld._val_dofs[irowq]; //Dirichlet bc
	   }
//============ END QTYZERO rhs ============

        for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
          const uint irowq = i+idim*BeOld._ndof[vb];
          currelem._KeM[vb](irowq,irowq) += (1-currelem._bc_eldofs[vb][irowq])*detb;
        }
                                           // end filling diagonal for Dirichlet bc

//============ QTYZERO x QTYZERO dofs matrix (A matrix) ============
        for (uint j=0; j<BeOld._ndof[vb]; j++) {
//============preparation for (j) ============
//here you put things that depend either on (j) or on (i,j)
          double                                  phij_g       =      currgp._phi_ndsQLVB_g[vb][BeOld._FEord][j];
          for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[vb][BeOld._FEord][j+idim*BeOld._ndof[vb]];

	  double Lap_g=_utils.dot(dphijdx_g,dphiidx_g,space_dim);
	  double lambdaDdphii_g=_utils.dot(VelAdj._val_g,dphiidx_g,space_dim);  //TODO why did i recompute it?
	  double lambdaDdphij_g=_utils.dot(VelAdj._val_g,dphijdx_g,space_dim);
//============ end preparation for (j) ============

          for (uint idim=0; idim<space_dim; idim++) {
            int irowq = i+idim*BeOld._ndof[vb];
            // diagonal blocks [1-5-9]
            currelem._KeM[vb](irowq,j+idim*BeOld._ndof[vb])
            += currelem._bc_eldofs[vb][irowq]*dtxJxW_g*(
                 NonStatMHDCONT*phij_g*phii_g/dt// time
                 + betaL2*phij_g*phii_g
                 + gammaLap*Lap_g
                 + S*phij_g*(lambdaDdphii_g  - VelAdj._val_g[idim]*dphiidx_g[idim] )
                 + S*phii_g*(lambdaDdphij_g  - VelAdj._val_g[idim]*dphijdx_g[idim] )
               );
            // block +1 [2-6-7]
            int idimp1=(idim+1)%space_dim;
            currelem._KeM[vb](irowq,j+idimp1*BeOld._ndof[vb])
            += currelem._bc_eldofs[vb][irowq]*dtxJxW_g*(
                 + S*phij_g*( -VelAdj._val_g[idim]*dphiidx_g[idimp1])
                 + S*phii_g*( -VelAdj._val_g[idimp1]*dphijdx_g[idim])
               );
#if (DIMENSION==3)
            // block +2 [3-4-8]
            int idimp2=(idim+2)%space_dim;
            currelem._KeM[vb](irowq,j+idimp2*BeOld._ndof[vb])
            += currelem._bc_eldofs[vb][irowq]*dtxJxW_g*(
                  + S*phij_g*( -VelAdj._val_g[idim]*dphiidx_g[idimp2]) 
                  + S*phii_g*( -VelAdj._val_g[idimp2]*dphijdx_g[idim])
               );
#endif
          }
                    
        } 
//============ END QTYZERO x QTYZERO dofs matrix (A matrix) ============

//============ QTYZERO x QTYONE dofs matrix (B^T matrix) // ( p*div(v) ) (NS eq) ============
       for (uint j=0; j < LagMultOld._ndof[vb]; j++) {
          const double psij_g = currgp._phi_ndsQLVB_g[vb][LagMultOld._FEord][j];
          const int jclml = j+space_dim*BeOld._ndof[vb];
          for (uint idim=0; idim<space_dim; idim++) {
            uint irowq=i+idim*BeOld._ndof[vb];
            currelem._KeM[vb](irowq,jclml) += currelem._bc_eldofs[vb][irowq]*dtxJxW_g*(-psij_g*dphiidx_g[idim]);
           }
        }
//============ END QTYZERO x QTYONE dofs matrix (B^T matrix) ============

//===================================================================
//============ END QTYZERO dofs (rows) ============
//===================================================================

//============================================
//============ QTYONE dofs (rows) ============
//===============================================
          if ( i < LagMultOld._ndof[vb] ) {
//============ preparation for (i) ============
          double psii_g = currgp._phi_ndsQLVB_g[vb][LagMultOld._FEord][i];
	  const uint irowl = i+space_dim*BeOld._ndof[vb];

//============ QTYONE rhs ============
          currelem._FeM[vb](irowl)=0.;
//============ END QTYONE rhs ============

//============ QTYONE x QTYONE ============
 //             currelem._KeM[vb](irowl,j+space_dim*el_ndof_q)  += dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt; // (KOMP dp/dt=rho*div) 
//============ END QTYONE x QTYONE ============

//============ QTYONE x QTYZERO (B matrix) q*div(u) ============
          for (uint j = 0; j < BeOld._ndof[vb]; j++) { // B element matrix 
//============ preparation for (j) ============
            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim]= currgp._dphidxyz_ndsQLVB_g[vb][BeOld._FEord][j+idim*BeOld._ndof[vb]];
            for (uint idim=0; idim<space_dim; idim++) currelem._KeM[vb](irowl,j+idim*BeOld._ndof[vb]) += -dtxJxW_g*psii_g*dphijdx_g[idim]; 
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
    _A[Level]->add_matrix(currelem._KeM[vb],currelem._el_dof_indices[vb]);
    _b[Level]->add_vector(currelem._FeM[vb],currelem._el_dof_indices[vb]);
    
  } 
  // end of element loop

  }//END VOLUME
  
    // *****************************************************************
    // *****************************************************************
    // *****************************************************************

    else if (vb==BB)  {//BEGIN BOUNDARY  // *****************************************************************
  

  for (uint iel=0;iel < (nel_e - nel_b) ; iel++) {

     currelem._KeM[vb].zero();
     currelem._FeM[vb].zero();

     currelem.get_el_nod_conn_lev_subd(vb,Level,_iproc,iel);
     currelem.get_el_DofObj_lev_subd(vb,Level,_iproc,iel);
     currelem.get_el_ctr(vb);
     
     currelem.ConvertElemCoordsToMappingOrd(vb,xyz);
    _mesh.TransformElemNodesToRef(vb,currelem._xx_nds[vb],xyz_refbox._val_dofs);    

     currelem.GetElDofsBc(vb,Level);
     
          BeOld.GetElDofsVect(vb,Level);
     LagMultOld.GetElDofsVect(vb,Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(vb,BeOld._FEord,currelem._bc_eldofs[vb]); //only the Quadratic Part is modified! /*OK DIR_PEN*/


//============ BC =======
       int     el_flag[NT] = {0,0}; //normal and tangential flag
#if DIMENSION==2
       double el_value[N1T1] = {0.,0.};  //1 normal and 1 tangential
#elif DIMENSION==3
       double el_value[N1T3] = {0.,0.,0.,0.}; //1 normal and 3 tangential
#endif
       double  dbl_pen[NT] = {0.,0.};  //normal and tangential penalty value
   
    
       Bc_GetElFlagValLevSubd(Level,_iproc,iel,el_flag,el_value);

if (_Dir_pen_fl == 1)  { 
       if (el_flag[NN] == 1) {   dbl_pen[NN]=penalty_val; } //normal dirichlet
       if (el_flag[TT] == 1) {   dbl_pen[TT]=penalty_val; } //tangential dirichlet
   }

//checks
//TODO here i should check that the nodal bc dirichlet i put correspond to the element NT flags

       uint press_fl=0;
       Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(vb,currelem._bc_eldofs[vb],BeOld,LagMultOld,press_fl); //compute the PRESSURE FLAG with the PRESSURE nodal bc flags
 //only the LINEAR PART is USED!!
       
//   if ( (1-el_flag[NN]) != press_fl)  {std::cout << "Sthg wrong with press elflags" << std::endl;abort();}

//========END BC============
//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
    for (uint qp=0; qp< el_ngauss; qp++) {

//======= "COMMON SHAPE PART"============================
for (uint fe = 0; fe < QL; fe++)     {        currgp.SetPhiElDofsFEVB_g (vb,fe,qp);  } //for velocity test functions AND for pressure shape functions
for (uint fe = 0; fe < QL; fe++)     {      currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp);   }

        const double det   = dt*currgp.JacVectBB_g(vb,xyz);
	const double dtxJxW_g = det*_eqnmap._qrule._weightVB[vb][qp];
//=======end "COMMON SHAPE PART"===================================

   xyz_refbox.val_g(vb);
      LagMultOld._qtyptr->Function_txyz(time,xyz_refbox._val_g/*xyz._val_g*/,LagMultOld._val_g);  //i prefer using the function instead of the p_old vector

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
//============ QTYZERO dofs (rows) ============
      for (uint i=0; i<BeOld._ndof[vb]; i++) {
//============ preparation for (i) ============
const double phii_g = currgp._phi_ndsQLVB_g[vb][BeOld._FEord][i];

//============ QTYZERO rhs ============
               for (uint idim=0; idim< space_dim; idim++)    {
             uint irowq=i+idim*BeOld._ndof[vb];
            currelem._FeM[vb](irowq)  += 
          currelem._bc_eldofs[vb][irowq]*           
           dtxJxW_g*(   -1.*/*press_fl*/(1-el_flag[NN])*LagMultOld._val_g[0]*currgp.get_normal_ptr()[idim]*phii_g  //  //OLD VALUES //AAA multiplying int times uint!!!

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

	   for (uint j=0; j<BeOld._ndof[vb]; j++) {
          const double phij_g = currgp._phi_ndsQLVB_g[vb][BeOld._FEord][j];

  currelem._KeM[vb](irowq,j+jdim*BeOld._ndof[vb]) +=                //projection over the physical (x,y,z) 
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
   
    _A[Level]->add_matrix(currelem._KeM[vb],currelem._el_dof_indices[vb]);//      std::cout << "KeM "<< vb << " " << currelem._KeM[vb].l1_norm() << std::endl;
    _b[Level]->add_vector(currelem._FeM[vb],currelem._el_dof_indices[vb]);//      std::cout << "FeM "<< vb << " " << currelem._FeM[vb].l2_norm() << std::endl;

  }
  //end bdry element loop
      
      
    }  
    
// END BOUNDARY  // **************************


  //TODO DO ALL THE DELETE***************
 //******************************DESTROY ALL THE Vect****************************  
// ===============domain ========
    delete [] xyz_refbox._val_g;
    delete [] xyz_refbox._val_dofs;
    
    delete [] xyz._val_g;
    delete [] xyz._val_dofs;
//================================= 
    
//=========Internal Quantities: no ifdef for them =========
   delete [] BeOld._val_g;
   delete [] BeOld._val_dofs;

   delete [] LagMultOld._val_g;
   delete [] LagMultOld._val_dofs;
//==============================================
  #if VELOCITY_QTY==1
        delete [] Vel._val_g;
        delete [] Vel._val_g3D;
	delete [] Vel._val_dofs;
	delete [] Vel._val_dofs3D;

#endif
  
        delete [] VelAdj._val_g;
        delete [] VelAdj._val_g3D;
	delete [] VelAdj._val_dofs;
	delete [] VelAdj._val_dofs3D;
  
// // //   for (uint k=0;k<(nvars_q);k++)    delete [] dxidx_g[k];
// // //   delete []dxidx_g;
// // //      
   
   
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << _eqname << ": assembled  Level " << Level
           << " with " << _A[Level]->m() << " dofs " << std::endl;
#endif     

  return;
}
