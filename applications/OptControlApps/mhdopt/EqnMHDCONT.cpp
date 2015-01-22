#include "EqnMHDCONT.hpp"

//library includes
#include "FemusDefault.hpp"

#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"

#include "Math.hpp"
#include "MeshTwo.hpp"
#include "GeomEl.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "FETypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "Quantity.hpp"
#include "QTYnumEnum.hpp"
#include "TimeLoop.hpp"
#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"


// application
#include "Opt_conf.hpp"
#include "EqnNS.hpp"
#include "EqnMHD.hpp"
#include "EqnNSAD.hpp"
#include "EqnMHDAD.hpp"

namespace femus {

//now, the idea is this:
//if you have specific data for this equation, then you may 
//define also specific member functions for it
  EqnMHDCONT::EqnMHDCONT( std::vector<Quantity*> int_map_in,
	           MultiLevelProblemTwo& equations_map_in,
                   std::string eqname_in,
                   std::string varname_in):
      SystemTwo(int_map_in,equations_map_in,eqname_in,varname_in)
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

//======================
  EqnMHDCONT::~EqnMHDCONT()  { }

  
/// This function assembles the matrix and the rhs:

  void EqnMHDCONT::GenMatRhs(const uint Level)  {

   const double time =  0.; //_eqnmap._timeloop._curr_time;
   
//======= TIME - STATIONARY OR NOT =======
  const int NonStatMHDCONT = (int) _phys.get("NonStatMHDCONT");
  const double        dt   = 1.; //_eqnmap._timeloop._timemap.get("dt");

//======== GEOMETRICAL ELEMENT =======
  const uint space_dim =       _mesh.get_dim();
  const uint mesh_ord  = (int) _mesh.GetRuntimeMap().get("mesh_ord");
  const uint meshql    = (int) _mesh.GetRuntimeMap().get("meshql"); //======== ELEMENT MAPPING =======

  //========= BCHandling =========
  const double penalty_val = _mesh.GetRuntimeMap().get("penalty_val");    
    //========= reference values =========
  
  //====== Physics
  double       IRem = 1./_phys.get("Rem");
  double          S = _phys.get("S");
  const double betaL2   = _phys.get("betaL2");
  const double gammaLap = _phys.get("gammaLap");
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

  {//BEGIN VOLUME

  const uint mesh_vb = VV;
  
    CurrentElem       currelem(VV,this,_mesh,_eqnmap._elem_type);
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
    
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity BeOld(currgp);
    BeOld._qtyptr   = _QtyInternalVector[QTYZERO];
    BeOld.VectWithQtyFillBasic();
    BeOld.Allocate();

    CurrentQuantity LagMultOld(currgp);
    LagMultOld._qtyptr   = _QtyInternalVector[QTYONE];
    LagMultOld.VectWithQtyFillBasic();
    LagMultOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();

    //==================Quadratic domain, auxiliary
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
  xyz_refbox.Allocate();
  
#if VELOCITY_QTY==1
    CurrentQuantity Vel(currgp);
    Vel._qtyptr      = _eqnmap._qtymap.get_qty("Qty_Velocity"); //an alternative cannot exist, because it is an Unknown of This Equation
    Vel.VectWithQtyFillBasic();
    Vel.Allocate();
#endif  

    //==================
    CurrentQuantity VelAdj(currgp);
    VelAdj._qtyptr      = _eqnmap._qtymap.get_qty("Qty_VelocityAdj");
    VelAdj.VectWithQtyFillBasic();
    VelAdj.Allocate();
  
    //==================
    CurrentQuantity Bhom(currgp); 
    Bhom._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldHom");
    Bhom.VectWithQtyFillBasic();
    Bhom.Allocate();
    
//===============
    CurrentQuantity BhomAdj(currgp); 
    BhomAdj._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldHomAdj"); 
    BhomAdj.VectWithQtyFillBasic();
    BhomAdj.Allocate();
  //=========END EXTERNAL QUANTITIES (couplings) =====

    const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level+1];
    const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level];


  for (uint iel=0; iel < (nel_e - nel_b); iel++) {
 
    currelem.Mat().zero();
    currelem.Rhs().zero(); 
     
    currelem.set_el_nod_conn_lev_subd(Level,_iproc,iel);
    currelem.SetMidpoint();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);
    _mesh.TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]); 

    currelem.SetElDofsBc(Level);
    
         BeOld.GetElemDofs(Level);  
    LagMultOld.GetElemDofs(Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),BeOld._FEord,currelem.GetBCDofFlag());  //only the Quadratic Part is modified!
 
    
    if ( Vel._eqnptr != NULL )        Vel.GetElemDofs(Level);
    else                              Vel._qtyptr->FunctionDof(Vel,time,&xyz_refbox._val_dofs[0]);
    if ( VelAdj._eqnptr != NULL )  VelAdj.GetElemDofs(Level);
    else                           VelAdj._qtyptr->FunctionDof(VelAdj,time,&xyz_refbox._val_dofs[0]);
    if ( Bhom._eqnptr != NULL )      Bhom.GetElemDofs(Level);
    else                             Bhom._qtyptr->FunctionDof(Bhom,time,&xyz_refbox._val_dofs[0]);
    if ( BhomAdj._eqnptr != NULL ) BhomAdj.GetElemDofs(Level);
    else                           BhomAdj._qtyptr->FunctionDof(BhomAdj,time,&xyz_refbox._val_dofs[0]);    
    

//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
   const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();
    
    for (uint qp = 0; qp < el_ngauss; qp++) {

//======= here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     {          
  currgp.SetPhiElDofsFEVB_g (fe,qp);
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}

 const double      det = dt*currgp.JacVectVV_g(xyz);   //InvJac: is unique!
 const double dtxJxW_g = det*_eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
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
    _A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    _b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());
    
  } 
  // end of element loop

  
  }//END VOLUME
  
    // *****************************************************************
    // *****************************************************************
    // *****************************************************************

 {//BEGIN BOUNDARY  // *****************************************************************
  
  const uint mesh_vb = BB;
  
    CurrentElem       currelem(BB,this,_mesh,_eqnmap._elem_type);
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
    
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity BeOld(currgp);
    BeOld._qtyptr   = _QtyInternalVector[QTYZERO];
    BeOld.VectWithQtyFillBasic();
    BeOld.Allocate();

    CurrentQuantity LagMultOld(currgp);
    LagMultOld._qtyptr   = _QtyInternalVector[QTYONE];
    LagMultOld.VectWithQtyFillBasic();
    LagMultOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();

    //==================Quadratic domain, auxiliary
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
  xyz_refbox.Allocate();
  
  //=========END EXTERNAL QUANTITIES (couplings) =====

    const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level+1];
    const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level];

    
  for (uint iel=0;iel < (nel_e - nel_b) ; iel++) {

     currelem.Mat().zero();
     currelem.Rhs().zero();

     currelem.set_el_nod_conn_lev_subd(Level,_iproc,iel);
     currelem.SetMidpoint();
     
     currelem.ConvertElemCoordsToMappingOrd(xyz);
    _mesh.TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]);    

     currelem.SetElDofsBc(Level);
     
          BeOld.GetElemDofs(Level);
     LagMultOld.GetElemDofs(Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),BeOld._FEord,currelem.GetBCDofFlag()); //only the Quadratic Part is modified! /*OK DIR_PEN*/


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
       Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(currelem.GetBCDofFlag(),BeOld,LagMultOld,press_fl); //compute the PRESSURE FLAG with the PRESSURE nodal bc flags
 //only the LINEAR PART is USED!!
       
//   if ( (1-el_flag[NN]) != press_fl)  {std::cout << "Sthg wrong with press elflags" << std::endl;abort();}

//========END BC============
//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
   const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();
    
    for (uint qp=0; qp< el_ngauss; qp++) {

//======= "COMMON SHAPE PART"============================
for (uint fe = 0; fe < QL; fe++)     {        
  currgp.SetPhiElDofsFEVB_g (fe,qp);
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);   
}

        const double det   = dt*currgp.JacVectBB_g(xyz);
	const double dtxJxW_g = det*_eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
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

	   for (uint j=0; j<BeOld._ndof; j++) {
          const double phij_g = currgp._phi_ndsQLVB_g[BeOld._FEord][j];

  currelem.Mat()(irowq,j+jdim*BeOld._ndof) +=                //projection over the physical (x,y,z) 
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
   
    _A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    _b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());

  }
  //end bdry element loop
      
  
    }  
    
// END BOUNDARY  // **************************

#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << _eqname << ": assembled  Level " << Level
           << " with " << _A[Level]->m() << " dofs " << std::endl;
#endif     

  return;
}


} //end namespace femus

