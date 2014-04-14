#include "EqnMHDAD.hpp"

#include "FemusDefault.hpp"

#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "LinearSolverM.hpp"

#include "Math.hpp"
#include "Physics.hpp"
#include "Quantity.hpp"
#include "MeshTwo.hpp"
#include "GeomEl.hpp"
#include "EquationsMap.hpp"
#include "FETypeEnum.hpp"
#include "FEElemBase.hpp"
#include "QRule.hpp"
#include "NormTangEnum.hpp"
#include "QTYnumEnum.hpp"
#include "TimeLoop.hpp"

#include "CurrGaussPoint.hpp"
#include "CurrElem.hpp"

#include "Opt_conf.hpp"
#include "EqnNS.hpp"
#include "EqnNSAD.hpp"
#include "EqnMHD.hpp"
#include "EqnMHDCONT.hpp"
#include "OptPhysics.hpp"

namespace femus {

/// Constructor.
  EqnMHDAD::EqnMHDAD(std::vector<Quantity*> int_map_in,
	           EquationsMap& mg_equations_map_in,
                   std::string eqname_in,
                   std::string varname_in):
      EqnBase(int_map_in,mg_equations_map_in,eqname_in,varname_in)      
      {

//=======  _var_names[]  ===========
    _var_names[0]="xix";
    _var_names[1]="xiy";
#if (DIMENSION==3)
    _var_names[2]="xiz";
#endif 
    _var_names[DIMENSION]="xip";
  
//=======  _refvalue[] ==============   
    _refvalue[0] = _QtyInternalVector[0]->_refvalue[0];
    _refvalue[1] = _QtyInternalVector[0]->_refvalue[1];
#if (DIMENSION==3)
    _refvalue[2] = _QtyInternalVector[0]->_refvalue[2];
#endif
    _refvalue[DIMENSION] = _QtyInternalVector[1]->_refvalue[0];

//========= MG solver ===================
   for(uint l=0;l<_NoLevels;l++)   _solver[l]->set_solver_type(SOLVERMHDAD);

//============= DIR PENALTY===============
   _Dir_pen_fl = MHDAD_DIR_PENALTY;   
   
    }


//==============
EqnMHDAD::~EqnMHDAD() {}


/// This function assembles the matrix and the rhs:


  void EqnMHDAD::GenMatRhsVB(const uint vb,const double time,const uint Level)  {

    CurrElem       currelem(*this,_eqnmap);
    CurrGaussPointBase & currgp = CurrGaussPointBase::build(_eqnmap, _mesh._dim);
   
    
//======= TIME - STATIONARY OR NOT =======
const int NonStatMHDAD = (int) _phys._physrtmap.get("NonStatMHDAD");
  const double   dt = _eqnmap._timeloop._timemap.get("dt");

//======== GEOMETRICAL ELEMENT =======
  const uint space_dim =       _mesh._dim;
  const uint  mesh_ord = (int) _mesh._mesh_rtmap.get("mesh_ord");
  const uint    meshql = (int) _mesh._mesh_rtmap.get("meshql");  //======== ELEMENT MAPPING =======

//========= BCHandling =========
  const double penalty_val =   _mesh._mesh_rtmap.get("penalty_val");    

//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
  
     //QTYZERO
    QuantityLocal BhomAdjOld(currgp,currelem);
    BhomAdjOld._qtyptr   = _QtyInternalVector[QTYZERO];
    BhomAdjOld.VectWithQtyFillBasic();
    BhomAdjOld._val_dofs = new double[BhomAdjOld._dim*BhomAdjOld._ndof[vb]];
    BhomAdjOld._val_g    = new double[BhomAdjOld._dim];
  
    //QTYONE
    QuantityLocal BhomLagMultAdjOld(currgp,currelem);
    BhomLagMultAdjOld._qtyptr   = _QtyInternalVector[QTYONE];
    BhomLagMultAdjOld.VectWithQtyFillBasic();
    BhomLagMultAdjOld._val_dofs = new double[BhomLagMultAdjOld._dim*BhomLagMultAdjOld._ndof[vb]];
    BhomLagMultAdjOld._val_g    = new double[BhomLagMultAdjOld._dim];  

//========= END INTERNAL QUANTITIES (unknowns of the equation) ================= 

//=========EXTERNAL QUANTITIES (couplings) =====

 //========= DOMAIN MAPPING
    QuantityLocal xyz(currgp,currelem);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof[VV] = _AbstractFE[xyz._FEord]->_ndof[VV];
    xyz._ndof[BB] = _AbstractFE[xyz._FEord]->_ndof[BB];
    xyz._val_dofs = new double[xyz._dim*xyz._ndof[vb]];
    xyz._val_g    = new double[xyz._dim];

//========== Quadratic domain, auxiliary  
  QuantityLocal xyz_refbox(currgp,currelem);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof[VV] = _mesh._GeomEl._elnds[VV][xyz_refbox._FEord];
  xyz_refbox._ndof[BB] = _mesh._GeomEl._elnds[BB][xyz_refbox._FEord];
  xyz_refbox._val_dofs = new double[xyz_refbox._dim*xyz_refbox._ndof[vb]]; 
  xyz_refbox._val_g    = new double[xyz_refbox._dim];   

  //==========     
    QuantityLocal Vel(currgp,currelem);
    Vel._qtyptr      = _eqnmap._qtymap.get_qty("Qty_Velocity");
    Vel.VectWithQtyFillBasic();
    Vel._val_dofs    = new double[Vel._dim*Vel._ndof[vb]];
    Vel._val_g       = new double[Vel._dim];
 
    //==========    
    QuantityLocal VelAdj(currgp,currelem);
    VelAdj._qtyptr      = _eqnmap._qtymap.get_qty("Qty_VelocityAdj");
    VelAdj.VectWithQtyFillBasic();
    VelAdj._val_dofs    = new double[VelAdj._dim*VelAdj._ndof[vb]];
    VelAdj._val_dofs3D  = new double[       3*VelAdj._ndof[vb]]; //needed for a vector product
    VelAdj._val_g       = new double[VelAdj._dim];
    VelAdj._val_g3D     = new double[3];   //needed for a vector product
   
    //==========    
    QuantityLocal Bhom(currgp,currelem);
    Bhom._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldHom");
    Bhom.VectWithQtyFillBasic();
    Bhom._val_dofs = new double[Bhom._dim*Bhom._ndof[vb]];

//=========
    QuantityLocal Bext(currgp,currelem);
    Bext._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldExt");
    Bext.VectWithQtyFillBasic();
    Bext._val_dofs = new double[Bext._dim*Bext._ndof[vb]];
  
//========= auxiliary, must be AFTER Bhom! //TODO this is an example of  Vect which is not associated to a Quantity
    QuantityLocal Bmag(currgp,currelem); //total
    Bmag._dim        = Bhom._dim;               //same as Bhom
    Bmag._FEord      = Bhom._FEord;             //same as Bhom
    Bmag._ndof[VV]   = _AbstractFE[Bmag._FEord]->_ndof[VV];
    Bmag._ndof[BB]   = _AbstractFE[Bmag._FEord]->_ndof[BB];
    Bmag._val_dofs   = new double[Bmag._dim*Bmag._ndof[vb]]; //we might use their own vector!
    Bmag._val_dofs3D = new double[        3*Bmag._ndof[vb]];
    Bmag._val_g      = new double[Bmag._dim];
    Bmag._curl_g3D   = new double[3];

    
//========= END EXTERNAL QUANTITIES =================
  //====== Physics
  OptPhysics *optphys; optphys = static_cast<OptPhysics*>(&_phys);
  //========= parameters
   double IRem =  1./optphys->_Rem;
   double S    = optphys->_S;
  
  //=========== Operators 

  double dphijdx_g[DIMENSION];
  double dphiidx_g[DIMENSION];
  double  curlBXlambda_g3D[3]; 

    const uint el_ngauss = _eqnmap._qrule._NoGaussVB[vb];
    
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
    
           BhomAdjOld.GetElDofsVect(vb,Level);
    BhomLagMultAdjOld.GetElDofsVect(vb,Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(vb,BhomAdjOld._FEord,currelem._bc_eldofs[vb]);  //only the Quadratic Part is modified!
    
    
     if ( Vel._eqnptr != NULL )      Vel.GetElDofsVect(vb,Level);
    else                             Vel._qtyptr->FunctionDof(vb,Vel,time,xyz_refbox._val_dofs);
    if ( VelAdj._eqnptr != NULL ) VelAdj.GetElDofsVect(vb,Level);
    else                          VelAdj._qtyptr->FunctionDof(vb,VelAdj,time,xyz_refbox._val_dofs);
    if ( Bhom._eqnptr != NULL )     Bhom.GetElDofsVect(vb,Level);
    else                            Bhom._qtyptr->FunctionDof(vb,Bhom,time,xyz_refbox._val_dofs);
    if ( Bext._eqnptr != NULL )     Bext.GetElDofsVect(vb,Level);
    else                            Bext._qtyptr->FunctionDof(vb,Bext,time,xyz_refbox._val_dofs);

//======SUM Bhom and Bext  //from now on, you'll only use Bmag //Bmag,Bext and Bhom must have the same orders!
    Math::zeroN(Bmag._val_dofs,Bmag._dim*Bmag._ndof[vb]);

    for (uint ivarq=0; ivarq < Bmag._dim; ivarq++)    { //ivarq is like idim
          for (uint d=0; d < Bmag._ndof[vb]; d++)    {
          const uint     indxq  =         d + ivarq*Bmag._ndof[vb];
          Bmag._val_dofs[indxq] = Bext._val_dofs[indxq] + Bhom._val_dofs[indxq];
	  }
    }    
//=======


    for (uint qp = 0; qp < el_ngauss; qp++) {
//=======here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     {          currgp.SetPhiElDofsFEVB_g (vb,fe,qp);  }
for (uint fe = 0; fe < QL; fe++)     {  currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp);  }  
	  
const double      det = dt*currgp.JacVectVV_g(vb,xyz);   //InvJac: is the same for both QQ and LL!
const double dtxJxW_g = det * _eqnmap._qrule._weightVB[vb][qp];
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { currgp.SetDPhiDxyzElDofsFEVB_g   (vb,fe,qp); }
for (uint fe = 0; fe < QL; fe++)     { currgp.ExtendDphiDxyzElDofsFEVB_g(vb,fe); }
//=======end of the "COMMON SHAPE PART"==================

      BhomAdjOld.val_g(vb);
      Bmag.curl_g(vb);
      Bmag.val_g(vb);
      Vel.val_g(vb);
      VelAdj.val_g(vb);

 //vector product
          Math::extend(VelAdj._val_g,VelAdj._val_g3D,space_dim);
          Math::cross(Bmag._curl_g3D,VelAdj._val_g3D,curlBXlambda_g3D);

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//============================================================== 

       for (uint i=0; i < BhomAdjOld._ndof[vb]; i++)     {
//======="COMMON tEST PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
        const double phii_g = currgp._phi_ndsQLVB_g[vb][BhomAdjOld._FEord][i];
        for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[vb][BhomAdjOld._FEord][i+idim*BhomAdjOld._ndof[vb]];
//======= END "COMMON tEST PART for QTYZERO" ==========

   	  double BDdphii_g      = Math::dot(  Bmag._val_g,dphiidx_g,space_dim);
	  double lambdaDdphii_g = Math::dot(VelAdj._val_g,dphiidx_g,space_dim);
	  
         for (uint idim=0; idim<space_dim; idim++) {
            const uint irowq = i+idim*BhomAdjOld._ndof[vb];
           currelem._FeM[vb](irowq) += currelem._bc_eldofs[vb][irowq]*dtxJxW_g*(
                             NonStatMHDAD*BhomAdjOld._val_g[idim]*phii_g/dt  //time
                            - S*curlBXlambda_g3D[idim]*phii_g                             //from NS
                            + S*(BDdphii_g*VelAdj._val_g[idim] - lambdaDdphii_g*Bmag._val_g[idim])     //from NS
                        )
                           + (1-currelem._bc_eldofs[vb][irowq])*detb*BhomAdjOld._val_dofs[irowq]; //Dirichlet bc
	   }

if (_Dir_pen_fl == 0)  {
        for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
          const uint irowq = i+idim*BhomAdjOld._ndof[vb];
          currelem._KeM[vb](irowq,irowq) += (1-currelem._bc_eldofs[vb][irowq])*detb;
        }        // end filling diagonal for Dirichlet bc
}
                                           
	 
        for (uint j=0; j<BhomAdjOld._ndof[vb]; j++) {// A element matrix
//======="COMMON SHAPE PART for QTYZERO": ==========
           double                                  phij_g       =      currgp._phi_ndsQLVB_g[vb][BhomAdjOld._FEord][j];
           for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[vb][BhomAdjOld._FEord][j+idim*BhomAdjOld._ndof[vb]];
//======= END "COMMON SHAPE PART for QTYZERO" ==========

           double Lap_g    = Math::dot(dphijdx_g,dphiidx_g,space_dim);
	  double Advphij_g = Math::dot(Vel._val_g,dphijdx_g,space_dim);
          
          for (uint idim=0; idim<space_dim; idim++) { //filled in as 1-2-3 // 4-5-6 // 7-8-9
            int irowq = i+idim*BhomAdjOld._ndof[vb];
            // diagonal blocks [1-5-9]
            currelem._KeM[vb](irowq,j+idim*BhomAdjOld._ndof[vb])
            += currelem._bc_eldofs[vb][irowq]*dtxJxW_g*(
                   NonStatMHDAD*phij_g*phii_g/dt// time
                 + LAP_MHD*IRem*(Lap_g)
                 + (1-LAP_MHD)*IRem*(   Lap_g - dphijdx_g[idim]* dphiidx_g[idim] )
                 - phii_g*(         Advphij_g - dphijdx_g[idim]*Vel._val_g[idim] )
               );
            // block +1 [2-6-7]
            int idimp1=(idim+1)%space_dim;
            currelem._KeM[vb](irowq,j+idimp1*BhomAdjOld._ndof[vb])
            += currelem._bc_eldofs[vb][irowq]*dtxJxW_g*(
                 + (1-LAP_MHD)*IRem*(  -dphijdx_g[idim]* dphiidx_g[idimp1] )
                 - phii_g*(            -dphijdx_g[idim]*Vel._val_g[idimp1] )
               );
#if (DIMENSION==3)
            // block +2 [3-4-8]
            int idimp2=(idim+2)%space_dim;
            currelem._KeM[vb](irowq,j+idimp2*BhomAdjOld._ndof[vb])
            += currelem._bc_eldofs[vb][irowq]*dtxJxW_g*(
                  + (1-LAP_MHD)*IRem*(-dphijdx_g[idim]* dphiidx_g[idimp2] )
                  - phii_g*(          -dphijdx_g[idim]*Vel._val_g[idimp2] )
                 );
#endif
          }
 
        } 
                                      // end A element matrix
      
       for (uint j=0; j < BhomLagMultAdjOld._ndof[vb]; j++) {// B^T element matrix ( p*div(v) )
          const double psij_g =  currgp._phi_ndsQLVB_g[vb][BhomLagMultAdjOld._FEord][j];
          const int jclml= j + space_dim*BhomAdjOld._ndof[vb];
          for (uint idim=0; idim<space_dim; idim++) {
            uint irowq = i+idim*BhomAdjOld._ndof[vb];
            currelem._KeM[vb](irowq,jclml) += currelem._bc_eldofs[vb][irowq]*dtxJxW_g*(-psij_g*dphiidx_g[idim]);
           }
        }
                                     // end B^T element matrix

          if (i<BhomLagMultAdjOld._ndof[vb]) {//  pressure equation (KOMP dp/dt=rho*div) 
          double psii_g = currgp._phi_ndsQLVB_g[vb][BhomLagMultAdjOld._FEord][i];
	  const uint irowl = i+space_dim*BhomAdjOld._ndof[vb];
          currelem._FeM[vb](irowl)=0.;  // rhs
 //             currelem._KeM[vb](irowl,j+space_dim*el_ndof_q)  += dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt;

          for (uint j=0; j<BhomAdjOld._ndof[vb]; j++) { // B element matrix q*div(u)
            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[vb][BhomAdjOld._FEord][j+idim*BhomAdjOld._ndof[vb]];
            for (uint idim=0; idim<space_dim; idim++) currelem._KeM[vb](irowl,j+idim*BhomAdjOld._ndof[vb]) += -dtxJxW_g*psii_g*dphijdx_g[idim]; 
                }
        }
                         // end pressure eq (cont)
      }
//===================================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) =====================
//===================================================================

    }
    // end element gaussian integration loop
    
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
     
            BhomAdjOld.GetElDofsVect(vb,Level);
     BhomLagMultAdjOld.GetElDofsVect(vb,Level);
   
     if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(vb,BhomAdjOld._FEord,currelem._bc_eldofs[vb]); //only the Quadratic Part is modified! /*OK DIR_PEN*/
       
  
    //============ BC =======
    
       int     el_flag[NT] = {0,0};
#if DIMENSION==2
       double el_value[N1T1] = {0.,0.};
#elif DIMENSION==3
       double el_value[N1T3] = {0.,0.,0.,0.}; 
#endif
       double  dbl_pen[NT] = {0.,0.};
   
    
       Bc_GetElFlagValLevSubd(Level,_iproc,iel,el_flag,el_value);

if (_Dir_pen_fl == 1)  { 
       if (el_flag[NN] == 1) {   dbl_pen[NN]=penalty_val; } //normal dirichlet
       if (el_flag[TT] == 1) {   dbl_pen[TT]=penalty_val; } //tangential dirichlet
   }
    
//========END BC============
     

    for (uint qp=0; qp< el_ngauss; qp++) {
//======= "COMMON SHAPE PART"============================
 for (uint fe = 0; fe < QL; fe++)     {
   currgp.SetPhiElDofsFEVB_g (vb,fe,qp); 
   currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp);   }

        const double det   = dt*currgp.JacVectBB_g(vb,xyz);
	const double dtxJxW_g = det * _eqnmap._qrule._weightVB[vb][qp];
//=======end "COMMON SHAPE PART"===================================   
      
      xyz_refbox.val_g(vb);
      BhomLagMultAdjOld._qtyptr->Function_txyz(time,xyz_refbox._val_g/*xyz._val_g*/,BhomLagMultAdjOld._val_g);  //i prefer using the function instead of the p_old vector
   
//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ===================
//==============================================================
           for (uint i=0; i < BhomAdjOld._ndof[vb]; i++) {
 
	const double phii_g = currgp._phi_ndsQLVB_g[vb][BhomAdjOld._FEord][i];

        for (uint idim=0; idim< space_dim; idim++)    {
             uint irowq=i+idim*BhomAdjOld._ndof[vb];
            currelem._FeM[vb](irowq)  += 
          currelem._bc_eldofs[vb][irowq]*           
           dtxJxW_g*(   -1.*/*press_fl*/(1-el_flag[NN])*BhomLagMultAdjOld._val_g[0]*currgp.get_normal_ptr()[idim]*phii_g  //  //OLD VALUES //AAA multiplying int times uint!!!

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

	   for (uint j=0; j< BhomAdjOld._ndof[vb]; j++) {
          const double phij_g = currgp._phi_ndsQLVB_g[vb][ BhomAdjOld._FEord][j];

  currelem._KeM[vb](irowq,j+jdim*BhomAdjOld._ndof[vb]) +=                //projection over the physical (x,y,z) 
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
//========= END FILLING ELEMENT MAT/RHS (i loop) ===============
//==============================================================      
      
    }  //gauss
    
    _A[Level]->add_matrix(currelem._KeM[vb],currelem._el_dof_indices[vb]);
    _b[Level]->add_vector(currelem._FeM[vb],currelem._el_dof_indices[vb]);

 }//elem loop

  }//END BOUNDARY ************************
  
  //******************************DESTROY ALL THE Vect****************************  TODO
  
  
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << _eqname << ": assembled  Level " << Level
           << " with " << _A[Level]->m() << " dofs " << std::endl;
#endif     

  return;
}



} //end namespace femus


