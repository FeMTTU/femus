#include "EqnNSAD.hpp"

#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"

#include "Math.hpp"
#include "MeshTwo.hpp"
#include "GeomEl.hpp"
#include "Domain.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "FETypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "Quantity.hpp"
#include "QTYnumEnum.hpp"
#include "TimeLoop.hpp"
#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"

#include "FemusDefault.hpp"

#include "Opt_conf.hpp"
#include "EqnNS.hpp"
#include "EqnMHD.hpp"
#include "EqnMHDAD.hpp"
#include "EqnMHDCONT.hpp"
#include "OptLoop.hpp"


namespace femus {

/// Constructor.
  EqnNSAD::EqnNSAD(std::vector<Quantity*> int_map_in,
	           MultiLevelProblemTwo& equations_map_in,
                   std::string eqname_in,
                   std::string varname_in):
      SystemTwo(int_map_in,equations_map_in,eqname_in,varname_in)      
      {

//=======  _var_names[]  ===========
    _var_names[0]="lambdax";
    _var_names[1]="lambday";
#if (DIMENSION==3)
    _var_names[2]="lambdaz";
#endif
    _var_names[DIMENSION]="lambdap";

//=======  _refvalue[] ==============   
     _refvalue[0] = _QtyInternalVector[0]->_refvalue[0];        /*1*/ 
     _refvalue[1] = _QtyInternalVector[0]->_refvalue[1];        /*1*/ 
#if (DIMENSION==3)
     _refvalue[2] =  _QtyInternalVector[0]->_refvalue[2];        /*1*/ 
#endif
     _refvalue[DIMENSION] = _QtyInternalVector[1]->_refvalue[0]; /*1*/
     

//============ solver ===============
   for(uint l=0;l<_NoLevels;l++)   _solver[l]->set_solver_type(SOLVERNSAD);
   
//============= DIR PENALTY ===============
   _Dir_pen_fl = NSAD_DIR_PENALTY;
    
    }



//================ DESTRUCTOR    
  EqnNSAD::~EqnNSAD()  {}



 void EqnNSAD::GenMatRhs(const uint Level)  {

    // ========= parameters
  const double alphaVel = _phys.get("alphaVel");
  const double IRe      = 1./_phys.get("Re");

  //=========== Operators 
  double dphijdx_g[DIMENSION];
  double dphiidx_g[DIMENSION];
    double    curlxiXB_g3D[3];

   const double time =  0.;  //_eqnmap._timeloop._curr_time;
   
//======= TIME - STATIONARY OR NOT =======
const int NonStatNSAD = (int) _phys.get("NonStatNSAD");
  const double   dt     = 1.;  //_eqnmap._timeloop._timemap.get("dt");

// //======== GEOMETRICAL ELEMENT =======
  const uint space_dim = _mesh.get_dim();
  const uint  mesh_ord = (int) _mesh.GetRuntimeMap().get("mesh_ord");
  const uint    meshql = (int) _mesh.GetRuntimeMap().get("meshql");  //======== ELEMENT MAPPING =======

  
//========= BCHandling =========
  const double penalty_val = _mesh.GetRuntimeMap().get("penalty_val");    
  
   {//BEGIN VOLUME    
   
   const uint mesh_vb = VV;
    
    CurrentElem       currelem(VV,this,_mesh,_eqnmap._elem_type);
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
   
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    //QTYZERO
    CurrentQuantity VelAdjOld(currgp);
    VelAdjOld._qtyptr   = _QtyInternalVector[QTYZERO];
    VelAdjOld.VectWithQtyFillBasic();
    VelAdjOld.Allocate();
  
    //QTYONE
    CurrentQuantity PressAdjOld(currgp);
    PressAdjOld._qtyptr   = _QtyInternalVector[QTYONE];
    PressAdjOld.VectWithQtyFillBasic();
    PressAdjOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================
  
//=========EXTERNAL QUANTITIES (couplings) =====

//========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();

//========== Quadratic domain, auxiliary  
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
  xyz_refbox.Allocate();
  
  CurrentQuantity Vel(currgp);
    Vel._qtyptr      = _eqnmap._qtymap.get_qty("Qty_Velocity");
    Vel.VectWithQtyFillBasic();
    Vel.Allocate();
  
  CurrentQuantity VelDes(currgp);
    VelDes._qtyptr   = _eqnmap._qtymap.get_qty("Qty_DesVelocity");
    VelDes.VectWithQtyFillBasic();
    VelDes.Allocate();

  CurrentQuantity Bhom(currgp);
    Bhom._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldHom");
    Bhom.VectWithQtyFillBasic();
    Bhom.Allocate();
 
  CurrentQuantity Bext(currgp);
    Bext._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldExt");
    Bext.VectWithQtyFillBasic();
    Bext.Allocate();

//========= auxiliary, must be AFTER Bhom!   //TODO this doesnt have any associated quantity!
  CurrentQuantity Bmag(currgp); //total
    Bmag._dim        = Bhom._dim;               //same as Bhom
    Bmag._FEord      = Bhom._FEord;             //same as Bhom
    Bmag._ndof       = _eqnmap._elem_type[currelem.GetDim()-1][Bmag._FEord]->GetNDofs();
    Bmag.Allocate();
    
//===============
  CurrentQuantity BhomAdj(currgp); 
    BhomAdj._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldHomAdj"); 
    BhomAdj.VectWithQtyFillBasic();
    BhomAdj.Allocate();
    
    //========= END EXTERNAL QUANTITIES =================


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
    
    VelAdjOld.GetElemDofs(Level);  
    PressAdjOld.GetElemDofs(Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),VelAdjOld._FEord,currelem.GetBCDofFlag());  //only the Quadratic Part is modified!
  
    
    if ( Vel._eqnptr != NULL )  Vel.GetElemDofs(Level);
    else                         Vel._qtyptr->FunctionDof(Vel,time,&xyz_refbox._val_dofs[0]);    //give the Hartmann flow, if not solving NS
    if ( Bhom._eqnptr != NULL )  Bhom.GetElemDofs(Level);
    else                         Bhom._qtyptr->FunctionDof(Bhom,time,&xyz_refbox._val_dofs[0]);
    if ( Bext._eqnptr != NULL )  Bext.GetElemDofs(Level);
    else                         Bext._qtyptr->FunctionDof(Bext,time,&xyz_refbox._val_dofs[0]);
    if ( BhomAdj._eqnptr != NULL )  BhomAdj.GetElemDofs(Level);
    else                            BhomAdj._qtyptr->FunctionDof(BhomAdj,time,&xyz_refbox._val_dofs[0]);    
    if ( VelDes._eqnptr != NULL )  VelDes.GetElemDofs(Level);
    else                           VelDes._qtyptr->FunctionDof(VelDes,time,&xyz_refbox._val_dofs[0]);    

 
//======SUM Bhom and Bext  //from now on, you'll only use Bmag //Bmag,Bext and Bhom must have the same orders!
    Math::zeroN(&Bmag._val_dofs[0],Bmag._dim*Bmag._ndof);

    for (uint ivarq=0; ivarq < Bmag._dim; ivarq++)    { //ivarq is like idim
          for (uint d=0; d < Bmag._ndof; d++)    {
          const uint     indxq  =         d + ivarq*Bmag._ndof;
          Bmag._val_dofs[indxq] = Bext._val_dofs[indxq] + Bhom._val_dofs[indxq];
	  }
    }    
//=======    
    
//=======    
///optimal control
    xyz_refbox.SetElemAverage();
  int el_flagdom = ElFlagControl(xyz_refbox._el_average,&_mesh);
//=======    


   const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();

   for (uint qp = 0; qp < el_ngauss; qp++) {
//=======here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     { 
  currgp.SetPhiElDofsFEVB_g (fe,qp);  
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  }  
	  
const double      det = dt*currgp.JacVectVV_g(xyz);   //InvJac: is the same for both QQ and LL!
const double dtxJxW_g = det*_eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { 
  currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp);
  currgp.ExtendDphiDxyzElDofsFEVB_g(fe); 
}
//=======end of the "COMMON SHAPE PART"==================

   VelAdjOld.val_g();
     BhomAdj.curl_g();
        Bmag.val_g();
         Vel.val_g();
      VelDes.val_g();
         Vel.grad_g();

//vector product
        Math::extend(&Bmag._val_g[0],&Bmag._val_g3D[0],space_dim);
        Math::cross(&BhomAdj._curl_g3D[0],&Bmag._val_g3D[0],curlxiXB_g3D);

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
       for (uint i = 0; i < VelAdjOld._ndof; i++)     {
//======="COMMON tEST PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
        const double phii_g = currgp._phi_ndsQLVB_g[VelAdjOld._FEord][i];
        for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim]= currgp._dphidxyz_ndsQLVB_g[VelAdjOld._FEord][i+idim*VelAdjOld._ndof];
//======= END "COMMON tEST PART for QTYZERO" ==========

         for (uint idim=0; idim<space_dim; idim++) {
            const uint irowq=i+idim*VelAdjOld._ndof;  //quadratic rows index
           currelem.Rhs()(irowq) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                            NonStatNSAD*VelAdjOld._val_g[idim]*phii_g/dt  //time
                           - curlxiXB_g3D[idim]*phii_g                    //this is due to the variation of velocity in the MAGNETIC ADVECTION, so it is due to a   NONLINEAR COUPLING "u times B", "MY_STATE times OTHER_STATE"
                           - alphaVel*el_flagdom*(Vel._val_g[idim] - VelDes._val_g[idim])*phii_g    //this is the dependence that counts
                           )
                           + (1-currelem.GetBCDofFlag()[irowq])*detb*VelAdjOld._val_dofs[irowq]; //Dirichlet bc
	   }

if (_Dir_pen_fl == 0)  { 
  for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
          const uint irowq = i+idim*VelAdjOld._ndof;
          currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
        }// end filling diagonal for Dirichlet bc
}
	 
        for (uint j=0; j<VelAdjOld._ndof; j++) {// A element matrix
//======="COMMON SHAPE PART for QTYZERO": ==========
           double                                  phij_g       =      currgp._phi_ndsQLVB_g[VelAdjOld._FEord][j];
           for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[VelAdjOld._FEord][j+idim*VelAdjOld._ndof];
//======= END "COMMON SHAPE PART for QTYZERO" ==========

          double     Lap_g = Math::dot(dphijdx_g,dphiidx_g,space_dim);
	  double Advphii_g = Math::dot(&Vel._val_g[0],dphiidx_g,space_dim);   //TODO can put it outside
          
          for (uint idim=0; idim<space_dim; idim++) { //filled in as 1-2-3 // 4-5-6 // 7-8-9
            int irowq = i+idim*VelAdjOld._ndof;
            // diagonal blocks [1-5-9]
            currelem.Mat()(irowq,j+idim*VelAdjOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                   NonStatNSAD* phij_g*phii_g/dt // time
                   + IRe*(dphijdx_g[idim]*dphiidx_g[idim] + Lap_g)      //Adjoint of D is D, Adjoint of Laplacian is Laplacian
                   + phij_g*phii_g*/*dveldx_g*/Vel._grad_g[idim][idim]  //Adjoint of Advection 1 delta(u) DOT grad(u): adj of nonlinear stuff has 2 TERMS (well, not always)
                   + phij_g*Advphii_g                                   //Adjoint of Advection 2 u DOT grad (delta(u)): adj of nonlinear stuff has 2 TERMS
               );
            // block +1 [2-6-7]
            int idimp1=(idim+1)%space_dim;
            currelem.Mat()(irowq,j+idimp1*VelAdjOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                   + IRe*(dphijdx_g[idim]*dphiidx_g[idimp1])
                   + phij_g*phii_g*/*dveldx_g*/Vel._grad_g[idimp1][idim]
               );
#if (DIMENSION==3)
            // block +2 [3-4-8]
            int idimp2=(idim+2)%space_dim;
            currelem.Mat()(irowq,j+idimp2*VelAdjOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                  + IRe*(dphijdx_g[idim]*dphiidx_g[idimp2])
                  + phij_g*phii_g*/*dveldx_g*/Vel._grad_g[idimp2][idim]
               );
#endif
          }
 
        } 
                                      // end A element matrix
      
       for (uint j=0; j<PressAdjOld._ndof; j++) {// B^T element matrix ( p*div(v) ) (NS eq)
          const double psij_g =  currgp._phi_ndsQLVB_g[PressAdjOld._FEord][j];
          const int jclml= j + space_dim*VelAdjOld._ndof;
          for (uint idim=0; idim<space_dim; idim++) {
            uint irowq = i+idim*VelAdjOld._ndof;
            currelem.Mat()(irowq,jclml) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(-psij_g*dphiidx_g[idim]);
           }
        }
                                     // end B^T element matrix

          if (i<PressAdjOld._ndof) {//  pressure equation (KOMP dp/dt=rho*div) 
          double psii_g = currgp._phi_ndsQLVB_g[PressAdjOld._FEord][i];
	  const uint irowl = i+space_dim*VelAdjOld._ndof;  //vertical offset
          currelem.Rhs()(irowl)=0.;  // rhs
 //             KeM(irowl,j+space_dim*el_ndof_q)  += dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt;

          for (uint j=0; j<VelAdjOld._ndof; j++) { // B element matrix q*div(u)
            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[VelAdjOld._FEord][j+idim*VelAdjOld._ndof];
            for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*VelAdjOld._ndof) += -dtxJxW_g*psii_g*dphijdx_g[idim]; 
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
    //QTYZERO
    CurrentQuantity VelAdjOld(currgp);
    VelAdjOld._qtyptr   = _QtyInternalVector[QTYZERO];
    VelAdjOld.VectWithQtyFillBasic();
    VelAdjOld.Allocate();
  
    //QTYONE
    CurrentQuantity PressAdjOld(currgp);
    PressAdjOld._qtyptr   = _QtyInternalVector[QTYONE];
    PressAdjOld.VectWithQtyFillBasic();
    PressAdjOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================
  
//=========EXTERNAL QUANTITIES (couplings) =====

//========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();

//========== Quadratic domain, auxiliary  
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
  xyz_refbox.Allocate();
  
//========= END EXTERNAL QUANTITIES =================


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

     VelAdjOld.GetElemDofs(Level);
     PressAdjOld.GetElemDofs(Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),VelAdjOld._FEord,currelem.GetBCDofFlag()); //only the Quadratic Part is modified! /*OK DIR_PEN*/
       

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

//checks
//TODO here i should check that the nodal bc dirichlet i put correspond to the element NT flags

       uint press_fl=0;
       Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(currelem.GetBCDofFlag(),VelAdjOld,PressAdjOld,press_fl); //compute the PRESSURE FLAG with the PRESSURE nodal bc flags
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

        const double det   = dt * currgp.JacVectBB_g(xyz);
	const double dtxJxW_g = det * _eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
//=======end "COMMON SHAPE PART"===================================   
      
   xyz_refbox.val_g(); // val_g(vb,xyz);   //CHECK the QUADRATICS!!!!!!!!!
      PressAdjOld._qtyptr->Function_txyz(time,&xyz_refbox._val_g[0]/*xyz._val_g*/,&PressAdjOld._val_g[0]);  //i prefer using the function instead of the p_old vector
      

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
      for (uint i=0; i < VelAdjOld._ndof; i++) {

	const double phii_g = currgp._phi_ndsQLVB_g[VelAdjOld._FEord][i];

         for (uint idim=0; idim< space_dim; idim++)    {
             uint irowq=i+idim*VelAdjOld._ndof;
          currelem.Rhs()(irowq)  += 
            currelem.GetBCDofFlag()[irowq]*           
           dtxJxW_g*(   -1.*/*press_fl*/(1-el_flag[NN])*PressAdjOld._val_g[0]*currgp.get_normal_ptr()[idim]*phii_g  //  //OLD VALUES //AAA multiplying int times uint!!!

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
if (_Dir_pen_fl == 1) {  //much faster than multiplying by _Dir_pen_fl=0 , and much better than removing the code with the #ifdef   
	   for (uint jdim=0; jdim< space_dim; jdim++)    {

	   for (uint j=0; j< VelAdjOld._ndof; j++) {
          const double phij_g = currgp._phi_ndsQLVB_g[ VelAdjOld._FEord][j];

  currelem.Mat()(irowq,j+jdim*VelAdjOld._ndof) +=                //projection over the physical (x,y,z) 
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
      
    }  //gauss
   
    _A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    _b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());

 }//elem loop
   
  }//END BOUNDARY ************************
  
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << _eqname << ": assembled  Level " << Level
           << " with " << _A[Level]->m() << " dofs" << std::endl;
#endif    

return;
}


} //end namespace femus

