#include "EqnMHDAD.hpp"

#include "FemusDefault.hpp"

#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"

#include "Math.hpp"
#include "Quantity.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "GeomEl.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "FETypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "QTYnumEnum.hpp"
#include "TimeLoop.hpp"

#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"

#include "Opt_conf.hpp"
#include "EqnNS.hpp"
#include "EqnNSAD.hpp"
#include "EqnMHD.hpp"
#include "EqnMHDCONT.hpp"

namespace femus {

/// Constructor.
  EqnMHDAD::EqnMHDAD(std::vector<Quantity*> int_map_in,
	           MultiLevelProblemTwo& mg_equations_map_in,
                   std::string eqname_in,
                   std::string varname_in):
      SystemTwo(int_map_in,mg_equations_map_in,eqname_in,varname_in)      
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


  void EqnMHDAD::GenMatRhs(const uint Level)  {

   const double time =  0.;  //_eqnmap._timeloop._curr_time;
   
  //========= parameters
   double IRem =  1./_phys.get("Rem");
   double S    = _phys.get("S");
  
  //=========== Operators 

  double dphijdx_g[DIMENSION];
  double dphiidx_g[DIMENSION];
  double  curlBXlambda_g3D[3]; 
    
//======= TIME - STATIONARY OR NOT =======
const int NonStatMHDAD = (int) _phys.get("NonStatMHDAD");
  const double   dt = 1.; //_eqnmap._timeloop._timemap.get("dt");

//======== GEOMETRICAL ELEMENT =======
  const uint space_dim =       _mesh.get_dim();
  const uint  mesh_ord = (int) _mesh.GetRuntimeMap().get("mesh_ord");
  const uint    meshql = (int) _mesh.GetRuntimeMap().get("meshql");  //======== ELEMENT MAPPING =======

//========= BCHandling =========
  const double penalty_val =   _mesh.GetRuntimeMap().get("penalty_val");    

   {//BEGIN VOLUME    

  const uint mesh_vb = VV;
  
    CurrentElem       currelem(VV,this,_mesh,_eqnmap._elem_type);
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
   
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
     //QTYZERO
    CurrentQuantity BhomAdjOld(currgp);
    BhomAdjOld._qtyptr   = _QtyInternalVector[QTYZERO];
    BhomAdjOld.VectWithQtyFillBasic();
    BhomAdjOld.Allocate();    
  
    //QTYONE
    CurrentQuantity BhomLagMultAdjOld(currgp);
    BhomLagMultAdjOld._qtyptr   = _QtyInternalVector[QTYONE];
    BhomLagMultAdjOld.VectWithQtyFillBasic();
    BhomLagMultAdjOld.Allocate();    
//========= END INTERNAL QUANTITIES (unknowns of the equation) ================= 

//=========EXTERNAL QUANTITIES (couplings) =====

 //========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();    

//========== Quadratic domain, auxiliary  
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
  xyz_refbox.Allocate();    

  //==========     
    CurrentQuantity Vel(currgp);
    Vel._qtyptr      = _eqnmap._qtymap.get_qty("Qty_Velocity");
    Vel.VectWithQtyFillBasic();
    Vel.Allocate();    
 
    //==========    
    CurrentQuantity VelAdj(currgp);
    VelAdj._qtyptr      = _eqnmap._qtymap.get_qty("Qty_VelocityAdj");
    VelAdj.VectWithQtyFillBasic();
    VelAdj.Allocate();    
   
    //==========    
    CurrentQuantity Bhom(currgp);
    Bhom._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldHom");
    Bhom.VectWithQtyFillBasic();
    Bhom.Allocate();    

//=========
    CurrentQuantity Bext(currgp);
    Bext._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldExt");
    Bext.VectWithQtyFillBasic();
    Bext.Allocate();    
  
//========= auxiliary, must be AFTER Bhom! //TODO this is an example of  Vect which is not associated to a Quantity
    CurrentQuantity Bmag(currgp); //total
    Bmag._dim        = Bhom._dim;               //same as Bhom
    Bmag._FEord      = Bhom._FEord;             //same as Bhom
    Bmag._ndof       = _eqnmap._elem_type[currelem.GetDim()-1][Bmag._FEord]->GetNDofs();
    Bmag.Allocate();    
    
//========= END EXTERNAL QUANTITIES =================

    const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level+1];
    const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level];

//=====================
//=====================    
    
  for (uint iel=0; iel < (nel_e - nel_b); iel++) {

    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.set_el_nod_conn_lev_subd(Level,_iproc,iel);
    currelem.SetMidpoint();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);    
    _mesh.TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]);

    currelem.SetElDofsBc(Level);
    
           BhomAdjOld.GetElemDofs(Level);
    BhomLagMultAdjOld.GetElemDofs(Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),BhomAdjOld._FEord,currelem.GetBCDofFlag());  //only the Quadratic Part is modified!
    
    
     if ( Vel._eqnptr != NULL )      Vel.GetElemDofs(Level);
    else                             Vel._qtyptr->FunctionDof(Vel,time,&xyz_refbox._val_dofs[0]);
    if ( VelAdj._eqnptr != NULL ) VelAdj.GetElemDofs(Level);
    else                          VelAdj._qtyptr->FunctionDof(VelAdj,time,&xyz_refbox._val_dofs[0]);
    if ( Bhom._eqnptr != NULL )     Bhom.GetElemDofs(Level);
    else                            Bhom._qtyptr->FunctionDof(Bhom,time,&xyz_refbox._val_dofs[0]);
    if ( Bext._eqnptr != NULL )     Bext.GetElemDofs(Level);
    else                            Bext._qtyptr->FunctionDof(Bext,time,&xyz_refbox._val_dofs[0]);

//======SUM Bhom and Bext  //from now on, you'll only use Bmag //Bmag,Bext and Bhom must have the same orders!
    Math::zeroN(&Bmag._val_dofs[0],Bmag._dim*Bmag._ndof);

    for (uint ivarq=0; ivarq < Bmag._dim; ivarq++)    { //ivarq is like idim
          for (uint d=0; d < Bmag._ndof; d++)    {
          const uint     indxq  =         d + ivarq*Bmag._ndof;
          Bmag._val_dofs[indxq] = Bext._val_dofs[indxq] + Bhom._val_dofs[indxq];
	  }
    }    
//=======


    const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();
    
    for (uint qp = 0; qp < el_ngauss; qp++) {
//=======here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     {          
  currgp.SetPhiElDofsFEVB_g (fe,qp); 
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}  
	  
const double      det = dt*currgp.JacVectVV_g(xyz);   //InvJac: is the same for both QQ and LL!
const double dtxJxW_g = det * _eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { 
  currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp);
  currgp.ExtendDphiDxyzElDofsFEVB_g(fe);
}
//=======end of the "COMMON SHAPE PART"==================

      BhomAdjOld.val_g();
      Bmag.curl_g();
      Bmag.val_g();
      Vel.val_g();
      VelAdj.val_g();

 //vector product
          Math::extend(&VelAdj._val_g[0],&VelAdj._val_g3D[0],space_dim);
          Math::cross(&Bmag._curl_g3D[0],&VelAdj._val_g3D[0],curlBXlambda_g3D);

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//============================================================== 

       for (uint i=0; i < BhomAdjOld._ndof; i++)     {
//======="COMMON tEST PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
        const double phii_g = currgp._phi_ndsQLVB_g[BhomAdjOld._FEord][i];
        for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BhomAdjOld._FEord][i+idim*BhomAdjOld._ndof];
//======= END "COMMON tEST PART for QTYZERO" ==========

   	  double BDdphii_g      = Math::dot(  &Bmag._val_g[0],dphiidx_g,space_dim);
	  double lambdaDdphii_g = Math::dot(&VelAdj._val_g[0],dphiidx_g,space_dim);
	  
         for (uint idim=0; idim<space_dim; idim++) {
            const uint irowq = i+idim*BhomAdjOld._ndof;
           currelem.Rhs()(irowq) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                             NonStatMHDAD*BhomAdjOld._val_g[idim]*phii_g/dt  //time
                            - S*curlBXlambda_g3D[idim]*phii_g                             //from NS
                            + S*(BDdphii_g*VelAdj._val_g[idim] - lambdaDdphii_g*Bmag._val_g[idim])     //from NS
                        )
                           + (1-currelem.GetBCDofFlag()[irowq])*detb*BhomAdjOld._val_dofs[irowq]; //Dirichlet bc
	   }

if (_Dir_pen_fl == 0)  {
        for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
          const uint irowq = i+idim*BhomAdjOld._ndof;
          currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
        }        // end filling diagonal for Dirichlet bc
}
                                           
	 
        for (uint j=0; j<BhomAdjOld._ndof; j++) {// A element matrix
//======="COMMON SHAPE PART for QTYZERO": ==========
           double                                  phij_g       =      currgp._phi_ndsQLVB_g[BhomAdjOld._FEord][j];
           for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BhomAdjOld._FEord][j+idim*BhomAdjOld._ndof];
//======= END "COMMON SHAPE PART for QTYZERO" ==========

           double Lap_g    = Math::dot(dphijdx_g,dphiidx_g,space_dim);
	  double Advphij_g = Math::dot(&Vel._val_g[0],dphijdx_g,space_dim);
          
          for (uint idim=0; idim<space_dim; idim++) { //filled in as 1-2-3 // 4-5-6 // 7-8-9
            int irowq = i+idim*BhomAdjOld._ndof;
            // diagonal blocks [1-5-9]
            currelem.Mat()(irowq,j+idim*BhomAdjOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                   NonStatMHDAD*phij_g*phii_g/dt// time
                 + LAP_MHD*IRem*(Lap_g)
                 + (1-LAP_MHD)*IRem*(   Lap_g - dphijdx_g[idim]* dphiidx_g[idim] )
                 - phii_g*(         Advphij_g - dphijdx_g[idim]*Vel._val_g[idim] )
               );
            // block +1 [2-6-7]
            int idimp1=(idim+1)%space_dim;
            currelem.Mat()(irowq,j+idimp1*BhomAdjOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                 + (1-LAP_MHD)*IRem*(  -dphijdx_g[idim]* dphiidx_g[idimp1] )
                 - phii_g*(            -dphijdx_g[idim]*Vel._val_g[idimp1] )
               );
#if (DIMENSION==3)
            // block +2 [3-4-8]
            int idimp2=(idim+2)%space_dim;
            currelem.Mat()(irowq,j+idimp2*BhomAdjOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                  + (1-LAP_MHD)*IRem*(-dphijdx_g[idim]* dphiidx_g[idimp2] )
                  - phii_g*(          -dphijdx_g[idim]*Vel._val_g[idimp2] )
                 );
#endif
          }
 
        } 
                                      // end A element matrix
      
       for (uint j=0; j < BhomLagMultAdjOld._ndof; j++) {// B^T element matrix ( p*div(v) )
          const double psij_g =  currgp._phi_ndsQLVB_g[BhomLagMultAdjOld._FEord][j];
          const int jclml= j + space_dim*BhomAdjOld._ndof;
          for (uint idim=0; idim<space_dim; idim++) {
            uint irowq = i+idim*BhomAdjOld._ndof;
            currelem.Mat()(irowq,jclml) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(-psij_g*dphiidx_g[idim]);
           }
        }
                                     // end B^T element matrix

          if (i<BhomLagMultAdjOld._ndof) {//  pressure equation (KOMP dp/dt=rho*div) 
          double psii_g = currgp._phi_ndsQLVB_g[BhomLagMultAdjOld._FEord][i];
	  const uint irowl = i+space_dim*BhomAdjOld._ndof;
          currelem.Rhs()(irowl)=0.;  // rhs
 //             currelem.Mat()(irowl,j+space_dim*el_ndof_q)  += dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt;

          for (uint j=0; j<BhomAdjOld._ndof; j++) { // B element matrix q*div(u)
            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[BhomAdjOld._FEord][j+idim*BhomAdjOld._ndof];
            for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*BhomAdjOld._ndof) += -dtxJxW_g*psii_g*dphijdx_g[idim]; 
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
    CurrentQuantity BhomAdjOld(currgp);
    BhomAdjOld._qtyptr   = _QtyInternalVector[QTYZERO];
    BhomAdjOld.VectWithQtyFillBasic();
    BhomAdjOld.Allocate();    
  
    //QTYONE
    CurrentQuantity BhomLagMultAdjOld(currgp);
    BhomLagMultAdjOld._qtyptr   = _QtyInternalVector[QTYONE];
    BhomLagMultAdjOld.VectWithQtyFillBasic();
    BhomLagMultAdjOld.Allocate();    
//========= END INTERNAL QUANTITIES (unknowns of the equation) ================= 

//=========EXTERNAL QUANTITIES (couplings) =====

 //========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();    

//========== Quadratic domain, auxiliary  
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
  xyz_refbox.Allocate();    
   
//========= END EXTERNAL QUANTITIES =================

    const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level+1];
    const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level];

//=====================
//=====================    
    
 for (uint iel=0;iel < (nel_e - nel_b) ; iel++) {

     currelem.Mat().zero();
     currelem.Rhs().zero();

     currelem.set_el_nod_conn_lev_subd(Level,_iproc,iel);
     currelem.SetMidpoint();
     
     currelem.ConvertElemCoordsToMappingOrd(xyz);
    _mesh.TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]);

     currelem.SetElDofsBc(Level);
     
            BhomAdjOld.GetElemDofs(Level);
     BhomLagMultAdjOld.GetElemDofs(Level);
   
     if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),BhomAdjOld._FEord,currelem.GetBCDofFlag()); //only the Quadratic Part is modified! /*OK DIR_PEN*/
       
  
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
     
    const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();

    for (uint qp=0; qp< el_ngauss; qp++) {
//======= "COMMON SHAPE PART"============================
 for (uint fe = 0; fe < QL; fe++)     {
   currgp.SetPhiElDofsFEVB_g (fe,qp); 
   currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);   }

        const double det   = dt*currgp.JacVectBB_g(xyz);
	const double dtxJxW_g = det * _eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
//=======end "COMMON SHAPE PART"===================================   
      
      xyz_refbox.val_g();
      BhomLagMultAdjOld._qtyptr->Function_txyz(time,&xyz_refbox._val_g[0]/*xyz._val_g*/,&BhomLagMultAdjOld._val_g[0]);  //i prefer using the function instead of the p_old vector
   
//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ===================
//==============================================================
           for (uint i=0; i < BhomAdjOld._ndof; i++) {
 
	const double phii_g = currgp._phi_ndsQLVB_g[BhomAdjOld._FEord][i];

        for (uint idim=0; idim< space_dim; idim++)    {
             uint irowq=i+idim*BhomAdjOld._ndof;
            currelem.Rhs()(irowq)  += 
          currelem.GetBCDofFlag()[irowq]*           
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

	   for (uint j=0; j< BhomAdjOld._ndof; j++) {
          const double phij_g = currgp._phi_ndsQLVB_g[ BhomAdjOld._FEord][j];

  currelem.Mat()(irowq,j+jdim*BhomAdjOld._ndof) +=                //projection over the physical (x,y,z) 
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
    
    _A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    _b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());

 }//elem loop

 
  }//END BOUNDARY ************************
  
  
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << _eqname << ": assembled  Level " << Level
           << " with " << _A[Level]->m() << " dofs " << std::endl;
#endif     

  return;
}



} //end namespace femus


