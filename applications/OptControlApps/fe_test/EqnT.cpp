#include "EqnT.hpp"

#include "FemusDefault.hpp"

#include "NumericVector.hpp"
#include "DenseVector.hpp"
#include "SparseMatrix.hpp"
#include "DenseMatrix.hpp"
#include "LinearSolverM.hpp"

#include "Files.hpp"
#include "Math.hpp"
#include "Physics.hpp"
#include "MeshTwo.hpp"
#include "GeomEl.hpp"
#include "EquationsMap.hpp"
#include "EqnBase.hpp"
#include "FETypeEnum.hpp"
#include "FEElemBase.hpp"
#include "NormTangEnum.hpp"
#include "VBTypeEnum.hpp"
#include "QTYnumEnum.hpp"
#include "Domain.hpp"
#include "TimeLoop.hpp"
#include "CurrGaussPoint.hpp"
#include "CurrElem.hpp"
#include "paral.hpp"

// application
#include "TempQuantities.hpp"
#include "TempPhysics.hpp"



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

//=======  _var_names: they are the names of the quantities which are unkwnowns to this equation  ===========
   for (uint i=0; i<int_map_in.size(); i++)  _var_names[i]=int_map_in[i]->_name;
  
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
  CurrGaussPointBase & currgp = CurrGaussPointBase::build(_eqnmap, _mesh.get_dim());
  
  TempPhysics* myphys; myphys = static_cast<TempPhysics*>(&_phys);

//========== PROCESSOR INDEX
  const uint myproc = _iproc;

//========= BCHandling =========
  const double penalty_val = _mesh._mesh_rtmap.get("penalty_val");    

  //======== ELEMENT MAPPING =======
  const uint space_dim =       _mesh.get_dim();
  const uint  meshql   = (int) _mesh._mesh_rtmap.get("meshql");
  const uint  mesh_ord = (int) _mesh._mesh_rtmap.get("mesh_ord");

  
//=========INTERNAL QUANTITIES (unknowns of the equation) =========     
    QuantityLocal Tempold(currgp,currelem);
    Tempold._qtyptr   = _QtyInternalVector[0]; 
    Tempold.VectWithQtyFillBasic();
    Tempold._val_dofs = new double[Tempold._dim*Tempold._ndof[vb]];
    Tempold._val_g    = new double[Tempold._dim];

// //=========INTERNAL QUANTITIES (unknowns of the equation) =========     
//     QuantityLocal Temp2(currgp,currelem);
//     Temp2._qtyptr   = _QtyInternalVector[1]; 
//     Temp2.VectWithQtyFillBasic();
//     Temp2._val_dofs = new double[Temp2._dim*Temp2._ndof[vb]];
//     Temp2._val_g    = new double[Temp2._dim];
// 
// //=========INTERNAL QUANTITIES (unknowns of the equation) =========     
//     QuantityLocal Temp3(currgp,currelem);
//     Temp3._qtyptr   = _QtyInternalVector[2]; 
//     Temp3.VectWithQtyFillBasic();
//     Temp3._val_dofs = new double[Temp3._dim*Temp3._ndof[vb]];
//     Temp3._val_g    = new double[Temp3._dim];
    
    
    //=========EXTERNAL QUANTITIES (couplings) =====
    //========= //DOMAIN MAPPING
    QuantityLocal xyz(currgp,currelem);  //no quantity
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof[VV] = _AbstractFE[VV][xyz._FEord]->_myelems[VV]->GetNDofs();
    xyz._ndof[BB] = _AbstractFE[VV][xyz._FEord]->_myelems[BB]->GetNDofs();
    xyz._val_dofs = new double[xyz._dim*xyz._ndof[vb]];
    xyz._val_g    = new double[xyz._dim];

    //==================Quadratic domain, auxiliary, must be QUADRATIC!!! ==========
  QuantityLocal xyz_refbox(currgp,currelem);  //no quantity
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof[VV] = _mesh._GeomEl[VV]._elnds[xyz_refbox._FEord];
  xyz_refbox._ndof[BB] = _mesh._GeomEl[BB]._elnds[xyz_refbox._FEord];
  xyz_refbox._val_dofs = new double[xyz_refbox._dim*xyz_refbox._ndof[vb]]; 
  xyz_refbox._val_g    = new double[xyz_refbox._dim];
  xyz_refbox._el_average.resize(VB);
  for (uint i=0; i<VB; i++)  xyz_refbox._el_average[i].resize(xyz_refbox._dim);
  //==================

  //==== AUXILIARY ==============
    double* dphijdx_g = new double[space_dim];
    double* dphiidx_g = new double[space_dim];
    double* dphijdx_gLL = new double[space_dim];
    double* dphiidx_gLL = new double[space_dim];
    double* dphijdx_gKK = new double[space_dim];
    double* dphiidx_gKK = new double[space_dim];
   //-----Nonhomogeneous Neumann------
 // Qflux = - k grad(T) by definition
//  QfluxDOTn>0: energy flows outside (cooling)  QfluxDOTn<0: energy flows inside (heating)
    double* Qflux_g = new double[space_dim];

  /// b) Element  Loop over the volume (n_elem)
   const uint el_ngauss = _eqnmap._qrule[vb].GetGaussPointsNumber();
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
//     Temp2.GetElDofsVect(vb,Level);
//     Temp3.GetElDofsVect(vb,Level);
    
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

    for (uint qp=0; qp< el_ngauss; qp++) {

//======= "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)   { 
  currgp.SetPhiElDofsFEVB_g (vb,fe,qp);
  currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp); 
}
	  
const double      det = currgp.JacVectVV_g(vb,xyz);
const double dtxJxW_g = det*_eqnmap._qrule[vb].GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { currgp.SetDPhiDxyzElDofsFEVB_g   (vb,fe,qp); }
for (uint fe = 0; fe < QL; fe++)     { currgp.ExtendDphiDxyzElDofsFEVB_g(vb,fe); }
//======= end of the "COMMON SHAPE PART"==================

 	Tempold.val_g(vb); 
//  	  Temp2.val_g(vb); 
//  	  Temp3.val_g(vb); 
	   
	   // always remember to get the dofs for the variables you use!
           // The point is that you fill the dofs with different functions...
           // you should need a flag to check if the dofs have been correctly filled
 
	  
      for (uint i=0; i < Tempold._ndof[vb]/*the maximum number is for biquadratic*/; i++)     {

        const double phii_g   = currgp._phi_ndsQLVB_g[vb][Tempold._FEord][i];
//         const double phii_gLL = currgp._phi_ndsQLVB_g[vb][Temp2._FEord][i];
//         const double phii_gKK = currgp._phi_ndsQLVB_g[vb][Temp3._FEord][i];

        for (uint idim = 0; idim < space_dim; idim++) dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[vb][Tempold._FEord][i+idim*Tempold._ndof[vb]];

//=========== FIRST ROW ===============
        currelem._FeM[vb](i) +=      
           currelem._bc_eldofs[vb][i]*dtxJxW_g*( 
                7.*phii_g
	  )
	   + (1-currelem._bc_eldofs[vb][i])*detb*(Tempold._val_dofs[i]);
        
        currelem._KeM[vb](i,i) +=  (1-currelem._bc_eldofs[vb][i])*detb;

// // // //========= SECOND ROW =====================
// // // 	 int ip1 = i + Tempold._ndof[vb]; 
// // // 	 
// // // 	if (i < _AbstractFE[ Temp2._FEord ]->_ndof[vb]) { 
// // // 	 currelem._FeM[vb](ip1) +=      
// // //            currelem._bc_eldofs[vb][ip1]*dtxJxW_g*( 
// // //                 0.07*phii_gLL
// // // 	  )
// // // 	   + (1-currelem._bc_eldofs[vb][ip1])*detb*(Temp2._val_dofs[i]);
// // //         
// // //          currelem._KeM[vb](ip1,ip1) +=  (1-currelem._bc_eldofs[vb][ip1])*detb;
// // // 	}
// // // 	
// // // //======= THIRD ROW ===================================
// // // 	 int ip2 = i + Tempold._ndof[vb] + Temp2._ndof[vb];
// // // 	 
// // // 	if (i < _AbstractFE[ Temp3._FEord ]->_ndof[vb]) { 
// // //            currelem._FeM[vb](ip2) +=      
// // //            currelem._bc_eldofs[vb][ip2]*dtxJxW_g*( 
// // //                 0.07*phii_gKK
// // // 	     )
// // // 	   + (1-currelem._bc_eldofs[vb][ip2])*detb*(Temp3._val_dofs[i]);
// // //         
// // //         currelem._KeM[vb](ip2,ip2) +=  (1-currelem._bc_eldofs[vb][ip2])*detb;
// // // 	}
	
	 // Matrix Assemblying ---------------------------
        for (uint j=0; j<Tempold._ndof[vb]; j++) {
          double phij_g   = currgp._phi_ndsQLVB_g[vb][Tempold._FEord][j];
//           double phij_gLL = currgp._phi_ndsQLVB_g[vb][Temp2._FEord][j];
//           double phij_gKK = currgp._phi_ndsQLVB_g[vb][Temp3._FEord][j];
	  
        for (uint idim = 0; idim < space_dim; idim++)   {
	  dphijdx_g  [idim] = currgp._dphidxyz_ndsQLVB_g[vb][Tempold._FEord][j+idim*Tempold._ndof[vb]]; 
// // // 	  dphijdx_gLL[idim] = currgp._dphidxyz_ndsQLVB_g[vb][Temp2._FEord]  [j+idim*Temp2._ndof[vb]]; 
// // // 	  dphijdx_gKK[idim] = currgp._dphidxyz_ndsQLVB_g[vb][Temp3._FEord]  [j+idim*Temp3._ndof[vb]]; 
          }
	  
	  
          double Lap_g   = Math::dot(dphijdx_g,dphiidx_g,space_dim);
          double Lap_gLL = Math::dot(dphijdx_gLL,dphiidx_gLL,space_dim);
          double Lap_gKK = Math::dot(dphijdx_gKK,dphiidx_gKK,space_dim);

	    int ip1 = i + Tempold._ndof[vb];
	    int jp1 = j + Tempold._ndof[vb];
// // // 	    int ip2 = i + Tempold._ndof[vb] + Temp2._ndof[vb];
// // // 	    int jp2 = j + Tempold._ndof[vb] + Temp2._ndof[vb];

 
//============ FIRST ROW state  delta T ===============
//======= DIAGONAL =============================
	   currelem._KeM[vb](i,j) +=        
            currelem._bc_eldofs[vb][i]*dtxJxW_g*( 
             Lap_g  
            );

// // // //=========== SECOND ROW  =============
// // // //===== DIAGONAL ===========================
// // //  	if (i < _AbstractFE[ Temp2._FEord ]->_ndof[vb])  { 
// // //   	if (j < _AbstractFE[ Temp2._FEord ]->_ndof[vb]) { 
// // //        currelem._KeM[vb](ip1,jp1) +=        
// // //             currelem._bc_eldofs[vb][ip1]*
// // //             dtxJxW_g*( 
// // //               Lap_gLL
// // //             ); 
// // // 	  }
// // // 	}
// // // //============= THIRD ROW  =============
// // // //======= DIAGONAL ==================
// // // 	if (i < _AbstractFE[ Temp3._FEord ]->_ndof[vb])  { 
// // //   	if (j < _AbstractFE[ Temp3._FEord ]->_ndof[vb]) { 
// // //           currelem._KeM[vb](ip2,jp2) +=        
// // //             currelem._bc_eldofs[vb][ip2]*
// // //               dtxJxW_g*( 
// // //               phij_gKK*phii_gKK
// // //             + Lap_gKK
// // //                
// // //             ); 
// // // 	  }
// // // 	}
	
	
	    
        }  //end j (col)
      }   //end i (row)
    } // end of the quadrature point qp-loop

    currelem._KeM[vb].print_scientific(std::cout);
    
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
// // //          Temp2.GetElDofsVect(vb,Level);
// // //          Temp3.GetElDofsVect(vb,Level);

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
        const double  det   = currgp.JacVectBB_g(vb,xyz);
        const double dtxJxW_g = det * _eqnmap._qrule[vb].GetGaussWeight(qp);
//=======end "COMMON SHAPE PART"===================================

       xyz.val_g(vb);
       
       static_cast<Temperature*>(_eqnmap._qtymap.get_qty("Qty_Temperature"))->heatflux_txyz(time,xyz._val_g,Qflux_g);
   
	Tempold.val_g(vb); //For the penalty Dirichlet //i need this for interpolating the old function at the gauss point

	   double QfluxDn_g=Math::dot( Qflux_g,currgp.get_normal_ptr(),space_dim );
	 
        for (uint i=0; i<Tempold._ndof[vb]; i++) {
   	const double phii_g =  currgp._phi_ndsQLVB_g[vb][Tempold._FEord][i]; 
	
       currelem._FeM[vb](i) +=
          0.*(currelem._bc_eldofs[vb][i]*
         el_Neum_flag*dtxJxW_g*(-QfluxDn_g)*phii_g    // beware of the sign  //this integral goes in the first equation
	 + el_penalty*dtxJxW_g*Tempold._val_g[0]*phii_g)  //clearly, if you continue using bc=0 for setting nodal Dirichlet, this must go outside
	 ; 
	 
         if (_Dir_pen_fl == 1) {
            for (uint j=0; j<Tempold._ndof[vb]; j++) {
               double phij_g = currgp._phi_ndsQLVB_g[vb][Tempold._FEord][j];
	       currelem._KeM[vb](i,j) += 0.*el_penalty*dtxJxW_g*phij_g*phii_g;
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
// // //   delete []  Temp3._val_g;    delete []  Temp3._val_dofs;
// // //   delete []  Temp2._val_g;    delete []  Temp2._val_dofs;
  delete []  Tempold._val_g;    delete []  Tempold._val_dofs;
  delete []  xyz._val_g;        delete []  xyz._val_dofs;

  delete [] dphijdx_g;
  delete [] dphiidx_g;
  delete [] dphijdx_gLL;
  delete [] dphiidx_gLL;
  delete [] dphijdx_gKK;
  delete [] dphiidx_gKK;


  
#ifdef DEFAULT_PRINT_INFO
  std::cout << " Matrix and RHS assembled for equation " << _eqname
            << " Level "<< Level << " dofs " << _A[Level]->n() << " vb = " << vb << std::endl;
#endif

  return;
}






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
    CurrGaussPointBase & currgp = CurrGaussPointBase::build(_eqnmap, _mesh.get_dim());

  //====== Physics cast
  TempPhysics *optphys; optphys = static_cast<TempPhysics*>(&_phys);

  //====== processor index
  const uint myproc = _iproc;
  const uint space_dim =      _mesh.get_dim();
  const uint mesh_ord = (int) _mesh._mesh_rtmap.get("mesh_ord");  
  const uint meshql   = (int) _mesh._mesh_rtmap.get("meshql");   //======== ELEMENT MAPPING =======
 
//========= DOMAIN MAPPING
    QuantityLocal xyz(currgp,currelem);
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof[VV] = _AbstractFE[VV][xyz._FEord]->_myelems[VV]->GetNDofs();
    xyz._ndof[BB] = _AbstractFE[VV][xyz._FEord]->_myelems[BB]->GetNDofs();
    xyz._val_dofs = new double[xyz._dim*xyz._ndof[vb]];
    xyz._val_g    = new double[xyz._dim];

//========== Quadratic domain, auxiliary  
  QuantityLocal xyz_refbox(currgp,currelem);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof[VV] = _mesh._GeomEl[VV]._elnds[xyz_refbox._FEord];
  xyz_refbox._ndof[BB] = _mesh._GeomEl[BB]._elnds[xyz_refbox._FEord];
  xyz_refbox._val_dofs = new double[xyz_refbox._dim*xyz_refbox._ndof[vb]]; 
  xyz_refbox._val_g    = new double[xyz_refbox._dim];
  xyz_refbox._el_average.resize(VB);
  for (uint i=0; i<VB; i++)  xyz_refbox._el_average[i].resize(xyz_refbox._dim);
  
    
   double integral = 0.;

      const uint el_ngauss = _eqnmap._qrule[vb].GetGaussPointsNumber();

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
 

    for (uint qp = 0; qp < el_ngauss; qp++) {

     for (uint fe = 0; fe < QL; fe++)     {  currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp);  }  
     
   const double  Jac_g = currgp.JacVectVV_g(vb,xyz);  //not xyz_refbox!      
   const double  wgt_g = _eqnmap._qrule[vb].GetGaussWeight(qp);

     for (uint fe = 0; fe < QL; fe++)     {          currgp.SetPhiElDofsFEVB_g (vb,fe,qp);  }

       integral += el_flagdom*wgt_g*Jac_g*1.;
   
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
    
    std::string intgr_fname = _eqnmap._files._app_path + "/" + DEFAULT_OUTPUTDIR + "/" + "alpha";
  
	std::ofstream intgr_fstream;
    
    if (paral::get_rank() ==0 ){ 
      intgr_fstream.open(intgr_fname.c_str(),ios_base::app); 
      intgr_fstream << _eqnmap._files._output_time << " " << optphys->_physrtmap.get("alphaT") << " " << optphys->_physrtmap.get("injsuc")<< " "  << J << " " << std::endl ; 
      intgr_fstream.close();  //you have to close to disassociate the file from the stream
}
 
    
  return J;  
    
}


/////////////////////////

double EqnT::ComputeNormControl (const uint vb, const uint Level, const uint reg_ord ) {

  //reg_ord = 0: L2
  //reg_ord = 1: H1

    CurrElem       currelem(*this,_eqnmap);  //TODO in these functions you only need the GEOMETRIC PART, not the DOFS PART
    CurrGaussPointBase & currgp = CurrGaussPointBase::build(_eqnmap, _mesh.get_dim());
  
  // processor index
  const uint myproc = _iproc;
  const uint space_dim =       _mesh.get_dim();
  const uint mesh_ord  = (int) _mesh._mesh_rtmap.get("mesh_ord");  
  const uint meshql    = (int) _mesh._mesh_rtmap.get("meshql");    //======== ELEMENT MAPPING =======

  
//======Functions in the integrand ============

//========= DOMAIN MAPPING
    QuantityLocal xyz(currgp,currelem);
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof[VV] = _AbstractFE[VV][xyz._FEord]->_myelems[VV]->GetNDofs();
    xyz._ndof[BB] = _AbstractFE[VV][xyz._FEord]->_myelems[BB]->GetNDofs();
    xyz._val_dofs = new double[xyz._dim*xyz._ndof[vb]];
    xyz._val_g    = new double[xyz._dim];

//========== Quadratic domain, auxiliary  
  QuantityLocal xyz_refbox(currgp,currelem);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof[VV] = _mesh._GeomEl[VV]._elnds[xyz_refbox._FEord];
  xyz_refbox._ndof[BB] = _mesh._GeomEl[BB]._elnds[xyz_refbox._FEord];
  xyz_refbox._val_dofs = new double[xyz_refbox._dim*xyz_refbox._ndof[vb]]; 
  xyz_refbox._val_g    = new double[xyz_refbox._dim];
  
    
   double integral = 0.;

//loop over the geom el types
      const uint el_ngauss = _eqnmap._qrule[vb].GetGaussPointsNumber();

//parallel sum
    const uint nel_e = _mesh._off_el[vb][_NoLevels*myproc+Level+1];
    const uint nel_b = _mesh._off_el[vb][_NoLevels*myproc+Level];
  
    for (int iel=0; iel < (nel_e - nel_b); iel++) {

      currelem.get_el_nod_conn_lev_subd(vb,Level,myproc,iel);
      currelem.get_el_DofObj_lev_subd(vb,Level,myproc,iel);
      currelem.get_el_ctr(vb);

      currelem.ConvertElemCoordsToMappingOrd(vb,xyz);
      _mesh.TransformElemNodesToRef(vb,currelem._xx_nds[vb],xyz_refbox._val_dofs);
     
  for (uint qp = 0; qp < el_ngauss; qp++) {

     for (uint fe = 0; fe < QL; fe++)     {
       currgp.SetPhiElDofsFEVB_g (vb,fe,qp);
       currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp);  
    }  
     
      const double  Jac_g = currgp.JacVectVV_g(vb,xyz);  //not xyz_refbox!      
      const double  wgt_g = _eqnmap._qrule[vb].GetGaussWeight(qp);


  double deltau_squarenorm_g = 0.;
   deltau_squarenorm_g += (1-reg_ord) * 1.; 

   for (uint idim = 0; idim < space_dim; idim++)  {   deltau_squarenorm_g  += reg_ord * 1. ;   }

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
