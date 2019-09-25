#include "CurrentQuantity.hpp"

#include "Quantity.hpp"
#include "ElemType.hpp"
#include "MultiLevelProblem.hpp"
#include "SystemTwo.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "CurrentElem.hpp"


#include "NumericVector.hpp"  //TODO this is needed for x_old


namespace femus {




     CurrentQuantity::CurrentQuantity(const CurrentElem<double> & curr_el_in)
     : _currEl(curr_el_in) 
     {  }

   
    //============================
//clearly, here we have the dphi at the dofs and the dof values
//so, we have to extend both the dphi and the dofs
//first we extend the dofs in here
//AAA: the dphi are assumed to be extended automatically before
//so, it must be called AFTER FILLING _dphidxyz_ndsQLVB_g3D
//instead,the myVect._val_dofs3D is filled here
//actually, that vector is only for service purposes 
//so it could be also translated into some temporary vector whenever needed
 /// @todo this function only works for a VECTOR quantity
 void CurrentQuantity::curl_g() {
   
  const uint       ord = _FEord;
  const uint el_nnodes = _ndof;
  const uint el_ndof_q = el_nnodes ;

//extend to 3D the dof values
  ExtendDofs();  
  
//set to zero  the gauss values
  for (uint i=0; i<3; i++) _curl_g3D[i]=0.;
  
//compute the gauss values
           for (uint eln=0; eln<el_nnodes; eln++)    { 

          //curl 
          for (uint idim=0; idim<3; idim++) {
	   uint idimp1 = (idim+1)%3;
	   uint idimp2 = (idim+2)%3;
// 	  _curl_g3D[idim] += 
// 	  (_currGP._dphidxyz_ndsQLVB_g3D[ord][eln+idimp1*el_nnodes] * _val_dofs3D[eln+idimp2*el_ndof_q] 
// 	 - _currGP._dphidxyz_ndsQLVB_g3D[ord][eln+idimp2*el_nnodes] * _val_dofs3D[eln+idimp1*el_ndof_q]); 
           }
	    //end curl

	   }


return;
 }
 
 
 const CurrentElem<double> &  CurrentQuantity::GetCurrentElem() const { 
      return _currEl;
 } 

///////////////
//here I should check that the    _grad_g double array has been allocated
//unfortunately there is no way to check this in C

 void CurrentQuantity::grad_g() {
   
        const uint ndim = _currEl.GetDim();
   const uint   el_ndof = _ndof;
   const uint     nvars = _dim;
   const uint  fe_order = _FEord;

   
   //set to zero
 for (uint ivar=0; ivar<nvars; ivar++) { for (uint idim=0; idim<ndim; idim++) {
                    _grad_g[ivar][idim]=0.;  }  }

//compute
       for (uint eln=0; eln<el_ndof; eln++)    { 
	  
          for (uint ivar=0; ivar<nvars; ivar++) {

	    const uint indxvar=eln+ivar*el_ndof;

	    for (uint idim=0; idim<ndim; idim++) {

	  const uint indxdim=eln+idim*el_ndof;

//               _grad_g[ivar][idim] += _currGP._dphidxyz_ndsQLVB_g[fe_order][indxdim]*_val_dofs[indxvar];
	       
	    }
	    }

      }
   
   return;
   
 }
  
 
//=================================================================== 
 void CurrentQuantity::val_g() {

const uint el_ndof = _ndof;
const uint nvars = _dim;
const uint FEord = _FEord;
//set to zero
 for (uint ivar=0; ivar<nvars; ivar++) _val_g[ivar]=0.;

//compute
       for (uint eln=0; eln<el_ndof; eln++)    { 
	  
          for (uint ivar=0; ivar<nvars; ivar++) {
	  
	       const uint indx=eln+ivar*el_ndof;
	  
//               _val_g[ivar] += _currGP._phi_ndsQLVB_g[FEord][eln]*_val_dofs[indx];
	       
	    }

      }
	 	 
 return;
 
 }

//nota che questa routine e' solo per i vect CON QUANTITY!!!
void CurrentQuantity::VectWithQtyFillBasic() {
  
  if ( _qtyptr == NULL ) {std::cout << " Vect must be endowed with a quantity for this to work " << std::endl; abort();}
    _eqnptr   = _qtyptr->_eqn; 
    _dim      = _qtyptr->_dim;
    _FEord    = _qtyptr->_FEord;
    _ndof     = _currEl.GetElemType(_FEord)->GetNDofs();
    
    
//   if (_SolName == "") { std::cout << " I need to have a solution " << std::endl; abort(); }
  
  
  

    return;
}

void CurrentQuantity::Allocate() {
  
    _val_g.resize(_dim);   
    _val_g3D.resize(3); 
    
    _val_dofs.resize(_dim*_ndof);
    _val_dofs3D.resize(3*_ndof);
    
    _grad_g.resize(_dim);
  for (uint i=0; i< _dim;i++) { _grad_g[i].resize(_currEl._mesh.get_dim()); }

    _grad_g3D.resize(_dim);
  for (uint i=0; i< _dim;i++) { _grad_g3D[i].resize(3); }

    _curl_g3D.resize(3);

    _el_average.resize(_dim);

  return; 
}


///copy the space_dim-sized dof vector into its 3D version
void CurrentQuantity::ExtendDofs() {
  
  //AAA: valid from ndim to 3

  const uint ndim = _currEl._mesh.get_dim();
  const uint el_ndofs = _ndof;
  //set to zero
  for (uint eln=0; eln<el_ndofs; eln++)  {
    for (uint i=0; i<3; i++) {
       _val_dofs3D[eln+i*el_ndofs] = 0.;
    }
  }
//extend
  for (uint eln=0; eln<el_ndofs; eln++)    {
    for (uint idim=0; idim<ndim; idim++) {
      _val_dofs3D[eln+idim*el_ndofs] = _val_dofs[eln+idim*el_ndofs];
    }
  }

  
  
  return;
}



//TODO we must change a little bit the way we do things
//First, when we pick the connectivity from the mesh, of course we know
// that if we have constant elements then we may not need the connectivity actually,
// but only the element number
//So, for every element we should pick the array of 
//"ALL THE GEOMETRICAL ENTITIES FROM WHICH WE EXTRACT ALL THE DOFS OF THE CURRENT ELEMENT"
//The numbering of the elements is simply based on the order of the connectivities.
//at present, we do not want to do a unique array, but we'll simply use 
// "el_conn" for the NODES
// "iel"     for the ELEMENTS

//Wait a minute... here is an issue. 
//We pick the dofs at the fine level, because it is our reference,
// it is the one that contains the VECTOR and not RESIDUAL, ERRORS, or whatever...
//but, together with that, we should pick things from the FINE ELEMENTS,
//but if we are in a coarse value, what element will we pick?
// We will pick the CHILDREN of the CURRENT ELEMENT, but we do not have it now...
//and WHICH CHILD?!?
// Or, if we go on _x_old[Level], will it be the correct VALUE, or is it the ERROR?
//we'll see how to do that, TODO TODO TODO, so far we go ahead
//so, instead of "iel" we should have "some child of iel".

// TODO we use the FINE LEVEL as a reference BOTH for the INITIAL CONDITION 
// and for the BOUNDARY CONDITION
// The point is that when we have elements and we are not on the FINE LEVEL,
// How do we know what to pick from the x_old FINE vector?
// I mean, PICKING is different from PRINTING:
//when you print you don't make any particular distinction,
// you just loop over ALL the ELEMENTS and you dont care about a PARTICULAR element

//this routine retrieves the dofs for the associated Quantity to this vector
//this quantity has an associated equation, which is the current one
//if the quantity has an associated equation, this is the routine, BELONGING TO THE EQUATION
// otherwise there is a routine BELONGING TO THE QUANTITY

//TODO wait a second: are we always sure that in our NONLINEAR ITERATIONS, or MULTIGRID ITERATIONS, or whatever,
//we ALWAYS have to pick the TRUE VALUES of a QUANTITY, or sometimes in the multigrid or nonlinear 
//algorithms we have to pick the ERRORS of that Quantity ?
//For example in the Full Multigrid Cycle...
// we need to give a meaning to the vectors we use in the multigrid cycles
//So far we have in the FINE x_old the solution of the problem,
// in all the others we have the ERROR epsilon.
// In the same way in the FINE RHS we have the true rhs,
// while in all the other rhs we have the RESIDUALS.

//TODO this routine must be fixed a lot
void CurrentQuantity::GetElemDofs()  {
  
  //we should put some try catch or something, to make sure that what we are calling here is already correctly filled as it should be
  //TODO FROM EQUATION HERE
  //TODO here the equation pointer is for sure different from zero, otherwise we cannot do anything
    if ( _eqnptr == NULL ) {std::cout << " We need the Equation here" << std::endl; abort();}

  
  const uint Lev_pick_dof = _eqnptr->GetGridn()-1;  //we use the FINE Level as reference //TODO wait, this is a mistake!!!
  // of course you want to take from the FINE LEVEL, but ONLY FOR THE NODES!!! FOR the ELEMENTS you have to take from EACH LEVEL!!!
  // or, you still have to take from the FINE provided that you give a map from element to fine!!!
  
  // AAAAAAAAA this is crucial here!!!
  //Ok, this routine has a problem. the point is when we do the multigrid cycle.
  //As we know, the multigrid cycle solves for the errors on all the grid, and eventually 
  //the only solution is the one we have at the FINE LEVEL, so only the FINE VECTOR is the vector from which to pick the dofs...
  //so  what if we are at some coarser level? what is the dof at the coarser level?
  //we have to take it from its FINEST CHILDS. These maybe 4,16,32, depending on the distance between current level and finest level
  
  
  
  
  const uint vect_ord = _FEord;
  int length_nodedof [QL];  
  length_nodedof[QQ] = _currEl._mesh._NoNodesXLev[_eqnptr->GetGridn() - 1];
  length_nodedof[LL] = _currEl._mesh._NoNodesXLev[_eqnptr->GetGridn() - 1];
  length_nodedof[KK] = _currEl._mesh._n_elements_vb_lev[VV][_currEl.GetLevel()];

   int off_total = 0;
   for (uint i = 0; i < _qtyptr->_pos; i++) off_total += _eqnptr->GetUnknownQuantitiesVector()[i]->_dim * _eqnptr->_dofmap._DofNumLevFE[ _currEl.GetLevel() ][ _eqnptr->GetUnknownQuantitiesVector()[i]->_FEord ];

   int DofObj = 0;

   for (uint ivar=0; ivar < _dim; ivar++)    {

         for (uint d = 0; d <  _ndof; d++)    {
               const uint     indx  = d + ivar * _ndof;

	          DofObj = _currEl.GetConn()[d];
	       
                 const uint dofkivar = _eqnptr->_dofmap.GetDofPosIn(Lev_pick_dof,DofObj + ivar*length_nodedof[vect_ord] + off_total);

	         _val_dofs[indx] =  ( *( _eqnptr->_LinSolver[Lev_pick_dof]->_EPSC ) )(dofkivar); 

	         }
	  }
  
 return; 
}


//this function fills the element dofs from a vector of CurrentQuantity objects
   void  CurrentQuantity::GetElemDofs(const std::vector<CurrentQuantity*> vec_in)  {
     
        for (uint ivar=0; ivar < _dim; ivar++)    {
              for (uint d = 0; d <  _ndof; d++)    {
               const uint     indx  = d + ivar * _ndof;
		_val_dofs[indx] =  vec_in[ivar]->_val_dofs[d];
		
	      }
	}
     
    return; 
   }



// void CurrentQuantity::GetElemDofs_two()  {
// 
//       unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex(_qtyptr->_name));  
//       unsigned nve2       = _currEl._mesh_new->el->GetElementDofNumber(kel,SolType2);
// 
//   
//       for (unsigned i=0;i<nve2;i++) {
// 	unsigned inode=_currEl._mesh_new->el->GetMeshDof(kel,i,SolType2);
// 	unsigned inode_Metis=_currEl._mesh_new->GetSolutionDof(inode,SolType2);
// 	for(int j=0; j<_currEl._mesh_new->GetDimension(); j++) {

// 	  Soli[indexVAR[j]][i] =  (*mysolution->_Sol[indVAR[j]])(inode_Metis);
// // 	  dofsVAR[j][i] = myLinEqSolver->GetSystemDof(indVAR[j],indexVAR[j],inode); 
// 	}
//       }




// clearly _el_average must be allocated already!!!

  void CurrentQuantity::SetElemAverage() {

       for (uint idim=0; idim< _dim; idim++)  _el_average[idim]=0.;

    for (uint idim=0; idim< _dim; idim++) {
       for (uint eln=0; eln< _ndof; eln++)    { 
        const uint indxn = eln+idim*_ndof;
	     _el_average[idim]   +=  _val_dofs[indxn];
       }
     }

  for (uint idim=0; idim< _dim; idim++)   _el_average[idim]= _el_average[idim]/_ndof;
  
   return; 
  }
  


} //end namespace femus


