/*=========================================================================

 Program: FEMUS
 Module: CurrentElem
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "CurrentElem.hpp"
#include "SystemTwo.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "Domain.hpp"
#include "ElemType.hpp"
#include "DofMap.hpp"

#include "CurrentQuantity.hpp"


namespace femus {

template < typename real_num_mov >
    CurrentElem<real_num_mov>::CurrentElem(const uint iel_in, const uint iproc_in, const uint level, const uint vb, const SystemTwo * eqn_in, const MultiLevelMeshTwo& mesh, const std::vector< std::vector<const elem_type*> >  & elem_type_in, const Mesh * mesh_new ):
    _eqn(eqn_in),
    _mesh(mesh),
    _dim(_mesh.get_dim()-vb),
    _elem_type(elem_type_in[mesh.get_dim()-vb -1]),
    _mesh_vb(vb),
    _Level(level),
    _iel(iel_in),
    _proc(iproc_in),
    _mesh_new(mesh_new)
    {
    
//========== Current "Geometric Element"  ========================
  uint elnodes = NVE[ _mesh._geomelem_flag[_dim-1] ][BIQUADR_FE];
  _el_conn.resize(elnodes);
  _el_conn_new.resize(elnodes);   
   _xx_nds.resize(_mesh.get_dim()*elnodes);
//========== Current "Geometric Element"  ========================

//========== Current "Equation Element"  ========================
  _el_n_dofs = 0;
     for (int fe = 0; fe < QL; fe++) {  _el_n_dofs += (_elem_type[fe]->GetNDofs() )*_eqn->_dofmap._nvars[fe]; }

  _el_dof_indices.resize(_el_n_dofs);
  _bc_eldofs.resize(_el_n_dofs);
  _KeM.resize(_el_n_dofs,_el_n_dofs);
  _FeM.resize(_el_n_dofs);
//========== Current "Equation Element"  ========================
      
      
  }
 

    
//Ok, in the new logic  the element and gauss points are going to be instantiations
// of the assemble routine. therefore they will be filled using the this-> pointer,
// because they need to access the data of the equation and turn those into data
//for each current element or current gauss point.



template < typename real_num_mov >
void CurrentElem<real_num_mov>::SetElDofsBc()  {

/*CHECK*/   if (_vol_iel_DofObj >= _mesh._n_elements_vb_lev[VV][_Level] ) { std::cout << "Out of the node_dof map FE KK range" << std::endl; abort();}

  const uint Lev_pick_bc_dof = _mesh._NoLevels -1;  //we use the FINE Level as reference
  
int off_local_el[QL];
off_local_el[QQ] = 0;
off_local_el[LL] = _eqn->_dofmap._nvars[QQ]*(_elem_type[QQ]->GetNDofs() );
off_local_el[KK] = _eqn->_dofmap._nvars[QQ]*(_elem_type[QQ]->GetNDofs() ) + _eqn->_dofmap._nvars[LL]*(_elem_type[LL]->GetNDofs() );
  

 int DofObj = 0;
for (int fe=0; fe < QL; fe++) {
for (uint ivar=0; ivar < _eqn->_dofmap._nvars[fe]; ivar++)    {
      for (uint d=0; d< _elem_type[fe]->GetNDofs(); d++)    {
	
	     if (fe < KK )       DofObj =        _el_conn[d];
	     else if (fe == KK)  DofObj = _vol_iel_DofObj;
	     
          const uint     indx  = d + ivar*_elem_type[fe]->GetNDofs() + off_local_el[fe];
	  _el_dof_indices[indx] = _eqn->_dofmap.GetDof(_Level,fe,ivar,DofObj);

	   const uint dofkivar = _eqn->_dofmap.GetDof(Lev_pick_bc_dof,fe,ivar,DofObj); 
             _bc_eldofs[indx] = _eqn->_bcond._bc[dofkivar]; 
	   
	 }
    } 
} // end fe
    
      return;

}


// ========================================================

/*
   HEX27:      7              18             6     			      
               o--------------o--------------o     			      
              /:             /              /|     			      
             / :            /              / |     			      
            /  :           /              /  |     			      
         19/   :        25/            17/   |     			      
          o--------------o--------------o    |     			      
         /     :        /              /|    |     			      
        /    15o       /    23o       / |  14o     
       /       :      /              /  |   /|     
     4/        :   16/             5/   |  / |     
     o--------------o--------------o    | /  |     			      
     |         :    |   26         |    |/   |     
     |  24o    :    |    o         |  22o    |     
     |         :    |       10     |   /|    |                                
     |        3o....|.........o....|../.|....o     
     |        .     |              | /  |   / 2      
     |       .    21|            13|/   |  /        
  12 o--------------o--------------o    | /         
     |     .        |              |    |/    zeta  
     |  11o         | 20o          |    o       ^    
     |   .          |              |   / 9      |   eta
     |  .           |              |  /         |  /
     | .            |              | /          | /
     |.             |              |/           |/ 
     o--------------o--------------o            -------> xi    
     0              8              1                
*/

///This function prints the element orientation
///A QUADRATIC MESH VECTOR is passed HERE  
///It is for debugging purposes
template < typename real_num_mov >
void CurrentElem<real_num_mov>::PrintOrientation() const {
  
      const uint mesh_dim = _mesh.get_dim();
      const uint el_nnodes   = _el_conn.size();

       std::vector<double>   xi(mesh_dim,0.);
       std::vector<double>  eta(mesh_dim,0.);
       std::vector<double> zeta(mesh_dim,0.);  //TODO this should only be instantiated  in the 3D case, in any case we avoid USING IT 

      for (uint idim=0; idim< mesh_dim; idim++) {
                              xi[idim] = _xx_nds[1+idim*el_nnodes]-_xx_nds[0+idim*el_nnodes]; //1-0 xi axis
                             eta[idim] = _xx_nds[3+idim*el_nnodes]-_xx_nds[0+idim*el_nnodes]; //3-0 eta axis
        if ( mesh_dim == 3 )   zeta[idim] = _xx_nds[4+idim*el_nnodes]-_xx_nds[0+idim*el_nnodes]; //4-0 zeta axis

     }

     std::cout << "asse xi   ";
       for (uint idim=0; idim< mesh_dim; idim++)  std::cout << xi[idim] << " ";
      std::cout <<std::endl;
     std::cout << "asse eta  ";
       for (uint idim=0; idim< mesh_dim; idim++)  std::cout << eta[idim] << " ";
      std::cout <<std::endl;
      if ( mesh_dim == 3 )  {
     std::cout << "asse zeta ";
       for (uint idim=0; idim< mesh_dim; idim++)  std::cout << zeta[idim] << " ";
      std::cout <<std::endl;
      }

      
      
      return; 
    }

 

 
 
// ========================================================
//   void CurrentElem<real_num_mov>::SetDofobjConnCoords(const uint Level,const uint isubd_in,const uint iel,
// 				uint el_conn[], double xx[]) const {///get the global node numbers for that element and their coordinates
///this routine does not yield the connectivity
//is this function called when the class members are already filled?
//i guess so
//oh maybe here's a reason why public protected and such things exist!
//only the class can call its members, and the class calls them 
//only when  she knows that she can...
//-------> PAY A LOT OF ATTENTION!!!
///this routine is "PARALLEL", in the sense that 
///it works only for the _iproc in which you are!
///pay attention to every routine that has _iproc!!!
///if you want to do a loop over ALL PROCESSES
/// you do not have to use iproc but you have to
/// PASS the SUBDOMAIN EXPLICITLY!!!

//TODO now we will adjust this routine so that it gives us, for the current level and the current processor,
// the vector of DOF OBJECTS that are basically the POSITIONS to be input in the _node_dof map.
//We should distinguish between NODE DOFObjects and ELEM DofObjects.,
// or DOFOBJECTS for QQ, LL and KK FE Families.
//well, here we are at the MESH LEVEL, so we will provide NODE and ELEM dofobjects.
// Then, each finite element family knows if she has to pick the NODES (AND WHICH NODES) or the ELEMENT
// We assume that the highest possible degree is QUADRATIC
// So, for instance, QUAD9 takes 9 NODES, QUAD4 takes 4 NODES, QUAD1 takes 1 ELEMENT, QUAD8 takes 8 NODES.
// The mesh is there to provide you with the list of ALL DOF OBJECTS for that element.




   // =====================================================================================
template < typename real_num_mov >
  void CurrentElem<real_num_mov>::SetDofobjConnCoords() {

    const uint mydim = _mesh.get_dim();
    const uint el_nnodes   = _el_conn.size();
          
   for (uint n=0; n<el_nnodes; n++)    {

     _el_conn[n] = _mesh._el_map[_mesh_vb][( _iel + _mesh._off_el[_mesh_vb][_mesh._NoLevels*_proc + _Level] )*el_nnodes+n];

      for (uint idim=0; idim < mydim; idim++) {
        const uint indxn = n+idim*el_nnodes;
          _xx_nds[indxn] = _mesh._xyz[_el_conn[n]+idim*_mesh._NoNodesXLev[_mesh._NoLevels-1]];
      }
   }
   
   
    int sum_elems_prev_sd_at_lev = 0;
      for (uint pr = 0; pr< _proc; pr++) { sum_elems_prev_sd_at_lev += _mesh._off_el[_mesh_vb][_mesh._NoLevels*pr + _Level + 1] - _mesh._off_el[_mesh_vb][ _mesh._NoLevels*pr + _Level]; }
    uint iel_DofObj = _iel + sum_elems_prev_sd_at_lev;
    if       (_mesh_vb == VV)  { _vol_iel_DofObj = iel_DofObj; }   
    else if  (_mesh_vb == BB)  { _vol_iel_DofObj = _mesh._el_bdry_to_vol[_Level][iel_DofObj]; }  
   

   
    return;
  }   

  
  
// ========================================================
  ///here, the original element order of the nodes is respected in xyz._val_dofs 
  //in practice this sets the dofs of xyz
  //maybe it belongs more to the Vect class
  //so far we leave it here, and we say that the current element delivers 
  //its data to the xyz Vect
  
  
template < typename real_num_mov >
void CurrentElem<real_num_mov>::ConvertElemCoordsToMappingOrd(CurrentQuantity& myvect) const {

  
  const uint  elndof = myvect._ndof;
  const uint vectdim = myvect._dim;
  const uint offset = _el_conn.size();
 
 //TODO ASSERT
 /* assert(*/ if (elndof > offset) { std::cout << "Quadratic transformation over linear mesh " << std::endl; abort(); }  /*);*/
  
  //information for passing from mesh to dofs
  
    for (uint d=0; d < elndof; d++) {
	                                        
      for (uint idim=0; idim < vectdim; idim++) {
          const uint     indxq  =    d + idim*elndof;

	  myvect._val_dofs[indxq] = _xx_nds[d+idim*offset];
      }
    }

  return;

 }   
   
//=================================
/// @deprecated
//This function is for NS type equations:
//it computes the flags for pressure and stress integrals 
//based on the pressure nodes
template < typename real_num_mov >
 int CurrentElem<real_num_mov>::Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(const uint ndof_in, const uint space_dim, const CurrentQuantity& press_in) const {
   
   int press_fl = 0;
   
	const uint el_ndof_p =  press_in._ndof;
	const uint el_ndof_u =  ndof_in;
	const uint   nvars_u =  space_dim;
        int press_sum=0;
           for (uint i=0; i< el_ndof_p; i++)   press_sum += _bc_eldofs[nvars_u*el_ndof_u + i]; //only one linear variable... pay attention when trying linear-linear
            if ( press_sum ==0 )                 {  press_fl = 1; }   //ALL zeros: ONLY PRESSURE

  
    return press_fl;
  }

  
// // // // ========================================================
// // // //This function transforms the node coordinates into the reference node coordinates
// // //   void CurrentElem<real_num_mov>::TransformElemNodesToRef(Domain* mydom, double* refbox_xyz) {
// // //    
// // //    std::vector<double>   x_in(_dim);
// // //    std::vector<double>   x_out(_dim);
// // //   const uint el_nds = _xx_nds.size()/_mesh.get_dim();
// // // 
// // //       for ( uint n=0; n < el_nds; n++ ) {
// // // 	
// // //    for ( uint idim=0; idim < _dim; idim++ )  x_in[idim] = _xx_nds[n + idim*el_nds];
// // //   
// // //   mydom->TransformPointToRef(&x_in[0],&x_out[0]);
// // // 
// // //    for ( uint idim=0; idim < _dim; idim++ )  refbox_xyz[n + idim*el_nds] = x_out[idim];
// // //    
// // //       }
// // //    
// // //   return; 
// // //  }

 
    // =====================================================================================
//   void CurrentElem<real_num_mov>::SetDofobjConn_twoCoordsDof() {

// for each scalar variable involved, you need to define
  // dofsVAR    # of solutions
  // indexVAR   # of solutions
  
  //indVAR   # of rows
  //SolType  # of columns
  
  //Now the thing is: we can have EQUATION BLOCKS that are NOT ASSOCIATED with any SOLUTION
  //So, not only you can have a solution without equation, you can also have an equation without solution

  //What is the difference between GetIndex and GetSolPdeIndex ??
  //  AX AY AZ are solutions or not? Wait... they only appear in the switch for SetBoundaryConditionComsol...
  
//       unsigned nel    = mymsh->GetNumberOfElements();

//          elem * myel = _mesh_new->el;
//         unsigned kel    = _mesh_new->IS_Mts2Gmt_elem[_iel]; 
//     short unsigned kelt = myel->GetElementType(kel);


//   unsigned order_ind = ml_sol->GetSolutionType(ml_sol->GetIndex("U")); 
//     unsigned nve        = myel->GetElementDofNumber(kel,order_ind);
//     
// 
// 
//     _el_conn_new.resize(nve);
//     
//     for (unsigned i=0;i<nve;i++) {
//       // gambit nodes
//            unsigned inode=(order_ind<3)?(myel->GetElementDofIndex(kel,i)-1u):(kel+i*nel);
//       unsigned inode=myel->GetElementDofIndex(kel,i)-1u;

//       // dof metis
//       /*metis_node2*/_el_conn_new[i] = mymsh->GetSolutionDof(inode,BIQUADR_FE);
//                       metis_node1[i] = mymsh->GetSolutionDof(inode,SolType[2*dim]);


//  	dofsVAR[j+dim][i]= mylsyspde->GetSystemDof(indVAR[j+dim],indexVAR[j+dim],inode);   

//     }
//     
//     
//     
//     return;
//   }   

 
    // =====================================================================================
//   void CurrentElem<real_num_mov>::SetCoords_two() {
// 
//      for (unsigned i=0;i<nve;i++) {
//       gambit nodes
//       unsigned inode=myel->GetElementDofIndex(kel,i)-1u;
//       dof metis
//       unsigned inode_Metis=mymsh->GetSolutionDof(inode,2);
//       metis_node2[i]=inode_Metis;
//       
//       unsigned inode_Metis=mymsh->GetSolutionDof(inode,2);
//       flag to know if the node "inode" lays on the fluid-solid interface
//       solidmark[i]=myel->GetNodeRegion(inode); // to check
//       for(int j=0; j<dim; j++) {
// 	Updated coordinates (Moving frame)
//         vx[j][i]= (*mymsh->_topology->_Sol[j])(inode_Metis) + (*mysolution->_Sol[indVAR[j]])(inode_Metis);
// 	Old coordinates (Moving frame)
//         vx_old[j][i]= (*mymsh->_topology->_Sol[j])(inode_Metis) + (*mysolution->_SolOld[indVAR[j]])(inode_Metis);
// 	Fixed coordinates (Reference frame)
// 	vx_hat[j][i]= (*mymsh->_topology->_Sol[j])(inode_Metis);  
// 	displacement dofs
// 	dofsVAR[j][i]= mylsyspde->GetSystemDof(indVAR[j],indexVAR[j],inode); 
// 	velocity dofs
// 	dofsVAR[j+dim][i]= mylsyspde->GetSystemDof(indVAR[j+dim],indexVAR[j+dim],inode);   
//       }
//     }
//  
//     
//     
//     
//     return;
//   }   

template class CurrentElem<double>;
 
} //end namespace femus


  
