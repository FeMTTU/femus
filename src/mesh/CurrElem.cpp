#include "CurrElem.hpp"


#include "EqnBase.hpp"
#include "EquationsMap.hpp"
#include "Utils.hpp"
#include "MeshTwo.hpp"
#include "GeomEl.hpp"

// #include "FEMap.h"
#include "FEElemBase.hpp"

#include "QuantityLocal.hpp"

    CurrElem::CurrElem(EqnBase & eqn_in, EquationsMap & e_map_in ):
    _eqn(eqn_in),
    _eqnmap(e_map_in)
    {
    
//========== ELEMENT: Current Geometric Element (SERVICE)  ========================
   const uint mesh_ord = (int) _eqnmap._mesh._mesh_rtmap.get("mesh_ord");
  _el_conn   = new uint*[VB];
  _xx_nds    = new double*[VB];
  _el_xm     = new double*[VB];
 for (int vb = 0; vb < VB; vb++) {
  uint elnodes = _eqnmap._mesh._GeomEl._elnds[vb][mesh_ord];     //TODO the mesh is quadratic
  _el_conn[vb] = new uint[ elnodes ];   
   _xx_nds[vb] = new double[_eqnmap._mesh._dim*elnodes ];
    _el_xm[vb] = new double[_eqnmap._mesh._dim];
 }  
//========== ELEMENT: Current Geometric Element (SERVICE)  ========================

for (int vb = 0; vb < VB; vb++) {
  _el_n_dofs[vb] = 0;
     for (int fe = 0; fe < QL; fe++) {  _el_n_dofs[vb] += (_eqnmap._AbstractFE[fe]->_ndof[vb])*_eqn._nvars[fe]; }    //EQUATION-RELATED
 }

for (int vb = 0; vb < VB; vb++) {
  _el_dof_indices[vb].resize(_el_n_dofs[vb]);                  //EQUATION-RELATED
  _bc_eldofs[vb] = new uint[_el_n_dofs[vb]];                   //EQUATION-RELATED
  _KeM[vb].resize(_el_n_dofs[vb],_el_n_dofs[vb]);              //EQUATION-RELATED
  _FeM[vb].resize(_el_n_dofs[vb]);                             //EQUATION-RELATED
}

      
      
  }
  
  

    CurrElem::~CurrElem() {
      
      
 for (int vb = 0; vb < VB; vb++) {
  delete [] _el_conn[vb];
  delete []  _xx_nds[vb];
  delete []   _el_xm[vb];  
}

 delete [] _el_conn;
 delete [] _xx_nds;
 delete []  _el_xm;     
      
      
      
      
  }
    
    
    
//Ok, in the new logic  the element and gauss points are going to be instantiations
// of the assemble routine. therefore they will be filled using the this-> pointer,
// because they need to access the data of the equation and turn those into data
//for each current element or current gauss point.



void CurrElem::GetElDofsBc(const uint vbfl, const uint Level)  {

/*CHECK*/   if (_vol_iel_DofObj[vbfl] >= _eqnmap._mesh._NoElements[VV][Level] ) { std::cout << "Out of the node_dof map FE KK range" << std::endl; abort();}

  const uint Lev_pick_bc_dof = _eqnmap._mesh._NoLevels -1;  //we use the FINE Level as reference
  
int off_local_el[QL];
off_local_el[QQ] = 0;
off_local_el[LL] = _eqn._nvars[QQ]*(_eqnmap._AbstractFE[QQ]->_ndof[vbfl]);
off_local_el[KK] = _eqn._nvars[QQ]*(_eqnmap._AbstractFE[QQ]->_ndof[vbfl]) + _eqn._nvars[LL]*(_eqnmap._AbstractFE[LL]->_ndof[vbfl]);
  

 int DofObj = 0;
for (int fe=0; fe < QL; fe++) {
for (uint ivar=0; ivar < _eqn._nvars[fe]; ivar++)    {
      for (uint d=0; d< _eqnmap._AbstractFE[fe]->_ndof[vbfl]; d++)    {
	
	     if (fe < KK )       DofObj =        _el_conn[vbfl][d];
	     else if (fe == KK)  DofObj = _vol_iel_DofObj[vbfl];
	     
          const uint     indx  = d + ivar*_eqnmap._AbstractFE[fe]->_ndof[vbfl] + off_local_el[fe];
	  _el_dof_indices[vbfl][indx] = _eqn._node_dof[Level][ DofObj + ivar*_eqn._DofNumLevFE[Level][fe] + _eqn._DofOffLevFE[Level][fe] ]; 

         if (fe < KK ) { const uint dofkivar = _eqn._node_dof[Lev_pick_bc_dof][ DofObj + ivar*_eqn._DofNumLevFE[Lev_pick_bc_dof][fe] + _eqn._DofOffLevFE[Lev_pick_bc_dof][fe] ]; 
             _bc_eldofs[vbfl][indx] = _eqn._bc[dofkivar]; }
         else if (fe == KK)    _bc_eldofs[vbfl][indx] = _eqn._bc_fe_kk[Level][ DofObj + ivar*_eqn._DofNumLevFE[Level][KK] ];
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
void CurrElem::get_el_orient(const uint vb) const {
  
      const uint mesh_ord = (int) _eqnmap._mesh._mesh_rtmap.get("mesh_ord");
      const uint el_nnodes   = _eqnmap._mesh._GeomEl._elnds[vb][mesh_ord];
      const uint mydim = _eqnmap._mesh._dim;

       std::vector<double>   xi(mydim,0.);
       std::vector<double>  eta(mydim,0.);
       std::vector<double> zeta(mydim,0.);  //TODO this should only be instantiated  in the 3D case, in any case we avoid USING IT 

      for (uint idim=0; idim< mydim; idim++) {
                              xi[idim] = _xx_nds[vb][1+idim*el_nnodes]-_xx_nds[vb][0+idim*el_nnodes]; //1-0 xi axis
                             eta[idim] = _xx_nds[vb][3+idim*el_nnodes]-_xx_nds[vb][0+idim*el_nnodes]; //3-0 eta axis
        if ( mydim == 3 )   zeta[idim] = _xx_nds[vb][4+idim*el_nnodes]-_xx_nds[vb][0+idim*el_nnodes]; //4-0 zeta axis

     }

     std::cout << "asse xi   ";
       for (uint idim=0; idim< mydim; idim++)  std::cout << xi[idim] << " ";
      std::cout <<std::endl;
     std::cout << "asse eta  ";
       for (uint idim=0; idim< mydim; idim++)  std::cout << eta[idim] << " ";
      std::cout <<std::endl;
      if ( mydim == 3 )  {
     std::cout << "asse zeta ";
       for (uint idim=0; idim< mydim; idim++)  std::cout << zeta[idim] << " ";
      std::cout <<std::endl;
      }

      
      
      return; 
    }

 

 
 
// ========================================================
//   void CurrElem::get_el_nod_conn_lev_subd(const uint vb,const uint Level,const uint isubd_in,const uint iel,
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


// ========================================================
///Compute the element center
///TODO must do this in the QUADRATIC case

  void CurrElem::get_el_ctr(const uint vb) const {

    const uint mesh_ord = (int) _eqnmap._mesh._mesh_rtmap.get("mesh_ord");    
    const uint el_nnodes   = _eqnmap._mesh._GeomEl._elnds[vb][mesh_ord];
    const uint mydim = _eqnmap._mesh._dim;

       for (uint idim=0; idim< mydim; idim++)  _el_xm[vb][idim]=0.;

    for (uint idim=0; idim< mydim; idim++) {
       for (uint eln=0; eln<el_nnodes; eln++)    { 
        const uint indxn = eln+idim*el_nnodes;
	     _el_xm[vb][idim]   +=   _xx_nds[vb][indxn];
       }
     }

  for (uint idim=0; idim< mydim; idim++)   _el_xm[vb][idim]= _el_xm[vb][idim]/el_nnodes;
  
   return; 
  }

   // =====================================================================================
  void CurrElem::get_el_nod_conn_lev_subd(const uint vb,const uint Level,const uint isubd_in,const uint iel) const {

    const uint mesh_ord = (int) _eqnmap._mesh._mesh_rtmap.get("mesh_ord");    
    const uint el_nnodes   = _eqnmap._mesh._GeomEl._elnds[vb][mesh_ord];
    const uint mydim = _eqnmap._mesh._dim;
          
   for (uint n=0; n<el_nnodes; n++)    {

     _el_conn[vb][n] = _eqnmap._mesh._el_map[vb][( iel + _eqnmap._mesh._off_el[vb][_eqnmap._mesh._NoLevels*isubd_in + Level] )*el_nnodes+n];

      for (uint idim=0; idim < mydim; idim++) {
        const uint indxn = n+idim*el_nnodes;
          _xx_nds[vb][indxn] = _eqnmap._mesh._xyz[_el_conn[vb][n]+idim*_eqnmap._mesh._NoNodesXLev[_eqnmap._mesh._NoLevels-1]];
      }
   }
    
    return;
  }   

 // =====================================================================================
void CurrElem::get_el_DofObj_lev_subd(const uint vb,const uint Level,const uint isubd_in,const uint iel) {

     int sum_elems_prev_sd_at_lev = 0;
      for (uint pr = 0; pr< isubd_in; pr++) { sum_elems_prev_sd_at_lev += _eqnmap._mesh._off_el[vb][_eqnmap._mesh._NoLevels*pr + Level + 1] - _eqnmap._mesh._off_el[vb][ _eqnmap._mesh._NoLevels*pr + Level]; }
    uint iel_DofObj = iel + sum_elems_prev_sd_at_lev;
    if       (vb == VV)  { _vol_iel_DofObj[vb] = iel_DofObj; }   
    else if  (vb == BB)  { _vol_iel_DofObj[vb] = _eqnmap._mesh._el_bdry_to_vol[Level][iel_DofObj]; }
          
    return; 
}

   
  
// ========================================================
  ///here, the original element order of the nodes is respected in xyz._val_dofs 
  //in practice this sets the dofs of xyz
  //maybe it belongs more to the Vect class
  //so far we leave it here, and we say that the current element delivers 
  //its data to the xyz Vect
  
  
void CurrElem::ConvertElemCoordsToMappingOrd(const uint vb,QuantityLocal& myvect) const {

  
  const uint  elndof = myvect._ndof[vb];
  const uint vectdim = myvect._dim;
  const uint mesh_ord = (int) _eqnmap._mesh._mesh_rtmap.get("mesh_ord");    
  const uint offset = _eqnmap._mesh._GeomEl._elnds[vb][mesh_ord];
 
 //TODO ASSERT
 /* assert(*/ if (elndof > offset) {std::cout << "Quadratic transformation over linear mesh " << std::endl;abort();}  /*);*/
  
  //information for passing from mesh to dofs
  
    for (uint d=0; d < elndof; d++) {
	                                        
      for (uint idim=0; idim < vectdim; idim++) {
          const uint     indxq  =    d + idim*elndof;

	  myvect._val_dofs[indxq] = _xx_nds[vb][d+idim*offset];
      }
    }

  return;

 }   
   
  