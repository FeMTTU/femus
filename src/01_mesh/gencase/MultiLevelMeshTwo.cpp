/*=========================================================================

 Program: FEMUS
 Module: MultiLevelMeshTwo
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "MultiLevelMeshTwo.hpp"

#include <sstream>
#include <vector>
#include <algorithm>

#include "Domain.hpp"
#include "Typedefs.hpp"
#include "XDMFWriter.hpp"
#include "Elem.hpp"
#include "GeomElTypeEnum.hpp"
#include "VBTypeEnum.hpp"
#include "FETypeEnum.hpp"

#include "paral.hpp"


namespace femus {




// ========================================================
MultiLevelMeshTwo::MultiLevelMeshTwo (const unsigned nolevels, const unsigned dim, const GeomElType geomel_type, const std::string mesh_file_in) :
         _dim(dim) {

    _eltype_flag[VV]= geomel_type;
	 
    _mesh_file.assign(mesh_file_in); 

    _geomelem_id.resize(_dim);
    _geomelem_flag.resize(_dim);
    
    if (_eltype_flag[VV] == HEX)  {
      _eltype_flag[BB] = QUAD;
      _geomelem_id[0] = "line";
      _geomelem_id[1] = "quad";
      _geomelem_id[2] = "hex";
      _geomelem_flag[0] = LINE;
      _geomelem_flag[1] = QUAD;
      _geomelem_flag[2] = HEX;
    }
    else if (_eltype_flag[VV] == TET)  {
      _eltype_flag[BB] = TRI;
      _geomelem_id[0] = "line";
      _geomelem_id[1] = "tri";
      _geomelem_id[2] = "tet";
      _geomelem_flag[0] = LINE;
      _geomelem_flag[1] = TRI;
      _geomelem_flag[2] = TET;
    }
    else if (_eltype_flag[VV] == WEDGE)  {
       std::cout << "Wedge not supported" << std::endl; abort(); 
    }
    else if (_eltype_flag[VV] == QUAD)  {
      _eltype_flag[BB] = LINE;
      _geomelem_id[0] = "line";
      _geomelem_id[1] = "quad";
      _geomelem_flag[0] = LINE;
      _geomelem_flag[1] = QUAD;
    }
    else if (_eltype_flag[VV] == TRI)  {
      _eltype_flag[BB] = LINE;
      _geomelem_id[0] = "line";
      _geomelem_id[1] = "tri";
      _geomelem_flag[0] = LINE;
      _geomelem_flag[1] = TRI;
    }
    else if (_eltype_flag[VV] == LINE)  {
      _geomelem_id[0] = "line";
      _geomelem_flag[0] = LINE;
      std::cout << "Geom Elem not supported" << std::endl; abort();
    }
    else  {  std::cout << "Geom Elem not supported" << std::endl; abort();   }
    
if ( _dim == 3  && (geomel_type != HEX && geomel_type != TET && geomel_type != WEDGE) )
        {  std::cout << "Inconsistent input file" << std::endl; abort();   }
    
if ( _dim == 2  && (geomel_type != QUAD && geomel_type != TRI ) )
        {  std::cout << "Inconsistent input file" << std::endl; abort();   }
    
if ( _dim == 1  && (geomel_type != LINE ) )
        {  std::cout << "Inconsistent input file" << std::endl; abort();   }

    
    _iproc    = paral::get_rank();
    _NoSubdom = paral::get_size();   
    _NoLevels = nolevels;
    
    if (MESH_ORDER != QQ) {
        std::cout << "Linear mesh not yet implemented" << std::endl;
        abort();
    }
    if (VB != 2) {
        std::cout << " NoFamFEM is not 2: to do! ";
        abort();
    }

    for (int vb=0;vb < VB; vb++) {
        _elnodes[vb][QQ] = NVE[ _geomelem_flag[_dim-1-vb] ][BIQUADR_FE];
        _elnodes[vb][LL] = NVE[ _geomelem_flag[_dim-1-vb] ][LINEAR_FE];
        _elnodes[vb][KK] = 1;
    }
    //i do not want to use the linear part actually!!

}

// ========================================================
// MultiLevelMeshTwo::~Mesh ()  {
//   
// }

// ========================================================
void MultiLevelMeshTwo::clear ()  {
 
   for (uint imesh =0; imesh<_NoLevels+1; imesh++) {
    delete []_Qnode_lev_Qnode_fine[imesh];
  }
 delete []_Qnode_lev_Qnode_fine;

   for (uint i = 0; i < QL_NODES; i++) {
       delete [] _off_nd[i];
   }

    for (uint imesh =0; imesh < VB; imesh++) {
    delete [] _el_map[imesh];
    delete [] _off_el[imesh];
    delete [] _n_elements_vb_lev[imesh];
  }


   for (uint lev=0; lev < _NoLevels; lev++) delete [] _el_bdry_to_vol[lev]; 
   delete []  _el_bdry_to_vol;
  
    delete[] _n_elements_vb_lev;
    delete[] _off_el;
    delete[] _el_map;
    delete[] _xyz;
    delete[] _NoNodesXLev;
    delete[] _off_nd;
    
  return;
}




// ========================================================
  void MultiLevelMeshTwo::SetDomain(Domain* domain_in)  {
    
    _domain = domain_in;
    
   return; 
  }
  
 // ========================================================
  Domain* MultiLevelMeshTwo::GetDomain() const {
    
   return _domain; 
   
  }










 
  
  










void MultiLevelMeshTwo::FillNodeSto() {

// ================================================
//    NODES (in nd_sto and nod_val)
// ================================================
    _nd_sto = new NodeSto*[_n_nodes];
    for (int i=0;i<_n_nodes;i++) {
        _nd_sto[i]= new NodeSto(100000,_NoLevels);    //initialization values are IMPORTANT!!!
    }

    // Storing the node information with element loop--------------  ==> we SUPERIMPOSE stuff, we fill more than once
    // looping over elements is essential if you want to separate the LINEAR from the QUADRATIC nodes
    for (int i=0;i<_n_elements_sum_levs[VV];i++) {
        for (uint k=0;k<_elnodes[VV][QQ];k++)   {

            int knode=_el_sto[i]->_elnds[k];
            _nd_sto[knode]->_id   = knode;             //libmesh node numbering
            if (_el_sto[i]->_subd > _nd_sto[knode]->_subd)   _nd_sto[knode]->_subd = _el_sto[i]->_subd; // HIGHEST subdomain number among the elements the node belongs to //libmesh does the MINIMUM: by definition, a Node's processor ID is the minimum processor ID for all of the elements which share the node.
            if (_el_sto[i]->_lev  < _nd_sto[knode]->_lev)    _nd_sto[knode]->_lev  = _el_sto[i]->_lev;  // COARSEST level among the elements the node belongs to
            if ( k <  _elnodes[VV][LL] &&
                    _el_sto[i]->_lev < _nd_sto[knode]->_levp)      _nd_sto[knode]->_levp = _el_sto[i]->_lev;  //COARSEST level also for linear
            _nd_sto[knode]->_var  = 1;
        }
    }

    return;

}


//=====================================
void MultiLevelMeshTwo::ComputeElemOffsetsBySubdLevel() {
    _off_el = new int*[VB];
    int sum_el_vb[VB];
    for (int vb=0;vb < VB; vb++) {
        _off_el[vb]=new int[_NoSubdom*_NoLevels+1];
        _off_el[vb][0]=0;
        sum_el_vb[vb]=0;
    }

    for (int vb=0;vb < VB; vb++) {
        for (int  isub=0;isub<_NoSubdom*_NoLevels; isub++)  {
            sum_el_vb[vb]         +=  _n_elements_sl_vb[vb][isub];
              _off_el[vb][isub+1]  =  sum_el_vb[vb];
        }
    }

    return;
}


//====================================
void MultiLevelMeshTwo::ComputeNodeOffsetsBySubdLevel()  {
    _off_nd = new int*[QL_NODES];
    int sum_nd_ql[QL_NODES];
    for (int fe=0;fe < QL_NODES; fe++) {
        _off_nd[fe]=new int[_NoSubdom*_NoLevels + 1];
        _off_nd[fe][0]=0;
        sum_nd_ql[fe]=0;
    }

    for (int fe=0;fe < QL_NODES; fe++) {
        for (int  isub=0; isub<_NoSubdom*_NoLevels; isub++)  {
            sum_nd_ql[fe]         +=  _n_nodes_sl_ql[fe][isub];
            _off_nd[fe][isub+1]  =  sum_nd_ql[fe];
        }
    }

    return;
}






// ==============================
// COMPUTE node map for "Extended Levels"
// ==============================
// this RECEIVES a NODE in FEMUS ordering
// and YIELDS the DOF NUMBERING of ONE SCALAR variable (sometimes QQ, sometimes LL, in MG...)
// ==================================================================
void MultiLevelMeshTwo::ComputeNodeMapExtLevels() {

    int NegativeOneFlag = -1;  
  
    _Qnode_fine_Qnode_lev=new int*[_NoLevels+1];

// quadratic levels
    for (int  ilev=0;ilev< _NoLevels;ilev++) {
        _Qnode_fine_Qnode_lev[ilev] = new int[_n_nodes];
        int count=0;
        for (int  inode=0;inode<_n_nodes;inode++)  _Qnode_fine_Qnode_lev[ilev][inode] = NegativeOneFlag;
        for (int iproc=0; iproc < _NoSubdom; iproc++) {
            for (int inode = _off_nd[QQ][iproc*_NoLevels];
                     inode < _off_nd[QQ][iproc*_NoLevels+ilev+1];inode++)
                _Qnode_fine_Qnode_lev[ilev][inode]=count++;//post increment
        }
    }

// linear coarse level
    _Qnode_fine_Qnode_lev[_NoLevels]=new int[_n_nodes];
    int count=0;
    for (int  inode=0;inode<_n_nodes;inode++)  _Qnode_fine_Qnode_lev[_NoLevels][inode] = NegativeOneFlag;
    for (int iproc=0; iproc < _NoSubdom; iproc++) {
        for (int  inode=0;inode <  _off_nd[LL][iproc*_NoLevels+1]
                - _off_nd[LL][iproc*_NoLevels];inode++)
            _Qnode_fine_Qnode_lev[_NoLevels][_off_nd[QQ][iproc*_NoLevels]+inode] = count++;
    }

    //After _Qnode_fine_Qnode_lev is filled, I can fill the inverse map... can I do it on the fly?!   //Let me do it here now
      _Qnode_lev_Qnode_fine = new uint *[_NoLevels+1];
     uint n_nodes_top=_NoNodesXLev[_NoLevels-1];

      for (uint ilev=0;ilev<=_NoLevels;ilev++) {
    _Qnode_lev_Qnode_fine[ilev] = new uint [_NoNodesXLev[ilev]];
     for (uint inode=0;inode<n_nodes_top;inode++) {
         uint val_lev = _Qnode_fine_Qnode_lev[ilev][inode];
      if ( val_lev != NegativeOneFlag ) _Qnode_lev_Qnode_fine[ilev][ val_lev ] = inode; 
       }
    }
    

    return;
}




} //end namespace femus


