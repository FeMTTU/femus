#include "MeshTwo.hpp"

#include <sstream>
#include <vector>
#include <algorithm>

#include "Domain.hpp"

#include "Typedefs.hpp"

#include "Files.hpp"
#include "IO.hpp"
#include "GeomEl.hpp"
#include "VBTypeEnum.hpp"
#include "FETypeEnum.hpp"

#include "paral.hpp"


namespace femus {




// ========================================================
Mesh::Mesh (const Files& files_in, const RunTimeMap<double>& map_in, const double Lref) :
         _files(files_in),
         _mesh_rtmap(map_in),
         _Lref(Lref),
         _GeomEl( (uint) map_in.get("dimension"), (uint) map_in.get("geomel_type") ),
         _n_GeomEl(map_in.get("numgeomels")),
         _dim(_mesh_rtmap.get("dimension"))
         {
	   
    _iproc    = paral::get_rank();
    _NoSubdom = paral::get_size();   
    _NoLevels = _mesh_rtmap.get("nolevels");
    
    const uint mesh_ord = (uint) _mesh_rtmap.get("mesh_ord");
    if (mesh_ord != 0) {
        std::cout << "Linear mesh not yet implemented" << std::endl;
        abort();
    }
    if (VB != 2) {
        std::cout << " NoFamFEM is not 2: to do! ";
        abort();
    }

    for (int vb=0;vb < VB; vb++) {
        _elnodes[vb][QQ]=_GeomEl._elnds[vb][mesh_ord];  //THE MESH ORD CAN ONLY BE QUADRATIC UP TO NOW!
        _elnodes[vb][LL]=_GeomEl._elnds[vb][LL];  //THE MESH ORD CAN ONLY BE QUADRATIC UP TO NOW!
        _elnodes[vb][KK]=1;
    }
    //i do not want to use the linear part actually!!

  
    _nodes_name = "/NODES";
    _elems_name = "/ELEMS";
//    _nd_coord_folder = "COORD";
//      _el_pid_name = "PID";
//     _nd_map_FineToLev = "MAP";


}

// ========================================================
Mesh::~Mesh ()  {
  
}

// ========================================================
void Mesh::clear ()  {
 
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
    delete[] _type_FEM;
    
  return;
}


// ========================================================
//This function transforms the node coordinates into the reference node coordinates
  void Mesh::TransformElemNodesToRef(const uint vb,const double* xx_qnds, double* refbox_xyz) const {
   
   double*   x_in = new double[_dim];
   double*   x_out = new double[_dim];
  const uint mesh_ord = (int) _mesh_rtmap.get("mesh_ord");
  const uint el_nds = _GeomEl._elnds[vb][mesh_ord];

      for (uint n=0;n < el_nds ;n++) {
	
   for (uint idim=0;idim < _dim; idim++)  x_in[idim] = xx_qnds[n + idim*el_nds];
  
  _domain->TransformPointToRef(x_in,x_out);

   for (uint idim=0;idim < _dim; idim++)  refbox_xyz[n + idim*el_nds] = x_out[idim];
   
      }
   
  delete[] x_in;
  delete[] x_out; 
  return; 
 }


// ========================================================
  void Mesh::SetDomain(Domain* domain_in)  {
    
    _domain = domain_in;
    
   return; 
  }
  
 // ========================================================
  Domain* Mesh::GetDomain() const {
    
   return _domain; 
   
  }


// ========================================================
/// Read mesh from hdf5 file (namefile) 
//is this function for ALL PROCESSORS 
// or only for PROC==0? Seems to be for all processors
// TODO do we need the leading "/" for opening a dataset?
void Mesh::ReadMeshFile()   {

  std::string    basemesh = DEFAULT_BASEMESH;
  std::string      ext_h5 = DEFAULT_EXT_H5;
  
  std::ostringstream meshname;
  meshname << _files._output_path << "/" << DEFAULT_CASEDIR << "/" << basemesh  << ext_h5;

//==================================
// OPEN FILE 
//==================================
  std::cout << " Reading mesh from= " <<  meshname.str() <<  std::endl;
  hid_t  file_id = H5Fopen(meshname.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
//   if (file_id < 0) {std::cout << "Mesh::read_c(): File Mesh input in data_in is missing"; abort();}
// he does not do this check, if things are wrong the H5Fopen function detects the error

//==================================
// DFLS (Dimension, VB, Levels, Subdomains)
// =====================
  uint topdata[4];
  IO::read_UIhdf5(file_id,"/DFLS",topdata);

//==================================
// CHECKS 
// ===================== 
 if (_NoLevels !=  topdata[2])  {std::cout << "Mesh::read_c. Mismatch: the number of mesh levels is " <<
   "different in the mesh file and in the configuration file" << std::endl;abort(); }


 if (_NoSubdom != topdata[3])  {std::cout << "Mesh::read_c. Mismatch: the number of mesh subdomains is " << _NoSubdom
                                   << " while the processor size of this run is " << paral::get_size()
                                   << ". Re-run gencase and your application with the same number of processors" << std::endl;abort(); }
  
//alright, there is a check to make: if you run gencase in 3D, and then you run main in 2D then things go wrong...
//this is because we read the dimension from 'dimension', we should read it from the mesh file in principle, 
//in fact it is that file that sets the space in which we are simulating...
//I'll put a check 

if (_dim != topdata[0] ) {std::cout << "Mesh::read_c. Mismatch: the mesh dimension is " << _dim
                                   << " while the dimension in the configuration file is " << _mesh_rtmap.get("dimension")
                                   << ". Recompile either gencase or your application appropriately" << std::endl;abort();}
//it seems like it doesn't print to file if I don't put the endline "<< std::endl".
//Also, "\n" seems to have no effect, "<< std::endl" must be used
//This fact doesn't seem to be related to PARALLEL processes that abort sooner than the others

if ( VB !=  topdata[1] )  {std::cout << "Mesh::read_c. Mismatch: the number of integration dimensions is " << _meshVB
                                   << " while we have VB= " << VB 
                                   << ". Re-run gencase and your application appropriately " << std::endl;abort(); }

//==================================
//======== mesh order
//==================================
//TODO should be known by the file, or hard coded there...
//the point is: what if you already have a mesh in external format?
//then you have to either read the mesh order from the external file, if such
//information can be retrieved, or provide it explicitly
    const uint mesh_ord = (int) _mesh_rtmap.get("mesh_ord");    

//==================================
// FEM element DoF number
// =====================
//Reading this is not very useful... well, it may be a check  
  _type_FEM=new uint[VB];
  IO::read_UIhdf5(file_id, "/ELNODES_VB",_type_FEM);

  for (int vb=0; vb<VB;vb++) {
if (_type_FEM[vb] !=  _GeomEl._elnds[vb][QQ] )  {std::cout << "Mesh::read_c. Mismatch: the element type of the mesh is" <<
   "different from the element type as given by the GeomEl" << std::endl; abort(); }
  }
  
// ===========================================
// ===========================================
//  NODES
// ===========================================
// ===========================================

 // ++++++++++++++++++++++++++++++++++++++++++++++++++
 // nodes X lev
 // ++++++++++++++++++++++++++++++++++++++++++++++++++
  _NoNodesXLev=new uint[_NoLevels+1];
  IO::read_UIhdf5(file_id, "/NODES/MAP/NDxLEV",_NoNodesXLev);

// ===========================================
//  COORDINATES  (COORD)
// ===========================================
  //in the mesh file now I have the coordinates for all levels, but let me read only those of the FINEST level
   double ILref = 1./_Lref;
   int lev_for_coords = _NoLevels-1;
  uint n_nodes =_NoNodesXLev[lev_for_coords];
  _xyz=new double[_dim*n_nodes];
  double *coord;coord=new double[n_nodes];
  for (uint kc=0;kc<_dim;kc++) {
    std::ostringstream Name; Name << "NODES/COORD/X" << kc+1 << "_L" << lev_for_coords;
    IO::read_Dhdf5(file_id,Name.str().c_str(),coord);
    for (uint inode=0;inode<n_nodes;inode++) _xyz[inode+kc*n_nodes]=coord[inode]*ILref; //NONdimensionalization!!!
  }
  delete []coord;

// ===================================================
//  OFF_ND
// ===================================================
_off_nd=new int*[QL_NODES];
	    
for (int fe=0;fe < QL_NODES; fe++)    {
  _off_nd[fe]=new int[_NoSubdom*_NoLevels+1];
  std::ostringstream namefe; namefe << "/NODES/MAP/OFF_ND" << "_F" << fe;
   IO::read_Ihdf5(file_id,namefe.str().c_str(),_off_nd[fe]);
  }
  
  
   // ====================================================
 // NODE MAP
 // ====================================================
  // node mapping
  //this map "mixes" linear and quadratic, in the sense that 
  // if you have levels 0 1 2, you can define 
  // 3 levels for the QUADRATIC DOFS, which are A=0 B=1 C=2,
  //and similarly you can define 3 levels for the linear dofs.
  //In order to be consistent with the levels you can define 
  // another auxiliary level 3 which is 
  // the coarse linear level
  //in each _node_map there is no "subdomain" stuff, only LEVEL
  // Each node_map contains the list of NODES of THAT LEVEL,
  //of course in the FINE NODE NUMBERING
  //in practice gives us, for every level, the FINE NODE INDICES corresponding to the DOFS
  //so it is like an INVERSE DOF MAP: from DOF TO NODE.
  //if you use levels 0 1 2, you do from QUADRATIC dof to FINE NODE
  //if you use levels 3 0 1, you do from LINEAR dof to FINE NODE
  
  uint n_nodes_top=_NoNodesXLev[_NoLevels-1];
  _Qnode_lev_Qnode_fine = new uint *[_NoLevels+1];
  _Qnode_fine_Qnode_lev = new int *[_NoLevels+1];

  for (uint ilev=0;ilev<=_NoLevels;ilev++) { //loop on Extended levels 
    _Qnode_lev_Qnode_fine[ilev] = new uint [_NoNodesXLev[ilev]];
    _Qnode_fine_Qnode_lev[ilev] = new int [n_nodes_top];  //THIS HAS TO BE INT because it has -1!!!
    std::ostringstream Name; Name << "/NODES/MAP/MAP" << "_XL" << ilev;
    IO::read_Ihdf5(file_id,Name.str().c_str(),_Qnode_fine_Qnode_lev[ilev]);
    for (uint inode=0;inode<n_nodes_top;inode++) {
         int val_lev = _Qnode_fine_Qnode_lev[ilev][inode];
      if ( val_lev != -1 ) _Qnode_lev_Qnode_fine[ilev][ val_lev ] = inode; //this doesnt have -1 numbers
      }
      //This is how you read the _node_map: 
      //- you remove the "-1" that come from the gencase
      //- and you dont put the content (we dont need it!) but the position,
      //   i.e. the node index in the fine numbering
   }
  
// ===========================================
//   /ELEMS
// ===========================================
   
// ===========================================
//   NUMBER EL
// ===========================================
  _n_elements_vb_lev=new uint*[VB];
  for (uint vb=0;vb< VB;vb++) {
    _n_elements_vb_lev[vb]=new uint[_NoLevels];
    std::ostringstream Name; Name << "/ELEMS/VB" << vb  <<"/NExLEV";
    IO::read_UIhdf5(file_id,Name.str().c_str(),_n_elements_vb_lev[vb]);
  }

// ===========================================
//   OFF_EL
// ===========================================
  _off_el=new int*[VB];

for (int vb=0; vb < VB; vb++)    {
  _off_el[vb] = new int [_NoSubdom*_NoLevels+1];
  std::ostringstream offname; offname << "/ELEMS/VB" << vb << "/OFF_EL";
    IO::read_Ihdf5(file_id,offname.str().c_str(),_off_el[vb]);
}

// ===========================================
//   CONNECTIVITY
// ===========================================
  _el_map=new uint*[VB];
for (int vb=0; vb < VB; vb++)    {
  _el_map[vb]=new uint [_off_el[vb][_NoSubdom*_NoLevels]*_GeomEl._elnds[vb][mesh_ord]];
  std::ostringstream elName; elName << "/ELEMS/VB" << vb  <<"/CONN";
  IO::read_UIhdf5(file_id,elName.str().c_str(),_el_map[vb]);
}

// ===========================================
//   BDRY EL TO VOL EL
// ===========================================
   _el_bdry_to_vol = new int*[_NoLevels];
  for (uint lev=0; lev < _NoLevels; lev++)    {
  _el_bdry_to_vol[lev] = new int[_n_elements_vb_lev[BB][lev]];
    std::ostringstream btov; btov << "/ELEMS/BDRY_TO_VOL_L" << lev;
  IO::read_Ihdf5(file_id, btov.str().c_str(),_el_bdry_to_vol[lev]);
  }
  
 // ===========================================
 //  CLOSE FILE 
 // ===========================================
  H5Fclose(file_id);

  return; 
   
}


// ===============================================================
//   print_mesh_h5(_nd_fm_libm, _nd_libm_fm, _el_fm_libm, _el_fm_libm_b);
void Mesh::PrintMeshHDF5() const  {

    std::ostringstream name;

    std::string basepath  = _files._app_path;
    std::string input_dir = DEFAULT_CASEDIR;
    std::string basemesh  = DEFAULT_BASEMESH;
    std::string ext_h5    = DEFAULT_EXT_H5;

    std::ostringstream inmesh;
    inmesh << basepath << "/" << input_dir << basemesh << ext_h5;

//==================================
// OPEN FILE
//==================================
    hid_t file = H5Fcreate(inmesh.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);

//==================================
// DFLS (Dimension, VB, Levels, Subdomains)
// =====================
    int *tdata;

    tdata    = new int[4];
    tdata[0] = _dim;
    tdata[1] = VB;
    tdata[2] = _NoLevels;
    tdata[3] = _NoSubdom;

    hsize_t dimsf[2];
    dimsf[0] = 4;
    dimsf[1] = 1;
    IO::print_Ihdf5(file,"DFLS", dimsf,tdata);
    delete [] tdata;

//==================================
// FEM element DoF number
// =====================
    int *ttype_FEM;
    ttype_FEM=new int[VB];

    for (uint vb=0; vb< VB;vb++)   ttype_FEM[vb]=_GeomEl._elnds[vb][QQ];

    dimsf[0] = VB;
    dimsf[1] = 1;
    IO::print_Ihdf5(file,"ELNODES_VB", dimsf,ttype_FEM);

// ===========================================
// ===========================================
//  NODES
// ===========================================
// ===========================================
    hid_t group_id = H5Gcreate(file, _nodes_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

// ++++ NODES/MAP ++++++++++++++++++++++++++++++++++++++++++++++
    std::string ndmap = _nodes_name + "/MAP";
    hid_t subgroup_id = H5Gcreate(file, ndmap.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
// nodes X lev
    std::string ndxlev =  ndmap + "/NDxLEV";
    dimsf[0] = _NoLevels+1;
    dimsf[1] = 1;
    IO::print_UIhdf5(file, ndxlev.c_str(), dimsf,_NoNodesXLev);

// ++++++++++++++++++++++++++++++++++++++++++++++++++
// node map (XL, extended levels)
// ++++++++++++++++++++++++++++++++++++++++++++++++++
    dimsf[0] = _n_nodes;
    dimsf[1] = 1;

    for (int ilev= 0;ilev < _NoLevels+1; ilev++) {
        name.str("");
        name << ndmap << "/MAP"<< "_XL" << ilev;
        IO::print_Ihdf5(file,name.str(),dimsf,_Qnode_fine_Qnode_lev[ilev]);
    }

// ++++++++++++++++++++++++++++++++++++++++++++++++++
//   OFF_ND: node offset quadratic and linear
// ++++++++++++++++++++++++++++++++++++++++++++++++++
    dimsf[0] = _NoSubdom*_NoLevels+1;
    dimsf[1] = 1;
    for (int fe=0;fe < QL_NODES; fe++) {
        std::ostringstream namefe;
        namefe <<  ndmap << "/OFF_ND" << "_F" << fe;
        IO::print_Ihdf5(file, namefe.str(),dimsf,_off_nd[fe]);
    }

    H5Gclose(subgroup_id);

// ===========================================
//  COORDINATES  (COORD)
// ===========================================
    //ok, we need to print the coordinates of the nodes for each LEVEL
    //The array _nod_coords holds the coordinates of the FINE Qnodes
    //how do we take the Qnodes of each level?
    //Well, I'd say we need to use  _off_nd for the quadratics
    //Ok, the map _nd_fm_libm goes from FINE FEMUS NODE ORDERING to FINE LIBMESH NODE ORDERING
    //I need to go from the QNODE NUMBER at LEVEL 
    //to the QNODE NUMBER at FINE LEVEL according to FEMUS,
    //and finally to the number at fine level according to LIBMESH.
    //The fine level is the only one where the Qnode numbering is CONTIGUOUS
    
    
    
    std::string ndcoords = _nodes_name + "/COORD";
    subgroup_id = H5Gcreate(file, ndcoords.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // nodes are PRINTED ACCORDING to FEMUS ordering, which is inode, i.e. the INVERSE of v[inode].second
    //you use this because you do DIRECTLY a NODE LOOP
    //now let us loop coordinates over ALL LEVELS
    for (int l=0; l<_NoLevels; l++)  {
    double * xcoord = new double[_NoNodesXLev[l]];
    for (int kc=0;kc<3;kc++) {
      
      int Qnode_lev=0; 
      for (uint isubdom=0; isubdom<_NoSubdom; isubdom++) {
            uint off_proc=isubdom*_NoLevels;
               for (int k1 = _off_nd[QQ][off_proc];
                        k1 < _off_nd[QQ][off_proc + l+1 ]; k1++) {
		 int Qnode_fine_fm = _Qnode_lev_Qnode_fine[l][Qnode_lev];
		 xcoord[Qnode_lev] = _nd_coords_libm[ _nd_fm_libm[Qnode_fine_fm].second + kc*_n_nodes ];
		  Qnode_lev++; 
		 }
	      } //end subdomain
	    
// old        for (int inode=0; inode <_n_nodes_lev[l];inode++)  xcoord[inode] = _nod_coords[_nd_fm_libm[inode].second+kc*_n_nodes]; //the offset is fine
        dimsf[0] = _NoNodesXLev[l];
        dimsf[1] = 1;
        name.str("");
        name << ndcoords << "/X" << kc+1<< "_L" << l;
        IO::print_Dhdf5(file,name.str(), dimsf,xcoord);
    }

    delete [] xcoord;
     }  //levels
    
    H5Gclose(subgroup_id);

    H5Gclose(group_id);

// ===========================================
// ===========================================
//   /ELEMS
// ===========================================
// ===========================================

    group_id = H5Gcreate(file, _elems_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    ElemStoBase** elsto_out;   //TODO delete it!
    elsto_out = new ElemStoBase*[_n_elements_sum_levs[VV]];
    for (int i=0;i<_n_elements_sum_levs[VV];i++) {
        elsto_out[i]= static_cast<ElemStoBase*>(_el_sto[i]);
    }

    ElemStoBase** elstob_out;
    elstob_out = new ElemStoBase*[_n_elements_sum_levs[BB]];
    for (int i=0; i<_n_elements_sum_levs[BB]; i++) {
        elstob_out[i]= static_cast<ElemStoBase*>(_el_sto_b[i]);
    }

    PrintElemVB(file,VV,_nd_libm_fm, elsto_out,_el_fm_libm);
    PrintElemVB(file,BB,_nd_libm_fm, elstob_out,_el_fm_libm_b);

    // ===============
    // print child to father map for all levels for BOUNDARY ELEMENTS
    // ===============

    for (int lev=0;lev<_NoLevels; lev++)  {
        std::ostringstream   bname;
        bname << _elems_name << "/BDRY_TO_VOL_L" << lev;
        dimsf[0] = _n_elements_vb_lev[BB][lev];
        dimsf[1] = 1;
        IO::print_Ihdf5(file,bname.str(), dimsf,_el_child_to_fath[lev]);
    }




//             std::cout <<  "==================" << std::endl;
//          for (int i=0;i<_n_elements_vb_lev[BB][lev];i++)  {
//              std::cout <<  _el_child_to_fath[lev][i] << std::endl;
//        }
// 	 }



// ok, so, i had to do two cast vectors so i could pass them both to the print_elem_vb routine
//now i have to be careful in destroying these, in relation with their fathers...

    //delete temp  //I AM HERE TODO
//     for (int i=0;i< _n_elements_sum_levs_vb[BB];i++) {   delete /*[]*/ elstob_out[i]; } // delete [] el_sto[i]; with [] it doesnt work
//    delete [] elstob_out;
//     for (int i=0;i< _n_elements_sum_levs_vb[VV];i++) {   delete /*[]*/ elsto_out[i]; } // delete [] el_sto[i]; with [] it doesnt work
//    delete [] elsto_out;
    H5Gclose(group_id);

// ===========================================
//   PID
// ===========================================
    group_id = H5Gcreate(file, "/PID", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (int  vb= 0; vb< VB;vb++) {
        for (int  ilev= 0;ilev< _NoLevels; ilev++)   PrintSubdomFlagOnQuadrCells(vb,ilev,inmesh.str().c_str());
    }

    H5Gclose(group_id);

// ===========================================
//  CLOSE FILE
// ===========================================
    H5Fclose(file);

    //============
    delete [] ttype_FEM;

    return;
}


// ========================================================
/// Write mesh to hdf5 file (namefile) 
///              as Mesh class (Mesh.h): 

void Mesh::PrintMeshFile (const std::string & /*namefile*/) const
{ 

  std::cout << "Implement it following the Gencase print mesh function" << std::endl; abort();
  
  return; 
}



// ==========================================================
//prints conn and stuff for either vol or bdry mesh
//When you have to construct the connectivity,
//you go back to the libmesh elem ordering,
//then you pick the nodes of that element in LIBMESH numbering,
//then you pick the nodes in femus NUMBERING,
//and that's it
void Mesh::PrintElemVB(hid_t file,const uint vb , int* nd_libm_fm , ElemStoBase** el_sto_in, std::vector<std::pair<int,int> > el_fm_libm_in ) const {

    std::ostringstream name;

    std::string auxvb[VB];
    auxvb[0]="0";
    auxvb[1]="1";
    std::string elems_fem = _elems_name;
    std::string elems_fem_vb = elems_fem + "/VB" + auxvb[vb];  //VV later

    hid_t subgroup_id = H5Gcreate(file, elems_fem_vb.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dimsf[2];
    dimsf[0] = 2;
    dimsf[1] = 1;
    int ndofm[2];
    ndofm[0]=_elnodes[vb][QQ];
    ndofm[1]=_elnodes[vb][LL];
    IO::print_Ihdf5(file,(elems_fem_vb + "/NDOF_FO_F1"), dimsf,ndofm);
    // NoElements ------------------------------------
    dimsf[0] = _NoLevels;
    dimsf[1] = 1;
    IO::print_UIhdf5(file,(elems_fem_vb + "/NExLEV"), dimsf,_n_elements_vb_lev[vb]);
    // offset
    dimsf[0] = _NoSubdom*_NoLevels+1;
    dimsf[1] = 1;
    IO::print_Ihdf5(file,(elems_fem_vb + "/OFF_EL"), dimsf,_off_el[vb]);

    // Mesh 1 Volume at all levels  packaging data  (volume) ----------
    //here you pick all the elements at all levels,
    //and you print their connectivities according to the libmesh ordering
    int *tempconn;
    tempconn=new int[_n_elements_sum_levs[vb]*_elnodes[vb][QQ]]; //connectivity of all levels
    for (int ielem=0;ielem<_n_elements_sum_levs[vb];ielem++) {
        for (uint inode=0;inode<_elnodes[vb][QQ];inode++) {
            int el_libm =   el_fm_libm_in[ielem].second;
            int nd_libm = el_sto_in[el_libm]->_elnds[inode];
            tempconn[inode+ielem*_elnodes[vb][QQ]] = nd_libm_fm[nd_libm];
        }
    }
    // global mesh hdf5 storage ------------------------------
    dimsf[0] = _n_elements_sum_levs[vb]*_elnodes[vb][QQ];
    dimsf[1] = 1;
    IO::print_Ihdf5(file,(elems_fem_vb + "/CONN"), dimsf,tempconn);

    // level connectivity ---------------------------------
    for (int ilev=0;ilev <_NoLevels; ilev++) {

        int *conn_lev=new int[_n_elements_vb_lev[vb][ilev]*_elnodes[vb][QQ]];  //connectivity of ilev

        int ltot=0;
        for (int iproc=0;iproc <_NoSubdom; iproc++) {
            for (int iel = _off_el[vb][iproc*_NoLevels+ilev];
                     iel < _off_el[vb][iproc*_NoLevels+ilev+1]; iel++) {
                for (uint inode=0;inode<_elnodes[vb][QQ];inode++) {
                    conn_lev[ltot*_elnodes[vb][QQ]+inode]=
                        tempconn[  iel*_elnodes[vb][QQ]+inode];
                }
                ltot++;
            }
        }
        // hdf5 storage     ----------------------------------
        dimsf[0] = _n_elements_vb_lev[vb][ilev]*_elnodes[vb][QQ];
        dimsf[1] = 1;

        name.str("");
        name << elems_fem_vb << "/CONN" << "_L" << ilev ;
        IO::print_Ihdf5(file,name.str(), dimsf,conn_lev);
        //clean
        delete []conn_lev;

    }


    delete []tempconn;

    H5Gclose(subgroup_id);

    return;


}

//=============================================
//this is different from the Mesh class function
//because we are using quadratic elements
void Mesh::PrintSubdomFlagOnQuadrCells(const int vb, const int Level, std::string filename) const {

    if (_iproc==0)   {

        //   const uint Level = /*_NoLevels*/_n_levels-1;
        const uint n_children = /*4*(_dim-1)*/1;  /*here we have quadratic cells*/

        uint      n_elements = _n_elements_vb_lev[vb][Level];
        int *ucoord;
        ucoord=new int[n_elements*n_children];
        int cel=0;
        for (int iproc=0; iproc<_NoSubdom; iproc++) {
            for (int iel = _off_el[vb][Level  + iproc*_NoLevels];
                     iel < _off_el[vb][Level+1 + iproc*_NoLevels]; iel++) {
                for (uint is=0; is< n_children; is++)
                    ucoord[cel*n_children + is]=iproc;
                cel++;
            }
        }

        hid_t file_id = H5Fopen(filename.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
        hsize_t dimsf[2];
        dimsf[0] = n_elements*n_children;
        dimsf[1] = 1;
        std::ostringstream name;
        name << "/PID/PID_VB" << vb << "_L" << Level;

        IO::print_Ihdf5(file_id,name.str(),dimsf,ucoord);

        H5Fclose(file_id);

    }

    return;
}



// ========================================================
/// It manages the printing in Xdmf format
void Mesh::PrintForVisualizationAllLEVAllVB()  const {
  
    const uint iproc=_iproc;
   if (iproc==0) {
       PrintConnLinAllLEVAllVB();
       PrintXDMFAllLEVAllVB();
    }
   
   return;
}


// ==================================================================
void Mesh::PrintMultimeshXdmf() const {

    std::string basepath  = _files._app_path;
    std::string input_dir = DEFAULT_CASEDIR;
    std::string multimesh = DEFAULT_MULTIMESH;
    std::string ext_xdmf  = DEFAULT_EXT_XDMF;
    std::string basemesh  = DEFAULT_BASEMESH;
    std::string ext_h5    = DEFAULT_EXT_H5;

    std::ostringstream inmesh_xmf;
    inmesh_xmf << basepath << "/" << input_dir << multimesh << ext_xdmf;
    std::ofstream out(inmesh_xmf.str().c_str());

    std::ostringstream top_file;
    top_file << basemesh << ext_h5;

    if (out.fail()) {
        std::cout << "Mesh::PrintMultimeshXdmf: The file is not open" << std::endl;
        abort();
    }

//it seems that there is no control on this, if the directory isn't there
//it doesnt give problems

    //strings for VB
    std::string  meshname[VB];
    meshname[VV]="VolumeMesh";
    meshname[BB]="BoundaryMesh";

    out << "<?xml version=\"1.0\" ?> \n";
    out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" \n";
//   out << "[ <!ENTITY HeavyData \"mesh.h5\"> ] ";
    out << "> \n";
    out << " \n";
    out << "<Xdmf> \n";
    out << "<Domain> \n";

    for (int vb=0;vb< VB;vb++) {

        for (int ilev=0; ilev<_NoLevels; ilev++) {
	  
            out << "<Grid Name=\"" << meshname[vb] << ilev <<"\"> \n";

            std::ostringstream hdf5_field;
            hdf5_field << _elems_name << "/VB" << vb << "/CONN" << "_L" << ilev;
            IO::PrintXDMFTopology(out,top_file.str(),hdf5_field.str(),_GeomEl.name[vb],_n_elements_vb_lev[vb][ilev],_n_elements_vb_lev[vb][ilev],_elnodes[vb][QQ]);
            std::ostringstream coord_lev; coord_lev << "_L" << ilev; 
	    IO::PrintXDMFGeometry(out,top_file.str(),_nodes_name+"/COORD/X",coord_lev.str(),"X_Y_Z","Float",_NoNodesXLev[ilev],1);
            std::ostringstream pid_field;
            pid_field << "PID/PID_VB"<< vb <<"_L"<< ilev;
            IO::PrintXDMFAttribute(out,top_file.str(),pid_field.str(),"PID","Scalar","Cell","Int",_n_elements_vb_lev[vb][ilev],1);

            out << "</Grid> \n";
	    
        }
    }

    out << "</Domain> \n";
    out << "</Xdmf> \n";
    out.close();

    return;
}

// ========================================================
/// It prints the volume/boundary Mesh (connectivity) in Xdmf format
 //this is where the file mesh.xmf goes: it is a good rule that the
 //file .xmf and the .h5 it needs should be in the same folder (so that READER and DATA are ALWAYS TOGETHER)
 //well, there can be MANY READERS for ONE DATA file (e.g. sol.xmf, case.xmf)
 // or MANY DATA files for  ONE READER,
 //but if you put ALL THE DATA in THE SAME FOLDER as the READER(s)
 //then they are always independent and the data will always be readable

 void Mesh::PrintXDMFAllLEVAllVB() const {

  std::string     basemesh = DEFAULT_BASEMESH;
  std::string     ext_xdmf = DEFAULT_EXT_XDMF;
  std::string       ext_h5 = DEFAULT_EXT_H5;
  std::string      connlin = DEFAULT_CONNLIN;

  std::ostringstream top_file; top_file << basemesh << connlin << ext_h5;
  std::ostringstream geom_file; geom_file << basemesh << ext_h5;

  std::ostringstream namefile;
  namefile << _files._output_path << "/" << basemesh << ext_xdmf;
 
  std::ofstream out (namefile.str().c_str());

  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" \n"; 
  out << " [ <!ENTITY HeavyData \"mesh.h5 \"> ] ";
  out << ">\n"; 
  out << " \n";
  out << "<Xdmf> \n" << "<Domain> \n";
  
          for(uint vb=0;vb< VB; vb++)  
	    for(uint l=0; l< _NoLevels; l++)
	      PrintXDMFGridVB(out,top_file,geom_file,l,vb);
   
   out << "</Domain> \n" << "</Xdmf> \n";
   out.close ();
   
   return;
}



// ========================================================
void Mesh::PrintXDMFGridVB(std::ofstream& out,
			      std::ostringstream& top_file,
			      std::ostringstream& geom_file, const uint Level, const uint vb) const  {

  std::string grid_mesh[VB];
  grid_mesh[VV]="Volume";
  grid_mesh[BB]="Boundary";
  
  std::ostringstream hdf_field; hdf_field << "MSHCONN_VB_" << vb << "_LEV_" << Level;
  
    uint nel = _n_elements_vb_lev[vb][Level];

    out << "<Grid Name=\"" << grid_mesh[vb].c_str() << "_L" << Level << "\"> \n";
    
   IO::PrintXDMFTopology(out,top_file.str(),hdf_field.str(),_GeomEl.pname[vb],nel*_GeomEl.n_se[vb],nel*_GeomEl.n_se[vb],_GeomEl._elnds[vb][LL]);    

   std::ostringstream coord_lev; coord_lev << "_L" << Level; 
   IO::PrintXDMFGeometry(out,geom_file.str(),"NODES/COORD/X",coord_lev.str(),"X_Y_Z","Float",_NoNodesXLev[Level],1);
    
    out << "</Grid> \n";  

   return;
}



// ========================================================
/// It prints the connectivity in hdf5 format
/// The changes are only for visualization of quadratic FEM

void Mesh::PrintConnLinAllLEVAllVB() const { 

    std::string    basemesh = DEFAULT_BASEMESH;
    std::string    ext_h5   = DEFAULT_EXT_H5;
    std::string    connlin  = DEFAULT_CONNLIN;

  std::ostringstream namefile;
  namefile << _files._output_path << "/" << basemesh << connlin << ext_h5; 

  std::cout << namefile.str() << std::endl;
  hid_t file = H5Fcreate (namefile.str().c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);  //TODO VALGRIND

//================================
// here we loop both over LEVELS and over VB, so everything is inside a UNIQUE FILE  
  for(uint l=0; l< _NoLevels; l++)
           for(uint vb=0;vb< VB; vb++)
	        PrintConnLinVB(file,l,vb);
//================================
	   
  H5Fclose(file); //TODO VALGRIND Invalid read of size 8: this is related to both H5Fcreate and H5Dcreate
  return;
}



void Mesh::PrintConnLinVB(hid_t file, const uint Level, const uint vb) const {
  
   int conn[8][8];   uint *gl_conn;
  
    uint icount=0;
    uint mode = _GeomEl._elnds[vb][QQ];
    uint n_elements = _n_elements_vb_lev[vb][Level];
    uint nsubel, nnodes;

    switch(_dim)   {
      case 2:  {
    switch(mode){
      // -----------------------------------
      case 9:	// Quad 9  0-4-1-5-2-6-3-7-8
      gl_conn=new uint[n_elements*4*4];
      conn[0][0] = 0;  conn[0][1] = 4;   conn[0][2] = 8; conn[0][3] = 7;// quad4  0-4-8-7
      conn[1][0] = 4;  conn[1][1] = 1;   conn[1][2] = 5; conn[1][3] = 8;// quad4  4-1-5-8
      conn[2][0] = 8;  conn[2][1] = 5;   conn[2][2] = 2; conn[2][3] = 6;// quad4  8-5-2-6
      conn[3][0] = 7;  conn[3][1] = 8;   conn[3][2] = 6; conn[3][3] = 3;// quad4  7-8-6-3
      nsubel=4;nnodes=4;
      break;
//=====================
      case 6:	// Quad 9  0-4-1-5-2-6-3-7-8
      gl_conn=new uint[n_elements*4*3];
      conn[0][0] = 0;  conn[0][1] = 3;   conn[0][2] = 5; // quad4  0-4-8-7
      conn[1][0] = 3;  conn[1][1] = 4;   conn[1][2] = 5; // quad4  4-1-5-8
      conn[2][0] = 3;  conn[2][1] = 1;   conn[2][2] = 4; // quad4  8-5-2-6
      conn[3][0] = 4;  conn[3][1] = 2;   conn[3][2] = 5; // quad4  7-8-6-3
      nsubel=4;nnodes=3;
      break;
//===================
     case 3: // boundary edge 2 linear 0-2-1
      gl_conn=new uint[n_elements*2*2];
      conn[0][0] = 0; conn[0][1] = 2;		// element 0-2
      conn[1][0] = 2; conn[1][1] = 1;		// element 1-2
      nsubel=2;nnodes=2;
     break;
      // -----------------------------------------
      default:   // interior 3D
      gl_conn = new uint[n_elements*mode];
      for (uint n=0;n<mode;n++) conn[0][n] = n;
      nsubel=1;nnodes=mode;
      break;
      } //end switch mode 1
 
 break; //end case 2
      }
      case 3:  {
         switch(mode){
  // ----------------------
      case  27: //  Hex 27 (8 Hex8)
      gl_conn=new uint[n_elements*8*8];
      conn[0][0] = 0;  conn[0][1] = 8;  conn[0][2] = 20; conn[0][3] = 11;
      conn[0][4] = 12; conn[0][5] = 21; conn[0][6] = 26; conn[0][7] = 24;
      conn[1][0] = 8;  conn[1][1] = 1;  conn[1][2] = 9;  conn[1][3] = 20;
      conn[1][4] = 21; conn[1][5] = 13; conn[1][6] = 22; conn[1][7] = 26;
      conn[2][0] = 11; conn[2][1] = 20; conn[2][2] = 10; conn[2][3] = 3;
      conn[2][4] = 24; conn[2][5] = 26; conn[2][6] = 23; conn[2][7] = 15;
      conn[3][0] = 20; conn[3][1] = 9;  conn[3][2] = 2;  conn[3][3] = 10;
      conn[3][4] = 26; conn[3][5] = 22; conn[3][6] = 14; conn[3][7] = 23;
      conn[4][0] = 12; conn[4][1] = 21;  conn[4][2] = 26; conn[4][3] = 24;
      conn[4][4] = 4;  conn[4][5] = 16;  conn[4][6] = 25; conn[4][7] = 19;
      conn[5][0] = 21;  conn[5][1] = 13; conn[5][2] = 22; conn[5][3] = 26;
      conn[5][4] = 16; conn[5][5] = 5;  conn[5][6] = 17; conn[5][7] = 25;
      conn[6][0] = 24; conn[6][1] = 26; conn[6][2] = 23; conn[6][3] = 15;
      conn[6][4] = 19; conn[6][5] = 25; conn[6][6] = 18; conn[6][7] = 7;
      conn[7][0] = 26; conn[7][1] = 22; conn[7][2] = 14; conn[7][3] = 23;
      conn[7][4] = 25; conn[7][5] = 17; conn[7][6] = 6;  conn[7][7] = 18;
      nsubel=8;nnodes=8;     
      break;
    // ---------------------------------------  
      case  10: // Tet10
      gl_conn=new uint[n_elements*8*4];
      conn[0][0] = 0; conn[0][1] = 4;  conn[0][2] = 6; conn[0][3] = 7;
      conn[1][0] = 4; conn[1][1] = 1;  conn[1][2] = 5; conn[1][3] = 8;
      conn[2][0] = 5; conn[2][1] = 2;  conn[2][2] = 6; conn[2][3] = 9;
      conn[3][0] = 7; conn[3][1] = 8;  conn[3][2] = 9; conn[3][3] = 3;
      conn[4][0] = 4; conn[4][1] = 8;  conn[4][2] = 6; conn[4][3] = 7;
      conn[5][0] = 4; conn[5][1] = 5;  conn[5][2] = 6; conn[5][3] = 8;
      conn[6][0] = 5; conn[6][1] = 9;  conn[6][2] = 6; conn[6][3] = 8;
      conn[7][0] = 7; conn[7][1] = 6;  conn[7][2] = 9; conn[7][3] = 8;
      nsubel=8;nnodes=4;
      break;
  // ---------------------------------------  
      case 6:	// Tri6  0-3-1-4-2-5
      gl_conn=new uint[n_elements*4*3];
      conn[0][0] = 0;  conn[0][1] = 3;   conn[0][2] = 5; // quad4  0-4-8-7
      conn[1][0] = 3;  conn[1][1] = 4;   conn[1][2] = 5; // quad4  4-1-5-8
      conn[2][0] = 3;  conn[2][1] = 1;   conn[2][2] = 4; // quad4  8-5-2-6
      conn[3][0] = 4;  conn[3][1] = 2;   conn[3][2] = 5; // quad4  7-8-6-3
      nsubel=4;nnodes=3;
      break;    
    // ---------------------------------------  
      case 9:  // Quad9 elements ( 4 Quad4)
      gl_conn=new uint[n_elements*4*4];
      conn[0][0] = 0; conn[0][1] =4; conn[0][2] = 8; conn[0][3] = 7;
      conn[1][0] = 4; conn[1][1] =1; conn[1][2] = 5; conn[1][3] = 8;
      conn[2][0]= 8;  conn[2][1] =5; conn[2][2] = 2; conn[2][3] = 6;
      conn[3][0] = 7; conn[3][1] = 8;conn[3][2] = 6; conn[3][3] = 3;
      nsubel=4;nnodes=4;
      break;
     // -----------------------------------------
    default:   // interior 3D
      gl_conn = new uint[n_elements*mode];
      for (uint n=0;n<mode;n++) conn[0][n] = n;
      nsubel=1;nnodes=mode;  
       break;
      } //end switch mode 2

      break;
      } //end case 3

       default:
      std::cout << "PrintConnLin not working" << std::endl; abort();
   break;
    }
    
      // mapping
    for (uint iproc = 0; iproc < _NoSubdom; iproc++) {
      for (int el = _off_el[vb][iproc*_NoLevels+Level];
           el <_off_el[vb][iproc*_NoLevels+Level+1]; el++) {
        for (uint se = 0; se < nsubel; se++) {
          for (uint i = 0; i < nnodes; i++) {
	    uint Qnode_fine = _el_map[vb][el*mode+conn[se][i]];
	    uint Qnode_lev = _Qnode_fine_Qnode_lev[Level][Qnode_fine];
            gl_conn[icount] = Qnode_lev;
	    icount++;
          }
        }
      }
    } 
    
     // Print mesh in hdf files
    hsize_t dimsf[2];  dimsf[0] =icount;  dimsf[1] = 1;
    hid_t dtsp = H5Screate_simple(2, dimsf, NULL);
    std::ostringstream Name; Name << "MSHCONN_VB_" << vb << "_LEV_" << Level;
    hid_t dtset = H5Dcreate(file,Name.str().c_str(),H5T_NATIVE_INT,dtsp,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);  //TODO VALGRIND
    H5Dwrite(dtset,H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT, gl_conn);
    H5Sclose(dtsp);
    H5Dclose(dtset);
    
    delete[] gl_conn;
  
    return;
}
  
  
  



// ===============================================================
/// this function is done only by _iproc == 0!
/// it prints the PID index on the cells of the linear mesh
void Mesh::PrintSubdomFlagOnLinCells(std::string filename) const {
  
 if (_iproc==0)   {

  hid_t file_id = H5Fopen(filename.c_str(),H5F_ACC_RDWR, H5P_DEFAULT); 

  for (uint l=0; l<_NoLevels; l++) {
    
  const uint n_children = 4*(_dim-1);
  uint n_elements = _n_elements_vb_lev[VV][l];
  int *ucoord;   ucoord=new int[n_elements*n_children];
  int cel=0;
  for (uint iproc=0; iproc < _NoSubdom; iproc++) {
    for (int iel =  _off_el[VV][ l   + iproc*_NoLevels];
              iel < _off_el[VV][ l+1 + iproc*_NoLevels]; iel++) {
      for (uint is=0; is< n_children; is++)      
	ucoord[cel*n_children + is]=iproc;
      cel++;
    }
  }
  
  hsize_t dimsf[2]; dimsf[0] = n_elements*n_children; dimsf[1] = 1;
  std::ostringstream pidname; pidname << "PID" << "_LEVEL" << l;
  
  IO::print_Ihdf5(file_id,pidname.str(),dimsf,ucoord);
  
     } //end levels
  
    H5Fclose(file_id);
    
  }
  
return;
}


// ==================================================================
//ok, this doesnt print the femus ordering like i wish, when i am in multiprocessor
//now, let us see.
//The point is that there is some mismatch in the construction of this function.
// In fact, this function receives a BOUNDARY ELEM in the LEVEL NUMERATION
// and it yields a VOLUME ELEMENT in the GLOBAL NUMERATION.
// Since I would prefer doing things by SEPARATE LEVEL, let us see.

// Now my point is: i want to convert from the ABSOLUTE VOLUME NUMBER to the "RELATIVE per LEVEL" VOLUME NUMBER.
//So for every level, I still want to pick the elements of only THAT level but STARTING from ZERO
// 
// 
     //ok now I modified this function consistently, in such a way that for every LEVEL i have an el_bdry_to_vol map
     //that deals with BOUNDARY ELEM numbers and VOLUME ELEM NUMBERS *at THAT LEVEL*
     // so now I already have the  element number IN THAT LEVEL, so that is what I want.

void Mesh::ElemChildToFather() {

    _el_child_to_fath = new int*[_NoLevels];
    for (int lev=0; lev<_NoLevels; lev++)  {
        std::cout << "Level=== " << lev << std::endl;
        _el_child_to_fath[lev] =  new int[_n_elements_vb_lev[BB][lev]];
        int ielem = 0;
        for (int isub=0;isub<_NoSubdom;isub++) {
            std::cout << "Subdomain= " << isub << std::endl;
            for (int iel = _off_el[BB][isub*_NoLevels+lev];
                     iel < _off_el[BB][isub*_NoLevels+lev+1]; iel++) {
                int   el_libm_b = _el_fm_libm_b[iel].second;
//                 std::cout << "boundary elem libmesh" << el_libm_b << std::endl;
                int father_id_libm = _el_sto_b[el_libm_b]->_vol_id;
//                 std::cout << "Volume elem libmesh " << father_id_libm << std::endl;
                int father_id_fm = _el_libm_fm[father_id_libm];
//                 std::cout << "Volume elem fm " << father_id_fm << std::endl;
                int boundary_proc = _el_sto_b[el_libm_b]->_subd;
//                 std::cout << "BOUNDARY PROC " << boundary_proc << std::endl;
                int volume_proc   = _el_sto[father_id_libm]->_subd;
		// father_id_fm is the ABSOLUTE number
		// I want to have the number SO THAT EACH LEVEL RESTARTS FROM ZERO
		// therefore, i want to remove 
		int el_sum_prev_sd_at_lev = 0;
		    for (int sd=0; sd < isub; sd++)  el_sum_prev_sd_at_lev += _off_el[VV][sd*_NoLevels + lev+1] - _off_el[VV][sd*_NoLevels + lev];
		int vol_elem_per_lev = father_id_fm -_off_el[VV][isub*_NoLevels + lev] + el_sum_prev_sd_at_lev;
//                 std::cout << "Volume PROC " << volume_proc  << std::endl;
                _el_child_to_fath[lev][ielem] = vol_elem_per_lev;  //father_id_fm;
//                 std::cout << "ielem pos " << ielem << std::endl;
                //TODO of course _el_sto and _el_sto_b have the Node and Elem numberings still in LIBMESH order
                //given the libmesh order, we must get the FEMuS order
                //TODO is it possible that when you have ONE PROC and ONE SUBDOMAIN
                // the LIBMESH and FEMUS orderings COINCIDE?
                //I guess I will build the INVERSE of the ELEMENT ORDERING
                /*CHECK*/
                if (ielem >= _n_elements_vb_lev[BB][lev]  ) {
                    std::cout << ielem << " Something wrong with number of elems" << std::endl;
                    abort();
                }

                ielem++;

            } //end iel
        } //end subdomain
    }//end lev


    return;
}



// ================================================
// REORDERING ELEMENT (lev>subdomain)
// ================================================
// the goal is to put the elements first by subdomain, and by level inside each subdomain
// in this way you can distinguish them easily if you want to do multigrid AND/OR domain decomposition
// WHAT IS THE WAY TO DISTINGUISH THEM? You can use the sum (L+P*N_L)
//rank_sl 0, rank_sl1: you regroup them: if they have the SAME SUBDOMAIN and SAME LEVEL, you have the SAME SUM (L+P*N_L)!
// Then you SORT them and so you have the FEMUS ordering! (This would be the ordering for CELL dofs!)
// in order to get back to libmesh ordering, just do el_fm_libm[kel].second
//the element info is always taken from the elem_sto file which is in libmesh ordering
//reordering the ELEMENTS is easier than reordering the NODES
// the point is that we have a list of elements of all levels, and a list of the FINE nodes.
// understanding what subd and lev the element have is immediate (->subd ->_lev), so you can group them.
// for the nodes you have to start from the ELEMENT info.
//for every node you have to figure out: what subdomain it belongs to, what level it belongs to, and also ...
// in practice for nodes you have to reason with an auxiliary level.
// if you have N mesh levels, i.e. N quadratic levels,
// then for the dofs you have to consider N levels for the LINEAR part and N levels for the quadratic part.
// Now, since QQ and LL can be interlaced, you can consider N+1 overall levels
// the coarsest one contains only the COARSE LINEAR nodes
//then you have the COARSE QUADRATIC NODES = FINER LINEAR NODES ... and so on
//so, for the ordering we could add a level with a +1...
// PARAVIEW ORDERING: it is the ordering we fix here. For 1 processor you first see the coarse linear nodes,
// but for more processors
//el_fm_libm is a vector of pairs
//we want to make it VB, the only point is El_Sto vs El_Sto_b
// ==================================================================
void Mesh::ReorderElementBySubdLev_VV() {

    _el_fm_libm.resize(_n_elements_sum_levs[VV]);  //every component of this vector is a pair of integers  ==> 3 things
    _el_libm_fm = new int[_n_elements_sum_levs[VV]];     //so far I need the INVERSE for the VOLUME ELEMENT LIST, not for the boundary so far
    //==== Volume ==========
    for (int kel=0;kel<_n_elements_sum_levs[VV];kel++) {
        // counting (level and subdom)
        int isubd = _el_sto[kel]->_subd;
        int  ilev = _el_sto[kel]->_lev;
        int rank_sl = ilev + isubd*_NoLevels;
        _n_elements_sl_vb[VV][rank_sl]++;
        _el_fm_libm[kel].first  = rank_sl ;
        _el_fm_libm[kel].second = kel;                                 //this is the libmesh ordering
    }

    // std reordering on pairs //pairs are reordered wrt their .FIRST element
    std::sort(_el_fm_libm.begin(),_el_fm_libm.end());

    //====== INVERSE Volume Elements ===========
    for (uint i = 0; i < _el_fm_libm.size(); i++)  _el_libm_fm[_el_fm_libm[i].second]=i;


    return;

}


// ==================================================================
void Mesh::ReorderElementBySubdLev_BB() {

    _el_fm_libm_b.resize(_n_elements_sum_levs[BB]);

    for (int kel=0;kel<_n_elements_sum_levs[BB];kel++) {
        int isubdom = _el_sto_b[kel]->_subd;  // _el_sto[_el_sto_b[kel]->_vol_id]->_subd;  //this is the only difference with VV: you get the SUBDOMAIN from the VOLUME!
        int ilev    = _el_sto_b[kel]->_lev;
        int rank_sl = ilev+isubdom*_NoLevels;
        _n_elements_sl_vb[BB][rank_sl]++;
        _el_fm_libm_b[kel].first  = rank_sl;
        _el_fm_libm_b[kel].second = kel;
    }
    // std reordering on pairs
    std::sort(_el_fm_libm_b.begin(),_el_fm_libm_b.end());

    return;

}




void Mesh::FillNodeSto() {

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


// ==================================================================
void Mesh::ReorderNodesBySubdLev() {

    _n_nodes_sl_ql     = new int*[QL_NODES];
    for (int fe=0;fe < QL_NODES; fe++) {
        int offdim2 = _NoSubdom*_NoLevels+1;
        _n_nodes_sl_ql[fe] = new int[offdim2];
        for (int i=0; i<offdim2; i++)      _n_nodes_sl_ql[fe][i]=0;
    }

    _nd_fm_libm.resize(_n_nodes);                        //every component of this vector is a pair of integers (some_rank_to_give_an_order, libmesh_index)
    _nd_libm_fm = new int[_n_nodes];                       //the new numbering in the desired P-L ordering

    // -----------------------------------------------
    // REORDERING NODES (subdom>lev>levP>var)
    // -----------------------------------------------
    //since the mesh is QUADRATIC and we also want to track COARSE LINEAR nodes, it's like thinking of an ADDITIONAL LEVEL
    // 1 + LP + (NL+1)(L + NL*SUBD)
    int n_variables=1;
    int n_nodes_0p=0;
    for (int knode=0;knode<_n_nodes;knode++) {
        _n_nodes_sl_ql[QQ][_NoLevels*_nd_sto[knode]->_subd + _nd_sto[knode]->_lev]++;   //order by quadratic rank
        if (_nd_sto[knode]->_levp <_NoLevels)    _n_nodes_sl_ql[LL][_NoLevels*_nd_sto[knode]->_subd + _nd_sto[knode]->_levp]++;
        if (_nd_sto[knode]->_levp == 0)    n_nodes_0p++;
        for (int klev=_nd_sto[knode]->_lev;klev<_NoLevels;klev++) _NoNodesXLev[klev]++;  // a node contributes to its level and all the finer
        _nd_fm_libm[knode].first= _nd_sto[knode]->_var
                                  + n_variables*(_nd_sto[knode]->_levp + (_NoLevels+1)*(_nd_sto[knode]->_lev + _NoLevels*_nd_sto[knode]->_subd)); //rank
        _nd_fm_libm[knode].second=knode;         //libmesh ordering value
    }
    _NoNodesXLev[_NoLevels]= n_nodes_0p;

    // std reordering on pairs
    sort(_nd_fm_libm.begin(),_nd_fm_libm.end());   //FROM NOW ON, nodes are ordered in this way
    //if i want to retrieve info about a node,  v[knode].second gives me the libmesh node ordering... BUT,
    //the point is that _nd_sto is not used later on
    //so, nd_libm_fm RECEIVES the libmesh ordering and YIELDS the FEMUS ordering
    //while for elements it is not important to store the femus index, with nodes we have to
    // nodes must be ordered some way, they are the dofs!
    //element connectivities can be order whatever you want

    for (uint i = 0; i < _nd_fm_libm.size(); i++)  _nd_libm_fm[_nd_fm_libm[i].second]=i;  //you start numbering from zero   //TODO we could overwrite v[i].first = i, without using nd_libm_fm, we dont need the rank later on


//given this node numbering, you can later fill the coords also

// el_fm_libm[kel].second receives the femus index and yields the libmesh index
// if you want the opposite just do a v_inv_el exactly like nd_libm_fm: this might be good for cell dofs!

//nd_libm_fm[i]       receives the libmesh index and yields the femus index
// v[i].second      receives the femus index and yields the libmesh index


    return;

}

//=====================================
void Mesh::ComputeElemOffsetsBySubdLevel() {
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
void Mesh::ComputeNodeOffsetsBySubdLevel()  {
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


// ============================================
// COMPUTE max number of elements on each quadratic node
// ============================================
// ==================================================================
void Mesh::ComputeMaxElXNode() {

    _elxnode = new int[_n_nodes];

    for (int inode=0;inode<_n_nodes;inode++) _elxnode[inode]=0;
    for (int isub=0;isub<_NoSubdom;isub++)
        for (int i = _off_el[VV][isub*_NoLevels];
                 i < _off_el[VV][isub*_NoLevels + 1];i++)
            for (uint k=0;k<_elnodes[VV][QQ];k++)   {
                _elxnode[_el_sto[_el_fm_libm[i].second]->_elnds[k]]++;
            }

    _maxelxnode=0;
    for (int inode=0;inode<_n_nodes;inode++) if (_elxnode[inode]>_maxelxnode) _maxelxnode=_elxnode[inode];
    std::cout<< " Max number of element on each node = " << _maxelxnode << " \n";

    delete []_elxnode;

    return;
}



// ==============================
// COMPUTE node map for "Extended Levels"
// ==============================
// this RECEIVES a NODE in FEMUS ordering
// and YIELDS the DOF NUMBERING of ONE SCALAR variable (sometimes QQ, sometimes LL, in MG...)
// ==================================================================
void Mesh::ComputeNodeMapExtLevels() {

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


