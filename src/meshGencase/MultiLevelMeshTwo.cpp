#include "MultiLevelMeshTwo.hpp"

#include <sstream>
#include <vector>
#include <algorithm>

#include "Domain.hpp"
#include "Typedefs.hpp"
#include "Files.hpp"
#include "IO.hpp"
#include "GeomEl.hpp"
#include "GeomElTypeEnum.hpp"
#include "VBTypeEnum.hpp"
#include "FETypeEnum.hpp"

#include "paral.hpp"


namespace femus {




// ========================================================
MultiLevelMeshTwo::MultiLevelMeshTwo (const Files& files_in, const FemusInputParser<double>& map_in, const std::string mesh_file_in) :
         _files(files_in),
         _mesh_rtmap(map_in),
         _dim(map_in.get("dimension")),
         _mesh_order(map_in.get("mesh_ord"))  {

    _mesh_file.assign(mesh_file_in); 

    std::vector <std::string>  geomelem; 
    geomelem.resize(_dim);
    
    if (map_in.get("geomel_type") == HEX)  {
      geomelem[0] = "line";
      geomelem[1] = "quad";
      geomelem[2] = "hex";
    }
    else if (map_in.get("geomel_type") == TET)  {
      geomelem[0] = "line";
      geomelem[1] = "tri";
      geomelem[2] = "tet";
    }
    else if (map_in.get("geomel_type") == WEDGE)  {
       std::cout << "Wedge not supported" << std::endl; abort(); 
    }
    else if (map_in.get("geomel_type") == QUAD)  {
      geomelem[0] = "line";
      geomelem[1] = "quad";
    }
    else if (map_in.get("geomel_type") == TRI)  {
      geomelem[0] = "line";
      geomelem[1] = "tri";
    }
    else if (map_in.get("geomel_type") == LINE)  {
      geomelem[0] = "line";
    }
    else  {  std::cout << "Geom Elem not supported" << std::endl; abort();   }
    
if ( _dim == 3  && (map_in.get("geomel_type") != HEX && map_in.get("geomel_type") != TET && map_in.get("geomel_type") != WEDGE) )
        {  std::cout << "Inconsistent input file" << std::endl; abort();   }
    
if ( _dim == 2  && (map_in.get("geomel_type") != QUAD && map_in.get("geomel_type") != TRI ) )
        {  std::cout << "Inconsistent input file" << std::endl; abort();   }
    
if ( _dim == 1  && (map_in.get("geomel_type") != LINE ) )
        {  std::cout << "Inconsistent input file" << std::endl; abort();   }

    
    //     How to initialize a std::vector of classes ==================
    _GeomEl.resize(_dim);
    for (int d=0; d < _dim; d++) _GeomEl[d].reserve(QL_NODES);

    for (int d=0; d < _dim; d++) { 
    for (int ql=0;ql < QL_NODES; ql++)  { 
          GeomEl geomel_temp( geomelem[d], ql ); 
         _GeomEl[d].push_back(geomel_temp); 
       }
    }
//    End How to initialize a std::vector of classes ==================

    _iproc    = paral::get_rank();
    _NoSubdom = paral::get_size();   
    _NoLevels = GetRuntimeMap().get("nolevels");
    
    const uint mesh_ord = (uint) map_in.get("mesh_ord");
    if (mesh_ord != 0) {
        std::cout << "Linear mesh not yet implemented" << std::endl;
        abort();
    }
    if (VB != 2) {
        std::cout << " NoFamFEM is not 2: to do! ";
        abort();
    }

    for (int vb=0;vb < VB; vb++) {
        _elnodes[vb][QQ] = GetGeomEl(_dim-1-vb,mesh_ord)._elnds;  //THE MESH ORD CAN ONLY BE QUADRATIC UP TO NOW!
        _elnodes[vb][LL] = GetGeomEl(_dim-1-vb,LL)._elnds;  //THE MESH ORD CAN ONLY BE QUADRATIC UP TO NOW!
        _elnodes[vb][KK] = 1;
    }
    //i do not want to use the linear part actually!!

  
    _nodes_name = "/NODES";
    _elems_name = "/ELEMS";
//    _nd_coord_folder = "COORD";
//      _el_pid_name = "PID";
//     _nd_map_FineToLev = "MAP";


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
    delete[] _type_FEM;
    
  return;
}


// ========================================================
//This function transforms the node coordinates into the reference node coordinates
  void MultiLevelMeshTwo::TransformElemNodesToRef(const uint elem_dim,const double* xx_qnds, double* refbox_xyz) const {
   
   double*   x_in = new double[_dim];
   double*   x_out = new double[_dim];
  const uint mesh_ord = (int) GetRuntimeMap().get("mesh_ord");
  const uint el_nds = GetGeomEl(elem_dim-1,mesh_ord)._elnds;

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
  void MultiLevelMeshTwo::SetDomain(Domain* domain_in)  {
    
    _domain = domain_in;
    
   return; 
  }
  
 // ========================================================
  Domain* MultiLevelMeshTwo::GetDomain() const {
    
   return _domain; 
   
  }


// ========================================================
/// Read mesh from hdf5 file (namefile) 
//is this function for ALL PROCESSORS 
// or only for PROC==0? Seems to be for all processors
// TODO do we need the leading "/" for opening a dataset?
// This routine reads the mesh file and also makes it NONDIMENSIONAL, so that everything is solved on a nondimensional mesh
void MultiLevelMeshTwo::ReadMeshFileAndNondimensionalize()   {

  std::string    basemesh = DEFAULT_BASEMESH;
  std::string      ext_h5 = DEFAULT_EXT_H5;
  
  std::ostringstream meshname;
  meshname << _files._output_path << "/" << basemesh  << ext_h5;

//==================================
// OPEN FILE 
//==================================
  std::cout << " Reading mesh from= " <<  meshname.str() <<  std::endl;
  hid_t  file_id = H5Fopen(meshname.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
//   if (file_id < 0) {std::cout << "MultiLevelMeshTwo::read_c(): File Mesh input in data_in is missing"; abort();}
// he does not do this check, if things are wrong the H5Fopen function detects the error

//==================================
// DFLS (Dimension, VB, Levels, Subdomains)
// =====================
  uint topdata[4];
  IO::read_UIhdf5(file_id,"/DFLS",topdata);

//==================================
// CHECKS 
// ===================== 
 if (_NoLevels !=  topdata[2])  {std::cout << "MultiLevelMeshTwo::read_c. Mismatch: the number of mesh levels is " <<
   "different in the mesh file and in the configuration file" << std::endl;abort(); }


 if (_NoSubdom != topdata[3])  {std::cout << "MultiLevelMeshTwo::read_c. Mismatch: the number of mesh subdomains is " << _NoSubdom
                                   << " while the processor size of this run is " << paral::get_size()
                                   << ". Re-run gencase and your application with the same number of processors" << std::endl;abort(); }
  
//alright, there is a check to make: if you run gencase in 3D, and then you run main in 2D then things go wrong...
//this is because we read the dimension from 'dimension', we should read it from the mesh file in principle, 
//in fact it is that file that sets the space in which we are simulating...
//I'll put a check 

if (_dim != topdata[0] ) {std::cout << "MultiLevelMeshTwo::read_c. Mismatch: the mesh dimension is " << _dim
                                   << " while the dimension in the configuration file is " << GetRuntimeMap().get("dimension")
                                   << ". Recompile either gencase or your application appropriately" << std::endl;abort();}
//it seems like it doesn't print to file if I don't put the endline "<< std::endl".
//Also, "\n" seems to have no effect, "<< std::endl" must be used
//This fact doesn't seem to be related to PARALLEL processes that abort sooner than the others

if ( VB !=  topdata[1] )  {std::cout << "MultiLevelMeshTwo::read_c. Mismatch: the number of integration dimensions is " << topdata[1]
                                   << " while we have VB= " << VB 
                                   << ". Re-run gencase and your application appropriately " << std::endl;abort(); }

//==================================
//======== mesh order
//==================================
//TODO should be known by the file, or hard coded there...
//the point is: what if you already have a mesh in external format?
//then you have to either read the mesh order from the external file, if such
//information can be retrieved, or provide it explicitly
    const uint mesh_ord = (int) GetRuntimeMap().get("mesh_ord");    

//==================================
// FEM element DoF number
// =====================
//Reading this is not very useful... well, it may be a check  
  _type_FEM=new uint[VB];
  IO::read_UIhdf5(file_id, "/ELNODES_VB",_type_FEM);

  for (int vb=0; vb<VB;vb++) {
if (_type_FEM[vb] !=  GetGeomEl(_dim-1-vb,QQ)._elnds )  {std::cout << "MultiLevelMeshTwo::read_c. Mismatch: the element type of the mesh is" <<
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
  _el_map[vb]=new uint [_off_el[vb][_NoSubdom*_NoLevels]*GetGeomEl(_dim-1-vb,mesh_ord)._elnds];
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



//=============================================
//this is different from the Mesh class function
//because we are using quadratic elements
void MultiLevelMeshTwo::PrintSubdomFlagOnQuadrCells(const int vb, const int Level, std::string filename) const {

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
void MultiLevelMeshTwo::PrintForVisualizationAllLEVAllVB()  const {
  
    const uint iproc=_iproc;
   if (iproc==0) {
       PrintConnLinAllLEVAllVB();
       PrintXDMFAllLEVAllVB();
    }
   
   return;
}


// ==================================================================
void MultiLevelMeshTwo::PrintMultimeshXdmf() const {

     if (_iproc==0) {
  
    std::string multimesh = DEFAULT_MULTIMESH;
    std::string ext_xdmf  = DEFAULT_EXT_XDMF;
    std::string basemesh  = DEFAULT_BASEMESH;
    std::string ext_h5    = DEFAULT_EXT_H5;

    std::ostringstream inmesh_xmf;
    inmesh_xmf << _files._output_path << "/" << multimesh << ext_xdmf;
    std::ofstream out(inmesh_xmf.str().c_str());

    std::ostringstream top_file;
    top_file << basemesh << ext_h5;

    if (out.fail()) {
        std::cout << "MultiLevelMeshTwo::PrintMultimeshXdmf: The file is not open" << std::endl;
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
            IO::PrintXDMFTopology(out,top_file.str(),hdf5_field.str(),GetGeomEl(_dim-1-vb,QQ).name,
				  _n_elements_vb_lev[vb][ilev],
				  _n_elements_vb_lev[vb][ilev],
				  _elnodes[vb][QQ]);
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

    }    //end iproc
    
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

 void MultiLevelMeshTwo::PrintXDMFAllLEVAllVB() const {

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
void MultiLevelMeshTwo::PrintXDMFGridVB(std::ofstream& out,
			      std::ostringstream& top_file,
			      std::ostringstream& geom_file, const uint Level, const uint vb) const  {

  std::string grid_mesh[VB];
  grid_mesh[VV]="Volume";
  grid_mesh[BB]="Boundary";
  
  std::ostringstream hdf_field; hdf_field << "MSHCONN_VB_" << vb << "_LEV_" << Level;
  
    uint nel = _n_elements_vb_lev[vb][Level];

    out << "<Grid Name=\"" << grid_mesh[vb].c_str() << "_L" << Level << "\"> \n";
    
   IO::PrintXDMFTopology(out,top_file.str(),hdf_field.str(),
			     GetGeomEl(_dim-1-vb,LL).name,
			 nel*GetGeomEl(_dim-1-vb,_mesh_order).n_se,
			 nel*GetGeomEl(_dim-1-vb,_mesh_order).n_se,
			     GetGeomEl(_dim-1-vb,LL)._elnds);    

   std::ostringstream coord_lev; coord_lev << "_L" << Level; 
   IO::PrintXDMFGeometry(out,geom_file.str(),"NODES/COORD/X",coord_lev.str(),"X_Y_Z","Float",_NoNodesXLev[Level],1);
    
    out << "</Grid> \n";  

   return;
}



// ========================================================
/// It prints the connectivity in hdf5 format
/// The changes are only for visualization of quadratic FEM

void MultiLevelMeshTwo::PrintConnLinAllLEVAllVB() const { 

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



void MultiLevelMeshTwo::PrintConnLinVB(hid_t file, const uint Level, const uint vb) const {
  
   int conn[8][8];   uint *gl_conn;
  
    uint icount=0;
    uint mode = GetGeomEl(_dim-1-vb,QQ)._elnds;
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
void MultiLevelMeshTwo::PrintSubdomFlagOnLinCells(std::string filename) const {
  
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


