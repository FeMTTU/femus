#include "EquationsMap.hpp"

// C++
#include <iomanip>
#include <sstream>

// local includes 
#include "FemusDefault.hpp"

#include "Files.hpp"
#include "IO.hpp"
#include "Physics.hpp"
#include "Quantity.hpp"
#include "MeshTwo.hpp"
#include "GeomEl.hpp"
#include "FEElemBase.hpp"
#include "QRule.hpp"
#include "TimeLoop.hpp"

#include "paral.hpp"

// ====================================================
/// This function constructs the equation map

EquationsMap::EquationsMap(Files& files_in,         // Utils pointer
                           Physics& mgphys_in,        // Physics pointer
                           QuantityMap& qtymap_in,
                           Mesh& mgmesh_in,           // Mesh pointer
                           std::vector<FEElemBase*>& absfe_in,
			   QRule&   qrule_in,
                           TimeLoop& timeloop_in ):
        _files(files_in),
        _phys(mgphys_in),
        _qtymap(qtymap_in),
        _mesh(mgmesh_in),
        _AbstractFE(absfe_in),
        _qrule(qrule_in),
        _timeloop(timeloop_in)  {}


// ====================================================
/// This function destroys the equations
void EquationsMap::clean() {
    for (EquationsMap::iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++) {
        delete eqn->second;
    }
}

// ====================================================
/// This sets dof initial and boundary conditions and sets the operators
//inside this functions there are a lot of new of class elements
//so, i must do a corresponding clean NOT in the DESTRUCTOR
//what is new'd in the constructor is deleted in the destructor
//what is new'd with this function is deleted with a corresponding clear function
//all these things are related to the Dofs, so the first thing is to settle them
// ===============================================================
/// This  function reads all the Operators from files
///  initialization of all levels: dofs and matrices;
/// initialize dof map
/// initialize BC
/// initialize MG Operators (TODO the values of the Restrictor Operators DEPEND on the boundary conditions)
/// initialize vectors (could do this even before the BC)
/// The INITIAL conditions can be done only after initVectors();

void  EquationsMap::setDofBcOpIc() {

    for (iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++) {
        EqnBase* mgsol = eqn->second;
        
#ifdef DEFAULT_PRINT_INFO
    std::cout << "\n Reading "  <<  mgsol -> _eqname << " Dof, Bc, Op, Ic \n";
#endif

//=====================
    mgsol -> ComputeMeshToDof();
//=====================
    mgsol -> GenBc();
    mgsol -> GenElBc();
//=====================
    mgsol -> initMGOps();
//=====================
    mgsol -> initVectors();     //TODO can I do it earlier than this position?
//=====================
    mgsol -> GenIc();              // initial solution

#ifdef DEFAULT_PRINT_INFO
    std::cout << " Dof, Bc, Op, Ic settled for"  <<  mgsol -> _eqname <<  "\n";
#endif
    }
    
    return;
}

// ==========================================================================================
/// This function performes all the Physics time step routines
void EquationsMap::OneTimestepEqnLoop(
    const double time,             // real time
    const uint delta_t_step_in     // integer time
) {
    // loop for time steps
    for (iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
        EqnBase* mgsol = eqn->second;
        mgsol -> MGTimeStep(time,delta_t_step_in);
    }
    return;
}

// =================================================================
/// This function prints xdmf and hdf5 file
void EquationsMap::PrintSol(const uint t_step, const double curr_time) {

    PrintSolHDF5(t_step);
    PrintSolXDMF(t_step,curr_time);

    return;
}

// =================================================================
/// This function prints the attributes into the corresponding hdf5 file
void EquationsMap::PrintSolHDF5(const uint t_flag ) {

    const uint    iproc =_mesh._iproc;
    if (iproc==0) {

        const uint ndigits  = _timeloop._timemap.get("ndigits");

        std::string    basepath = _files.get_basepath();
        std::string output_dir  = DEFAULT_OUTPUTDIR;
        std::string outtime_dir = _files.get_frtmap().get("OUTTIME_DIR");
        std::string    basesol  = DEFAULT_BASESOL;
        std::string     ext_h5  = DEFAULT_EXT_H5;
        std::ostringstream filename;
        filename << basepath << "/" << output_dir << outtime_dir
        << basesol << "." << setw(ndigits) << setfill('0') << t_flag << ext_h5;

        hid_t   file= H5Fcreate(filename.str().c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
        H5Fclose(file);

        EquationsMap::const_iterator pos=_equations.begin();
        EquationsMap::const_iterator pos_e=_equations.end();
        for (;pos!=pos_e;pos++)    {
            EqnBase *mgsol=pos->second;
            mgsol->PrintVector(filename.str());
        }

    } //end print iproc

    return;
}

// ===================================================================
/// It prints the attributes in  Xdmf format for one time step
// Qui dobbiamo separare le NODE variables dalle CELL variables!
//le node variables sono le _nvars[QQ] e _nvars[LL], le Cell sono _nvars[KK]
//The SCALAR variables are ordered FIRST QUADRATIC, then LINEAR, then CONSTANT

//Ok, there is a problem when I want to print on coarser grids,
// it seems that it needs the coarse coordinates
//Ok i think i solved it

//Now the thing is: i need to print from FINEST LEVEL, so that i can see in paraview the FINEST sequence!!!
//Ok, now I want to print an XDMF for every level, because I've seen that paraview has quite a little bit
// of problems in reading files with multiple grids, so I'll put everything in separate files per level.
//The file name depends on the level so i will put it inside the level loop

// If you want to put all the grids in the SAME XDMF file, then it is better to put the FINEST first, and then all the others...
// This is due to the Paraview TIME Files that would read the FIRST GRID THEY MEET in the XDMF FILE!
//So you should do  in reverse order:  for (int l=NoLevels-1; l >= 0; l--) 
// I DO NOT WANT TO SEPARATE the HDF5 file, because it works fine as a single file with no problems

void EquationsMap::PrintSolXDMF(const uint t_step,const double curr_time) {

    const uint    iproc =_mesh._iproc;
    if (iproc==0) {

      const uint ndigits  = _timeloop._timemap.get("ndigits");

      const uint NoLevels = _mesh._NoLevels;

        //FE print
        std::string DofType[QL];
        DofType[QQ] = "Node";
        DofType[LL] = "Node";
        DofType[KK] = "Cell";

        std::string    basepath = _files.get_basepath();
        std::string output_dir  = DEFAULT_OUTPUTDIR;
        std::string outtime_dir = _files.get_frtmap().get("OUTTIME_DIR");
        std::string basesol     = DEFAULT_BASESOL;
        std::string basemesh    = DEFAULT_BASEMESH;
        std::string aux_xdmf    = DEFAULT_AUX_XDMF;
        std::string connlin     = DEFAULT_CONNLIN;
        std::string     ext_h5  = DEFAULT_EXT_H5;
        std::string    ext_xdmf = DEFAULT_EXT_XDMF;

// =================================
// ============= LEVELS ============
// =================================
	
	for (uint l=0; l < NoLevels; l++) {
	
	std::ostringstream filename_xdmf;
        filename_xdmf << basepath << "/" << output_dir << outtime_dir
        << basesol << "." << setw(ndigits) << setfill('0') << t_step << "_l" << l << ext_xdmf;
        std::ostringstream hdf_file; hdf_file << basesol << "." << setw(ndigits) << setfill('0') << t_step << ext_h5;

	std::ofstream out(filename_xdmf.str().c_str());
        if (out.fail()) {
            std::cout << "EquationsMap::print_soln_xmf: cannot print " << filename_xdmf.str().c_str() << std::endl;
            abort();
        }

        //  ++++++++++++ Header ++++++++++++++
        out << "<?xml version=\"1.0\" ?> \n";
        out << "<!DOCTYPE Xdmf SYSTEM ";
        out <<  "\"" << basepath << "/" << output_dir << outtime_dir << "/" << aux_xdmf << "\" \n";
//   out << " [ <!ENTITY HeavyData \"\"> ] ";
        out << ">\n";
        out << "<Xdmf> \n" << "<Domain> \n";
	

        int NGeomObjOnWhichToPrint[QL];
        NGeomObjOnWhichToPrint[QQ] = _mesh._NoNodesXLev[l];
        NGeomObjOnWhichToPrint[LL] = _mesh._NoNodesXLev[l];
        NGeomObjOnWhichToPrint[KK] = _mesh._NoElements[VV][l]*_mesh._GeomEl.n_se[VV];
	  
	out << "<Grid Name=\"Volume_L" << l << "\"> \n";

	out << "<Time Value =\"" << curr_time << "\" /> \n";

	PrintXDMFTopologyGeometry(out,l,VV);

	EquationsMap::const_iterator pos1   = _equations.begin();
        EquationsMap::const_iterator pos1_e = _equations.end();
        for (;pos1!=pos1_e;pos1++)   {
            EqnBase *mgsol=pos1->second;
            int OffVarNames[QL];
            OffVarNames[QQ] = 0;
            OffVarNames[LL] = mgsol->_nvars[QQ];
            OffVarNames[KK] = mgsol->_nvars[QQ] + mgsol->_nvars[LL];
            for (int fe=0; fe<QL; fe++)  {
                for (uint ivar=0;ivar< mgsol->_nvars[fe]; ivar++)   {
                   std::ostringstream var_name; var_name << mgsol->_var_names[ OffVarNames[fe] + ivar] << "_LEVEL" << l;
                   IO::PrintXDMFAttribute(out,hdf_file.str(),var_name.str(),var_name.str(),"Scalar",DofType[fe],"Float",NGeomObjOnWhichToPrint[fe],1);
                }
            } // end fe
        }

        out << "</Grid>\n";
	
        // ============= END GRID ===========	

	out << "</Domain> \n" << "</Xdmf> \n";
        out.close();
	
	}//end levels
	
    } //end iproc==0

    return;
}

// ========================================================================
/// This function read the solution form all the system (restart)
void EquationsMap::ReadSol(const uint t_step, double& time_out) {

    const uint ndigits      = _timeloop._timemap.get("ndigits");
    std::string    basepath = _files.get_basepath();
    std::string output_dir  = DEFAULT_OUTPUTDIR;
    std::string outtime_dir = _files.get_frtmap().get("OUTTIME_DIR");
    std::string    basesol  = DEFAULT_BASESOL;
    std::string   ext_xdmf  = DEFAULT_EXT_XDMF;
    std::string     ext_h5  = DEFAULT_EXT_H5;
// ---------------------------------------------------
    // reading time from from sol.N.xmf file
    // ---------------------------------------------------
    // open file -----------------------------
    std::ostringstream namefile;
    namefile << basepath << "/" << output_dir << outtime_dir
    << basesol << "." << setw(ndigits) << setfill('0') << t_step << "_l" << (_mesh._NoLevels - 1) << ext_xdmf;  //TODO here we should avoid doing this process TWICE because we already do it in the TransientSetup calling function

#ifdef DEFAULT_PRINT_INFO // --------  info ------------------ 
    std::cout << "\n EquationsMap::read_soln: Reading time  from "
              << namefile.str().c_str();
#endif  // -------------------------------------------
    std::ifstream in ;
    in.open(namefile.str().c_str());  //associate the file stream with the name of the file
    if (!in.is_open()) {
        std::cout << " read_soln: restart .xmf file not found "  << std::endl;
        abort();
    }

    // reading time from xmf file --------------
    std::string buf="";
    while (buf != "<Time") in >> buf;
    in >> buf >> buf;
    buf=buf.substr(2,buf.size()-3);
//create an istringstream from a string
    std::istringstream buffer(buf);
    double restart_time;
    buffer >> restart_time;

    //pass  the time value to the calling function
    time_out = restart_time;

//add parameter to system dont need that now
//   _utils.set_par("restartime",restart_time);

    // ---------------------------------------------------
    // reading data from  sol.N.h5
    // ---------------------------------------------------
    // file name -----------------------------------------
    namefile.str("");  //empty string
    namefile << basepath << "/" << output_dir << outtime_dir
    << basesol << "." << setw(ndigits) << setfill('0') << t_step << ext_h5;
    //if i put the path of this file to be relative, will the read depend on where I launched the executable...
    // or where the executable is I think... no, the path is given by where the executable is LAUNCHED

#ifdef DEFAULT_PRINT_INFO  // --------------- info ---------------
    std::cout << "\n EquationsMap::read_soln: Reading from file "
              << namefile.str().c_str() << std::endl;
#endif // ---------------------------------------------
    // loop reading over the variables ---------------------
    for (EquationsMap::const_iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
        EqnBase *mgsol=eqn->second;
        mgsol->ReadVector(namefile.str());
    } //  loop --------------------------------------------------------

    return;
}

// ====================================================
/// This function prints initial and boundary data in xdmf+hdf5 format
// of course whenever you change the fields printed in the case h5 file
// then you need to change also the xdmf file
// we should do a routine that for a given field prints both the hdf5 dataset
// and the  xdmf tag... well it's not so automatic, because you need to know
// what is the grid on which to print, bla bla bla
void EquationsMap::PrintCase(const uint t_init) {
  
     PrintCaseHDF5(t_init);
     PrintCaseXDMF(t_init);

    return;
} 


// =============================================================================
/// This function prints initial and boundary data in hdf5 fromat
/// in the file case.h5
void EquationsMap::PrintCaseHDF5(const uint t_init) {

    const uint    iproc =_mesh._iproc;
    if (iproc==0) {

        const uint ndigits      = _timeloop._timemap.get("ndigits");
        std::string    basepath = _files.get_basepath();
        std::string output_dir  = DEFAULT_OUTPUTDIR;
        std::string outtime_dir = _files.get_frtmap().get("OUTTIME_DIR");
        std::string    basecase = DEFAULT_BASECASE;
        std::string   ext_xdmf  = DEFAULT_EXT_XDMF;
        std::string     ext_h5  = DEFAULT_EXT_H5;

        std::ostringstream filename;
        filename << basepath << "/" << output_dir << outtime_dir
        << basecase << "." << setw(ndigits) << setfill('0') << t_init << ext_h5;

        hid_t file = H5Fcreate(filename.str().c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        _mesh.PrintSubdomFlagOnLinCells(filename.str());
        H5Fclose(file);

        EquationsMap::const_iterator pos   = _equations.begin();
        EquationsMap::const_iterator pos_e = _equations.end();
        for (;pos!=pos_e;pos++) {
            EqnBase *mgsol=pos->second;
            mgsol->PrintVector(filename.str());          // initial solution
            mgsol->PrintBc(filename.str());            // boundary condition
        }

    } //end iproc

    return;
}



// ====================================================================
/// It prints the Xdmf file to read the initial and boundary conditions
//now, this function does exactly the same as print sol;
//more over, it prints the PID, and the BOUNDARY CONDITIONS
//so it prints, INITIAL CONDITIONS, BOUNDARY CONDITIONS, PID
//let us split so that we can have a unique function

//clearly, you need to know WHERE to print this file.
//so you need the absolute paths
//but, inside the lines of this file, you dont need to put the absolute paths,
//because you already know you'll not separate .xmf and .h5

void EquationsMap::PrintCaseXDMF(const uint t_init) {

    const uint    iproc =_mesh._iproc;
    if (iproc==0) {

        const uint NoLevels = _mesh._NoLevels;
        const uint ndigits  = _timeloop._timemap.get("ndigits");

        std::string     basepath = _files.get_basepath();
        std::string    input_dir = DEFAULT_CASEDIR;
        std::string   output_dir = DEFAULT_OUTPUTDIR;
        std::string  outtime_dir = _files.get_frtmap().get("OUTTIME_DIR");
        std::string     basecase = DEFAULT_BASECASE;
        std::string     basemesh = DEFAULT_BASEMESH;
        std::string       ext_h5 = DEFAULT_EXT_H5;
        std::string     ext_xdmf = DEFAULT_EXT_XDMF;
        std::string     aux_xdmf = DEFAULT_AUX_XDMF;
        std::string      connlin = DEFAULT_CONNLIN;
        std::string  bdry_suffix = DEFAULT_BDRY_SUFFIX;

        //FE print
        std::string DofType[QL];
        DofType[QQ] = "Node";
        DofType[LL] = "Node";
        DofType[KK] = "Cell";

        std::string var_name[VB];
        std::string var_type[VB];
	
// =================================
// ============= LEVELS ============
// =================================
	
	for (uint l=0; l < NoLevels; l++) {
	  
        std::ostringstream filename_xdmf;
        filename_xdmf << basepath << "/" << output_dir << outtime_dir
        << basecase  << "." << setw(ndigits) << setfill('0') << t_init << "_l" << l << ext_xdmf;
        std::ostringstream hdf_file;
        hdf_file <<  basecase << "." << setw(ndigits) << setfill('0') << t_init << ext_h5;

        std::ofstream out(filename_xdmf.str().c_str());
        if (out.fail()) {
            std::cout << "EquationsMap::print_case_xmf: cannot print " << filename_xdmf.str().c_str() << std::endl;
            abort();
        }

        // BEGIN XDMF =======
        out << "<?xml version=\"1.0\" ?> \n";
        out << "<!DOCTYPE Xdmf SYSTEM ";
        out <<  "\"" << basepath << output_dir << outtime_dir << "/" << aux_xdmf << "\" \n";
//    out << " [ <!ENTITY HeavyData \"\"> ] ";
        out << ">\n";
        out << "<Xdmf> \n" << "<Domain> \n";
	

        int NGeomObjOnWhichToPrint[QL];
        NGeomObjOnWhichToPrint[QQ] = _mesh._NoNodesXLev[l];
        NGeomObjOnWhichToPrint[LL] = _mesh._NoNodesXLev[l];
        NGeomObjOnWhichToPrint[KK] = _mesh._NoElements[VV][l]*_mesh._GeomEl.n_se[VV];

	out << "<Grid Name=\"Volume_L" << l << "\"> \n";

        // TOPOLOGY GEOMETRY ===========
        PrintXDMFTopologyGeometry(out,l,VV);

	// ===== PID ======
        std::ostringstream  pid_name; pid_name << "PID" << "_LEVEL" << l;
	IO::PrintXDMFAttribute(out,hdf_file.str(),pid_name.str(),pid_name.str(),"Scalar",DofType[KK],"Int",NGeomObjOnWhichToPrint[KK],1);

        // ATTRIBUTES FOR EACH SYSTEM ===========
        EquationsMap::const_iterator pos1=_equations.begin();
        EquationsMap::const_iterator pos1_e=_equations.end();
        for (;pos1!=pos1_e;pos1++)   {
            EqnBase *mgsol=pos1->second;
            int OffVarNames[QL];
            OffVarNames[QQ] = 0;
            OffVarNames[LL] = mgsol->_nvars[QQ];
            OffVarNames[KK] = mgsol->_nvars[QQ] + mgsol->_nvars[LL];
            for (int fe=0; fe<QL; fe++)  {
                for (uint ivar=0; ivar < mgsol->_nvars[fe]; ivar++)     {
		    std::ostringstream  varstream; varstream << mgsol->_var_names[OffVarNames[fe] + ivar] << "_LEVEL" << l;
                    var_name[VV] = varstream.str();
                    var_type[VV] = "Float";
                    var_name[BB] = var_name[VV] + bdry_suffix;
                    var_type[BB] = "Int";
                    for (int vb=0;vb<VB; vb++) {
                        IO::PrintXDMFAttribute(out,hdf_file.str(),var_name[vb],var_name[vb],"Scalar",DofType[fe],var_type[vb],NGeomObjOnWhichToPrint[fe],1);
                    }
                }
            } //end fe

        }  //end eqn

        out << "</Grid>\n";

	out << "</Domain> \n" << "</Xdmf> \n";
        out.close();
	
	} //end levels

    } //if iproc==0
    
    return;
}


//print topology and geometry, useful for both case.xmf and sol.xmf
void EquationsMap::PrintXDMFTopologyGeometry(std::ofstream& out, const uint Level, const uint vb) {

    //Mesh
    uint n_elements = _mesh._NoElements[vb][Level];

    std::string basemesh   = DEFAULT_BASEMESH;
    std::string connlin    = DEFAULT_CONNLIN;
    std::string     ext_h5 = DEFAULT_EXT_H5;
    
    //connectivity
    std::ostringstream connfile; connfile <<  basemesh << connlin <<  ext_h5;
    std::ostringstream hdf5_field; hdf5_field << "MSHCONN_VB_" << vb << "_LEV_" << Level;
    //coordinates
    std::ostringstream coord_file; coord_file <<  basemesh <<  ext_h5;
    
    IO::PrintXDMFTopology(out,connfile.str(),hdf5_field.str(),_mesh._GeomEl.pname[vb],n_elements*_mesh._GeomEl.n_se[vb],n_elements*_mesh._GeomEl.n_se[vb],_mesh._GeomEl._elnds[vb][LL]);
    std::ostringstream coord_lev; coord_lev << "_L" << Level; 
    IO::PrintXDMFGeometry(out,coord_file.str(),"NODES/COORD/X",coord_lev.str(),"X_Y_Z","Float",_mesh._NoNodesXLev[Level],1);

    return;
}

////////////////////////////////
////////////////////////////////

void EquationsMap::TransientSetup()  {

    const uint restart  = _timeloop._timemap.get("initial_step");
    const uint ndigits  = _timeloop._timemap.get("ndigits");

    std::string    basepath = _files.get_basepath();
    std::string  output_dir = DEFAULT_OUTPUTDIR;
    std::string outtime_dir = _files.get_frtmap().get("OUTTIME_DIR");
    std::string   lastrun_f = DEFAULT_LAST_RUN;
    std::string     basesol = DEFAULT_BASESOL;
    std::string    ext_xdmf = DEFAULT_EXT_XDMF;
    std::string      ext_h5 = DEFAULT_EXT_H5;

    std::string    basecase = DEFAULT_BASECASE;
    std::string    basemesh = DEFAULT_BASEMESH;

    std::string  aux_xdmf   = DEFAULT_AUX_XDMF;
    std::string  connlin    = DEFAULT_CONNLIN;


    std::string lastrun_str;
    lastrun_str = basepath + "/" + output_dir + "/" + lastrun_f;

//now, every run, restart or not, has a new output dir.
//So, if you restart, you have to copy sol.N.h5 and sol.N.xmf
//to the new output directory
//just a raw copy, nothing more, because the file paths in sol.N.xmf
// DO NOT DEPEND ON THE OUTPUTDIR, only on the INPUT_DIR for now.
//the problem is that you have to know the PREVIOUS output_dir
// to automatically copy the files...
//let us just try by hand now...no we cant, because we will not know
// the NEW output_dir...
//We have to find a way to keep track of the PREVIOUS output dir.
//--- use a shell variable. $femus_last_run
//SORRY, BUT YOU CANT DO THAT!
// A process cannot export to parent processes.
// So the solution will be to print a very small file, called
// femus_last_run,that contains the name of the last output dir.

//We know the NEW one at this point because we constructed it.
// Well, we put them in the main output and read from there.
//So you see once more that it is important to make the difference
// between the BASE_OUTPUT and the OVERALL_OUTPUT paths
//We will distinguish them later.

    //------initial data
    if (restart) {
        _timeloop._t_idx_in  = restart;
        std::cout << "We wish to restart from time step " << _timeloop._t_idx_in << std::endl;
        std::cout << "\n *+*+* TimeLoop::transient_setup: RESTART  " << std::endl;


        if (paral::get_rank() == 0) {

            std::cout << " Reading the  output dir to from run_to_restart_from file (fill it with what you want)" << std::endl;
            //check if last_run is there
            std::string lastone;
            std::ifstream last_run;
            last_run.open(lastrun_str.c_str());
            if (!last_run.is_open()) {
                std::cout << "There is no last_run file" << std::endl;
                abort();
            }
            //read from last_run
            last_run >> lastone >> lastone;  //"run_to_restart_from" is the STRING i use in the file, so there is one intermediate "buffer"...
            //AAA output is there
            stringstream tidxin;
            tidxin << setw(ndigits) << setfill('0') << _timeloop._t_idx_in;
            std::cout << " Restarting from run: " << lastone << std::endl;
	    std::ostringstream cp_src_xmf_stream;
	    cp_src_xmf_stream << basepath << "/" << output_dir << "/" << lastone << "/" << basesol << "." << tidxin.str() << "_l" << (_mesh._NoLevels - 1) << ext_xdmf;
            std::string cp_src_xmf = cp_src_xmf_stream.str();
            fstream file_cp_xmf(cp_src_xmf.c_str());
            if (!file_cp_xmf.is_open()) {
                std::cout <<"No xmf file"<< std::endl;
                abort();
            }

	    std::ostringstream cp_src_h5_stream;
	    cp_src_h5_stream << basepath << "/" << output_dir << "/" << lastone << "/" << basesol << "." << tidxin.str() /*<< "_l" << (_mesh._NoLevels - 1)*/ << ext_h5;
            std::string cp_src_h5 = cp_src_h5_stream.str();
            fstream file_cp_h5(cp_src_h5.c_str());
            if (!file_cp_h5.is_open()) {
                std::cout <<"No h5 file"<< std::endl;
                abort();
            }

            std::string cp_dest_dir =  basepath + "/" + output_dir + "/" +  outtime_dir;

            std::string cp_cmd_xmf = "cp " + cp_src_xmf  + " " +  cp_dest_dir;
            std::string cp_cmd_h5  = "cp " + cp_src_h5   + " " +  cp_dest_dir;

            std::cout << "Copying the two restart files to the NEW output directory created before" << std::endl;
            //you should first check that sol.N.xmf and sol.N.h5 are there

            std::cout <<cp_cmd_xmf	 << std::endl;
            std::cout <<cp_cmd_h5  << std::endl;
            system(cp_cmd_xmf.c_str() );
            system(cp_cmd_h5.c_str() );

        }

//here you should wait so that you are sure that sol.N.xmf and sol.N.h5 have been copied to the new directory
        std::cout << "***** Barrier so that all the processors have the sol.N.xmf and sol.N.h5 to read from" << std::endl;
#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif

//now the question arises: if to restart,you need to re-read the sol. files,
//you will also need to
// - COPY the mesh.h5 and mesh_conn_lin.h5 (you can avoid regenerating them) (and also mesh.xmf to be able to read them)
// - COPY the MG OPERATORS (if you have them in your run) if you keep the same number of levels
//The rule is: what is connected to the PROBLEM to be solved (DOMAIN, BC's) should be:
        // LEVEL INDEPENDENT
        //PROCESSOR INDEPENDENT
//An appropriate RESTART of the solution of a problem should be
//INDEPENDENT of the METHOD to SOLVE the PROBLEM (multigrid or not, parallel computing or not...)
//So the files for the DOMAIN and the BOUNDARY CONDITIONS should reflect that.
// e.g. consider mesh.h5 and case.h5
//given some FINE DISCRETIZATION (I dont consider the case of RESTARTING with a DIFFERENT FINE discretization...
//then during the run one might do adaptive mesh, but this is another thing... ),
//some parts of them are good for a SPECIFIC NOLEVELS or a SPECIFIC NO_PROCESSORS
//while others, IN PARTICULAR THOSE THAT ARE NECESSARY FOR RESTART, must be INDEPENDENT OF THAT.

        ReadSol(_timeloop._t_idx_in,_timeloop._time_in); //read  sol.N.h5 and sol.N.xmf
        //AAA: here _t_idx_in is intent in, _time_in is intent out
        //reading files can be done in parallel with no problems
        //well, not really... reading can be done if you are sure that the file is there at the moment you call to read it
        //so let's put this inside procid=0 .... NO, THIS READ IS DONE IN PARALLEL!

    }


    else {
        //Set up initial time step index and time value
        _timeloop._t_idx_in = 0;                            //time step index
        _timeloop._time_in = 0.;                           //time absolute value

        PrintSol(_timeloop._t_idx_in,_timeloop._time_in);  //print sol.0.h5 and sol.0.xmf
        //AAA: here _t_idx_in is intent-in, and also _time_in is intent-in
    }

    std::cout << "\nInitial time index: " << _timeloop._t_idx_in
              << ", Initial time value: " << _timeloop._time_in << std::endl;

    //-------final data
    const int nsteps       = _timeloop._timemap.get("nsteps");
    const double dt        = _timeloop._timemap.get("dt");

    _timeloop._t_idx_final = _timeloop._t_idx_in + nsteps;
    _timeloop._time_final  = _timeloop._time_in + nsteps*dt;

    std::cout << "\nFinal time index: " << _timeloop._t_idx_final
              << ", Final time value: " << _timeloop._time_final
              << std::endl;

    //now you can update last_run with new_run for a following run
    //well,actually before putting the last_run you should be sure that this run was completely finished.
    //That is why I'd better put this call at the end of the main program
//   _utils._files.PrintRunForRestart(DEFAULT_LAST_RUN);

//------- print
    //this happens when the output dir is already set
    //at this point this is already true
    PrintCase(_timeloop._t_idx_in);       //print caseN.xmf&h5 = IC + BC flags
    _timeloop.transient_print_xmf(_timeloop._t_idx_in,_timeloop._t_idx_final,_mesh._NoLevels); //print timeN.xmf

    return;
}

///////////////////////////////////////////////////////
/*standard time loop over the map equations.
 The equations are solved in alphabetical order given by the map*/
void EquationsMap::TransientLoop()  {

    //  parameters
    double         dt = _timeloop._timemap.get("dt");
    int print_step =    _timeloop._timemap.get("printstep");

    double curr_time = _timeloop._time_in;  //initialize current time

    for (uint curr_step = _timeloop._t_idx_in + 1; curr_step <= _timeloop._t_idx_final; curr_step++) {

        curr_time += dt;

#if DEFAULT_PRINT_TIME==1 // only for cpu time check --------
        std::clock_t  start_time=std::clock();
#endif // ------------------------------------------- 

        // set up the time step
        std::cout << "\n  ** Solving time step " << curr_step
                  << ", time = "                 << curr_time  << " ***" << std::endl;

        const uint delta_t_step = curr_step -_timeloop._t_idx_in;


        //  time step for each system, without printing (good)
        OneTimestepEqnLoop(curr_time, delta_t_step);

#if DEFAULT_PRINT_TIME==1 // only for cpu time check --------
        std::clock_t    end_time=std::clock();
#endif  // ------------------------------------------

        // print solution
        if (delta_t_step%print_step == 0) PrintSol(curr_step,curr_time);   //print sol.N.h5 and sol.N.xmf


#if DEFAULT_PRINT_TIME==1 // only for cpu time check --------
        std::clock_t    end_time2=std::clock();
        std::cout <<" Time solver ----->= "   << double(end_time- start_time)/ CLOCKS_PER_SEC
                  <<" Time printing ----->= " << double(end_time2- end_time) / CLOCKS_PER_SEC <<
                  std::endl;
#endif  // ------------------------------------------


    }   // end time loop

    return;
}

