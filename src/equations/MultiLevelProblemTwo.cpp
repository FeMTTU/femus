/*=========================================================================

 Program: FEMUS
 Module: MultiLevelProblemTwo
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "MultiLevelProblemTwo.hpp"

// C++
#include <iomanip>
#include <sstream>

// local includes 
#include "FemusDefault.hpp"

#include "Files.hpp"
#include "XDMFWriter.hpp"
#include "XDMFWriter.hpp"
#include "Quantity.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "TimeLoop.hpp"

#include "paral.hpp"


namespace femus {



// ====================================================
/// This function constructs the equation map

MultiLevelProblemTwo::MultiLevelProblemTwo(Files& files_in,
                           FemusInputParser<double> & phys_in,
                           QuantityMap& qtymap_in,
                           MultiLevelMeshTwo& mesh_in,
                           std::vector< std::vector<elem_type*> >  & elem_type_in,
			   std::vector<Gauss>   qrule_in ):
        _files(files_in),
        _phys(phys_in),
        _qtymap(qtymap_in),
        _mesh(mesh_in),
        _elem_type(elem_type_in),
        _qrule(qrule_in)  {}


// ====================================================
/// This function destroys the equations
void MultiLevelProblemTwo::clean() {
    for (MultiLevelProblemTwo::iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++) {
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

void  MultiLevelProblemTwo::setDofBcOpIc() {

    for (iterator eqn = _equations.begin(); eqn != _equations.end(); eqn++) {
        SystemTwo* mgsol = eqn->second;
        
#ifdef DEFAULT_PRINT_INFO
    std::cout << "\n Reading "  <<  mgsol -> _eqname << " Dof, Bc, Op, Ic \n";
#endif

//=====================
    mgsol -> _dofmap.ComputeMeshToDof();
//=====================
    mgsol -> GenerateBdc();
    mgsol -> GenerateBdcElem();
//=====================
    mgsol -> ReadMGOps();
//=====================
    mgsol -> initVectors();     //TODO can I do it earlier than this position?
//=====================
    mgsol -> Initialize();              // initial solution

#ifdef DEFAULT_PRINT_INFO
    std::cout << " Dof, Bc, Op, Ic settled for"  <<  mgsol -> _eqname <<  "\n";
#endif
    }
    
    return;
}






// ========================================================================
/// This function read the solution form all the system (restart)
void MultiLevelProblemTwo::ReadSol(const uint t_step, double& time_out) const {

    const uint ndigits      = DEFAULT_NDIGITS;
    std::string    basesol  = DEFAULT_BASESOL;
    std::string   ext_xdmf  = DEFAULT_EXT_XDMF;
    std::string     ext_h5  = DEFAULT_EXT_H5;
// ---------------------------------------------------
    // reading time from from sol.N.xmf file
    // ---------------------------------------------------
    // open file -----------------------------
    std::ostringstream namefile;
    namefile << _files.GetOutputPath() << "/" 
    << basesol << "." << setw(ndigits) << setfill('0') << t_step << "_l" << (_mesh._NoLevels - 1) << ext_xdmf;  //TODO here we should avoid doing this process TWICE because we already do it in the TransientSetup calling function

#ifdef DEFAULT_PRINT_INFO // --------  info ------------------ 
    std::cout << "\n MultiLevelProblemTwo::read_soln: Reading time  from "
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
    namefile << _files.GetOutputPath() << "/"
    << basesol << "." << setw(ndigits) << setfill('0') << t_step << ext_h5;
    //if i put the path of this file to be relative, will the read depend on where I launched the executable...
    // or where the executable is I think... no, the path is given by where the executable is LAUNCHED

#ifdef DEFAULT_PRINT_INFO  // --------------- info ---------------
    std::cout << "\n MultiLevelProblemTwo::read_soln: Reading from file "
              << namefile.str().c_str() << std::endl;
#endif // ---------------------------------------------
    // loop reading over the variables ---------------------
    for (MultiLevelProblemTwo::const_iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
        SystemTwo *mgsol=eqn->second;
    } //  loop --------------------------------------------------------

    return;
}


} //end namespace femus
