#include "EqnBase.hpp"

// C++ 
#include <sstream>
#include <limits>

#include "FEMTTUConfig.h"
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "FemusDefault.hpp"

#include "NormTangEnum.hpp"
#include "QTYnumEnum.hpp"
#include "Quantity.hpp"
#include "Utils.hpp"
#include "Physics.hpp"
#include "MeshTwo.hpp"
#include "GeomEl.hpp"
#include "EquationsMap.hpp"
#include "QRule.hpp"
#include "FEElemBase.hpp"
#include "Files.hpp"
#include "CurrElem.hpp"

#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "LinearSolverM.hpp"
#include "DenseMatrix.hpp"



//the most important things for an EqnBase are:
//the number of variables
//the names
//other stuff but let us stop here now
EqnBase::EqnBase(std::vector<Quantity*> int_map_in,
                 EquationsMap& e_map_in,
                 std::string eqname_in,
                 std::string varname_in):
        _files(e_map_in._files),
        _phys(e_map_in._phys),
        _mesh(e_map_in._mesh),
        _AbstractFE(e_map_in._AbstractFE),
        _eqnmap(e_map_in),
        //=============
        _eqname(eqname_in),
        _NoLevels(e_map_in._mesh._NoLevels), //you can do that
        _var_names(NULL),
        _refvalue(NULL) {

//============= init Quantities ================
//internal std::vector
//the equal puts the two equal to each other
//as an alternative, since std::vector is a class, i could have used the initialization list above
//like i do with REFERENCES. In that case the copy constructor would be called.
    _QtyInternalVector = int_map_in;

//============= init n_vars================
    initNVars();
//========== varnames and refvalue ==============
    initVarNamesRefValues(varname_in);

    //========== processor number ==============
    _iproc=_mesh._iproc;

// ========= PENALTY DIRICHLET FLAG ==============
//put a default to zero, then every Eqn will OVERRIDE it
    _Dir_pen_fl = 0;

// ========= ELEM BC AUX ==============
 _number_tang_comps[0] = 0;
 _number_tang_comps[1] = 1;
 _number_tang_comps[2] = 3;
    
    //========= solver package ===========
    _solver = new LinearSolverM*[_NoLevels];     //well, we clearly use the same package for all levels...
    for (uint l=0;l<_NoLevels;l++) _solver[l] = LinearSolverM::build().release();
    

}



// ===================================================
/// This function  is the EqnBase destructor.
//the important thing is that these destructions occur AFTER
//the destructions of the levels inside
//these things were allocated and filled in various init functions
//either in the Base or in the DA
//here, we destroy them in the place where they belong
// pay attention to the fact that a lot of these delete are ok only if the respective function
// where the new is is called!
EqnBase::~EqnBase() {

 //========= MGOps  ===========================
    for (uint Level =0; Level<_NoLevels; Level++) {
        delete _A[Level];
        if (Level < _NoLevels - 1) delete _Rst[Level];
        if (Level > 0)             delete _Prl[Level];
    }

    _A.clear();
    _Rst.clear();
    _Prl.clear();

 //======== Vectors ==========================
    for (uint Level =0; Level<_NoLevels; Level++) {
        delete _b[Level];
        delete _res[Level];
        delete _x[Level];
        delete _x_old[Level];
        delete _x_oold[Level];
        delete _x_tmp[Level];
    }

         _b.clear();
       _res.clear();
         _x.clear();
     _x_old.clear();
    _x_oold.clear();
     _x_tmp.clear();

 //========= solver
    for (uint l=0;l<_NoLevels;l++)  delete _solver[l];
    delete []_solver;

 //========= DOF MAP ==============================
    for (uint Level = 0; Level < _NoLevels; Level++)  delete [] _node_dof[ Level];
    delete [] _node_dof;
    for (uint Level = 0; Level < _NoLevels; Level++)  { delete [] _DofNumLevFE[Level]; delete [] _DofOffLevFE[Level]; } 
    delete [] _DofNumLevFE;
    delete [] _DofOffLevFE;
    delete [] _Dim;           // dimension system Ax=b
    
 //========= refvalue and varnames =============
    delete [] _refvalue;
    delete [] _var_names;

 //=========== BOUNDARY CONDITIONS =================
 //===nodal
    delete[] _bc;                                                                   // boundary condition flag
 //===constant
     for (uint l=0; l < _NoLevels; l++)   delete [] _bc_fe_kk[l];
     delete [] _bc_fe_kk;
 //===penalty
    clearElBc();   /*if (_Dir_pen_fl==1)*/ //DO IT ALWAYS!

    
}


//=========================
//based on the Qty Internal Vector, i can define HOW MANY variables this equation has
//then compute the total number of variables
void EqnBase::initNVars()  {

    for (uint fe = 0; fe < QL; fe++)  _nvars[fe]= 0;

    for (uint i=0;i< _QtyInternalVector.size(); i++) {
        for (uint fe = 0; fe < QL; fe++)  if ( _QtyInternalVector[i]->_FEord == fe) {
                _nvars[fe] +=  _QtyInternalVector[i]->_dim;
            }
        }

    
    _n_vars = 0;
    for (uint fe = 0; fe < QL; fe++) _n_vars += _nvars[fe];

    
  for (uint fe = 0; fe < QL; fe++) {
    _VarOff[fe] = 0; 
  for (uint fe2 = 0; fe2 < fe; fe2++) _VarOff[fe] += _nvars[fe2];
  }
  
    
//CHECK LINEAR MESH WITH QUADRATIC SHAPES
    const uint mesh_ord = (int) _mesh._mesh_rtmap.get("mesh_ord");
    if (mesh_ord == 1 && _nvars[QQ]>0) {
        std::cout << "Can't handle linear mesh with quadratic variables, do a quadratic mesh" << std::endl;
        abort();
    }

    return;
}


//====================
// by default, all the reference values are initialized to 1.
//This is a function that doesnt make distinction BETWEEN various FE,
// it treats the variables in the same manner
void EqnBase::initVarNamesRefValues(std::string varname_in) {

    assert(_n_vars > 0);

    _var_names = new std::string[_n_vars];       // names
    _refvalue  = new      double[_n_vars];           // refvalues

    std::ostringstream name;
    for (uint i=0;i< _n_vars; i++) { // variable name
        name.str("");
        name << varname_in << i+1;
        _var_names[i] = name.str();
         _refvalue[i] = 1.;
    }

    return;
}



// ============================================================
/// This function initializes the system degrees of freedom (dof)
/// for every Level you have a different dof map,
/// but this is allocated as the fine level, therefore
/// except the fine level the _nodedof is "full of holes"
/// Why dont we allocate the dof map "without holes" ?
/// because we want to refer to the FINE MESH node numbering
/// which is unique.
/// In fact, off_nd_q and off_nd_l give us the RANGES in the FINE NODES
/// The dof map IS NOT THE IDENTITY !
/// node numbering is what is given by the GENCASE;
/// dof numbering instead starts from ZERO for every level!
/// i think it is the identity only on the COARSE LEVEL
/// For every level, you construct a _node_dof[Level]
/// here you start giving an index to the dofs
/// starting SUBDOMAIN BY SUBDOMAIN
/// For each subdomain, you pick FIRST the QUADRATIC then the LINEAR dofs
// void EqnBase::InitMeshToDof(const uint Level) {
//let us do the level loop inside here,
//so we can allocate and fill the nodedof
//clearly, if we put a new inside here, we cannot call this function twice,
//because we have only one delete.
//so, we should either call this function in the base constructor or something
//the point is that here the loop is over the SCALAR components already,
//but first it should be over the ordered Quantities
//Pick a quantity, pick its FE order, its number of scalar components
// _off_nd[QQ] gives the node indices corresponding to the quadratic dofs for that lev and subd
//_off_nd[LL] is different: it gives how many linear nodes you have after one quadratic node
//so, to get the node index, you have to start from the quadratic node and add the offset of the linear
// so, _off_nd[QQ] gives ABSOLUTE NODE values
// instead, _off_nd[LL] gives "RELATIVE" node values wrt the quadratic nodes:
// in order to obtain ABSOLUTE nodes you have to START FROM THE ABSOLUTE QUADRATIC and move
// according to the relative offset
// "quadratic nodes" = "nodes which correspond to quadratic dofs"
// "linear nodes" = "nodes which correspond to linear dofs"
// so those that are ABSOLUTE may be used also as RELATIVE when you do DIFFERENCES
// but those that are RELATIVE may only be used as RELATIVE
// notice that for each level, you have to SUBDIVIDE into SUBDOMAINS,
//and within each subdomain you'll have to put ux,uy,uz and p
// la versione "ASSOLUTA" e' in pratica la _node_map
// l'unica differenza e' che la node_map e' su UN LIVELLO IN PIU'
//noi praticamente dobbiamo prendere DUE LIVELLI CONTIGUI e metterli UNO IN FILA ALL'ALTRO
//infatti io guarderei anche alla node_map per riempire questa node_dof,
//perche' secondo me sono la stessa cosa!

//here we have something to consider carefully
//we must understand how to build the nodedof map
//here, we first put all the QUADRATIC VARIABLES
// and then all the LINEAR variables
// - so we see that things here are based on nvars
// all divided. Suppose that you have quadratic
//linear and quadratic, it first puts the
// two quadratic and then the linear one
//this seems to be alright

//Now we have to add the elements as dofs
//InitMeshToDof basically fills the _node_dof map
//now we are gonna do an _elem_dof thing

//From here we see that QUADRATIC DOFS are put first,
// then LINEAR DOFS,
//and then we'll add CONSTANT DOFS.

//Devo pensare se fare elem_dof oppure una unica node_dof,
//che diventera' piuttosto "_geometrical_entity_to_dof"

//count starts from zero for every level and then it grows.

//node dof e' una mappa che tiene conto di TUTTI i PROCESSORI

//Now the point is: Once you introduce ELEMENTS,
//then you're gonna have NODES to which correspond
// BOTH QUADRATIC DOFS
// AND LINEAR DOFS
// and then you have ELEMENTS to which correspond only CONSTANT DOFS.
//So, do we need to INCREASE the OFFSET for the "GEOMETRICAL ENTITIES"?!
//We'll see that.
//Now, actually the biggest number you could have of GEOMETRICAL ENTITIES
//is the NUMBER OF QUADRATIC NODES.
//But let us not think of QUADRATIC NODES as "quadratic".
//Let us think of NODES as MESH GEOMETRICAL ENTITIES.
//We have NODES, we have ELEMENTS,
//those are GEOMETRICAL ENTITIES.
//On the other hand, we have FINITE ELEMENT FAMILIES,
//each of which is characterized by a certain SET OF DEGREES OF FREEDOM,
//so, HOW DO WE MATCH THE TWO?

//===========
// ogni processore inizia ad occupare la node_dof per un pezzo ben preciso,
//un blocco contiguo che pero' e' diverso per ogni processore
//     OFF_ND
//     || ---|---|---  || ---|---|---  ||
//       <--SUBD0----->  <---SUBD1----->
//       <--->           <--->              Level0 - NODES
//       <------->       <------>           Level1 - NODES
//       <------------>  <------------>     Level2 - NODES

//     This was for the NODES, for the ELEMENTS we have sthg like
//     OFF_EL
//     || ---|---|---  || ---|---|---  ||
//       <--SUBD0----->  <---SUBD1---->
//       <--->           <--->              Level0 - ELEMENTS
//           <-->             <-->          Level1 - ELEMENTS
//                <--->           <--->     Level2 - ELEMENTS

//============
// linear -----------------------------------
//why do we put the linears in those positions?
//maybe because we want some correspondence with the mesh
//anyway, the point is that this node_dof is actually full of -1,
//at all levels, both QQ and LL
//the only level that is completely WITHOUT -1 is the QUADRATIC FINEST level.
//now, i dont understand why we don't do something like
//       for (int k1= 0;
//                k1< _mesh._off_nd[LL][off_proc + Level+1]  - _mesh._off_nd[LL][off_proc]; k1++)
//So, why don't we start from the ZERO position? Because we are supposed to put the dof number
//at the CORRECT NODE POSITION; so, if it is a quadratic node, we put it at the quadratic position,
//if it is a linear node, we put it at the linear position,
//if it is an element dof, we put it at the element position, and so on and so on.
//Now, for the elements things may be simpler, because one element belongs to ONLY ONE Level,
//as a matter of fact, you dont really have an ABSOLUTE NUMBERING of ELEMENTS.
//Well, the connectivity gives the numbering of the Elements at each level
//Since we do not have to rely on some kind of Absolute Numbering for the Elements,
//while we have to consider that for the NODES (because we have a well established list of nodes
//on which we have to print all our fields),
//then for the elements we will simply start from zero for every level
// (of course we will go ahead after the previous quadratic and linear dofs)

//So, remember that you have to put the DOF in the position of the corresponding NODE,
//so that you can use this map!
//So since the positions are always with respect to the QUADRATIC MESH,
// you have to set the quadratic nodes, and as a matter of fact, linear nodes 
// are a subset of quadratic nodes

//
// From the fact that we add this guy for the linear dofs _mesh._off_nd[QQ][off_proc]
//I think that it means that the NODES are numbered FIRST QUADRATIC, THEN LINEAR... is that true?
//Wait, some quadratic nodes are also linear already...

//The count variable will change from Level to Level.
//We have that "count" must be equal to _Dim[Level]. 
//If that is not the case, we are missing something


//Now, the point is: Every subdomain fills different pieces of the _node_dof map

//TODO I want to write the three loops for QUADRATIC, LINEAR and CONSTANT ALTOGETHER

//Now, what happens here... 
//NODES are ordered in a unique set, separated by subdomain and by Level
//Somewhat it makes sense to have nodes in a unique set because the same node belongs to DIFFERENT LEVELS.
//It happens that elements are ordered as a unique set too,
// but in this case each Element belongs to ONE and ONLY ONE LEVEL.
// So there is a unique list of element, which we will call ABSOLUTE list,
// and the elements in there have an ABSOLUTE number.
// The ABSOLUTE element number is used to retrieve the connectivity
// Now, for every Level, we want to have the list of elements 
// starting from ZERO. In this way at each level we can fill the node_dof appropriately.
// Of course, you want each subdomain to start at a DIFFERENT POINT,
//So you dont want k1 to start from zero for all processors
// Maybe the best way to do that is:
// 1) to SUM the elements of the previous domains,
// and put that as an OFFSET for starting;
// 2) also, it could be good to read the elem connectivities at SEPARATE levels.
// so, for every level we have to subtract the sum of the elements at the previous levels!

// So, so far, we decide we'll explore the ELEM CONNECTIVITY with the ABSOLUTE Element number,
// and the node_dof map with the "LEVEL-BASED" Element number.
//so we must subtract the sum of the previous levels, 
//which is a function of the current subdomain.

//TODO we must make a CLEAR DISTINCTION between DOF_OBJECTS and DOF_INDICES!!!
// The DOF_OBJECTS are the INPUTS in the node_dof map
// The DOF_INDICES are the OUTPUTS
// Every GEOMETRICAL ENTITY (NODE or ELEMENT) is characterized by MORE DOFS;
// so, in order to make a FUNCTION, we must do the association
//    GEOMETRICAL ENTITY --> DOF_OBJECTS --> DOF_INDICES
// Now the idea is this: for the NODES, everything was built in such a way that
// the DOF_OBJECTS are the SAME AT ALL LEVELS.
//So, if you have a dof_object of some node at some level L, that is going to be the SAME dof_object 
// also at the FINE LEVEL.
// So, dof_object means "INPUT POSITION in the node_dof map"
// Now, for ELEMENTS the situation is different.
// We can choose two possibilities:
// either build a _bc array FOR EVERY LEVEL for the element;
// or building a _bc array ONLY AT THE FINE LEVEL,
// and then know how to convert from the fine level to the CURRENT LEVEL.
// In this second case you would say: for every level, what is the boundary condition
// at my level given the boundary conditions on the fine level?

// ok, i need to print some FIELD at ALL LEVELS
// quindi mi ci vuole una routine che mi LINEARIZZA le CONNETTIVITA' a TUTTI I LIVELLI,
// in modo che io posso sempre vedere, siccome paraview non mi stampa quadratico in 3D.
// dopodiche' la stampa verra' fatta a ciascun livello su ciascun mesh LINEARIZZATO.
// Allora, il punto e' questo: se usi un multigriglia e fai su e giu', allora i vettori intermedi
//non hanno significato perche' contengono valori temporanei, residui, ecc...
//ma se fai la stessa risoluzione a diversi livelli, allora ti puo' interessare stampare indipendentemente 
// a diversi livelli! sia per le bc sia per x_old. Facciamolo.

//=====================================
// The DofObjects here, i.e. the "dof carriers", are NODES and ELEMENTS.
// Assume a QUADRATIC MESH.
// - For QQ FE, all Nodes are dof carriers
// - For LL FE, SOME Nodes are dof carriers (and the others are -1)
// - For KK FE, all Elements are dof carriers

// - For the NODES, we always refer to the FINE numbering (Then, you may always want to use the maps 
// to obtain the LEVEL numberings)
// - For the ELEMENTS, it is different, because any Element belongs to ONLY ONE level.
// So, what you do is you count the number of elements of the previous processors at THAT level,
// and then you go ahead up to the following processor
//======================================

// When do you have a -1 in the node_dof? For the LINEAR nodes, and also for the QUADRATIC nodes at the non-fine level.

//TODO this function is based on the assumption that we order
// FIRST QUADRATIC, THEN LINEAR, THEN CONSTANT DOFS
//I would like to see what happens if we remove this ordering
// (of course all the matrices and MG operators would have consistent ordering)
// also, it affects the way the ELEM MATRIX and RHS are FILLED

// The node dof is used by MATRICES, RECTANGULAR or NOT, and VECTORS,
//so it is basically the fundamental part of the equation

void EqnBase::ComputeMeshToDof() {
  
  
    assert(  _n_vars > 0);
    assert(_NoLevels > 0);

    _Dim      = new uint [_NoLevels];

    for (uint Level=0; Level < _NoLevels; Level++)   {
      
      int ndof_onevar[QL];
      ndof_onevar[QQ] = _mesh._NoNodesXLev[Level];
      ndof_onevar[LL] = _mesh._NoNodesXLev[_NoLevels];
      if (Level>0) ndof_onevar[LL] = _mesh._NoNodesXLev[Level-1];
      ndof_onevar[KK] =  _mesh._NoElements[VV][Level];

      _Dim[Level] = 0;
      for (uint fe=0;fe<QL;fe++) _Dim[Level] += _nvars[fe]*ndof_onevar[fe];

        std::cout << " Level "<< Level
                  << ", Mesh Quadratic Level Nodes " << _mesh._NoNodesXLev[Level]
                  << ", Quadratic DOFs " << _nvars[QQ]*ndof_onevar[QQ]
                  << ", Linear DOFs "    << _nvars[LL]*ndof_onevar[LL]
                  << ", Constant DOFs "  << _nvars[KK]*ndof_onevar[KK]  << std::endl;
		  
    }  
  
  
  _DofLocLevProcFE = new uint**[_NoLevels];     
   for (uint Level = 0; Level < _NoLevels; Level++) {
       _DofLocLevProcFE[Level] = new uint*[_mesh._NoSubdom];
        for (uint isubdom=0; isubdom<_mesh._NoSubdom; isubdom++) {
            _DofLocLevProcFE[Level][isubdom] = new uint[QL];
            _DofLocLevProcFE[Level][isubdom][QQ] = _mesh._off_nd[QQ][isubdom*_NoLevels+Level+1] - _mesh._off_nd[QQ][isubdom*_NoLevels];
            _DofLocLevProcFE[Level][isubdom][LL] = _mesh._off_nd[LL][isubdom*_NoLevels+Level+1] - _mesh._off_nd[LL][isubdom*_NoLevels];
            _DofLocLevProcFE[Level][isubdom][KK] = _mesh._off_el[VV][isubdom*_NoLevels+Level+1] - _mesh._off_el[VV][isubdom*_NoLevels+Level];
       }
    }
  
  //SETUP DOF OFFSETS for NODE DOF MAP ********************
  //I need to set up the dof offsets for every level
  //Since for nodes they do not depend on fe, but for elements they do, in general they depend on fe
  _DofNumLevFE = new uint*[_NoLevels];
  for (uint Level = 0; Level < _NoLevels; Level++) {
  _DofNumLevFE[Level] = new uint[QL];
  _DofNumLevFE[Level][QQ] = _mesh._NoNodesXLev[_NoLevels-1];
  _DofNumLevFE[Level][LL] = _mesh._NoNodesXLev[_NoLevels-1];
  _DofNumLevFE[Level][KK] = _mesh._NoElements[VV][Level];
   }
      
  _DofOffLevFE = new uint*[_NoLevels];
  for (uint Level = 0; Level < _NoLevels; Level++) {
  _DofOffLevFE[Level] = new uint[QL];
  
 for (int fe=0; fe<QL; fe++) { 
  _DofOffLevFE[Level][fe] = 0;
  for (int fe2=0; fe2 < fe; fe2++) _DofOffLevFE[Level][fe] += _nvars[fe2]*_DofNumLevFE[Level][fe2];
    }
  }

  //fill node dof ******************************************
    _node_dof = new int*[_NoLevels];

    for (uint Level = 0; Level < _NoLevels; Level++) {
      
       int nodedof_size = 0;
            for (uint fe=0;fe<QL;fe++)      nodedof_size += _nvars[fe]*_DofNumLevFE[Level][fe];
   
        _node_dof[Level] = new int[nodedof_size];
            for (int k1=0;k1< nodedof_size;k1++) _node_dof[Level][k1] = -1;

        uint dof_count_lev=0;
        for (uint isubdom=0; isubdom<_mesh._NoSubdom; isubdom++) {
            uint off_proc=isubdom*_NoLevels;

            // quadratic -----------------------------------
            for (uint ivar=0;ivar<_nvars[QQ];ivar++) {
                for (int fine_node = _mesh._off_nd[QQ][off_proc];
                         fine_node < _mesh._off_nd[QQ][off_proc+Level+1];fine_node++) {
 #ifndef NDEBUG
 		  if ( fine_node >= (int) _DofNumLevFE[Level][QQ] ) {  std::cout << " Wrong QQ " << std::endl; abort(); }
 #endif	    
                    _node_dof[Level][ fine_node + ivar*_DofNumLevFE[Level][QQ] + _DofOffLevFE[Level][QQ] ] = dof_count_lev;
                    dof_count_lev++;
                }
            }

            // linear -----------------------------------
            for (uint ivar=0;ivar<_nvars[LL];ivar++) {
                for (int fine_node = _mesh._off_nd[QQ][off_proc];
                         fine_node < _mesh._off_nd[QQ][off_proc]
                        + _mesh._off_nd[LL][off_proc + Level+1]
                        - _mesh._off_nd[LL][off_proc]; fine_node++) {
#ifndef NDEBUG
  		  if ( fine_node >= (int) _DofNumLevFE[Level][LL] ) {  std::cout << " Wrong LL " << std::endl; abort(); }
#endif
                   _node_dof[Level][ fine_node + ivar*_DofNumLevFE[Level][LL] + _DofOffLevFE[Level][LL] ] = dof_count_lev;
                    dof_count_lev++;
                }
            }
           
            // constant -----------------------------------
            int sum_elems_prev_sd_at_lev = 0;
	    for (uint pr = 0; pr < isubdom; pr++) { sum_elems_prev_sd_at_lev += _mesh._off_el[VV][pr*_NoLevels + Level + 1] - _mesh._off_el[VV][pr*_NoLevels + Level]; }
           
            for (uint ivar = 0; ivar<_nvars[KK]; ivar++) {
                for (uint elem =  sum_elems_prev_sd_at_lev;
		                 elem <  sum_elems_prev_sd_at_lev
		                     + _mesh._off_el[VV][off_proc + Level+1]
                                     - _mesh._off_el[VV][off_proc + Level]; elem++) {
#ifndef NDEBUG
		  if ( elem >= _DofNumLevFE[Level][KK] ) {  std::cout << " Wrong KK " << std::endl; abort(); }
#endif
                    _node_dof[Level][elem + ivar*_DofNumLevFE[Level][KK] + _DofOffLevFE[Level][KK] ] = dof_count_lev;
                    dof_count_lev++;
                }
            }

        }  //end subdomain
        
#ifndef NDEBUG
       if ( dof_count_lev  != _Dim[Level] ) { std::cout << "There is a mismatch between count and _Dim[Level] " << dof_count_lev << " " << _Dim[Level] << std::endl; abort(); }
#endif

        
#ifdef DEFAULT_PRINT_INFO
        std::cout << "EqnBase::InitMeshToDof(D)   Level= " << Level  << std::endl;
#endif

    } //end Level
    
    PrintMeshToDof();

    return;
}



// For every Level, we loop over the DOF FE FAMILIES.
// For every DOF FE FAMILY, we associate the GEOMETRICAL ENTITIES on top of which that FE FAMILY is built.
// So, QUAD9 FE is associated to QUADRATIC QUAD9 NODES
// You could also associate QUAD8 FE to QUADRATIC QUAD9 NODES,

//On the other hand, we could think first of the GEOMETRICAL ENTITIES: "QUADRATIC NODES", "LINEAR NODES", "ELEMENTS",
//And loop over the GEOMETRICAL ENTITIES, and ask: WHAT ARE THE DOFS ASSOCIATED TO THESE ENTITIES?
//So, I would better just consider "QUADRATIC NODES" and "ELEMENTS", or "QUAD9 NODES" and "ELEMENTS";
// or, if you have a "QUAD8 MESH", you would have "QUAD8 NODES" and "ELEMENTS".
//Now, what do you do if you have a QUAD9 FE Family? No problem, in your computations you use another DOF,
//which is not directly associated to any Geometrical Entity (in the "computer science" sense that the 9th node does not exists in the file..).
//But, eventually, when you PRINT the QUAD9 ONTO the QUAD8 MESH, you somehow should take into account the presence of the 9TH DOF,
// WHICH SHOULD CONTRIBUTE TO the DOF COEFFICIENTS in the PROJECTED SPACE.

//So, putting the GEOMETRICAL ENTITIES FIRST, "NODES" and "ELEMENTS",
// you may have:
//NODES with QUADRATIC AND LINEAR FE DOFS
//NODES with QUADRATIC FE DOFS
//ELEMENTS with CONSTANT DOFS

//What if you have "GEOMETRICAL ENTITIES" WITH "NO ASSOCIATED DOFS"? They just do not contribute.
// Or "DOFS" with "NO GEOMETRICAL ENTITY", viceversa? You just PROJECT their contribution to the OTHER EXISTING GEOMETRICAL ENTITIES that HAVE DOFS of THAT SAME FE FAMILY.
//So, another example, if you want to PRINT a QUAD4 FE SPACE onto a QUAD9 NODE GEOMETRY, 
//you would be going from a SMALLER to a LARGER space,
// so you would do some PROLONGATION.
//If you went from a LARGER to a SMALLER space,
// you would do some sort of PROJECTION ONTO SUBSPACE.
// Somehow you may see the MG operators in the same way:

// Restriction  = FROM LARGER TO SMALLER SPACE (PROJECTION ONTO SUBSPACE)
// Prolongation = FROM SMALLER TO LARGER SPACE

// Now, let us go back, and let us consider being in the case of looping over the DOF FE FAMILIES, 
// and for each FE Family we pick the CORRESPONDING "DOF GENERATING" GEOMETRICAL ENTITIES.
// The point is: the existence of the DOF should actually INTRINSIC, in the sense that it does not depend
// on the UNDERLYING "DOF GENERATING" GEOMETRICAL ENTITIES.
// You might wanna have QUAD9 FE on a linear mesh, for instance, and do all the computations as quadratic
// and only at the end you project to the LINEAR MESH.

// "GEOMETRICAL ENTITY" (NODES or ELEMENTS) may be also called "DOF GENERATING GEOMETRICAL OBJECT" or "DOF OBJECT"

//The only thing that we are not printing here is the SUBDOMAIN SUBDIVISION

//You should print this map while you are building it, and show precisely
//where you SEPARATE things between one SUBDOMAIN and the other

void EqnBase::PrintMeshToDof() const {

  
      for (uint Level = 0; Level< _NoLevels; Level++) { 
	
	std::cout << "========== Level " << Level << " =========="  << std::endl;

	for (int fe=0; fe<QL; fe++) {
	  
	  	std::cout << "========== FE " << fe << " =========="  << std::endl;
	  
	for (uint ivar=0; ivar<_nvars[fe]; ivar++) {
	  
  	  	std::cout << "========== VAR " << ivar << " =========="  << std::endl;
		
	     for (uint i=0; i < _DofNumLevFE[Level][fe]; i++) {
  
	  std::cout <<  _node_dof[Level][ i + ivar*_DofNumLevFE[Level][fe] + _DofOffLevFE[Level][fe] ]   << std::endl;
	  
	  //inside here I can also find a way to print at the end of each subdomain, depending on "i"
	  
	     } //end i
	     
	  }
	}
	
      }
  
  
 return; 
}




//==================
//this function DEPENDS in _iproc!!!
void EqnBase::initVectors() {

    //allocation
         _x.resize(_NoLevels);
     _x_old.resize(_NoLevels);
    _x_oold.resize(_NoLevels);
     _x_tmp.resize(_NoLevels);
         _b.resize(_NoLevels);
       _res.resize(_NoLevels);

    for (uint Level = 0; Level< _NoLevels; Level++) {

    uint ml[QL];    for (int fe=0; fe<QL; fe++) ml[fe] = _DofLocLevProcFE[Level][_iproc][fe];
    uint m_l = 0;
    for (int fe=0; fe<QL; fe++)  m_l +=  ml[fe]*_nvars[fe];
    
        _b[Level] = NumericVector::build().release();
        _b[Level]->init(_Dim[Level],m_l,false,AUTOMATIC);
        _res[Level] = NumericVector::build().release();
        _res[Level]->init(_Dim[Level],m_l,false,AUTOMATIC);
        _x[Level] = NumericVector::build().release();
        _x[Level]->init(_Dim[Level],m_l,false,AUTOMATIC);
        _x_old[Level] = NumericVector::build().release();
        _x_old[Level]->init(_Dim[Level],false, SERIAL);
        _x_oold[Level] = NumericVector::build().release();
        _x_oold[Level]->init(_Dim[Level],false, SERIAL);
        _x_tmp[Level] = NumericVector::build().release();
        _x_tmp[Level]->init(_Dim[Level],false, SERIAL);

    } //end level loop

    return;
}

 //TODO every processor does ALL because bc is SERIAL
 //Now here we have to think how to impose the boundary conditions for KK
 //Now, the bc_read function orders the boundary condition flags of a QUADRATIC or LINEAR DOF OBJECT.
 //TODO The point is that you already have to know the ORDER in which the UNKNOWNS are settled in your system,
 // FIRST QUADRATIC, THEN LINEAR, THEN CONSTANT.
 // as a matter of fact, the bc_field should be ordered like a DOF vector,
 //so have first quadratic, then linear, then constant variables.
 //The fact is that by now we only have ONE bc FIELD, ordered like dofs yes, but ONLY FOR THE *FINE* LEVEL.
 // This is because every coarser NODE is always also FINE, and so that is ok.
 //But, for constant elements I want to make a LEVEL DEPENDENT bc_map,
 //  whose length changes with Level, just like the node_dof map.
 //Then you will have to loop over the BOUNDARY ELEMENT at EACH LEVEL,
 // and assign a FLAG to the VOLUME DOF ELEMENTS.
 // climbing back from Boundary to Volume (we know how to do that)
 // In order to pick the boundary element, we pick its MIDDLE POINT.
 
 //TODO why don't we better impose the boundary conditions
 // based on QUANTITIES and then SCALAR COMPONENTS of QUANTITIES?
 // Then depending on the type of dof you would have to 
 // loop over the DOF OBJECTS generating those dofs!
 // So you would loop either over NODES or over ELEMENTS and so on.
 
 //Wait, we did two functions for QQ-LL and KK, but actually one function seems to be enough!!!
 //The only thing is that you need to distinguish between passing the NODE COORDINATE or the MIDDLE POINT of an ELEMENT,
 //but we can MODIFY the bc_read INTERFACE and MAKE A UNIQUE FUNCTION that ENCOMPASSES EVERYTHING.
 // Yes we can do that.
 
 //Now, can we make a UNIQUE Function for filling the _bc fields, both for Nodes and for Elements?
 //One of the two is defined only on the FINE level, the other on ALL LEVELS
 
 //BC is a property of the DOFS. Now, since to every DOF at any level has a corresponding FINE DOF,
 // and for consistency the bc_flag must be THE SAME at ALL LEVELS, then for the nodes
 // we have a UNIQUE FIELD for all levels.
 // Now, of course the DOF NUMBER CHANGES from one level to another. 
 // But, the POSITION of the DofObject (= "source" of the dof number) in the _node_dof[Level]
 // is THE SAME!
 //That is why we have a lot of (-1) in the _node_dof map.
 
 //Ok, now i did a UNIQUE FUNCTION GenBc()
 // The point here is that for every equation I have to set TWO FUNCTIONS, bc_read and BCElemKKRead.
 //One is based on the NODE COORDINATE, one on the MIDDLE POINT of the BOUNDARY ELEMENT
 // So, for every equation i should implement two separate functions bc_read and BCElemKKRead,
 // which i don't really like
 // The problem is that bc_read is called only on the fine level... no problem, I can call it also at the other levels.
 // but, in one case i will only use xp, in the other case i will use only el_xm, so i should pass fake variables in one case or the other
 //As far as i use two SEPARATE FIELDS for the boundary conditions, _bc AND _bc_fe_kk,
 // it is obvious that I have to use two separate functions
 //Then the problem is when you try to switch one Quantity from a NODE BASED FE Family to an ELEM BASED FE FAMILY
 //You have to change the way you  ENFORCE the BOUNDARY CONDITIONS.
 //TODO the problem is that here we should write things better.
 // We should write things in such a way that if we set an equation WITHOUT VARIABLES the code still runs. 
 // Here this does not happen.
 // if you allocate bc_flag with zero components, then you pass the pointer to the bc_read function which will set the components.
 // So you have to call the bc_read function only if the number of variables is not zero, and in the same way you have to call the function
 // BCElemKKRead only if the number of ElemBased variables is not zero!
 
 
 //SWITCHING one VARIABLE from KK to LL.
 // since the length of bc_value and bc_value_kk change, in the routine i cannot leave both selections fixed!
 
 //Now the point is: I am sorry but we have to modify the structure of the boundary conditions...
 //I think we have to enforce the boundary conditions SEPARATELY for EACH QUANTITY
 //If the quantity is NODE-BASED, then we will use bc_read;
 //If the quantity is ELEMENT BASED, we will use bcelemKKread;
 //Now the point is only this: we want to do this in a UNIQUE routine,
 //so that we write the if's of the geometry only once. 
 //I guess first of all I should base everything on the BOUNDARY MIDDLE POINT.
 
//As a matter of fact I am taking the nodes but these are taken in an element loop,
// so actually i am not looping straight away on the nodes, so i have superimpositions anyway.
//so YES, i'll base anything on the ELEMENT MIDDLE POINT!
//Ok, I did that. Now the list of the bc for each scalar variable is only in one place, under bc_read
 
//Now that routine is done in such a way that you have to order FIRST QQ, THEN LL, THEN KK;
//In this way you will read correctly

//TODO one thing that may be optimized here is that there are computations that could be avoided when you DO NOT HAVE QQ variables, or LL variables, or KK variables.

// The good thing here would be to have a unique bc array at EACH LEVEL. 
// The distinction between level is especially good for the elements, because the elements 
//only belong to ONE level, it's not like the nodes...

//TODO now i am looping over elements first, and nodes inside each element next.
//therefore, if an element puts a zero and a following element puts a one, the one wins.
//so we must find a way which is independent of the element order
//if you loop over the nodes you involve them only once...
//well, there are two reasons for which you superimpose elements:
// 1 - because you loop over the boundary elements
// 2 - because you loop over all faces which have "intersections"

// I think we should do two routines for the flags: one that loops 
// over the BOUNDARY ELEMENTS,
// and one that loops over the VOLUME ELEMENTS.

//So far, let us concentrate on the boundary loop. What is the way to check 
//not to superimpose the ones on the zeros?
//just put a check before checking

//TODO ricorda che i bc sono praticamente dei dof fields,
// quindi in principio possono essere fatti in PARALLELO!
// anche x_old e' un vettore di dof. Alla fine non sono cose 
//molto diverse in fondo!

//Quindi bisogna pensare ad una cosa piu' generale.
//Faro' il bc come un Numeric Vector di INTERI,
//cosi' come x_old e' un NumericVector di DOUBLE.

//in questo modo il loop per riempire le bc 
// o il loop per riempire x_old
// non saranno molto diversi, anzi molto simili!

//per come e' fatto il mio sistema,
//direi che la cosa giusta e' fare le BC a TUTTI I LIVELLI.
// in questo modo e' tutto molto piu' chiaro.
//la cosa che devo capire un attimo e' come avere 
// il valore fisico A TUTTI I LIVELLI.
//Per i nodi e' immediato, perche' un nodo coarser e' anche 
// un nodo fine.
// Per un elemento coarser, bisogna fare cosi': ottieni la soluzione 
// FINE per gli elementi. Poi la ridistribuisci sugli elementi coarse,
// a tutti i livelli.

//Se volessimo fare le cose ottimizzate, 
//dovremmo trattare cosi' lo STORAGE della SOLUZIONE:
//Per i dof il cui dof object e' un NODO, 
// mettiamo in un vettore SOLO FINE.

//per i dof il cui dof object e' un ELEMENTO,
//mettiamo in un vettore DI TUTTI I LIVELLI.

//in questo modo minimizzeremmo lo spazio di memoria, che mi sembra 
//una cosa giusta, anche se si devono splittare un attimo queste cose.

//il punto e' che se magari facciamo prima QQ, poi KK, poi LL,
// avremmo dei pezzi separati. 
//Ma non c'e' problema, perche' noi proprio all'atto della stampa
//facciamo la separazione!
//ora, dobbiamo farla anche durante l'ALGORITMO DI SOLUZIONE.
//Pertanto, definiro' nell'equazione 
// un VETTORE di STORAGE NODALE per le VARIABILI NODALI
// e un VETTORE di STORAGE di ELEMENTO per le VARIABILI DI ELEMENTO.

//Questo vettore servira' per RETRIEVE THE DOFS FOR THE TRUE VALUES FOR ANY EQUATION,
// and to PRINT ALL THE VARIABLES, node or cell.



void EqnBase::GenBc() {
  
    std::string     basepath = _files.get_basepath();
    std::string    input_dir = _files.get_frtmap().get("INPUT_DIR");
    std::string          ibc = DEFAULT_IBC;
    std::string     ext_xdmf = DEFAULT_EXT_XDMF;
    std::string       ext_h5 = DEFAULT_EXT_H5;
    std::string  bdry_suffix = DEFAULT_BDRY_SUFFIX;
    
    std::ostringstream ibc_fileh5;
    ibc_fileh5  << basepath << "/" << input_dir << "/" << ibc << ext_h5;
// TODO actually, we should first COPY this file in the outtime dir, then READ it from there!
    std::ifstream in(ibc_fileh5.str().c_str());

    const uint Lev_pick_bc_NODE_dof = _NoLevels-1;  //we use the FINE Level as reference
    
 //************************************************   
 //******** NODE BASED ****************************   
    const uint mesh_ord    = (int) _mesh._mesh_rtmap.get("mesh_ord");
    const uint offset      = _mesh._NoNodesXLev[_NoLevels-1];
    const uint el_nnodes_b = _mesh._GeomEl._elnds[BB][mesh_ord];
    double* normal = new double[_mesh._dim];  //TODO remove this, it is useless

    int  *bc_flag = new int[_n_vars];

    _bc = new int[_Dim[Lev_pick_bc_NODE_dof]];
    for (uint i1=0;i1< _Dim[Lev_pick_bc_NODE_dof];i1++) _bc[i1] = DEFAULT_BC_FLAG;    // set 1 all the points for  bc (boundary condition)

    //both here and in bc_read you have to put the value bc=0
    //so that you get the IDENTITY OPERATOR and you only have to provide the
    //function for the RHS
 //******** NODE BASED ****************************   
 //***************************************************   

 //**************************************************   
 //******** ELEM BASED ******************************  
     CurrElem       currelem(*this,_eqnmap);

    _bc_fe_kk             =  new int*[_NoLevels];
    int* DofOff_Lev_kk    =  new int[_NoLevels];
    
    for (uint Level=0; Level <_NoLevels; Level++)   { //loop over the levels

          DofOff_Lev_kk[Level] = _nvars[KK]*_DofNumLevFE[Level][KK];
              _bc_fe_kk[Level] = new int[DofOff_Lev_kk[Level]];

        for (int i=0; i < DofOff_Lev_kk[Level]; i++)    _bc_fe_kk[Level][i] = DEFAULT_BC_FLAG;
    }
 //******** ELEM BASED ******************************   
 //************************************************   

 // Both the NODE BASED and the ELEM BASED are initialized to 1 over THE WHOLE VOLUME,
// so over ALL THE DOFS.
// Then you call the BOUNDARY FUNCTIONS, and these act only on boundary nodes or boundary elements,
// but then of course their contribution goes into the volume,
// because every Boundary NODE belongs to the VOLUME
// and every Boundary ELEMENT communicates its value ("as a TRACE") to the corresponding VOLUME ELEM DOF.

    if (!in) {   // -------- reading bc from function
 
    for (uint Level=0; Level <_NoLevels;Level++)   { //loop over the levels
	
        for (uint isubd=0; isubd<_mesh._NoSubdom; ++isubd) {
            uint iel0 = _mesh._off_el[BB][ _NoLevels*isubd + Level];
            uint ielf = _mesh._off_el[BB][ _NoLevels*isubd + Level+1];
            for (uint iel=0; iel < (ielf-iel0); iel++) {

	        currelem.get_el_nod_conn_lev_subd(BB,Level,isubd,iel);
                currelem.get_el_ctr(BB);
		
 	    for (uint ivar=0; ivar< _n_vars; ivar++)  bc_flag[ivar] = DEFAULT_BC_FLAG; //this is necessary here to re-clean!

                 if (_n_vars >0) bc_read(currelem._el_xm[BB],normal,bc_flag);

  //******************* ONLY FINE LEVEL, NODE VARS ***************** 
   if (Level == Lev_pick_bc_NODE_dof)  { 
                for (uint i=0; i<  el_nnodes_b; i++)  {
                        const uint fine_node = _mesh._el_map[BB][(iel+iel0)*el_nnodes_b+i];

                    //Set the quadratic fields
                    if (i<_AbstractFE[QQ]->_ndof[BB])
		      for (uint ivar=0; ivar<_nvars[QQ]; ivar++) {
                            int kdof = _node_dof[Lev_pick_bc_NODE_dof][ fine_node + ivar*_DofNumLevFE[Lev_pick_bc_NODE_dof][QQ] + _DofOffLevFE[Lev_pick_bc_NODE_dof][QQ] ];
                           if (_bc[kdof] != 0) _bc[kdof] = bc_flag[ ivar + _VarOff[QQ]];
                        }
                    // Set the linear fields
                    if (i<_AbstractFE[LL]->_ndof[BB]) {
                        for (uint ivar = 0; ivar < _nvars[LL]; ivar++) {
                            int kdof = _node_dof[Lev_pick_bc_NODE_dof][ fine_node + ivar*_DofNumLevFE[Lev_pick_bc_NODE_dof][LL] + _DofOffLevFE[Lev_pick_bc_NODE_dof][LL] ];
                           if (_bc[kdof] != 0) _bc[kdof] = bc_flag[ ivar + _VarOff[LL]];
                        }
                    }
                } // i
                
           }
   //******************* END ONLY FINE LEVEL ***************** 

 //******************* ALL LEVELS, ELEM VARS ***************** 
		 int sum_elems_prev_sd_at_lev = 0;
                 for (uint pr = 0; pr < isubd; pr++) { sum_elems_prev_sd_at_lev += _mesh._off_el[BB][_NoLevels*pr + Level + 1] - _mesh._off_el[BB][ _NoLevels*pr + Level]; }
		 
		 int bdry_iel_lev =  iel + sum_elems_prev_sd_at_lev;
		 int vol_iel =  _mesh._el_bdry_to_vol[Level][bdry_iel_lev];
		 
 	        for (uint ivar= 0; ivar < _nvars[KK];  ivar++)  { 
		  int dof_kk_pos_lev =  vol_iel + ivar*_DofNumLevFE[Level][KK];
	            if (_bc_fe_kk[Level][ dof_kk_pos_lev ] != 0)  _bc_fe_kk[Level][ dof_kk_pos_lev ] = bc_flag[ ivar + _VarOff[KK] ];
		}
//******************* END ALL LEVELS, ELEM VARS **************   
                 
            } // end of element loop
        }  // end subdomain
        
     }//end nolevels
       
        std::cout << "\n GenBc(DA): boundary conditions defined by  bc_read" << "\n \n";
    }
    
    // *********  reading from file ***************************
    else {     // boundary from the top grid
        std::cout << "\n GenBc: READING FROM FILE STILL WITHOUT KK ELEMENTS, CHECK new and delete in here! ************** " << ibc_fileh5.str() << "\n \n";
        abort();
	
        // temporary vector
        uint  n_nodes_q=_mesh._NoNodesXLev[_NoLevels-1];
        int *data=new int[ n_nodes_q ];

        // open file
        std::ostringstream tabname;
        hid_t  file = H5Fopen(ibc_fileh5.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

        // variable loop
        for (uint ivar=0;ivar<_nvars[0];ivar++) { // ----------------------------
            // nodes for boundary elements
            int el_nodes=el_nnodes_b; /*  if (ivar >= _nvars[0]) el_nodes=NDOF_PB;*/ /*TODO check here*/
            // read
            IO::read_Ihdf5(file,"/"+_var_names[ivar] + bdry_suffix,data);// from hdf5 file

            // storage boundary conditions <- data <- hdf5 file
            for (uint isubd=0;isubd<_mesh._NoSubdom;++isubd) {
                int iel0=_mesh._off_el[BB][_NoLevels-1+_NoLevels*isubd];
                int ielf=_mesh._off_el[BB][_NoLevels-1+_NoLevels*isubd+1];
                for (int iel=0;iel <ielf-iel0; iel++) {

                    // Element node loop
                    for (int i=0; i< el_nodes; i++)     {// element loop
                        const uint k=_mesh._el_map[BB][(iel+iel0)*el_nnodes_b+i]; // global node
                        int kdof= _node_dof[_NoLevels-1][k+ivar* offset];
                        _bc[kdof]=data[k];
                    }
                }// iel
            }//proc
        }// ivar ------------------------------------------------------------
        // clean
        delete []data;
        H5Fclose(file);
        std::cout << "\n GenBc(DA): boundary conditions defined from file " << ibc_fileh5.str() << "\n \n";
    }
    
    delete []bc_flag;
    delete []DofOff_Lev_kk;
    
    delete [] normal;
    
    return;
}



/// boundary conditions  from function
void EqnBase::bc_read(double /*xp*/[],double /*normal*/[],int bc[]) {
    for (uint ivar=0;ivar<_n_vars;ivar++) bc[ivar]=1;
    return;
}


// =====================================================================
/// This function  defines the boundary conditions for  DA systems:
// how does this function behave when e.g. both velocity and pressure are linear?
//this function runs over ALL the domain, so it's not only on a portion.
//ok we must think of this function in terms of the quantities
//we must also keep in mind that the DofMap is built by following a GIVEN ORDER;
// so, if you have 2 quadratic and 1 linear, you choose to put first the quadratic then the linear
// here, we loop over BOUNDARY ELEMENTS
//for every boundary element, loop over its nodes
// pick the coordinate of that node
//for every node, pick all its degrees of freedom from the dofmap
//at THAT COORDINATE, the DoF loop corresponds to a Quantity loop



// ================================================================
/// This function prints the boundary conditions: for quad and linear fem
//ok, this routine prints FLAGS, not VALUES.
//therefore, I'm not interested in INTERPOLATING FLAGS.
//Flags are either one or zero, no intermediate values are requested
//if all the nodes are 0, then it means do the pressure integral ---> put a zero
//if at least one node is 1 ---> put a 1. (anyway, the bc_p are not used in "volume" manner)
//the problem is always that we are printing LINEAR VARIABLES on a QUADRATIC mesh,
// so actually the interpolated values are there just for printing...
//

//TODO I must do in such a way as to print bc just like x_old,
// it is supposed to be exactly the same routine

//TODO I should print these bc flags, as well as the solution vectors, for ALL LEVELS!
// in this way I could see what happens at every level
// I already have the connectivities at all levels, thanks to gencase
// The only problem is that in 3D I do not see the connectivities, so I should convert all the levels to LINEARIZED REFINED

//also, apart from very little things, the routine is very similar to printing any field defined on each level
// And it's an integer and when you do the average you need to use the ceil() function...

//ok, when you print a field on a coarser grid, you need to specify a coarser topology,
// but also a coarser number of points...
// but, when I print only BOUNDARY MESH or VOLUME MESH at all levels,
// I put the FINE numbers and it's ok...
// so how does it pick the coordinates?


// A = NUMBER SEEN in PARAVIEW
// B = NUMBER in HDF5

//Ok, when i load coarser meshes, the number of cells changes but not the number of points...
// but the thing is: it seems that the NODES are ORDERED in SUCH A WAY THAT 
// you have first COARSER then 

//AAA: the equation A=B only holds with the WHOLE MESH!!!
// When you add some FILTER in paraview, for instance if you do a CLIP,
// Then the numbers DO NOT CORRESPOND ANYMORE!!!

//Ok, if i have to print a data attribute on a grid, i need to provide a SMALLER NUMBER of COORDINATES,
// or maybe make a LARGER DATA VECTOR...
//Due to the way we order the nodes, I guess we can say in the XDMF that the array CAN BE CUT! Let's see...
//You cannot put in the XDMF a SMALLER NUMBER than the dimension of the HDF5 ARRAY,
//otherwise it gives SEGMENTATION FAULT!
//So I need to have HDF5 fields for the coordinates of EACH LEVEL!
//I dont wanna put useless values, so I'll just print the coordinates for each level...

//Now there is a problem with the TIME printing... I am including all the solution files but it is picking up
// the COARSE LEVEL instead of the FINEST ONE

//TODO the problem is that now every solution file contains DIFFERENT GRIDS,
// and when I load the TIME file then it loads THE FIRST ONE IT MEETS in the SOL FILES!!!

// AAA, ecco l'errore! la connettivita' che percorriamo per interpolare i valori lineari 
// e' quella FINE, ma noi ora dobbiamo prendere quella DI CIASCUN LIVELLO SEPARATAMENTE!


void EqnBase::PrintBc(std::string namefile) {

    hid_t file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

    std::string  bdry_suffix = DEFAULT_BDRY_SUFFIX;

    const uint Lev_pick_bc_NODE_dof = _NoLevels-1;  //we use the FINE Level as reference
    
    // ==========================================
    // =========== FOR ALL LEVELS ===============
    // ==========================================
    for (uint Level = 0; Level < _NoLevels; Level++)  {
    
      std::ostringstream grname; grname << "LEVEL" << Level;
  
   int NGeomObjOnWhichToPrint[QL];
    NGeomObjOnWhichToPrint[QQ] = _mesh._NoNodesXLev[Level];
    NGeomObjOnWhichToPrint[LL] = _mesh._NoNodesXLev[Level];
    NGeomObjOnWhichToPrint[KK] = _mesh._NoElements[VV][Level]*_mesh._GeomEl.n_se[VV];
  
    const uint n_nodes_lev = _mesh._NoNodesXLev[Level];
    int* sol_on_Qnodes = new int[n_nodes_lev];  //this vector will contain the values of ONE variable on ALL the QUADRATIC nodes
    //for the quadratic variables it'll be just a copy, for the linear also interpolation
    
    // ===================================
    // ========= QUADRATIC ================
    // ===================================
    for (uint ivar=0;ivar<_nvars[QQ]; ivar++)        {
      
      for (uint isubdom=0; isubdom<_mesh._NoSubdom; isubdom++) {
            uint off_proc=isubdom*_NoLevels;
     
          for (int fine_node = _mesh._off_nd[QQ][off_proc];
	           fine_node < _mesh._off_nd[QQ][off_proc+Level+1]; fine_node ++) {

	int pos_in_sol_vec_lev = _node_dof[Lev_pick_bc_NODE_dof][ fine_node + ivar*_DofNumLevFE[ Lev_pick_bc_NODE_dof ][QQ] + _DofOffLevFE[ Lev_pick_bc_NODE_dof ][QQ] ];
	int pos_on_Qnodes_lev  = _mesh._Qnode_fine_Qnode_lev[Level][ fine_node ]; 

	sol_on_Qnodes[ pos_on_Qnodes_lev ] = _bc[pos_in_sol_vec_lev];
          }
       }  //end subd

       std::ostringstream var_name; var_name << _var_names[ ivar + _VarOff[QQ] ] << "_" << grname.str() << bdry_suffix;
       hsize_t dimsf[2];  dimsf[0] = NGeomObjOnWhichToPrint[QQ];  dimsf[1] = 1;
       IO::print_Ihdf5(file_id,var_name.str(),dimsf,sol_on_Qnodes);
    }

    // ===================================
    // ========= LINEAR ==================
    // ===================================
    uint elnds[QL_NODES];
    elnds[QQ] =_mesh._GeomEl._elnds[VV][QQ];
    elnds[LL] =_mesh._GeomEl._elnds[VV][LL];
    double *elsol_c=new double[elnds[LL]];

    for (uint ivar=0; ivar<_nvars[LL]; ivar++)   {

	for (uint isubdom=0; isubdom<_mesh._NoSubdom; isubdom++) {
	              uint off_proc=isubdom*_NoLevels;

            for (int fine_node =   _mesh._off_nd[QQ][off_proc];
                     fine_node <   _mesh._off_nd[QQ][off_proc]
                    + _mesh._off_nd[LL][off_proc+Level+1]
                    - _mesh._off_nd[LL][off_proc]; fine_node++) {
	      
            int pos_in_sol_vec_lev = _node_dof[Lev_pick_bc_NODE_dof][ fine_node  + ivar*_DofNumLevFE[ Lev_pick_bc_NODE_dof ][LL] + _DofOffLevFE[ Lev_pick_bc_NODE_dof ][LL] ]; 
 	    int pos_on_Qnodes_lev  = _mesh._Qnode_fine_Qnode_lev[Level][ fine_node ];

	    sol_on_Qnodes[ pos_on_Qnodes_lev ]=  _bc[pos_in_sol_vec_lev]; 
          }
        }
 
        //  2bB element interpolation over the fine mesh 
        for (uint iproc=0; iproc<_mesh._NoSubdom; iproc++) {
               uint off_proc = iproc*_NoLevels;
	       int iel_b = _mesh._off_el[VV][off_proc + Level];
	       int iel_e = _mesh._off_el[VV][off_proc + Level + 1];
            for (int iel = 0; iel < (iel_e-iel_b); iel++) {

                for (uint in=0; in < elnds[LL]; in++)  {
 		  int pos_Qnode_fine = _mesh._el_map[VV][ (iel+iel_b)*elnds[QQ]+in ];
		  int pos_Qnode_lev  = _mesh._Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];
        	  elsol_c[in]= sol_on_Qnodes[ pos_Qnode_lev ];
		}
		
		  for (uint in=0; in < elnds[QQ]; in++) { // mid-points
                    double sum=0;
                    for (uint jn=0; jn<elnds[LL]; jn++) {
                        sum += _AbstractFE[LL]->get_prol(in*elnds[LL]+jn)*elsol_c[jn];
                    }
                    
                    int pos_Qnode_fine = _mesh._el_map[VV][ (iel+iel_b)*elnds[QQ]+in ];
                    int pos_Qnode_lev  = _mesh._Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];
              
                    sol_on_Qnodes[ pos_Qnode_lev ] = ceil(sum);   //the ceiling, because you're putting double over int!
                }
            }
        } // 2bB end interpolation over the fine mesh

       std::ostringstream var_name; var_name << _var_names[ ivar + _VarOff[LL] ] << "_" << grname.str() << bdry_suffix;
       hsize_t  dimsf[2];  dimsf[0] = NGeomObjOnWhichToPrint[LL];  dimsf[1] = 1;
       IO::print_Ihdf5(file_id,var_name.str(),dimsf,sol_on_Qnodes);
    } // ivar
    
    delete [] elsol_c;
    delete [] sol_on_Qnodes;

     // ===================================
     // ========= CONSTANT ================
     // ===================================
     for (uint ivar=0; ivar<_nvars[KK]; ivar++)        {
      
  int *sol_on_cells;   sol_on_cells = new int[ NGeomObjOnWhichToPrint[KK] ];
  
  int cel=0;
  for (uint iproc=0; iproc<_mesh._NoSubdom; iproc++) {
               uint off_proc = iproc*_NoLevels;
	       
            int sum_elems_prev_sd_at_lev = 0;
	    for (uint pr = 0; pr < iproc; pr++) { sum_elems_prev_sd_at_lev += _mesh._off_el[VV][pr*_NoLevels + Level + 1] - _mesh._off_el[VV][pr*_NoLevels + Level]; }

	    for (int iel = 0;
              iel <    _mesh._off_el[VV][off_proc + Level+1]
                     - _mesh._off_el[VV][off_proc + Level]; iel++) {
      for (uint is=0; is< _mesh._GeomEl.n_se[VV]; is++) {      
	sol_on_cells[cel*_mesh._GeomEl.n_se[VV] + is] = _bc_fe_kk[Level][iel + sum_elems_prev_sd_at_lev + ivar*_mesh._NoElements[VV][Level]]; //this depends on level!
      }
      cel++;
    }
  }
  
  std::ostringstream var_name; var_name << _var_names[ ivar + _VarOff[KK] ] << "_" << grname.str() << bdry_suffix;
  hsize_t dimsf[2]; dimsf[0] = NGeomObjOnWhichToPrint[KK]; dimsf[1] = 1;
  IO::print_Ihdf5(file_id,var_name.str(),dimsf,sol_on_cells);   
      
      delete [] sol_on_cells;
      
      } //end KK   
    
    } //end Level
    

     H5Fclose(file_id);

    return;
}




//===========================================
/////////////////// ELEM BC ////////////////////
///////////////////////////////
///    el_flag[0] = pressure integral (yes=1 or no=0)
///    el_flag[1] = stress integral
///    el_flag[2] = normal component of velocity (fixed=yes=1)
///   el_flag[3] = tangential component(-s)

//elembc
//this must only pick my _iproc and my Level
//so what is the length of this per proc and per level?
//clearly, here we do not know in which level we are, so we must
//compute the lengths for every processor and every level
//well, when we first allocate the vector, since it is serial,
//we allocate it to the maximum length
//then, it is when we FILL it that we have to take care of the correct range
//we'll have to do a map for each LEVEL, how do we do that?

//this map only takes into account the fine level, but how do i consider all the levels?
//for the dofs, you just do it for the FINE LEVEL, because then a FINE node is also a coarser node
//what could we do here?Once we pick the el_xm we must know at which level we are.
//Well,actually we can pick the boundary element by its iel, like one does for the connectivity
//if the iel is in that range, compute the middle point and fill the _elem_bc

//---------> WHAT IS THE PROBLEM ABOUT THIS FUNCTION?
//---- Why does the get_el_nod_conn(BB,..) work in GenMatRhsVB
//but HERE it doesnt?
//I think because here we are looping over ALL LEVELS and ALL SUBDOMAINS (procs)
//while in the GenMatRhsVB every processor
//does HIS OWN SUBDOMAIN (_iproc)
// and ONE LEVEL AT A TIME (Level)!
//instead, here every single processor does all of this!
//it is a SERIAL routine
//instead, in GenMatRhs every processor does HIS OWN PART, BOTH AT VOLUME and AT BOUNDARY
//do we have to loop over

//each processor does the loop over ALL the PROCESSORS, because this routine is serial

//compute the elem middle point
//we are on the boundary, so we are doing Quad9
//why was it always working with MORE THAN ONE LEVEL?
//it seems like the boundary connectivity map with ONE level and MORE THAN TWO PROCS doesnt work...

//here, I distinguished the levels.
//can i think of something like the dofs, where you can do all the things only at the FINE LEVEL?
//if a dof is fine, it is also coarse,
//if an element is fine, is it also coarse? no
//if an element is coarse, is it also fine? no
//if a flag is set on a fine element, then can i obtain the corresponding COARSE flag?
//you have to know the SONS of your FATHER element
//let us do a map with THREE indices:
//first: processor
//second: level
//third: iel
//well, the el_connectivity goes first by proc then by level,
//so we decide to follow the el_connectivity
//no, i think that i should put first the small range indices, then the larger ones
//i wont use more than 4-5 level
//i may want to use a lot of processors
//every processor at every level may have a lot of elements

//this routine is GENERAL. The only thing that is NOT GENERAL is the use of 4 element boundary flags
// and of ONE normal value and ONE tangential VALUE
//for now, we will leave things like this
void EqnBase::GenElBc()  {

     CurrElem       currelem(*this,_eqnmap);  
  
      uint space_dim = _mesh._dim;

    _elem_bc       =  new int**[_NoLevels];
    _elem_val_norm =  new double**[_NoLevels];
    _elem_val_tg   =  new double**[_NoLevels];

    for (uint Level=0; Level <_NoLevels;Level++)   { //loop over the levels

        _elem_bc[Level]       = new int*[_mesh._NoSubdom];//4*_mesh._NoElements[BB][Level] this was wrong, actually it is L + P*NoLevels, i wasnt considering the others but it was working! For instance in two D the numbers for the two procs are the same, here's why!
        _elem_val_norm[Level] = new double*[_mesh._NoSubdom];//_mesh._NoElements[BB][Level]
        _elem_val_tg[Level]   = new double*[_mesh._NoSubdom];//_mesh._NoElements[BB][Level]

        for (uint isubd=0;isubd<_mesh._NoSubdom;++isubd) {  //loop over the subdomains

            uint iel_b=_mesh._off_el[BB][ _NoLevels*isubd +  Level];
            uint iel_e=_mesh._off_el[BB][ _NoLevels*isubd +  Level+1];

            _elem_bc[Level][isubd]  = new int[/*4*/2*(iel_e-iel_b)]; /*normal and tangential*/ //4*_mesh._NoElements[BB][Level+isubd*_NoLevels] that was not correct
            _elem_val_norm[Level][isubd]  = new double[ 1*(iel_e-iel_b) ];  //_mesh._NoElements[BB][Level+isubd*_NoLevels]

            _elem_val_tg[Level][isubd]  = new double[ _number_tang_comps[space_dim - 1]*(iel_e-iel_b) ];   //_mesh._NoElements[BB][Level+isubd*_NoLevels]

            for (uint iel=0;iel < (iel_e - iel_b); iel++) {  //loop over the elems of that level&subdomain

                int surf_id=0;  //could do this outside also
                int     el_flag[2] = {0,0};
                std::vector<double> el_value(1 + _number_tang_comps[space_dim - 1],0.); //1 normal and 1 tangential or 1 normal and 3 tangential

                currelem.get_el_nod_conn_lev_subd(BB,Level,isubd,iel);
                currelem.get_el_ctr(BB);

                //read the bc's //the read forgets all levels and subdomains, it is only based on the MIDDLE POINT
                elem_bc_read(currelem._el_xm[BB],surf_id,&el_value[0],el_flag);

                std::cout << "Bdry " << surf_id << " normal: " <<  el_flag[0] << " tang: " <<  el_flag[1] << std::endl;

                //store the bc's
                _elem_bc[Level][isubd][/*4*/2*iel]   =  el_flag[0];
                _elem_bc[Level][isubd][/*4*/2*iel+1] =  el_flag[1];
                /*              _elem_bc[Level][isubd][4*iel+2] =  el_flag[2];
                              _elem_bc[Level][isubd][4*iel+3] =  el_flag[3];*/
                _elem_val_norm[Level][isubd][iel]   = el_value[0];
 
		for (uint i=0; i< _number_tang_comps[space_dim - 1]; i++) {
		  _elem_val_tg[Level][isubd][ _number_tang_comps[space_dim-1]*iel + i ] = el_value[1+i];
		}
// // //   WAS
// // // #if dimension==2
// // //                 _elem_val_tg[Level][isubd][iel]   = el_value[1];
// // // #elif dimension==3
// // //                 _elem_val_tg[Level][isubd][3*iel]   = el_value[1];
// // //                 _elem_val_tg[Level][isubd][3*iel+1]   = el_value[2];
// // //                 _elem_val_tg[Level][isubd][3*iel+2]   = el_value[3];
// // // 
// // // #endif

            }
        }

    }//end nolevels


    return;
}


void EqnBase::elem_bc_read(double */*el_xm*/,int& surf_id, double *value,int* el_flag) {

  std::cout << "AAAAAAAAAA you are using the default elem_bc_read which is basically MEANINGLESS" << std::endl;

  uint space_dim = _mesh._dim;

  surf_id=0;  

  el_flag[NN]=1;    //yes normal component
  el_flag[TT]=1;    //yes tang component
     
  value[NN]=0.;  //value of the normal 

   for (uint i=0; i < _number_tang_comps[space_dim - 1]; i++) value[1+i] = 0.;

// // //   WAS
// // // #if dimension==2
// // //   value[1]=0.;  //value of the tangential
// // // #elif dimension==3
// // //   value[1]=0.;  //value of the tangential
// // //   value[2]=0.;  //value of the  tangential
// // //   value[3]=0.;  //value of the tangential
// // // #endif

  return;
}


///////////////////////////////
void EqnBase::Bc_GetElFlagValLevSubd(const uint Level,const uint isubd,const uint iel,int* el_flag,double* el_value ) const {

    el_flag[0] =        _elem_bc[Level][isubd][/*4*/2*iel]   ;
    el_flag[1] =        _elem_bc[Level][isubd][/*4*/2*iel+1] ;
    /*    el_flag[2] =        _elem_bc[Level][isubd][4*iel+2] ;
        el_flag[3] =        _elem_bc[Level][isubd][4*iel+3] ;*/
    el_value[0] =  _elem_val_norm[Level][isubd][iel]   ;

  uint space_dim = _mesh._dim;

  for (uint i=0; i < _number_tang_comps[space_dim - 1]; i++) el_value[1+i] = _elem_val_tg[Level][isubd][ _number_tang_comps[space_dim - 1]*iel + i ];
    
// // //   WAS
// #if dimension==2
//     el_value[1] =    _elem_val_tg[Level][isubd][iel]   ;
// #elif dimension==3
//     el_value[1] =    _elem_val_tg[Level][isubd][3*iel]   ;
//     el_value[2] =    _elem_val_tg[Level][isubd][3*iel+1]   ;
//     el_value[3] =    _elem_val_tg[Level][isubd][3*iel+2]   ;
// #endif

    return;
}



////////////////////////////////
void EqnBase::clearElBc() {

    for (uint Level=0; Level <_NoLevels;Level++)   {

        for (uint isubd=0;isubd<_mesh._NoSubdom;++isubd) {
            delete [] _elem_bc[Level][isubd];
            delete [] _elem_val_norm[Level][isubd];
            delete [] _elem_val_tg[Level][isubd];
        }

        delete [] _elem_bc[Level];
        delete [] _elem_val_norm[Level];
        delete [] _elem_val_tg[Level];

    }

    delete [] _elem_bc;
    delete [] _elem_val_norm;
    delete [] _elem_val_tg;

    return;
}

/////////////////// ELEM BC ////////////////////




// ============================================================================
// we must set initial conditions for all levels and all types of dofs
// ok ic_read takes care of the initial conditions based on the coordinates.
// the point is that we have to locate the NODES for writing the FUNCTIONS,
// but also we need to locate the ELEMENTS for the cells.
// So I guess the ic_read will receive also the MIDDLE POINT of the ELEMENT
//let me add it to the arguments of ic_read, as last argument
// Then when you have a node based dof, you will use the node,
// when you have the elem based dof, you will use the cell center
// Of course, for your problem, you need to know THE ORDER of the SCALAR VARIABLES:
// what QQ you have, then what LL, then what KK...
// maybe we should find a way to do this automatically: define the list of 
//scalar variables and to each of them give a NUMBER which represents the ORDER.
// Of course it can only be up to the user to put the right function at the right place!

//TODO AAA: oh this is very interesting... the initial conditions are set IN PARALLEL,
// only by the current processor!

//TODO ok the point is this: we are looping over elements, and then over nodes in the element,
//so when we pick the element values we pick them more than once actually...
//but that's ok, we do like this so far, otherwise you would need two separate routines,
// one for the NODE dofs and one for the ELEM dofs, and if you want to switch a quantity
// from NODE to ELEM that may be better to have things together... well, of course you have to
// distinguish 
//So of course when you change FE you would need to use x_node instead of x_elem...
// That's why the best thing would be to PROVIDE, for ANY DOF, the ASSOCIATED DOF CARRIER coordinate,
// which could be either a node or an element!
//but anyway, let us go ahead like this now, it is fair already

//In order to set initial conditions, we may think in two ways:
//either loop over dofs, and for each dof find its associated DofObject,
// or loop over DofObjects, and find the dofs associated to them
//so far we are looping over DofObj

//Then, the number of dofs associated to a certain dof obj is variable... you may have
//only the quadratic, only the linear, both,... for element DofObj you only have the constants of course...

//here we are not looping over processors, we are INSIDE ONLY ONE PROCESSOR,
// so we are looping on the ELEMENTS of THIS PROCESSOR (there is only a local element loop)

//Basically in the PRINTING ROUTINES we have loops where the VARIABLES are EXTERNAL
// and inside each variables you loop over processors and then elements
//here instead each processor loops on its own elements, and inside each element
// we loop over all the types of variables

//For the printing we needed to do like this, because we need to print an array
// for EACH VARIABLE, so the single variable is the main actor
//For this routine the main goal is to fill all of x_old in PARALLEL manner,
// so the main actor must be the PROCESSOR and the ELEMENTS of the processor
//Probably if you wanted to do Printing with PARALLEL file system 
// then you would need to give priority to processors first...

//TODO: so far we have a node_dof that is explored based on the FINEST LEVEL for the NODES,
// and based on EACH LEVEL for the elements. 
//So, if we want to make it more general, we just put the offsets to be dependent on LEVEL,
// and for the nodes they are simply NOT DEPENDENT!

//TODO this function's goal is basically just to fill 
// the x vector in parallel, and then put it into x_old.
//Of course it does this in a GEOMETRIC WAY.
// xold is a vector of DOFS, so what we need are the 
// GEOMETRIC COORDINATES of the DOF OBJECTS (Node, Elem)
// associated to these dofs. In this way we can 


/// This function generates the initial conditions:
void EqnBase::GenIc() {

    const uint mesh_ord    = (int) _mesh._mesh_rtmap.get("mesh_ord");
  
    std::string   basepath = _files.get_basepath();
    std::string  input_dir = _files.get_frtmap().get("INPUT_DIR");
    std::string        ibc = DEFAULT_IBC;
    std::string     ext_h5 = DEFAULT_EXT_H5;
    std::ostringstream ibc_filexmf;
    ibc_filexmf << basepath << "/"<< input_dir << ibc << ext_h5;
    std::ifstream in(ibc_filexmf.str().c_str());

    if (!in) {

        CurrElem       currelem(*this,_eqnmap);  
     
        const uint  coords_fine_offset = _mesh._NoNodesXLev[_NoLevels-1];
        const uint  el_nnodes = _mesh._GeomEl._elnds[mesh_ord][VV];
        double*      xp = new double[_mesh._dim];
        double* u_value = new double[_n_vars];

       std::cout << "\n====================== GenIc:  Now we are setting them for all levels! ========================" << "\n \n";

    for (uint Level = 0; Level< _NoLevels; Level++) {

            uint iel_b = _mesh._off_el[VV][ _iproc*_NoLevels + Level ];
            uint iel_e = _mesh._off_el[VV][ _iproc*_NoLevels + Level + 1];

	    for (uint iel=0; iel < (iel_e - iel_b); iel++) {
	  
	        currelem.get_el_nod_conn_lev_subd(VV,Level,_iproc,iel);
                currelem.get_el_ctr(VV);
	
            // we are looping over the mesh nodes, but it's a fake loop because we do not depend on "i" for the elements
            for (uint i=0; i < el_nnodes ; i++) {
                int fine_node = _mesh._el_map[VV][ i + ( iel + iel_b )*el_nnodes ];
                for (uint idim=0; idim<_mesh._dim; idim++) xp[idim] = _mesh._xyz[ fine_node + idim*coords_fine_offset ];
                for (uint ivar=0; ivar<_n_vars; ivar++) u_value[ivar] = 0.;

                // ===================================================
                // user definition reading function ----------------
                ic_read(xp,u_value,currelem._el_xm[VV]);
                // -------------------------------------------------
                // ===================================================

                // Set the quadratic fields
                if ( i < _AbstractFE[QQ]->_ndof[VV] ) {
                    for (uint ivar=0; ivar<_nvars[QQ]; ivar++) {
		      int dof_pos_lev = _node_dof[Level][ fine_node + ivar*_DofNumLevFE[Level][QQ] + _DofOffLevFE[Level][QQ] ];
                        _x[Level]->set( dof_pos_lev, u_value[ivar + _VarOff[QQ]] );
		    }
                }
                // Set the linear fields
                if ( i < _AbstractFE[LL]->_ndof[VV] ) {
                    for (uint ivar=0; ivar<_nvars[LL]; ivar++) {
		      int dof_pos_lev = _node_dof[Level][ fine_node + ivar*_DofNumLevFE[Level][LL] + _DofOffLevFE[Level][LL] ];
                        _x[Level]->set( dof_pos_lev, u_value[ivar + _VarOff[LL]] );
                        }
		    }
		    
                if ( i < _AbstractFE[KK]->_ndof[VV] ) {
		  
		int sum_elems_prev_sd_at_lev = 0;
	    for (uint pr = 0; pr < _iproc; pr++) { sum_elems_prev_sd_at_lev += _mesh._off_el[VV][pr*_NoLevels + Level + 1] - _mesh._off_el[VV][pr*_NoLevels + Level]; }
	    
		  //elem my level
                    for (uint ivar=0; ivar<_nvars[KK]; ivar++) {
		      int elem_lev = iel + sum_elems_prev_sd_at_lev;
		      int dof_pos_lev = _node_dof[Level][ elem_lev + ivar*_DofNumLevFE[Level][KK] + _DofOffLevFE[Level][KK] ];
                        _x[Level]->set( dof_pos_lev, u_value[ivar + _VarOff[KK]] );
                         }
		    }
            }  //end i loop
            
        } // end of element loop

        _x[Level]->localize(*_x_old[Level]);
        _x_old[Level]->close();
	
    } //end Level
    
        delete [] u_value;
        delete [] xp;

#ifdef DEFAULT_PRINT_INFO
        std::cout << "\n GenIc(Base): Initial solution defined by ic_read" << "\n \n";
#endif
    }
    else {// -------------------- file reading --> data_in/case.h5
        std::cout << "^^^^^ WE HAVE TO CHECK BECAUSE WE ADDED CONSTANT ELEMENTS ^^^^^^" << std::endl; abort();

        ReadVector(ibc_filexmf.str());
        std::cout << "\n GenIc(Base): Initial solution defined by " <<  ibc_filexmf.str().c_str() << "\n \n";
    }

    in.close();

    return;
}




/// initial conditions  from function
void EqnBase::ic_read(double /*xp*/[],double ic[], double /*el_xm*/[]) {
    for (uint ivar=0;ivar<_n_vars;ivar++) ic[ivar]=0.;
    return;
}



// ======================================================
/// This function controls the time step operations:
//first, prepare the A matrices for all levels
// and the b rhs for the fine level
//then call the multigrid solver
//eventually, update the old solution

//this is the NONLINEAR loop
//it is settled as a multigrid, but for 0 levels it is a normal single grid
//it returns the l2 norm of the difference between two nonlinear steps, x - x_old

//time is passed only because the matrices and rhs can depend on time
//   its value is given by TimeLoop

/// AAAAAAA TODO I have to use this for the PENALTY METHOD!!!
/// In general, it is safer to call both of them because
/// you may fill both Ke and Fe in both and they are needed at all levels
//     GenMatRhs(time,Level,0); // matrix
//     GenRhsB(time,Level);
//yes with NESTED! //with nonnested the rhs is overwritten,you're solving for another problem in the other levels
//with non-nested you do not fill the rhs (mode=0) nor you fill with the boundary integrals (GenRhsB commented)
//pay attention to this!
//i dont wanna use mode any longer!
//mode is used to avoid filling the RHS
//because the multigrid has another RHS at the non-fine levels.
//But actually the RHS is AUTOMATICALLY OVERWRITTEN for the non-fine levels,
//because it is always the RESTRICTION of the RESIDUAL from the UPPER LEVEL

//TODO Pay attention here, we are filling BOTH A and B at ALL LEVELS,
// But then in the multigrid algorithm b for some levels will be OVERWRITTEN!!!
//The point is that you still have to fill A, also with the boundary conditions if necessary (SEE PENALTY DIRICHLET BC)
// Of course, the boundary conditions must be satisfied at all levels

/// Basically, the single time step consists in  ASSEMBLING + SOLVING

double EqnBase::MGTimeStep(const double time, const uint iter) {

    std::cout  << std::endl << " Solving " << _eqname << " , step " << iter << std::endl;

    ///A0) Put x_old into x_oold
    *_x_oold[_NoLevels-1] = *_x_old[_NoLevels-1];

    /// A) Assemblying 
#if  DEFAULT_PRINT_TIME==1
    std::clock_t start_time=std::clock();
#endif

    for (uint Level = 0 ; Level < _NoLevels; Level++) {

        _A[Level]->zero();
        _b[Level]->zero();

        GenMatRhsVB(VV,time,Level);
        GenMatRhsVB(BB,time,Level);

#ifdef DEFAULT_PRINT_INFO
        _A[Level]->close();
        double ANorm=_A[Level]->l1_norm();
	
//	_A[Level]->print_graphic(true); TODO should pass this true or false as a parameter
	
        std::cout << " ANorm l1 " << Level << " "  << ANorm  << std::endl;
#endif
    }

#if    DEFAULT_PRINT_TIME==1
    std::clock_t end_time=std::clock();
    std::cout << " ================ Assembly time = " << double(end_time- start_time) / CLOCKS_PER_SEC
              << " s " << std::endl;
#endif

///std::cout << " $$$$$$$$$ Prepared A for all levels and b for the fine level $$$$$$" << std::endl;
// The matrices R and P for all levels were already prepared at read-time
// and they do not depend on  time (because no adaptive refinement is performed)
// They only depend on the nodes (dofs!)

/// C) Reinitialization of  the ILU (if ILUINIT defined),
// #ifdef ILUINIT
//   if (iter%ILUINIT==0) {
//     std::cout << " ILU Reinitialized: " << std::endl;
//     for (int Level = NoLevels - 1; Level > 0; Level--) *L[Level].ILUExists = _LPFalse;
//   }
// #endif

    /// D) Solution of the linear MGsystem
#ifdef DEFAULT_PRINT_TIME
        std::clock_t start_time_sol = std::clock();
#endif
      
        MGSolve(DEFAULT_EPS_LSOLV, DEFAULT_MAXITS_LSOLV);
	
// //     for (uint Level = 0 ; Level < _NoLevels; Level++)  { 
// //       _solver[Level]->solve(*_A[Level],*_x[Level],*_b[Level],1.e-6,40);
// //            _x[Level]->localize(*_x_old[Level]);   // x_old = x
// //     }

#if    DEFAULT_PRINT_TIME==1
    std::clock_t end_time_sol = std::clock();
    std::cout << " ================ Solver time = " << double(end_time_sol- start_time_sol) / CLOCKS_PER_SEC
              << " s "<< std::endl;
#endif

      
/// std::cout << "$$$$$$$$$ Computed the x with the MG method $$$$$$$" << std::endl;

    /// E) Update of the old solution at the top Level
    _x[_NoLevels-1]->localize(*_x_old[_NoLevels-1]);   // x_old = x
#ifdef DEFAULT_PRINT_INFO
    std::cout << "$$$$$$$$$ Updated the x_old solution $$$$$$$$$" << std::endl;
#endif
/// std::cout << "$$$$$$$$$ Check the convergence $$$$$$$" << std::endl;

    _x_tmp[_NoLevels-1]->zero();
    _x_tmp[_NoLevels-1]->add(+1.,*_x_oold[_NoLevels-1]);
    _x_tmp[_NoLevels-1]->add(-1.,*_x_old[_NoLevels-1]);
    // x_oold -x_old =actually= (x_old - x)
    //(x must not be touched, as you print from it)
    //x_oold must not be touched ! Because it's used later for UPDATING Becont!
    //so you must create a temporary vector necessarily.

    _x_tmp[_NoLevels-1]->close();
    double deltax_norm=_x_tmp[_NoLevels-1]->l2_norm();
    std::cout << " $$$$$$ " << _eqname << " error l2 " << deltax_norm << std::endl;
    std::cout << " $$$$$$ " << _eqname << " error linfty " << _x_tmp[_NoLevels-1]->linfty_norm() << std::endl;
//AAA when the vectors have nan's, the norm becomes zero!
//when the residual norm in pre and post smoothing is too big,
//then it doesnt do any iterations, the solver doesnt solve anymore, so the solution remains frozen

    return deltax_norm;  //TODO do we have to be based on l2norm or linfty norm???
}



// =================================================================
//if the residual norm is small enough,exit the cycle, and so also the MGSolve
//this is also the check for the single level solver
//what is the meaning of having multiple cycles for single-grid solver?
//they are all like just a single linear solver loop, where convergence has already been reached,
// and you do another check after the previous one in the linear solver loop
//the big question is:
//TODO why dont you do "res_fine < Eps1"
// instead of  "res_fine < Eps1*(1.+ bNorm_fine)" ???
//because one is for the absolute error and another one is for the relative error
/// This function solves the discrete problem with multigrid solver
void EqnBase::MGSolve(double Eps1,          // tolerance for the linear solver
                      int MaxIter,           // n iterations
                      const uint Gamma,     // Control V W cycle
                      const uint Nc_pre,    // n pre-smoothing cycles
                      const uint Nc_coarse, // n coarse cycles
                      const uint Nc_post    // n post-smoothing cycles
                     ) {

#ifdef DEFAULT_PRINT_INFO
    std::cout << "######### BEGIN MG SOLVE ########" << std::endl;
#endif
    double res_fine;

    _b[_NoLevels-1]->close();
    double bNorm_fine =     _b[_NoLevels-1]->l2_norm();
    _x_old[_NoLevels-1]->close();
    double x_old_fine = _x_old[_NoLevels-1]->l2_norm();

#ifdef DEFAULT_PRINT_INFO
    std::cout << " bNorm_fine l2 "     <<  bNorm_fine                     << std::endl;
    std::cout << " bNorm_fine linfty " << _b[_NoLevels-1]->linfty_norm()  << std::endl;
    std::cout << " xold_fine l2 "      <<  x_old_fine                     << std::endl;
#endif

    // FAS Multigrid (Nested) ---------
    bool NestedMG=false;
    if (NestedMG) {
        _x[0]->zero();
        MGStep(0,1.e-20,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post);

        //smooth on the coarse level WITH PHYSICAL b !
        //and compute the residual

        for (uint Level = 1; Level < _NoLevels; Level++) {

            _x[Level]->matrix_mult(*_x[Level-1],*_Prl[Level]);  //**** project the solution

            res_fine = MGStep(Level,Eps1,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post);

        }
    } // NestedMG

    // V or W cycle
    int cycle = 0;
    bool exit_mg = false;

    while (!exit_mg && cycle<MaxIter) {

///std::cout << "@@@@@@@@@@ BEGIN MG CYCLE @@@@@@@@"<< std::endl;

///std::cout << "@@@@@@@@@@ start on the finest level @@@@@@@@"<< std::endl;

        res_fine = MGStep(_NoLevels-1,Eps1,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post);
        // MGCheck(_NoLevels-1); // check projection-restriction

///std::cout << "@@@@@@@@@@ back to the finest level @@@@@@@@"<< std::endl;

///std::cout << "@@@@@@@@@@ END MG CYCLE @@@@@@@@"<< std::endl;

        std::cout << "@@@@@@@@@@ CHECK THE RESIDUAL NORM OF THE FINEST LEVEL @@@@@@@@"<< std::endl;

        std::cout << "res_fine: " << res_fine << std::endl;
        std::cout << "bNorm_fine: " << bNorm_fine << std::endl;

        if (res_fine < Eps1*(1. + bNorm_fine)) exit_mg = true;

        cycle++;

#ifdef DEFAULT_PRINT_INFO
        std::cout << " cycle= " << cycle   << " residual= " << res_fine << " \n";
#endif

    }

#ifdef DEFAULT_PRINT_INFO
    std::cout << "######### END MG SOLVE #######"<< std::endl;
#endif
    return;
}

// ====================================================================
/// This function does one multigrid step
// solve Ax=b
// compute the residual res=b- Ax
// restrict the residual R*res
//notice that the A and b for the POST-smoothing are the same
//as for the pre-smoothing


double EqnBase::MGStep(int Level,            // Level
                       double Eps1,          // Tolerance
                       int MaxIter,          // n iterations - number of mg cycles
                       const uint Gamma,     // Control V W cycle
                       const uint Nc_pre,    // n pre-smoothing smoother iterations
                       const uint Nc_coarse, // n coarse smoother iterations
                       const uint Nc_post    // n post-smoothing smoother iterations
                      ) {


    std::pair<uint,double> rest;

    if (Level == 0) {
///  std::cout << "************ REACHED THE BOTTOM *****************"<< std::endl;

#ifdef DEFAULT_PRINT_CONV
        _x[Level]->close();
        double xNorm0=_x[Level]->linfty_norm();
        _b[Level]->close();
        double bNorm0=_b[Level]->linfty_norm();
        _A[Level]->close();
        double ANorm0=_A[Level]->l1_norm();
        std::cout << "Level " << Level << " ANorm l1 " << ANorm0 << " bNorm linfty " << bNorm0  << " xNormINITIAL linfty " << xNorm0 << std::endl;
#endif

        rest = _solver[Level]->solve(*_A[Level],*_x[Level],*_b[Level],DEFAULT_EPS_LSOLV_C,Nc_coarse);  //****** smooth on the coarsest level

#ifdef DEFAULT_PRINT_CONV
        std::cout << " Coarse sol : res-norm: " << rest.second << " n-its: " << rest.first << std::endl;
        _x[Level]->close();
        std::cout << " Norm of x after the coarse solution " << _x[Level]->linfty_norm() << std::endl;
#endif

        _res[Level]->resid(*_b[Level],*_x[Level],*_A[Level]);      //************ compute the coarse residual

        _res[Level]->close();
        std::cout << "COARSE Level " << Level << " res linfty " << _res[Level]->linfty_norm() << " res l2 " << _res[Level]->l2_norm() << std::endl;

    }

    else {

///  std::cout << "************ BEGIN ONE PRE-SMOOTHING *****************"<< std::endl;
#ifdef DEFAULT_PRINT_TIME
        std::clock_t start_time=std::clock();
#endif
#ifdef DEFAULT_PRINT_CONV
        _x[Level]->close();
        double xNormpre=_x[Level]->linfty_norm();
        _b[Level]->close();
        double bNormpre=_b[Level]->linfty_norm();
        _A[Level]->close();
        double ANormpre=_A[Level]->l1_norm();
        std::cout << "Level " << Level << " ANorm l1 " << ANormpre << " bNorm linfty " << bNormpre  << " xNormINITIAL linfty " << xNormpre << std::endl;
#endif

        rest=_solver[Level]->solve(*_A[Level],*_x[Level],*_b[Level],DEFAULT_EPS_PREPOST, Nc_pre); //****** smooth on the finer level

#ifdef DEFAULT_PRINT_CONV
        std::cout << " Pre Lev: " << Level << ", res-norm: " << rest.second << " n-its: " << rest.first << std::endl;
#endif
#ifdef DEFAULT_PRINT_TIME
        std::clock_t end_time=std::clock();
        std::cout << " time ="<< double(end_time- start_time) / CLOCKS_PER_SEC << std::endl;
#endif

        _res[Level]->resid(*_b[Level],*_x[Level],*_A[Level]);//********** compute the residual

///    std::cout << "************ END ONE PRE-SMOOTHING *****************"<< std::endl;

///    std::cout << ">>>>>>>> BEGIN ONE DESCENT >>>>>>>>>>"<< std::endl;

        _b[Level-1]->matrix_mult(*_res[Level],*_Rst[Level-1]);//****** restrict the residual from the finer grid ( new rhs )
        _x[Level-1]->close();                                  //initial value of x for the presmoothing iterations
        _x[Level-1]->zero();                                  //initial value of x for the presmoothing iterations

        for (uint g=1; g <= Gamma; g++)
            MGStep(Level-1,Eps1,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post); //***** call MGStep for another possible descent

//at this point you have certainly reached the COARSE level

///      std::cout << ">>>>>>>> BEGIN ONE ASCENT >>>>>>>>"<< std::endl;
#ifdef DEFAULT_PRINT_CONV
        _res[Level-1]->close();
        std::cout << "BEFORE PROL Level " << Level << " res linfty " << _res[Level-1]->linfty_norm() << " res l2 " << _res[Level-1]->l2_norm() << std::endl;
#endif

        _res[Level]->matrix_mult(*_x[Level-1],*_Prl[Level]);//******** project the dx from the coarser grid
#ifdef DEFAULT_PRINT_CONV
        _res[Level]->close();
        //here, the _res contains the prolongation of dx, so the new x is x_old + P dx
        std::cout << "AFTER PROL Level " << Level << " res linfty " << _res[Level]->linfty_norm() << " res l2 " << _res[Level]->l2_norm() << std::endl;
#endif
        _x[Level]->add(*_res[Level]);// adding the coarser residual to x
        //initial value of x for the post-smoothing iterations
        //_b is the same as before
///   std::cout << "************ BEGIN ONE POST-SMOOTHING *****************"<< std::endl;
        // postsmooting (Nc_post)
#ifdef DEFAULT_PRINT_TIME
        start_time=std::clock();
#endif
#ifdef DEFAULT_PRINT_CONV
        _x[Level]->close();
        double xNormpost=_x[Level]->linfty_norm();
        _b[Level]->close();
        double bNormpost=_b[Level]->linfty_norm();
        _A[Level]->close();
        double ANormpost=_A[Level]->l1_norm();
        std::cout << "Level " << Level << " ANorm l1 " << ANormpost << " bNorm linfty " << bNormpost << " xNormINITIAL linfty " << xNormpost << std::endl;
#endif

        rest=_solver[Level]->solve(*_A[Level],*_x[Level],*_b[Level],DEFAULT_EPS_PREPOST,Nc_post);  //***** smooth on the coarser level

#ifdef DEFAULT_PRINT_CONV 
        std::cout<<" Post Lev: " << Level << ", res-norm: " << rest.second << " n-its: " << rest.first << std::endl;
#endif
#ifdef DEFAULT_PRINT_TIME
        end_time=std::clock();
        std::cout<< " time ="<< double(end_time- start_time) / CLOCKS_PER_SEC << std::endl;
#endif

        _res[Level]->resid(*_b[Level],*_x[Level],*_A[Level]);   //*******  compute the residual

///    std::cout << "************ END ONE POST-SMOOTHING *****************"<< std::endl;

    }

    _res[Level]->close();

    return  rest.second;  //it returns the residual norm of whatever level you are in
    //TODO if this is the l2_norm then also the nonlinear solver is computed in the l2 norm
    //WHAT NORM is THIS?!? l2, but PRECONDITIONED!!!

}

// =============================================
/// Check for Prolong and Restr Operators
void EqnBase::MGCheck(int Level) const {

    _x[Level-1]->matrix_mult(*_x[Level],*_Rst[Level-1]);
    _x[Level]  ->matrix_mult(*_x[Level-1],*_Prl[Level]);
    return;
}


// =========================================
/// This function is the main initialisation function:
//here, the vectors involved in the system are initialized per level
//a lot of things are performed here
//let us split this into two parts
//the multigrid part, which depends on the files also
//and the MatVec part
//actually, also the matrix depends on the file due to the sparsity pattern
//clearly, this function is called after InitMeshToDof

//==================
///this function DEPENDS in _iproc!!!
/// Every process has a portion of the structures of every level,
///but of course he knows also what is the GLOBAL dimension.
//Here the matrix is initialized, so you give the dimensions,
//but i would do altogether in the ReadMatrix files


void EqnBase::initMGOps() {

    std::string     basepath = _files.get_basepath();
    std::string   output_dir = _files.get_frtmap().get("OUTPUT_DIR");
    std::string  outtime_dir = _files.get_frtmap().get("OUTTIME_DIR");
    std::string     f_matrix = DEFAULT_F_MATRIX;
    std::string       f_rest = DEFAULT_F_REST;
    std::string       f_prol = DEFAULT_F_PROL;
    std::string       ext_h5 = DEFAULT_EXT_H5;

    std::ostringstream filename;
    std::string filename_base;
    filename_base =  basepath + "/" + output_dir + "/" + outtime_dir + "/";
    
        filename.str("");     filename << filename_base << f_rest << ext_h5;
        ReadRest(filename.str());
    
        filename.str("");   filename << filename_base << f_matrix << ext_h5;
        ReadMatrix(filename.str());
  
        filename.str("");    filename << filename_base << f_prol << ext_h5;
        ReadProl(filename.str());

    return;
}



//=============================
//This function depends on _iproc
//remember that a new does not vanish when inside a for loop,
//because the variable was declared outside
// it is NOT like a DECLARATION INSIDE

//The point is that now we must be ready to READ the matrix,
//also the part with CONSTANT ELEMENTS.
//so we have to figure out how to put all the things together.
//basically we need to remember that we have a bunch of
// quantities that may be QQ, LL, KK,
// so eventually we'll have a certain number of QQ, LL and KK
//now, we simply have to decide HOW to ADD THEM in the MATRIX.
// so far, I think we can decide to put
//FIRST all the QQ, then LL, then KK.

// The DofMap was already initialized much earlier than now.
// As a matter of fact, it was needed for setting the BC flags
// It would be nice not to have a lot of places where _Dim is recomputed...

// initDofDim
// initDof
// initVectors

//The GRAPH, or SPARSITY PATTERN, with PETSC functions can be filled JUST WITH ZEROS;
//but with LASPACK  you need do put the DOF VALUES INSIDE.

// so, I got an error that seems to tell me that I basically didn't fill the graph
//correctly. It seems like one processor fills it and the other ones do not,
//half of the lines seem to be unfilled.

//TODO the std::cerr is not redirected to file, think of redirecting also it to file.
// The idea is that every output to terminal must be redirected to file.
// so printf must not be used 

// Every processor will call this function

//Now, the dimension of "graph" is global,
// so it is the dimension of the whole

//TODO: idea about SYNCHRONIZING OUTPUTS.
//I would really like all the outputs OF ALL THE LIBRARIES (HFD5, PETSC, MPI, LIBMESH) involved in FEMuS
//to be REDIRECTED to run_log, altogether.
//What is the best way to make it possible?
//I guess from shell, in principle. The point is how to generate the string for the time for the folder!
// I should use bash functions, make an environment variable, and let all the processes read that environment variable.
// And also, i would like to be able to choose whether to DUPLICATE them, 
// so both to file and to std output, or only to file, or the other cases.

//I need to have the TIME STRING. Bash has the date() function, I want to have it my way
//OK, THAT'S IT:  date +%Y_%m_%d-%H_%M_%S


  //TODO what is this mlinit used for?!?
                      // it is used to determine with WHAT DOF NUMBER the CURRENT PROCESSOR _iproc BEGINS.
                      // NOTICE THAT THE _NODE_DOF MAP IS FILLED COMPLETELY BY ALL PROCESSORS, 
                      // BUT HERE IT IS READ IN PARALLEL BY EACH SEPARATE PROCESSOR !!!
                      //THEREFORE EACH PROCESSOR WILL START WITH ITS OWN INITIAL_DOF AND FINISH WITH ITS OWN FINAL_DOF.
                      //The difference between NODES and ELEMENTS is that for nodes
                      // we keep the offset as the FINEST NODES, because somehow 
                      //we have to keep the correspondence with the NODES OF THE MESH, 
                      //especially when we have to PRINT ON THE MESH
                      //as a matter of fact, we have to bring everything to the FINE MESH
                      //so we can view what happens on the FINE MESH
                      //every processor will fill QUADRATIC, LINEAR and CONSTANT DOFS
                      //being that these dofs belong to a given processor, they will be CONTIGUOUS in the matrix
                      // So you'll have
                      
                      // QQ LL KK || QQ LL KK || QQ LL KK 
                      //< proc0  >||< proc1  >||< proc2 > 
                      
                      // Now, consider the level to be FROZEN, we concentrate on the SUBDOMAINS
                      // Every processor will begin filling the rows at some point.
                      // Every processor has a range of rows, which is mrow_lev_proc_t
                      // Within this range you'll have to put QQ, LL and KK dofs
                      //so within this range you'll have a number of dofs for every fe given by mrow_lev_proc[r],
                      // of course each of them being multiplied by the number of variables of that type
                      // Now, you have to compare how you EXPLORE the NODE_DOF MAP and the MATRIX ROWS
                      //the node_dof map is basically separated into OFFSETS of GEOMETRICAL ENTITIES,
                      //so these ranges are BASED ON THE MESH GEOMETRICAL ENTITIES.
                      // every GEOMETRICAL ENTITY (NODE or ELEMENT) can have ONE OR MORE ASSOCIATED DOFS.
                      //NODES CAN HAVE QUADRATIC DOFS or LINEAR DOFS associated
                      //ELEMENTS CAN HAVE CONSTANT DOFS associated.
                      //so, the node dof is constructed by picking
                      // FIRST NODE QUADRATIC DOFS, for all the quadratic unknowns of the system
                      //  THEN NODE LINEAR DOFS, for all the linear unknowns of the system
                      //  THEN ELEMENT CONSTANT DOFS, for all the constant unknowns of the system
                      // so "QUADRATIC NODES" can be defined as "NODES to which QUADRATIC VARIABLES are associated",
                      //  "LINEAR NODES" means                  "NODES to which LINEAR variables are associated"
                      // ELEMENTS are "GEOMETRICAL ENTITIES to which CONSTANT variables are associated"
                      // So the node dof is divided by 
                      // "GEOMETRICAL ENTITIES for QUADRATIC DOFS" = "QUADRATIC NODES", the length of which must be multiplied by the "NUMBER OF QUADRATIC DOF VARIABLES"
                      // "GEOMETRICAL ENTITIES for LINEAR DOFS"    = "LINEAR NODES",    the length of which must be multiplied by the "NUMBER OF LINEAR DOF VARIABLES"
                      // "GEOMETRICAL ENTITIES for CONSTANT DOFS"  = "ELEMENTS"         the length of which must be multiplied by the "NUMBER OF CONSTANT DOF VARIABLES"
                      //This of course assumes that every GEOMETRICAL ENTITY has the SAME NUMBER of ASSOCIATED DOFS for every FE type.  
                      //So you may say that the node_dof is ordered by GEOMETRICAL ENTITIES first, and WITHIN EACH GEOMETRICAL ENTITY you find ALL THE VARIABLES of THAT FE TYPE.
                      // This is because actually we have a ONE-TO-ONE CORRESPONDENCE between GEOMETRICAL ENTITY and ASSOCIATED DOF.
                      // In a sense, we start with a given DOF family, and, going backwards, we ask ourselves "WHAT IS THE SET OF GEOMETRICAL ENTITIES NATURALLY ASSOCIATED to THIS DOF FE Family?"
                      //For instance, we might have QUAD8 and QUAD9, and associate BOTH OF THEM to a "QUAD9 NODES".In that case we would have EMPTY points,
                      // and it is actually the same that happens when we have the LINEAR DOF FAMILY and we decide to associate the "QUADRATIC MESH NODES" to it for various reasons: multigrid, parallel computing, etc...)
                      //So our current situation is:
                      // QUADRATIC FE DOFS ---> QUADRATIC MESH NODES (at the FINE level)
                      // LINEAR FE DOFS    ---> QUADRATIC MESH NODES (at the FINE level)
                      // CONSTANT FE DOFS  ---> MESH ELEMENTS  (at the CURRENT level)
                      //So, we see that for various reasons we may actually have EMPTY SPACES for two reasons:
                      // 1) BECAUSE we are at some level which is NOT THE FINE
                      // 2) even at the FINE Level, we have LINEAR FINE DOFS built on top of  QUADRATIC FINE NODES
                      
                      //so, _nvars[QQ] is the number of QUADRATIC FE variables
                      
  //NOW, THE POINT IS THIS: FOR EVERY SUBDOMAIN, Do we want to have ALL DOFS CONTIGUOUS or NOT?
  //Yes, of course we want, so in the same subdomain we would have QUADRATIC DOFS, THEN LINEAR DOFS, THEN CONSTANT DOFS.
  //Otherwise, we would have that dofs of the SAME SUBDOMAIN are DISJOINT.
  //Keeping all the dofs of the same subdomain CONTIGUOUS, the equations do not change, they will be the same, 
  // but in this way we will have AS FEW OFF-PROCESSOR COLUMNS as POSSIBLE!
  // So, we do need to BUILD the NODE-DOF MAP APPROPRIATELY for the PARALLEL CASE!
  //Of course, in a SERIAL ENVIRONMENT everything works fine
  
//=======================
//The dimension of the Graph is the total number of rows
// This Graph is the TOTAL GRAPH for ALL THE MATRIX
//The Graph is basically a std::VECTOR of std::VECTORS,
//   therefore you can use the resize() function both on graph and on graph[].
//Also, by the way, with std::vectors you have the overloading of the operator[]
//=====================

// Now, the point is: What is mlinit[]?
// Now, based on the mesh Node and Element numbering, we built the node_dof map.
//Now we have to Read the ONE-VARIABLE-MATRICES for EVERY LEVEL, 
//and USE THEM to BUILD THE WHOLE MATRIX WITH ALL THE VARIABLES and SO ON.
//Now, we already know how the dof map is. So we gave a number to every dof.
// So after doing the NODE ORDERING and ELEMENT ORDERING (VV and BB) in the GENCASE,
// now it is time to use those orderings consequently and have the 
//DOF ORDERING.
//The node map was built in such a way that ALL the DOFS associated to DofObjects (Elems or Nodes)
//of a given subdomain be CONTIGUOUS in the DOF ROW INDEXING.
//Now, given this DOF ORDERING, every PROCESSOR WILL TAKE CARE of ITS OWN DOFS,
// which are contiguous by construction of the node_dof map;
// therefore every processor will pick that block of contiguous row numbers
// and to each of them it will associate the corresponding LENGTH and OFFLENGTH

// Notice how GRAPH has GLOBAL dimension, but SINCE THIS READ MATRIX is PARALLEL,
//Then each processor will fill only PART OF IT. So, actually it makes no sense to make graph be GLOBAL
// and then have each processor only fill one part of it.

//ReadMatrix NEEDS the NODE_DOF to be filled BEFORE
  //Now, we have to understand these steps:
  //1) How we GENERATE the SINGLE DOF VARIABLE sparsity pattern in Gencase;
  //2) How we BUILD the NODE_DOF MAP in PARALLEL based on these SINGLE VARIABLE SPARSITY PATTERNS;
  //3) How we READ the MATRIX SPARSITY PATTERNS and COMBINE THEM WITH the NODE_DOF.
  
  //Basically, you have to COMBINE the "DOFS" (Matrix.h5) with the "MESH" (mesh.h5).
                      
// The node_dof of course takes into account the NUMBER of VARIABLES for each FE type
// So, if the node_dof was not filled correctly, then when you 


void EqnBase::ReadMatrix(const  std::string& namefile) {

  _A.resize(_NoLevels);

    for (uint Level = 0; Level< _NoLevels; Level++) {

      std::ostringstream groupname_lev; groupname_lev <<  "LEVEL" << Level;
  
    uint*** rowcln;    //number of rows and columns for one variable couple
    uint*** length_row; //for every FE row and FE column
    uint*** length_offrow;
    uint*** pos_row;

    rowcln        = new uint**[QL];
    length_row    = new uint**[QL];
    length_offrow = new uint**[QL];
    pos_row       = new uint**[QL];

    hid_t  file = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

    for (int r=0;r<QL;r++) {  //row FE type

        rowcln[r]        = new uint*[QL];
        length_row[r]    = new uint*[QL];
        length_offrow[r] = new uint*[QL];
        pos_row[r]       = new uint*[QL];

        for (int c=0;c<QL;c++) {   //col FE type

            std::ostringstream fe_couple;
            fe_couple << "_F" << r << "_F" << c;

            // dimension
            std::ostringstream dim_name;
            dim_name << groupname_lev.str() << "/" << "DIM" << fe_couple.str();
            rowcln[r][c] = new uint[2];
            hid_t dataset=H5Dopen(file, dim_name.str().c_str(), H5P_DEFAULT);
            H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,rowcln[r][c]);

            //row length
            std::ostringstream len_name;
            len_name << groupname_lev.str() << "/" << "LEN" << fe_couple.str();
            length_row[r][c]=new uint[ rowcln[r][c][0]+1 ];
            dataset=H5Dopen(file,len_name.str().c_str(), H5P_DEFAULT);
            H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,length_row[r][c]);

            // matrix off diagonal
            std::ostringstream offlen_name;
            offlen_name << groupname_lev.str() << "/" << "OFFLEN" << fe_couple.str();
            length_offrow[r][c]=new uint[ rowcln[r][c][0]+1 ];
            dataset=H5Dopen(file,offlen_name.str().c_str(), H5P_DEFAULT);
            H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,length_offrow[r][c]);

            // matrix pos //must stay AFTER reading length_offrow
            std::ostringstream pos_name;
            pos_name << groupname_lev.str() << "/" << "POS" << fe_couple.str();
            pos_row[r][c]=new uint[ length_row[r][c][rowcln[r][c][0]] ];
            dataset=H5Dopen(file,pos_name.str().c_str(), H5P_DEFAULT);
            H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,pos_row[r][c]);

        } //end col
    } //end row
    
    H5Fclose(file);

//============================================================================
//============ compute things for the sparsity pattern =======================
//============================================================================
    
    int NoLevels  = _mesh._NoLevels;

    uint mrow_glob_t = 0;
    for (int fe=0; fe<QL; fe++) mrow_glob_t += _nvars[fe]*rowcln[fe][fe][0];
    uint ncol_glob_t = mrow_glob_t;

    uint mrow_lev_proc_t  = 0;
    for (int fe=0; fe<QL; fe++)  mrow_lev_proc_t +=  _DofLocLevProcFE[Level][_iproc][fe]*_nvars[fe];
    uint ncol_lev_proc_t  = mrow_lev_proc_t;

    uint DofObjInit_lev_PrevProcs[QL];  //TODO what is this? it is the ROW INDEX at which to begin for every processor
                      
     for (int r=0; r<QL; r++)     DofObjInit_lev_PrevProcs[r] = 0;
         
    for (uint isubd=0; isubd<_iproc; isubd++) {
        DofObjInit_lev_PrevProcs[QQ] += _DofLocLevProcFE[Level][isubd][QQ];
        DofObjInit_lev_PrevProcs[LL] += _DofLocLevProcFE[Level][isubd][LL];
        DofObjInit_lev_PrevProcs[KK] += _DofLocLevProcFE[Level][isubd][KK];
    }    
    
    
    _A[Level] = SparseMatrix::build().release();
// //     _A[Level]->init(_Dim[Level],_Dim[Level], mrow_lev_proc_t, mrow_lev_proc_t); //TODO BACK TO a REASONABLE INIT

    Graph graph;
    graph.resize(mrow_glob_t);
    graph._m  = mrow_glob_t;
    graph._n  = ncol_glob_t;
    graph._ml = mrow_lev_proc_t;
    graph._nl = ncol_lev_proc_t;
    graph._ml_start = _node_dof[Level][ _mesh._off_nd[QQ][_iproc*NoLevels] ];
    // TODO is this used? I guess it is used by update_sparsity_pattern !
    // Every subdomain has a local set of dofs, and these dofs start at a specific point.
    // Now, remember that _mesh._off_nd[QQ] should only be used for computing offsets, so differences.
    // Here, it is used ALONE, because it gives you the NODE (in femus ordering) AT WHICH THE CURRENT SUBDOMAIN BEGINS,
    //and then from _node_dof[Level] (which was already constructed) you get THE LOCAL DOF AT THE CURRENT LEVEL TO START FROM.
    // Clearly, pay attention when you add elements, because in that case you would need to REDO the _node_dof map !!!
    //TODO: also, what happens if you have a system with ONLY ELEMENT BASED DOFS?!?
    
    
    int FELevel[QL];
    FELevel[QQ] = Level;
    FELevel[LL] = (Level+_NoLevels)%(_NoLevels+1); //This is the map for the level of the LINEAR DOFS
    FELevel[KK] = Level;

    int off_onevar[QL];
    off_onevar[QQ] = _mesh._NoNodesXLev[_NoLevels-1];
    off_onevar[LL] = _mesh._NoNodesXLev[_NoLevels-1];
    off_onevar[KK] = _mesh._NoElements[VV][Level];
    
    uint  off_EachFEFromStart[QL];
    off_EachFEFromStart[QQ] = 0;
    off_EachFEFromStart[LL] = _nvars[QQ]*off_onevar[QQ];
    off_EachFEFromStart[KK] = _nvars[QQ]*off_onevar[QQ]+_nvars[LL]*off_onevar[LL];

 //==============   
     for (int r=0;r<QL;r++) {
      for (uint ivar=0; ivar<_nvars[r]; ivar++) {

        for (uint DofObj_lev = DofObjInit_lev_PrevProcs[r]; DofObj_lev < DofObjInit_lev_PrevProcs[r] + _DofLocLevProcFE[Level][_iproc][r]; DofObj_lev++) {

            int dof_pos, irow;
	         if  (r<KK) {  dof_pos = _mesh._Qnode_lev_Qnode_fine[FELevel[r]][ DofObj_lev ];  }
            else if (r==KK) {  dof_pos = DofObj_lev; }
                    irow = _node_dof[Level][ dof_pos + ivar*off_onevar[r] + off_EachFEFromStart[r] ]; 

	    int len[QL];   for (int c=0;c<QL;c++) len[c] = 0;
	    for (int c=0;c<QL;c++) len[c] = (length_row[r][c][ DofObj_lev+1 ]-length_row[r][c][ DofObj_lev ]);
            int rowsize = 0; 
	    for (int c=0;c<QL;c++) rowsize +=_nvars[c]*len[c];
	      graph[irow].resize(rowsize + 1);  //There is a +1 because in the last position you memorize the number of offset dofs in that row

// // // #ifdef FEMUS_HAVE_LASPACK
// // // 
// // //             for (uint jvar=0; jvar<_nvars[QQ]; jvar++) {
// // //             // quadratic-quadratic
// // //                 for (int j=0; j<len[QQ]; j++) {
// // //                     graph[irow][j+jvar*len[QQ]] = _node_dof[Level][ _mesh._node_map[FELevel[QQ]][pos_row[QQ][QQ][j+length_row[QQ][QQ][DofObj_lev]]]+jvar*_mesh._NoNodes[_NoLevels-1]];
// // //                 }
// // // 	    }
// // //                 // quadratic-linear 
// // //                 for (uint jvar=0; jvar<_nvars[LL]; jvar++) {
// // //                     for (int j=0; j<len[LL]; j++) {
// // //                         graph[irow][j+jvar*len[LL]+_nvars[QQ]*len[QQ]] = _node_dof[Level][_mesh._node_map[FELevel[LL]][pos_row[QQ][LL][j+length_row[QQ][LL][DofObj_lev]]]+(jvar+_nvars[QQ])*_mesh._NoNodes[_NoLevels-1]];
// // //                     }
// // //                 }
// // // 
// // //             
// // //             for (uint jvar=0; jvar<_nvars[QQ]; jvar++) {
// // //                 for (int j=0; j<len[QQ]; j++) {
// // //                     graph[irow][j+jvar*len[QQ]] = _node_dof[Level][_mesh._node_map[FELevel[QQ]][pos_row[LL][QQ][j+length_row[LL][QQ][DofObj_lev]]]+jvar*offset];
// // //                 }
// // //             }
// // //             
// // //             for (uint jvar=0; jvar<_nvars[LL]; jvar++) {
// // //                 for (int j=0; j<len[LL]; j++) {
// // //                     graph[irow][j+jvar*len[LL]+_nvars[QQ]*len[QQ]] = _node_dof[Level][ _mesh._node_map[FELevel[LL]][pos_row[LL][LL][j+length_row[LL][LL][DofObj_lev]]]+(jvar+_nvars[QQ])*offset];
// // //                 }
// // //             }
// // //            
// // #endif
 
            int lenoff[QL];  for (int c=0;c<QL;c++) lenoff[c] = 0;
	    for (int c=0;c<QL;c++)  lenoff[c] = length_offrow[r][c][DofObj_lev+1] - length_offrow[r][c][ DofObj_lev ];
	    int lenoff_size = 0;
	    for (int c=0;c<QL;c++) lenoff_size += _nvars[c]*lenoff[c];
            graph[irow][rowsize] = lenoff_size;      // last stored value is the number of in-matrix nonzero off-diagonal values

        }
    }
  } //end r
     
//===========================
        std::cout << " Matrix \n";
	graph.print();
//===========================

    _A[Level]->update_sparsity_pattern_old(graph);  //TODO see how it works

    //  clean ===============
    graph.clear();

    for (int r=0;r<QL;r++) {  //row FE type

        for (int c=0;c<QL;c++) {   //col

            delete []  rowcln[r][c];
            delete []  length_row[r][c];
            delete []  length_offrow[r][c];
            delete []  pos_row[r][c];

        }

        delete []  rowcln[r];
        delete []  length_row[r];
        delete []  length_offrow[r];
        delete []  pos_row[r];

    }

    delete  []  rowcln;
    delete  []  length_row;
    delete  []  length_offrow;
    delete  []  pos_row;

    
    } //end levels   

#ifdef DEFAULT_PRINT_INFO
    std::cout << " ReadMatrix: matrix reading "  << std::endl;
#endif

    return;
}



//=============================
//This function depends on _iproc
void EqnBase::ReadProl(const std::string& name) {

    _Prl.resize(_NoLevels);  //TODO one place is left empty in practice, we can optimize this!!!

    for (uint Level = 1; Level< _NoLevels; Level++) {
  
    uint Lev_c = Level-1;
    uint Lev_f = Level;    

        int FEXLevel_c[QL];
        FEXLevel_c[QQ] = Level-1;                                 //COARSE Level for QUADRATIC
        FEXLevel_c[LL] = (Level-1+_NoLevels)%(_NoLevels+1);
        FEXLevel_c[KK] = Level-1;                                 //COARSE Level for CONSTANT //TODO is this used?
        int FEXLevel_f[QL];
        FEXLevel_f[QQ] = Level;                                  //FINE Level for QUADRATIC
        FEXLevel_f[LL] = Level-1;                                // AAA look at the symmetry, this is exactly (_n_levels + Level1 + 1)%(_n_levels + 1); ! //FINE Level for LINEAR:   Level1=0 means coarse linear, a finer linear is the first coarse quadratic, and so on and so on
        FEXLevel_f[KK] = Level;                                  //FINE Level for CONSTANT //TODO is this used?

    std::ostringstream groupname_lev; groupname_lev <<  "LEVEL" << Lev_c << "_" << Lev_f;
  
    int** rowcln;
    int** len;
    int** lenoff;
    int** Prol_pos;
    double** Prol_val;

    rowcln     = new int*[QL];
    len        = new int*[QL];
    lenoff     = new int*[QL];
    Prol_pos   = new int*[QL];
    Prol_val   = new double*[QL];

    hid_t  file = H5Fopen(name.c_str(),H5F_ACC_RDWR, H5P_DEFAULT); //TODO do I need to open it here for every level?!?

    for (int fe=0; fe<QL; fe++) {

        std::ostringstream fe_family;
        fe_family <<  "_F" << fe;

        //==== DIM ========
        std::ostringstream name0;
        name0 << groupname_lev.str() << "/" << "DIM" << fe_family.str();
        rowcln[fe]= new int[2];
        hid_t dataset = H5Dopen(file, name0.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,rowcln[fe]);
        int nrow=rowcln[fe][0];
        //===== LEN =======
        std::ostringstream name2;
        name2 << groupname_lev.str() << "/" << "LEN" << fe_family.str();
        len[fe] = new int[nrow+1];
        dataset = H5Dopen(file,name2.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,len[fe]);
        //==== OFFLEN ========
        std::ostringstream name3;
        name3 << groupname_lev.str() << "/" << "OFFLEN" << fe_family.str();
        lenoff[fe] = new int[nrow+1];
        dataset = H5Dopen(file,name3.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,lenoff[fe]);
        //===== POS =======
        std::ostringstream name1;
        name1 << groupname_lev.str() << "/" << "POS" << fe_family.str();
        uint count3 = len[fe][nrow];
        Prol_pos[fe] = new int[count3];
        dataset = H5Dopen(file,name1.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,Prol_pos[fe]);
        //===== VAL =======
        std::ostringstream name1b;
        name1b << groupname_lev.str() << "/" << "VAL" << fe_family.str();
        Prol_val[fe] = new double[count3];
        dataset = H5Dopen(file,name1b.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,Prol_val[fe]);

    }
   
    H5Fclose(file);

//=================== end reading =====
    
//======= From here on the EQUATION comes into play, because we have to take into account
// the number of variables of every FE type
//Level goes from 1 to NoLevels-1
   
    uint off_proc = _iproc*_NoLevels;

    _Prl[ Lev_f ] = SparseMatrix::build().release();
// // //     _Prl[ Lev_f ]->init(0,0,0,0); //TODO BACK TO A REASONABLE INIT

    // local matrix dimension
    uint ml[QL]; uint nl[QL];
     for (int fe=0; fe<QL; fe++) { 
       if (fe < KK) {
       ml[fe] = _mesh._off_nd[fe][off_proc + Lev_f +1] - _mesh._off_nd[fe][off_proc];    //  local quadratic    //COARSE (rows)
       nl[fe] = _mesh._off_nd[fe][off_proc + Lev_c +1] - _mesh._off_nd[fe][off_proc];    // global quadratic  //FINE QUADRATIC (cols)
       }
       else if (fe == KK) { 
       ml[fe] = _mesh._off_el[VV][off_proc + Lev_f +1] - _mesh._off_el[VV][off_proc + Lev_f];
       nl[fe] = _mesh._off_el[VV][off_proc + Lev_c +1] - _mesh._off_el[VV][off_proc + Lev_c];
      }
     }
    // pattern dimension
        int nrowt=0;int nclnt=0;
        for (int fe=0;fe<QL;fe++) {
	  nrowt += _nvars[fe]*rowcln[fe][0];
          nclnt += _nvars[fe]*rowcln[fe][1];
	}
    
    Graph pattern;
    pattern.resize(nrowt);
    pattern._m=nrowt;
    pattern._n=nclnt;
    pattern._ml = 0;
    pattern._nl = 0;
     for (int fe=0;fe<QL;fe++) { 
        pattern._ml += _nvars[fe]*ml[fe];  //  local _m
        pattern._nl += _nvars[fe]*nl[fe];  //  local _n
     }
    uint ml_start = _node_dof[Level][_mesh._off_nd[QQ][off_proc]];
    pattern._ml_start = ml_start;

    uint ml_init[QL]; //up to the current processor
      for (int fe=0;fe<QL;fe++) { 
           ml_init[fe]=0;
        for (uint isubd=0;isubd<_iproc; isubd++) {
       if (fe < KK)       ml_init[fe] += _mesh._off_nd[fe][isubd*_NoLevels + Lev_f +1] - _mesh._off_nd[fe][isubd*_NoLevels];
       else if (fe == KK) ml_init[fe] += _mesh._off_el[VV][isubd*_NoLevels + Lev_f +1] - _mesh._off_el[VV][isubd*_NoLevels + Lev_f];
	}
      }
  
    //============= POSITION =========
   for (int fe=0;fe<QL;fe++) {

     for (uint ivar=0;ivar<_nvars[fe];ivar++) {
        for (unsigned int i = ml_init[fe]; i < ml_init[fe]+ml[fe]; i++) {
          int dof_pos_f;
          if (fe < KK)         dof_pos_f = _mesh._Qnode_lev_Qnode_fine[ FEXLevel_f[fe] ][ i ];  //end fe < ql
          else if (fe == KK)   dof_pos_f = i;
          
            int irow  = _node_dof[ Lev_f ][ dof_pos_f + ivar*_DofNumLevFE[ Lev_f ][fe] + _DofOffLevFE[ Lev_f ][fe] ];

	    uint ncol =    len[fe][i+1] -    len[fe][i];
            uint noff = lenoff[fe][i+1] - lenoff[fe][i];
            pattern[irow].resize(ncol+1);
            pattern[irow][ncol] = noff;
// #ifdef FEMUS_HAVE_LASPACK
            for (uint j=0; j<ncol; j++) {
	      int dof_pos_lev_c = Prol_pos[fe][j+len[fe][i]];
	      int dof_pos_c;
	      if      (fe  < KK) dof_pos_c = _mesh._Qnode_lev_Qnode_fine[ FEXLevel_c[fe] ][ dof_pos_lev_c ];
              else if (fe == KK) dof_pos_c = dof_pos_lev_c; 
	      
	      pattern[irow][j] = _node_dof[ Lev_c ][ dof_pos_c + ivar*_DofNumLevFE[ Lev_c ][fe] + _DofOffLevFE[ Lev_c ][fe] ];
            }
// #endif

        }
	    }//end fe < KK
   } //end fe
   

    std::cout << "Printing Prolongator ===========" << std::endl;
    pattern.print();
    _Prl[ Lev_f ]->update_sparsity_pattern_old(pattern);  //TODO shall we make the two operations of updating sparsity pattern and setting values together?

//=========== VALUES ===================
    DenseMatrix *valmat;
    std::vector<uint> tmp(1);
   for (int fe=0; fe<QL; fe++) { 
       for (uint ivar=0;ivar<_nvars[fe];ivar++) {
          for (unsigned int i=ml_init[fe]; i<ml_init[fe] + ml[fe]; i++) {
          int dof_pos_f;
          if (fe < KK)        dof_pos_f = _mesh._Qnode_lev_Qnode_fine[ FEXLevel_f[fe] ][ i ];  //end fe < ql
          else if (fe == KK)  dof_pos_f = i;

          int irow  = _node_dof[ Lev_f ][ dof_pos_f + ivar*_DofNumLevFE[ Lev_f ][fe] + _DofOffLevFE[ Lev_f ][fe] ];

            uint ncol = len[fe][i+1]-len[fe][i];
            tmp[0] = irow;
            std::vector< uint> ind(pattern[irow].size()-1);
            for (uint j=0; j<ind.size(); j++) ind[j] = pattern[irow][j];
            valmat = new DenseMatrix(1,ncol);
            for (uint j=0; j<ncol; j++)(*valmat)(0,j) = Prol_val[fe][j+len[fe][i]];
            _Prl[ Lev_f ]->add_matrix(*valmat,tmp,ind);
            delete  valmat;
        }
     }
   }  //end fe


    for (int fe=0;fe<QL;fe++) {
        delete [] Prol_val[fe];
        delete [] Prol_pos[fe];
        delete [] len[fe];
        delete [] lenoff[fe];
    }

    delete [] Prol_val;
    delete [] Prol_pos;
    delete [] len;
    delete [] lenoff;

    pattern.clear();

    _Prl[  Lev_f ]->close();  //TODO do we need this?
//     if (_iproc==0) _Prl[  Lev_f ]->print_personal();
//     _Prl[  Lev_f ]->print_graphic(false); //TODO should pass this true or false as a parameter
   } //end levels
    
#ifdef DEFAULT_PRINT_INFO
    std::cout << " ReadProl(B): read Op " << name.c_str() << std::endl;
#endif

    return;
}



// =================================================================
//This function depends on _iproc
//TODO AAA This operator here DEPENDS on the BOUNDARY CONDITIONS
// while the prolongator does not depend on the bc's... CHECK THAT
//Also notice that the Matrix, Restriction and Prolongation sparsity patterns
//are assembled by looping over NODES.
// Then, the values are set with NODE LOOPS for Rest and Prol,
// and with ELEMENT LOOP for the MATRIX (the Assemble Function)
//Level goes from 0 to < _NoLevels - 1 ==> Level is COARSE here
  
//    uint Lev_c = Level;
//   uint Lev_f = Level+1;
	//with these you explore arrays that go from 0  to _NoLevels - 1
                           //so where the distinction between QQ and LL is already made
                           // with the EXTENDED levels you explore things that have an additional level,
                           // and so can work both with QQ and with LL
    //the point is: there are parts where you cannot use extended levels, and parts where you can
    //for instance, in this routine the FINE LEVELS and the FINE EXTENDED LEVELS will both be ok,
    //so we can use them in both cases, but we cannot say the same for the COARSE levels and ext levels
    
    //AAA fai molta attenzione: per esplorare la node_dof devi usare Lev_c e Lev_f,
    //perche' sono legati ai DOF (devi pensare che la questione del mesh e' gia' risolta)
void EqnBase::ReadRest(const std::string& name) {
 
  _Rst.resize(_NoLevels);  //TODO why do it bigger?
  
  for (uint Level = 0; Level< _NoLevels - 1; Level++) {
    
    uint Lev_c = Level;
    uint Lev_f = Level+1;
    
    std::ostringstream groupname_lev; groupname_lev <<  "LEVEL" << Lev_f << "_" << Lev_c;
   
    int** rowcln;
    int** len;
    int** lenoff;
    int** Rest_pos;
    double** Rest_val;

    rowcln   = new int*[QL];
    len      = new int*[QL];
    lenoff   = new int*[QL];
    Rest_pos = new int*[QL];
    Rest_val = new double*[QL];

        hid_t  file = H5Fopen(name.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

	for (int fe=0;fe<QL;fe++) {

        std::ostringstream fe_family;
        fe_family <<  "_F" << fe;

        //==== DIM ========
        std::ostringstream name0;
        name0 << groupname_lev.str() << "/" << "DIM" << fe_family.str();
        rowcln[fe]= new int[2];
        hid_t dataset = H5Dopen(file, name0.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,rowcln[fe]);
        int nrow=rowcln[fe][0];
        //===== LEN =======
        std::ostringstream name2;
        name2 << groupname_lev.str() << "/" << "LEN" << fe_family.str();
        len[fe] = new int[nrow+1];
        dataset = H5Dopen(file,name2.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,len[fe]);
        //==== OFFLEN ========
        std::ostringstream name3;
        name3 << groupname_lev.str() << "/" << "OFFLEN" << fe_family.str();
        lenoff[fe] = new int[nrow+1];
        dataset = H5Dopen(file,name3.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,lenoff[fe]);
        //===== POS =======
        std::ostringstream name1;
        name1 << groupname_lev.str() << "/" << "POS" << fe_family.str();
        uint count3 = len[fe][nrow];
        Rest_pos[fe] = new int[count3];
        dataset = H5Dopen(file,name1.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,Rest_pos[fe]);
        //===== VAL =======
        std::ostringstream name1b;
        name1b << groupname_lev.str() << "/" << "VAL" << fe_family.str();
        Rest_val[fe] = new double[count3];
        dataset = H5Dopen(file,name1b.str().c_str(), H5P_DEFAULT);
        H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,Rest_val[fe]);

    }
    
        H5Fclose(file);

//======= From here on the EQUATION comes into play, because we have to take into account
// the number of variables of every FE type
//TODO if you watch carefully, you see that the COARSE levels are EXTENDED,
// while the FINE LEVELS never risk to be "extended", the maximum for them is (_NoLevels-1) !
    
        int FEXLevel_c[QL];
        FEXLevel_c[QQ] = Level;                                 //COARSE Level for QUADRATIC
        FEXLevel_c[LL] = (_NoLevels + Level)%(_NoLevels + 1);   //COARSE Level for LINEAR: Level1=0 means coarse linear, a finer linear is the first coarse quadratic, and so on and so on
        FEXLevel_c[KK] = Level;                                 //COARSE Level for CONSTANT
        int FEXLevel_f[QL];
        FEXLevel_f[QQ] = Level+1;                                  //FINE Level for QUADRATIC
        FEXLevel_f[LL] = Level;                                // AAA look at the symmetry, this is exactly (_n_levels + Level1 + 1)%(_n_levels + 1); ! //FINE Level for LINEAR:   Level1=0 means coarse linear, a finer linear is the first coarse quadratic, and so on and so on
        FEXLevel_f[KK] = Level+1;                                  //FINE Level for CONSTANT

    uint off_proc=_NoLevels*_iproc;

    _Rst[Lev_c] = SparseMatrix::build().release();
// // //     _Rst[Lev_c]->init(0,0,0,0);   //TODO BACK TO A REASONABLE INIT  //we have to do this before appropriately!!!

    int nrowt=0;int nclnt=0;
        for (int fe=0;fe<QL;fe++) {
	  nrowt += _nvars[fe]*rowcln[fe][0];
          nclnt += _nvars[fe]*rowcln[fe][1];
	}


    Graph pattern;
    pattern.resize(nrowt);
    pattern._m  = nrowt;
    pattern._n  = nclnt;        // global dim _m x _n
    pattern._ml = 0;            //  local _m
    pattern._nl = 0;            //  local _n
     for (int fe=0;fe<QL;fe++) { 
        pattern._ml += _nvars[fe]*_DofLocLevProcFE[Lev_c][_iproc][fe]; 
        pattern._nl += _nvars[fe]*_DofLocLevProcFE[Lev_f][_iproc][fe];
     } 
     
    // starting indices for local matrix
   uint ml_start = _node_dof[Lev_c][_mesh._off_nd[QQ][off_proc]];     //  offset proc nodes /*TODO is this correct?!? MAybe this is what is giving me trouble*/
   pattern._ml_start = ml_start;
   uint DofObjInit_lev_PrevProcs_c[QL];
        for (int fe=0;fe<QL;fe++) { DofObjInit_lev_PrevProcs_c[fe] = 0;  }
        for (int fe=0;fe<QL;fe++) { 
    for (uint isubd=0;isubd<_iproc; isubd++) { //up to the current processor
       if (fe < KK)       /*mlinit*/ DofObjInit_lev_PrevProcs_c[fe] += _mesh._off_nd[fe][ isubd*_NoLevels + Lev_c +1 ]  - _mesh._off_nd[fe][isubd*_NoLevels];  
       else if (fe == KK) /*mlinit*/ DofObjInit_lev_PrevProcs_c[fe] += _mesh._off_el[VV][ isubd*_NoLevels + Lev_c +1 ]  - _mesh._off_el[VV][isubd*_NoLevels + Lev_c];   
         }
     }

    //============= POSITION =========
        for (int fe=0; fe<QL; fe++) { 
    for (uint ivar=0;ivar<_nvars[fe];ivar++) {
        for (unsigned int i = DofObjInit_lev_PrevProcs_c[fe]; i< DofObjInit_lev_PrevProcs_c[fe] + _DofLocLevProcFE[Lev_c][_iproc][fe]; i++) {
	  int dof_pos_c;
          if (fe < KK)         dof_pos_c = _mesh._Qnode_lev_Qnode_fine[ FEXLevel_c[fe] ][ i ];
          else if (fe == KK)   dof_pos_c = i;
	  
            int irow  = _node_dof[ Lev_c ][ dof_pos_c + ivar*_DofNumLevFE[ Lev_c ][fe] + _DofOffLevFE[ Lev_c ][fe] ];

	    uint ncol = len[fe][i+1]-len[fe][i];
            uint noff = lenoff[fe][i+1] - lenoff[fe][i];
            pattern[irow].resize(ncol+1);  //when you do resize, it puts a zero in all positions
	    pattern[irow][ncol] = noff;
// pattern structure (was for laspack only)
            for (uint j=0; j<ncol; j++) {
	      int dof_pos_lev_f = Rest_pos[fe][ j+len[fe][i] ];
	      int dof_pos_f;
    
	      if      (fe  < KK)  dof_pos_f = _mesh._Qnode_lev_Qnode_fine[ FEXLevel_f[fe] ][ dof_pos_lev_f ];
              else if (fe == KK)  dof_pos_f = dof_pos_lev_f;
              pattern[irow][j] = _node_dof[ Lev_f ][ dof_pos_f + ivar*_DofNumLevFE[ Lev_f ][fe] + _DofOffLevFE[ Lev_f ][fe] ];
                }

              }
           }
	} //end fe loop


    std::cout << "Printing Restrictor ===========" << std::endl;
    pattern.print();
    _Rst[Lev_c]->update_sparsity_pattern_old(pattern);  //TODO see 
//         _Rst[Lev_c]->close();
//     if (_iproc==0) _Rst[Lev_c]->print_personal(); //there is no print function for rectangular matrices, and print_personal doesnt seem to be working...
// la print stampa il contenuto, ma io voglio solo stampare lo sparsity pattern!
     //Allora cosa faccio: riempio di zeri e poi la stampo! No, di zeri no!!! devi riempirla con qualcos'altro!
//TODO how can I print the sparsity pattern in Petsc BEFORE FILLING the MATRIX ?!?!
    //==============================
//========= SET VALUES =========
//==============================

    DenseMatrix *valmat;
    std::vector<uint> tmp(1);
        for (int fe=0;fe<QL;fe++) {

    for (uint ivar=0;ivar<_nvars[fe];ivar++) {
        for (unsigned int i = DofObjInit_lev_PrevProcs_c[fe]; i< DofObjInit_lev_PrevProcs_c[fe] + _DofLocLevProcFE[Lev_c][_iproc][fe]; i++) {
	  int dof_pos_c;
          if (fe < KK)         dof_pos_c = _mesh._Qnode_lev_Qnode_fine[ FEXLevel_c[fe] ][ i ];
          else if (fe == KK)   dof_pos_c = i;
            int irow     = _node_dof[      Lev_c][ dof_pos_c + ivar*_DofNumLevFE[ Lev_c ][fe] + _DofOffLevFE[ Lev_c ][fe] ];
            int irow_top = _node_dof[_NoLevels-1][ dof_pos_c + ivar*_DofNumLevFE[ _NoLevels-1 ][fe] + _DofOffLevFE[ _NoLevels-1 ][fe] ];
            uint ncol = len[fe][i+1]-len[fe][i];
            tmp[0]=irow;
            std::vector< uint> ind(pattern[irow].size()-1);
// 	    std::cout << "\n ==== " << irow << ": ";
            for (uint i1=0;i1<ind.size();i1++) { ind[i1] = pattern[irow][i1]; /*std::cout << " " << ind[i1] << " ";*/}
            valmat = new DenseMatrix(1,ncol);  //TODO add a matrix row by row...
            for (uint j=0; j<ncol; j++) (*valmat)(0,j) = _bc[irow_top]*Rest_val[fe][ j+len[fe][i] ];
            _Rst[Lev_c]->add_matrix(*valmat,tmp,ind);
            delete  valmat;
            }// end dof loop
         } // end var loop
       } //end fe

    for (int fe=0;fe<QL;fe++) {
        delete [] Rest_val[fe];
        delete [] Rest_pos[fe];
        delete [] len[fe];
        delete [] lenoff[fe];
    }

    delete [] Rest_val;
    delete [] Rest_pos;
    delete [] len;
    delete [] lenoff;

    pattern.clear();

    _Rst[Lev_c]->close();   //TODO Do we really need this?
//     if (_iproc==0)  _Rst[Lev_c]->print_personal(std::cout);
//     _Rst[Lev_c]->print_graphic(false); // TODO should pass this true or false as a parameter

  } //end levels
  
#ifdef DEFAULT_PRINT_INFO
    std::cout << " ReadRest(B): read Op " << name.c_str() << std::endl;
#endif
    return;
}


// ============================================================
/// This function prints the solution: for quad and linear fem
/// the mesh is assumed to be quadratic, so we must print on quadratic nodes
/// first we print QUADRATIC variables over quadratic nodes, straightforward
/// then, LINEAR variables over quadratic nodes --> we interpolate
/// The NODE numbering of our mesh is so that LINEAR NODES COME FIRST
/// So, to interpolate we do the following:
/// First, we put the values of the LINEAR (coarse) MESH in the FIRST POSITIONS of the sol vector
/// Then, we loop over the elements:
///     we pick the element LINEAR values FROM THE FIRST POSITIONS of the SOL VECTOR itself,
///     we multiply them by the PROLONGATOR,
///     we put the result in the sol vector
/// Therefore the sol values are REWRITTEN every time and OVERWRITTEN because of adjacent elements.
// Ok, allora dobbiamo fare in modo che questa routine stampi Quadratici, Lineari e Costanti
// Ricordiamoci che tutto viene fatto in un mesh "RAFFINATO UNA VOLTA IN PIU' RISPETTO AL LIVELLO PIU' FINE"
// per quanto riguarda gli  elementi, noi abbiamo un valore per ogni elemento,
//e poi dovremo trasferirlo ai suoi figli originati dal raffinamento
// quindi dovremo avere una mappa che per ogni elemento ci da' i suoi FIGLI originati da un raffinamento,
// il tutto ovviamente seguendo gli ordinamenti delle connettivita' di FEMuS.
// Allora, per quanto riguarda gli elementi, viene fatta una mesh_conn_lin che si occupa di mostrare gli 
// elementi del mesh linearizzato
//questa connettivita' linearizzata riguarda soltanto il livello FINE
// finora SIA le VARIABILI QUADRATICHE sia le VARIABILI LINEARI 
//sono stampate sul "MESH FINE QUADRATICO, reso come MESH LINEARIZZATO"
// per le variabili COSTANTI, esse le stampero' sulla CELLA, che ovviamente sara' una CELLA LINEARIZZATA,
// la stessa che uso per il resto, il mio file soluzione .xmf credo che possa utilizzare UNA SOLA TIPOLOGIA.
// Quindi dobbiamo appoggiarci alla stessa topologia "raffinata linearizzata".
// Per fare questo penso che possiamo copiare quello che fa il file case.h5!
    
//The printing of the constant will be very similar to the printing of the PID for the case;
// so we could make a common routine which means "print_cell_property_on_linearized_mesh"

//By the way, observe that with multigrid and with the problems about printing in XDMF format with Paraview
//we are in practice having ONE COARSER LEVEL, for LINEAR COARSE NODES,
// and ONE FINER LEVEL for LINEAR NODES, given by this "auxiliary further mesh refinement"

//TODO Parallel: is this function called only by one proc I hope...
//TODO That information should stay INSIDE a FUNCTION, not outside

//TODO this function works only with Level=NoLevels - 1 !
//facciamola funzionare per ogni livello.
//ora che stampo le connettivita' a tutti i livelli per il mesh LINEARIZZATO, posso fare come voglio.
//Quindi stampo su MSH_CONN di un certo livello
//Ok, dobbiamo distinguere come si esplora la node_dof map e la lunghezza di ogni livello
//
    // ===============================================
//An idea could be: can we make the projection of any FE solution onto this "REFINED LINEARIZED" MESH
// an automatic process, by using the Prol operators?

//Now, the behaviour of this function should be somewhat "parallel" to the construction of the node_dof.

//TODO the group must be REOPENED by EACH EQUATION.. forse Gopen anziche' Gcreate
// I need to have a group that is first created, then closed and reopened without being emptied...
// it should do the same as a File does

//ok, so far we will print the variables with _LEVEL... now we have to fix the wrapper as well

//Ok, now that we make it print things for ALL variables at ALL levels, we need to JUMP on the _node_dof map.
// TODO the NODE DOF NUMBERING is "CONSECUTIVE" ONLY ON THE FINE LEVEL... and only for QUADRATIC variables, 	for (uint i=0;i< n_nodes;i++) 
// which is the same order as the mesh!
//So only on the FINE level you can loop over NODES, in all the other cases 
// you need to loop over ELEMENTS
//For all the other levels we have to JUMP... so we will jump the same way as when we build the node dof!
//Now for all variables, loop over all subdomains, collect all the values, and print them...
// I cannot print them all together

      //TODO TODO TODO here it is more complicated... pos_sol_ivar will sum up to the QUADRATIC NODES,
      //but here we pick the quadratic nodes MORE THAN ONCE, so how can we count them EXACTLY?!!
      //the problem is that when you are not on the fine level you have jumps, so, in order to count the LINEAR nodes,
      //you need a flag for every quadratic node that tells you if it is also linear or not...
      //we should either build a flag field in advance, so that we don't do it now, or maybe later is ok
      //since we don't keep any extra info about each node (we do not have a Node class, or an Elem class),
      //now it's time we have some flags...
      //i guess i should loop over all the quadratic nodes, and let the flag for the linear only
      // ok so when i pick the linear nodes i may set the flag is linear there
      //flag = 1 means: it is linear, otherwise 0

// // //       for (int i = _mesh._off_nd[QQ][off_proc];
// // // 	       i < _mesh._off_nd[QQ][off_proc+Level+1];i++) {
// // //       if (i< _mesh._off_nd[QQ][off_proc]
// // // 	   + _mesh._off_nd[LL][off_proc + Level+1 ]
// // // 	   - _mesh._off_nd[LL][off_proc])
// // //         
// // // 	 }
  //now, in the quadratic nodes, what are the positions of the linear nodes?    
  //if I am not wrong, the first nodes in the COARSE LEVEL in every processor
  //are the LINEAR ONES, that's why the offsets are what they are...
  //in fact, we are looping on the GEomElObjects of each level.

//     int* flag_is_linear = new int[n_nodes]; // we are isolating MESH NODES correspoding to LINEAR DOFS
    //what, but we already know from the node numbering what nodes are linear, because we distinguished them, right?!
//     for (uint i=0;i< n_nodes;i++) flag_is_linear[i] = 0;

//LINEAR VARIABLES
//these are re-dimensionalized //so you dont need to multiply below! //the offset for _node_dof is always quadratic, clearly, because the mesh is quadratic
//        //I am filling two QUADRATIC arrays first by the LINEAR POSITIONS
	  //so pay attention that if you are not setting to ZERO for the linear case, but exactly replacing at the required points

//Always remember that in order to pick the DofValues you have to provide the DofObjects in the correct manner
	  
// QUADRATIC VARIABLES  
// pos_in_mesh_obj gives me the position of the GEomElObject: in fact the quadratic dofs are built on the quadratic GeomEls exactly in this order
//here we are picking the NODES per subd and level, so we are sure don't pass MORE TIMES on the SAME NODE

// PRoblem with the linear variables in PrintVector and PrintBc. 
//There is one line that brings to mistake, but TWO different mistakes.
// in PrintBc it seems to be related to HDF5;
// in PrintVector it seems to concern PETSC
//so that is the wrong line, if i comment it everything seems to work for any processor.
//with two and three processors it seems to give even different errors...
//when there is an error related to HDF5, the stack has the _start thing...
// with two procs there is an error related to Petsc,
// with three procs there is an error related to HDF5...

//If I only use PrintBc and not PrintVector, 
//both with 2 and 3 processors the errors are related to HDF5... this is so absolutely weird...

//Now it seems like I am restricted to that line. That line is responsible for the error.
//Then, what is wrong with that? I would say: I am doing sthg wrong on the array positions.
//Let me put a control on the position. maybe the _el_map is somehow ruined

//Quando hai un segmentation fault, devi concentrarti non tanto sui VALORI ASSEGNATI
// quanto sugli INDICI DEGLI ARRAY !!!

        // to do the interpolation over the fine mesh you loop over the ELEMENTS
        //So of course you will pass through most nodes MORE THAN ONCE
        //after setting correctly the linear nodes, then the quadratic ones are done by projection, without any problem
        //for every element, I take the linear nodes.
        //Then I loop over the quadratic nodes, and the value at each of them is the prolongation
        // of the linear values.
        //So, for every quadratic node, I compute the value and set the value 
        //now, the connectivity is that of the quadratic mesh, but with respect
        //to the FINE node numbering
        //now we have to convert from the fine node numbering to the node numbering at Level!!!
        //do we already have something to do this in the MESH?!?
        //given a quadratic node in FINE NUMBERING, obtained by an element AT LEVEL L,
        //can we obtain the position of that node AT LEVEL L
        //don't we have the connectivities AT LEVEL L,
        //where the numbering goes from ZERO to n_nodes_level?
        //I would say it is exactly the _node_dof FOR THE FIRST QUADRATIC VARIABLE!
        //TODO OF COURSE THAT WOULD IMPLY THAT AT LEAST ONE QUADRATIC VARIABLE is BUILT
        //TODO Also, I have to check that COARSER GEOMETRIES are ASSOCIATED to COARSER TOPOLOGIES!
        // _node_dof in serial is the IDENTITY, BUT NOT IN PARALLEL!!
        //Per passare dal Qnode in FINE NUMBERING al Qnode in SERIAL NUMBERING
        //non basta passarlo alla _node_dof del livello, perche' quando sei in parallelo
        //lui conta prima i dof quadratici, poi quelli lineari, poi quelli costanti, e bla bla bla...
        //e quindi bisogna cambiare il modo di pigliarli!
        //Si' bisogna usare qualcosa di SIMILE alla NODE_DOF, ma NON la node_dof,
        //perche' quella, per dato livello, conta i dof QQ,LL,KK del proc0, poi QQ,LL,KK del proc1, and so on...
        //OK OK OK! Direi che quello che dobbiamo usare e' proprio la node_map, quella che leggiamo dal mesh.h5!!!
        //Quindi _node_dof e' per TUTTI i DOF,
        // mentre _node_map e' solo per i NODI GEOMETRICI!
        //mi sa pero' che la _node_map fa l'opposto di quello che vogliamo noi, cioe'
        //dato il nodo di un certo livello ti restituisce il nodo FINE
        //Noi invece abbiamo il NODO FINE, e vogliamo avere il nodo di un certo livello.
        //Questo si otterrebbe leggendo TUTTA la map e non comprimendola come facciamo...
        //oppure io la ricostruisco rifacendo il loop della node_dof ma solo per UNA variabile quadratica!
        //praticamente, MA NON DEL TUTTO!, la _node_map e' l'inverso della _node_dof!!!
        //TODO ma non c'e' nessun altro posto nel codice in cui devi passare 
        //da QNODI FINI a QNODI di LIVELLO?
        // sembra di no, piu' che altro devi passare da QNODI FINI a DOF di LIVELLO!
        //ora qui, essendo che dobbiamo stampare su griglie di diverso livello,
        //e siccome le CONNETTIVITA' di TUTTI I LIVELLI SONO DATE RISPETTO all'unico FINE NUMBERING
        // (TODO beh, volendo potrei utilizzare le connettivita' "di livello" che stampo nel file msh_conn_lin...
        //Il fatto e' che quelle non sono variabile di classe (in Mesh) e quindi mi limito a calcolare, stampare e distruggere...)
        // ALLORA DEVO FARE IL PASSAGGIO da QNODI FINI a QNODI DI LIVELLO.
        //Questo me lo costruisco io domattina.
        //Siccome e' gia' stato calcolato nel gencase, allora mi conviene evitare di fare questo calcolo.
        //Quando leggo il vettore dal gencase, posso costruire un vettore di PAIRS...
        //ovviamente e' piu' lento, perche' deve cercare un elemento con degli if...
        //allora faccio un array 2x, per cui calcolo io l'indice da estrarre...

// We have to find the LINEAR NODES positions in the QUADRATIC LIST
// the QUADRATIC list at each level is based on 
// the "LINEARIZED" CONNECTIVITIES of that level
//    plus the LIST OF COORDINATES of that level
	    
//Now there is another mistake, in printing the LINEAR CONNECTIVITIES!
//The fact is always that _el_map gives the Qnodes in FINE NUMBERING,
//while I want the Qnode numbering AT EACH LEVEL.
//because at coarser levels I have fewer nodes, so the list of coords is shorter!

//Ok, the printing of the linearized connectivities is WRONG at NON-FINE LEVEL
        
//TODO AAAAAAAAAAAAAAAAA: Well, the thing is this: the connectivities of any level must be expressed 
// in terms of the node numbering OF THAT LEVEL!
//if the connectivities are expressed with the FINE NODE NUMBERING, 
//then we always have to convert to the LEVEL NODE NUMBERING!!!
//I guess the COORDINATES are WRONG...
//Ok, first of all we already have the connectivities at all levels
// of quadratic elements in mesh.h5, in FINE NODE NUMBERING.
//We only have to convert them to LEVEL NODE NUMBERING.

//Ok, I fixed the COORDINATES of the Qnodes at EVERY LEVEL.
//Plus, I have the "LINEARIZED" CONNECTIVITIES at EVERY LEVEL.
//Therefore, it seems like I can print any vector at EVERY LEVEL.
//Now, I still have a problem for the LINEAR variables.
//I guess I'm putting them in the WRONG PLACES before the interpolations.

//allora, quando faccio i loop con _off_nd, quei numeri li' corrispondono alle posizioni nel vettore sol?
//TODO chi da' l'ordine dei nodi del mesh A CIASCUN LIVELLO?
// Le COORDINATE DEL MESH A CIASCUN LIVELLO, quello e' l'ordine di riferimento!

//ok, the position i corresponds to the FINE MESH... then you must translate it for the LEVEL mesh

//For the linear variables we have TWO PROC LOOPS

         //Now I guess I have to pick the position of the linear nodes on the mesh by using the ExtendedLevel on the map...
         //In the following interpolation the loop is on the elements. From the elements you pick the connectivities,
//          which yield you the FINE QUADRATIC NODES, from which you pick the node numbers at level


//  PRINT OF SOL0000
// the first sol, and the case, are printed BEFORE the INITIAL CONDITIONS are set,
// or the initial conditions are set only at the FINE LEVEL?
// I would say the first one, otherwise I would expect different solutions...
// The first sol and the case are supposed to be equal...
//No wait, what i am saying is not true, because the first sol and the case 
// at the FINE LEVEL are printed as they should,
// after the initial conditions,
//but not on COARSER LEVELS!
// Ok now the initial conditions are set only AT THE FINE LEVEL.
// Depending on the kind of multigrid cycle, you may want to set the initial conditions
// ONLY AT THE FINE LEVEL or AT ALL LEVELS...
// For us, let us just do the GenIc function in such a way that all the levels can be treated separately.
// then, if we need it, we call it for ALL LEVELS, or we call it for ONLY THE FINE, or ONLY THE COARSE, or whatever...

//TODO ok, we have to remember one basic principle about our multigrid algorithm:
//the true solution is contained only at the FINE LEVEL!
//at all the other levels we have the DELTAS
//so, the level where to pick the values is the FINE LEVEL.
//Then, of course, we will print at all levels, but the VALUES from x_old are taken from the FINE LEVEL
//no... but wait a second... i am printing at all levels, so that's fine! I wanna print the RESIDUAL for all levels,
//except for the fine level where i print the true solution

// This prints All Variables of One Equation    
void EqnBase::PrintVector(std::string namefile) {

    hid_t file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
   
    // ==========================================
    // =========== FOR ALL LEVELS ===============
    // ==========================================
    for (uint Level = 0; Level < _NoLevels; Level++)  {

      std::ostringstream grname; grname << "LEVEL" << Level;
//      hid_t group_id = H5Gcreate(file_id, grname.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//      hid_t group_id = H5Gopen(file_id, grname.str().c_str(),H5P_DEFAULT);
    
    int NGeomObjOnWhichToPrint[QL];
    NGeomObjOnWhichToPrint[QQ] = _mesh._NoNodesXLev[Level];
    NGeomObjOnWhichToPrint[LL] = _mesh._NoNodesXLev[Level];
    NGeomObjOnWhichToPrint[KK] = _mesh._NoElements[VV][Level]*_mesh._GeomEl.n_se[VV];
    
    const uint n_nodes_lev = _mesh._NoNodesXLev[Level];
    double* sol_on_Qnodes  = new double[n_nodes_lev];  //TODO VALGRIND //this is QUADRATIC because it has to hold  either quadratic or linear variables and print them on a QUADRATIC mesh
    
    // ===================================
    // ========= QUADRATIC ===============
    // ===================================
    for (uint ivar=0; ivar<_nvars[QQ]; ivar++)        {
      
      int pos_in_mesh_obj = 0;   
         for (uint isubdom=0; isubdom<_mesh._NoSubdom; isubdom++) {
            uint off_proc=isubdom*_NoLevels;
     
            for (int fine_node = _mesh._off_nd[QQ][off_proc];
                     fine_node < _mesh._off_nd[QQ][off_proc+Level+1]; fine_node++) {
      
  	int pos_in_sol_vec_lev = _node_dof[Level][fine_node + ivar*_DofNumLevFE[ Level ][QQ] + _DofOffLevFE[ Level ][QQ] ];
	int pos_on_Qnodes_lev = _mesh._Qnode_fine_Qnode_lev[Level][ fine_node ]; 

#ifndef NDEBUG
         if ( pos_on_Qnodes_lev >= (int) n_nodes_lev ) { std::cout << "^^^^^^^OUT OF THE ARRAY ^^^^^^" << std::endl; abort(); }
#endif
        sol_on_Qnodes[ pos_on_Qnodes_lev/* pos_in_mesh_obj*/ ] = (*_x_old[Level])(pos_in_sol_vec_lev) * _refvalue[ ivar + _VarOff[QQ] ];
	pos_in_mesh_obj++;
	  }
       }  //end subd
       
#ifndef NDEBUG
	 if (pos_in_mesh_obj != NGeomObjOnWhichToPrint[QQ]) { std::cout << "Wrong counting of quadratic nodes" << std::endl; abort(); }
#endif

     std::ostringstream var_name;  var_name << _var_names[ ivar + _VarOff[QQ] ] << "_" << grname.str(); 	 //         std::string var_name = grname.str() + "/" + _var_names[ivar];
     hsize_t  dimsf[2];  dimsf[0] = NGeomObjOnWhichToPrint[QQ];  dimsf[1] = 1;
     IO::print_Dhdf5(file_id,var_name.str(),dimsf,sol_on_Qnodes);   //TODO VALGRIND

     }

    // =================================
    // ========= LINEAR ================
    // =================================
    uint elnds[QL_NODES];
    elnds[QQ] = _mesh._GeomEl._elnds[VV][QQ];
    elnds[LL] = _mesh._GeomEl._elnds[VV][LL];
    double* elsol_c = new double[elnds[LL]];
    
    for (uint ivar=0; ivar<_nvars[LL]; ivar++)        {
      
//               for (uint i=0; i< n_nodes_lev; i++) { sol_on_Qnodes[i] = 0.; }
    
	for (uint isubdom=0; isubdom<_mesh._NoSubdom; isubdom++) {
	     uint off_proc=isubdom*_NoLevels;
            for (int fine_node = _mesh._off_nd[QQ][off_proc];
                     fine_node < _mesh._off_nd[QQ][off_proc]+
                         _mesh._off_nd[LL][off_proc + Level+1 ]
                       - _mesh._off_nd[LL][off_proc]; fine_node++) {
	      
	    int pos_in_sol_vec_lev = _node_dof[Level][ fine_node + ivar*_DofNumLevFE[ Level ][LL] + _DofOffLevFE[ Level ][LL] ];
 	    int pos_on_Qnodes_lev = _mesh._Qnode_fine_Qnode_lev[Level][ fine_node ];

#ifndef NDEBUG
	 if ( pos_in_sol_vec_lev == -1 ) { std::cout << "Not correct DOF number at required level" << std::endl; abort(); }
         if ( pos_on_Qnodes_lev >= (int) n_nodes_lev ) { std::cout << "^^^^^^^OUT OF THE ARRAY ^^^^^^" << std::endl; abort(); }
#endif

         sol_on_Qnodes[ pos_on_Qnodes_lev ] = (*_x_old[Level])(pos_in_sol_vec_lev) * _refvalue[ ivar + _VarOff[LL] ];
	 
            }
        }

        //  2bB element interpolation over the fine mesh -----------------------
        // the way you filled linear positions before completely affects what happens next, which is only geometric
        for (uint iproc=0; iproc<_mesh._NoSubdom; iproc++) {
               uint off_proc = iproc*_NoLevels;
	       int iel_b = _mesh._off_el[VV][off_proc + Level];
	       int iel_e = _mesh._off_el[VV][off_proc + Level + 1];
	       
            for (int iel = 0; iel < (iel_e-iel_b); iel++) {
      
                for (uint in=0; in < elnds[LL]; in++) {
		  int pos_Qnode_fine = _mesh._el_map[VV][ (iel+iel_b)*elnds[QQ]+in ];
		  int pos_Qnode_lev  = _mesh._Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];
		  elsol_c[in] = sol_on_Qnodes[ pos_Qnode_lev ];   /**_refvalue[ivar]*/ //Do not multiply here!
		}

                for (uint in=0; in < elnds[QQ]; in++) { //TODO this loop can be done from elnds[LL] instead of from 0
                    double sum=0.;
                    for (uint jn=0; jn<elnds[LL]; jn++) {
                        sum += _AbstractFE[LL]->get_prol(in*elnds[LL]+jn)*elsol_c[jn];
                    }
                    
                    int pos_Qnode_fine = _mesh._el_map[VV][ (iel+iel_b)*elnds[QQ]+in ];       //Qnode in FINE NUMBERING
                    int pos_Qnode_lev  = _mesh._Qnode_fine_Qnode_lev[Level][pos_Qnode_fine];  //Qnode in Level NUMBERING

#ifndef NDEBUG
                    if ( pos_Qnode_lev == -1 ) { std::cout << "Not correct node number at required level" << std::endl; abort(); }
 		    if ( pos_Qnode_lev >= (int) n_nodes_lev ) { std::cout << "^^^^^^^OUT OF THE ARRAY ^^^^^^" << std::endl; abort(); }
#endif

 		    sol_on_Qnodes[ pos_Qnode_lev ] = sum;
                 }
              }
          } // 2bB end interpolation over the fine mesh --------
        
     std::ostringstream var_name; var_name << _var_names[ ivar + _VarOff[LL] ] << "_" << grname.str();
     hsize_t  dimsf[2]; dimsf[0] = NGeomObjOnWhichToPrint[LL];  dimsf[1] = 1;
     IO::print_Dhdf5(file_id,var_name.str(),dimsf,sol_on_Qnodes);
     
    } // ivar linear

      delete []elsol_c;
      delete []sol_on_Qnodes;

     // ===================================
     // ========= CONSTANT ================
     // ===================================
  double *sol_on_cells;   sol_on_cells = new double[ NGeomObjOnWhichToPrint[KK] ];

  for (uint ivar=0; ivar<_nvars[KK]; ivar++)        {
      
  int cel=0;
  for (uint iproc=0; iproc<_mesh._NoSubdom; iproc++) {
               uint off_proc = iproc*_NoLevels;
   
            int sum_elems_prev_sd_at_lev = 0;
	    for (uint pr = 0; pr < iproc; pr++) { sum_elems_prev_sd_at_lev += _mesh._off_el[VV][pr*_NoLevels + Level + 1] - _mesh._off_el[VV][pr*_NoLevels + Level]; }

	    for (int iel = 0;
              iel <    _mesh._off_el[VV][off_proc + Level+1]
                      - _mesh._off_el[VV][off_proc + Level]; iel++) {
             int elem_lev = iel + sum_elems_prev_sd_at_lev;
	  int dof_pos_lev = _node_dof[Level][ elem_lev + ivar*_DofNumLevFE[ Level ][KK] + _DofOffLevFE[ Level ][KK] ];   
      for (uint is=0; is< _mesh._GeomEl.n_se[VV]; is++) {      
	   sol_on_cells[cel*_mesh._GeomEl.n_se[VV] + is] = (*_x_old[Level])(dof_pos_lev) * _refvalue[ ivar + _VarOff[KK] ];
      }
      cel++;
    }
  }
  
  std::ostringstream varname; varname << _var_names[ ivar + _VarOff[KK] ] << "_" << grname.str();         //   std::string varname = grname.str() + "/" + _var_names[_nvars[QQ]+_nvars[LL]+ivar];
  hsize_t dimsf[2]; dimsf[0] = NGeomObjOnWhichToPrint[KK]; dimsf[1] = 1;
  IO::print_Dhdf5(file_id,varname.str(),dimsf,sol_on_cells);   
      
    } //end KK

     delete [] sol_on_cells;
  
//         H5Gclose(group_id);
	
    } //end Level
    
    H5Fclose(file_id);   //TODO VALGRIND

    return;
}


// ===================================================
/// This function reads the system solution from namefile.h5
//TODO this must be modified in order to take into account KK element dofs
void EqnBase::ReadVector(std::string namefile) {
//this is done in parallel

  std::cout << "ReadVector still has to be written for CONSTANT elements, BEWARE!!! ==============================  " << std::endl;
  
    const uint Level = _NoLevels-1;
    
    const uint mesh_ord = (int) _mesh._mesh_rtmap.get("mesh_ord");
    const uint offset   =       _mesh._NoNodesXLev[_NoLevels-1];

    // file to read
    double *sol=new double[offset]; // temporary vector
    hid_t  file_id = H5Fopen(namefile.c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

    // reading loop over system varables
    for (uint ivar=0;ivar< _nvars[LL]+_nvars[QQ]; ivar++) {
        uint el_nds=_mesh._GeomEl._elnds[VV][QQ];
        if (ivar >= _nvars[QQ]) el_nds = _mesh._GeomEl._elnds[VV][LL]; // quad and linear
        // reading ivar param
       std::ostringstream grname; grname << _var_names[ivar] << "_" << "LEVEL" << Level;
        IO::read_Dhdf5(file_id,grname.str(),sol);
        double Irefval = 1./_refvalue[ivar]; // units

        // storing  ivar variables (in parallell)
        for (int iel=0;iel <  _mesh._off_el[0][_iproc*_NoLevels+_NoLevels]
                -_mesh._off_el[0][_iproc*_NoLevels+_NoLevels-1]; iel++) {
            uint elem_gidx=(iel+_mesh._off_el[0][_iproc*_NoLevels+_NoLevels-1])*_mesh._GeomEl._elnds[VV][mesh_ord];
            for (uint i=0; i<el_nds; i++) { // linear and quad
                int k=_mesh._el_map[0][elem_gidx+i];   // the global node
                _x[_NoLevels-1]->set(_node_dof[_NoLevels-1][k+ivar*offset],sol[k]*Irefval); // set the field
            }
        }
    }

    _x[_NoLevels-1]->localize(*_x_old[_NoLevels-1]);
    // clean
    H5Fclose(file_id);
    delete []sol;
    
    return;
}






//here, we could pass a std::VECTOR of the UNKNOWN INTERNAL QUANTITIES,
//not a map because we want things to be ordered.
//this vector must be passed to the EqnNS constructor so that it can be passed
//to the DA. 
//so we must do things OUTSIDE the NS constructor, i.e. in the EquationsMap.
//We cannot do INSIDE the NS constructor because it is called AFTER the DA constructor.
//where Velocity is in the first place and pressure is in the second
//and later one can add OTHER Quantities with
//at this point things are done for only two types,
//quadratic and linear
//in general, suppose that you want linear velocity/linear pressure,
//then you do not use the quadratics
//BUT, we decided to have every shape funcs and stuff
//both for QQ and for LL
//So, we'll leave these like that, 
//so that linear temperature can get velocity values from QUADRATIC velocity
//but, actually, we should do that only the things for the INVOLVED ORDERS
// (Internally INVOLVED and Externally INVOLVED) are first allocated HERE
// and THEN FILLED in the GenMatRhsVB
//this would mean doing a LOOP of "new" and correspondingly a LOOP of "delete"
//for the involved quantities
//clearly, it is important that these "new" and delete things are not done INSIDE
//the Gauss loop but OUTSIDE, a priori
//the problem is that in order to know these we must do a loop 
//not only over the Internal Quantities but also the EXTERNAL ones  
//because they need to be interpolated as well!
//also, once we know which orders are involved, it is not absolutely immediate
// what DERIVATIVE orders need to be actually filled: phi,dphidxezeta, dphidxyz,dphidxyz3D
//While the FE order depends on the Involved QUANTITIES,
//the DERIVATIVE order depends on the Involved OPERATORS.
//Involved Quants or Ops means BOTH Internally AND Externally!
//therefore, we cannot know.
//One Day, we'll do an Equation as a list of Quantities and a list of Operators acting
//on each of them.
//So the rule will be: we FILL all that we ALLOCATE here.
//ok, so we'll do a vector of internal, and a map of external for now...
//in this way we may write in more general form all the things

//===============related to the SINGLE EQUATION
//=========number of variables
//the number of variables of my own equation must be deduced from the INTERNAL QUANTITIES
//now let us consider TWO internal quantities 
//the InternalVector actually belongs to the BASE class,
//but n_vars[] belongs to THIS class
//so, we must consider that we have to turn things into QUADRATIC and LINEAR,
//no matter how the Quantities are.
//Suppose that you have FOUR Quantities, some of them are scalar, some vectors,
//some are linear, some are quadratic, we must turn them into
//nvars[QQ] and nvars[LL].
//Then, when trying to pick the DOFS, we only need to know the offsets
//for passing from the FIRST Quantity (some components and some order)
//to the Second Quantity. Alright.
//so now we''ll just do a TWO by TWO. From TWO Quantities to TWO different FE.
//---> the important thing is not to use nvars_q_in and nvars_l_in any longer!

//---> now, we must understand if the QQ and LL variables used later are
//used in the sense "quadratic or linear - PUBLIC DOMAIN"
//or Quadratic for the FIRST Quantity and LINEAR for the SECOND Quantity
//----so, we have to understand where they are called.
//_fe[0] must be called for being the _fe of quadratic,
//not for being the fe of velocity...
//well, actually i think that all the variables are treated here for being 
//quadratic or linear variables










 
 



//=======================================

// one day this will be   GetElDofs(myvect,iel);
  //the external quantities might go ahead on their own,
//the only thing we need is the Gauss point value
//The thing is simple: I need the value of the external quantity at a gauss point,
//FOR THE ELEMENT I AM IN
//So, both the SEARCH FOR DOFS and the INTERPOLATION AT THE GAUSS POINT
//can be done INDEPENDENTLY BY EACH QUANTITY!
  //for a given element,
   //every EQUATION and EVERY VARIABLE IN THE EQUATION
   //has ITS OWN _node_dof map!!!
//considering how the dof map is filled now, it works only for QUADRATIC and LINEAR FINITE ELEMENTS
//we should REWRITE the _node_dof map FOR EACH FINITE ELEMENT ASSOCIATED TO THE VARIABLES
//for now, we use the connectivity of el_conn which works for BOTH QUADRATIC AND LINEAR.

// TODO THAT THIS FUNCTION is NOT COMPLETE! If you want to pick the PRESSURE from NS,
// you cant do like that!
// you need to know how the dofmap behaves.
// For now, and throughout the code, there is no dofmap class to handle this thing.
//But, the dofmap is always like this: First the quadratic variables, then the linear variables
//also, el_conn is always quadratic and trhis is not a problem, since 
//you only have to loop over the vertices which are the FIRST in the connectivity order.
//so, (idim+ var_off)*dof_off_gl 
//myvect can be either quadratic or linear
//if myvect == QQ, varoff=0
//if myvect == LL, varoff=nvars_q. So this depends on the DA equation in which I find myself.
// This is ok, because I have an equation pointer, so i know where the variables lies on.
// TODO For now, we cannot pick only one component, but we can pick
//either ALL the QUADRATIC, or ALL the LINEAR components
//for now all our equations behave in this manner, but later we'll change this thing
//this is something to consider.
//What happens if I have both velocity and pressure that are linear,
//and i want to pick only velocity or Only pressure?
//i want to pick only a PART of the linear variables
//we need this distinction because the _node_dof map is divided into the quadratic and the linear part
//this is because we would have TWO DISTINCT PHYSICAL TYPES of LINEAR VARIABLES
//so, if both variables are linear, you have that lin_off = 0, because nodedof starts with linear stuff;
//but then, while Velocity is the First and it starts right at the beginning, Pressure
//starts right after velocity.
//So, for every DA equation, instead of distinguishing between 
//the QUADRATIC variables and the LINEAR variables
//We distinguish between a FIRST Quantity and a SECOND Quantity
//(later we will do just a LIST of MORE Ordered Quantities)
//based on them we have the number of variables nvars[FIRST] and nvars[SECOND] and 
//so we can turn them into nvars[QQ] and nvars[LL],
//on which all the loops are based
//So, we'll have that for the FIRST Quantity the offset is ZERO,
//for the second quantity the offset is the number of variables of the First Quantity:
//lin_off = nvars[FIRST]
//all the loops like multigrid operators or stuff are based only on _nvars[0] and _nvars[1],
//so they're not aware if the first is velocity and the second is pressure,
//they just dont care, they work BLINDLY for an ARBITRARY QUANTITY of QUADRATIC and LINEAR variables







//**************************************************
//here, we should explicitly pass the THREE QUANTITIES,
// or better a vector of Quantities;
//the structures of nodedof 
// and bc_eldofs
// are both ordered 
// first by QUADRATIC and 
// then by LINEAR variables.
//But, the Vect we provide are ordered 
// by the way we provide them...
//so they may be first quadratic then linear then quadratic...
//we should remember the passage from the initial order 
//to the QL order...
//Since later you fill the matrix based on the 
//QL order, so you explicitly do what the initvectors
//does at the beginning,
//the point is that you must pay attention to
//the order given at the BEGINNING
//So we have a hybrid reasoning: 
// by List of Quantities 
// and 
// by List of Quadratic and Linear Variables.
//here we are picking things from structures 
//that are all ordered by QL variables
//but some of those are distinguished by Quantity...




//now, we'll consider that qty zero and qty one are ordered
//this must be done because in the el_dof_indices and bc_eldofs the nodes are ordered, first quadratic then linear,
//or first linear then linear...
//The KeM and FeM are filled as ORDERED, FIRST QTYZERO("velocity") then QTYONE("pressure")
//so the bc_eldofs and el_dof_indices are correspondingly ORDERED
//therefore also the corresponding global structures are ORDERED
//if you want to write FIRST the PRESSURE then the VELOCITY equations you must "REPOSITION" the element matrix.
//if you do linear linear how can you be sure that velocity comes first?
//I must assure that also the GLOBAL bc and _node_dof start from the qty_zero (either quadratic or linear) and then qty_one


//I'll do the version for a vector of Vect
//when you explore the global structures, 
//they all have the structure with quadratic and linear
//distinguished, and ordered, first quadratic then linear.
//So we have to explore these with nvars[QQ] and nvars[LL]
//Now, the point is: when we want to explore the xold 
//to retrieve the dof values, we have to do the inverse 
//path that allowed us to 
//the first quadratic goes to the first quantity
//the second quadratic goes to the second quantity
//the third quadratic goes to the third quantity... and so on...
//in that case you have to consider that you may have 
//more variables for a single quantity
//we can say:we loop over the quantities,
// in the order we gave
// if they are quadratic we pick the dofs from the quadratic part, in the order we need;
//if they are linear, we pick the dofs from the linear part.
//then, we must keep in mind the current cursor for both the quadratic and the linear parts
//E' possibile fare in modo che tutto sia costruito sulle quantities?
//beh, in qualche modo, quando riempi i pezzi della matrice di elemento, 
// ti servira' una distinzione per componenti... anche se sicuramente avere la quantita'
//tutta insieme ti puo' permettere di scrivere in modo compatto
//certi operatori specie su strutture vettoriali
//direi che e' il caso di separare quello che e' suddiviso per QL
// da quello che e' suddiviso per Quantity

//finora l'idea e': inizializzo per quantita' fisiche,
//poi converto tutto subito in QL,
//poi costruisco tutto per QL

// ad esempio, le condizioni al contorno sono messe tutte insieme,
// e anche le condizioni iniziali sono tutte insieme.
// Se io volessi un modo per riempire le condizioni al contorno
// non per QL ma per Quantita',
//dovrei fare il loop sulle quantita' nell'ordine che e' stato dato,
//poi fare il loop sul numero di componenti di ciascuna quantita',
//e in base a quello mettere una incognita oppure un'altra

//TODO: Now, we are MULTILEVEL. Nevertheless, the boundary conditions are imposed
//only on the FINE LEVEL. Now, the point is that the DofObj for 
//Now, the boundary conditions are imposed by NODE COORDINATES,
// lopping over ALL THE NODES, also the VOLUME NODES.
// Then, you distinguish the BOUNDARY NODES from the VOlUME nodes 
// by the COORDINATES
// In the same way, we will distinguish the VOLUME ELEMENTS from the 
// BOUNDARY elements by the COORDINATES of the CENTER.
// the point is that we are talking about the VOLUME ELEMENTS,
// so we do not know how big they are for putting ifs on the coordinates
// Then, we will loop over volume elements;
// if a volume element HAS a FACE ON THE BOUNDARY, 
// and if that face satisfies the "LINE CONDITION" on the boundary,
// then SET the corresponding BC.
// The point is that we need the list of faces for each element, 
// and then the coordinates of the nodes of each face.
// so, so far, let me avoid all this stuff and put the boundary conditions
// only on the QQ and LL
// The important thing to remember is that the bc_flag 
// is a flag that is associated to EACH DOF,
// in particular to EACH ROW.

//We need to understand what this routine receives and what it gives.
//This routine receives DOF CARRIERS, in terms of NODE DOF CARRIERS and ELEM DOF CARRIERS.
//What does it give? it gives DOF POSITIONS, and DOF BC FLAGS.
//Of course for any DOF CARRIER there may be ONE OR MORE DOFS,
//so one or more DOF POSITIONS, and one or more DOF BC FLAGS
//So, given a "dof carrier" we do not know how many dofs are associated to it.
//Therefore, an alternative would be: starting from the dofs, can we retrieve
//what their dof carriers are?
//Well, the point is that for purpose of domain decomposition
// it is natural to loop over elements, and so it is natural 
// to say: what are the DOF CARRIERS in this element?
// Well, on the other hand, one might say: what are the DOFS in this element?
//But, since you are in an element, to retrieve the dofs you necessarily start from 
//the GEOMETRICAL ENTITIES that are available to you in there.
//So, from "DOF CARRIER" to DOF is a "GEOMETRICAL approach"
//Instead, from DOF to "DOF CARRIER" would be an "ALGEBRAIC approach"...,
//but I dont think that anybody uses that

//So the idea is that I always provide all the possible DOF CARRIERS,
//and then they may be used or not depending on the existence of variables of each FE type





//you see that this function is const, i.e. it doesnt modify any member of the class
//but actually, IT DOES!
//Simply, the class is made of POINTERS. Then, I've created a new pointer
//which is EQUAL to the pointer of the class, and then I MODIFY what is inside that!
//So, many functions here can be called CONST even if they are not!
  void EqnBase::Bc_ConvertToDirichletPenalty(const uint vb, const uint ql, uint* bc) const {

    const uint ndof  = _AbstractFE[ql]->_ndof[vb];
    const uint nvars = _nvars[ql];

    for (uint ivarq=0; ivarq < nvars; ivarq++) {
           for (uint d=0; d< ndof; d++)    { 
          const uint     indxq  =         d + ivarq*ndof;
	       bc[indxq] = 1 ;     }
    }
          
    return;
  }

//=================================
//This function is for NS type equations:
//it computes the flags for pressure and stress integrals 
//based on the pressure nodes
 void EqnBase::Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(const uint vb,
	    const uint *bc_eldofs,const QuantityLocal & Velold_in,const QuantityLocal& press_in,uint &press_fl) const {

	const uint el_ndof_p =  press_in._ndof[vb];
	const uint el_ndof_u =  Velold_in._ndof[vb];
	const uint   nvars_u =  Velold_in._dim;
        int press_sum=0;
           for (uint i=0; i< el_ndof_p; i++)   press_sum += bc_eldofs[nvars_u*el_ndof_u + i]; //only one linear variable... pay attention when trying linear-linear
            if ( press_sum ==0 )                 {  press_fl = 1; }   //ALL zeros: ONLY PRESSURE

  
    return;
  }

  
  
  

//============================
//glnode may be a multiple of the node number,
//since node->dof is a non-function map
//either because you have one node - one scalar dof, but vector variables
//or because you have one node - many dofs



///this function scales the passed dof vector and takes into account the boundary conditions as well
///a function like that can be useful also for multiplying/dividing by reference values or the like

void EqnBase::Bc_ScaleDofVec(NumericVector* myvec,  double ScaleFac /*, dimension */ ) {
//only works with Level = NoLevels - 1 , because bc is only on the finest level
  
  //pass the pointer to the array,the dimension of the array, and the scale factor  
//remember that these vectors we pick are made by quadratic and linear parts,
//so we actually know what kind of vector these are => the dimension can be taken from _Dim

//if you have to scale, scale  BOTH vector components, otherwise you could have a non div free vector after the scaling,
//because you scaled only one component.
//Even if you're controlling only with one component at the boundary, both components of the lifting function
//will change due to the variation of a single component at the boundary


for (uint i=0; i < _Dim[_NoLevels-1]; i++) { //loop over all the dofs, both quadratic and linear
  
  if (_bc[i] == 1 ) {  //if the dofs are not fixed, scale them
  
    myvec->set( i, (*myvec)(i)*ScaleFac );

  }

}
  
  
  
 return; 
}

//add only where boundary conditions are not fixed
void EqnBase::Bc_AddDofVec(NumericVector* vec_in,NumericVector* vec_out ) {

for (uint i=0; i < _Dim[_NoLevels-1]; i++) { 

    if (_bc[i] == 1 ) {
  
      vec_out->add (i,(*vec_in)(i));

    }
  
}

return;

}

void EqnBase::Bc_AddScaleDofVec(NumericVector* vec_in,NumericVector* vec_out,const double ScaleFac ) {
//add a vector multiplied by a constant (only where it is not fixed)
  
for (uint i=0; i < _Dim[_NoLevels-1]; i++) { 

    if (_bc[i] == 1 ) {
  
      vec_out->add (i,(*vec_in)(i)*ScaleFac);

    }
  
}

return;

}



// // // // void EqnBase::func_xyz() {

// // // // //   double xyz[dimension];
// // // //   double xez[dimension];
// // // //   
// // // //   //Give the point xyz in the 'dimensional' domain
// // // //   
// // // //   //NONdimensionalize THE COORDINATES,
// // // //   //but i think here they are already dimensionalized
// // // //   //do a code that depends EXPLICITLY from what was done before.
// // // //   //A good practice would be to check if something had been done before...
// // // //   
// // // //   
// // // //   // find in which element my point is 
// // // //   
// // // //   //get the dofs of that element, in some order
// // // // 
// // // //   //compute the coordinate transformation
// // // //   //as we have to transform the real coordinates into canonical ones,
// // // //   //we have to invert x = x(\xi)
// // // //   //so we have to compute \xi as a function of x
// // // // 
// // // // 
// // // // xez[0]=0.16;
// // // // xez[1]=0.4867;
// // // // 
// // // // double* lin_shapes = new double[_eldof[0][1]];
// // // // 
// // // // //void FE::shapes_xez(const uint ql, const uint el_ndof,const double* xez, double* shapes) {
// // // // //for now i'll take the linear from the quadratics...
// // // // _fe[0]->shapes_xez(1,_eldof[0][1],xez,lin_shapes);
// // // // 
// // // // 
// // // // //transform the nondimensional physical coordinates to the 
// // // //   //coordinates in the reference element (\xi \eta)
// // // //   
// // // //   //compute the shape functions at that reference elements
// // // // //first the domain was linear, now pick the quadratic shapes for velocity i.e.
// // // // double* quad_shapes = new double[_eldof[0][0]];
// // // // 
// // // // _fe[0]->shapes_xez(0,_eldof[0][0],xez,quad_shapes);
// // // // 
// // // // double* quad_dofs = new double[_eldof[0][0]];
// // // // for (uint i=0;i<_eldof[0][0];i++ ) quad_dofs[i]= 87.;
// // // // //clearly the order in which these dofs are given must be the same as the connectivity order!
// // // // //
// // // // 
// // // //   //do the sum and you get the value...
// // // //   //scalar product in R^n, not R^dim
// // // // 
// // // // double resu = _utils.dotN(quad_shapes,quad_dofs,_eldof[0][0]);
// // // // 
// // // // std::cout << "The final result is " <<  resu <<  std::endl;
// // // //   
// // // //   //given the computed vector of dofs
// // // //   //given the real point xyz
// // // //   //compute the function at that point
// // // //   //find which element I am in
// // // //   //retrieve from the global dof only the dofs of THAT element
// // // // //   for every dof, transform to the canonical element
// // // // //retrieve the shape fncs for every local dof
// // // // //translate the global coordinate into the local coordinate
// // // // //compute all the functions at the local coordinate
// // // // //multiply and you're done
// // // // //to test things you may first do them ALL INSIDE A FUNCTION,
// // // // // without calling from outside, then you put it in a 
// // // // //MODULAR FASHION.
// // // // 



  
  
//  return; 
// // // // }


#include "CurrElem.hpp"
#include "CurrGaussPoint.hpp"
//=================================================================
//=================================================================
//this function computes the integral of some function
//i would wanna do it without the quantity interface
//also do it with all the procs inside
//remember that the mesh is nondimensional here
// at the end you have to multiply again

//i do not want to need no Equation nor Quantity
//clearly in this way the FE order is not given,
//but if you use the function at gauss coordinates you do not need that basically
//clearly since i do not use any Quantity, i cannot consider the framework of the Domain
//and the refbox coordinates, this is just a simple function, for testing
//i think i cannot pass a function pointer to a member function,
//maybe i can do that only if the function is STATIC
//since here the order is not related to any equation, i can pick the  quadrature rule
//as INDEPENDENT of the ORDER
//in practice, here we are picking a 

double EqnBase::FunctionIntegral (const uint vb, double (*pt2func)(double, const double* ) ) {

const uint Level = _NoLevels - 1;
const uint myproc= _iproc;
  double time=0.;
  
    CurrElem       currelem(*this,_eqnmap);
//     CurrGaussPoint   currgp(_eqnmap);    
    CurrGaussPointBase & currgp = CurrGaussPointBase::build(_eqnmap, _mesh._dim);

  //======== ELEMENT MAPPING =======
  const uint meshql = (int) _mesh._mesh_rtmap.get("meshql");  
 
//========= DOMAIN MAPPING
    QuantityLocal xyz(currgp,currelem);
    xyz._dim      = _mesh._dim;
    xyz._FEord    = meshql;
    xyz._ndof[VV] = _AbstractFE[xyz._FEord]->_ndof[VV];
    xyz._ndof[BB] = _AbstractFE[xyz._FEord]->_ndof[BB];
    xyz._val_dofs = new double[xyz._dim*xyz._ndof[vb]];
    xyz._val_g    = new double[xyz._dim];

  double integral = 0.;
  
//loop over the geom el types
      const uint el_ngauss = _eqnmap._qrule._NoGaussVB[vb];   //elem gauss points  

//parallel sum
    const uint nel_e = _mesh._off_el[vb][_NoLevels*myproc+Level+1];
    const uint nel_b = _mesh._off_el[vb][_NoLevels*myproc+Level];
  
    for (uint iel=0; iel < (nel_e - nel_b); iel++) {
  
    currelem.get_el_nod_conn_lev_subd(vb,Level,myproc,iel);
    currelem.get_el_ctr(vb); 
    
    currelem.ConvertElemCoordsToMappingOrd(vb,xyz);

     
    for (uint qp = 0; qp < el_ngauss; qp++) {

for (uint fe = 0; fe < QL; fe++)     {  currgp.SetDPhiDxezetaElDofsFEVB_g (vb,fe,qp);  }  
     
double  Jac_g=0.;
          if (vb==0)   Jac_g = currgp.JacVectVV_g(vb,xyz);  //not xyz_refbox!      
     else if (vb==1)   Jac_g = currgp.JacVectBB_g(vb,xyz);  //not xyz_refbox!      

   const double  wgt_g = _eqnmap._qrule._weightVB[vb][qp];

     for (uint fe = 0; fe < QL; fe++)     {          currgp.SetPhiElDofsFEVB_g (vb,fe,qp);  }

 xyz.val_g(vb);
double myval_g = pt2func(time,xyz._val_g); 

 
  integral += wgt_g*Jac_g*myval_g;

   
    }//gauss loop
     
    }//element loop
    
         std::cout << std::endl << "vb = " << vb << " ^^^^^^^^^^^^^^^^^L'integrale sul processore "<< myproc << " vale: " << integral << std::endl;

    double weights_sum = 0.;
    for (uint qp = 0; qp < el_ngauss; qp++)  weights_sum += _eqnmap._qrule._weightVB[vb][qp];
       std::cout << std::endl << "vb = " << vb << " ^^^^^^^^^^^^^^^^^ La somma dei pesi  vale: " << weights_sum << std::endl;

       double J=0.;
   #ifdef HAVE_MPI
//       MPI_Reduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );  //This one gives J only to processor 0 !
      MPI_Allreduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  //THIS IS THE RIGHT ONE!!
   #else
   J = integral;
   #endif
    
     std::cout << std::endl << "vb = " << vb << " ^^^^^^^^^^^^^^^^^L'integrale totale vale: " << J << std::endl;

   
  return J;  
  
}
