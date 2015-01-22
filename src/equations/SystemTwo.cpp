/*=========================================================================

 Program: FEMUS
 Module: SystemTwo
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "SystemTwo.hpp"

// C++ 
#include <sstream>
#include <limits>

#include "FEMTTUConfig.h"
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "FemusDefault.hpp"

#include "MgSmootherEnum.hpp"
#include "NormTangEnum.hpp"
#include "QTYnumEnum.hpp"
#include "Quantity.hpp"
#include "IO.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "GeomEl.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "Files.hpp"
#include "CurrentElem.hpp"
#include "CurrentGaussPoint.hpp"

#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"
#include "DenseMatrix.hpp"


namespace femus {





//the most important things for an SystemTwo are:
//the number of variables
//the names
//other stuff but let us stop here now
SystemTwo::SystemTwo(std::vector<Quantity*> int_map_in,
                 MultiLevelProblemTwo& e_map_in,
                 std::string eqname_in,
                 std::string varname_in):
        _files(e_map_in._files),
        _phys(e_map_in._phys),
        _mesh(e_map_in._mesh),
        _eqnmap(e_map_in),
        //=============
        _eqname(eqname_in),
        _NoLevels(e_map_in._mesh._NoLevels), //you can do that
        _var_names(NULL),
        _refvalue(NULL),
        _dofmap(*this,e_map_in._mesh) {

//============= init Quantities ================
//internal std::vector
//the equal puts the two equal to each other
//as an alternative, since std::vector is a class, i could have used the initialization list above
//like i do with REFERENCES. In that case the copy constructor would be called.
    _QtyInternalVector = int_map_in;

//============= init n_vars================
    _dofmap.initNVars();
//========== varnames ==============
    initVarNames(varname_in);
//========== RefValues ==============
    initRefValues();
   
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
    _solver = new LinearEquationSolver*[_NoLevels];     //well, we clearly use the same package for all levels...
    for (uint l=0;l<_NoLevels;l++) _solver[l] = LinearEquationSolver::build(0,NULL,NO_SMOOTHER).release();
    

}


//====================
// by default, all the reference values are initialized to 1.
//This is a function that doesnt make distinction BETWEEN various FE,
// it treats the variables in the same manner
void SystemTwo::initVarNames(std::string varname_in) {

    assert(_dofmap._n_vars > 0);

    _var_names = new std::string[_dofmap._n_vars];       // names

    std::ostringstream name;
    for (uint i=0;i< _dofmap._n_vars; i++) { // variable name
        name.str("");
        name << varname_in << i+1;
        _var_names[i] = name.str();
    }

    return;
}



void SystemTwo::initRefValues() {

    assert(_dofmap._n_vars > 0);

    _refvalue  = new      double[_dofmap._n_vars];           // refvalues

    for (uint i=0;i< _dofmap._n_vars; i++) _refvalue[i] = 1.;

    return;
}



// ===================================================
/// This function  is the SystemTwo destructor.
//the important thing is that these destructions occur AFTER
//the destructions of the levels inside
//these things were allocated and filled in various init functions
//either in the Base or in the DA
//here, we destroy them in the place where they belong
// pay attention to the fact that a lot of these delete are ok only if the respective function
// where the new is is called!
SystemTwo::~SystemTwo() {

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




//==================
//this function DEPENDS in _iproc!!!
void SystemTwo::initVectors() {

    //allocation
         _x.resize(_NoLevels);
     _x_old.resize(_NoLevels);
    _x_oold.resize(_NoLevels);
     _x_tmp.resize(_NoLevels);
         _b.resize(_NoLevels);
       _res.resize(_NoLevels);

    for (uint Level = 0; Level< _NoLevels; Level++) {

    uint ml[QL];    for (int fe=0; fe<QL; fe++) ml[fe] = _dofmap._DofLocLevProcFE[Level][_iproc][fe];
    uint m_l = 0;
    for (int fe=0; fe<QL; fe++)  m_l +=  ml[fe]*_dofmap._nvars[fe];
    
        _b[Level] = NumericVector::build().release();
        _b[Level]->init(_dofmap._Dim[Level],m_l,false,AUTOMATIC);
        _res[Level] = NumericVector::build().release();
        _res[Level]->init(_dofmap._Dim[Level],m_l,false,AUTOMATIC);
        _x[Level] = NumericVector::build().release();
        _x[Level]->init(_dofmap._Dim[Level],m_l,false,AUTOMATIC);
        _x_old[Level] = NumericVector::build().release();
        _x_old[Level]->init(_dofmap._Dim[Level],false, SERIAL);
        _x_oold[Level] = NumericVector::build().release();
        _x_oold[Level]->init(_dofmap._Dim[Level],false, SERIAL);
        _x_tmp[Level] = NumericVector::build().release();
        _x_tmp[Level]->init(_dofmap._Dim[Level],false, SERIAL);

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

//

void SystemTwo::GenBc() {
  
    std::string    input_dir = DEFAULT_CONFIGDIR;
    std::string          ibc = DEFAULT_IBC;
    std::string     ext_xdmf = DEFAULT_EXT_XDMF;
    std::string       ext_h5 = DEFAULT_EXT_H5;
    std::string  bdry_suffix = DEFAULT_BDRY_SUFFIX;
    
    std::ostringstream ibc_fileh5;
    ibc_fileh5  << "./" << input_dir << "/" << ibc << ext_h5;
// TODO actually, we should first COPY this file in the outtime dir, then READ it from there!
    std::ifstream in(ibc_fileh5.str().c_str());

    const uint Lev_pick_bc_NODE_dof = _NoLevels-1;  //we use the FINE Level as reference
    
 //************************************************   
 //******** NODE BASED ****************************   
    const uint mesh_ord    = (int) _mesh.GetRuntimeMap().get("mesh_ord");
    const uint offset      = _mesh._NoNodesXLev[_NoLevels-1];
    const uint el_nnodes_b = _mesh.GetGeomEl(_mesh.get_dim()-1-BB,mesh_ord)._elnds;
    double* normal = new double[_mesh.get_dim()];  //TODO remove this, it is useless

    int  *bc_flag = new int[_dofmap._n_vars];

    _bc = new int[_dofmap._Dim[Lev_pick_bc_NODE_dof]];
    for (uint i1=0;i1< _dofmap._Dim[Lev_pick_bc_NODE_dof];i1++) _bc[i1] = DEFAULT_BC_FLAG;    // set 1 all the points for  bc (boundary condition)

    //both here and in bc_read you have to put the value bc=0
    //so that you get the IDENTITY OPERATOR and you only have to provide the
    //function for the RHS
 //******** NODE BASED ****************************   
 //***************************************************   

 //**************************************************   
 //******** ELEM BASED ******************************  
     CurrentElem       currelem(BB,this,_mesh,_eqnmap._elem_type);

    _bc_fe_kk             =  new int*[_NoLevels];
    int* DofOff_Lev_kk    =  new int[_NoLevels];
    
    for (uint Level=0; Level <_NoLevels; Level++)   { //loop over the levels

          DofOff_Lev_kk[Level] = _dofmap._nvars[KK]*_dofmap._DofNumLevFE[Level][KK];
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
            uint iel_b = _mesh._off_el[BB][ _NoLevels*isubd + Level];
            uint iel_e = _mesh._off_el[BB][ _NoLevels*isubd + Level+1];
            for (uint iel=0; iel < (iel_e - iel_b); iel++) {

	        currelem.set_el_nod_conn_lev_subd(Level,isubd,iel);
                currelem.SetMidpoint();
		
 	    for (uint ivar=0; ivar< _dofmap._n_vars; ivar++)  bc_flag[ivar] = DEFAULT_BC_FLAG; //this is necessary here to re-clean!

                 if (_dofmap._n_vars >0) bc_read(currelem.GetMidpoint(),normal,bc_flag);

  //******************* ONLY FINE LEVEL, NODE VARS ***************** 
   if (Level == Lev_pick_bc_NODE_dof)  { 
                for (uint i=0; i<  el_nnodes_b; i++)  {
                        const uint fine_node = _mesh._el_map[BB][(iel+iel_b)*el_nnodes_b+i];

                    //Set the quadratic fields
                    if ( i < _eqnmap._elem_type[_mesh.get_dim()-1-BB][QQ]->GetNDofs() )
		      for (uint ivar=0; ivar<_dofmap._nvars[QQ]; ivar++) {
                            int kdof = _dofmap.GetDof(Lev_pick_bc_NODE_dof,QQ,ivar,fine_node);
                           if (_bc[kdof] != 0) _bc[kdof] = bc_flag[ ivar + _dofmap._VarOff[QQ]];
                        }
                    // Set the linear fields
                    if ( i < _eqnmap._elem_type[_mesh.get_dim()-1-BB][LL]->GetNDofs() ) {
                        for (uint ivar = 0; ivar < _dofmap._nvars[LL]; ivar++) {
                            int kdof = _dofmap.GetDof(Lev_pick_bc_NODE_dof,LL,ivar,fine_node);
                           if (_bc[kdof] != 0) _bc[kdof] = bc_flag[ ivar + _dofmap._VarOff[LL]];
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
		 
 	        for (uint ivar= 0; ivar < _dofmap._nvars[KK];  ivar++)  { 
		  int dof_kk_pos_lev =  vol_iel + ivar*_dofmap._DofNumLevFE[Level][KK];
	            if (_bc_fe_kk[Level][ dof_kk_pos_lev ] != 0)  _bc_fe_kk[Level][ dof_kk_pos_lev ] = bc_flag[ ivar + _dofmap._VarOff[KK] ];
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
        for (uint ivar=0;ivar< _dofmap._nvars[QQ];ivar++) { // ----------------------------
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
                        int kdof= _dofmap.GetDof(_NoLevels-1,QQ,ivar,k);
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
//---- Why does the get_el_nod_conn(BB,..) work in GenMatRhs
//but HERE it doesnt?
//I think because here we are looping over ALL LEVELS and ALL SUBDOMAINS (procs)
//while in the GenMatRhs every processor
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
void SystemTwo::GenElBc()  {

     CurrentElem       currelem(BB,this,_mesh,_eqnmap._elem_type);  
  
      uint space_dim = _mesh.get_dim();

    _elem_bc       =  new int**[_NoLevels];
    _elem_val_norm =  new double**[_NoLevels];
    _elem_val_tg   =  new double**[_NoLevels];

    for (uint Level=0; Level <_NoLevels;Level++)   { //loop over the levels

        _elem_bc[Level]       = new int*[_mesh._NoSubdom];//4*_mesh._n_elements_vb_lev[BB][Level] this was wrong, actually it is L + P*NoLevels, i wasnt considering the others but it was working! For instance in two D the numbers for the two procs are the same, here's why!
        _elem_val_norm[Level] = new double*[_mesh._NoSubdom];//_mesh._n_elements_vb_lev[BB][Level]
        _elem_val_tg[Level]   = new double*[_mesh._NoSubdom];//_mesh._n_elements_vb_lev[BB][Level]

        for (uint isubd=0;isubd<_mesh._NoSubdom;++isubd) {  //loop over the subdomains

            uint iel_b=_mesh._off_el[BB][ _NoLevels*isubd +  Level];
            uint iel_e=_mesh._off_el[BB][ _NoLevels*isubd +  Level+1];

            _elem_bc[Level][isubd]  = new int[/*4*/2*(iel_e-iel_b)]; /*normal and tangential*/ //4*_mesh._n_elements_vb_lev[BB][Level+isubd*_NoLevels] that was not correct
            _elem_val_norm[Level][isubd]  = new double[ 1*(iel_e-iel_b) ];  //_mesh._n_elements_vb_lev[BB][Level+isubd*_NoLevels]

            _elem_val_tg[Level][isubd]  = new double[ _number_tang_comps[space_dim - 1]*(iel_e-iel_b) ];   //_mesh._n_elements_vb_lev[BB][Level+isubd*_NoLevels]

            for (uint iel=0;iel < (iel_e - iel_b); iel++) {  //loop over the elems of that level&subdomain

                int surf_id=0;  //could do this outside also
                int     el_flag[2] = {0,0};
                std::vector<double> el_value(1 + _number_tang_comps[space_dim - 1],0.); //1 normal and 1 tangential or 1 normal and 3 tangential

                currelem.set_el_nod_conn_lev_subd(Level,isubd,iel);
                currelem.SetMidpoint();

                //read the bc's //the read forgets all levels and subdomains, it is only based on the MIDDLE POINT
                elem_bc_read(currelem.GetMidpoint(),surf_id,&el_value[0],el_flag);

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


///////////////////////////////
void SystemTwo::Bc_GetElFlagValLevSubd(const uint Level,const uint isubd,const uint iel,int* el_flag,double* el_value ) const {

    el_flag[0] =        _elem_bc[Level][isubd][/*4*/2*iel]   ;
    el_flag[1] =        _elem_bc[Level][isubd][/*4*/2*iel+1] ;
    /*    el_flag[2] =        _elem_bc[Level][isubd][4*iel+2] ;
        el_flag[3] =        _elem_bc[Level][isubd][4*iel+3] ;*/
    el_value[0] =  _elem_val_norm[Level][isubd][iel]   ;

  uint space_dim = _mesh.get_dim();

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
void SystemTwo::clearElBc() {

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
void SystemTwo::GenIc() {

    const uint mesh_ord    = (int) _mesh.GetRuntimeMap().get("mesh_ord");
  
    std::string  input_dir = DEFAULT_CONFIGDIR;
    std::string        ibc = DEFAULT_IBC;
    std::string     ext_h5 = DEFAULT_EXT_H5;
    std::ostringstream ibc_filexmf;
    ibc_filexmf << "./"<< input_dir << ibc << ext_h5;
    std::ifstream in(ibc_filexmf.str().c_str());

    if (!in) {

        CurrentElem       currelem(VV,this,_mesh,_eqnmap._elem_type);  
     
        const uint  coords_fine_offset = _mesh._NoNodesXLev[_NoLevels-1];
        const uint  el_nnodes = _mesh.GetGeomEl(_mesh.get_dim()-1-VV,mesh_ord)._elnds;
        double*      xp = new double[_mesh.get_dim()];
        double* u_value = new double[_dofmap._n_vars];

       std::cout << "\n====================== GenIc:  Now we are setting them for all levels! ========================" << "\n \n";

    for (uint Level = 0; Level< _NoLevels; Level++) {

            uint iel_b = _mesh._off_el[VV][ _iproc*_NoLevels + Level ];
            uint iel_e = _mesh._off_el[VV][ _iproc*_NoLevels + Level + 1];

	    for (uint iel=0; iel < (iel_e - iel_b); iel++) {
	  
	        currelem.set_el_nod_conn_lev_subd(Level,_iproc,iel);
                currelem.SetMidpoint();
	
            // we are looping over the mesh nodes, but it's a fake loop because we do not depend on "i" for the elements
            for (uint i=0; i < el_nnodes ; i++) {
                int fine_node = _mesh._el_map[VV][ i + ( iel + iel_b )*el_nnodes ];
                for (uint idim=0; idim<_mesh.get_dim(); idim++) xp[idim] = _mesh._xyz[ fine_node + idim*coords_fine_offset ];
                for (uint ivar=0; ivar < _dofmap._n_vars; ivar++) u_value[ivar] = 0.;

                // ===================================================
                // user definition reading function ----------------
                ic_read(xp,u_value,currelem.GetMidpoint());
                // -------------------------------------------------
                // ===================================================

                // Set the quadratic fields
                if ( i < _eqnmap._elem_type[_mesh.get_dim()-1-VV][QQ]->GetNDofs() ) {
                    for (uint ivar=0; ivar < _dofmap._nvars[QQ]; ivar++) {
		      int dof_pos_lev = _dofmap.GetDof(Level,QQ,ivar,fine_node);
                        _x[Level]->set( dof_pos_lev, u_value[ivar + _dofmap._VarOff[QQ]] );
		    }
                }
                // Set the linear fields
                if ( i <  _eqnmap._elem_type[_mesh.get_dim()-1-VV][LL]->GetNDofs() ) {
                    for (uint ivar=0; ivar<_dofmap._nvars[LL]; ivar++) {
		      int dof_pos_lev = _dofmap.GetDof(Level,LL,ivar,fine_node);
                        _x[Level]->set( dof_pos_lev, u_value[ivar + _dofmap._VarOff[LL]] );
                        }
		    }
		    
                if ( i < _eqnmap._elem_type[_mesh.get_dim()-1-VV][KK]->GetNDofs() ) {
		  
		int sum_elems_prev_sd_at_lev = 0;
	    for (uint pr = 0; pr < _iproc; pr++) { sum_elems_prev_sd_at_lev += _mesh._off_el[VV][pr*_NoLevels + Level + 1] - _mesh._off_el[VV][pr*_NoLevels + Level]; }
	    
		  //elem my level
                    for (uint ivar=0; ivar<_dofmap._nvars[KK]; ivar++) {
		      int elem_lev = iel + sum_elems_prev_sd_at_lev;
		      int dof_pos_lev = _dofmap.GetDof(Level,KK,ivar,elem_lev);
                        _x[Level]->set( dof_pos_lev, u_value[ivar + _dofmap._VarOff[KK]] );
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

        IO::read_system_solutions(ibc_filexmf.str(),&_mesh,&_dofmap,this);
        std::cout << "\n GenIc(Base): Initial solution defined by " <<  ibc_filexmf.str().c_str() << "\n \n";
    }

    in.close();

    return;
}




/// initial conditions  from function
// void SystemTwo::ic_read(const double * /*xp*/, double * ic, const double * /*el_xm*/) const {
//     for (uint ivar=0;ivar<_n_vars;ivar++) ic[ivar]=0.;
//     return;
// }






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
void SystemTwo::MGSolve(double Eps1,          // tolerance for the linear solver
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


double SystemTwo::MGStep(int Level,            // Level
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
void SystemTwo::MGCheck(int Level) const {

    _x[Level-1]->matrix_mult(*_x[Level],*_Rst[Level-1]);
    _x[Level]  ->matrix_mult(*_x[Level-1],*_Prl[Level]);
    return;
}




// =========================================
void SystemTwo::ReadMGOps() {

    std::string     f_matrix = DEFAULT_F_MATRIX;
    std::string       f_rest = DEFAULT_F_REST;
    std::string       f_prol = DEFAULT_F_PROL;
    std::string       ext_h5 = DEFAULT_EXT_H5;

    std::ostringstream filename;
    std::string filename_base;
    filename_base = _files._output_path + "/";
    
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


void SystemTwo::ReadMatrix(const  std::string& namefile) {

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
            IO::read_UIhdf5(file,dim_name.str().c_str(),rowcln[r][c]);

            //row length
            std::ostringstream len_name;
            len_name << groupname_lev.str() << "/" << "LEN" << fe_couple.str();
            length_row[r][c]=new uint[ rowcln[r][c][0]+1 ];
            IO::read_UIhdf5(file,len_name.str().c_str(),length_row[r][c]);

            // matrix off diagonal
            std::ostringstream offlen_name;
            offlen_name << groupname_lev.str() << "/" << "OFFLEN" << fe_couple.str();
            length_offrow[r][c]=new uint[ rowcln[r][c][0]+1 ];
            IO::read_UIhdf5(file,offlen_name.str().c_str(),length_offrow[r][c]);

            // matrix pos //must stay AFTER reading length_offrow
            std::ostringstream pos_name;
            pos_name << groupname_lev.str() << "/" << "POS" << fe_couple.str();
            pos_row[r][c]=new uint[ length_row[r][c][rowcln[r][c][0]] ];
            IO::read_UIhdf5(file,pos_name.str().c_str(),pos_row[r][c]);

        } //end col
    } //end row
    
    H5Fclose(file);

//============================================================================
//============ compute things for the sparsity pattern =======================
//============================================================================
    
    int NoLevels  = _mesh._NoLevels;
    uint off_proc = NoLevels*_iproc;

    uint mrow_glob_t = 0;
    for (int fe=0; fe<QL; fe++) mrow_glob_t += _dofmap._nvars[fe]*rowcln[fe][fe][0];
    uint ncol_glob_t = mrow_glob_t;

    uint mrow_lev_proc_t  = 0;
    for (int fe=0; fe<QL; fe++)  mrow_lev_proc_t +=  _dofmap._DofLocLevProcFE[Level][_iproc][fe]*_dofmap._nvars[fe];
    uint ncol_lev_proc_t  = mrow_lev_proc_t;

    uint DofObjInit_lev_PrevProcs[QL];  //TODO what is this? it is the ROW INDEX at which to begin for every processor
                      
     for (int r=0; r<QL; r++)     DofObjInit_lev_PrevProcs[r] = 0;
         
    for (uint isubd=0; isubd<_iproc; isubd++) {
        DofObjInit_lev_PrevProcs[QQ] += _dofmap._DofLocLevProcFE[Level][isubd][QQ];
        DofObjInit_lev_PrevProcs[LL] += _dofmap._DofLocLevProcFE[Level][isubd][LL];
        DofObjInit_lev_PrevProcs[KK] += _dofmap._DofLocLevProcFE[Level][isubd][KK];
    }    
    
    
    _A[Level] = SparseMatrix::build().release();
// //     _A[Level]->init(_Dim[Level],_Dim[Level], mrow_lev_proc_t, mrow_lev_proc_t); //TODO BACK TO a REASONABLE INIT

    Graph graph;
    graph.resize(mrow_glob_t);
    graph._m  = mrow_glob_t;
    graph._n  = ncol_glob_t;
    graph._ml = mrow_lev_proc_t;
    graph._nl = ncol_lev_proc_t;
    graph._ml_start = _dofmap.GetStartDof(Level,off_proc);
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
    off_onevar[KK] = _mesh._n_elements_vb_lev[VV][Level];
    
    uint  off_EachFEFromStart[QL];
    off_EachFEFromStart[QQ] = 0;
    off_EachFEFromStart[LL] = _dofmap._nvars[QQ]*off_onevar[QQ];
    off_EachFEFromStart[KK] = _dofmap._nvars[QQ]*off_onevar[QQ] + _dofmap._nvars[LL]*off_onevar[LL];

 //==============   
     for (int r=0;r<QL;r++) {
      for (uint ivar=0; ivar < _dofmap._nvars[r]; ivar++) {

        for (uint DofObj_lev = DofObjInit_lev_PrevProcs[r]; DofObj_lev < DofObjInit_lev_PrevProcs[r] + _dofmap._DofLocLevProcFE[Level][_iproc][r]; DofObj_lev++) {

            int dof_pos, irow;
	         if  (r<KK) {  dof_pos = _mesh._Qnode_lev_Qnode_fine[FELevel[r]][ DofObj_lev ];  }
            else if (r==KK) {  dof_pos = DofObj_lev; }
                    irow = _dofmap.GetDof(Level,r,ivar,dof_pos); 

	    int len[QL];   for (int c=0;c<QL;c++) len[c] = 0;
	    for (int c=0;c<QL;c++) len[c] = (length_row[r][c][ DofObj_lev+1 ]-length_row[r][c][ DofObj_lev ]);
            int rowsize = 0; 
	    for (int c=0;c<QL;c++) rowsize +=_dofmap._nvars[c]*len[c];
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
	    for (int c=0;c<QL;c++) lenoff_size += _dofmap._nvars[c]*lenoff[c];
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
void SystemTwo::ReadProl(const std::string& name) {

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
	  nrowt += _dofmap._nvars[fe]*rowcln[fe][0];
          nclnt += _dofmap._nvars[fe]*rowcln[fe][1];
	}
    
    Graph pattern;
    pattern.resize(nrowt);
    pattern._m=nrowt;
    pattern._n=nclnt;
    pattern._ml = 0;
    pattern._nl = 0;
     for (int fe=0;fe<QL;fe++) { 
        pattern._ml += _dofmap._nvars[fe]*ml[fe];  //  local _m
        pattern._nl += _dofmap._nvars[fe]*nl[fe];  //  local _n
     }
    uint ml_start = _dofmap.GetStartDof(Level,off_proc);
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

     for (uint ivar=0;ivar<_dofmap._nvars[fe];ivar++) {
        for (unsigned int i = ml_init[fe]; i < ml_init[fe]+ml[fe]; i++) {
          int dof_pos_f;
          if (fe < KK)         dof_pos_f = _mesh._Qnode_lev_Qnode_fine[ FEXLevel_f[fe] ][ i ];  //end fe < ql
          else if (fe == KK)   dof_pos_f = i;
          
            int irow  = _dofmap.GetDof(Lev_f,fe,ivar,dof_pos_f);

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
	      
	      pattern[irow][j] = _dofmap.GetDof(Lev_c,fe,ivar,dof_pos_c);
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
       for (uint ivar=0;ivar < _dofmap._nvars[fe];ivar++) {
          for (unsigned int i=ml_init[fe]; i<ml_init[fe] + ml[fe]; i++) {
          int dof_pos_f;
          if (fe < KK)        dof_pos_f = _mesh._Qnode_lev_Qnode_fine[ FEXLevel_f[fe] ][ i ];  //end fe < ql
          else if (fe == KK)  dof_pos_f = i;

          int irow  = _dofmap.GetDof(Lev_f,fe,ivar,dof_pos_f);

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
void SystemTwo::ReadRest(const std::string& name) {
 
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
	  nrowt += _dofmap._nvars[fe]*rowcln[fe][0];
          nclnt += _dofmap._nvars[fe]*rowcln[fe][1];
	}


    Graph pattern;
    pattern.resize(nrowt);
    pattern._m  = nrowt;
    pattern._n  = nclnt;        // global dim _m x _n
    pattern._ml = 0;            //  local _m
    pattern._nl = 0;            //  local _n
     for (int fe=0;fe<QL;fe++) { 
        pattern._ml += _dofmap._nvars[fe]*_dofmap._DofLocLevProcFE[Lev_c][_iproc][fe]; 
        pattern._nl += _dofmap._nvars[fe]*_dofmap._DofLocLevProcFE[Lev_f][_iproc][fe];
     } 
     
    // starting indices for local matrix
   uint ml_start = _dofmap.GetStartDof(Lev_c,off_proc);    //  offset proc nodes
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
    for (uint ivar=0;ivar<_dofmap._nvars[fe];ivar++) {
        for (unsigned int i = DofObjInit_lev_PrevProcs_c[fe]; i< DofObjInit_lev_PrevProcs_c[fe] + _dofmap._DofLocLevProcFE[Lev_c][_iproc][fe]; i++) {
	  int dof_pos_c;
          if (fe < KK)         dof_pos_c = _mesh._Qnode_lev_Qnode_fine[ FEXLevel_c[fe] ][ i ];
          else if (fe == KK)   dof_pos_c = i;
	  
            int irow  = _dofmap.GetDof(Lev_c,fe,ivar,dof_pos_c);

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
              pattern[irow][j] = _dofmap.GetDof(Lev_f,fe,ivar,dof_pos_f);
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

    for (uint ivar=0;ivar<_dofmap._nvars[fe];ivar++) {
        for (unsigned int i = DofObjInit_lev_PrevProcs_c[fe]; i< DofObjInit_lev_PrevProcs_c[fe] + _dofmap._DofLocLevProcFE[Lev_c][_iproc][fe]; i++) {
	  int dof_pos_c;
          if (fe < KK)         dof_pos_c = _mesh._Qnode_lev_Qnode_fine[ FEXLevel_c[fe] ][ i ];
          else if (fe == KK)   dof_pos_c = i;
            int irow     = _dofmap.GetDof(Lev_c,fe,ivar,dof_pos_c);
            int irow_top = _dofmap.GetDof(_NoLevels-1,fe,ivar,dof_pos_c);
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









//here, we could pass a std::VECTOR of the UNKNOWN INTERNAL QUANTITIES,
//not a map because we want things to be ordered.
//this vector must be passed to the EqnNS constructor so that it can be passed
//to the DA. 
//so we must do things OUTSIDE the NS constructor, i.e. in the MultiLevelProblemTwo.
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
// and THEN FILLED in the GenMatRhs
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
  void SystemTwo::Bc_ConvertToDirichletPenalty(const uint elem_dim, const uint ql, uint* bc) const {

    const uint ndof  = _eqnmap._elem_type[elem_dim-1][ql]->GetNDofs();
    const uint nvars = _dofmap._nvars[ql];

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
 void SystemTwo::Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(
	    const uint *bc_eldofs,const CurrentQuantity & Velold_in,const CurrentQuantity& press_in,uint &press_fl) const {

	const uint el_ndof_p =  press_in._ndof;
	const uint el_ndof_u =  Velold_in._ndof;
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

void SystemTwo::Bc_ScaleDofVec(NumericVector* myvec,  double ScaleFac /*, dimension */ ) {
//only works with Level = NoLevels - 1 , because bc is only on the finest level
  
  //pass the pointer to the array,the dimension of the array, and the scale factor  
//remember that these vectors we pick are made by quadratic and linear parts,
//so we actually know what kind of vector these are => the dimension can be taken from _Dim

//if you have to scale, scale  BOTH vector components, otherwise you could have a non div free vector after the scaling,
//because you scaled only one component.
//Even if you're controlling only with one component at the boundary, both components of the lifting function
//will change due to the variation of a single component at the boundary


for (uint i=0; i < _dofmap._Dim[_NoLevels-1]; i++) { //loop over all the dofs, both quadratic and linear
  
  if (_bc[i] == 1 ) {  //if the dofs are not fixed, scale them
  
    myvec->set( i, (*myvec)(i)*ScaleFac );

  }

}
  
  
  
 return; 
}

//add only where boundary conditions are not fixed
void SystemTwo::Bc_AddDofVec(NumericVector* vec_in,NumericVector* vec_out ) {

for (uint i=0; i < _dofmap._Dim[_NoLevels-1]; i++) { 

    if (_bc[i] == 1 ) {
  
      vec_out->add (i,(*vec_in)(i));

    }
  
}

return;

}

void SystemTwo::Bc_AddScaleDofVec(NumericVector* vec_in,NumericVector* vec_out,const double ScaleFac ) {
//add a vector multiplied by a constant (only where it is not fixed)
  
for (uint i=0; i < _dofmap._Dim[_NoLevels-1]; i++) { 

    if (_bc[i] == 1 ) {
  
      vec_out->add (i,(*vec_in)(i)*ScaleFac);

    }
  
}

return;

}



// // // // void SystemTwo::func_xyz() {

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




} //end namespace femus
