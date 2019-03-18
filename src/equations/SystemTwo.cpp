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
#include <cassert>

#include "FemusConfig.hpp"
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "FemusDefault.hpp"

#include "LinearEquationSolverEnum.hpp"
#include "NormTangEnum.hpp"
#include "XDMFWriter.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "MultiLevelProblem.hpp"
#include "CurrentElem.hpp"
#include "CurrentGaussPoint.hpp"

#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"
#include "DenseMatrix.hpp"




namespace femus {


//the most important things for a SystemTwo are:
//the number of variables
//the names
//other stuff but let us stop here now
SystemTwo::SystemTwo(MultiLevelProblem& e_map_in, const std::string & eqname_in, const unsigned int number, const LinearEquationSolverType & smoother_type):
        _dofmap(this,e_map_in.GetMeshTwo()),
        _bcond(&_dofmap),
        NonLinearImplicitSystem(e_map_in,eqname_in,number,smoother_type) { }


void SystemTwo::init_unknown_vars() {

//============= init n_vars================
    _dofmap.initNVars();
//========== varnames ==============
    initVarNames();
//========== RefValues ==============
    initRefValues();

    return;
 }

//====================
// by default, all the reference values are initialized to 1.
//This is a function that doesnt make distinction BETWEEN various FE,
// it treats the variables in the same manner
void SystemTwo::initVarNames() {

    assert(_dofmap._n_vars > 0);

    _var_names.resize(_dofmap._n_vars);       // names
    
//=======  _var_names: they are the names of the quantities which are unkwnowns to this equation  ===========
    std::ostringstream name;
    uint count = 0;
   for (uint i=0; i< _UnknownQuantitiesVector.size(); i++) {
        for (uint j=0; j < _UnknownQuantitiesVector[i]->_dim; j++) {
          name.str("");
          name << _UnknownQuantitiesVector[i]->_name << j;
	  _var_names[count] = name.str();
	  count++;
	}
     }
   
    return;
}



void SystemTwo::initRefValues() {

    assert(_dofmap._n_vars > 0);

    _refvalue.resize(_dofmap._n_vars);

    uint count = 0;
   for (uint i=0; i< _UnknownQuantitiesVector.size(); i++) {
        for (uint j=0; j < _UnknownQuantitiesVector[i]->_dim; j++) {
	  _refvalue[count] = _UnknownQuantitiesVector[i]->_refvalue[j];
	  count++;
	}
   }
    
    return;
}








//==================
//this function DEPENDS in _iproc!!!
void SystemTwo::initVectors() {


    for (uint Level = 0; Level< GetGridn(); Level++) {

    uint ml[QL];    for (int fe=0; fe<QL; fe++) ml[fe] = _dofmap._DofLocLevProcFE[Level][GetMLProb().GetMeshTwo()._iproc][fe];
    uint m_l = 0;
    for (int fe=0; fe<QL; fe++)  m_l +=  ml[fe]*_dofmap._nvars[fe];
    
        _LinSolver[Level]->_RESC = NumericVector::build().release();
        _LinSolver[Level]->_RESC->init(_dofmap._Dim[Level],m_l,false,AUTOMATIC);
        _LinSolver[Level]->_RES = NumericVector::build().release();
        _LinSolver[Level]->_RES->init(_dofmap._Dim[Level],m_l,false,AUTOMATIC);
        _LinSolver[Level]->_EPS = NumericVector::build().release();
        _LinSolver[Level]->_EPS->init(_dofmap._Dim[Level],m_l,false,AUTOMATIC);
        _LinSolver[Level]->_EPSC = NumericVector::build().release();
        _LinSolver[Level]->_EPSC->init(_dofmap._Dim[Level],false, SERIAL);

    } //end level loop
    

    return;
}

 



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

//TTOODDOO AAA: oh this is very interesting... the initial conditions are set IN PARALLEL,
// only by the current processor!

//TTOODDOO ok the point is this: we are looping over elements, and then over nodes in the element,
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

//TTOODDOO: so far we have a node_dof that is explored based on the FINEST LEVEL for the NODES,
// and based on EACH LEVEL for the elements. 
//So, if we want to make it more general, we just put the offsets to be dependent on LEVEL,
// and for the nodes they are simply NOT DEPENDENT!

//TTOODDOO this function's goal is basically just to fill 
// the x vector in parallel, and then put it into x_old.
//Of course it does this in a GEOMETRIC WAY.
// xold is a vector of DOFS, so what we need are the 
// GEOMETRIC COORDINATES of the DOF OBJECTS (Node, Elem)
// associated to these dofs. In this way we can 

// I want to change the loop for initializing a vector
//

/// This function generates the initial conditions:
void SystemTwo::Initialize() {

      
        const uint  coords_fine_offset = GetMLProb().GetMeshTwo()._NoNodesXLev[GetGridn()-1];
        std::vector<double>      xp(GetMLProb().GetMeshTwo().get_dim());
        std::vector<double> u_value(_dofmap._n_vars);

       std::cout << "\n====================== Initialize:  Now we are setting them for all levels! ========================" << "\n \n";

    for (uint Level = 0; Level< GetGridn(); Level++) {
      
            Mesh	*mymsh	=  GetMLProb()._ml_msh->GetLevel(Level);
            const unsigned myproc  = mymsh->processor_id();

            uint iel_b = GetMLProb().GetMeshTwo()._off_el[VV][ GetMLProb().GetMeshTwo()._iproc*GetGridn() + Level ];
            uint iel_e = GetMLProb().GetMeshTwo()._off_el[VV][ GetMLProb().GetMeshTwo()._iproc*GetGridn() + Level + 1];

	    for (uint iel=0; iel < (iel_e - iel_b); iel++) {
	  
                CurrentElem       currelem(iel,myproc,Level,VV,this,GetMLProb().GetMeshTwo(),GetMLProb().GetElemType(),mymsh);  
	
	        currelem.SetDofobjConnCoords();
		
                const uint  el_dof_objs = NVE[ GetMLProb().GetMeshTwo()._geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];

            for (uint q=0; q < _UnknownQuantitiesVector.size() ; q++) {
		      
	  std::vector<double>  value(_UnknownQuantitiesVector[q]->_dim,0.);
            //the fact is that THERE ARE DIFFERENT DOF OBJECTS for DIFFERENT FE families
	    //for each family we should only pick the dof objects that are needed
	    //what changes between the FE families is the DOF OBJECT YOU PROVIDE: it could be a NODE or a CELL
	    // Notice that for some elements you don't have the midpoint of the element!
       
        for (uint ivar=0; ivar < _UnknownQuantitiesVector[q]->_dim; ivar++) {
       
          for (uint k=0; k < currelem.GetElemType(_UnknownQuantitiesVector[q]->_FEord)->GetNDofs() ; k++) {
	    
                const int fine_node = GetMLProb().GetMeshTwo()._el_map[VV][ k + ( iel + iel_b )*el_dof_objs ];
                for (uint idim = 0; idim < GetMLProb().GetMeshTwo().get_dim(); idim++) xp[idim] = GetMLProb().GetMeshTwo()._xyz[ fine_node + idim*coords_fine_offset ];

	  _UnknownQuantitiesVector[q]->initialize_xyz(&xp[0],value);
		    const int dof_pos_lev = _dofmap.GetDofQuantityComponent(Level,_UnknownQuantitiesVector[q],ivar,fine_node);
		    _LinSolver[Level]->_EPS->set( dof_pos_lev, value[ivar] );
		      
	    }  //dof objects 
         }

              } //qty
        
        
        } // end of element loop

        _LinSolver[Level]->_EPS->localize(*_LinSolver[Level]->_EPSC);
        _LinSolver[Level]->_EPSC->close();
	
    } //end Level
    

#ifdef DEFAULT_PRINT_INFO
        std::cout << "\n Initialize(Base): Initial solution defined by ic_read" << "\n \n";
#endif

    return;
}







} //end namespace femus
