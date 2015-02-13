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

#include "FEMTTUConfig.h"
#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "FemusDefault.hpp"

#include "MgSmootherEnum.hpp"
#include "NormTangEnum.hpp"
#include "QTYnumEnum.hpp"
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
SystemTwo::SystemTwo(MultiLevelProblem& e_map_in, const std::string & eqname_in, const unsigned int number, const MgSmoother & smoother_type):
        _dofmap(this,e_map_in.GetMeshTwo()),
        _bcond(&_dofmap),
        LinearImplicitSystem(e_map_in,eqname_in,number,smoother_type) { }


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
    for (uint Level =0; Level< GetGridn(); Level++) {
        delete _A[Level];
        if (Level < GetGridn() - 1) delete _Rst[Level];
        if (Level > 0)             delete _Prl[Level];
    }

    _A.clear();
    _Rst.clear();
    _Prl.clear();

 //======== Vectors ==========================
    for (uint Level =0; Level<GetGridn(); Level++) {
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

   
}



void SystemTwo::init_sys() {

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

    //allocation
         _x.resize(GetGridn());
     _x_old.resize(GetGridn());
    _x_oold.resize(GetGridn());
     _x_tmp.resize(GetGridn());
         _b.resize(GetGridn());
       _res.resize(GetGridn());

    for (uint Level = 0; Level< GetGridn(); Level++) {

    uint ml[QL];    for (int fe=0; fe<QL; fe++) ml[fe] = _dofmap._DofLocLevProcFE[Level][GetMLProb().GetMeshTwo()._iproc][fe];
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

// I want to change the loop for initializing a vector
//

/// This function generates the initial conditions:
void SystemTwo::Initialize() {

      
        const uint  coords_fine_offset = GetMLProb().GetMeshTwo()._NoNodesXLev[GetGridn()-1];
        std::vector<double>      xp(GetMLProb().GetMeshTwo().get_dim());
        std::vector<double> u_value(_dofmap._n_vars);

       std::cout << "\n====================== Initialize:  Now we are setting them for all levels! ========================" << "\n \n";

    for (uint Level = 0; Level< GetGridn(); Level++) {
      
       CurrentElem       currelem(Level,VV,this,GetMLProb().GetMeshTwo(),GetMLProb().GetElemType());  
        const uint  el_dof_objs = NVE[ GetMLProb().GetMeshTwo()._geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];

            uint iel_b = GetMLProb().GetMeshTwo()._off_el[VV][ GetMLProb().GetMeshTwo()._iproc*GetGridn() + Level ];
            uint iel_e = GetMLProb().GetMeshTwo()._off_el[VV][ GetMLProb().GetMeshTwo()._iproc*GetGridn() + Level + 1];

	    for (uint iel=0; iel < (iel_e - iel_b); iel++) {
	  
	        currelem.SetDofobjConnCoords(GetMLProb().GetMeshTwo()._iproc,iel);
                currelem.SetMidpoint();

            for (uint q=0; q < _UnknownQuantitiesVector.size() ; q++) {
		      
	  std::vector<double>  value(_UnknownQuantitiesVector[q]->_dim,0.);
            //the fact is that THERE ARE DIFFERENT DOF OBJECTS for DIFFERENT FE families
	    //for each family we should only pick the dof objects that are needed
	    //what changes between the FE families is the DOF OBJECT YOU PROVIDE: it could be a NODE or a CELL
	    // Notice that for some elements you don't have the midpoint of the element!
     if (_UnknownQuantitiesVector[q]->_FEord < KK) {
       
        for (uint ivar=0; ivar < _UnknownQuantitiesVector[q]->_dim; ivar++) {
       
          for (uint k=0; k < currelem.GetElemType(_UnknownQuantitiesVector[q]->_FEord)->GetNDofs() ; k++) {
	    
                const int fine_node = GetMLProb().GetMeshTwo()._el_map[VV][ k + ( iel + iel_b )*el_dof_objs ];
                for (uint idim = 0; idim < GetMLProb().GetMeshTwo().get_dim(); idim++) xp[idim] = GetMLProb().GetMeshTwo()._xyz[ fine_node + idim*coords_fine_offset ];

	  _UnknownQuantitiesVector[q]->initialize_xyz(&xp[0],value);
		    const int dof_pos_lev = _dofmap.GetDofQuantityComponent(Level,_UnknownQuantitiesVector[q],ivar,fine_node);
		    _x[Level]->set( dof_pos_lev, value[ivar] );
		      
	    }  //dof objects 
         }

     }
     else if (_UnknownQuantitiesVector[q]->_FEord == KK) { 

	    for (uint ivar=0; ivar < _UnknownQuantitiesVector[q]->_dim; ivar++) {
	    
                for (uint k=0; k < currelem.GetElemType(_UnknownQuantitiesVector[q]->_FEord)->GetNDofs() ; k++) { //only 1
		  
       int sum_elems_prev_sd_at_lev = 0;
	  for (uint pr = 0; pr < GetMLProb().GetMeshTwo()._iproc; pr++) { sum_elems_prev_sd_at_lev += GetMLProb().GetMeshTwo()._off_el[VV][pr*GetGridn() + Level + 1] - GetMLProb().GetMeshTwo()._off_el[VV][pr*GetGridn() + Level]; }
	  
          currelem.GetMidpoint();
	  
	  _UnknownQuantitiesVector[q]->initialize_xyz(&currelem.GetMidpoint()[0],value);

	      const int elem_lev = iel + sum_elems_prev_sd_at_lev;
	      const int dof_pos_lev = _dofmap.GetDofQuantityComponent(Level,_UnknownQuantitiesVector[q],ivar,elem_lev);
              _x[Level]->set( dof_pos_lev, value[ivar] );
	    
 	                       }  //k
	                  } //ivar
       
                     }     //end KK
	    
              } //qty
        
        
        } // end of element loop

        _x[Level]->localize(*_x_old[Level]);
        _x_old[Level]->close();
	
    } //end Level
    

#ifdef DEFAULT_PRINT_INFO
        std::cout << "\n Initialize(Base): Initial solution defined by ic_read" << "\n \n";
#endif

    return;
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

    _b[GetGridn()-1]->close();
    double bNorm_fine =     _b[GetGridn()-1]->l2_norm();
    _x_old[GetGridn()-1]->close();
    double x_old_fine = _x_old[GetGridn()-1]->l2_norm();

#ifdef DEFAULT_PRINT_INFO
    std::cout << " bNorm_fine l2 "     <<  bNorm_fine                     << std::endl;
    std::cout << " bNorm_fine linfty " << _b[GetGridn()-1]->linfty_norm()  << std::endl;
    std::cout << " xold_fine l2 "      <<  x_old_fine                     << std::endl;
#endif

    // FAS Multigrid (Nested) ---------
    bool NestedMG=false;
    if (NestedMG) {
        _x[0]->zero();
        MGStep(0,1.e-20,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post);

        //smooth on the coarse level WITH PHYSICAL b !
        //and compute the residual

        for (uint Level = 1; Level < GetGridn(); Level++) {

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

        res_fine = MGStep(GetGridn()-1,Eps1,MaxIter,Gamma,Nc_pre,Nc_coarse,Nc_post);

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

        rest = _LinSolver[Level]->solve(*_A[Level],*_A[Level],*_x[Level],*_b[Level],DEFAULT_EPS_LSOLV_C,Nc_coarse);  //****** smooth on the coarsest level

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

        rest = _LinSolver[Level]->solve(*_A[Level],*_A[Level],*_x[Level],*_b[Level],DEFAULT_EPS_PREPOST, Nc_pre); //****** smooth on the finer level

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

        rest = _LinSolver[Level]->solve(*_A[Level],*_A[Level],*_x[Level],*_b[Level],DEFAULT_EPS_PREPOST,Nc_post);  //***** smooth on the coarser level

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

// =========================================
void SystemTwo::ReadMGOps(const std::string output_path) {

    std::string     f_matrix = DEFAULT_F_MATRIX;
    std::string       f_rest = DEFAULT_F_REST;
    std::string       f_prol = DEFAULT_F_PROL;
    std::string       ext_h5 = DEFAULT_EXT_H5;

    std::ostringstream filename;
    std::string filename_base;
    filename_base = output_path + "/";
    
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

  _A.resize(GetGridn());

    for (uint Level = 0; Level< GetGridn(); Level++) {

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
            XDMFWriter::read_UIhdf5(file,dim_name.str().c_str(),rowcln[r][c]);

            //row length
            std::ostringstream len_name;
            len_name << groupname_lev.str() << "/" << "LEN" << fe_couple.str();
            length_row[r][c]=new uint[ rowcln[r][c][0]+1 ];
            XDMFWriter::read_UIhdf5(file,len_name.str().c_str(),length_row[r][c]);

            // matrix off diagonal
            std::ostringstream offlen_name;
            offlen_name << groupname_lev.str() << "/" << "OFFLEN" << fe_couple.str();
            length_offrow[r][c]=new uint[ rowcln[r][c][0]+1 ];
            XDMFWriter::read_UIhdf5(file,offlen_name.str().c_str(),length_offrow[r][c]);

            // matrix pos //must stay AFTER reading length_offrow
            std::ostringstream pos_name;
            pos_name << groupname_lev.str() << "/" << "POS" << fe_couple.str();
            pos_row[r][c]=new uint[ length_row[r][c][rowcln[r][c][0]] ];
            XDMFWriter::read_UIhdf5(file,pos_name.str().c_str(),pos_row[r][c]);

        } //end col
    } //end row
    
    H5Fclose(file);

//============================================================================
//============ compute things for the sparsity pattern =======================
//============================================================================
    
    int NoLevels  = GetMLProb().GetMeshTwo()._NoLevels;
    uint off_proc = NoLevels*GetMLProb().GetMeshTwo()._iproc;

    uint mrow_glob_t = 0;
    for (int fe=0; fe<QL; fe++) mrow_glob_t += _dofmap._nvars[fe]*rowcln[fe][fe][0];
    uint ncol_glob_t = mrow_glob_t;

    uint mrow_lev_proc_t  = 0;
    for (int fe=0; fe<QL; fe++)  mrow_lev_proc_t +=  _dofmap._DofLocLevProcFE[Level][GetMLProb().GetMeshTwo()._iproc][fe]*_dofmap._nvars[fe];
    uint ncol_lev_proc_t  = mrow_lev_proc_t;

    uint DofObjInit_lev_PrevProcs[QL];  //TODO what is this? it is the ROW INDEX at which to begin for every processor
                      
     for (int r=0; r<QL; r++)     DofObjInit_lev_PrevProcs[r] = 0;
         
    for (uint isubd=0; isubd<GetMLProb().GetMeshTwo()._iproc; isubd++) {
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
    // Now, remember that GetMLProb().GetMeshTwo()._off_nd[QQ] should only be used for computing offsets, so differences.
    // Here, it is used ALONE, because it gives you the NODE (in femus ordering) AT WHICH THE CURRENT SUBDOMAIN BEGINS,
    //and then from _node_dof[Level] (which was already constructed) you get THE LOCAL DOF AT THE CURRENT LEVEL TO START FROM.
    // Clearly, pay attention when you add elements, because in that case you would need to REDO the _node_dof map !!!
    //TODO: also, what happens if you have a system with ONLY ELEMENT BASED DOFS?!?
    
    
    int FELevel[QL];
    FELevel[QQ] = Level;
    FELevel[LL] = (Level+GetGridn())%(GetGridn()+1); //This is the map for the level of the LINEAR DOFS
    FELevel[KK] = Level;

    int off_onevar[QL];
    off_onevar[QQ] = GetMLProb().GetMeshTwo()._NoNodesXLev[GetGridn()-1];
    off_onevar[LL] = GetMLProb().GetMeshTwo()._NoNodesXLev[GetGridn()-1];
    off_onevar[KK] = GetMLProb().GetMeshTwo()._n_elements_vb_lev[VV][Level];
    
    uint  off_EachFEFromStart[QL];
    off_EachFEFromStart[QQ] = 0;
    off_EachFEFromStart[LL] = _dofmap._nvars[QQ]*off_onevar[QQ];
    off_EachFEFromStart[KK] = _dofmap._nvars[QQ]*off_onevar[QQ] + _dofmap._nvars[LL]*off_onevar[LL];

 //==============   
     for (int r=0;r<QL;r++) {
      for (uint ivar=0; ivar < _dofmap._nvars[r]; ivar++) {

        for (uint DofObj_lev = DofObjInit_lev_PrevProcs[r]; DofObj_lev < DofObjInit_lev_PrevProcs[r] + _dofmap._DofLocLevProcFE[Level][GetMLProb().GetMeshTwo()._iproc][r]; DofObj_lev++) {

            int dof_pos, irow;
	         if  (r<KK) {  dof_pos = GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[FELevel[r]][ DofObj_lev ];  }
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
// // //                     graph[irow][j+jvar*len[QQ]] = _node_dof[Level][ GetMLProb().GetMeshTwo()._node_map[FELevel[QQ]][pos_row[QQ][QQ][j+length_row[QQ][QQ][DofObj_lev]]]+jvar*GetMLProb().GetMeshTwo()._NoNodes[GetGridn()-1]];
// // //                 }
// // // 	    }
// // //                 // quadratic-linear 
// // //                 for (uint jvar=0; jvar<_nvars[LL]; jvar++) {
// // //                     for (int j=0; j<len[LL]; j++) {
// // //                         graph[irow][j+jvar*len[LL]+_nvars[QQ]*len[QQ]] = _node_dof[Level][GetMLProb().GetMeshTwo()._node_map[FELevel[LL]][pos_row[QQ][LL][j+length_row[QQ][LL][DofObj_lev]]]+(jvar+_nvars[QQ])*GetMLProb().GetMeshTwo()._NoNodes[GetGridn()-1]];
// // //                     }
// // //                 }
// // // 
// // //             
// // //             for (uint jvar=0; jvar<_nvars[QQ]; jvar++) {
// // //                 for (int j=0; j<len[QQ]; j++) {
// // //                     graph[irow][j+jvar*len[QQ]] = _node_dof[Level][GetMLProb().GetMeshTwo()._node_map[FELevel[QQ]][pos_row[LL][QQ][j+length_row[LL][QQ][DofObj_lev]]]+jvar*offset];
// // //                 }
// // //             }
// // //             
// // //             for (uint jvar=0; jvar<_nvars[LL]; jvar++) {
// // //                 for (int j=0; j<len[LL]; j++) {
// // //                     graph[irow][j+jvar*len[LL]+_nvars[QQ]*len[QQ]] = _node_dof[Level][ GetMLProb().GetMeshTwo()._node_map[FELevel[LL]][pos_row[LL][LL][j+length_row[LL][LL][DofObj_lev]]]+(jvar+_nvars[QQ])*offset];
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

    _Prl.resize(GetGridn());  //TODO one place is left empty in practice, we can optimize this!!!

    for (uint Level = 1; Level< GetGridn(); Level++) {
  
    uint Lev_c = Level-1;
    uint Lev_f = Level;    

        int FEXLevel_c[QL];
        FEXLevel_c[QQ] = Level-1;                                 //COARSE Level for QUADRATIC
        FEXLevel_c[LL] = (Level-1+GetGridn())%(GetGridn()+1);
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
   
    uint off_proc = GetMLProb().GetMeshTwo()._iproc*GetGridn();

    _Prl[ Lev_f ] = SparseMatrix::build().release();
// // //     _Prl[ Lev_f ]->init(0,0,0,0); //TODO BACK TO A REASONABLE INIT

    // local matrix dimension
    uint ml[QL]; uint nl[QL];
     for (int fe=0; fe<QL; fe++) { 
       if (fe < KK) {
       ml[fe] = GetMLProb().GetMeshTwo()._off_nd[fe][off_proc + Lev_f +1] - GetMLProb().GetMeshTwo()._off_nd[fe][off_proc];    //  local quadratic    //COARSE (rows)
       nl[fe] = GetMLProb().GetMeshTwo()._off_nd[fe][off_proc + Lev_c +1] - GetMLProb().GetMeshTwo()._off_nd[fe][off_proc];    // global quadratic  //FINE QUADRATIC (cols)
       }
       else if (fe == KK) { 
       ml[fe] = GetMLProb().GetMeshTwo()._off_el[VV][off_proc + Lev_f +1] - GetMLProb().GetMeshTwo()._off_el[VV][off_proc + Lev_f];
       nl[fe] = GetMLProb().GetMeshTwo()._off_el[VV][off_proc + Lev_c +1] - GetMLProb().GetMeshTwo()._off_el[VV][off_proc + Lev_c];
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
        for (uint isubd=0;isubd<GetMLProb().GetMeshTwo()._iproc; isubd++) {
       if (fe < KK)       ml_init[fe] += GetMLProb().GetMeshTwo()._off_nd[fe][isubd*GetGridn() + Lev_f +1] - GetMLProb().GetMeshTwo()._off_nd[fe][isubd*GetGridn()];
       else if (fe == KK) ml_init[fe] += GetMLProb().GetMeshTwo()._off_el[VV][isubd*GetGridn() + Lev_f +1] - GetMLProb().GetMeshTwo()._off_el[VV][isubd*GetGridn() + Lev_f];
	}
      }
  
    //============= POSITION =========
   for (int fe=0;fe<QL;fe++) {

     for (uint ivar=0;ivar<_dofmap._nvars[fe];ivar++) {
        for (unsigned int i = ml_init[fe]; i < ml_init[fe]+ml[fe]; i++) {
          int dof_pos_f;
          if (fe < KK)         dof_pos_f = GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_f[fe] ][ i ];  //end fe < ql
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
	      if      (fe  < KK) dof_pos_c = GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_c[fe] ][ dof_pos_lev_c ];
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
          if (fe < KK)        dof_pos_f = GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_f[fe] ][ i ];  //end fe < ql
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
//     if (GetMLProb().GetMeshTwo()._iproc==0) _Prl[  Lev_f ]->print_personal();
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
//Level goes from 0 to < GetGridn() - 1 ==> Level is COARSE here
  
//    uint Lev_c = Level;
//   uint Lev_f = Level+1;
	//with these you explore arrays that go from 0  to GetGridn() - 1
                           //so where the distinction between QQ and LL is already made
                           // with the EXTENDED levels you explore things that have an additional level,
                           // and so can work both with QQ and with LL
    //the point is: there are parts where you cannot use extended levels, and parts where you can
    //for instance, in this routine the FINE LEVELS and the FINE EXTENDED LEVELS will both be ok,
    //so we can use them in both cases, but we cannot say the same for the COARSE levels and ext levels
    
    //AAA fai molta attenzione: per esplorare la node_dof devi usare Lev_c e Lev_f,
    //perche' sono legati ai DOF (devi pensare che la questione del mesh e' gia' risolta)
void SystemTwo::ReadRest(const std::string& name) {
 
  _Rst.resize(GetGridn());  //TODO why do it bigger?
  
  for (uint Level = 0; Level< GetGridn() - 1; Level++) {
    
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
// while the FINE LEVELS never risk to be "extended", the maximum for them is (GetGridn()-1) !
    
        int FEXLevel_c[QL];
        FEXLevel_c[QQ] = Level;                                 //COARSE Level for QUADRATIC
        FEXLevel_c[LL] = (GetGridn() + Level)%(GetGridn() + 1);   //COARSE Level for LINEAR: Level1=0 means coarse linear, a finer linear is the first coarse quadratic, and so on and so on
        FEXLevel_c[KK] = Level;                                 //COARSE Level for CONSTANT
        int FEXLevel_f[QL];
        FEXLevel_f[QQ] = Level+1;                                  //FINE Level for QUADRATIC
        FEXLevel_f[LL] = Level;                                // AAA look at the symmetry, this is exactly (_n_levels + Level1 + 1)%(_n_levels + 1); ! //FINE Level for LINEAR:   Level1=0 means coarse linear, a finer linear is the first coarse quadratic, and so on and so on
        FEXLevel_f[KK] = Level+1;                                  //FINE Level for CONSTANT

    uint off_proc=GetGridn()*GetMLProb().GetMeshTwo()._iproc;

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
        pattern._ml += _dofmap._nvars[fe]*_dofmap._DofLocLevProcFE[Lev_c][GetMLProb().GetMeshTwo()._iproc][fe]; 
        pattern._nl += _dofmap._nvars[fe]*_dofmap._DofLocLevProcFE[Lev_f][GetMLProb().GetMeshTwo()._iproc][fe];
     } 
     
    // starting indices for local matrix
   uint ml_start = _dofmap.GetStartDof(Lev_c,off_proc);    //  offset proc nodes
   pattern._ml_start = ml_start;
   uint DofObjInit_lev_PrevProcs_c[QL];
        for (int fe=0;fe<QL;fe++) { DofObjInit_lev_PrevProcs_c[fe] = 0;  }
        for (int fe=0;fe<QL;fe++) { 
    for (uint isubd=0;isubd<GetMLProb().GetMeshTwo()._iproc; isubd++) { //up to the current processor
       if (fe < KK)       /*mlinit*/ DofObjInit_lev_PrevProcs_c[fe] += GetMLProb().GetMeshTwo()._off_nd[fe][ isubd*GetGridn() + Lev_c +1 ]  - GetMLProb().GetMeshTwo()._off_nd[fe][isubd*GetGridn()];  
       else if (fe == KK) /*mlinit*/ DofObjInit_lev_PrevProcs_c[fe] += GetMLProb().GetMeshTwo()._off_el[VV][ isubd*GetGridn() + Lev_c +1 ]  - GetMLProb().GetMeshTwo()._off_el[VV][isubd*GetGridn() + Lev_c];   
         }
     }

    //============= POSITION =========
        for (int fe=0; fe<QL; fe++) { 
    for (uint ivar=0;ivar<_dofmap._nvars[fe];ivar++) {
        for (unsigned int i = DofObjInit_lev_PrevProcs_c[fe]; i< DofObjInit_lev_PrevProcs_c[fe] + _dofmap._DofLocLevProcFE[Lev_c][GetMLProb().GetMeshTwo()._iproc][fe]; i++) {
	  int dof_pos_c;
          if (fe < KK)         dof_pos_c = GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_c[fe] ][ i ];
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
    
	      if      (fe  < KK)  dof_pos_f = GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_f[fe] ][ dof_pos_lev_f ];
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
//     if (GetMLProb().GetMeshTwo()._iproc==0) _Rst[Lev_c]->print_personal(); //there is no print function for rectangular matrices, and print_personal doesnt seem to be working...
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
        for (unsigned int i = DofObjInit_lev_PrevProcs_c[fe]; i< DofObjInit_lev_PrevProcs_c[fe] + _dofmap._DofLocLevProcFE[Lev_c][GetMLProb().GetMeshTwo()._iproc][fe]; i++) {
	  int dof_pos_c;
          if (fe < KK)         dof_pos_c = GetMLProb().GetMeshTwo()._Qnode_lev_Qnode_fine[ FEXLevel_c[fe] ][ i ];
          else if (fe == KK)   dof_pos_c = i;
            int irow     = _dofmap.GetDof(Lev_c,fe,ivar,dof_pos_c);
            int irow_top = _dofmap.GetDof(GetGridn()-1,fe,ivar,dof_pos_c);
            uint ncol = len[fe][i+1]-len[fe][i];
            tmp[0]=irow;
            std::vector< uint> ind(pattern[irow].size()-1);
// 	    std::cout << "\n ==== " << irow << ": ";
            for (uint i1=0;i1<ind.size();i1++) { ind[i1] = pattern[irow][i1]; /*std::cout << " " << ind[i1] << " ";*/}
            valmat = new DenseMatrix(1,ncol);  //TODO add a matrix row by row...
            for (uint j=0; j<ncol; j++) (*valmat)(0,j) = _bcond._bc[irow_top]*Rest_val[fe][ j+len[fe][i] ];
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
//     if (GetMLProb().GetMeshTwo()._iproc==0)  _Rst[Lev_c]->print_personal(std::cout);
//     _Rst[Lev_c]->print_graphic(false); // TODO should pass this true or false as a parameter

  } //end levels
  
#ifdef DEFAULT_PRINT_INFO
    std::cout << " ReadRest(B): read Op " << name.c_str() << std::endl;
#endif
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
