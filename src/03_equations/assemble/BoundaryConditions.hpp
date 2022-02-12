/*=========================================================================

 Program: FEMUS
 Module: BoundaryConditions
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_BoundaryConditions_hpp__
#define __femus_equations_BoundaryConditions_hpp__

#include "Typedefs.hpp"



namespace femus {

 class NumericVector;
 class DofMap;
 
class BoundaryConditions {


public:

     BoundaryConditions(const DofMap* dofmap_in);
    ~BoundaryConditions();

//=======================================================================
//==== BOUNDARY CONDITIONS of the equation ========= (procs,levels) ==     //MultilevelSolution
//=======================================================================
    int   *_bc;         //==== NODAL DIRICHLET ======== ///< boundary conditions map (top level)  // POINTWISE(NODAL) FLAG for the BOUNDARY DOFS = FLAG for the tEST FUNCTIONS //TODO this should be PrintNumericVector of the equation, integer instead of double! do it when you make it parallel especially! //Later on I will do a bc for every level, considering the ELEMENT DOFS
//     int  **_bc_fe_kk;   //==== CELL DIRICHLET ========
          void    GenerateBdc();          //MultilevelSolution
//========= treating NumericVectors, related to Dirichlet Boundary Conditions! =======
          void Bc_ScaleDofVec(NumericVector * myvec,  double ScaleFac);
          void Bc_AddDofVec(NumericVector* myvec, NumericVector* myvec2 );
          void Bc_AddScaleDofVec(NumericVector* vec_in,NumericVector* vec_out,const double ScaleFac );
	  
	  

 //====PENALTY DIRICHLET ======Elem BC=====================
//    uint  _Dir_pen_fl;         ///flag for penalty with Dirichlet (0=no penalty, 1=yes penalty) //this penalty is for ALL the QUADRATIC variables
//  static const int _number_tang_comps[3];  //  {0,1,3} Number of tangential components in 1D,2D,3D (see BC in general domains) : use it only for 2D and 3D
//     int  ***_elem_bc;        ///[LEVELS][IPROC][2xELEMENTSxLEV&PROC]
//  double  ***_elem_val_norm;  ///[LEVELS][IPROC][1xELEMENTSxLEV&PROC]
//  double  ***_elem_val_tg;    ///[LEVELS][IPROC][(1or3)xELEMENTSxLEV&PROC]
//           void    GenerateBdcElem();
//           void elem_bc_read(const double * el_xm, int& surf_id, double * value,int * el_flag) const;
//           void  clearElBc();
//           void Bc_GetElFlagValLevSubd(const uint Level,const uint isubd,const uint iel,int* el_flag,double* el_value ) const;
	  
    
  private:
    
    const DofMap*  _dofmap;
  

};



} //end namespace femus





#endif