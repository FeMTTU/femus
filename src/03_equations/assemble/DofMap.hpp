/*=========================================================================

 Program: FEMUS
 Module: DofMap
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_DofMap_hpp__
#define __femus_equations_DofMap_hpp__


#include <string>



#include "FETypeEnum.hpp"
#include "VBTypeEnum.hpp"
#include "Typedefs.hpp"
#include "MultiLevelMeshTwo.hpp"
// #include "SystemTwo.hpp"  //TODO opening this creates a conflict, had to move one inline function to the .cpp

namespace femus {


  class SystemTwo;
  class Quantity;
  
//=======================================================================
//==== DOF MAP of the equation ============ (procs,levels) ==============   //// LinearEquation (each level)
//=======================================================================
  
class DofMap {


public:

     DofMap(const SystemTwo * eqn_in, const MultiLevelMeshTwo& mesh);
    ~DofMap();

    const SystemTwo * _eqn;
    const MultiLevelMeshTwo & _mesh;
    
//====== data =======
  uint  *    _Dim;            //number of dofs per level
  uint **    _DofNumLevFE;     
  uint **    _DofOffLevFE;     
  uint ***   _DofLocLevProcFE;     
  uint       _nvars[QL];      ///  number of SCALAR variables          
  uint       _VarOff[QL];
  uint       _n_vars;           ///< number of SCALAR variables
  
//====== functions =======
          void initNVars();
          void ComputeMeshToDof();
          void PrintMeshToDof() const;

  inline  int GetDof(const uint Level,const uint fe,const uint ivar,const uint i) const;
  /*inline*/  int GetDofQuantityComponent(const uint Level, const Quantity* quantity_in, const uint quantity_var,const uint dofobj) const;  //TODO this should be in the .hpp instead
  inline  int GetDofPosIn(const uint Level,const uint pos_in) const;
  inline  int GetStartDof(const uint Level, const uint offproc) const;  //TODO understand this, possibly remove it

  
  private:
    
    int    **  _node_dof; ///< dof map
  

};




  int DofMap::GetDofPosIn(const uint Level,const uint pos_in) const {
         return _node_dof[Level][pos_in];
    }

    /** at level Level, for a given FE, take the dof of the ivar_fe variable of that FE*/
  int DofMap::GetDof(const uint Level,const uint fe,const uint ivar_fe,const uint dofobj) const {
         return _node_dof[Level][ dofobj + ivar_fe*_DofNumLevFE[Level][fe] + _DofOffLevFE[Level][fe] ];
    }
    
  int DofMap::GetStartDof(const uint Level, const uint offproc) const {
          return _node_dof[Level][  _mesh._off_nd[QQ][offproc] ];
       }
 

} //end namespace femus





#endif