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

#ifndef __dofmap_hpp__
#define __dofmap_hpp__


#include <string>



#include "FETypeEnum.hpp"
#include "VBTypeEnum.hpp"
#include "Typedefs.hpp"
#include "MultiLevelMeshTwo.hpp"


namespace femus {


  class SystemTwo;
  
//=======================================================================
//==== DOF MAP of the equation ============ (procs,levels) ==============   //// LinearEquation (each level)
//=======================================================================
class DofMap {


public:

     DofMap(const SystemTwo& eqn_in, const MultiLevelMeshTwo& mesh);
    ~DofMap();

    const SystemTwo & _eqn;
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
  inline  int GetStartDof(const uint Level, const uint offproc) const;  //TODO understand this, possibly remove it

  
  private:
    
    int    **  _node_dof; ///< dof map
  

};


    int DofMap::GetDof(const uint Level,const uint fe,const uint ivar,const uint dofobj) const {
         return _node_dof[Level][ dofobj + ivar*_DofNumLevFE[Level][fe] + _DofOffLevFE[Level][fe] ];
    }
    
     inline  int DofMap::GetStartDof(const uint Level, const uint offproc) const {
          return _node_dof[Level][  _mesh._off_nd[QQ][offproc] ];
       }
 

} //end namespace femus





#endif