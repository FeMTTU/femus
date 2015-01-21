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



namespace femus {


  class EqnBase;
  class MeshTwo;
  
//=======================================================================
//==== DOF MAP of the equation ============ (procs,levels) ==============   //// LinearEquation (each level)
//=======================================================================
class DofMap {


public:

     DofMap(const EqnBase& eqn_in, const MeshTwo& mesh);
    ~DofMap();

    const EqnBase & _eqn;
    const MeshTwo & _mesh;
    
//====== data =======
  int    **  _node_dof; ///< dof map
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







};




} //end namespace femus





#endif