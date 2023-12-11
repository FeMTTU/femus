/*=========================================================================

 Program: FEMUS
 Module: Writer
 Authors: Eugenio Aulisa, Simone Bnà, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

#include "Writer.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"
#include "MultiLevelSolution.hpp"
#include "ElemType.hpp"
#include "FElemTypeEnum_list.hpp"




namespace femus {


  Writer::Writer (const MultiLevelSolution* ml_sol) : 
  Writer_one_level(),
    _ml_sol (ml_sol),
    _ml_mesh (ml_sol->GetMLMesh())
{
      
    _gridn = _ml_mesh->GetNumberOfLevels();
        
  }
  
  

  Writer::Writer (const MultiLevelMesh* ml_mesh) :
  Writer_one_level(),
    _ml_sol (NULL), 
    _ml_mesh (ml_mesh) {
      
    _gridn = _ml_mesh->GetNumberOfLevels();
    
      
  }

  


  
  std::unique_ptr<Writer> Writer::build (const WriterEnum format, const MultiLevelSolution * ml_sol)  {

    switch (format) {
      case VTK: {
        std::unique_ptr<Writer>   ap (new VTKWriter (ml_sol));
        return ap;
      }
      case GMV: {
        std::unique_ptr<Writer>   ap (new GMVWriter (ml_sol));
        return ap;
      }
#ifdef HAVE_HDF5
      case XDMF: {
        std::unique_ptr<Writer>   ap (new XDMFWriter (ml_sol));
        return ap;
      }
#endif
      default: {
        std::cout << "Format not supported" << std::endl;
        abort();
      }
    } //end switch

  }

  std::unique_ptr<Writer> Writer::build (const WriterEnum format, const MultiLevelMesh * ml_mesh)  {

    switch (format) {
      case VTK: {
        std::unique_ptr<Writer>   ap (new VTKWriter (ml_mesh));
        return ap;
      }
      case GMV: {
        std::unique_ptr<Writer>   ap (new GMVWriter (ml_mesh));
        return ap;
      }
#ifdef HAVE_HDF5
      case XDMF: {
        std::unique_ptr<Writer>   ap (new XDMFWriter (ml_mesh));
        return ap;
      }
#endif
      default: {
        std::cout << "Format not supported" << std::endl;
        abort();
      }
    } //end switch

  }


    const Solution * Writer::get_solution(const unsigned level_in) const {
      return  ( _ml_sol != NULL ) ? _ml_sol->GetSolutionLevel( level_in - 1 )  :  NULL;
    }  
  
  
} //end namespace femus


