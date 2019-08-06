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

#include "mpi.h"

#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "SparseMatrix.hpp"
#include "ElemType.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"



namespace femus {


  const unsigned Writer::FemusToVTKorToXDMFConn[27] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 23, 21, 20, 22, 24, 25, 26};

  Writer::Writer (MultiLevelSolution* ml_sol) :
    _ml_sol (ml_sol), _ml_mesh (ml_sol->_mlMesh) {
    _gridn = _ml_mesh->GetNumberOfLevels();
    _moving_mesh = 0;
    _graph = false;
    _surface = false;
  }

  Writer::Writer (MultiLevelMesh* ml_mesh) :
    _ml_sol (NULL), _ml_mesh (ml_mesh) {
    _gridn = _ml_mesh->GetNumberOfLevels();
    _moving_mesh = 0;
    _graph = false;
    _surface = false;
  }

  Writer::~Writer() { }


  std::unique_ptr<Writer> Writer::build (const WriterEnum format, MultiLevelSolution * ml_sol)  {

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

  std::unique_ptr<Writer> Writer::build (const WriterEnum format, MultiLevelMesh * ml_mesh)  {

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

  void Writer::SetMovingMesh (std::vector<std::string>& movvars_in) {
    _moving_mesh = 1;
    _moving_vars = movvars_in;
  }

  void Writer::SetGraphVariable (const std::string &graphVaraible) {
    _graph = true;
    _surface = false;
    _graphVariable = graphVaraible;
  }

  void Writer::SetSurfaceVariables (std::vector < std::string > &surfaceVariable) {
    _surface = true;
    _graph = false;
    _surfaceVariables = surfaceVariable;
  }


} //end namespace femus


