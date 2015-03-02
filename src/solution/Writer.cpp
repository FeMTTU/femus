/*=========================================================================

 Program: FEMUS
 Module: Writer
 Authors: Eugenio Aulisa, Simone Bnà
 
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
#include "SparseMatrix.hpp"
#include "ElemType.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"



namespace femus {


  Writer::Writer(MultiLevelSolution& ml_sol):
    _ml_sol(ml_sol)
  {
    _gridn = ml_sol._ml_msh->GetNumberOfLevels();
    _gridr = ml_sol._ml_msh->GetNumberOfGridTotallyRefined();
    _moving_mesh = 0;
  }

  Writer::~Writer() { }


  std::auto_ptr<Writer> Writer::build(const WriterEnum format, MultiLevelSolution * ml_sol)  {
     
    switch (format) {
      case VTK: {
	std::auto_ptr<Writer>   ap(new VTKWriter(*ml_sol)); return ap;
      }
      case GMV: {
	std::auto_ptr<Writer>   ap(new GMVWriter(*ml_sol)); return ap;
      }
      case XDMF: {
	std::auto_ptr<Writer>   ap(new XDMFWriter(*ml_sol)); return ap;
      }
      default: {
	std::cout << "Format not supported" << std::endl; 
	abort(); 
      }
    } //end switch
  }


  void Writer::SetMovingMesh(std::vector<std::string>& movvars_in){
    _moving_mesh = 1;
    _moving_vars = movvars_in;
  }

  
} //end namespace femus


