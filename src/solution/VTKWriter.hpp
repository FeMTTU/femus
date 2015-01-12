/*=========================================================================

 Program: FEMUS
 Module: VTKWriter
 Authors: Eugenio Aulisa, Simone Bnà

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __vtkwriter_h_
#define __vtkwriter_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Writer.hpp"


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class MultiLevelProblem;


class VTKWriter : public Writer {

public:

    /** Constructor. */
    VTKWriter(MultiLevelSolution& ml_sol);

    /** Destructor */
    virtual ~VTKWriter();

    /** write output function */
    virtual void write_system_solutions(const char order[], std::vector<std::string>& vars, const unsigned time_step=0);
    
private:
  
    /** femus to vtk cell type map */
    static short unsigned int femusToVtkCellType[3][6];

};

/**    
 * VTK Cell type: from http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
 * 
 * VTK_EMPTY_CELL                         = 0, 
 * VTK_VERTEX                             = 1, 
 * VTK_POLY_VERTEX                        = 2, 
 * VTK_LINE                               = 3,
 * VTK_POLY_LINE                          = 4, 
 * VTK_TRIANGLE                           = 5, 
 * VTK_TRIANGLE_STRIP                     = 6, 
 * VTK_POLYGON                            = 7,
 * VTK_PIXEL                              = 8, 
 * VTK_QUAD                               = 9, 
 * VTK_TETRA                              = 10, 
 * VTK_VOXEL                              = 11,
 * VTK_HEXAHEDRON                         = 12, 
 * VTK_WEDGE                              = 13, 
 * VTK_PYRAMID                            = 14, 
 * VTK_PENTAGONAL_PRISM                   = 15,
 * VTK_HEXAGONAL_PRISM                    = 16, 
 * VTK_QUADRATIC_EDGE                     = 21, 
 * VTK_QUADRATIC_TRIANGLE                 = 22, 
 * VTK_QUADRATIC_QUAD                     = 23,
 * VTK_QUADRATIC_POLYGON                  = 36, 
 * VTK_QUADRATIC_TETRA                    = 24, 
 * VTK_QUADRATIC_HEXAHEDRON               = 25, 
 * VTK_QUADRATIC_WEDGE                    = 26,
 * VTK_QUADRATIC_PYRAMID                  = 27, 
 * VTK_BIQUADRATIC_QUAD                   = 28, 
 * VTK_TRIQUADRATIC_HEXAHEDRON            = 29, 
 * VTK_QUADRATIC_LINEAR_QUAD              = 30,
 * VTK_QUADRATIC_LINEAR_WEDGE             = 31, 
 * VTK_BIQUADRATIC_QUADRATIC_WEDGE        = 32,
 * VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON   = 33,
 * VTK_BIQUADRATIC_TRIANGLE               = 34,
 * VTK_CUBIC_LINE                         = 35, 
 * VTK_CONVEX_POINT_SET                   = 41, 
 * VTK_POLYHEDRON                         = 42, 
 * VTK_PARAMETRIC_CURVE                   = 51,
 * VTK_PARAMETRIC_SURFACE                 = 52, 
 * VTK_PARAMETRIC_TRI_SURFACE             = 53, 
 * VTK_PARAMETRIC_QUAD_SURFACE            = 54, 
 * VTK_PARAMETRIC_TETRA_REGION            = 55,
 * VTK_PARAMETRIC_HEX_REGION              = 56, 
 * VTK_HIGHER_ORDER_EDGE                  = 60, 
 * VTK_HIGHER_ORDER_TRIANGLE              = 61, 
 * VTK_HIGHER_ORDER_QUAD                  = 62,
 * VTK_HIGHER_ORDER_POLYGON               = 63, 
 * VTK_HIGHER_ORDER_TETRAHEDRON           = 64, 
 * VTK_HIGHER_ORDER_WEDGE                 = 65, 
 * VTK_HIGHER_ORDER_PYRAMID               = 66,
 * VTK_HIGHER_ORDER_HEXAHEDRON            = 67, 
 */
  
 short unsigned int VTKWriter::femusToVtkCellType[3][6]= {{12,10,13,9,5,3},{25,24,26,23,22,21},{29,24,32,28,22,21}};


} //end namespace femus



#endif
