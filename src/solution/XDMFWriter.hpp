/*=========================================================================

 Program: FEMUS
 Module: XDMFWriter
 Authors: Eugenio Aulisa, Simone Bnà, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __xdmfwriter_h_
#define __xdmfwriter_h_

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Writer.hpp"
#include "MultiLevelMeshTwo.hpp"

namespace femus {



class XDMFWriter : public Writer {

public:

    /** Constructor. */
    XDMFWriter(MultiLevelSolution& ml_sol);
    
    /** Destructor */
    virtual ~XDMFWriter();

    /** write output function */
    virtual void write_system_solutions(const char order[], std::vector<std::string>& vars, const unsigned time_step = 0);

    /** write a wrapper file for paraview to open all the files of an history toghether */
    void write_solution_wrapper(const char type[]) const;

    /** */
    static void PrintXDMFTopologyGeometry(std::ofstream& out,const unsigned Level, const unsigned vb, const MultiLevelMeshTwo& mesh);

};


} //end namespace femus



#endif
