/*=========================================================================

 Program: FEMUS
 Module: MeshMetisPartitioning
 Authors: Simone Bnà, Eugenio Aulisa

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MeshMetisPartitioning.hpp"
#include "Mesh.hpp"
#include "FemusConfig.hpp"

#ifdef HAVE_METIS
#include "metis.h"
#endif

//C++ include
#include <iostream>


namespace femus {



  MeshMetisPartitioning::MeshMetisPartitioning(const Mesh& mesh) : MeshPartitioning(mesh) {

  }


//------------------------------------------------------------------------------------------------------
  void MeshMetisPartitioning::DoPartition(std::vector <unsigned>& partition, const bool& AMR) {

    int nnodes = _mesh.GetNumberOfNodes();
    int nelem = _mesh.GetNumberOfElements();

    std::vector <int> epart(nelem);

    if(_nprocs == 1) {
      //serial computation
      for(unsigned i = 0; i < nelem; i++) {
        epart[i] = 0;
      }
    }
    else if(_nprocs > nelem) {
      std::cout << "Error In MeshMetis::DoPartition, the number of processes " << _nprocs
                << " is greater than the number of elements " << nelem << std::endl;
      abort();
    }
    else if(_nprocs == nelem) {
      for(unsigned i = 0; i < nelem; i++) {
        epart[i] = i;
      }
    }
    else {

#ifndef HAVE_METIS
      std::cerr << "Fatal error: Metis library was not found. Metis partioning algorithm cannot be called!" << std::endl;
      exit(1);
#endif

      unsigned eind_size =   _mesh.el->GetElementNumber(geom_elems[HEX]) * _mesh.el->GetNVE(0, CONTINUOUS_BIQUADRATIC)
                           + _mesh.el->GetElementNumber(geom_elems[TET]) * _mesh.el->GetNVE(1, CONTINUOUS_BIQUADRATIC)
                           + _mesh.el->GetElementNumber(geom_elems[WEDGE]) *  _mesh.el->GetNVE(2, CONTINUOUS_BIQUADRATIC)
                           + _mesh.el->GetElementNumber(geom_elems[QUAD]) * _mesh.el->GetNVE(3, CONTINUOUS_BIQUADRATIC)
                           + _mesh.el->GetElementNumber(geom_elems[TRI]) * _mesh.el->GetNVE(4, CONTINUOUS_BIQUADRATIC)
                           + _mesh.el->GetElementNumber(geom_elems[LINE]) * _mesh.el->GetNVE(5, CONTINUOUS_BIQUADRATIC);

      std::vector < idx_t > eptr(nelem + 1);
      std::vector < idx_t > eind(eind_size);

      std::vector < int > npart(nnodes);

      idx_t objval;
      idx_t options[METIS_NOPTIONS];

      METIS_SetDefaultOptions(options);

      options[METIS_OPTION_NUMBERING] = 0;
      options[METIS_OPTION_DBGLVL]   = 0;
      options[METIS_OPTION_CTYPE]    = METIS_CTYPE_SHEM;
      options[METIS_OPTION_PTYPE]    = METIS_PTYPE_KWAY;
      options[METIS_OPTION_IPTYPE]   = METIS_IPTYPE_RANDOM;
      options[METIS_OPTION_CONTIG]   = 0;
      options[METIS_OPTION_MINCONN]  = 1;
      options[METIS_OPTION_NITER]    = 10;
      options[METIS_OPTION_UFACTOR]  = 100;

      eptr[0] = 0;
      unsigned counter = 0;
      for(unsigned iel = 0; iel < nelem; iel++) {
        unsigned ndofs = _mesh.el->GetElementDofNumber(iel, 2);
        eptr[iel + 1] = eptr[iel] + ndofs;
        for(unsigned inode = 0; inode < ndofs; inode++) {
          eind[counter] = _mesh.el->GetElementDofIndex(iel, inode);
          counter++;
        }
      }


      int ncommon = (AMR || _mesh.GetDimension() == 1) ? 1 : _mesh.GetDimension() + 1;

      //I call the Mesh partioning function of Metis library (output is epart(own elem) and npart (own nodes))
      int err = METIS_PartMeshDual(&nelem, &nnodes, &eptr[0], &eind[0], NULL, NULL, &ncommon, &_nprocs, NULL, options, &objval, &epart[0], &npart[0]);

      if(err == METIS_OK) {
        std::cout << " METIS PARTITIONING IS OK " << std::endl;
      }
      else if(err == METIS_ERROR_INPUT) {
        std::cout << " METIS_ERROR_INPUT " << std::endl;
        exit(1);
      }
      else if(err == METIS_ERROR_MEMORY) {
        std::cout << " METIS_ERROR_MEMORY " << std::endl;
        exit(2);
      }
      else {
        std::cout << " METIS_GENERIC_ERROR " << std::endl;
        exit(3);
      }



    }

    partition.resize(epart.size());
    for(unsigned i = 0; i < epart.size(); i++) {
      partition[i] = epart[i];
    }

    return;
  }

  void MeshMetisPartitioning::DoPartition(std::vector <unsigned>& partition, const Mesh& meshc) {
      
    partition.resize(_mesh.GetNumberOfElements());
    const unsigned refIndex = _mesh.GetMeshElements()->GetRefIndex( _mesh.GetDimension() );
    for(int isdom = 0; isdom < _nprocs; isdom++) {
      for(unsigned iel = meshc.GetElementOffset(isdom); iel < meshc.GetElementOffset(isdom + 1); iel++) {
        for(unsigned j = 0; j < refIndex; j++) {
          partition[ iel * refIndex + j ] = isdom;
        }
      }
    }
    
  }
  

}
