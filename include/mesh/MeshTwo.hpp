#ifndef __mesh123D_h
#define __mesh123D_h

//C++ 
#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "hdf5.h"

#include "Typedefs.hpp"
#include "FEMTTUConfig.h"
#include "RunTimeMap.hpp"
#include "GeomEl.hpp"
class Files;
class QuantityLocal;
class Domain;

class Mesh  {

public:
  
    Files& _files; 
    RunTimeMap<double>  _mesh_rtmap;
    uint _dim;               ///< spatial dimension
    uint _meshVB;          /// Number of FEM "manifold" Families (VB)  F

// ===== ABSTRACT GEOMEL(S) =====
    GeomEl    _GeomEl;
    const uint _n_GeomEl;         //number of geom elem types
    uint*      _type_FEM;         //just for check

// ====== DOMAIN SHAPE (TODO optional => pointer) ----- //if I put it as reference I'd have to initialize it
    Domain* _domain;
    Domain* GetDomain();
    void TransformElemNodesToRef(const uint vb,const double* xx_qnds,double* refbox_xyz) const;

// ==== PARALLEL ===
    uint _iproc;       /// current subdomain
    uint _NoSubdom;    /// Number of subdomains (subdomain P)

// ==== MULTIGRID ====
    uint _NoLevels;          /// Number of Levels (L)

// ==== ELEMENTS ====
    uint** _NoElements;        /// Number of Elements        [VB][L]                                     
    int**  _off_el;            /// offset of element numbers [VB][L+P*_NoLevels]         
    uint** _el_map;            /// Connectivity map         [VB][L+P*_NoLevels]
    int**  _el_bdry_to_vol;

// ===== NODES ======
    uint*   _NoNodesXLev;            ///VOLUME, at the EXTENDED LEVELS
    double* _xyz;                    ///< node coordinates
    uint**  _Qnode_lev_Qnode_fine;   ///ONLY USED FOR THE EQUATION /// FROM DOF TO NODE: Q or L dofs by using an add aux level [_NoLevels+1][NoNodes[L]] //ONLY VOLUME  //NO PROC info, these are DOFS in the "SERIALIZED NUMBERING"
    uint**  _Qnode_fine_Qnode_lev;   ///Un nodo fine puo' o non appartenere ad un certo livello... quindi devo continuare a mantenere i -1, altrimenti dovrei fare un search che fa perdere tempo
    int**   _off_nd;                 ///ONLY USED FOR THE EQUATION /// FROM DOF TO NODE: offsets for dofs  [QL][L+P*_NoLevels]                     //ONLY VOLUME
                                     // on the fine node numbering, the nodes corresponding to linear dofs are numbered FIRST... or not?

//===== Constructors/ Destructor ===========
     Mesh (Files& files_in, RunTimeMap<double>& map_in, const double Lref, Domain* domain_in=NULL);
    ~Mesh ();			  ///<  Level Mesh Destructor
    void clear ();		  ///<  substructure Destructor

    //======= Print/read functions =======
    void PrintForVisualizationAllLEVAllVB() const;
    void PrintSubdomFlagOnLinCells(std::string filename) const;
    void PrintXDMFAllLEVAllVB() const;
    void PrintConnLinAllLEVAllVB() const;
    void PrintXDMFGridVB(std::ofstream& out, std::ostringstream& top_file, std::ostringstream& geom_file,const uint Level, const uint vb) const;
    void PrintConnLinVB(hid_t file, const uint Level, const uint vb) const; 
    void PrintMeshFile(const std::string & namefile) const;
    void ReadMeshFile();
    
protected:

    const double _Lref;          ///Reference length for non-dimensionalization

};

#endif
