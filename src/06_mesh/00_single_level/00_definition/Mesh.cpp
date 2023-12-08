/*=========================================================================

 Program: FEMUS
 Module: Mesh
 Authors: Eugenio Aulisa

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "Mesh.hpp"
#include "MeshGeneration.hpp"
#include "GambitIO.hpp"
#include "MED_IO.hpp"
#include "MeshMetisPartitioning.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"

// C++ includes
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>


namespace femus {

  
bool (* Mesh::_SetRefinementFlag)(const std::vector < double >& x, const int &ElemGroupNumber,const int &level) = NULL; 


  bool Mesh::_IsUserRefinementFunctionDefined = false;

  unsigned Mesh::_dimension = 2;                                    ///@todo I don't like the default dimension to be 2
  unsigned Mesh::_ref_index = 4; // 8*DIM[2]+4*DIM[1]+2*DIM[0];     ///@todo I don't like the default dimension to be 2
  unsigned Mesh::_ref_face_index = 2; // 4*DIM[2]+2*DIM[1]+1*DIM[0];    ///@todo I don't like the default dimension to be 2


// === Constructors / Destructor - BEGIN =================

//------------------------------------------------------------------------------------------------------
  Mesh::Mesh() {

    _coarseMsh = NULL;

    for(int i = 0; i < NFE_FAMS; i++) {
      _ProjCoarseToFine[i] = NULL;
    }

  }


  Mesh::~Mesh() {
      
    delete el;
    
    _topology->FreeSolutionVectors();
    delete _topology;

    for(unsigned i = 0; i < NFE_FAMS; i++) {
      if(_ProjCoarseToFine[i]) {
        delete _ProjCoarseToFine[i];
        _ProjCoarseToFine[i] = NULL;
      }
    }
    
  }


// === Constructors / Destructor - END =================
  
  

/// print Mesh info
  void Mesh::PrintInfo() const {

      PrintInfoLevel();
      PrintInfoElements();
      PrintInfoNodes();
      std::cout << std::endl;

  }
  
  
  void Mesh::PrintInfoLevel() const {
    std::cout << " Mesh Level                  : " <<  GetLevel()  << std::endl;
  }
  
  
  void Mesh::PrintInfoElements() const {
    std::cout << " Number of elements          : " << _nelem  << std::endl;
  }
  
  void Mesh::PrintInfoNodes() const {
    std::cout << " Number of linear nodes      : " << _dofOffset[0][_nprocs] << std::endl;
    std::cout << " Number of quadratic nodes   : " << _dofOffset[1][_nprocs] << std::endl;
    std::cout << " Number of biquadratic nodes : " << _dofOffset[2][_nprocs] << std::endl;
  }
  
  
  
  const unsigned Mesh::_numberOfMissedBiquadraticNodes[N_GEOM_ELS] = {0, 5, 3, 0, 1, 0};
  const double Mesh::_baricentricWeight[N_GEOM_ELS][NFE_FAMS][MAXIMUM_NUMBER_OF_NON_BIQUADRATIC_NODES_TO_USE_IN_WEIGHTING] = {
    {},
    {
      { -1. / 9., -1. / 9., -1. / 9.,  0    , 4. / 9., 4. / 9., 4. / 9., 0.   , 0.   , 0.   },
      { -1. / 9., -1. / 9.,  0.   , -1. / 9., 4. / 9., 0.   , 0.   , 4. / 9., 4. / 9., 0.   },
      {  0.   , -1. / 9., -1. / 9., -1. / 9., 0.   , 4. / 9., 0.   , 0.   , 4. / 9., 4. / 9.},
      { -1. / 9.,  0.    , -1. / 9., -1. / 9., 0.   , 0.   , 4. / 9., 4. / 9., 0.   , 4. / 9.},
      { -1. / 8., -1. / 8., -1. / 8., -1. / 8., 1. / 4., 1. / 4., 1. / 4., 1. / 4., 1. / 4., 1. / 4.}
    },
    {
      { -1. / 9., -1. / 9., -1. / 9., 0.   ,  0.   ,  0.   , 4. / 9., 4. / 9., 4. / 9.},
      {  0.   ,  0.   ,  0.   , -1. / 9., -1. / 9., -1. / 9., 0.   , 0.   , 0.   , 4. / 9., 4. / 9., 4. / 9.},
      {
        0.   ,  0.   ,  0.   , 0.   ,  0.   ,  0.   , 0.   , 0.   , 0.   , 0.   , 0.   , 0.   ,
        -1. / 9., -1. / 9., -1. / 9., 4. / 9.,  4. / 9.,  4. / 9.
      },
    },
    {},
    {{ -1. / 9., -1. / 9., -1. / 9., 4. / 9., 4. / 9., 4. / 9.}}
  };


  
  std::vector < unsigned > Mesh::PartitionForElements_refinement(const bool AMR, const Mesh* mshc) const {
  
    std::vector < unsigned > partition;

    partition.resize(/*_mesh.*/GetNumberOfElements());

    MeshMetisPartitioning meshMetisPartitioning(*this/*_mesh*/);

    if(AMR == true) {
      meshMetisPartitioning.DoPartition(partition, AMR);
    }
    else {
      meshMetisPartitioning.DoPartition(partition, *mshc);
    }
  
    return partition;

  }
  
  

  std::vector < unsigned >  Mesh::PartitionForElements() const {

    std::vector < unsigned > partition;
    
    partition.resize(GetNumberOfElements());
    
    const bool flag_for_ncommon_in_metis = false;
    
    MeshMetisPartitioning meshMetisPartitioning(*this);
    
    meshMetisPartitioning.DoPartition(partition, flag_for_ncommon_in_metis);

    return partition;
    
  }
  
  
  void Mesh::PartitionElements_and_FillDofMapAllFEFamilies() {

    std::vector < unsigned > partition = PartitionForElements();
    
    FillISvectorDofMapAllFEFamilies(partition);
    
    std::vector<unsigned> ().swap(partition);

  }


  /**
   *  This function generates the coarse Mesh level, $l_0$, from an input Mesh file
   **/
  void Mesh::ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool>& type_elem_flag) {


    const bool read_groups = true; //by default groups are read
    const bool read_boundary_groups = true; //by default boundary groups are read

    ReadCoarseMesh(name, Lref, type_elem_flag, read_groups, read_boundary_groups);

  }


  void Mesh::ReadCoarseMeshFile(const std::string& name, const double Lref, std::vector<bool>& type_elem_flag, const bool read_groups, const bool read_boundary_groups) {


    if(name.rfind(".neu") < name.size()) {
      GambitIO(*this).read(name, _coords, Lref, type_elem_flag, read_groups, read_boundary_groups);
    }

    // else if (name.rfind (".obj") < name.size()) {
    //   obj_io (*this).read (name, _coords, Lref, type_elem_flag);
    // }

#ifdef HAVE_HDF5
    else if(name.rfind(".med") < name.size()) {
      MED_IO(*this).read(name, _coords, Lref, type_elem_flag, read_groups, read_boundary_groups);
    }
#endif
    else {
      std::cerr << " ERROR: Unrecognized file extension: " << name
                << "\n   I understand the following:\n\n"
                << "     *.neu -- Gambit Neutral File\n"
                << "     *.med -- MED File\n"
                << std::endl;
      abort();
    }

  }
  
  

  void Mesh::ReadCoarseMeshBeforePartitioning(const std::string& name, const double Lref, std::vector<bool> & type_elem_flag, const bool read_groups, const bool read_boundary_groups) {

//==== Level - BEGIN ==============================
    SetLevel(0);
//==== Level - END ==============================

//==== AMR - BEGIN ==============================
    SetIfHomogeneous(true);
//==== AMR - END ==============================

//==== Coords, coarse, then to topology - BEGIN ==============================
    _coords.resize(3);
//==== Coords, coarse, then to topology - END ==============================

    

    ReadCoarseMeshFile(name, Lref, type_elem_flag, read_groups, read_boundary_groups);



//==== Only in Coarse generations - BEGIN ==============================
    AddBiquadraticNodesNotInMeshFile();

    el->ShrinkToFitElementDof();
    el->ShrinkToFitElementNearFace();
//==== Only in Coarse generations - END ==============================


  }  
  

  /**
   *  This function generates the coarse Mesh level, $l_0$, from an input Mesh file
   **/
  void Mesh::ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool>& type_elem_flag, const bool read_groups, const bool read_boundary_groups) {
      

    ReadCoarseMeshBeforePartitioning(name, Lref, type_elem_flag, read_groups, read_boundary_groups);
    
    PartitionElements_and_FillDofMapAllFEFamilies();

    BuildElementAndNodeStructures();
   
//====  CharacteristicLength ======== 
    SetCharacteristicLengthOfCoarsestLevel();  //doesn't need dofmap

//====  Print Info ======== 
    PrintInfo();  //needs dofmap
    
    
  }
  
  

  void Mesh::BuildElementAndNodeStructures() {
      
    GetMeshElements()->BuildElem_NearFace_NearElem_using_NearVertex();  //must stay here, cannot be anticipated. Does it need dofmap already? I don't think so, but it needs the elem reordering and maybe also the node reordering
    
    GetMeshElements()->ScatterElement_Level_Type_Group_Material___NearFace();
    
    GetMeshElements()->ScatterElementDof();
    
    
    //------

    
    //needs dofmap
//====  Topology, Coordinates - BEGIN ======== 
    Topology_InitializeCoordinates();
    Topology_FillCoordinates();
//====  Topology, Coordinates - END ======== 


//====  Topology, AMR - BEGIN ======== 
    Topology_InitializeAMR();  // needs to be here because the order of AddSolution is that 
//====  Topology, AMR - END ======== 
    

//====  Topology, Solidnodeflag - BEGIN ======== 
    Topology_InitializeSolidNodeFlag();
    Topology_FillSolidNodeFlag();
//====  Topology, Solidnodeflag - END ======== 


    InitializeAndPossiblyFillAmrRestriction(false);       /** @todo is it really needed here? */
    
  }

  
  void Mesh::InitializeAndPossiblyFillAmrRestriction(const bool amr) {

    if(amr) {
      
      GetMeshElements()->GetAMRRestriction(this, this->GetMeshElements());
//       for(unsigned soltype = 0; soltype < NFE_FAMS_C_ZERO_LAGRANGE; soltype++) {
//         std::cout << "solution type = " << soltype << std::endl;
//         for(std::map<unsigned, std::map<unsigned, double> >::iterator it1 = _mesh.GetAmrRestrictionMap()[soltype].begin(); it1 != _mesh.GetAmrRestrictionMap()[soltype].end(); it1++) {
//           std::cout << it1->first << "\t";
//           for(std::map<unsigned, double> ::iterator it2 = _mesh.GetAmrRestrictionMap()[soltype][it1->first].begin(); it2 != _mesh.GetAmrRestrictionMap()[soltype][it1->first].end(); it2++) {
//             std::cout << it2->first << " (" << it2->second << ")  ";
// 
//           }
//           std::cout << std::endl;
//         }
//       }
      
    }
    else{

     _amrRestriction.resize(NFE_FAMS_C_ZERO_LAGRANGE);
     
    }
    

  }
  



 /// this needs all the dof maps, for the continuous Lagrange elements
  void Mesh::Topology_InitializeCoordinates() {
      
    _topology = new Solution(this);

    _topology->AddSolution("X", LAGRANGE, SECOND, 1, 0);
    _topology->AddSolution("Y", LAGRANGE, SECOND, 1, 0);
    _topology->AddSolution("Z", LAGRANGE, SECOND, 1, 0);

    _topology->ResizeSolutionVector("X");  //needs dofmap
    _topology->ResizeSolutionVector("Y");  //needs dofmap
    _topology->ResizeSolutionVector("Z");  //needs dofmap

  }
  

  void Mesh::Topology_FillCoordinates() {
    
    //set coordinates -----------
    _topology->GetSolutionByName("X") = _coords[0];
    _topology->GetSolutionByName("Y") = _coords[1];
    _topology->GetSolutionByName("Z") = _coords[2];
    //set coordinates -----------

  }
  
  
  
  
/// This needs the dof maps, for the discontinuous Lagrange elements
  void Mesh::Topology_InitializeAMR() {
      
    _topology->AddSolution("AMR", DISCONTINUOUS_POLYNOMIAL, ZERO, 1, 0);

    _topology->ResizeSolutionVector("AMR");

  }
  

/// This needs the dof maps, for the continuous Lagrange elements
  void Mesh::Topology_InitializeSolidNodeFlag() {
      
    _topology->AddSolution("solidMrk", LAGRANGE, SECOND, 1, 0);
    
    _topology->ResizeSolutionVector("solidMrk");

  }
  

  
  void Mesh::SetCharacteristicLengthOfCoarsestLevel() {
      
    //compute max and min coords -----------
    std::vector < double > xMax(3, 0.);
    std::vector < double > xMin(3, 0.);
    for(unsigned i = 0; i < _coords[0].size(); i++) {
      for(unsigned k = 0; k < 3; k++) {
        if(xMax[k] < _coords[k][i]) xMax[k] = _coords[k][i];
        if(xMin[k] > _coords[k][i]) xMin[k] = _coords[k][i];
      }
    }
    //compute max and min coords - end -----------
    
    _cLength = sqrt(pow(xMax[0] - xMin[0], 2) + pow(xMax[1] - xMin[1], 2) + pow(xMax[2] - xMin[2], 2));
    
  }
  
  
  /**
   *  This function generates the coarse Box Mesh level using the built-in generator
   *   ///@todo seems like GenerateCoarseBoxMesh doesn't assign flags to faces correctly, need to check that
   **/
  void Mesh::GenerateCoarseBoxMesh(
    const unsigned int nx, const unsigned int ny, const unsigned int nz,
    const double xmin, const double xmax,
    const double ymin, const double ymax,
    const double zmin, const double zmax,
    const ElemType elemType, std::vector<bool>& type_elem_flag) {
    

//==== Level - BEGIN ==============================
    SetLevel(0);
//==== Level - END ==============================

//==== AMR - BEGIN ==============================
    SetIfHomogeneous(true);
//==== AMR - END ==============================

//==== Coords, coarse, then to topology - BEGIN ==============================
    _coords.resize(3);
//==== Coords, coarse, then to topology - END ==============================
    

    
    MeshTools::Generation::BuildBox(*this, _coords, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, elemType, type_elem_flag);


    
//==== Only in Coarse generations - BEGIN ==============================
    AddBiquadraticNodesNotInMeshFile();

    el->ShrinkToFitElementDof();
    el->ShrinkToFitElementNearFace();
//==== Only in Coarse generations - END ==============================


    PartitionElements_and_FillDofMapAllFEFamilies();

    BuildElementAndNodeStructures();
    
//====  CharacteristicLength ======== 
    SetCharacteristicLengthOfCoarsestLevel();  //doesn't need dofmap

//====  Print Info ======== 
    PrintInfo();  //needs dofmap
    
  }


  void Mesh::Topology_FillSolidNodeFlag() {


    NumericVector& NodeMaterial =  _topology->GetSolutionByName("solidMrk");

    NodeMaterial.zero();

    for(int iel = _elementOffset[_iproc]; iel < _elementOffset[_iproc + 1]; iel++) {
      int flag_mat = GetElementMaterial(iel);

      if(flag_mat == 4) { ///@todo Where on Earth do we say that 4 is a special flag for the solid
        unsigned elementType = GetElementType(iel);
        unsigned nve = el->GetNVE(elementType, CONTINUOUS_BIQUADRATIC);

        for(unsigned i = 0; i < nve; i++) {
          unsigned inode = GetSolutionDof(i, iel, CONTINUOUS_BIQUADRATIC);
          NodeMaterial.set(inode, 1);
        }
      }
    }

    NodeMaterial.close();

  }


  
  void Mesh::SetFiniteElementPtr(/*const*/ elem_type* OtherFiniteElement[N_GEOM_ELS][NFE_FAMS]) {
      
    for(int i = 0; i < N_GEOM_ELS; i++)
      for(int j = 0; j < NFE_FAMS; j++)
        _finiteElement[i][j] = OtherFiniteElement[i][j];
      
  }

  
 
  

  std::vector <unsigned> Mesh::dofmap_Node_based_dof_offsets_Compute_Node_mapping_and_Node_ownSize() {
  // at this point the elements have been reordered, but not the nodes. The new node numbering starting from the med node numbering is happening here
      

     std::vector<unsigned>  partition(GetNumberOfNodes(), _nprocs);
    std::vector <unsigned>  mapping_from_mesh_file_to_femus(GetNumberOfNodes(), 0);

    for(unsigned k = 0; k < NFE_FAMS_C_ZERO_LAGRANGE; k++) {
      _ownSize[k].assign(_nprocs, 0);
    }

    unsigned counter = 0;

    for(int isdom = 0; isdom < _nprocs; isdom++) {
        
      for(unsigned k = 0; k < NFE_FAMS_C_ZERO_LAGRANGE; k++) {
        for(unsigned iel = _elementOffset[isdom]; iel < _elementOffset[isdom + 1]; iel++) {
          unsigned nodeStart = (k == 0) ? 0 : el->GetElementDofNumber(iel, k - 1);
          unsigned nodeEnd = el->GetElementDofNumber(iel, k);

          for(unsigned inode = nodeStart; inode < nodeEnd; inode++) {
            unsigned ii = el->GetElementDofIndex(iel, inode); //at this point the elements were reordered, but not the nodes

            if(partition[ii] > isdom) {
              partition[ii] = isdom;
              mapping_from_mesh_file_to_femus[ii] = counter;
              counter++;

              for(int j = k; j < NFE_FAMS_C_ZERO_LAGRANGE; j++) {
                _ownSize[j][isdom]++;
              }
            }
          }
        }
      }
      
    }
  
    
    
    return mapping_from_mesh_file_to_femus;
  
  }
// *******************************************************


  void Mesh::initialize_elem_offsets() {
      
    _elementOffset.resize(_nprocs + 1);
    _elementOffset[0] = 0;
    
  }
  

  void Mesh::dofmap_all_fe_families_initialize() {
      
    //BEGIN Initialization for k = 0,1,2,3,4
    for(int k = 0; k < NFE_FAMS; k++) {
      _dofOffset[k].resize(_nprocs + 1);
      _dofOffset[k][0] = 0;
    }
    
    for(int k = 0; k < NFE_FAMS; k++) {
              _ownSize[k].assign(_nprocs, 0); 
            _ghostDofs[k].resize(_nprocs);
    }
    //END Initialization for k = 0,1,2,3,4

    
  }

  
   void Mesh::build_elem_offsets(const std::vector <unsigned> & partition)  {
       
          //BEGIN building the  metis2mesh_file element list 
    std::vector <unsigned>  element_mapping(GetNumberOfElements());

    unsigned counter = 0;

    for(int isdom = 0; isdom < _nprocs; isdom++) {  // isdom = iprocess
      for(unsigned iel = 0; iel < GetNumberOfElements(); iel++) {
        if(partition[iel] == isdom) {
          //filling the Metis to Mesh element mapping
          element_mapping[ iel ] = counter;
          counter++;
          _elementOffset[isdom + 1] = counter;
        }
      }
    }

    el->ReorderMeshElement_Type_Level_Group_Material___NearFace_rows_ChildElem_columns(element_mapping);  ///this is needed because later there will be another reordering based on Group and Material

    
    
    //Dof
    el->ReorderMeshElement_Dof_stuff(element_mapping);
    
    
    
   }
   
   
   
  
   void Mesh::mesh_reorder_elem_quantities()  {
 

//======== Inverse and direct element mapping - BEGIN =============
     
    std::vector < unsigned > inverse_element_mapping(GetNumberOfElements());

    for(unsigned iel = 0; iel < GetNumberOfElements(); iel++) {
      inverse_element_mapping[iel] = iel;
    }
    

    for(int isdom = 0; isdom < _nprocs; isdom++) {

// Old, much slower (while below is better) **********
//       for (unsigned i = _elementOffset[isdom]; i < _elementOffset[isdom + 1] - 1; i++) {
//         unsigned iel = inverse_element_mapping[i];
//         unsigned ielMat = el->GetElementMaterial (iel);
//         unsigned ielGroup = el->GetElementGroup (iel);
//         for (unsigned j = i + 1; j < _elementOffset[isdom + 1]; j++) {
//           unsigned jel = inverse_element_mapping[j];
//           unsigned jelMat = el->GetElementMaterial (jel);
//           unsigned jelGroup = el->GetElementGroup (jel);
//           if (jelMat < ielMat || (jelMat == ielMat && jelGroup < ielGroup || (jelGroup == ielGroup && iel > jel))) {
//             inverse_element_mapping[i] = jel;
//             inverse_element_mapping[j] = iel;
//             iel = jel;
//             ielMat = jelMat;
//             ielGroup = jelGroup;
//           }
//         }
//       }

      unsigned jel, iel;
      short unsigned jelMat, jelGroup, ielMat, ielGroup;

      unsigned n = _elementOffset[isdom + 1u] - _elementOffset[isdom];
      
      while(n > 1) {
        unsigned newN = 0u;
        for(unsigned j = _elementOffset[isdom] + 1u; j < _elementOffset[isdom] + n ; j++) {
          jel = inverse_element_mapping[j];
          jelMat = el->GetElementMaterial(jel);
          jelGroup = el->GetElementGroup(jel);

          iel = inverse_element_mapping[j - 1];
          ielMat = el->GetElementMaterial(iel);
          ielGroup = el->GetElementGroup(iel);

          if( jelMat < ielMat || (jelMat == ielMat && (jelGroup < ielGroup || (jelGroup == ielGroup && jel < iel) ) ) ) {
            inverse_element_mapping[j - 1] = jel;
            inverse_element_mapping[j] = iel;
            newN = j - _elementOffset[isdom];
          }
        }
        n = newN;
      }

    }
    

    std::vector <unsigned>  element_mapping(GetNumberOfElements());
        
    for(unsigned i = 0; i < GetNumberOfElements(); i++) {
      element_mapping[inverse_element_mapping[i]] = i;
    }
    
//======== Inverse and direct element mapping - END =============


    el->ReorderMeshElement_Type_Level_Group_Material___NearFace_rows_ChildElem_columns(element_mapping);
    
    
    
    //Dof
    el->ReorderMeshElement_Dof_stuff(element_mapping);
    
    
    //END building the  metis2mesh_file element list 

    
}



   void Mesh::dofmap_Element_based_dof_offsets_build()  {
       
    //BEGIN building element based dofs -  k = 3,4 

    // ghost vs owned nodes: 3 and 4 have no ghost nodes
    for(unsigned k = NFE_FAMS_C_ZERO_LAGRANGE; k < NFE_FAMS; k++) {
      _ownSize[k].assign(_nprocs, 0);
    }

    for(int isdom = 0; isdom < _nprocs; isdom++) {
      _ownSize[3][isdom] = _elementOffset[isdom + 1] - _elementOffset[isdom];
      _ownSize[4][isdom] = (_elementOffset[isdom + 1] - _elementOffset[isdom]) * (_dimension + 1);
    }

    for(int k = NFE_FAMS_C_ZERO_LAGRANGE; k < NFE_FAMS; k++) {
      _ghostDofs[k].resize(_nprocs);

      for(int isdom = 0; isdom < _nprocs; isdom++) {
        _dofOffset[k][isdom + 1] = _dofOffset[k][isdom] + _ownSize[k][isdom];
        _ghostDofs[k][isdom].resize(0);
      }
    }

    //END building element based dofs -  k = 3,4

   }
   
   

   void Mesh::dofmap_Node_based_dof_offsets_Continue_biquadratic()  {

    for(int i = 1 ; i <= _nprocs; i++) {
      _dofOffset[2][i] = _dofOffset[2][i - 1] + _ownSize[2][i - 1];
    }

   }


   void Mesh::mesh_reorder_node_quantities(const std::vector <unsigned> & mapping)  {
       
     
    el->ReorderElementDof_columns_Using_node_mapping(mapping);

    //reorder coordinate vector at coarse level - BEGIN ----
    if(GetLevel() == 0) {
      std::vector <double> coords_temp;

      for(int i = 0; i < 3; i++) {
        coords_temp = _coords[i];

        for(unsigned j = 0; j < GetNumberOfNodes(); j++) {
          _coords[i][mapping[j]] = coords_temp[j];
        }
      }
    }
    //reorder coordinate vector at coarse level - END ----


   }
   

    void Mesh::dofmap_Node_based_dof_offsets_Ghost_nodes_search_Complete_biquadratic() {
 
    for(int k = 0; k < NFE_FAMS_C_ZERO_LAGRANGE; k++) {
      _ghostDofs[k].resize(_nprocs);

      for(int isdom = 0; isdom < _nprocs; isdom++) {
        std::map < unsigned, bool > ghostMap;

        for(unsigned iel = _elementOffset[isdom]; iel < _elementOffset[isdom + 1]; iel++) {
          for(unsigned inode = 0; inode < el->GetElementDofNumber(iel, k); inode++) {
            unsigned ii = el->GetElementDofIndex(iel, inode);

            if(ii < _dofOffset[2][isdom]) {
              ghostMap[ii] = true;
            }
          }
        }

        _ghostDofs[k][isdom].resize(ghostMap.size());
        unsigned counter = 0;

        for(std::map < unsigned, bool >::iterator it = ghostMap.begin(); it != ghostMap.end(); it++) {
          _ghostDofs[k][isdom][counter] = it->first;
          counter++;
        }
      }
    }

  }
  

  
    void Mesh::dofmap_Node_based_dof_offsets_Complete_linear_quadratic() {

    //BEGIN completing k = 0, 1

    for(unsigned k = 0; k < 2; k++) {

      std::vector < unsigned > ownedGhostCounter(_nprocs , 0);
      unsigned counter = 0;

      _originalOwnSize[k].resize(_nprocs);

      for(int isdom = 0; isdom < _nprocs; isdom++) {

        //owned nodes
        for(unsigned inode = _dofOffset[2][isdom]; inode < _ownSize[k][isdom] + _dofOffset[2][isdom]; inode++) {
          counter++;
        }

        for(unsigned inode = 0; inode < _ghostDofs[k][isdom].size(); inode++) {
          unsigned ghostNode = _ghostDofs[k][isdom][inode];

          unsigned ksdom = BisectionSearch_find_processor_of_dof(ghostNode, 2);

          int upperBound = _dofOffset[2][ksdom] + _ownSize[k][ksdom];

          if(ghostNode < upperBound) {
            _ghostDofs[k][isdom][inode] =  ghostNode  - _dofOffset[2][ksdom] + _dofOffset[k][ksdom];
          }
          else if(_ownedGhostMap[k].find(ghostNode) != _ownedGhostMap[k].end()) {
            _ghostDofs[k][isdom][inode] =  _ownedGhostMap[k][ghostNode];
          }
          else { // owned ghost nodes
            _ownedGhostMap[k][ ghostNode ] = counter;
            counter++;
            ownedGhostCounter[isdom]++;

            for(unsigned jnode = inode; jnode < _ghostDofs[k][isdom].size() - 1; jnode++) {
              _ghostDofs[k][isdom][jnode] = _ghostDofs[k][isdom][jnode + 1];
            }

            _ghostDofs[k][isdom].resize(_ghostDofs[k][isdom].size() - 1);
            inode--;
          }
        }

        _originalOwnSize[k][isdom] = _ownSize[k][isdom];
        _ownSize[k][isdom] += ownedGhostCounter[isdom];
        _dofOffset[k][isdom + 1] = _dofOffset[k][isdom] + _ownSize[k][isdom];
      }
    }

    //END completing for k = 0, 1
    

  }
  
  
  void Mesh::dofmap_all_fe_families_clear_ghost_dof_list_for_other_procs() {
      
    //delete ghost dof list all but _iproc
    for(int isdom = 0; isdom < _nprocs; isdom++) {
      if(isdom != _iproc)
        for(int k = 0; k < NFE_FAMS; k++) {
          _ghostDofs[k][isdom].resize(0);
        }
    }
      

  }
  
 
  void Mesh::dofmap_build_all_fe_families_and_elem_and_node_structures() {
      
            dofmap_all_fe_families_initialize();

           std::vector < unsigned > elem_partition_from_mesh_file_to_new = elem_offsets();   ///@todo these are actually unused returns
  
           std::vector < unsigned > node_mapping_from_mesh_file_to_new = node_offsets();     ///@todo these are actually unused returns
           
           dofmap_all_fe_families_clear_ghost_dof_list_for_other_procs();
       
           BuildElementAndNodeStructures();

//====  CharacteristicLength ======== 
           SetCharacteristicLengthOfCoarsestLevel();  //doesn't need dofmap

//====  Print Info ======== 
           PrintInfo();  //needs dofmap
    
  }

  
  
  void Mesh::set_node_counts() {
      
    SetNumberOfNodes(_dofOffset[ CONTINUOUS_BIQUADRATIC ][_nprocs]);
    el->SetNodeNumber(_dofOffset[CONTINUOUS_BIQUADRATIC ][_nprocs]);

  }
  

  void Mesh::set_elem_counts_per_subdomain() {
      
    el->SetElementOffsets(_elementOffset, _iproc, _nprocs);

  }
  
  
  void Mesh::deallocate_node_mapping(std::vector < unsigned > & node_mapping) const {
      
    std::vector<unsigned> ().swap(node_mapping);    //     node_mapping.resize(0);  //resize DOES NOT FREE memory!!!
      
  }

  
  
  // // // ======== ELEM OFFSETS =========================================================  
  void Mesh::FillISvectorElemOffsets(std::vector < unsigned >& partition) {
      
   
   initialize_elem_offsets();
   
   build_elem_offsets(partition);
   
   mesh_reorder_elem_quantities();
   
   
   set_elem_counts_per_subdomain();
   
      
  }

  std::vector < unsigned >  Mesh::elem_offsets() {

      std::vector < unsigned > elem_partition_from_mesh_file_to_new = PartitionForElements(); 
           
      FillISvectorElemOffsets(elem_partition_from_mesh_file_to_new);
           
      dofmap_Element_based_dof_offsets_build();
      
      return elem_partition_from_mesh_file_to_new;
      
}
  
  
  // // // ======== NODE OFFSETS =========================================================  
  void Mesh::FillISvectorNodeOffsets() {
      
    std::vector < unsigned > node_mapping =  dofmap_Node_based_dof_offsets_Compute_Node_mapping_and_Node_ownSize();
    mesh_reorder_node_quantities(node_mapping);
    
    deallocate_node_mapping(node_mapping);
      
  }
  
  std::vector < unsigned >  Mesh::node_offsets() {
      
           std::vector < unsigned > node_mapping_from_mesh_file_to_new = dofmap_Node_based_dof_offsets_Compute_Node_mapping_and_Node_ownSize();
           
           mesh_reorder_node_quantities(node_mapping_from_mesh_file_to_new);
           
           dofmap_Node_based_dof_offsets_Continue_biquadratic();
           dofmap_Node_based_dof_offsets_Ghost_nodes_search_Complete_biquadratic();
           dofmap_Node_based_dof_offsets_Complete_linear_quadratic();
           
           set_node_counts();
  
      return node_mapping_from_mesh_file_to_new;
}

 /**
  * dof map: piecewise liner 0, quadratic 1, bi-quadratic 2, piecewise constant 3, piecewise linear discontinuous 4
  */
  void Mesh::FillISvectorDofMapAllFEFamilies(std::vector < unsigned >& partition) {

        dofmap_all_fe_families_initialize();
     
      //BEGIN  k = 3, 4
     FillISvectorElemOffsets(partition);
     
        dofmap_Element_based_dof_offsets_build();  
      //END  k = 3, 4
        

      //BEGIN building for k = 0,1,2
     FillISvectorNodeOffsets();
    
      dofmap_Node_based_dof_offsets_Continue_biquadratic();
      
        //BEGIN ghost nodes search k = 0, 1, 2
      dofmap_Node_based_dof_offsets_Ghost_nodes_search_Complete_biquadratic();
        //END ghost nodes search k = 0, 1, 2
      
      // --END building for k = 2, but incomplete for k = 0, 1
      
      //BEGIN completing for k = 0,1
      dofmap_Node_based_dof_offsets_Complete_linear_quadratic();
      //END completing for k = 0,1
    
    set_node_counts(); //also, it shouldn't use dofOffset
    
      //END building for k = 0,1,2

     
     
        dofmap_all_fe_families_clear_ghost_dof_list_for_other_procs();
    
  }


// *******************************************************
  unsigned Mesh::BisectionSearch_find_processor_of_dof(const unsigned& dof, const short unsigned& solType) const {

    unsigned isdom0 = 0;
    unsigned isdom1 = _nprocs ;
    unsigned isdom = _iproc;

    while(dof < _dofOffset[solType][isdom] || dof >= _dofOffset[solType][isdom + 1]) {
      if(dof < _dofOffset[solType][isdom]) isdom1 = isdom;
      else isdom0 = isdom + 1;

      isdom = (isdom0 + isdom1) / 2;
    }

    return isdom;
  }
// *******************************************************

  unsigned Mesh::GetSolutionDof(const unsigned& i, const unsigned& iel, const short unsigned& solType) const {

    unsigned dof;

    switch(solType) {
        
      case CONTINUOUS_LINEAR: { // linear Lagrange
        unsigned iNode = el->GetElementDofIndex(iel, i);  //GetMeshDof(iel, i, solType);
        unsigned isdom = BisectionSearch_find_processor_of_dof(iNode, CONTINUOUS_BIQUADRATIC);

        if(iNode < _dofOffset[CONTINUOUS_BIQUADRATIC][isdom] + _originalOwnSize[CONTINUOUS_LINEAR][isdom]) {
          dof = (iNode - _dofOffset[CONTINUOUS_BIQUADRATIC][isdom]) + _dofOffset[CONTINUOUS_LINEAR][isdom];
        }
        else {
          dof = _ownedGhostMap[CONTINUOUS_LINEAR].find(iNode)->second;
        }
      }
      break;

      case CONTINUOUS_SERENDIPITY: { // quadratic Lagrange
        unsigned iNode = el->GetElementDofIndex(iel, i);  //GetMeshDof(iel, i, solType);
        unsigned isdom = BisectionSearch_find_processor_of_dof(iNode, CONTINUOUS_BIQUADRATIC);

        if(iNode < _dofOffset[CONTINUOUS_BIQUADRATIC][isdom] + _originalOwnSize[CONTINUOUS_SERENDIPITY][isdom]) {
          dof = (iNode - _dofOffset[CONTINUOUS_BIQUADRATIC][isdom]) + _dofOffset[CONTINUOUS_SERENDIPITY][isdom];
        }
        else {
          dof = _ownedGhostMap[CONTINUOUS_SERENDIPITY].find(iNode)->second;
        }
      }
      break;

      case CONTINUOUS_BIQUADRATIC: // bi-quadratic Lagrange
        dof = el->GetElementDofIndex(iel, i);  //GetMeshDof(iel, i, solType);
        break;

      case DISCONTINUOUS_CONSTANT: // piecewise constant
        // in this case use i=0
        dof = iel;
        break;

      case DISCONTINUOUS_LINEAR: // piecewise linear discontinuous
        unsigned isdom = BisectionSearch_find_processor_of_dof(iel, DISCONTINUOUS_CONSTANT);
        unsigned offset = _elementOffset[isdom];
        unsigned offsetp1 = _elementOffset[isdom + 1];
        unsigned ownSize = offsetp1 - offset;
        unsigned offsetPWLD = offset * (_dimension + 1);
        unsigned locIel = iel - offset;
        dof = offsetPWLD + (i * ownSize) + locIel;
        break;
    }

    return dof;
  }

// *******************************************************

  unsigned Mesh::GetSolutionDof(const unsigned& ielc, const unsigned& i0, const unsigned& i1, const short unsigned& solType, const Mesh* mshc) const {

    unsigned dof;

    switch(solType) {
        
      case CONTINUOUS_LINEAR: { // linear Lagrange
        unsigned iNode = mshc->el->GetChildElementDof(ielc, i0, i1);
        unsigned isdom = BisectionSearch_find_processor_of_dof(iNode, CONTINUOUS_BIQUADRATIC);

        if(iNode < _dofOffset[CONTINUOUS_BIQUADRATIC][isdom] + _originalOwnSize[CONTINUOUS_LINEAR][isdom]) {
          dof = (iNode - _dofOffset[CONTINUOUS_BIQUADRATIC][isdom]) + _dofOffset[CONTINUOUS_LINEAR][isdom];
        }
        else {
          dof = _ownedGhostMap[CONTINUOUS_LINEAR].find(iNode)->second;
        }
      }
      break;

      case CONTINUOUS_SERENDIPITY: { // quadratic Lagrange
        unsigned iNode = mshc->el->GetChildElementDof(ielc, i0, i1);
        unsigned isdom = BisectionSearch_find_processor_of_dof(iNode, CONTINUOUS_BIQUADRATIC);

        if(iNode < _dofOffset[CONTINUOUS_BIQUADRATIC][isdom] + _originalOwnSize[CONTINUOUS_SERENDIPITY][isdom]) {
          dof = (iNode - _dofOffset[CONTINUOUS_BIQUADRATIC][isdom]) + _dofOffset[CONTINUOUS_SERENDIPITY][isdom];
        }
        else {
          dof = _ownedGhostMap[CONTINUOUS_SERENDIPITY].find(iNode)->second;
        }
      }
      break;

      case CONTINUOUS_BIQUADRATIC: // bi-quadratic Lagrange
        dof = mshc->el->GetChildElementDof(ielc, i0, i1);
        break;

      case DISCONTINUOUS_CONSTANT: // piecewise constant
        // in this case use i=0
        dof = mshc->el->GetChildElement(ielc, i0);
        break;

      case DISCONTINUOUS_LINEAR: // piecewise linear discontinuous
        unsigned iel = mshc->el->GetChildElement(ielc, i0);
        unsigned isdom = BisectionSearch_find_processor_of_dof(iel, DISCONTINUOUS_CONSTANT);
        unsigned offset = _elementOffset[isdom];
        unsigned offsetp1 = _elementOffset[isdom + 1];
        unsigned ownSize = offsetp1 - offset;
        unsigned offsetPWLD = offset * (_dimension + 1);
        unsigned locIel = iel - offset;
        dof = offsetPWLD + (i1 * ownSize) + locIel;
        break;
    }

    return dof;
  }

// *******************************************************

  
  
  
  
  SparseMatrix* Mesh::GetCoarseToFineProjectionRestrictionOnCoarse(const unsigned& solType) {

    if(solType >= NFE_FAMS) {
      std::cout << "Wrong argument range in function \"GetCoarseToFineProjection\": "
                << "solType is greater then SolTypeMax" << std::endl;
      abort();
    }

    if(_ProjCoarseToFine[solType])
      BuildCoarseToFineProjection(solType, "coarse");

    return _ProjCoarseToFine[solType];
  }


  SparseMatrix* Mesh::GetCoarseToFineProjection(const unsigned& solType) {

    if(solType >= NFE_FAMS) {
      std::cout << "Wrong argument range in function \"GetCoarseToFineProjection\": "
                << "solType is greater then SolTypeMax" << std::endl;
      abort();
    }

    if(!_ProjCoarseToFine[solType])
      BuildCoarseToFineProjection(solType, "fine");

    return _ProjCoarseToFine[solType];
  }



  void Mesh::BuildCoarseToFineProjection(const unsigned& solType, const char el_dofs[]) {

    if(!_coarseMsh) {
      std::cout << "Error! In function \"BuildCoarseToFineProjection\": the coarse mesh has not been set" << std::endl;
      abort();
    }

    if(!_ProjCoarseToFine[solType]) {

    // ------------------- Sparsity pattern size - BEGIN
      const int nf     = _dofOffset[solType][_nprocs];
      const int nc     = _coarseMsh->_dofOffset[solType][_nprocs];
      const int nf_loc = _ownSize[solType][_iproc];
      const int nc_loc = _coarseMsh->_ownSize[solType][_iproc];

      //build matrix sparsity pattern size
      NumericVector* NNZ_d = NumericVector::build().release();

      if(n_processors() == 1) {  // IF SERIAL
        NNZ_d->init(nf, nf_loc, false, SERIAL);
      }
      else { // IF PARALLEL
        if(solType < NFE_FAMS_C_ZERO_LAGRANGE) {  // GHOST nodes only for Lagrange FE families
          NNZ_d->init(nf, nf_loc, _ghostDofs[solType][processor_id()], false, GHOSTED);
        }
        else { //piecewise discontinuous variables have no ghost nodes
          NNZ_d->init(nf, nf_loc, false, PARALLEL);
        }
      }

      NNZ_d->zero();

      NumericVector* NNZ_o = NumericVector::build().release();
      NNZ_o->init(*NNZ_d);
      NNZ_o->zero();

      for(int isdom = _iproc; isdom < _iproc + 1; isdom++) {
        for(int iel = _coarseMsh->_elementOffset[isdom]; iel < _coarseMsh->_elementOffset[isdom + 1]; iel++) {
          const short unsigned ielt = _coarseMsh->GetElementType(iel);
            Get_Prolongation_SparsityPatternSize_OneElement_OneFEFamily(*this, *_coarseMsh, iel, NNZ_d, NNZ_o, el_dofs, GetFiniteElement(ielt, solType) );
        }
      }

      NNZ_d->close();
      NNZ_o->close();

      const unsigned offset = _dofOffset[solType][_iproc];
      std::vector <int> nnz_d(nf_loc);
      std::vector <int> nnz_o(nf_loc);

      for(int i = 0; i < nf_loc; i++) {
        nnz_d[i] = static_cast <int>(floor((*NNZ_d)(offset + i) + 0.5));
        nnz_o[i] = static_cast <int>(floor((*NNZ_o)(offset + i) + 0.5));
      }

      delete NNZ_d;
      delete NNZ_o;
    // ------------------- Sparsity pattern size - END


      
    // ------------------- Prolongator - BEGIN
      //build matrix
      _ProjCoarseToFine[solType] = SparseMatrix::build().release();
      _ProjCoarseToFine[solType]->init(nf, nc, nf_loc, nc_loc, nnz_d, nnz_o);

      // loop on the coarse grid
      for(int isdom = _iproc; isdom < _iproc + 1; isdom++) {
        for(int iel = _coarseMsh->_elementOffset[isdom]; iel < _coarseMsh->_elementOffset[isdom + 1]; iel++) {
          const short unsigned ielt = _coarseMsh->GetElementType(iel);
            Build_Prolongation_OneElement_OneFEFamily(*this, *_coarseMsh, iel, _ProjCoarseToFine[solType], el_dofs, GetFiniteElement(ielt, solType) );
        }
      }

      _ProjCoarseToFine[solType]->close();
    // ------------------- Prolongator - END
      
    }

  }

  
  
//----------------------------------------------------------------------------------------------------
//BEGIN  build matrix sparsity pattern size and build prolongator matrix for single solution
//-----------------------------------------------------------------------------------------------------

  void Mesh::Get_Prolongation_SparsityPatternSize_OneElement_OneFEFamily(const Mesh& meshf,
                                         const Mesh& meshc,
                                         const int& ielc,
                                         NumericVector* NNZ_d,
                                         NumericVector* NNZ_o,
                                         const char is_fine_or_coarse[],
                                         const elem_type * elem_type_in) const
  {

    const unsigned soltype_in = elem_type_in->GetSolType();
    const unsigned      ndofs = elem_type_in->GetNDofs();
    const unsigned      ndofs_fine = elem_type_in->GetNDofsFine();
    
      unsigned n_elemdofs = 0;
      if ( !strcmp(is_fine_or_coarse, "fine") )        n_elemdofs = ndofs_fine;
      else if ( !strcmp(is_fine_or_coarse, "coarse") ) n_elemdofs = ndofs;
      
    if(meshc.GetRefinedElementIndex(ielc)) {  // coarse2fine prolongation
      
      for(int i = 0; i < n_elemdofs ; i++) {
        
        const std::pair<int, int> id_0_1 = elem_type_in->GetKVERT_IND(i);
        
        const int irow = meshf.GetSolutionDof(ielc, id_0_1.first, id_0_1.second, soltype_in, &meshc);

        const int iproc = meshf.BisectionSearch_find_processor_of_dof(irow, soltype_in);
        
        const int ncols = elem_type_in->Get_Prolongator_Num_Columns(i);
        
        unsigned counter_off_diag = 0;

        for(int k = 0; k < ncols; k++) {
          int j = elem_type_in->Get_Prolongator_Index(i, k);
          int jcolumn = meshc.GetSolutionDof(j, ielc, soltype_in);

          if(jcolumn <  meshc.dofmap_get_dof_offset(soltype_in, iproc)     ||
             jcolumn >= meshc.dofmap_get_dof_offset(soltype_in, iproc + 1)   ) counter_off_diag++;
        }

        NNZ_d->set(irow, ncols - counter_off_diag);
        NNZ_o->set(irow, counter_off_diag);
      }
      
    }
    else { // coarse2coarse prolongation
      
      for(int i = 0; i < ndofs; i++) {
        
        int irow = meshf.GetSolutionDof(ielc, 0, i , soltype_in, &meshc);

        int iproc = meshf.BisectionSearch_find_processor_of_dof(irow, soltype_in);
        int jcolumn = meshc.GetSolutionDof(i, ielc, soltype_in);

        if( jcolumn <  meshc.dofmap_get_dof_offset(soltype_in, iproc)      || 
            jcolumn >= meshc.dofmap_get_dof_offset(soltype_in, iproc + 1)     ) {
          NNZ_o->set(irow, 1);
        }
        else {
          NNZ_d->set(irow, 1);
        }
      }
      
    }
    
    
  }
  
  

  void Mesh::Build_Prolongation_OneElement_OneFEFamily(const Mesh& meshf,
                                    const Mesh& meshc, 
                                    const int& ielc,
                                    SparseMatrix* Projmat, 
                                    const char is_fine_or_coarse[],
                                  const elem_type * elem_type_in) const
  {
 
      const unsigned soltype_in = elem_type_in->GetSolType();
      const unsigned      ndofs = elem_type_in->GetNDofs();
      const unsigned      ndofs_fine = elem_type_in->GetNDofsFine();
    
      unsigned n_elemdofs = 0;
      if ( !strcmp(is_fine_or_coarse, "fine") )        n_elemdofs = ndofs_fine;
      else if ( !strcmp(is_fine_or_coarse, "coarse") ) n_elemdofs = ndofs;
      
      if(meshc.GetRefinedElementIndex(ielc)) {  // coarse2fine prolongation

      std::vector<int> jcols( ndofs );

      for(int i = 0; i < n_elemdofs /*_nf*/; i++) {
        
        const std::pair<int, int> id_0_1 = elem_type_in->GetKVERT_IND(i);
                
        const int irow = meshf.GetSolutionDof(ielc, id_0_1.first, id_0_1.second, soltype_in, &meshc);
        const int ncols  =  elem_type_in->Get_Prolongator_Num_Columns(i);
        
        jcols.assign(ncols, 0);

        for(int k = 0; k < ncols; k++) {
          int j = elem_type_in->Get_Prolongator_Index(i, k);
          int jcolumn = meshc.GetSolutionDof(j, ielc, soltype_in);
          jcols[k] = jcolumn;
        }

        Projmat->insert_row(irow, ncols, jcols, elem_type_in->Get_Prolongator_Values_Row(i) );
      }
      
    }
    else { // coarse2coarse prolongation
      
      std::vector <int> jcol(1);
      double one = 1.;

      for(int i = 0; i < ndofs; i++) {
        
        const int irow = meshf.GetSolutionDof(ielc, 0, i , soltype_in, &meshc);
        jcol[0] = meshc.GetSolutionDof(i, ielc, soltype_in);
        Projmat->insert_row(irow, 1, jcol, &one);
        
      }
      
    }
    
    
  }

//----------------------------------------------------------------------------------------------------
//END  build matrix sparsity pattern size and build prolongator matrix for single solution
//-----------------------------------------------------------------------------------------------------

  
  

  short unsigned Mesh::GetRefinedElementIndex(const unsigned& iel) const {
    return static_cast <short unsigned>((*_topology->_Sol[_amrIndex])(iel) + 0.25);
  }

  short unsigned Mesh::GetElementGroup(const unsigned int& iel) const {
    return el->GetElementGroup(iel);
  }

  short unsigned Mesh::GetElementMaterial(const unsigned int& iel) const {
    return el->GetElementMaterial(iel);
  }

  short unsigned Mesh::GetElementType(const unsigned int& iel) const {
    return el->GetElementType(iel);
  }

  bool Mesh::GetSolidMark(const unsigned int& inode) const {
    return static_cast <short unsigned>((*_topology->_Sol[_solidMarkIndex])(inode) + 0.25);
  }


  /** Only for parallel */
  unsigned Mesh::GetElementDofNumber(const unsigned& iel, const unsigned& type) const {
    return el->GetNVE(GetElementType(iel), type);
  }

  /** Only for parallel */
  const unsigned Mesh::GetElementFaceType(const unsigned& kel, const unsigned& jface) const {
    unsigned kelt = GetElementType(kel);
    const unsigned FELT[N_GEOM_ELS][2] = {{3, 3}, {4, 4}, {3, 4}, {5, 5}, {5, 5}, {6, 6}};
    const unsigned felt = FELT[kelt][jface >= GetElementFaceNumber(kel, 0)];
    return felt;
  }

  /** Only for parallel */
  unsigned Mesh::GetLocalFaceVertexIndex(const unsigned& iel, const unsigned& iface, const unsigned& jnode) const {
    return el->GetIG(GetElementType(iel), iface, jnode);
  }

  /** this needs the abstract elem type */
  unsigned Mesh::GetLocalFaceVertexIndex_PassElemType(const short unsigned & el_type, const unsigned& iface, const unsigned& jnode) const {
    return el->GetIG(el_type, iface, jnode);
  }


  /** Only for parallel */
  unsigned Mesh::GetElementFaceDofNumber(const unsigned& iel, const unsigned jface, const unsigned& type) const {
    assert(type < NFE_FAMS_C_ZERO_LAGRANGE);   ///@todo relax this
    return el->GetNFACENODES(GetElementType(iel), jface, type);
  }

  /** Only for parallel */
  unsigned Mesh::GetElementFaceNumber(const unsigned& iel, const unsigned& type) const {
    return el->GetNFC(GetElementType(iel), type);
  }
  
  /** this needs the abstract elem type */
  unsigned Mesh::GetElementFaceNumber_PassElemType(const short unsigned & el_type, const unsigned& type) const {
    return el->GetNFC(el_type, type);
  }

// *******************************************************

  void Mesh::AddBiquadraticNodesNotInMeshFile() {

    unsigned int nnodes = GetNumberOfNodes();

    //intialize to UINT_MAX
    for(unsigned iel = 0; iel < el->GetElementNumber(); iel++) {
      const unsigned elementType = el->GetElementType(iel);

      if(elementType == 1 || elementType == 2) {
        for(unsigned inode = el->GetElementDofNumber(iel, 2) - _numberOfMissedBiquadraticNodes[elementType];
            inode < el->GetElementDofNumber(iel, 2); inode++) {
          el->SetElementDofIndex(iel, inode, UINT_MAX);
        }
      }
    }

    // generate face dofs for tet and wedge elements
    for(unsigned iel = 0; iel < el->GetElementNumber(); iel++) {
      const unsigned elementType = el->GetElementType(iel);

      if(elementType == 1 || elementType == 2) {
        for(unsigned iface = el->GetElementFaceNumber(iel, 0); iface < el->GetElementFaceNumber(iel, 1); iface++) {       //on all the faces that are triangles
          const unsigned inode = el->GetElementDofNumber(iel, 1) + iface;

          if(UINT_MAX == el->GetElementDofIndex(iel, inode)) {
            el->SetElementDofIndex(iel, inode, nnodes);
            const unsigned i1 = el->GetFaceVertexIndex(iel, iface, 0);
            const unsigned i2 = el->GetFaceVertexIndex(iel, iface, 1);
            const unsigned i3 = el->GetFaceVertexIndex(iel, iface, 2);
            bool faceHasBeenFound = false;

            for(unsigned jel = iel + 1; jel < el->GetElementNumber(); jel++) {
              for(unsigned jface = el->GetElementFaceNumber(jel, 0); jface < el->GetElementFaceNumber(jel, 1); jface++) {
                unsigned jnode = el->GetElementDofNumber(jel, 1) + jface;

                if(UINT_MAX == el->GetElementDofIndex(jel, jnode)) {
                  const unsigned j1 = el->GetFaceVertexIndex(jel, jface, 0);
                  const unsigned j2 = el->GetFaceVertexIndex(jel, jface, 1);
                  const unsigned j3 = el->GetFaceVertexIndex(jel, jface, 2);

                  if((i1 == j1 || i1 == j2 || i1 == j3) &&
                      (i2 == j1 || i2 == j2 || i2 == j3) &&
                      (i3 == j1 || i3 == j2 || i3 == j3)) {
                    el->SetElementDofIndex(jel, jnode, nnodes);
                    faceHasBeenFound = true;
                    break;
                  }
                }
              }

              if(faceHasBeenFound) {
                break;
              }
            }

            ++nnodes;
          }
        }
      }
    }

    // generates element dofs for tet, wedge and triangle elements
    for(unsigned iel = 0; iel < el->GetElementNumber(); iel++) {
      if(1 == el->GetElementType(iel)) {     //tet
        el->SetElementDofIndex(iel, 14, nnodes);
        ++nnodes;
      }
      else if(2 == el->GetElementType(iel)) {     //wedge
        el->SetElementDofIndex(iel, 20, nnodes);
        ++nnodes;
      }
      else if(4 == el->GetElementType(iel)) {     //triangle
        el->SetElementDofIndex(iel, 6, nnodes);
        ++nnodes;
      }
    }

    el->SetNodeNumber(nnodes);
    SetNumberOfNodes(nnodes);
//     std::cout <<"nnodes after="<< nnodes << std::endl;
    
    

    // add the coordinates of the biquadratic nodes not included in gambit
    _coords[0].resize(nnodes);
    _coords[1].resize(nnodes);
    _coords[2].resize(nnodes);

    for(int iel = 0; iel < el->GetElementNumber(); iel++) {
      unsigned elementType = el->GetElementType(iel);

      unsigned ndof = el->GetElementDofNumber(iel, 2);
      unsigned jstart = ndof - _numberOfMissedBiquadraticNodes[elementType];

      for(unsigned j = jstart; j < ndof; j++) {

        unsigned jnode = el->GetElementDofIndex(iel, j);

        _coords[0][jnode] = 0.;
        _coords[1][jnode] = 0.;
        _coords[2][jnode] = 0.;

        for(int i = 0; i < jstart; i++) {
          unsigned inode = el->GetElementDofIndex(iel, i);

          for(int k = 0; k < GetDimension(); k++) {
            _coords[k][jnode] += _coords[k][inode] * _baricentricWeight[ elementType ][j - jstart][i];
          }
        }
      }
    }
    
    
  }
  
  

  basis* Mesh::GetBasis(const short unsigned& ielType, const short unsigned& solType)  const {
    return _finiteElement[ielType][solType]->GetBasis();
  }

  
  void Mesh::GetElementNodeCoordinates(std::vector < std::vector <double > > &xv, const unsigned &iel, const unsigned &solType) const {
      
    xv.resize(_dimension);
    unsigned ndofs = el->GetElementDofNumber(iel, solType);
    for(int d = 0; d < _dimension; d++) {
      xv[d].resize(ndofs);
    }
    for(unsigned j = 0; j < ndofs; j++) {
      unsigned xdof  = GetSolutionDof(j, iel, solType);
      for(int d = 0; d < _dimension; d++) {
        xv[d][j] = (*_topology->_Sol[d])(xdof);
      }
    }
    
  }

  
  
  
} //end namespace femus
