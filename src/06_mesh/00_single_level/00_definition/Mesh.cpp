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


#include "PolynomialBases.hpp"


// C++ includes
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <algorithm>


namespace femus {

  
bool (* Mesh::_SetRefinementFlag)(const std::vector < double >& x, const int &ElemGroupNumber,const int &level) = NULL; 

    const std::string Mesh::_x_name = "X";
    const std::string Mesh::_y_name = "Y";  
    const std::string Mesh::_z_name = "Z";  
    const std::string Mesh::_amrIndex_name = "AMR";
    const std::string Mesh::_solidMark_name = "solidMrk";
    
  bool Mesh::_IsUserRefinementFunctionDefined = false;


// === Constructors / Destructor - BEGIN =================

//------------------------------------------------------------------------------------------------------
  Mesh::Mesh() {

  }


  Mesh::~Mesh() {
      
    delete el;
    
    _topology->FreeSolutionVectors();
    delete _topology;


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


    Initialize_Level_AMR_Coords();


    ReadCoarseMeshFile(name, Lref, type_elem_flag, read_groups, read_boundary_groups);


    AddNodes_in_CoarseGen();


  }  
  

  /**
   *  This function generates the coarse Mesh level, $l_0$, from an input Mesh file
   **/
  void Mesh::ReadCoarseMesh(const std::string& name, const double Lref, std::vector<bool>& type_elem_flag, const bool read_groups, const bool read_boundary_groups) {
      

    ReadCoarseMeshBeforePartitioning(name, Lref, type_elem_flag, read_groups, read_boundary_groups);
    
    PartitionElements_and_FillDofMapAllFEFamilies();

    BuildElementAndNodeStructures();
   
    
    
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


//=========


//====  CharacteristicLength ========
    ComputeCharacteristicLength();  //doesn't need dofmap

//====  Print Info ========
    PrintInfo();  //needs dofmap

    
  }

  
  void Mesh::InitializeAndPossiblyFillAmrRestriction(const bool amr) {

    if(amr) {
      
      GetAMRRestrictionAndAMRSolidMark(this, GetMeshElements() );
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

    _topology->AddSolution(_x_name, LAGRANGE, SECOND, 1, 0);
    _topology->AddSolution(_y_name, LAGRANGE, SECOND, 1, 0);
    _topology->AddSolution(_z_name, LAGRANGE, SECOND, 1, 0);

    _topology->ResizeSolutionVector(_x_name);  //needs dofmap
    _topology->ResizeSolutionVector(_y_name);  //needs dofmap
    _topology->ResizeSolutionVector(_z_name);  //needs dofmap

  }
  

  void Mesh::Topology_FillCoordinates() {
    
    //set coordinates -----------
    _topology->GetSolutionByName(_x_name) = _coords[0];
    _topology->GetSolutionByName(_y_name) = _coords[1];
    _topology->GetSolutionByName(_z_name) = _coords[2];
    //set coordinates -----------

  }
  
  
  
  
/// This needs the dof maps, for the discontinuous Lagrange elements
  void Mesh::Topology_InitializeAMR() {
      
    _topology->AddSolution(_amrIndex_name, DISCONTINUOUS_POLYNOMIAL, ZERO, 1, 0);

    _topology->ResizeSolutionVector(_amrIndex_name);

  }
  

/// This needs the dof maps, for the continuous Lagrange elements
  void Mesh::Topology_InitializeSolidNodeFlag() {
      
    _topology->AddSolution(_solidMark_name, LAGRANGE, SECOND, 1, 0);
    
    _topology->ResizeSolutionVector(_solidMark_name);

  }
  

  
  void Mesh::ComputeCharacteristicLength() {
      
    //compute max and min coords - BEGIN -----------
    std::vector < double > xMax(3, 0.);
    std::vector < double > xMin(3, 0.);
    for(unsigned i = 0; i < _coords[0].size(); i++) {
      for(unsigned k = 0; k < 3; k++) {
        if(xMax[k] < _coords[k][i]) xMax[k] = _coords[k][i];
        if(xMin[k] > _coords[k][i]) xMin[k] = _coords[k][i];
      }
    }
    //compute max and min coords - END -----------
    
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


    GenerateCoarseBoxMeshBeforePartitioning(nx, ny, nz,
    xmin, xmax,
    ymin, ymax,
    zmin, zmax,
    elemType, type_elem_flag);



    PartitionElements_and_FillDofMapAllFEFamilies();

    BuildElementAndNodeStructures();
    
  }



  void Mesh::Initialize_Level_AMR_Coords() {

//==== Level - BEGIN ==============================
    SetLevel(0);
//==== Level - END ==============================

//==== AMR - BEGIN ==============================
    SetIfHomogeneous(true);
//==== AMR - END ==============================

//==== Coords, coarse, then to topology - BEGIN ==============================
    _coords.resize(3);
//==== Coords, coarse, then to topology - END ==============================

  }



  void Mesh::GenerateCoarseBoxMeshBeforePartitioning(
    const unsigned int nx, const unsigned int ny, const unsigned int nz,
    const double xmin, const double xmax,
    const double ymin, const double ymax,
    const double zmin, const double zmax,
    const ElemType elemType, std::vector<bool>& type_elem_flag) {


    Initialize_Level_AMR_Coords();


    MeshTools::Generation::BuildBox(*this, _coords, nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, elemType, type_elem_flag);


    AddNodes_in_CoarseGen();


  }



  void Mesh::AddNodes_in_CoarseGen() {

//==== Only in Coarse generations - BEGIN ==============================
    AddBiquadraticNodesNotInMeshFile();

    el->ShrinkToFitElementDof();
    el->ShrinkToFitElementNearFace();
//==== Only in Coarse generations - END ==============================

  }



  void Mesh::Topology_FillSolidNodeFlag() {


    NumericVector& NodeMaterial =  _topology->GetSolutionByName(_solidMark_name);

    NodeMaterial.zero();

    for(int iel = GetElementOffset(_iproc); iel < GetElementOffset(_iproc + 1); iel++) {
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
        for(unsigned iel = GetElementOffset(isdom); iel < GetElementOffset(isdom + 1); iel++) {
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
//       for (unsigned i = GetElementOffset(isdom); i < GetElementOffset(isdom + 1) - 1; i++) {
//         unsigned iel = inverse_element_mapping[i];
//         unsigned ielMat = el->GetElementMaterial (iel);
//         unsigned ielGroup = el->GetElementGroup (iel);
//         for (unsigned j = i + 1; j < GetElementOffset(isdom + 1); j++) {
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

      unsigned n = _elementOffset[isdom + 1u] - GetElementOffset(isdom);
      
      while(n > 1) {
        unsigned newN = 0u;
        for(unsigned j = GetElementOffset(isdom) + 1u; j < GetElementOffset(isdom) + n ; j++) {
          jel = inverse_element_mapping[j];
          jelMat = el->GetElementMaterial(jel);
          jelGroup = el->GetElementGroup(jel);

          iel = inverse_element_mapping[j - 1];
          ielMat = el->GetElementMaterial(iel);
          ielGroup = el->GetElementGroup(iel);

          if( jelMat < ielMat || (jelMat == ielMat && (jelGroup < ielGroup || (jelGroup == ielGroup && jel < iel) ) ) ) {
            inverse_element_mapping[j - 1] = jel;
            inverse_element_mapping[j] = iel;
            newN = j - GetElementOffset(isdom);
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
      _ownSize[3][isdom] = GetElementOffset(isdom + 1) - GetElementOffset(isdom);
      _ownSize[4][isdom] = (GetElementOffset(isdom + 1) - GetElementOffset(isdom)) * (_dimension + 1);
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

        for(unsigned iel = GetElementOffset(isdom); iel < GetElementOffset(isdom + 1); iel++) {
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
        const unsigned iNode = el->GetElementDofIndex(iel, i);  //GetMeshDof(iel, i, solType);
        const unsigned isdom = BisectionSearch_find_processor_of_dof(iNode, CONTINUOUS_BIQUADRATIC);

        if(iNode < _dofOffset[CONTINUOUS_BIQUADRATIC][isdom] + _originalOwnSize[CONTINUOUS_LINEAR][isdom]) {
          dof = (iNode - _dofOffset[CONTINUOUS_BIQUADRATIC][isdom]) + _dofOffset[CONTINUOUS_LINEAR][isdom];
        }
        else {
          dof = _ownedGhostMap[CONTINUOUS_LINEAR].find(iNode)->second;
        }
      }
      break;

      case CONTINUOUS_SERENDIPITY: { // quadratic Lagrange
        const unsigned iNode = el->GetElementDofIndex(iel, i);  //GetMeshDof(iel, i, solType);
        const unsigned isdom = BisectionSearch_find_processor_of_dof(iNode, CONTINUOUS_BIQUADRATIC);

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
        const unsigned isdom = BisectionSearch_find_processor_of_dof(iel, DISCONTINUOUS_CONSTANT);
        const unsigned offset = GetElementOffset(isdom);
        const unsigned offsetp1 = GetElementOffset(isdom + 1);
        const unsigned ownSize = offsetp1 - offset;
        const unsigned offsetPWLD = offset * (_dimension + 1);
        const unsigned locIel = iel - offset;
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
        unsigned iNode = mshc->GetMeshElements()->GetChildElementDof(ielc, i0, i1);
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
        unsigned iNode = mshc->GetMeshElements()->GetChildElementDof(ielc, i0, i1);
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
        dof = mshc->GetMeshElements()->GetChildElementDof(ielc, i0, i1);
        break;

      case DISCONTINUOUS_CONSTANT: // piecewise constant
        // in this case use i=0
        dof = mshc->GetMeshElements()->GetChildElement(ielc, i0);
        break;

      case DISCONTINUOUS_LINEAR: // piecewise linear discontinuous
        unsigned iel = mshc->GetMeshElements()->GetChildElement(ielc, i0);
        unsigned isdom = BisectionSearch_find_processor_of_dof(iel, DISCONTINUOUS_CONSTANT);
        unsigned offset = GetElementOffset(isdom);
        unsigned offsetp1 = GetElementOffset(isdom + 1);
        unsigned ownSize = offsetp1 - offset;
        unsigned offsetPWLD = offset * (_dimension + 1);
        unsigned locIel = iel - offset;
        dof = offsetPWLD + (i1 * ownSize) + locIel;
        break;
    }

    return dof;
  }

// *******************************************************

  
  
  


  
  

  short unsigned Mesh::GetRefinedElementIndex(const unsigned& iel) const {
    return static_cast <short unsigned>((*GetTopology()->_Sol[_amrIndex])(iel) + 0.25);
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
    return static_cast <short unsigned>((*GetTopology()->_Sol[_solidMarkIndex])(inode) + 0.25);
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

    // Nodes, initialize - BEGIN
    unsigned int nnodes = GetNumberOfNodes();
    // Nodes, initialize - END


    //intialize to UINT_MAX - BEGIN
    for(unsigned iel = 0; iel < el->GetElementNumber(); iel++) {
      const unsigned elementType = el->GetElementType(iel);

      if(elementType == 1 || elementType == 2) {
        for(unsigned inode = el->GetElementDofNumber(iel, 2) - _numberOfMissedBiquadraticNodes[elementType];
            inode < el->GetElementDofNumber(iel, 2); inode++) {
          el->SetElementDofIndex(iel, inode, UINT_MAX);
        }
      }
    }
    //intialize to UINT_MAX - END


    // generate face dofs for tet and wedge elements - BEGIN
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
    // generate face dofs for tet and wedge elements - END

    // generates element dofs for tet, wedge and triangle elements - BEGIN
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
    // generates element dofs for tet, wedge and triangle elements - END


    // Nodes - BEGIN
    el->SetNodeNumber(nnodes);
    SetNumberOfNodes(nnodes);
    // Nodes - END

    

    // add the coordinates of the biquadratic nodes not included in the mesh file - BEGIN
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
    // add the coordinates of the biquadratic nodes not included in the mesh file - END

    
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
        xv[d][j] = (*GetTopology()->_Sol[d])(xdof);
      }
    }
    
  }

  

  
  void Mesh::GetAMRRestrictionAndAMRSolidMark(Mesh *msh, const elem *el_in)
  {

    std::vector < std::map < unsigned,  std::map < unsigned, double  > > > & restriction = msh->GetAmrRestrictionMap();
    restriction.resize(NFE_FAMS_C_ZERO_LAGRANGE);

    std::vector < std::map < unsigned, bool > > & interfaceSolidMark = msh->GetAmrSolidMark();
    interfaceSolidMark.resize(NFE_FAMS_C_ZERO_LAGRANGE);

    std::vector < MyVector<unsigned> > interfaceElement;
    std::vector < MyMatrix<unsigned> > interfaceLocalDof;
    std::vector < std::vector < MyMatrix<unsigned> > > interfaceDof;
    std::vector < std::vector < MyMatrix<unsigned> > > levelInterfaceSolidMark;
    std::vector < std::vector < MyMatrix< double > > > interfaceNodeCoordinates;

    const unsigned level_of_list = el_in->GetLevelOfRefinementForList();
    const unsigned curr_proc = msh->processor_id();
    const unsigned n_procs = msh->n_processors();
    
    interfaceElement.resize(level_of_list + 1);
    interfaceLocalDof.resize(level_of_list + 1);
    interfaceDof.resize(NFE_FAMS_C_ZERO_LAGRANGE);
    levelInterfaceSolidMark.resize(NFE_FAMS_C_ZERO_LAGRANGE);
    for (unsigned i = 0; i < NFE_FAMS_C_ZERO_LAGRANGE; i++) {
      interfaceDof[i].resize(level_of_list + 1);
      levelInterfaceSolidMark[i].resize(level_of_list + 1);
    }
    interfaceNodeCoordinates.resize(level_of_list + 1);
    const unsigned dim = msh->GetDimension();


    for (unsigned ilevel = 0; ilevel <= level_of_list; ilevel++) {
      
      //BEGIN interface element search
      interfaceElement[ilevel] = MyVector <unsigned> ( el_in->get_n_elements_owned() );
      unsigned counter = 0;
      for (unsigned i = el_in->GetElementLevelArray().begin(); i < el_in->GetElementLevelArray().end(); i++) {
        if (ilevel == el_in->GetElementLevel(i) ) {
          for (unsigned j = el_in->GetElementNearFaceArray().begin(i); j < el_in->GetElementNearFaceArray().end(i); j++) {
            if (-1 == el_in->GetElementNearFaceArray()[i][j] ) {
              interfaceElement[ilevel][counter] = i;
              counter++;
              break;
            }
          }
        }
      }
      interfaceElement[ilevel].resize(counter);
      interfaceElement[ilevel].stack();
      //END interface element search

      //BEGIN interface node search
      std::vector< unsigned > offset = interfaceElement[ilevel].getOffset();
      interfaceLocalDof[ilevel] = MyMatrix <unsigned>(offset, el_in->GetNVE(0, 2), UINT_MAX);
      for (unsigned i = interfaceElement[ilevel].begin(); i < interfaceElement[ilevel].end(); i++) {
        unsigned iel =  interfaceElement[ilevel][i];
        std::map <unsigned, bool> ldofs;
        for (unsigned jface = el_in->GetElementNearFaceArray().begin(iel); jface < el_in->GetElementNearFaceArray().end(iel); jface++) {
          if (-1 == el_in->GetElementNearFaceArray()[iel][jface]) {
            for (unsigned k = 0; k < el_in->GetNFACENODES( el_in->GetElementType(iel), jface, 2); k++) {
              const unsigned index = el_in->GetIG( el_in->GetElementType(iel), jface, k);
              ldofs[index] = true;
            }
          }
        }
        unsigned j = 0;
        for (std::map<unsigned, bool>::iterator it = ldofs.begin(); it != ldofs.end(); it++) {
          interfaceLocalDof[ilevel][i][j] = it->first;
          j++;
        }
      }
      interfaceLocalDof[ilevel].shrinkToFit(UINT_MAX);
      //END interface node search

      //BEGIN interface node dof global search, one for each soltype
      MyVector <unsigned> rowSize = interfaceLocalDof[ilevel].getRowSize();
      for (unsigned soltype = 0; soltype < NFE_FAMS_C_ZERO_LAGRANGE; soltype++) {
        interfaceDof[soltype][ilevel] = MyMatrix< unsigned > (rowSize, UINT_MAX);
        levelInterfaceSolidMark[soltype][ilevel] = MyMatrix< unsigned > (rowSize, UINT_MAX);
        for (unsigned i = interfaceLocalDof[ilevel].begin(); i < interfaceLocalDof[ilevel].end(); i++) {
          unsigned iel = interfaceElement[ilevel][i];
          unsigned counter = 0;
          for (unsigned j = interfaceLocalDof[ilevel].begin(i); j < interfaceLocalDof[ilevel].end(i); j++) {
            unsigned jloc = interfaceLocalDof[ilevel][i][j];
            if (jloc < el_in->GetElementDofNumber(iel, soltype)) {
              unsigned jdof  = msh->GetSolutionDof(jloc, iel, soltype);
              interfaceDof[soltype][ilevel][i][counter] = jdof;
              unsigned jdof2  = msh->GetSolutionDof(jloc, iel, 2);
              levelInterfaceSolidMark[soltype][ilevel][i][counter] = msh->GetSolidMark(jdof2);
              counter++;
            }
            else {
              break;
            }
          }
        }
        interfaceDof[soltype][ilevel].shrinkToFit(UINT_MAX);
        levelInterfaceSolidMark[soltype][ilevel].shrinkToFit(UINT_MAX);
      }
      //END interface node dof global search, one for each soltype

      //BEGIN interface node coordinates
      interfaceNodeCoordinates[ilevel].resize(dim);
      for (unsigned k = 0; k < dim; k++) {
        interfaceNodeCoordinates[ilevel][k] = MyMatrix< double > (rowSize, 0.);
      }
      for (unsigned i = interfaceLocalDof[ilevel].begin(); i < interfaceLocalDof[ilevel].end(); i++) {
        unsigned iel = interfaceElement[ilevel][i];
        for (unsigned j = interfaceLocalDof[ilevel].begin(i); j < interfaceLocalDof[ilevel].end(i); j++) {
          unsigned jnode = interfaceLocalDof[ilevel][i][j];
          unsigned xDof  = msh->GetSolutionDof(jnode, iel, 2);
          for (unsigned k = 0; k < dim; k++) {
            interfaceNodeCoordinates[ilevel][k][i][j] = (*msh->GetTopology()->_Sol[k])(xDof);
          }
        }
      }
      //END interface node coordinates

    }


    for (unsigned soltype = 0; soltype < NFE_FAMS_C_ZERO_LAGRANGE; soltype++) {
      for (int ilevel = 0; ilevel < level_of_list; ilevel++) {
        for (int jlevel = ilevel + 1; jlevel <= level_of_list; jlevel++) {
          for (unsigned lproc = 0; lproc < n_procs; lproc++) {
            interfaceDof[soltype][jlevel].broadcast(lproc);
            levelInterfaceSolidMark[soltype][jlevel].broadcast(lproc);
            for (unsigned d = 0; d < dim; d++) {
              interfaceNodeCoordinates[jlevel][d].broadcast(lproc);
            }
            std::map< unsigned, bool> candidateNodes;
            std::map< unsigned, bool> elementNodes;

            for (unsigned i = interfaceDof[soltype][ilevel].begin(); i < interfaceDof[soltype][ilevel].end(); i++) {

              candidateNodes.clear();

              std::vector < std::vector < std::vector <double > > > aP(NFE_FAMS_C_ZERO_LAGRANGE);
              bool aPIsInitialized = false;

              unsigned iel = interfaceElement[ilevel][i];
              short unsigned ielType = el_in->GetElementType(iel);

              elementNodes.clear();
              for (unsigned j = 0; j < el_in->GetElementDofNumber(iel, soltype); j++) {
                unsigned jdof  = msh->GetSolutionDof(j, iel, soltype);
                elementNodes[jdof] = true;
              }

              std::vector < std::vector <double > > xv;
              msh->GetElementNodeCoordinates(xv, iel, CONTINUOUS_BIQUADRATIC);
              unsigned ndofs = xv[0].size();

              double r;
              std::vector <double> xc;
              GetConvexHullSphere(xv, xc, r, 0.01);
              double r2 = r * r;

              std::vector < std::vector< double > > xe;
              GetBoundingBox(xv, xe, 0.01);


              for (unsigned k = interfaceDof[soltype][jlevel].begin(); k < interfaceDof[soltype][jlevel].end(); k++) {
                for (unsigned l = interfaceDof[soltype][jlevel].begin(k); l < interfaceDof[soltype][jlevel].end(k); l++) {
                  unsigned ldof = interfaceDof[soltype][jlevel][k][l];
                  if (candidateNodes.find(ldof) == candidateNodes.end() || candidateNodes[ldof] != false) {
                    double d2 = 0.;
                    std::vector<double> xl(dim);
                    for (int d = 0; d < dim; d++) {
                      xl[d] = interfaceNodeCoordinates[jlevel][d][k][l];
                      d2 += (xl[d] - xc[d]) * (xl[d] - xc[d]);
                    }
                    bool insideHull = true;
                    if (d2 > r2) {
                      insideHull = false;
                    }
                    for (unsigned d = 0; d < dim; d++) {
                      if (xl[d] < xe[d][0] || xl[d] > xe[d][1]) {
                        insideHull = false;
                      }
                    }
                    if (insideHull) {
                      if (elementNodes.find(ldof) == elementNodes.end()) {

                        if (!aPIsInitialized) {
                          aPIsInitialized = true;
                          std::vector < std::vector <double> > x1(dim);
                          for (unsigned jtype = 0; jtype < NFE_FAMS_C_ZERO_LAGRANGE; jtype++) {
                            ProjectNodalToPolynomialCoefficients(aP[jtype], xv, ielType, jtype) ;
                          }
                        }

                        std::vector <double> xi;
                        GetClosestPointInReferenceElement(xv, xl, ielType, xi);
                        GetInverseMapping(2, ielType, aP, xl, xi);

                        bool insideDomain = CheckIfPointIsInsideReferenceDomain(xi, ielType, 0.0001);
                        if (insideDomain) {
                          for (unsigned j = interfaceDof[soltype][ilevel].begin(i); j < interfaceDof[soltype][ilevel].end(i); j++) {
                            unsigned jloc = interfaceLocalDof[ilevel][i][j];

                            const basis* base = msh->GetFiniteElement(ielType, soltype)->GetBasis();
                            double value = base->eval_phi(jloc, xi);

                            if (fabs(value) >= 1.0e-10) {
                              unsigned jdof = interfaceDof[soltype][ilevel][i][j];
                              if (restriction[soltype][jdof].find(jdof) == restriction[soltype][jdof].end()) {
                                restriction[soltype][jdof][jdof] = 1.;
                                unsigned jdof2  = msh->GetSolutionDof(jloc, iel, 2);
                                interfaceSolidMark[soltype][jdof] = levelInterfaceSolidMark[soltype][ilevel][i][j];
                              }
                              restriction[soltype][jdof][ldof] = value;
                              restriction[soltype][ldof][ldof] = 10.;
                              interfaceSolidMark[soltype][ldof] = levelInterfaceSolidMark[soltype][jlevel][k][l];
                              candidateNodes[ldof] = true;
                            }
                          }
                        }
                        else {
                          candidateNodes[ldof] = false;
                        }
                      }
                      else {
                        candidateNodes[ldof] = false;
                      }
                    }
                  }
                }
              }
            }
            interfaceDof[soltype][jlevel].clearBroadcast();
            levelInterfaceSolidMark[soltype][jlevel].clearBroadcast();
            for (unsigned d = 0; d < dim; d++) {
              interfaceNodeCoordinates[jlevel][d].clearBroadcast();
            }
          }
        }
      }




      NumericVector* pvector;
      pvector = NumericVector::build().release();
      pvector->init(n_procs, 1 , false, AUTOMATIC);

      unsigned counter = 1;
      while (counter != 0) {
        counter = 0;

        //BEGIN  saving the restriction object in parallel vectors and matrices

        MyVector <unsigned> rowSize(restriction[soltype].size(), 0);
        unsigned cnt1 = 0;
        for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
          rowSize[cnt1] = restriction[soltype][it1->first].size();
          cnt1++;
        }
        rowSize.stack();

        std::vector< unsigned > offset = rowSize.getOffset();

        MyVector <unsigned> masterNode(offset);
        MyMatrix <unsigned> slaveNodes(rowSize);
        MyMatrix <double> slaveNodesValues(rowSize);

        cnt1 = 0;
        for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
          masterNode[offset[curr_proc] + cnt1] = it1->first;
          unsigned cnt2 = 0;
          for (std::map<unsigned, double> ::iterator it2 = restriction[soltype][it1->first].begin(); it2 != restriction[soltype][it1->first].end(); it2++) {
            slaveNodes[offset[curr_proc] + cnt1][cnt2] = it2->first;
            slaveNodesValues[offset[curr_proc] + cnt1][cnt2] = it2->second;
            cnt2++;
          }
          cnt1++;
        }
        //END  saving the restriction object in parallel vectors and matrices


        //BEGIN filling the restriction object with infos coming form the parallel vectors and matrices
        unsigned solutionOffset   = msh->dofmap_get_dof_offset(soltype, curr_proc);
        unsigned solutionOffsetp1 = msh->dofmap_get_dof_offset(soltype, curr_proc + 1);
        for (unsigned lproc = 0; lproc < n_procs; lproc++) {
          masterNode.broadcast(lproc);
          slaveNodes.broadcast(lproc);
          slaveNodesValues.broadcast(lproc);
          for (unsigned i = slaveNodes.begin(); i < slaveNodes.end(); i++) {
            unsigned inode = masterNode[i];
            if (inode >= solutionOffset && inode < solutionOffsetp1 && // inode belongs to curr_proc
                restriction[soltype].find(inode) == restriction[soltype].end()) { // but inode is not set as master node of curr_proc
              counter++;
              for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) { //copy information for lproc to curr_proc
                unsigned jnode = slaveNodes[i][j];
                restriction[soltype][inode][jnode] = slaveNodesValues[i][j];
              }
            }
            else { // either inode does not belong to curr_proc or it was already defined as a master for curr_proc
              if (restriction[soltype].find(inode) != restriction[soltype].end()) { // inode is already defined as master node in curr_proc (either it does or does not belong to curr_proc)
                for (unsigned j = slaveNodes.begin(i); j < slaveNodes.end(i); j++) { // loop on all the columns of restriction[lproc][inode]
                  unsigned jnode = slaveNodes[i][j];
                  double value = slaveNodesValues[i][j];
                  if (inode != jnode || value > 5.) { // if off-diagonal or hanging node for lproc
                    restriction[soltype][inode][jnode] =  value;
                  }
                  if (restriction[soltype].find(jnode) == restriction[soltype].end()) { // if jnode is not yet a master node for curr_proc
                    counter++;
                    for (unsigned k = masterNode.begin(); k < masterNode.end(); k++) {
                      if (masterNode[k] == jnode) { // and if jnode is also a master node for lproc
                        for (unsigned l = slaveNodes.begin(k); l < slaveNodes.end(k); l++) {
                          unsigned lnode = slaveNodes[k][l];
                          restriction[soltype][jnode][lnode] = slaveNodesValues[k][l]; //copy the rule of lptoc into curr_proc
                        }
                        break;
                      }
                    }
                  }
                }
              }
            }
          }
          masterNode.clearBroadcast();
          slaveNodes.clearBroadcast();
          slaveNodesValues.clearBroadcast();
        }
        //END filling the restriction object with infos coming form the parallel vectors and matrices


        pvector->set(curr_proc, counter);
        pvector->close();
        counter = static_cast <unsigned>(floor(pvector->l1_norm() + 0.5));
      }
      delete pvector;

//       for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restriction[soltype].begin(); it1 != restriction[soltype].end(); it1++) {
//         unsigned inode = it1->first;
//         if (restriction[soltype][inode][inode] > 5.) {
//           if (restriction[soltype][inode].size() > 1) {
//             for (std::map<unsigned, std::map<unsigned, double> >::iterator it2 = restriction[soltype].begin(); it2 != restriction[soltype].end(); it2++) {
//               unsigned jnode = it2->first;
//               if (jnode != inode && restriction[soltype][jnode].find(inode) != restriction[soltype][jnode].end()) {
//                 double value =  restriction[soltype][jnode][inode];
//                 for (std::map<unsigned, double> ::iterator it3 = restriction[soltype][inode].begin(); it3 != restriction[soltype][inode].end(); it3++) {
//                   unsigned knode = it3->first;
//                   if (knode != inode) {
//                     restriction[soltype][jnode][knode] = it3->second * value;
//                   }
//                 }
//               }
//             }
//           }
//           restriction[soltype][inode].clear();
//           restriction[soltype][inode][inode] = 0.;
//         }
//       }



      std::vector<std::vector < unsigned > > genealogy;
      std::vector<std::vector < double > > heredity;
      std::vector< unsigned > index;
      std::map < unsigned,  std::map < unsigned, double  > >  restrictionCopy = restriction[soltype];

      for (std::map<unsigned, std::map<unsigned, double> >::iterator it1 = restrictionCopy.begin(); it1 != restrictionCopy.end(); it1++) { // loop all over master, hanging and master+hanging nodes

        genealogy.resize(1);
        heredity.resize(1);
        index.resize(1);

        genealogy[0].resize(1);
        heredity[0].resize(1);

        unsigned inode = it1->first;

        if (restrictionCopy[inode][inode] < 5.) { // only if a real master node

          // initialize master node genealogy and heredity at level 0
          genealogy[0][0] = inode;
          heredity[0][0] = 1.;
          index[0] = 0;

          restriction[soltype][inode].clear();
          restriction[soltype][inode][inode] = 1.;

          unsigned level = 1;
          while ( level > 0 ) {

            // initialize master node genealogy and heredity at genemeric level

            unsigned father = genealogy[level - 1][index[level - 1]];

            genealogy.resize(level + 1);
            genealogy[level].reserve( restrictionCopy[ father ].size() - 1 );
            genealogy[level].resize(0);

            heredity.resize(level + 1);
            heredity[level].reserve( restrictionCopy[ father ].size() - 1 );
            heredity[level].resize(0);

            index.resize(level + 1);
            index[level] = 0;

            unsigned cnt  = 0;
            for (std::map <unsigned, double>::iterator it2 = restrictionCopy[ father ].begin(); it2 != restrictionCopy[ father ].end(); it2++) { // loop on all the father sons
              unsigned son = it2->first;
              bool alreadyFound = false;
              for (unsigned klevel = 0; klevel < level; klevel++) { // check if the son is in the previous genealogy
                for (unsigned k = 0; k < genealogy[klevel].size(); k++) {
                  if ( genealogy[klevel][k] == son ) alreadyFound = true;
                }
              }
              if (!alreadyFound ) { // if never found add the the restionction value in the master node line and zero the hanging node line
                genealogy[level].resize( genealogy[level].size() + 1);
                heredity[level].resize(heredity[level].size() + 1);

                genealogy[level][cnt] = son;
                heredity[level][cnt] = it2->second * heredity[level - 1][index[level - 1]];

                restriction[soltype][inode][son] += heredity[level][cnt];
                cnt++;

                restriction[soltype][son].clear();
                restriction[soltype][son][son] = 0.;
              }
            }

            if ( cnt > 0) {
              level++;
            }
            else {
              bool test = true;
              while ( test && level > 0) {
                index[level - 1]++;
                test = false;
                if ( index[level - 1] == genealogy[level - 1].size() ) {
                  level--;
                  test = true;
                }
              }
            }
          }
        }
	else{
	  restriction[soltype][inode].clear();
	  restriction[soltype][inode][inode] = 0.;
	}
      }

      MyVector <unsigned> InterfaceSolidMarkNode(interfaceSolidMark[soltype].size());
      MyVector <short unsigned> InterfaceSolidMarkValue(interfaceSolidMark[soltype].size());

      unsigned cnt = 0;
      for (std::map<unsigned, bool >::iterator it = interfaceSolidMark[soltype].begin(); it != interfaceSolidMark[soltype].end(); it++) {
        InterfaceSolidMarkNode[cnt] = it->first;
        InterfaceSolidMarkValue[cnt] = it->second;
        cnt++;
      }
      InterfaceSolidMarkNode.stack();
      InterfaceSolidMarkValue.stack();

      for (unsigned lproc = 0; lproc < n_procs; lproc++) {
        InterfaceSolidMarkNode.broadcast(lproc);
        InterfaceSolidMarkValue.broadcast(lproc);
        for (unsigned i = InterfaceSolidMarkNode.begin(); i < InterfaceSolidMarkNode.end(); i++) {
          unsigned jnode = InterfaceSolidMarkNode[i];
          if ( restriction[soltype].find(jnode) != restriction[soltype].end()) {
            interfaceSolidMark[soltype][jnode] = InterfaceSolidMarkValue[i];
          }
        }
        InterfaceSolidMarkNode.clearBroadcast();
        InterfaceSolidMarkValue.clearBroadcast();
      }
    }
  }


  
  
  
  
  
  
} //end namespace femus
