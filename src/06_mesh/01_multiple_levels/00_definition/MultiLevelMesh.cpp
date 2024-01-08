/*=========================================================================

 Program: FEMUS
 Module: MultiLevelMesh
 Authors: Simone Bnà, Eugenio Aulisa, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "MultiLevelMesh.hpp"
#include "Mesh.hpp"
#include "ElemType.hpp"
#include "Elem.hpp"
#include "SparseMatrix.hpp"
#include "NumericVector.hpp"
#include "FemusConfig.hpp"
#include "MeshRefinement.hpp"
#include "Domain.hpp"

#include "Files.hpp"
#include "MED_IO.hpp"

//C++ include
#include <iostream>


namespace femus {


//==== Constructors / Destructor - BEGIN ======== 

MultiLevelMesh::~MultiLevelMesh() {  
    ///@todo these delete's should not be in the Destructor, because the new's are not in the constructor, 
    // so you create problems when you try to use the default copy constructor and copy pointers!

    
    DeleteLevelsZero();
    
    
    DeleteWriter();
    

    DeleteFETypesForExistingGeomElements();
    
    
}

//---------------------------------------------------------------------------------------------------
MultiLevelMesh::MultiLevelMesh() : _gridn0(0)
  {

  InitializeGeomElemFlag();
  
  
  InitializeFETypes();
  
  
  InitializeWriter();

}



//---------------------------------------------------------------------------------------------------
MultiLevelMesh::MultiLevelMesh(const unsigned short &igridn,
                               const unsigned short &igridr,
                               const char mesh_file[], 
                               const char GaussOrder[], 
                               const double Lref,
                               bool (* SetRefinementFlag)(const std::vector < double > &x, const int &ElemGroupNumber,const int &level) )  :
    _gridn0(igridn)
    {
        

    InitializeGeomElemFlag();
  
    InitializeFETypes();
  
    InitializeWriter();
    

    InitializeLevelsZeroAndAllocateCoarse(_gridn0);

    std::cout << " Reading corse mesh from file: " << mesh_file << std::endl;
    
    _level0[0]->ReadCoarseMesh(mesh_file, Lref, _finiteElementGeometryFlag);

    BuildFETypesBasedOnExistingCoarseMeshGeomElements(GaussOrder);

    
    //finer meshes **************
    RefineMesh(igridn, igridr, SetRefinementFlag);
  

}


//==== Constructors / Destructor - END ======== 




//==== Basic - BEGIN  ========


/** Get the dimension of the problem (1D, 2D, 3D) from one Mesh (level 0 always exists, after initialization) */
const unsigned MultiLevelMesh::GetDimension() const {
      return _level0[LEVEL_AT_WHICH_YOU_PICK_THE_DIM]->GetDimension();
    }
    
    
//---------------------------------------------------------------------------------------------

void MultiLevelMesh::PrintInfo() const {
    
    std::cout << " Number of uniform mesh refinement: " << _gridn << std::endl;
    for(int i = 0; i < _gridn; i++) {
      _level[i]->PrintInfo();
    }
    
}



//==== Basic - END  ========



//==== Coarse level - BEGIN  ========


void MultiLevelMesh::ReadCoarseMesh(const char mesh_file[], const char GaussOrder[], const double Lref)
{
    
  const bool read_groups = true; //by default groups are read
  const bool read_boundary_groups = true; //by default boundary groups are read
  
    ReadCoarseMesh(mesh_file, GaussOrder, Lref, read_groups, read_boundary_groups);
   
}

void MultiLevelMesh::ReadCoarseMeshFileReadingBeforePartitioning(const std::string mesh_file, const double Lref, const bool read_groups, const bool read_boundary_groups)
{
    
    _gridn0 = 1;

    InitializeLevelsZeroAndAllocateCoarse(_gridn0);
    
    std::cout << " Reading corse mesh from file: " << mesh_file << std::endl;
    
    _level0[0]->ReadCoarseMeshBeforePartitioning(mesh_file, Lref, _finiteElementGeometryFlag, read_groups, read_boundary_groups);

}

void MultiLevelMesh::ReadCoarseMeshOnlyFileReading(const char mesh_file[], const double Lref, const bool read_groups, const bool read_boundary_groups)
{
    
    _gridn0 = 1;

    InitializeLevelsZeroAndAllocateCoarse(_gridn0);
    
    std::cout << " Reading corse mesh from file: " << mesh_file << std::endl;
    
    _level0[0]->ReadCoarseMesh(mesh_file, Lref, _finiteElementGeometryFlag, read_groups, read_boundary_groups);

}


void MultiLevelMesh::ReadCoarseMesh(const char mesh_file[], const char GaussOrder[], const double Lref, const bool read_groups, const bool read_boundary_groups)
{
    
   ReadCoarseMeshOnlyFileReading(mesh_file, Lref, read_groups, read_boundary_groups);

    BuildFETypesBasedOnExistingCoarseMeshGeomElements(GaussOrder);

    PrepareNewLevelsForRefinement();
    
}


void MultiLevelMesh::ReadCoarseMesh(const std::string mesh_file, const double Lref, const bool read_groups, const bool read_boundary_groups) {
  
  ReadCoarseMeshFileReadingBeforePartitioning(mesh_file, Lref, read_groups, read_boundary_groups);
    
  GetLevelZero(0)->dofmap_build_all_fe_families_and_elem_and_node_structures();
 
  BuildFETypesBasedOnExistingCoarseMeshGeomElements();
  
  PrepareNewLevelsForRefinement();

  }
  

void MultiLevelMesh::ReadCoarseMesh(const std::string mesh_file) {

  const double Lref = 1.; //by default we do not scale the coordinates
  const bool read_groups = true; //by default groups are read
  const bool read_boundary_groups = true; //by default boundary groups are read

  ReadCoarseMesh(mesh_file, Lref, read_groups, read_boundary_groups);

}



void MultiLevelMesh::GenerateCoarseBoxMesh(
        const unsigned int nx, const unsigned int ny, const unsigned int nz,
        const double xmin, const double xmax,
        const double ymin, const double ymax,
        const double zmin, const double zmax,
        const ElemType type, const char GaussOrder[])
{
    
    _gridn0 = 1;

    InitializeLevelsZeroAndAllocateCoarse(_gridn0);
    
    std::cout << " Building brick mesh using the built-in mesh generator" << std::endl;

    _level0[0]->GenerateCoarseBoxMesh(nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,type, _finiteElementGeometryFlag);

    BuildFETypesBasedOnExistingCoarseMeshGeomElements(GaussOrder);

    PrepareNewLevelsForRefinement();
    
}


//==== Coarse level - END  ========


//==== Coarse level, Geom Elem - BEGIN ========

 void MultiLevelMesh::InitializeGeomElemFlag() {
     
    _finiteElementGeometryFlag.resize(N_GEOM_ELS, false);
     
 }
 

//==== Coarse level, Geom Elem - END ======== 



//==== Coarse level, Group Info - BEGIN  ========

   void  MultiLevelMesh::set_group_info(const std::string relative_path_to_input_folder) {
      
  
  MED_IO  med_io(*GetLevel(0)/* coarse or whatever?*/);
  
  const std::string mesh_file_location = Files::get_input_file_with_prefix(get_mesh_filename(), relative_path_to_input_folder)/*relative_path_to_input_folder + Files::_application_input_directory + ml_sol->GetMLMesh()->get_mesh_filename() */;
  ;
  
    hid_t  file_id = med_io.open_mesh_file( mesh_file_location );

    const std::vector< std::string > mesh_menus = med_io.get_mesh_names(file_id);
         
    
        std::vector< std::vector< GroupInfo > >  group_info_all_meshes(mesh_menus.size());  

    for (unsigned m = 0; m < group_info_all_meshes.size(); m++)   {
        
            group_info_all_meshes[m] = med_io.get_all_groups_per_mesh(file_id, mesh_menus[m]);
     }    
 
     _group_info_all_meshes = group_info_all_meshes;
 
  }

//==== Coarse level, Group Info - END  ========



//==== Multilevel - BEGIN ======== 


void MultiLevelMesh::RefineMesh( const unsigned short &igridn, 
                                 const unsigned short &igridr,
                                 bool (* SetRefinementFlag)(const std::vector < double >& x, const int &ElemGroupNumber, const int &level) )
{

    _gridn0 = igridn;

    _level0.resize(_gridn0);

    //totally refined meshes **************
    RefineMeshesTotally(igridr);

    
    //partially refined meshes **************
    RefineMeshesPartially(igridr, SetRefinementFlag);
    

   // copy _level0 into the new levels _level **************
    CopyLevelsZeroIntoNewLevels();

}

  
  
  void MultiLevelMesh::RefineMeshesTotally(const unsigned short &igridr) {
      
    for (unsigned i = 1; i < igridr; i++) {
      MeshRefinement meshcoarser(*_level0[i-1u]);
      meshcoarser.FlagAllElementsToBeRefined();

      _level0[i] = new Mesh();
      MeshRefinement meshfiner(*_level0[i]);
      meshfiner.RefineMesh(i, _level0[i-1], _finiteElement);
    }
      
      
  }
  

   void MultiLevelMesh::RefineMeshesPartially(const unsigned short &igridr,
                                              bool (* SetRefinementFlag)(const std::vector < double > &x, const int &ElemGroupNumber,const int &level) ) {
       
 
      if(SetRefinementFlag == NULL) {

    }
    else{
      Mesh::_SetRefinementFlag = SetRefinementFlag;
      Mesh::set_is_refinement_function_defined(true);
    }
    
    //Here there are two alternatives:
    //When igridn is STRICTLY LARGER THAN igridr, one must provide a refinement flag function, ... 
    //When igridn is SMALLER THAN OR EQUAL igridr, this loop is ignored
    for (unsigned i = igridr; i < _gridn0; i++) {
        
      if(Mesh::get_is_refinement_function_defined() == false) {
        std::cout << "Set Refinement Region flag is not defined! " << std::endl;
        exit(1);
      }
      else {
	MeshRefinement meshcoarser(*_level0[i-1u]);
        meshcoarser.FlagElementsToBeRefined();
	//meshcoarser.FlagOnlyEvenElementsToBeRefined();
      }
      
      _level0[i] = new Mesh();
      MeshRefinement meshfiner(*_level0[i]);
      meshfiner.RefineMesh(i,_level0[i-1],_finiteElement);
    
    }
    
    

   }
   
   
   void MultiLevelMesh::CopyLevelsZeroIntoNewLevels() {
       
   _gridn = _gridn0;
    _level.resize(_gridn);
    for(int i = 0; i < _gridn; i++) {
        _level[i] = _level0[i];
     }
    
   }
   
   
   
    void MultiLevelMesh::InitializeLevelsZeroAndAllocateCoarse(const unsigned short gridn) {
  
    _level0.resize(gridn);

    //coarse mesh
    _level0[0] = new Mesh();
        
    }
    
    
    
   void MultiLevelMesh::DeleteLevelsZero() {
       
    for (unsigned i = 0; i < _level0.size(); i++) {
        if(_level0[i] != NULL)  delete _level0[i];
    }
    
   }

   

  
  
void MultiLevelMesh::PrepareNewLevelsForRefinement() {
    
    _gridn = _gridn0;
    _level.resize(_gridn);
    _level[0] = _level0[0];

}

//---------------------------------------------------------------------------------------------

void MultiLevelMesh::EraseCoarseLevels(unsigned levels_to_be_erased) {
    
    _gridn -= levels_to_be_erased;
    
    for(int i=0; i<_gridn; i++) {
      _level[i]=_level0[i+levels_to_be_erased];
      _level[i]->SetLevel(i);
    }
    
}



void MultiLevelMesh::AddAMRMeshLevel()
{

  //AMR refine mesh
   _level0.resize(_gridn0+1u);

  MeshRefinement meshcoarser(*_level0[_gridn0-1u]);
  meshcoarser.FlagElementsToBeRefined();

  _level0[_gridn0] = new Mesh();
  MeshRefinement meshfiner(*_level0[_gridn0]);
  meshfiner.RefineMesh(_gridn0,_level0[_gridn0-1u],_finiteElement);

  _level.resize(_gridn+1u);
  _level[_gridn]=_level0[_gridn0];

  _gridn0++;
  _gridn++;
}


//==== Multilevel - END ======== 



//==== Multilevel - File output - BEGIN  ========

 void MultiLevelMesh::InitializeWriter() {
     
   _writer = NULL;
    
 }

 
 void MultiLevelMesh::DeleteWriter() {
     
       if(_writer != NULL) delete _writer;
    
 } 

//==== Multilevel - File output - END ========



//==== FE - BEGIN ======== 


 void MultiLevelMesh::InitializeFETypes() {
     
   for(int i = 0; i < N_GEOM_ELS; i++) {
     for(int j = 0; j < NFE_FAMS; j++) {
      _finiteElement[i][j] = NULL;
     }
   }
 
 }
 


 void MultiLevelMesh::DeleteFETypesForExistingGeomElements() {
     
        for(unsigned i = 0; i < N_GEOM_ELS; i++){
        
      if( _finiteElementGeometryFlag[i]) {
      for(unsigned j = 0; j < NFE_FAMS; j++){
        if (_finiteElement[i][j] != NULL)  { delete _finiteElement[i][j]; }
        }
      }
      
    }

 }


 
 
  void MultiLevelMesh::BuildFETypesBasedOnExistingCoarseMeshGeomElements(const char GaussOrder[]) {
      
    if(_finiteElementGeometryFlag[0]) {
      _finiteElement[0][0] = new elem_type_3D("hex",  fe_fams[CONTINUOUS_LINEAR].c_str()      , GaussOrder);
      _finiteElement[0][1] = new elem_type_3D("hex",  fe_fams[CONTINUOUS_SERENDIPITY].c_str() , GaussOrder);
      _finiteElement[0][2] = new elem_type_3D("hex",  fe_fams[CONTINUOUS_BIQUADRATIC].c_str() , GaussOrder);
      _finiteElement[0][3] = new elem_type_3D("hex",  fe_fams[DISCONTINUOUS_CONSTANT].c_str() , GaussOrder);
      _finiteElement[0][4] = new elem_type_3D("hex",  fe_fams[DISCONTINUOUS_LINEAR].c_str()   , GaussOrder);
    }
    if(_finiteElementGeometryFlag[1]) {
      _finiteElement[1][0] = new elem_type_3D("tet",  fe_fams[CONTINUOUS_LINEAR].c_str()     , GaussOrder);
      _finiteElement[1][1] = new elem_type_3D("tet",  fe_fams[CONTINUOUS_SERENDIPITY].c_str(), GaussOrder);
      _finiteElement[1][2] = new elem_type_3D("tet",  fe_fams[CONTINUOUS_BIQUADRATIC].c_str(), GaussOrder);
      _finiteElement[1][3] = new elem_type_3D("tet",  fe_fams[DISCONTINUOUS_CONSTANT].c_str(), GaussOrder);
      _finiteElement[1][4] = new elem_type_3D("tet",  fe_fams[DISCONTINUOUS_LINEAR].c_str()  , GaussOrder);
    }
    if(_finiteElementGeometryFlag[2]) {
      _finiteElement[2][0] = new elem_type_3D("wedge",  fe_fams[CONTINUOUS_LINEAR].c_str()     , GaussOrder);
      _finiteElement[2][1] = new elem_type_3D("wedge",  fe_fams[CONTINUOUS_SERENDIPITY].c_str(), GaussOrder);
      _finiteElement[2][2] = new elem_type_3D("wedge",  fe_fams[CONTINUOUS_BIQUADRATIC].c_str(), GaussOrder);
      _finiteElement[2][3] = new elem_type_3D("wedge",  fe_fams[DISCONTINUOUS_CONSTANT].c_str(), GaussOrder);
      _finiteElement[2][4] = new elem_type_3D("wedge",  fe_fams[DISCONTINUOUS_LINEAR].c_str()  , GaussOrder);
    }
    if(_finiteElementGeometryFlag[3]) {
      _finiteElement[3][0] = new elem_type_2D("quad",  fe_fams[CONTINUOUS_LINEAR].c_str()      , GaussOrder);
      _finiteElement[3][1] = new elem_type_2D("quad",  fe_fams[CONTINUOUS_SERENDIPITY].c_str() , GaussOrder);
      _finiteElement[3][2] = new elem_type_2D("quad",  fe_fams[CONTINUOUS_BIQUADRATIC].c_str() , GaussOrder);
      _finiteElement[3][3] = new elem_type_2D("quad",  fe_fams[DISCONTINUOUS_CONSTANT].c_str() , GaussOrder);
      _finiteElement[3][4] = new elem_type_2D("quad",  fe_fams[DISCONTINUOUS_LINEAR].c_str()   , GaussOrder);
    }
    if(_finiteElementGeometryFlag[4]) {
      _finiteElement[4][0] = new elem_type_2D("tri",   fe_fams[CONTINUOUS_LINEAR].c_str()     , GaussOrder);
      _finiteElement[4][1] = new elem_type_2D("tri",   fe_fams[CONTINUOUS_SERENDIPITY].c_str(), GaussOrder);
      _finiteElement[4][2] = new elem_type_2D("tri",   fe_fams[CONTINUOUS_BIQUADRATIC].c_str(), GaussOrder);
      _finiteElement[4][3] = new elem_type_2D("tri",   fe_fams[DISCONTINUOUS_CONSTANT].c_str(), GaussOrder);
      _finiteElement[4][4] = new elem_type_2D("tri",   fe_fams[DISCONTINUOUS_LINEAR].c_str()  , GaussOrder);
    }
    _finiteElementGeometryFlag[5]=1;
    _finiteElement[5][0] = new elem_type_1D("line",   fe_fams[CONTINUOUS_LINEAR].c_str()     , GaussOrder);
    _finiteElement[5][1] = new elem_type_1D("line",   fe_fams[CONTINUOUS_SERENDIPITY].c_str(), GaussOrder);
    _finiteElement[5][2] = new elem_type_1D("line",   fe_fams[CONTINUOUS_BIQUADRATIC].c_str(), GaussOrder);
    _finiteElement[5][3] = new elem_type_1D("line",   fe_fams[DISCONTINUOUS_CONSTANT].c_str(), GaussOrder);
    _finiteElement[5][4] = new elem_type_1D("line",   fe_fams[DISCONTINUOUS_LINEAR].c_str()  , GaussOrder);
    
    
    SetFiniteElementPtrOnCoarseMesh();
    
  }

  
  
  void MultiLevelMesh::BuildFETypesBasedOnExistingCoarseMeshGeomElements() {
      
    if(_finiteElementGeometryFlag[0]) {
      _finiteElement[0][0] = new elem_type_3D("hex",  fe_fams[CONTINUOUS_LINEAR].c_str()     );
      _finiteElement[0][1] = new elem_type_3D("hex",  fe_fams[CONTINUOUS_SERENDIPITY].c_str());
      _finiteElement[0][2] = new elem_type_3D("hex",  fe_fams[CONTINUOUS_BIQUADRATIC].c_str());
      _finiteElement[0][3] = new elem_type_3D("hex",  fe_fams[DISCONTINUOUS_CONSTANT].c_str());
      _finiteElement[0][4] = new elem_type_3D("hex",  fe_fams[DISCONTINUOUS_LINEAR].c_str()  );
    }
    if(_finiteElementGeometryFlag[1]) {
      _finiteElement[1][0] = new elem_type_3D("tet",  fe_fams[CONTINUOUS_LINEAR].c_str()     );
      _finiteElement[1][1] = new elem_type_3D("tet",  fe_fams[CONTINUOUS_SERENDIPITY].c_str());
      _finiteElement[1][2] = new elem_type_3D("tet",  fe_fams[CONTINUOUS_BIQUADRATIC].c_str());
      _finiteElement[1][3] = new elem_type_3D("tet",  fe_fams[DISCONTINUOUS_CONSTANT].c_str());
      _finiteElement[1][4] = new elem_type_3D("tet",  fe_fams[DISCONTINUOUS_LINEAR].c_str()  );
    }
    if(_finiteElementGeometryFlag[2]) {
      _finiteElement[2][0] = new elem_type_3D("wedge",  fe_fams[CONTINUOUS_LINEAR].c_str()     );
      _finiteElement[2][1] = new elem_type_3D("wedge",  fe_fams[CONTINUOUS_SERENDIPITY].c_str());
      _finiteElement[2][2] = new elem_type_3D("wedge",  fe_fams[CONTINUOUS_BIQUADRATIC].c_str());
      _finiteElement[2][3] = new elem_type_3D("wedge",  fe_fams[DISCONTINUOUS_CONSTANT].c_str());
      _finiteElement[2][4] = new elem_type_3D("wedge",  fe_fams[DISCONTINUOUS_LINEAR].c_str()  );
    }
    if(_finiteElementGeometryFlag[3]) {
      _finiteElement[3][0] = new elem_type_2D("quad",  fe_fams[CONTINUOUS_LINEAR].c_str()     );
      _finiteElement[3][1] = new elem_type_2D("quad",  fe_fams[CONTINUOUS_SERENDIPITY].c_str());
      _finiteElement[3][2] = new elem_type_2D("quad",  fe_fams[CONTINUOUS_BIQUADRATIC].c_str());
      _finiteElement[3][3] = new elem_type_2D("quad",  fe_fams[DISCONTINUOUS_CONSTANT].c_str());
      _finiteElement[3][4] = new elem_type_2D("quad",  fe_fams[DISCONTINUOUS_LINEAR].c_str()  );
    }
    if(_finiteElementGeometryFlag[4]) {
      _finiteElement[4][0] = new elem_type_2D("tri",  fe_fams[CONTINUOUS_LINEAR].c_str()     );
      _finiteElement[4][1] = new elem_type_2D("tri",  fe_fams[CONTINUOUS_SERENDIPITY].c_str());
      _finiteElement[4][2] = new elem_type_2D("tri",  fe_fams[CONTINUOUS_BIQUADRATIC].c_str());
      _finiteElement[4][3] = new elem_type_2D("tri",  fe_fams[DISCONTINUOUS_CONSTANT].c_str());
      _finiteElement[4][4] = new elem_type_2D("tri",  fe_fams[DISCONTINUOUS_LINEAR].c_str()  );
    }
    _finiteElementGeometryFlag[5]=1;                                                         
    _finiteElement[5][0] = new elem_type_1D("line",  fe_fams[CONTINUOUS_LINEAR].c_str()      );
    _finiteElement[5][1] = new elem_type_1D("line",  fe_fams[CONTINUOUS_SERENDIPITY].c_str() );
    _finiteElement[5][2] = new elem_type_1D("line",  fe_fams[CONTINUOUS_BIQUADRATIC].c_str() );
    _finiteElement[5][3] = new elem_type_1D("line",  fe_fams[DISCONTINUOUS_CONSTANT].c_str() );
    _finiteElement[5][4] = new elem_type_1D("line",  fe_fams[DISCONTINUOUS_LINEAR].c_str()   );
    
    
    SetFiniteElementPtrOnCoarseMesh();
    
  }

  
  
  void MultiLevelMesh::InitializeQuadratureWithFEEvalsOnExistingCoarseMeshGeomElements(const char * GaussOrder) {
      
    if(_finiteElementGeometryFlag[0]) {
      _finiteElement[0][0]->initialize_quadrature_with_fe_evals_from_child("hex", GaussOrder);
      _finiteElement[0][1]->initialize_quadrature_with_fe_evals_from_child("hex", GaussOrder);
      _finiteElement[0][2]->initialize_quadrature_with_fe_evals_from_child("hex", GaussOrder);
      _finiteElement[0][3]->initialize_quadrature_with_fe_evals_from_child("hex", GaussOrder);
      _finiteElement[0][4]->initialize_quadrature_with_fe_evals_from_child("hex", GaussOrder);
    }
    if(_finiteElementGeometryFlag[1]) {
      _finiteElement[1][0]->initialize_quadrature_with_fe_evals_from_child("tet", GaussOrder); 
      _finiteElement[1][1]->initialize_quadrature_with_fe_evals_from_child("tet", GaussOrder);
      _finiteElement[1][2]->initialize_quadrature_with_fe_evals_from_child("tet", GaussOrder);
      _finiteElement[1][3]->initialize_quadrature_with_fe_evals_from_child("tet", GaussOrder);
      _finiteElement[1][4]->initialize_quadrature_with_fe_evals_from_child("tet", GaussOrder);
    }
    if(_finiteElementGeometryFlag[2]) {
      _finiteElement[2][0]->initialize_quadrature_with_fe_evals_from_child("wedge", GaussOrder);
      _finiteElement[2][1]->initialize_quadrature_with_fe_evals_from_child("wedge", GaussOrder);
      _finiteElement[2][2]->initialize_quadrature_with_fe_evals_from_child("wedge", GaussOrder);
      _finiteElement[2][3]->initialize_quadrature_with_fe_evals_from_child("wedge", GaussOrder);
      _finiteElement[2][4]->initialize_quadrature_with_fe_evals_from_child("wedge", GaussOrder);
    }
    if(_finiteElementGeometryFlag[3]) {
      _finiteElement[3][0]->initialize_quadrature_with_fe_evals_from_child("quad", GaussOrder);
      _finiteElement[3][1]->initialize_quadrature_with_fe_evals_from_child("quad", GaussOrder);
      _finiteElement[3][2]->initialize_quadrature_with_fe_evals_from_child("quad", GaussOrder);
      _finiteElement[3][3]->initialize_quadrature_with_fe_evals_from_child("quad", GaussOrder);
      _finiteElement[3][4]->initialize_quadrature_with_fe_evals_from_child("quad", GaussOrder);
    }
    if(_finiteElementGeometryFlag[4]) {
      _finiteElement[4][0]->initialize_quadrature_with_fe_evals_from_child("tri", GaussOrder);
      _finiteElement[4][1]->initialize_quadrature_with_fe_evals_from_child("tri", GaussOrder);
      _finiteElement[4][2]->initialize_quadrature_with_fe_evals_from_child("tri", GaussOrder);
      _finiteElement[4][3]->initialize_quadrature_with_fe_evals_from_child("tri", GaussOrder);
      _finiteElement[4][4]->initialize_quadrature_with_fe_evals_from_child("tri", GaussOrder);
    }
    _finiteElementGeometryFlag[5] = 1;
    _finiteElement[5][0]->initialize_quadrature_with_fe_evals_from_child("line", GaussOrder);
    _finiteElement[5][1]->initialize_quadrature_with_fe_evals_from_child("line", GaussOrder);
    _finiteElement[5][2]->initialize_quadrature_with_fe_evals_from_child("line", GaussOrder);
    _finiteElement[5][3]->initialize_quadrature_with_fe_evals_from_child("line", GaussOrder);
    _finiteElement[5][4]->initialize_quadrature_with_fe_evals_from_child("line", GaussOrder);
    
    
  }
  

  void MultiLevelMesh::SetFiniteElementPtrOnCoarseMesh() {
      
    _level0[0]->SetFiniteElementPtr(_finiteElement);
    
  }
  

  
//==== FE - END ======== 

  
  

//==== Domain (optional) - BEGIN ======== 


// ========================================================
  void MultiLevelMesh::SetDomain(Domain* domain_in)  {

    _domain = domain_in;

   return;
  }

 // ========================================================
  Domain* MultiLevelMesh::GetDomain() const {

   return _domain;

  }

//==== Domain (optional) - END ======== 



} //end namespace femus


