/*=========================================================================

 Program: FEMUS
 Module: CurrentElem
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_CurrentElem_hpp__
#define __femus_equations_CurrentElem_hpp__


#include "FETypeEnum.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "DenseVector.hpp"
#include "DenseMatrix.hpp"
#include "ElemType.hpp"


namespace femus {



class SystemTwo;
class CurrentQuantity;



template < typename real_num_mov >
class CurrentElem {
    

  public:
      
    CurrentElem(const unsigned dim, const Mesh * msh_in);

    CurrentElem(const uint iel_in, const uint iproc_in, const uint level, const uint vb, const SystemTwo*,  const MultiLevelMeshTwo& mesh, const std::vector< std::vector<const elem_type*> >  & elem_type, const Mesh * mesh_new);

    inline const uint  GetVb() const {
      return _mesh.get_dim() - _dim;
    }
    
    inline const uint  GetDim() const {
      return _dim;
    }
    

    
    inline const uint*  GetConn() const {
      return &_el_conn[0];
    }
    
    inline uint  GetVolIel() const {
      return _vol_iel_DofObj;
    }
    
    inline const double*  GetNodeCoords() const {
      return &_xx_nds[0];
    }
    
    inline std::vector<uint>  GetDofIndices() const {
      return _el_dof_indices;
    }
    
    //TODO here i have to return non-const because of a function that i will change...
    inline uint*  GetBCDofFlag() {
      return &_bc_eldofs[0];
    }
    
    /** Returns a reference to the element rhs */
    inline DenseVector & Rhs() {
      return _FeM;
    }
    
    /** Returns a reference to the element matrix */
    inline DenseMatrix & Mat() {
      return _KeM;
    }
    
    void  SetDofobjConnCoords();
    
    
    void  PrintOrientation() const;
    
    void  ConvertElemCoordsToMappingOrd(CurrentQuantity& myvect) const;
    
    /** needs the EQUATION basically */
    void  SetElDofsBc();
 
    /** */
    inline const elem_type* GetElemType(const uint fe) const { return  _elem_type[fe]; }
    
    inline const std::vector<const elem_type*> &  GetElemTypeVectorFE() const { return _elem_type; }
     
    /**  */  
    int Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(const uint ndof_in, const uint space_dim, const CurrentQuantity& press_in) const;
   
//     void TransformElemNodesToRef(Domain* mydom, double* refbox_xyz);
    
    const uint GetLevel() const {return _Level;}
    
   //TODO make these private
  //========== Equation-related ========================               
  const SystemTwo * _eqn;  //con questo puoi accedere a dati e funzioni DEL PADRE, NON al FIGLIO
  const MultiLevelMeshTwo & _mesh;


  
  // ========= NEW ===============================================================================
   inline const std::vector< real_num_mov > & get_elem_center_bdry() const {    return _elem_center_bdry_3d;  }
   inline const std::vector< real_num_mov > & get_elem_center()      const {    return _elem_center_3d;  }
    
   void  set_elem_center(const unsigned int iel, const unsigned int xType);
   void  set_elem_center_bdry_3d();
    
   inline short unsigned int geom_type() const { return geom_elem_type; }
  
   void set_coords_at_dofs_and_geom_type(const unsigned int dim,  const unsigned int xType);
   
   vector < vector < real_num_mov > > & get_coords_at_dofs() {  return _coords_at_dofs; }
   vector < vector < real_num_mov > > & get_coords_at_dofs_3d() {  return _coords_at_dofs_3d; }
   const vector < vector < real_num_mov > > & get_coords_at_dofs() const {  return _coords_at_dofs; }
   const vector < vector < real_num_mov > > & get_coords_at_dofs_3d() const {  return _coords_at_dofs_3d; }
   const vector < vector < real_num_mov > > & get_coords_at_dofs_bdry_3d() const {  return _coords_at_dofs_bdry_3d; }
   
   const real_num_mov & get_coords_at_dofs(const unsigned int idim,  const unsigned int idof) const {  return _coords_at_dofs[idim][idof]; }
   const real_num_mov & get_coords_at_dofs_3d(const unsigned int idim,  const unsigned int idof) const {  return _coords_at_dofs_3d[idim][idof]; }
 
  
  void   set_coords_at_dofs_bdry_3d(const unsigned int iel, const unsigned int jface, const unsigned int solType_coords);

  private:

  std::vector < std::vector < real_num_mov > >  _coords_at_dofs;         // must be adept if the domain is moving, otherwise double
  std::vector < std::vector < real_num_mov > >  _coords_at_dofs_3d;
  std::vector< real_num_mov > _elem_center_3d;                           // element center point
  
  std::vector < std::vector < real_num_mov > >  _coords_at_dofs_bdry_3d;
  std::vector< real_num_mov > _elem_center_bdry_3d;

  const uint _dim;         //spatial dimension of the current element (can be different from the mesh dimension!)
  /*const */unsigned _max_size_elem_dofs;                   // conservative: based on line3, quad9, hex27
  const Mesh * _mesh_new;
  short unsigned int geom_elem_type;
  
  static constexpr unsigned int _space_dim = 3;
  
// private functions
  void   set_geom_type(const unsigned int iel) {    geom_elem_type = _mesh_new->GetElementType(iel); }

  void   set_coords_at_dofs(const unsigned int iel, 
                            const unsigned int xType,
                            std::vector < std::vector < real_num_mov > > & coords_at_dofs_in, 
                            const unsigned int space_dim);

// === OLD =====================================================================================
  const std::vector<const elem_type*>  &  _elem_type;
    
// ========================================================================================
//========== Current "EQUATION" Element (ql are TOGETHER ): needs the EQUATION ========================               
  uint                   _el_n_dofs;
  DenseMatrix                  _KeM; 
  DenseVector                  _FeM;
  std::vector<uint> _el_dof_indices;                  // this must become a vect of vect
  std::vector<uint>      _bc_eldofs;                  // this must become a vect of vect
  
// ========================================================================================
//==========  Current Geometric Element:  needs the MESH  ========================
   std::vector<uint>   _el_conn;             /// vector of the global nodes for that element         [NNDS];
   uint    _vol_iel_DofObj;     /// i need to put the element also.
   std::vector<uint>   _el_conn_new;
   std::vector<double> _xx_nds;              /// vector of the node coordinates for that element     [_spacedimension*NNDS];  // this must become a vect of vect
   const uint _mesh_vb;     //index for the mesh

   const uint _Level;  //the level to which the element belongs
      
   const uint _iel;  //the index of the element (input for the parallel partition)
   const uint _proc;
   
  
  };



template < typename real_num_mov >
    CurrentElem<real_num_mov>::CurrentElem(const unsigned dim, const Mesh * msh_in) :
    _max_size_elem_dofs(static_cast< unsigned >(ceil(pow(3, dim)))),
    _dim(dim),
    _mesh_new(msh_in),
//    rest to be thrown away 
    _eqn(NULL),
    _mesh(MultiLevelMeshTwo()),
    _elem_type(std::vector<const elem_type*>()),
    _mesh_vb(1u),
    _Level(1u),
    _iel(1),
    _proc(1)  {
        
    constexpr unsigned int space_dim = 3;    
        
  _coords_at_dofs.resize(dim);
  for (unsigned i = 0; i < dim; i++)  _coords_at_dofs[i].reserve(_max_size_elem_dofs);
        
  _coords_at_dofs_3d.resize(space_dim);
  for (unsigned i = 0; i < space_dim; i++)  _coords_at_dofs_3d[i].reserve(_max_size_elem_dofs);
  
  _elem_center_3d.resize(space_dim);

  //bdry
  _coords_at_dofs_bdry_3d.resize(space_dim);
  for (unsigned i = 0; i < space_dim; i++)  _coords_at_dofs_bdry_3d[i].reserve(_max_size_elem_dofs);
  
    _elem_center_bdry_3d.resize(space_dim);

  
    }
    
    
template < typename real_num_mov >
 void    CurrentElem<real_num_mov>::set_coords_at_dofs_and_geom_type(const unsigned int iel, const unsigned int xType) {
     
     
    set_coords_at_dofs(iel, xType, _coords_at_dofs, _dim);
    
    set_coords_at_dofs(iel, xType, _coords_at_dofs_3d, 3);
   
    set_geom_type(iel);
    
    
  }
    

template < typename real_num_mov >
 void    CurrentElem<real_num_mov>::set_coords_at_dofs(const unsigned int iel, 
                                                                     const unsigned int xType,
                                                                     std::vector < std::vector < real_num_mov > > & coords_at_dofs_in, 
                                                                     const unsigned int space_dim) {
         
  unsigned nDofx = _mesh_new->GetElementDofNumber(iel, xType);
      
   for (int idim = 0; idim < space_dim; idim++)    coords_at_dofs_in[idim].resize(nDofx);
   
    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = _mesh_new->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < space_dim; jdim++) {
        coords_at_dofs_in[jdim][i] = (*_mesh_new->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
      }
    }
    
    
  }
  
  
   
  template < typename real_num_mov >
  void CurrentElem<real_num_mov>::set_coords_at_dofs_bdry_3d(const unsigned int iel, const unsigned int jface, const unsigned int solType_coords) {

      constexpr unsigned int space_dim = 3;
      
    		unsigned nve_bdry = _mesh_new->GetElementFaceDofNumber(iel, jface, solType_coords);
            
	        for (unsigned idim = 0; idim < space_dim; idim++) {  _coords_at_dofs_bdry_3d[idim].resize(nve_bdry); }
		const unsigned felt_bdry = _mesh_new->GetElementFaceType(iel, jface);    
		for(unsigned i_bdry=0; i_bdry < nve_bdry; i_bdry++) {
		  unsigned int i_vol = _mesh_new->GetLocalFaceVertexIndex(iel, jface, i_bdry);
                  unsigned iDof = _mesh_new->GetSolutionDof(i_vol, iel, solType_coords);
		  for(unsigned idim=0; idim < space_dim; idim++) {
		      _coords_at_dofs_bdry_3d[idim][i_bdry] = (*_mesh_new->_topology->_Sol[idim])(iDof);
		  }
		}  
      
      
  }   
  
  
// ========================================================
///Compute the element center

template < typename real_num_mov >
  void CurrentElem<real_num_mov>::set_elem_center(const unsigned int iel, const unsigned int solType_coords) {

  unsigned nDofx = _mesh_new->GetElementDofNumber(iel, solType_coords);

  std::fill(_elem_center_3d.begin(), _elem_center_3d.end(), 0.);
    
  for (unsigned j = 0; j < 3; j++) {  
      for (unsigned i = 0; i < nDofx; i++) {
         _elem_center_3d[j] += _coords_at_dofs_3d[j][i];
       }
    }
    
   for (unsigned j = 0; j < 3; j++) { _elem_center_3d[j] = _elem_center_3d[j]/nDofx; } 
      
  }    

  
// since the Mesh is  Bi-TriQuadratic order, here I need the formula for the center of an Edge, a Triangle, a Quadrilateral
// On the other hand, instead of COMPUTING this quantity, I'd rather RETRIEVE IT from the LAST NODE!
template < typename real_num_mov >
 void CurrentElem<real_num_mov>::set_elem_center_bdry_3d() {
     
            
           std::fill(_elem_center_bdry_3d.begin(), _elem_center_bdry_3d.end(), 0.);
            
            for (unsigned j = 0; j < _space_dim; j++) {
                for (unsigned i = 0; i <       _coords_at_dofs_bdry_3d[j].size(); i++) {
                    _elem_center_bdry_3d[j] += _coords_at_dofs_bdry_3d[j][i];
                }
            }
            
            for (unsigned j = 0; j < _space_dim; j++) {
                _elem_center_bdry_3d[j] = _elem_center_bdry_3d[j]/_coords_at_dofs_bdry_3d[j].size();
            }
                        
 }

} //end namespace femus

// Ok, in this case if we want to reach all the other classes we just 
// use the eq pointer, through which we have the eqn map, through which 
// we have everything...
//Through an equation we get to the equationsmap which gives us access to everything we need
//Ok, I need to pass the Equation, and also I will pass the ref to EqMap, 
// which is protected in the SystemTwo class and I want to leave it protected

//The geometry does not depend on the specific equation
//The rest yes

// CurrentElement and CurrentGaussPoint are supposed to be "External containers"
// to hold the informations needed for those objects

//Since we structured the GenMatRhs in terms of VB, then the current elem and gauss point
//are going to have all the VB structure... maybe we could do an object for only VV
// and one for only BB, maybe templatized. So far we go ahead like this.

//here there should be two pointers:
// to the ABSTRACT GEometrical Element
// and to the ABSTRACT FE Element

//Basically earlier we were using the Mesh and the Equation to handle what was needed by the current element.
//Now we are using the current element and current gauss to call from the big classess,
// Mesh and Equation, what is needed for the assemblying

 //alright, here I need a Current Element, to perform the loop
 //the current element only needs the eqn map, so we can use it everywhere
 //TODO here we only need the GEOMETRIC PART of the CurrentElem, so maybe we will split 
 //between the CurrGeomEl and the CurrFEEl
//in fact this is used for many element loop, but just to retrieve the geometrical properties
//like coords, middle point, etc.

//=======================
//ok here we would need a "REFERENCE REAL ELEMENT"
//the current element contains the absolute coordinates, 
//but not in the reference frame
//ok now i want to set the element center but NOT BASED ON THE MESH

//So, so far we have the CurrentElem class is 
//split into a CURR GEOMETRIC and a CURR EQUATION part.
//The curr geometric is basically filled with the MESH class
//The curr fe is basically filled with the EQUATION class

//====================================================
// The differences with the other applications are that:
// the element matrix is split into blocks, while here it is only one
// Therefore, also the el_dof_indices and the bc_eldof flags are split as vector of vectors

#endif
