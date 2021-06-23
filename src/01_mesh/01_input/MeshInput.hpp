/*=========================================================================

 Program: FEMUS
 Module: MeshInput
 Authors: Simone Bnà
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_mesh_MeshInput_hpp__
#define __femus_mesh_MeshInput_hpp__


// Local includes
#include "Mesh.hpp"

// C++ includes
#include <string>
#include <vector>

namespace femus
{


/**
 * This class defines an abstract interface for \p mesh input.
 * Specific classes derived from this class actually implement
 * reading various mesh formats.
 */

// ------------------------------------------------------------
// MeshInput class definition
template <class MT>
class MeshInput
{
 protected:

  /**
   * Default constructor. Will set the _obj to NULL, effectively
   * rendering this object useless.
   */
  explicit
  MeshInput (bool is_parallel_format = false);

  /**
   * Constructor.  Takes a writeable reference to an object.
   * This is the constructor required to read an object.
   */
  explicit
  MeshInput (MT&, const bool is_parallel_format = false);

 public:

  /**
   * Destructor.
   */
  virtual ~MeshInput ();

  /**
   * This method implements reading a mesh from a specified file.
   */
  virtual void read (const std::string& name, vector < vector < double> > &vt, const double Lref, std::vector<bool> &type_elem_flag, const bool read_groups, const bool read_boundary_groups) = 0;


 protected:

  /**
   * Returns the object as a writeable reference.
   */
  MT& GetMesh ();

  const MT & GetMesh () const;

private:

  /**
   * A pointer to a non-const object object.
   * This allows us to read the object from file.
   */
  MT* _obj;

  /**
   * Flag specifying whether this format is parallel-capable.
   * If this is false (default) I/O is only permitted when the mesh
   * has been serialized.
   */
  const bool _is_parallel_format;
};



// ------------------------------------------------------------
// MeshInput inline members
template <class MT>
inline
MeshInput<MT>::MeshInput (const bool is_parallel_format) :
  _obj (NULL),
  _is_parallel_format(is_parallel_format)
{
}



template <class MT>
inline
MeshInput<MT>::MeshInput (MT& obj, const bool is_parallel_format) :
  _obj (&obj),
  _is_parallel_format(is_parallel_format)
{
//    if (!_is_parallel_format && !this->GetMesh().is_serial())
//      {
//        if (this->GetMesh().processor_id() == 0)
//  	{
//  	  std::cout << "Warning: This MeshOutput subclass only supports meshes which have been serialized!" << std::endl;
//         }
//      }
}



template <class MT>
inline
MeshInput<MT>::~MeshInput ()
{
}



template <class MT>
inline
MT& MeshInput<MT>::GetMesh ()
{
  if (_obj == NULL) exit(1);
  return *_obj;
}


template <class MT>
inline
const MT & MeshInput<MT>::GetMesh() const 
{
  if (_obj == NULL) exit(1);
  return *_obj;
}



} // namespace femus


#endif 
