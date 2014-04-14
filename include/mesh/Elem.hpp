/*=========================================================================

 Program: FEMUS
 Module: Elem
 Authors: Eugenio Aulisa
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __elem_hpp__
#define __elem_hpp__


namespace femus {


//vertexes,edges,faces,interior,element,element+derivatives
const unsigned NVE[6][6]= {{8,20,26,27,1,4},  //hex
  {4,10,10,10,1,4},   //tet
  {6,15,18,18,1,4},   //wedge
  {4, 8, 8, 9,1,3},   //quad
  {3, 6, 6, 6,1,3},   //tri
  {2, 3, 3, 3,1,2}
};  //line

/**
 * The elem class
*/

class elem {

private:
  int **kel;
  unsigned *elr;
  unsigned *_child_elem_memory;
  unsigned **_child_elem;
  bool _child_elem_flag;

  unsigned **kvtel;
  unsigned *kvtel_memory;
  unsigned *nve;

  unsigned *kvert_memory;
  int *kel_memory;

  unsigned nvt,nv0,nv1,nv2;
  unsigned nel,nelt[6];
  unsigned nelr,nelrt[6];
  unsigned ngroup;

  short unsigned *elt,*elg,*elmat;
  bool *_node_region;
  bool  _node_region_flag;
  unsigned *elf;
  unsigned nelf;
  unsigned **kvert;
public:
// ******************************
  elem(const unsigned & other_nel);
  elem(const elem *elc, const unsigned refindex);
  ~elem() {
    delete [] kvert_memory;
    delete [] kvert;
    delete [] kel_memory;
    delete [] kel;
    delete [] elt;
    delete [] elf;
    delete [] elg;
    delete [] elmat;
    delete [] elr;
    delete [] kvtel_memory;
    delete [] kvtel;
    delete [] nve;
    if(_node_region_flag) delete [] _node_region;
    if(_child_elem_flag){
      delete [] _child_elem_memory;
      delete [] _child_elem;
    }
  }
  // ******************************
  unsigned GetDof(const unsigned iel,const unsigned &inode,const unsigned &type)const;
  unsigned GetDofCoarse(const unsigned iel,const unsigned &inode,const unsigned &type)const;
  unsigned GetElementDofNumber(const unsigned &iel,const unsigned &type=3) const;
  /**
  * Return the local->global node number
  **/
  unsigned GetElementVertexIndex(const unsigned &iel,const unsigned &inode)const {
    return kvert[iel][inode];
  };
  const unsigned* GetElementVertexAddress(const unsigned &iel,const unsigned &inode)const;
  void SetElementVertexIndex(const unsigned &iel,const unsigned &inode, const unsigned &value);
  // ******************************
  unsigned GetFaceVertexIndex(const unsigned &iel,const unsigned &iface, const unsigned &inode) const;
  // ******************************
  short unsigned GetElementType(const unsigned &iel) const;
  void SetElementType(const unsigned &iel, const short unsigned &value);
  // ******************************
  short unsigned GetElementGroup(const unsigned &iel) const;
  void SetElementGroup(const unsigned &iel, const short unsigned &value);
  void SetElementMaterial(const unsigned &iel, const short unsigned &value);
  short unsigned GetElementMaterial(const unsigned &iel) const;
  unsigned GetElementGroupNumber() const;
  void SetElementGroupNumber(const unsigned &value);
  // ******************************
  int GetFaceElementIndex(const unsigned &iel,const unsigned &iface) const;
  void SetFaceElementIndex(const unsigned &iel,const unsigned &iface, const int &value);
  // ******************************
  unsigned GetIndex(const char name[]) const;
  // ******************************
  unsigned GetElementNumber(const char* name="All") const;
  void AddToElementNumber(const unsigned &value, const char name[]);
  void AddToElementNumber(const unsigned &value, short unsigned ielt);
  // ******************************
  unsigned GetRefinedElementNumber(const char name[]="All") const;
  unsigned GetRefinedElementNumber(short unsigned ielt) const;
  void AddToRefinedElementNumber(const unsigned &value, const char name[]="All");
  void AddToRefinedElementNumber(const unsigned &value, short unsigned ielt);
  void InitRefinedToZero();
  unsigned GetRefinedElementIndex(const unsigned &iel) const;
  void SetRefinedElementIndex(const unsigned &iel, const unsigned &value);
  // ******************************
  unsigned GetNodeNumber()const;
  unsigned GetVertexNodeNumber()const;
  unsigned GetMidpointNodeNumber()const;
  unsigned GetCentralNodeNumber()const;
  void SetNodeNumber(const unsigned &value);
  void SetVertexNodeNumber(const unsigned &value);
  void SetMidpointNodeNumber(const unsigned &value);
  void SetCentralNodeNumber(const unsigned &value);
  // ******************************
  unsigned GetElementFaceNumber(const unsigned &iel,const unsigned &type=1)const;
  unsigned GetElementSquareFaceNumber(const unsigned &iel)const;
  unsigned GetElementTriangleFaceNumber(const unsigned &iel)const;
  // ******************************
  void AllocateVertexElementMemory();
  unsigned GetVertexElementNumber(const unsigned &inode)const;
  unsigned GetVertexElementIndex(const unsigned &inode,const unsigned &jnode)const;
  const unsigned *GetVertexElementAddress(const unsigned &inode,const unsigned &jnode)const;
  void SetVertexElementIndex(const unsigned &inode,const unsigned &jnode, const unsigned &value);
  // ******************************
  unsigned GetElementFather(const unsigned &iel) const;
  void SetElementFather(const unsigned &iel, const unsigned &value);
  void SetNumberElementFather(const unsigned &value);
  // **********************************
  bool GetNodeRegion(const unsigned &jnode) const;
  void SetNodeRegion(const unsigned &jnode, const bool &value);
  void AllocateNodeRegion();
  void AllocateChildrenElement(const unsigned &ref_index);
  void SetChildElement(const unsigned &iel,const unsigned &json, const unsigned &value);
  unsigned GetChildElement(const unsigned &iel,const unsigned &json) const;
  
};


} //end namespace femus


#endif


// ********************  class solver**************************

/*

//         7------14-------6
//        /|              /|
//       / |             / |
//     15  |   25      13  |
//     /  19      22   /  18
//    /    |          /    |
//   4------12-------5     |  
//   | 23  |   26    | 21  |
//   |     3------10-|-----2
//   |    /          |    / 
//  16   /  20      17   /
//   | 11      24    |  9
//   | /             | /
//   |/              |/
//   0-------8-------1


//           5
//          /|\
//         / | \
//        /  |  \
//      11   |  10
//      /   14    \
//     /     |     \
//    /      |      \
//   3--------9------4
//   |  17   |  16   |   
//   |       2       |
//   |      / \      |
//   |     /   \     |
//  12    / 15  \   13 
//   |   8       7   |
//   |  /         \  | 
//   | /           \ | 
//   |/             \|   
//   0-------6-------1

 
//            3 
//           /|\
//          / | \
//         /  |  \
//        9   |   8
//       /    |    \  
//      /     |     \
//     /      7      \ 
//    2-------|5------1
//     \      |      /    
//      \     |     /  
//       \    |    /
//        6   |   4  
//         \  |  /    
//          \ | /    
//           \|/   
//            0                


//      3-----6-----2
//      |           |
//      |           |
//      7     8     5
//      |           |
//      |           |
//      0-----4-----1


//      2 
//      | \
//      |   \
//      5     4
//      |       \
//      |         \
//      0-----3----1


//
//	0-----2-----1
//
*/