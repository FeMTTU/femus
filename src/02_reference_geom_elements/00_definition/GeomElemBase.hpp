#ifndef __femus_meshGencase_GeomElemBase_hpp__
#define __femus_meshGencase_GeomElemBase_hpp__


#include <string>
#include <vector>
#include <iostream>



namespace femus {



class GeomElemBase  {


// ===  Constructors / Destructor - BEGIN =================
public:
    ///runtime selection of Geom Elem
    static  GeomElemBase* build(const std::string geomel_id_in, const uint fe_family);
    
// ===  Constructors / Destructor - END =================


// Geom, Independent of Number of Nodes - BEGIN ===
    virtual unsigned int get_dimension() const = 0;
    virtual unsigned int n_nodes_linear() const = 0;
// Geom, Independent of Number of Nodes - END ===
    
    
// Geom, Dependent of Number of Nodes - BEGIN ===
    virtual unsigned int n_nodes()       const  = 0;
    
    virtual std::vector<unsigned> get_nodes_of_face(const unsigned f) const  = 0;
// Geom, Dependent of Number of Nodes - END ===
    

    
// Refinement - BEGIN ===
public:
    
    virtual float get_embedding_matrix(const uint,const uint,const uint) = 0; //[/*NCHILDS*/][/*NNDS*/][/*NNDS*/]
    virtual double get_prol(const uint) = 0;
    
// Refinement - END ===



// File names - BEGIN ===
public:
    
    virtual std::string  get_name_med()  const  = 0;
    virtual std::string  get_name_xdmf() const  = 0;
// File names - END===
    
    
   
};


} //end namespace femus



// ******************** GEOMETRIC ELEMENTS, NODE NUMBERING - BEGIN **************************

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
//      |   6   \
//      |         \
//      0-----3----1


//
//	0-----2-----1
//
*/
// ******************** GEOMETRIC ELEMENTS, NODE NUMBERING - END **************************



#endif
