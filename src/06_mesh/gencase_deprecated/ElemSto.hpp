#ifndef __femus_meshGencase_ElemSto_hpp__
#define __femus_meshGencase_ElemSto_hpp__

// was like this ================================================
//    ELEMENTS (in elem_sto and bd_elem_sto)
// ---------------------------------------------
//    elem_sto[Id(0), nodes (NDOF_FEM),
//           lev(NDOF_FEM+1),
//            pr(NDOF_FEM+2),
//       parent (NDOF_FEM+3),
//           Nchildren(NDOF_FEM+4),
//          CHILDREN]
// --------------------------------------------------
//    bd_elem_sto[Id(0),
//               elem_id (1),
//             lev(2),
//            side(3),
//           n_nodes (4),
//           side-nodes]
// ================================================

#include "Typedefs.hpp"


namespace femus {



//===========================
class ElemStoBase  {

public:

    ElemStoBase(int nnds_in, uint spacedim_in);
    ~ElemStoBase();

    int    _id;
    int    _lev;
    int    _subd;
    int    _nnds;
    int*   _elnds;
    uint _spacedim;


};


//============================
class ElemStoVol : public ElemStoBase {

public:

    ElemStoVol(int nnds_in, uint spacedim);
    ~ElemStoVol();


    int _par;
    int _nch;
    int* _elchs;  //number of children per element, any geom


};



//========================
class ElemStoBdry : public ElemStoBase {

public:

    ElemStoBdry(int nnds_in, uint spacedim);
    ~ElemStoBdry();

    int _vol_id;  // id of the volume element it belongs to
    int _nside;   //number of the side (0-1-2-3..)


};


// ================================================
//    NODES (in nd_sto and nod_val)
// ---------------------------------------------
//    nd_sto[Index(0), pr(1),lev (2),levp(3),var(4)]
// --------------------------------------------------
//    _nod_coords[x,y,z]
// ================================================


//========================
class NodeSto {

public:

    NodeSto(int nl1_in,int n_levs_in);
    ~NodeSto();

    int _id;
    int _lev;
    int _subd;
    int _levp;  //for linear nodes
    int _var;


};



//ok, now you can do a vectors of objects, each of which is an object, not a pointer to an object
//STATIC or DINAMIC allocation is sthg that doesnt depend only on the exterior. An object of a
//  class might be instantiated "statically" but some of its elements are instantiated "dynamically"...
//then what happens?where is the object, partly in the heap and partly in the stack?
//new ElemSto: the class goes all in the heap
//ElemSto with new inside: the class goes partly in the stack and partly in the heap?!?

//templates
//ok, it is the compiler that decides whether he can determine the template parameters
//its not a question about "const" or "not const"; he sees if he can add things.

//if it is a const int, you can do it; if it is a reference to sthg else, you CANT do.
//if it is jut an int, you cant do. In fact, the compiler does not execute the program
//so he cant tell you whether the variable can be modified
//but, in the case of const int, the compiler CAN RETRIEVE the NUMBER
// a "chain of const int" ALSO WORKS !!! See this:
//  const int gino=8;
// const int pippo=gino;
//
// ElemSto<pippo> gigio;

//  So the point is: WHEN CAN THE COMPILER RETRIEVE values WITHOUT EXECUTING?
// Whenever he can, you can instantiate template classes.
//I guess with static stuff he can


//with explicit instantiation the code is already generated for a class template.

//can we make a TEMPLATE struct wrt n_dofs
//until you dont include a header file in any source, it is NOT COMPILED clearly!
//the problem is always how to know things at run time
// the solution is template with a static build function

// interesting: when a template is NOT instantiated, its content is COMPLETELY IGNORED
//by the compiler! So, even if you write errors in it, the compiler DOESNT SEE THEM!!!

//how do i do a static array of ElemSto objects that i construct with a parameter?
//i guess you cant
// doing an array of Classes would mean CALLING a CONSTRUCTOR for each of its elements.
//it looks like this is not allowed
//Well, IT DEPENDS. If you define a DEFAULT CONSTRUCTOR (without arguments) you may do... TODO study later

// you can a static array of 100 POINTERS to elem sto. in this way you DONT CALL ANY CONSTRUCTOR!
//then LATER you call the constructors with "new"
//BOTH for double and for Classes, variable length arrays are forbidden by ISO C++!!!




} //end namespace femus



#endif
