#ifndef __mggeomel_h__
#define __mggeomel_h__


#include <string>

#include "Typedefs_conf.hpp"
#include "VBTypeEnum.hpp"
#include "FETypeEnum.hpp"


//here, we must do in such a way that we dont need the "false dynamic" initializer
//Prol is the vector that holds the matrix for passing 
//   from LINEAR to QUADRATIC NODES,
//or from LAGRANGE LINEAR to LAGRANGE QUADRATIC DOFS
// The embedding matrices have different meaning, they are used for REFINEMENT 
// but of the SAME ELEMENT TYPE: so the childs are of the same type as the father
// Prol is for CHANGING the ORDER.
//So, we can distiguish between PROLONGATION (ONE ELEMENT, from linear to quadratic )
//                          and REFINEMENT (FROM ONE element of some type TO MORE elements of the SAME TYPE) 
//static is a STORAGE CLASS: the possibilities are static, auto, extern, register
// Every variable and function has one of these storage classes; if a declaration does not specify
// the storage class, a context-dependent default is used (eg., extern for all top-level declarations
// in a source file; auto for variables declared in function bodies).

// Storage class 	Lifetime 	Linkage
// extern 	static 	external (whole program)
// static 	static 	internal (translation unit only)
// auto, register 	function call 	(none)
// static member variables: in C++, member variables declared as static inside class definitions 
// are class variables (shared between all class instances, as opposed to instance variables).


//This geomel contains both the QUADRATIC and LINEAR parts, 
//and both the VOLUME and BOUNDARY Parts, for something, and other for other things...

//We should retain the QL thing in a single class, but templatize it 
//with respect to the SHAPE
//why doesnt libmesh do that?
//Well,libmesh distinguishes the elements by their SHAPE,
//but in VB they are NOT GROUPED.
//so an Hex and Quad are not grouped,
//and also Tet and Tri are not grouped.
//It is more important to group
//     Hex and Tet as Cells
//and Quad and Tri as Faces
//Also, here the "multigrid" part concerns only the VOLUME
//because the operators act on the whole volume DOFS, which
//clearly involve the boundary dofs as well

//Well, actually here we have a hybrid dependency,
//because every geomel has both a QL part and a VB part
//The VB part is automatic
//The QL part is unfortunately EXPLICIT
//Well, also the QL part must be explicitly specified
//The idea is doing templates for QUADRILATERAL and TRIANGLE,
//and doing them as EXPLICITLY INSTANTIATED, so that i have 
//both of them available
//Then, based on the ELTYPE at RUN-TIME, i can decide what to instantiate.
//for now, we should leave
//the problem about GeomEl is that its allocation depends
// on the type
//can i pass nontype parameters and use them for the TEMPLATE?
//what if the nontype parameter is a CLASS instantiation or a struct instantiation?
//then i'll pass it as a NON_TYPE template parameter
//when i pass non-type template parameters, may they be known at run time also?
//no, they must be known at compile time, because the compiler needs to know what 
//code to generate
//so you may do sort of "STATIC" ALLOCATIONS with TEMPLATES instead of DYNAMIC ALLOCATIONS
//by PASSING parameters to the constructor at run-time.
//as an alternative you generate code by EXPLICIT INSTANTIATION,
//in which case you do not need to generate the code of the functions/class types you call,
//it was already generated before.

//The possible sons could be based on 'space_dimension' and ELTYPE. Every of them has the QUADRATIC and LINEAR parts:
//Hex278Quad94
//Hex278Quad94
//We must clarify the connection between VB and QL
//VB is used a lot in all the equations
//QL is in practice almost never used, in fact the point is that
//it is EXPLICIT. It could be used as LINEAR only if the mesh is LINEAR,
//otherwise it is always intended that it is always QUADRATIC
//if we remove that dependency the code will not work most likely.
//We can turn that dependency into a template nontype parameter.
//This one will turn the dependency on QL as templatizedd, while now 
//it is explicit.
//so, if i go back and make that dependency explicit, i will not easily
//switch between linear mesh and quadratic mesh.
//Well, we have to think that this is MORE than a GEOM ELEMENT,
//it is a MIXEDD ELEMENT, because it contains q and l.
//For the mesh we will only give EITHER a LIN or a QUAD or a CONST order mesh.
//This Mixed is more tightly related to the EQUATION
//So this is more QUADRATIC and LINEAR at the same time,
//it is more like a FE element.
//Later I will do a MeshGeomEl<VB> that has only one order;
//this GeomEl<VB> that has 2 orders and is used
// for multigrid purposes.
//We simply have to distinguish where we need 2 ORDERS
// for REFINEMENT Purposes
// for MULTIGRID Purposes
// for QUAD and LINEAR unknowns.

//for multigrid purposes, you remain at the same order and just do CHILDS of the SAME ORDER.
//So, you dont need Quadratic and LINEAR for MULTIGRID,
//you need a SEPARATE Quadratic<VB> and a SEPARATE Linear<VB>.

//Alright so _Prol is only used for the MESH printing of a LINEAR variable.
//So if we do a MeshGeomEl that has the ability to prolongate the linear variables
//into quadratic. That is associated to the MESH.

//Then, the Mesh will be constructed based on that MeshGeomEl.

//On the other hand, for the multigrid part, i will need TWO Geometric Elements.
// One that is ONLY LINEAR, for LINEAR unknowns, LINEARGeomel<VB>
// Another one that is ONLY QUADRATIC, for Quadratic unknowns QuadraticGeomel<VB>
//Actually these are better FEElements.

//When i have gencase, i'll pass BOTH the MeshGeomEl, because the mesh creation is based on that,
//and the LINEAR and QUADRATIC FEElem, because they are used for the dofs


//Each of them will have ONLY the QUADRATIC or the LINEAR part

//The MeshGeomEl will contain also the Mesh Order, read through the utils
//I'll keep the VB association
//so: MeshGeomEl <VB> that contains the two orders
// FE <VB> quadratic
//FE <VB> linear
//No the point is that you should make them TEMPLATED

//So for FE we will do
//Hex27,Hex8,Quad9,Quad4
//Tet10,Tet4,Tri6,Tri3

//For the GeomEl we will do
// LagrHex27,LagrQuad9
//LagrTet10, LagrTri6
//linear mesh elements are not supported yet

//notice that here the embedding matrices are used as properties of the FE,
//not of the GEOMETRIC elements.

//now, we'll have to do things separately for every Lagrange 
// Then, we'll end up with a linear lagrange and a quadratic lagrange,
//and we'll have to pass these two to the MULTIGRID Part of our GENCASE
//now, in order to have a GENERIC LINEAR LAGRANGE 
//                  and a GENERIC QUADRATIC LAGRANGE,
//either we do two Fathers and convert/use them as fathers
//or try with some pre-compiled templating
//templating is quite diffficult in this case, as 
//you would need a pre-filled struct to be used as a non-type template parameter.
// First, can we pass structs as non-type template parameters?
// Second, how does the compiler act to understand if something is known at compile time?

// I think i'll stick to the Father thing for now (unless i see templates for static data members...)

//====== STATIC vs TEMPLATES
// With TEMPLATES you usually have to put everything in the header
// With STATIC members in C++ instead you have to put the instance in the SOURCE file.!
//Basically, static is something that belongs MORE to the CLASS TYPE than to a specific Class Instantiation.

//So suppose that you have a template for which you need the implementation in the header file
//because it may be implicitly instantiated. Then what if you have static members,
//whose implementations are constrained to stay in a SOURCE file?


//Ok, now i have split all the FE
//For the GeomEl i want to do it ONLY QUADRATIC for now
//So, i only want the geometric element to be EITHER QUADRATIC  or LINEAR, alright.


class GeomEl  {

public:

     GeomEl(const uint dim_in,const uint geomel_type);
    ~GeomEl();
    
    const uint _geomel_type;
    const uint _dim;
    uint _elnds[VB][QL_NODES]; // I SWITCHED THIS ONE: VB first, QL next!!!   //number of nodes of one element [VB][QL]
    std::string name[VB];             ///< element name (volume+surface)
    std::string pname[VB];            ///< print element name (volume+surface) for XDMF print (linear)

 //===== Multigrid   
    uint n_se[VB];                    ///< number of linear subelements (volume+surface)
    
    //from linear to quadratic
//     static const double _Prol[NNDS*NNDSL];  //actually NNDSL should be NDOF_P: one part is geometric, one part is mathematic
    //======embedding matrices
//     static const float _embedding_matrix_q[NCHILDS][NNDS][NNDS];
//     static const float _embedding_matrix_l[NCHILDS][NNDSL][NNDSL];

};

#endif

//every FE must have its own prol like it has its own embedding matrix. Then the print routine asks to each variable how to be printed
//(whether to be prolongated or not)
//if in the mesh object i instantiate in the constructor the correct thing, 
// with which

//concerning the Prol, every FE must have a Prol, in the sense that it must have 
//a way to be printed on every GEOM ELEM type.
//So,for now we'll give the prol only to the linear elements,
//in order to allow them to have the possibility of being projected on the GEOM el
//e should have a prol for FEELemQuad4 to GeomElemQuad8
// and another prol for  FEELemQuad4 to GeomElemQuad9
//the question is still: do we have to associate it to the MATH part
//or to the GEOM part
//We give it to the math part and we say: for a given geom el part, what do you give us, my dear math part?