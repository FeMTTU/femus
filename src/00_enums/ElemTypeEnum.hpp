#ifndef __femus_enums_ElemTypeEnum_hpp__
#define __femus_enums_ElemTypeEnum_hpp__

/** Defines an \p enum for geometric element types. */

enum ElemType {
    EDGE2=0,    // 0
    EDGE3,      // 1
    EDGE4,      // 2

    TRI3,       // 3
    TRI6,       // 4

    QUAD4,      // 5
    QUAD8,      // 6
    QUAD9,      // 7

    TET4,       // 8
    TET10,      // 9

    HEX8,       // 10
    HEX20,      // 11
    HEX27,      // 12

    PRISM6,     // 13
    PRISM15,    // 14
    PRISM18,    // 15

    PYRAMID5,   // 16

    INFEDGE2,   // 17

    INFQUAD4,   // 18
    INFQUAD6,   // 19

    INFHEX8,    // 20
    INFHEX16,   // 21
    INFHEX18,   // 22

    INFPRISM6,  // 23
    INFPRISM12, // 24

    NODEELEM,   // 25

    REMOTEELEM,   // 26

    INVALID_ELEM
};  // 27 - should always be last

#endif
