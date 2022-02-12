/*=========================================================================

Program: FEMUS
Module: HtmlInputParser
Authors: Simone Bnà

Copyright (c) FEMUS
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_utils_HtmlInputParser_hpp__
#define __femus_utils_HtmlInputParser_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "InputParser.hpp"
#include "iostream"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------

namespace femus {

/**
* This class is the concrete class for Html format input parsing. It's not yet implemented.
*/

class HtmlInputParser : public InputParser {

public:

    /** Constructor */
    HtmlInputParser(const std::string& infile);

    /** Destructor. */
    ~HtmlInputParser ();

    /** See father class */
    const bool isTrue(const std::string& base, const std::string& type) { return true; }; // added only for watning pourposes in compile

    /** See father class */
    const unsigned int getSize(const std::string& base) { return 0u; };// added only for watning pourposes in compile

    /** See father class */
    const int getValue(const std::string& fullpath, const int defaultvalue) { return 0; };// added only for watning pourposes in compile

    /** See father class */
    const double getValue(const std::string& fullpath, const double defaultvalue) { return 0.; };// added only for watning pourposes in compile

    /** See father class */
    const ElemType getValue(const std::string& fullpath, const ElemType defaultvalue) { return EDGE2; };// added only for watning pourposes in compile

    /** See father class */
    const std::string getValue(const std::string& fullpath, const std::string defaultvalue) { return "a"; };// added only for watning pourposes in compile

    /** See father class */
    const FEOrder getValue(const std::string& fullpath, const FEOrder defaultvalue) { return ZERO; };// added only for watning pourposes in compile

    /** See father class */
    const MgType getValue(const std::string& fullpath, const MgType defaultvalue) { return F_CYCLE; };// added only for watning pourposes in compile

    /** See father class */
    const std::string getValueFromArray(const std::string& basepath, const unsigned int index,
                                        const std::string entry, const std::string defaultvalue) { return "a"; };// added only for watning pourposes in compile

    /** See father class */
    const BDCType getValueFromArray(const std::string& basepath, const unsigned int index,
                                    const std::string entry, const BDCType defaultvalue) {return DIRICHLET;};// added only for watning pourposes in compile

};

inline HtmlInputParser::HtmlInputParser(const std::string& infile) : InputParser(infile) {

}

inline HtmlInputParser::~HtmlInputParser() {

}

} //end namespace femus


#endif
