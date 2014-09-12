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

#ifndef __htmlinputparser_hpp__
#define __htmlinputparser_hpp__

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
    const bool isTrue(const std::string& base, const std::string& type) {};

    /** See father class */
    const unsigned int getSize(const std::string& base) {};

    /** See father class */
    const int getValue(const std::string& fullpath, const int defaultvalue) {};

    /** See father class */
    const double getValue(const std::string& fullpath, const double defaultvalue) {};

    /** See father class */
    const ElemType getValue(const std::string& fullpath, const ElemType defaultvalue) {};

    /** See father class */
    const std::string getValue(const std::string& fullpath, const std::string defaultvalue) {};

    /** See father class */
    const FEOrder getValue(const std::string& fullpath, const FEOrder defaultvalue) {};

    /** See father class */
    const MgType getValue(const std::string& fullpath, const MgType defaultvalue) {};

    /** See father class */
    const std::string getValueFromArray(const std::string& basepath, const unsigned int index,
                                        const std::string entry, const std::string defaultvalue) {};

    /** See father class */
    const BDCType getValueFromArray(const std::string& basepath, const unsigned int index,
                                    const std::string entry, const BDCType defaultvalue) {};

};

inline HtmlInputParser::HtmlInputParser(const std::string& infile) : InputParser(infile) {

}

inline HtmlInputParser::~HtmlInputParser() {

}

} //end namespace femus


#endif
