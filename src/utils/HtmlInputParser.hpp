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
* This class is the concrete class for Html format input parsing.
*/

class HtmlInputParser : public InputParser {

public:

    /** Constructor */
    HtmlInputParser(const std::string& infile);
    
    /** Destructor. */
    ~HtmlInputParser ();

    void printInfo() {std::cout << "I'm Html" << std::endl;};
    
    bool isTrue(const std::string& base, const std::string& type) {};
    
    const int getValue(const std::string& fullpath, const int defaultvalue) {};
    
    const double getValue(const std::string& fullpath, const double defaultvalue) {};
};


/**
 * -------------------- inline functions ---------------------
 */

inline HtmlInputParser::HtmlInputParser(const std::string& infile) : InputParser(infile) {

}

inline HtmlInputParser::~HtmlInputParser() {

}


} //end namespace femus


#endif