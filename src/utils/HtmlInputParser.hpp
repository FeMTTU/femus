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
    
    const bool isTrue(const std::string& base, const std::string& type) {};
    
    const unsigned int getSize(const std::string& base) {};
    
    const int getValue(const std::string& fullpath, const int defaultvalue) {};
    
    const double getValue(const std::string& fullpath, const double defaultvalue) {};
    
    const ElemType getValue(const std::string& fullpath, const ElemType defaultvalue) {};
    
    const std::string getValue(const std::string& fullpath, const std::string defaultvalue) {};
    
    const FEOrder getValue(const std::string& fullpath, const FEOrder defaultvalue) {};
    
    const MgType getValue(const std::string& fullpath, const MgType defaultvalue) {};
    
    const std::string getValueFromArray(const std::string& basepath, const unsigned int index, 
					const std::string entry, const std::string defaultvalue) {};
					
    const BDCType getValueFromArray(const std::string& basepath, const unsigned int index, 
					 const std::string entry, const BDCType defaultvalue) {};
					 

					
  
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