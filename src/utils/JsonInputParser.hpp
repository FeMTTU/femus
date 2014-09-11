/*=========================================================================

Program: FEMUS
Module: JsonInputParser
Authors: Simone Bnà

Copyright (c) FEMUS
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __jsoninputparser_hpp__
#define __jsoninputparser_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "InputParser.hpp"
#include <json/json.h>
#include <json/value.h>  

#include "iostream"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------

namespace femus {

/**
* This class is the concrete class for Json format input parsing.
*/

class JsonInputParser : public InputParser {

public:

    /** Constructor */
    JsonInputParser(const std::string& infile);
    
    /** Destructor. */
    ~JsonInputParser ();

    void printInfo() {std::cout << "I'm Json" << std::endl;};
    
    const bool isTrue(const std::string& base, const std::string& type);
    
    const unsigned int getSize(const std::string& base);
    
    const int getValue(const std::string& fullpath, const int defaultvalue);
    
    const double getValue(const std::string& fullpath, const double defaultvalue);
    
    const ElemType getValue(const std::string& fullpath, const ElemType defaultvalue);
    
    const std::string getValue(const std::string& fullpath, const std::string defaultvalue);
    
    const FEOrder getValue(const std::string& fullpath, const FEOrder defaultvalue);
    
    const MgType getValue(const std::string& fullpath, const MgType defaultvalue);
    
    const std::string getValueFromArray(const std::string& basepath, const unsigned int index, 
					const std::string entry, const std::string defaultvalue);
    
    const BDCType getValueFromArray(const std::string& basepath, const unsigned int index, 
			         const std::string entry, const BDCType defaultvalue);
    
private:    
    
  Json::Value _root;   //< will contains the root value after parsing.
  Json::Reader _reader;
    
};

} //end namespace femus


#endif