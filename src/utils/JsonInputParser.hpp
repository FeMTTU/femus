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

#ifndef __femus_utils_JsonInputParser_hpp__
#define __femus_utils_JsonInputParser_hpp__

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

    /** See father class */
    const bool isTrue(const std::string& base, const std::string& type);

    /** See father class */
    const unsigned int getSize(const std::string& base);

    /** See father class */
    const int getValue(const std::string& fullpath, const int defaultvalue);

    /** See father class */
    const double getValue(const std::string& fullpath, const double defaultvalue);

    /** See father class */
    const ElemType getValue(const std::string& fullpath, const ElemType defaultvalue);

    /** See father class */
    const std::string getValue(const std::string& fullpath, const std::string defaultvalue);

    /** See father class */
    const FEOrder getValue(const std::string& fullpath, const FEOrder defaultvalue);

    /** See father class */
    const MgType getValue(const std::string& fullpath, const MgType defaultvalue);

    /** See father class */
    const std::string getValueFromArray(const std::string& basepath, const unsigned int index,
                                        const std::string entry, const std::string defaultvalue);

    /** See father class */
    const BDCType getValueFromArray(const std::string& basepath, const unsigned int index,
                                    const std::string entry, const BDCType defaultvalue);

private:

    /** To be added */
    template <typename T_sys>
    const Json::Value getJsonValue(const std::string& basepath, const T_sys defaultvalue);

    Json::Value _root;     //< will contains the root value after parsing.
    Json::Reader _reader;  //< will contain the json file reader

};

} //end namespace femus


#endif
