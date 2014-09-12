/*=========================================================================

  Program: FEMUS
  Module: JsonInputParser
  Authors: Simone Bnà

  Copyright (c) FEMTTU
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/


//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include "JsonInputParser.hpp"
#include "cmath"

// C++ includes


namespace femus {

//-----------------------------------------------------------------------------------------
JsonInputParser::JsonInputParser(const std::string& infile) : InputParser(infile) {

    bool parsingSuccessful = _reader.parse(_infile, _root );
    if ( !parsingSuccessful )
    {
        // report to the user the failure and their locations in the document.
        std::cout  << "Failed to parse configuration\n" << _reader.getFormatedErrorMessages();
        exit(1);
    }

}

//-----------------------------------------------------------------------------------------
JsonInputParser::~JsonInputParser() {

}

//-----------------------------------------------------------------------------------------
const bool JsonInputParser::isTrue(const std::string& base, const std::string& type) {

    Json::Value jsonvalue = getJsonValue(base, "");

    bool boolflag = jsonvalue.get(type, false).asBool();

    return boolflag;

}

//-----------------------------------------------------------------------------------------
const unsigned int JsonInputParser::getSize(const std::string& fullpath) {

    unsigned int size;
    Json::Value jsonvalue = getJsonValue(fullpath, 0);

    size = jsonvalue.size();
    return size;
}

//-----------------------------------------------------------------------------------------
const int JsonInputParser::getValue(const std::string& fullpath, const int defaultvalue) {

    int value = defaultvalue;
    Json::Value jsonvalue = getJsonValue(fullpath, defaultvalue);

    value = jsonvalue.asInt();

    return value;
}

//-----------------------------------------------------------------------------------------
const double JsonInputParser::getValue(const std::string& fullpath, const double defaultvalue) {

    double value = defaultvalue;
    Json::Value jsonvalue = getJsonValue(fullpath, defaultvalue);

    value = jsonvalue.asDouble();

    return value;
}

//-----------------------------------------------------------------------------------------
const ElemType JsonInputParser::getValue(const std::string& fullpath, const ElemType defaultvalue) {

    ElemType value = defaultvalue;
    std::string valuestr = "";
    Json::Value jsonvalue = getJsonValue(fullpath, "");

    valuestr = jsonvalue.asString();

    if(valuestr == "Quad9")
    {
        value = QUAD9;
    }
    else if(valuestr == "Tri6")
    {
        value = TRI6;
    }
    else if(valuestr == "Edge3")
    {
        value = EDGE3;
    }
    else if(valuestr == "Hex27")
    {
        value = HEX27;
    }
    else
    {
        value = INVALID_ELEM;
    }

    return value;

}

//-----------------------------------------------------------------------------------------
const std::string JsonInputParser::getValue(const std::string& fullpath, const std::string defaultvalue) {

    std::string value = "";
    Json::Value jsonvalue = getJsonValue(fullpath, defaultvalue);

    value = jsonvalue.asString();

    return value;

}

//-----------------------------------------------------------------------------------------
const FEOrder JsonInputParser::getValue(const std::string& fullpath, const FEOrder defaultvalue) {

    FEOrder value = defaultvalue;
    std::string valuestr = "";
    Json::Value jsonvalue = getJsonValue(fullpath, "");

    valuestr = jsonvalue.asString();

    if (valuestr == "first")
    {
        value = FIRST;
    }
    else if (valuestr == "serendipity")
    {
        value = SERENDIPITY;
    }
    else if (valuestr == "second")
    {
        value = SECOND;
    }
    else
    {
        std::cerr << " Error: finite element order not supported!" << std::endl;
        exit(1);
    }

    return value;
}

//-----------------------------------------------------------------------------------------
const MgType JsonInputParser::getValue(const std::string& fullpath, const MgType defaultvalue) {

    MgType value = defaultvalue;
    std::string valuestr = "";
    Json::Value jsonvalue = getJsonValue(fullpath, "");

    valuestr = jsonvalue.asString();

    if ( valuestr == "V_cycle" )
    {
        value = V_CYCLE;
    }
    else if( valuestr == "F_cycle" )
    {
        value = F_CYCLE;
    }
    else {
        std::cout << "The selected MG cycle does not exist!" << std::endl;
        exit(1);
    }

    return value;

}

//-----------------------------------------------------------------------------------------
const std::string JsonInputParser::getValueFromArray(const std::string& basepath, const unsigned int index,
        const std::string entry, const std::string defaultvalue) {

    std::string value = defaultvalue;
    Json::Value jsonvalue = getJsonValue(basepath, "");

    jsonvalue = jsonvalue[index].get(entry,defaultvalue);
    value = jsonvalue.asString();

    return value;
}

//-----------------------------------------------------------------------------------------
const BDCType JsonInputParser::getValueFromArray(const std::string& basepath, const unsigned int index,
        const std::string entry, const BDCType defaultvalue) {

    BDCType value = defaultvalue;
    std::string valuestr = "";
    Json::Value jsonvalue = getJsonValue(basepath, "");

    jsonvalue = jsonvalue[index].get(entry,defaultvalue);
    valuestr = jsonvalue.asString();

    if ( valuestr == "dirichlet") {
        value = DIRICHLET;
    }
    else if( valuestr == "neumann") {
        value = NEUMANN;
    }
    else {
        std::cout << "Boundary condition not implemented!" << std::endl;
        exit(1);
    }

    return value;

}

//-----------------------------------------------------------------------------------------
template <typename T_sys>
const Json::Value JsonInputParser::getJsonValue(const std::string& basepath, const T_sys defaultvalue) {

    std::string str1 = extractWord(basepath, 1);
    std::string str_i;
    int numWords = howManyWords(basepath);

    Json::Value jsonvalue = _root[str1];

    for(int i=2; i<numWords; i++) {

        str_i = extractWord(basepath, i);

        if(jsonvalue.get(str_i,false).asBool()) {
            jsonvalue = jsonvalue.get(str_i,str_i);
        }
        else
        {
            std::cerr << "Error: json input parsing failed" << std::endl;
            std::cerr << "The entry path is: " << basepath << std::endl;
            std::cerr << "The parser didn't find the word: " << str_i << std::endl;
            std::cerr << "Check the input file for misprinting error!" << std::endl;
            exit(1);
        }
    }

    str_i = extractWord(basepath, numWords);
    jsonvalue = jsonvalue.get(str_i,defaultvalue);

    return jsonvalue;
}


}
