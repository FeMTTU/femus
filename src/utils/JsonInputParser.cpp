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

JsonInputParser::JsonInputParser(const std::string& infile) : InputParser(infile) {

    bool parsingSuccessful = _reader.parse(_infile, _root );
    if ( !parsingSuccessful )
    {
        // report to the user the failure and their locations in the document.
        std::cout  << "Failed to parse configuration\n" << _reader.getFormatedErrorMessages();
        exit(1);
    }

}

JsonInputParser::~JsonInputParser() {

}


const bool JsonInputParser::isTrue(const std::string& base, const std::string& type) {

    std::string str1 = extractWord(base, 1);
    int numWords = howManyWords(base);

    Json::Value pippo = _root[str1];

    for(int i=2; i<=numWords; i++) {

        std::string str_i = extractWord(base, i);

        if(pippo.get(str_i,false).asBool()) {
            pippo = pippo.get(str_i,str_i);
        }
        else
        {
            std::cerr << "Error: json input parsing failed" << std::endl;
            std::cerr << "The entry path is: " << base << std::endl;
            std::cerr << "The parser didn't find the word: " << str_i << std::endl;
            std::cerr << "Check the input file for misprinting error!" << std::endl;
            exit(1);
        }
    }

    bool boolflag = pippo.get(type, false).asBool();

    return boolflag;

}


const unsigned int JsonInputParser::getSize(const std::string& fullpath) {
  
    unsigned int size;
    std::string str1 = extractWord(fullpath, 1);
    std::string str_i;
    int numWords = howManyWords(fullpath);

    Json::Value pippo = _root[str1];

    for(int i=2; i<numWords; i++) {

        str_i = extractWord(fullpath, i);

        if(pippo.get(str_i,false).asBool()) {
            pippo = pippo.get(str_i,str_i);
        }
        else
        {
            std::cerr << "Error: json input parsing failed" << std::endl;
            std::cerr << "The entry path is: " << fullpath << std::endl;
            std::cerr << "The parser didn't find the word: " << str_i << std::endl;
            std::cerr << "Check the input file for misprinting error!" << std::endl;
            exit(1);
        }
    }

    str_i = extractWord(fullpath, numWords);
    pippo = pippo.get(str_i,"");
    
    size = pippo.size();
    return size;
}



const int JsonInputParser::getValue(const std::string& fullpath, const int defaultvalue) {

    int value = defaultvalue;

    std::string str1 = extractWord(fullpath, 1);
    std::string str_i;
    int numWords = howManyWords(fullpath);

    Json::Value pippo = _root[str1];

    for(int i=2; i<numWords; i++) {

        str_i = extractWord(fullpath, i);

        if(pippo.get(str_i,false).asBool()) {
            pippo = pippo.get(str_i,str_i);
        }
        else
        {
            std::cerr << "Error: json input parsing failed" << std::endl;
            std::cerr << "The entry path is: " << fullpath << std::endl;
            std::cerr << "The parser didn't find the word: " << str_i << std::endl;
            std::cerr << "Check the input file for misprinting error!" << std::endl;
            exit(1);
        }
    }

    str_i = extractWord(fullpath, numWords);
    pippo = pippo.get(str_i,defaultvalue);
    
    value = pippo.asInt();

    return value;
}

const double JsonInputParser::getValue(const std::string& fullpath, const double defaultvalue) {

    double value = defaultvalue;

    std::string str1 = extractWord(fullpath, 1);
    std::string str_i;
    int numWords = howManyWords(fullpath);

    Json::Value pippo = _root[str1];

    for(int i=2; i<numWords; i++) {

        str_i = extractWord(fullpath, i);

        if(pippo.get(str_i,false).asBool()) {
            pippo = pippo.get(str_i,str_i);
        }
        else
        {
            std::cerr << "Error: json input parsing failed" << std::endl;
            std::cerr << "The entry path is: " << fullpath << std::endl;
            std::cerr << "The parser didn't find the word: " << str_i << std::endl;
            std::cerr << "Check the input file for misprinting error!" << std::endl;
            exit(1);
        }
    }

    str_i = extractWord(fullpath, numWords);
    pippo = pippo.get(str_i,defaultvalue);
    
    value = pippo.asDouble();

    return value;
}


const ElemType JsonInputParser::getValue(const std::string& fullpath, const ElemType defaultvalue) {
  
    ElemType value = defaultvalue;
    std::string valuestr = "";

    std::string str1 = extractWord(fullpath, 1);
    std::string str_i;
    int numWords = howManyWords(fullpath);

    Json::Value pippo = _root[str1];

    for(int i=2; i<numWords; i++) {

        str_i = extractWord(fullpath, i);

        if(pippo.get(str_i,false).asBool()) {
            pippo = pippo.get(str_i,str_i);
        }
        else
        {
            std::cerr << "Error: json input parsing failed" << std::endl;
            std::cerr << "The entry path is: " << fullpath << std::endl;
            std::cerr << "The parser didn't find the word: " << str_i << std::endl;
            std::cerr << "Check the input file for misprinting error!" << std::endl;
            exit(1);
        }
    }

    str_i = extractWord(fullpath, numWords);
    pippo = pippo.get(str_i,"Quad9");

    valuestr = pippo.asString();

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


const std::string JsonInputParser::getValue(const std::string& fullpath, const std::string defaultvalue) {
    std::string value = "";

    std::string str1 = extractWord(fullpath, 1);
    std::string str_i;
    int numWords = howManyWords(fullpath);

    Json::Value pippo = _root[str1];

    for(int i=2; i<numWords; i++) {

        str_i = extractWord(fullpath, i);

        if(pippo.get(str_i,false).asBool()) {
            pippo = pippo.get(str_i,str_i);
        }
        else
        {
            std::cerr << "Error: json input parsing failed" << std::endl;
            std::cerr << "The entry path is: " << fullpath << std::endl;
            std::cerr << "The parser didn't find the word: " << str_i << std::endl;
            std::cerr << "Check the input file for misprinting error!" << std::endl;
            exit(1);
        }
    }

    str_i = extractWord(fullpath, numWords);
    pippo = pippo.get(str_i,defaultvalue);
    value = pippo.asString();
   
    return value;
   
}


const FEOrder JsonInputParser::getValue(const std::string& fullpath, const FEOrder defaultvalue) {
  
    FEOrder value = defaultvalue;
    std::string valuestr = "";

    std::string str1 = extractWord(fullpath, 1);
    std::string str_i;
    int numWords = howManyWords(fullpath);

    Json::Value pippo = _root[str1];

    for(int i=2; i<numWords; i++) {

        str_i = extractWord(fullpath, i);

        if(pippo.get(str_i,false).asBool()) {
            pippo = pippo.get(str_i,str_i);
        }
        else
        {
            std::cerr << "Error: json input parsing failed" << std::endl;
            std::cerr << "The entry path is: " << fullpath << std::endl;
            std::cerr << "The parser didn't find the word: " << str_i << std::endl;
            std::cerr << "Check the input file for misprinting error!" << std::endl;
            exit(1);
        }
    }

    str_i = extractWord(fullpath, numWords);
    pippo = pippo.get(str_i,"first");

    valuestr = pippo.asString();
    
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


const MgType JsonInputParser::getValue(const std::string& fullpath, const MgType defaultvalue) {
  
    MgType value = defaultvalue;
    std::string valuestr = "";

    std::string str1 = extractWord(fullpath, 1);
    std::string str_i;
    int numWords = howManyWords(fullpath);

    Json::Value pippo = _root[str1];

    for(int i=2; i<numWords; i++) {

        str_i = extractWord(fullpath, i);

        if(pippo.get(str_i,false).asBool()) {
            pippo = pippo.get(str_i,str_i);
        }
        else
        {
            std::cerr << "Error: json input parsing failed" << std::endl;
            std::cerr << "The entry path is: " << fullpath << std::endl;
            std::cerr << "The parser didn't find the word: " << str_i << std::endl;
            std::cerr << "Check the input file for misprinting error!" << std::endl;
            exit(1);
        }
    }

    str_i = extractWord(fullpath, numWords);
    pippo = pippo.get(str_i,"V_cycle");

    valuestr = pippo.asString();
  

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

const std::string JsonInputParser::getValueFromArray(const std::string& basepath, const unsigned int index, 
					             const std::string entry, const std::string defaultvalue) {
      
    std::string value = defaultvalue;  
    std::string str1 = extractWord(basepath, 1);
    std::string str_i;
    int numWords = howManyWords(basepath);

    Json::Value pippo = _root[str1];

    for(int i=2; i<numWords; i++) {

        str_i = extractWord(basepath, i);

        if(pippo.get(str_i,false).asBool()) {
            pippo = pippo.get(str_i,str_i);
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
    pippo = pippo.get(str_i,"");
  
    pippo = pippo[index].get(entry,defaultvalue);
    value = pippo.asString();
  
    return value;
}

const BDCType JsonInputParser::getValueFromArray(const std::string& basepath, const unsigned int index, 
					 const std::string entry, const BDCType defaultvalue) {
  
    BDCType value = defaultvalue;  
    std::string valuestr = "";
    std::string str1 = extractWord(basepath, 1);
    std::string str_i;
    int numWords = howManyWords(basepath);

    Json::Value pippo = _root[str1];

    for(int i=2; i<numWords; i++) {

        str_i = extractWord(basepath, i);

        if(pippo.get(str_i,false).asBool()) {
            pippo = pippo.get(str_i,str_i);
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
    pippo = pippo.get(str_i,"");
  
    pippo = pippo[index].get(entry,defaultvalue);
    valuestr = pippo.asString();
  
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


}
