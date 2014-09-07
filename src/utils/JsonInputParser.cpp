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
 

bool JsonInputParser::isTrue(const std::string& base, const std::string& type) {
   
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
  if(pippo == defaultvalue) {
    std::cout << "warning: " << fullpath << " has no value or the variable name " 
      << str_i << " is not present! default value " << defaultvalue << " will be used " << std::endl;  
  }
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
  if(pippo == defaultvalue) {
    std::cout << "warning: " << fullpath << " has no value or the variable name " 
      << str_i << " is not present! default value " << defaultvalue << " will be used " << std::endl;  
  }
  value = pippo.asDouble();
  
  return value;
}
 
 
 
}