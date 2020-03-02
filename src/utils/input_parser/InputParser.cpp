/*=========================================================================

  Program: FEMUS
  Module: InputParser
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
#include "InputParser.hpp"
#include "JsonInputParser.hpp"
#include "HtmlInputParser.hpp"
#include <iostream>
#include <cstdlib>
#include <fstream>

namespace femus {
  
InputParser::InputParser(const std::string& infile) {
  
    _infile = readInputTestFile( infile.c_str() );
    if ( _infile.empty() )
    {
        std::cout << "Failed to read input or empty input: " << infile.c_str() << std::endl;
        exit(1);
    }

}

InputParser::~InputParser() {

}  

//----------------------------------------------------------------------------------------------
std::unique_ptr<InputParser> InputParser::build(const std::string& infile) {
   
   int format = 0;
   
   if(infile.rfind(".json") < infile.size()) {
     format = 0;
   }
   else if (infile.rfind(".html") < infile.size()) {
     format = 1; 
   }
   else {
    std::cerr << " ERROR: Unrecognized file extension: " << infile
	      << "\n   I understand the following:\n\n"
	      << "     *.json -- JSon File\n"
	      << "     *.html -- HTML File (not yet implemented)\n"
              << std::endl;
     abort();
   }
  
   switch (format)  {

    case 0:  {
      std::unique_ptr<InputParser> ap(new JsonInputParser(infile));
      return ap;
    }
    case 1:  {
      std::cerr << "Error: HTMLInput parser not yet implemented! " << std::endl;
      abort();
    } 
    default:
      std::cerr << "ERROR:  format not implemented: " << std::endl;
      abort();
    }
  
  }

//----------------------------------------------------------------------------------------------
std::string InputParser::readInputTestFile(const char *path)
{
    FILE *file = fopen( path, "rb" );
    if ( !file )
        return std::string("");
    fseek( file, 0, SEEK_END );
    long size = ftell( file );
    fseek( file, 0, SEEK_SET );
    std::string text;
    char *buffer = new char[size+1];
    buffer[size] = 0;
    if ( fread( buffer, 1, size, file ) == (unsigned long)size )
        text = buffer;
    fclose( file );
    delete[] buffer;
    return text;
}

unsigned int InputParser::howManyWords(const std::string& input) const {
  
  const char* it = input.c_str();
  int numWords = 1;
  int length = input.length();
  
  for(int i=0; i<length; i++) {
    if(it[i] == '.') {
      numWords++;
    }
  }
  
  return numWords;
  
}


const std::string InputParser::extractWord(const std::string& input, const int iter) const {
    
  std::string strword = "";
  int length = input.length();
  const char* it = input.c_str();
  int numWords = 1;
  
  
  for(int i=0; i<length; i++) {
    if(numWords == iter && it[i] != '.') {
      strword.push_back(it[i]);    
    }
    if(it[i] == '.') {
      numWords++;
    }
  }
 
  return strword;
}




}
