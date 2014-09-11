/*=========================================================================

Program: FEMUS
Module: InputParser
Authors: Simone Bnà

Copyright (c) FEMUS
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __inputparser_hpp__
#define __inputparser_hpp__


//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <memory>
#include <string>
#include <iostream>
#include <cstdlib>
#include "ElemTypeEnum.hpp"
#include "FElemTypeEnum.hpp"
#include "MgTypeEnum.hpp"
#include "BDCTypeEnum.hpp"

//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------

namespace femus {
  
/**
* This class is an abstract class for input parsing.
*/

class InputParser {

public:

    /** Constructor */
    InputParser(const std::string& infile);
    
    /** Destructor */
    virtual ~InputParser();
    
    /** Builds a \p a using the b in \p c */
    static std::auto_ptr<InputParser> build(const std::string& infile, const int format);
    
    virtual void printInfo() = 0;
    
    virtual const bool isTrue(const std::string& basepath, const std::string& type ) = 0;
    
    virtual const unsigned int getSize(const std::string& fullpath) = 0;
    
    virtual const int getValue(const std::string& fullpath, const int defaultvalue) = 0;
    
    virtual const double getValue(const std::string& fullpath, const double defaultvalue) = 0;
    
    virtual const ElemType getValue(const std::string& fullpath, const ElemType defaultvalue) = 0;
    
    virtual const FEOrder getValue(const std::string& fullpath, const FEOrder defaultvalue) = 0;
    
    virtual const std::string getValue(const std::string& fullpath, const std::string defaultvalue) = 0;
    
    virtual const MgType getValue(const std::string& fullpath, const MgType defaultvalue) = 0;
    
    virtual const std::string getValueFromArray(const std::string& basepath, const unsigned int index, 
					        const std::string entry, const std::string defaultvalue) = 0;
					   
    virtual const BDCType getValueFromArray(const std::string& basepath, const unsigned int index, 
 					 const std::string entry, const BDCType defaultvalue) = 0;
					 
protected:
  
    std::string readInputTestFile(const char *path);    
    
    unsigned int howManyWords(const std::string& input) const; 
    
    const std::string extractWord(const std::string& input, const int iter) const;
    
    std::string _infile;
    
    
};

} //end namespace femus


#endif