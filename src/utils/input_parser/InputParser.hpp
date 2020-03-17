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

#ifndef __femus_utils_InputParser_hpp__
#define __femus_utils_InputParser_hpp__


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
* This class is an abstract class for input parsing. It implements the Factory Method pattern. Using the static \p build function
* you can istantiate the concrete classes HTMLInputParser (not yet implemented) or JsonInputParser.
* 
* How to use it:
* the path of the variable in the json input file is expressed as a sequence of words separated by a point
* e.g. :
* "mesh.first.type.box"
* 
* use the getValue function to read and get the value
* use the getSize function to get the size of a vector
* use the isTrue function to check if the basepath.type value in the file is equal or not to the type value you ask to   
*/

class InputParser {

public:

    /** Constructor */
    InputParser(const std::string& infile);

    /** Destructor */
    virtual ~InputParser();

    /** Build the concrete class using the Factory Method pattern */
    static std::unique_ptr<InputParser> build(const std::string& infile);

    /** return true if the argument type is equal to the basepath.type value in the json input file */
    virtual const bool isTrue(const std::string& basepath, const std::string& type ) = 0;

    /** return the size of the vector at the \p fullpath position in the json input file */
    virtual const unsigned int getSize(const std::string& fullpath) = 0;

    /** return the int value at the \p fullpath position in the json input file */
    virtual const int getValue(const std::string& fullpath, const int defaultvalue) = 0;

    /** return the double value at the \p fullpath position in the json input file */
    virtual const double getValue(const std::string& fullpath, const double defaultvalue) = 0;

    /** return the ElemType value at the \p fullpath position in the json input file */
    virtual const ElemType getValue(const std::string& fullpath, const ElemType defaultvalue) = 0;

    /** return the FEOrder value at the \p fullpath position in the json input file */
    virtual const FEOrder getValue(const std::string& fullpath, const FEOrder defaultvalue) = 0;

    /** return the string value at the \p fullpath position in the json input file */
    virtual const std::string getValue(const std::string& fullpath, const std::string defaultvalue) = 0;

    /** return the MgType value at the \p fullpath position in the json input file */
    virtual const MgType getValue(const std::string& fullpath, const MgType defaultvalue) = 0;

    /** return the string value at the \p fullpath position in the json input file */
    virtual const std::string getValueFromArray(const std::string& basepath, const unsigned int index,
            const std::string entry, const std::string defaultvalue) = 0;

    /** return the BDCType value at the \p fullpath position in the json input file */
    virtual const BDCType getValueFromArray(const std::string& basepath, const unsigned int index,
                                            const std::string entry, const BDCType defaultvalue) = 0;

protected:

    /** return a string that stands for the input file for the Jsoncpp library*/
    std::string readInputTestFile(const char *path);

    /** return how many words are in the input string */
    unsigned int howManyWords(const std::string& input) const;

    /** return the word in the input string at the \p iter position */
    const std::string extractWord(const std::string& input, const int iter) const;

    std::string _infile; //< the input file


};

} //end namespace femus


#endif
