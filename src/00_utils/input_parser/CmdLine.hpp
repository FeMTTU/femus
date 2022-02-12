#ifndef __femus_utils_CmdLine_hpp__
#define __femus_utils_CmdLine_hpp__

#include <iostream>
#include <utility>
#include <cstdlib>
#include <cstring>
#include <sstream>


namespace femus {



//the advantage of adding a namespace is that you dont have to INSTANTIATE OBJECTS.
//So, if you need to pass data or functions to some class, you dont have to INSTANTIATE an EXTERNAL object
//containing those data or those functions and then PASS that EXTERNAL OBJECT to the target object:
//those data and functions are like GLOBAL, you can pass them just with an include

// TODO this class should be with a map of STRINGS which then are turned into DOUBLES or INTEGERS according to the NEEDS!!!
// OK I CHANGED IT, it might have comments for the DOUBLES but it is with STRINGS NOW!
// We should make a TEMPLATE CLASS or TEMPLATE FUNCTION

namespace {
  
  std::map<std::string,std::string> _cmdline;///< Parameter map  //
//in this position it is like in a class, namely, a place where i can "drop" (once,or more) and then "drag"
//also, it is as if it were PRIVATE because you cannot access it from other files than this
//=> you necessarily need a get() like for protected data in classes
// how about const attributes for functions? It's only for classes:
// error: non-member function ‘double CmdLine::get(std::string)’ cannot have cv-qualifier
//to protect this datum I could use const... but in that way i couldnt fill it after initialization, so no...
//since this map is in anonymous namespace, it has file scope, so it can be seen only in this file

//what if I do a map of strings and strings and then I transform them into int or double on the basis of what i need?

}

namespace CmdLine {
  
/////////
inline void parse(int argc, char** argv) {
//if you have a class with one function and one datum, then put that datum in that function and avoid doing a class. 
//In that class data are not so important, they are not modified... well, that datum is important because you 
// use it once and then you store things in it without refilling every time.
//But the datum is not so important in the sense that you dont have DIFFERENT INSTANTIATION WITH CHANGING DATA,
//it's STATIC, you dont fill it in different ways, you fill it once and for all and then you only USE it.
//so, we choose to put it in an anonymous namespace whose scope is this file

  for (int i=0;i<argc;i++) {

    std::cout << argv[i] << std::endl; 

    //if the first two characters are --, then set that pair in the map
    
   if   (!strncmp(argv[i],"--",2) )      { 

     //here you should first check that argv[i+1] is a valid double, OR integer
     
     //at this point i am sure that argv[i] is a char while argv[i+1] is a double
     //but both of them are char* in principle
     //so I have to find a way to convert them to their respective type in the map
     //char* to string
     //char* to double
     
//argv[i] == "--dt" The overloading of this operator does not exist!
// /home/femus/Software/femus/src/CmdLineParse.C:25:24: warning: comparison with string literal results in unspecified behaviour
//     _cmdline.insert( std::make_pair(argv[i],/*strtod(*/argv[i+1]/*,NULL)*/) );  //this gives instantiated from here!
//     double temp  = strtod(argv[i+1],NULL);  gives a problem in parallel

     std::string l_val(argv[i]);

     //     double rightval = strtod(argv[i+1],NULL);
     std::stringstream ss_right; 
     ss_right << argv[i+1];
//      double r_val = 0;
     std::string r_val("");
     ss_right >> r_val;
       
      _cmdline.insert( std::make_pair(l_val,r_val ) );
      
      std::cout << "****Command line " << l_val << " " << _cmdline.find(l_val)->second << std::endl ;   
      
   }

  }

return;
}
  
  
  
////////
std::string get(std::string name) {
  
  if (_cmdline.find(name) == _cmdline.end()) { std::cout << "Flag " << name << " missing, please provide it" << std::endl; abort(); }
  
 return _cmdline.find(name)->second; 

}
  
  
  
  
}//end namespace







} //end namespace femus



#endif