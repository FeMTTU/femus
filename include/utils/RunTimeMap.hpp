#ifndef _rtimemap_h_
#define _rtimemap_h_


#include <cstdlib>
#include <sstream>
#include <fstream>  
#include <iostream>
#include <iomanip>
#include <map>
#include <string>

#include "FemusDefault.hpp"

//can i use a template string as a nontype template parameter? No:
//  http://stackoverflow.com/questions/5547852/string-literals-not-allowed-as-non-type-template-parameters

template <typename T>
class RunTimeMap   {

public:
  
  std::map<std::string, T > _rtmap;
  const std::string               _tag_name;  //for the tag
  const std::string               _basepath;


   RunTimeMap(const std::string class_name, const std::string basepath_in );
  ~RunTimeMap();

  void read();         ///< Reading file names from file
  void print() const;

  void   set(const std::string & name, T & value);
  T      get(const std::string & name) const;

};


//=============================================
//=============================================

//=============
template <typename T>
RunTimeMap<T>::RunTimeMap(const std::string class_name_in, const std::string basepath_in) : 
_tag_name(class_name_in), 
_basepath(basepath_in)
{ 
  //_tag_name = class_name; 
  //fill it in the initialization list instead of here;
  //    you may also put it here because it is not a reference nor a const;
  //    but in the initialization list you avoid calling for an empty constructor, IIRC
  
}

template <typename T>
RunTimeMap<T>::~RunTimeMap() {}

//==================================================
//this is the syntax for templated functions
//of the class template
//if the template is not needed, then the functions
// do not need this syntax, do they? //TODO
//i'm afraid you have to do it always, even if you dont use that template
//in that function (or you do not SEEM to use that template...)
template <typename T>
 void   RunTimeMap<T>::set(const std::string & name, T & value)  {
    _rtmap.insert(make_pair(name,value));
  }

//=========================================  
template <typename T>
    T RunTimeMap<T>::get(const std::string & name) const {
    return _rtmap.find(name)->second;
  }
  

//==========================================================  
//remember that you have to pass the tag name WITHOUT < > !
//this class does not contain the FILE informations, so we must 
//pass that explicitly to the read function
//if one day we remove the passing of the Utils to everything, then we still have the Files...
// Files should be sthg like global variables..

template <typename T>
  void RunTimeMap<T>::read() {

    std::ostringstream filename;
    filename << _basepath <<  "/" << DEFAULT_CONFIGDIR << "/" << DEFAULT_RUNTIMECONF;  


std::ifstream fin(filename.str().c_str());
  std::string buf="";
  T value;

  std::string cl_begin = "<" +  _tag_name + ">";
  std::string cl_end   = "</" + _tag_name + ">";
  
  
  if (fin != NULL) {
    while (!fin.eof() && buf != cl_begin ) {
      fin >> buf;
      } 

     if (buf == cl_begin ) {

      do {
           fin >> buf;
	   
                  if (buf == "#") getline(fin, buf); // comment line  //fin.ignore(200,'\n'); TODO what is the diff
                else if (buf == cl_end) {break;} //exit the INNER nearest enclosing loop 
                else { fin >> value; set(buf,value); } // set new parameter
   
	 } while (buf != cl_end);

      } 
      
    else{ std::cout << _tag_name << " ::read_par (RunTimeMap): no " << cl_begin << " field found" << std::endl; abort();}

  }
  else {std::cout <<  _tag_name << " ::read_par (RunTimeMap): no parameter file found" << std::endl; abort();}
 
   fin.close();
 
return;  

  }





// ============================================================
/// This function prints all the  parameters to stdout

//need typename because  of dependent scope... TODO
template <typename T>
void RunTimeMap<T>::print() const {
  
#ifdef DEFAULT_PRINT_INFO

  std::cout << "\n ============================= \n";
  std::cout << _tag_name << "   RunTimeMap: " << (int)_rtmap.size() << " parameters: \n";
  
  typename std::map<std::string,T>::const_iterator   pos=_rtmap.begin();
  typename std::map<std::string,T>::const_iterator pos_e=_rtmap.end();
  for (;pos!=pos_e;pos++){
    std::string name=pos->first; // get name
    T value=pos->second;    // get value
    std::cout << "  " << std::left << std::setw(15) << name << " = " 
	      << std::setprecision(12) << value << std::endl;
  }
  
#endif  
  
  return;
}


#endif
