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


namespace femus {



template <typename T>
class RunTimeMap   {

public:

   RunTimeMap(const std::string class_name, const std::string basepath_in );
  ~RunTimeMap();

  void read();
  void print() const;

  void   set(const std::string & name, const T & value);
  const T      get(const std::string & name) const;  //TODO do const or not for the return... sembra che l'errore in debug mode fosse proprio qui!!!

  const std::string  get_rbasepath() const {return _basepath;}
  
  std::string               _basepath; // questo e' il basepath da cui parte la COPIA dei file di input e cosi' via
                                       //TODO let me put it public just because I don't want to do the SET function right now...
                                       //TODO siccome lo posso MODIFICARE in base al RESTART non posso piu' metterlo CONST
  std::map<std::string, T > _rtmap;

private:

  const std::string               _tag_name;
  
};


//=============================================
//=============================================

//=============
template <typename T>
RunTimeMap<T>::RunTimeMap(const std::string class_name_in, const std::string basepath_in) : 
_basepath(basepath_in),
_tag_name(class_name_in)
{
  
     read();
     print();

  // I AM ALMOST SURE THAT THERE IS SOME PROBLEM WITH THE CONSTRUCTION of the _rtmap object...
  // it sounds like it's not well initialized
  //clearly the object is not const, so we need to find some other way to "sort of initialize" it
//    _rtmap.clear();

// IF I UNDERSTOOD CORRECTLY FROM THE DEBUGGER, the DEFAULT CONSTRUCTOR for STD::MAP is called here !!!  

  
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
 void   RunTimeMap<T>::set(const std::string & name, const T & value)  {
    _rtmap.insert(make_pair(name,value));
  }

//=========================================  
template <typename T>
  const  T RunTimeMap<T>::get(const std::string & name) const {
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



} //end namespace femus



#endif
