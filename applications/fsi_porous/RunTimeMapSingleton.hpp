#ifndef _rtimemap_h_
#define _rtimemap_h_

//============================================================
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <map>
#include <string>

// ===============================================
//due to the implementation of the getInstance function, this class is instantiated only once, so the idea
//is to put all the initializations in that function (doing const static would be even better and safer)
// chiaramente, per avere un'UNICA istanziazione devo mettere tutti i dati STATICI,
// il costruttore privato,
// e una funzione pubblica chiamata getInstance();
// e' chiaro che se ho anche dei dati o dei metodi non-statici allora le singole istanziazioni cambieranno,
// ma la parte statica rimane tale.
// Ora, ogni istanziazione viene costruita dentro la funzione getInstance(), quindi tale funzione avra'
// almeno tutti gli argomenti del costruttore

//e' chiaro che essendo un'implementazione singola non puoi avere valori specifici


// // // class Singleton
// // // {
// // // private:
// // //   
// // //     static unsigned int counter; //static because shared
// // //     static bool instanceFlag;  //this must be static (doesn't mean that it's not gonna be changed), because it must be SHARED
// // //     static Singleton *single;  //also, this is static because it must be SHARED
// // //     Singleton()
// // //     {
// // //         //private constructor, so that the class can only be instantiated with getInstance()
// // //     }
// // // public:
// // //     static Singleton* getInstance();
// // //     void method();
// // //     ~Singleton()
// // //     {
// // //         instanceFlag = false;
// // //     }
// // // };


//can i use a template string as a nontype template parameter? No:
//  http://stackoverflow.com/questions/5547852/string-literals-not-allowed-as-non-type-template-parameters

template <typename T>
class RunTimeMap   {

private:
  
  static bool instanceFlag;  
  static RunTimeMap<T> * single;

  static std::map<std::string, T > _rtmap;
  static std::string _file_name;  // Files&   _files;    //for the paths of the files
  std::string _tag_name;          // for the tag  -   instantiation specific => non static... mmmh, i don't think i can do something like that

   RunTimeMap( std::string class_name, std::string file_name );
   
public:
    static RunTimeMap<T> * getInstance(std::string class_name, std::string file_name );
  
  ~RunTimeMap();

  void read();
  void print() const;

  void   set(const std::string & name, T & value);
  T      get(const std::string & name) const;

};


//static member with template... cannot put it in a .cpp file, must leave everything in the header...
//so i will leave this in the header
//INITIALIZATION of static data members ===========
template <typename T>  bool RunTimeMap<T>::instanceFlag = false;
template <typename T>  RunTimeMap<T>* RunTimeMap<T>::single = NULL;
template <typename T>  std::map<std::string, T > RunTimeMap<T>::_rtmap;
template <typename T>  std::string  RunTimeMap<T>::_file_name;
// template <typename T>  std::string  RunTimeMap<T>::_tag_name;



template <typename T>
 RunTimeMap<T> * RunTimeMap<T>::getInstance(std::string class_name, std::string file_name ) {
   
      if(! instanceFlag)
    {
        instanceFlag = true;
        single = new RunTimeMap<T> (class_name,file_name);
        return single;
    }
    else
    {
        return single;
    }
    
   
 }





//=============
template <typename T>
RunTimeMap<T>::RunTimeMap(std::string class_name_in,std::string file_name_in)  : 
_tag_name(class_name_in)
{ 

  _file_name = file_name_in; 
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
    std::string abs_path = _file_name; // _files.get_femus_dir() + "/applications/" + _files.get_myapp_name() + "/";  
    filename << abs_path; //<< _files.get("CONFIG_DIR") << "/" << _files.get("BASEPARUTILS");  

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
	   
                  if (buf == "#") getline(fin, buf); // comment line
                else if (buf == cl_end) {break;} //exit the INNER nearest enclosing loop 
                else { fin >> value; set(buf,value); } // set new parameter
   
	 } while (buf != cl_end);

      } 
      
    else{ std::cout << _tag_name << " ::read_par (RunTimeMap): no " << cl_begin << " field found" << std::endl; abort();}

  }
  else {std::cout <<  _tag_name << " ::read_par (RunTimeMap): no parameter file found" << std::endl; abort();}
 
return;  
 
    
  }





// ============================================================
/// This function prints all the  parameters to stdout

//need typename because  of dependent scope... TODO
template <typename T>
void RunTimeMap<T>::print() const {
  
// #ifdef PRINT_INFO

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
  
// #endif  
  
  return;
}


#endif
