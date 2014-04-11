#ifndef _mg_files_h_
#define _mg_files_h_

// C++
#include <map>
#include <string>
#include <fstream>

#include "FemusDefault.hpp"
#include "RunTimeMap.hpp"

// =======================================
//    Files Class
// This class handles all the interactions of an application with the FILESYSTEM.
// It takes care of the mesh filename
// It takes care of generating the output directory for every instant of time
// It takes care of the RESTART procedure
// 
// =======================================
class Files {

public:

  std::string  _input_path; //this is where the input files are located BEFORE YOU COPY THEM to the OUTTIME DIR!!!!! it has to alternatives in case of restart or not
  std::string _output_path; //this is the OUTTIME DIR!!! it's always the same
  std::string    _app_path; //path of the application
  bool _restart_flag; 
  
   Files(const std::string &/*  = DEFAULT_BASEPATH*/);  //TODO seems like it doesn't work with ONE DEFAULT PARAMETER
  ~Files();


// RunTimeMap ===================
  inline       RunTimeMap<std::string> * get_frtmap_ptr()       {return  &_frtmap;} 
  inline       RunTimeMap<std::string> & get_frtmap()       {return  _frtmap;}   //trying to return the pointer instead of the reference //non-const version //WHO DECIDES whether to use THIS FUNCTION or the OTHER ONE? TODO OVERLOADING
  inline const RunTimeMap<std::string> & get_frtmap() const {return  _frtmap;}    // I WANT THIS TO RETURN a REFERENCE, because this is going to call the READ FUNCTION which modifies the object //THIS FUNCTION WILL BE CALLED by functions inside the class basically, which have the "const this" as well
  inline const std::string             get_basepath() const {return  _frtmap.get_rbasepath(); }  //THIS RETURNS a COPY

// Directory management =========
         void CheckIODirectories();
  static void CheckDir(const std::string& dir_name_in, const std::string& my_name_in);
  static void CheckDirOrAbort(const std::string& dir_name_in, const std::string& my_name_in);
  
// Copy ========================= 
  void CopyInputFiles() const;

// Restart ======================
  void ConfigureRestart();
  void PrintRunForRestart(const std::string run_name_in) const;
  
// LOG ==========================
// Stream redirect to file ======
         void RedirectCout(std::streambuf* sbuf,  std::ofstream& file_in) const;
  static void RedirectCoutFinalize(std::streambuf* sbuf);

private:
  
  RunTimeMap<std::string>           _frtmap;   //this map cannot be declared as CONST because at some point it is FILLED and it is not in the initialization in the constructor!!!
  
// Directory management
         void ComposeOutdirName();
  static void CheckDirOrMake(const std::string& dir_name_in, const std::string& my_name_in);

// Copy 
  void CopyFile(std::string  f_in,std::string  f_out) const;

  
};


#endif