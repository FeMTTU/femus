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
  std::string _output_path;
  std::string _output_time; //this is the OUTTIME DIR!!!
  std::string    _app_path; //path of the application

  bool _restart_flag; 
  
   Files(const std::string &/*  = DEFAULT_BASEPATH*/);  //TODO seems like it doesn't work with ONE DEFAULT PARAMETER
  ~Files();

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
         void log_petsc() const;

private:
  
// Directory management
         void ComposeOutdirName();
  static void CheckDirOrMake(const std::string& dir_name_in, const std::string& my_name_in);

// Copy 
  void CopyFile(std::string  f_in,std::string  f_out) const;
  
};


#endif