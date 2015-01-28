/*=========================================================================

 Program: FEMUS
 Module: Files
 Authors: Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __files_hpp__
#define __files_hpp__

// C++
#include <map>
#include <string>
#include <fstream>

#include "FemusDefault.hpp"
#include "FemusInputParser.hpp"


namespace femus {



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
  
   Files();
  ~Files();

// Directory management =========
         void CheckIODirectories();
  
// Copy ========================= 
  void CopyInputFiles() const;

// Restart ======================
  void ConfigureRestart();
  void PrintRunForRestart(const std::string run_name_in) const;
  
// LOG ==========================
// Stream redirect to file ======
         void RedirectCout() const;
  static void RedirectCoutFinalize(std::streambuf* sbuf);
         void log_petsc() const;

// get=============
  std::string  GetOutputPath() const {
    return _output_path;
  }
  
  std::string  GetInputPath() const {
    return _input_path;
  }
  
  std::string  GetOutputTime() const {
    return _output_time;
  }
  
  bool  GetRestartFlag() const {
    return _restart_flag;
  }
 
private:
  
  static std::ofstream file_sbuf;  //needed for I/O purposes
  std::string  _input_path; //this is where the input files are located BEFORE YOU COPY THEM to the OUTTIME DIR!!!!! it has to alternatives in case of restart or not
  std::string _output_path;
  std::string _output_time; //this is the OUTTIME DIR!!!

  bool _restart_flag; 
  
  
// Directory management
         void ComposeOutdirName();
  static void CheckDirOrMake(const std::string& dir_name_in, const std::string& my_name_in);
  static void CheckDir(const std::string& dir_name_in, const std::string& my_name_in);
  static void CheckDirOrAbort(const std::string& dir_name_in, const std::string& my_name_in);

// Copy 
  void CopyFile(std::string  f_in,std::string  f_out) const;
  
};



} //end namespace femus



#endif