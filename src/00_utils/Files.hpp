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

#ifndef __femus_utils_Files_hpp__
#define __femus_utils_Files_hpp__

#include "FemusDefault.hpp"
#include "FemusInputParser.hpp"

// C++
#include <map>
#include <string>
#include <fstream>


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
         void CheckIODirectories(const bool use_output_time_folder);
  
// Copy ========================= 
  void CopyInputFiles() const;

// Restart ======================
  void ConfigureRestart();
  void PrintRunForRestart(const std::string run_name_in) const;
  
// LOG ==========================
// Stream redirect to file ======
         void RedirectCout(const bool redirect_cout_to_file) const;
  static void RedirectCoutFinalize(std::streambuf* sbuf);
         void log_petsc() const;

// get=============
  std::string  GetOutputPath() const {
    return _output_path;
  }
  
  std::string  GetInputPath() const {
    return _input_path;
  }

static  std::string get_input_file_with_prefix(const std::string input_file, const std::string relative_location_of_input_folder)  {


      std::ostringstream mystream; mystream << relative_location_of_input_folder  /*"./"*/ << DEFAULT_INPUTDIR << "/" << input_file;
      const std::string infile = mystream.str();

      return infile;

   }

  std::string  GetOutputTime() const {
    return _output_time;
  }
  
  bool  GetRestartFlag() const {
    return _restart_flag;
  }
  
  static void CheckDir(const std::string& dir_name_in, const std::string& my_name_in);
  
private:
  
  static std::ofstream file_sbuf;  //needed for I/O purposes
  std::string  _input_path; //this is where the input files are located BEFORE YOU COPY THEM to the OUTTIME DIR!!!!! it has to alternatives in case of restart or not
  std::string _output_path;
  std::string _output_time; //this is the OUTTIME DIR!!!

  bool _restart_flag; 
  
  
// Directory management
         void ComposeOutdirName(const bool use_output_time_folder);
  static void CheckDirOrMake(const std::string& dir_name_in, const std::string& my_name_in);
//   static void CheckDir(const std::string& dir_name_in, const std::string& my_name_in);
  static void CheckDirOrAbort(const std::string& dir_name_in, const std::string& my_name_in);

// Copy 
  void CopyFile(std::string  f_in,std::string  f_out) const;
  
};



} //end namespace femus



#endif
