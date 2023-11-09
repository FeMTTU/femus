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



#include <string>
#include <fstream>
#include <sstream>


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

  
// === Constructors / Destructor  - BEGIN =================
public:
  
   Files();
   
  ~Files();
// === Constructors / Destructor  - END =================

// Directory management - BEGIN =========
public:
  
  void CheckIODirectories(const bool use_output_time_folder);

  static void CheckDir(const std::string& dir_name_in, const std::string& my_name_in);
  
private:
  
         void ComposeOutdirName(const bool use_output_time_folder);
  static void CheckDirOrMake(const std::string& dir_name_in, const std::string& my_name_in);
//   static void CheckDir(const std::string& dir_name_in, const std::string& my_name_in);
  static void CheckDirOrAbort(const std::string& dir_name_in, const std::string& my_name_in);
// Directory management - END =========

  


// Input - BEGIN ======================
public:
  
  static const std::string runtime_config_filename() { return "femus_conf.in"; };  //how to have it in the header!
  
  static const std::string  _application_input_directory;
  
  std::string  GetInputPath() const {
    return _input_path;
  }
  
  
static  std::string get_input_file_with_prefix(const std::string input_file, const std::string relative_location_of_input_folder)  {

      std::ostringstream mystream; mystream << relative_location_of_input_folder  /*"./"*/ << Files::_application_input_directory << "/" << input_file;
      const std::string infile = mystream.str();

      return infile;

   }
   
  private:
    
  std::string  _input_path; //this is where the input files are located BEFORE YOU COPY THEM to the OUTTIME DIR!!!!! it has to alternatives in case of restart or not
// Input - END ======================


// Input, Copy - BEGIN ========================= 
 public:
  
  void CopyInputFiles() const;
  
 private:

  void CopyFile(std::string  f_in,std::string  f_out) const;
// Input, Copy - END ========================= 

  
// Output - BEGIN ======================
public:
  
  static const std::string  _application_output_directory;
  
  std::string  GetOutputPath() const {
    return _output_path;
  }
  
  std::string  GetOutputTime() const {
    return _output_time;
  }
  
  
 private:
   
  std::string _output_path;
  std::string _output_time; //this is the OUTTIME DIR!!!
 
// Output - END ======================


// Output, LOG - BEGIN ==========================
 public:

// Stream redirect to file ======
         void RedirectCout(const bool redirect_cout_to_file) const;
  static void RedirectCoutFinalize(std::streambuf* sbuf);
         void log_petsc() const;

 private:
   
  static const std::string _run_log_basename;
  static const std::string _run_log_extension;
         
// Output, LOG - END ==========================


// Restart - BEGIN ======================
 public:
   
 // # if a run reaches the end, then we write it as a "default restart" run
  static const std::string run_to_restart_from_string() { return "run_to_restart_from"; };  //how to have it in the header!

  void ConfigureRestart();
  void PrintRunForRestart(const std::string run_name_in) const;

  bool  GetRestartFlag() const {
    return _restart_flag;
  }
  
 private:

  bool _restart_flag; 
// Restart - END ======================

  
// Mesh files database - BEGIN ======================
 public:
   
  static const std::string mesh_folder_path() {  return "src/06_mesh/00_single_level/01_input/00_mesh_files/";  }

// Mesh files database - END ======================
  
private:
  
  static std::ofstream file_sbuf;  //needed for I/O purposes
  

};


/// @todo  FORWARD SLASHES and BACKSLASHES //perhaps one day we should consider to treat these slashes appropriately, in case we compile on Windows etc.




} //end namespace femus



#endif
