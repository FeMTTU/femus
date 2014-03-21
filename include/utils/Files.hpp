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
// =======================================
class Files {

protected:
  
  std::string                       _basepath;
  RunTimeMap<std::string>           _frtmap;   //this map cannot be declared as CONST because at some point it is FILLED and it is not in the initialization in the constructor!!!

public:
  
  std::ofstream _case_data;

  // Constructor-Destructor ------------------------
   Files(std::string = DEFAULT_BASEPATH);  ///< Constructor //TODO seems like it doesn't work with ONE DEFAULT PARAMETER
  ~Files(); ///< Destructor

  // Return functions ---------------------------------
  inline RunTimeMap<std::string>   get_frtmap() const {return  _frtmap;}    //TODO this is allowed!!!
  inline RunTimeMap<std::string> & get_frtmap()       {return  _frtmap;}         // I WANT THIS TO RETURN a REFERENCE, because this is going to call the READ FUNCTION which modifies the object 
  inline            std::string  get_basepath() const {return _basepath; }  //THIS RETURNS a COPY


  void PrintRun(const std::string run_name_in) const;   //for restart

  void ComposeOutdirName();
  void InitCaseData();
  void CloseCaseData();

  void CheckDir(const std::string& dir_name_in, const std::string& my_name_in) const;
  void CheckDirOrMake(const std::string& dir_name_in, const std::string& my_name_in) const;  //TODO isn't this very similar to the previous one?
  void CheckDirOrAbort(const std::string& dir_name_in, const std::string& my_name_in) const;

  void CheckIODirectories();
  void RedirectCout(std::streambuf* sbuf,  std::ofstream& file_in);

  void CopyFile(std::string  f_in,std::string  f_out) const;
 
  void CopyGencaseFiles();

};


#endif

