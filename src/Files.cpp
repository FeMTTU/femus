// c/c++
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>  //for mkdir
#include <ctime>        //for output dir name
#include <cstring>      //for strcpy

// class
#include "Files.hpp"

// configuration includes
#include "FemusDefault.hpp"
#include "paral.hpp"//to get iproc HAVE_MPI is inside here



// =================================================
//         Files Class functions
// =================================================

  Files::Files(std::string  string_in) : 
        _basepath(string_in),
        _frtmap("Files",_basepath)   { 
	  
	    if (_basepath == "")  { std::cout << " Set the basepath in the command line" << std::endl;    abort(); }

	}  ///< Constructor

	
  Files::~Files() { } ///< Destructor


//=============================================
void Files::InitCaseData() {
  
  //construct the common output file stream for the data
  //first,build its name
   std::ostringstream data_name;
   data_name << get_basepath() << "/" << get_frtmap().get("OUTPUT_DIR") << get_frtmap().get("OUTTIME_DIR") << get_frtmap().get("CASE_DATA");

   //open the common output file stream
     if (paral::get_rank() == 0 )   {
   _case_data.open(data_name.str().c_str());  //I'll close it in the Files destructor

     }
     //Is this open in parallel?
   //It seems like it's working even without proc==0. I'll put it to be safe.
  
  
  
  return;
}

void Files::CloseCaseData() { 
  _case_data.close();
   return;
  }


/// This function just checks if the directory is there
/// The first argument reads the absolute path of the parent directory,
/// the second argument reads the name (no absolute path) of the directory to be checked
///The sum of the first + the second is absolute, as it should be in order to do the check
void Files::CheckDirOrAbort(const std::string& dir_name_in, const std::string& my_name_in) const {
//can be done by all processes
  
  // input directory ---------------------------------
  std::ostringstream dirname;
  dirname << dir_name_in << "/" << my_name_in ;

  std::ifstream in;  in.open(dirname.str().c_str()); 

  if (!in.is_open()) {
    std::cout << " No " << dirname.str()  << " directory: Abort "  << std::endl;
    abort();
  }  
  
  
  return;
}

///This function checks one directory
///the first argument is th absoltue path of the paretn directory
///the second argument is the name of the dierctory
///This function gives error if the directory is already there,
///but for the gencase I dont have any problems
///So I'll make a function that creates the directory if it isnt there,
///but does not complain in the other case

  //this check_dir is designed to work onnly with ABSOULUTE PATHS
  //you pass to this function the PARENT STRING and YOUR DESIRED NAME you make the new directory
  
  //---WHEN YOU CALL THIS FUNCTION, YOU HAVE TO BE SURE THAT THE PARENT DIR ALREADY EXISTS! put a check
  //---ALSO, YOU HAVE TO PASS THE ABSOLUTE PATH of the PARENT DIR, and the RELATIVE PATH of the NEW_DIR
  
//for what it does, this function must be run by only one processor
// So either you call it with an external if, or internally.
// Since external people may not know that, you must put it INTERNALLY
//on the other hand, if only one processor does that, the other processes must be synchronized,
//otherwise the check may not be 

void Files::CheckDirOrMake(const std::string& dir_name_in, const std::string& my_name_in) const {
  
  if (paral::get_rank() == 0 )   {
  
  // input directory ---------------------------------
  std::ostringstream abs_dirname;
  abs_dirname << dir_name_in << "/" << my_name_in ;

  std::ifstream in;  in.open(abs_dirname.str().c_str()); 

  if (!in.is_open()) {
    std::cout << " No " << abs_dirname.str()  << " directory: I'll create it. "  << std::endl;

    //AAA if you do mkdir A/B/C/D, A,B and C must already exist!
//in shell language you do mkdir -p A/B/C/D
//how do you do that in C?
//I should do it my own, or see the contents of the shell mkdir source file.
//clearly, the C version of typical shell functions are not as complete as those.
//- I could do mkdir and cd inside here... getting output/ and the date separately.
//-I could try mkdir at,but what  file descriptor do I have to pass?  mkdirat(1,"gigio",0777); exit(30);
//Or, I can try with boost library

//For now,I'll just do the brute force way

//   std::ostringstream base_output;
//   base_output << _femus_dir << "/" << "output" ;
//   int status2 = mkdir(base_output.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
//      if (status2 != 0) { std::cout << "No problem for now" << std::endl;}

   int status = mkdir(abs_dirname.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
      if (status != 0) {std::cout << "MKDIR error: " << status << std::endl;abort();}

//       status = mkdir("cassa/",S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
//instead of checking if the file is there, I can do it directly with mkdir
//but this is not the shell function,so it wont result in an abort of my program unless you explicitly put it.
//this is one reason why a function that was born for the shell is not so well done in C

  }
  else {std::cout<< "Some conflicting run at exactly the same time, unlikely but check it"<<std::endl;abort();}

in.close();

  } //end proc==0
  
//after proc==0 put a barrier, so that all the threads will have that directory from this point on
#ifdef HAVE_MPI
MPI_Barrier(MPI_COMM_WORLD);
#endif

  
return;
}


///if the directory isnt there, create it
///if the directory is there, nothing wrong
void Files::CheckDir(const std::string& dir_name_in, const std::string& my_name_in) const {
  
  if (paral::get_rank() == 0 )   {
  
  // input directory ---------------------------------
  std::ostringstream dirname;
  dirname <<  dir_name_in << "/" << my_name_in ;

  std::ifstream in;  in.open(dirname.str().c_str()); 

  if (!in.is_open()) {
    std::cout << " No " << dirname.str()  << " directory: I'll create it. "  << std::endl;

   int status = mkdir(dirname.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
      if (status != 0) {std::cout << "MKDIR error: " << status << std::endl; abort();}

  }
  else {std::cout<< "That's alright, " << my_name_in << " is already there." << std::endl;}

in.close();

  } //end proc==0
  
//after proc==0 put a barrier, so that all the threads will have that directory from this point on
#ifdef HAVE_MPI
MPI_Barrier(MPI_COMM_WORLD);
#endif

  
return;
}




  //the goal of this routine is to print a file with the current OUTTIME_DIR in it
  //it should be called "print outtime dir"
    //here you have to pass the RELATIVE_NAME of the file to be printed, 
  //it composes it wrt the OUTPUT_DIR
  //call this function when the OUTPUT DIR is already created
  //and also after OUTTIME_DIR is composed
  
   //print this string to /base_output/new_run
  //Files does not have Utils yet

void Files::PrintRun(const std::string run_name_in) const {


   if (paral::get_rank() == 0 )   {

   std::string run("");
   run = get_basepath() + "/" + get_frtmap().get("OUTPUT_DIR") + "/" + run_name_in; //AAA BASE OUTPUT
   std::cout << "Print the run " << run << "to file" << std::endl;

   std::ofstream run_file; run_file.open(run.c_str());
   run_file << get_frtmap().get("OUTTIME_DIR");
  
   run_file.close();
  
   }
  
  return;
}

///This function composes the outdir name
//the format i want is sthg like out_YMD_hms
//to avoid any possible lack of synchronism,
//only processor 0 will compose this name
//and broadcast it to all the other processors
//how do I broadcast a string?
//but then give it to all processors, i.e. to all instantiations of the class Files
//then the others will wait until this name has been composed

//The good thing to do is to SYNC all processor CLOCKS at the beginning, to be sure 
// that they get exactly the same time
//can this be done in OpenMPI?
//MPI_WTIME_IS_GLOBAL may not be changed... it is set to 3. see man MPI_Comm_set_attr
//in OpenMPI the clocks are not guaranteed to be synchronized.
//so, I cannot expect that

////////////    
/////////// Test
//  to show that the processors may have different times
// if (paral::get_rank() == 0) {system("sleep 5");}
// and uncomment the if rank below
/////////// END test
////////////

	//set the new output dir in the map
//only processor 0 will compose this name
//then it will broadcast it to all the others 

 /* Write formatted members of sys_time into the sysdate_buffer. */
//   #define NCHAR_SYSDATE_BUFFER 256
//   char sysdate_buffer[NCHAR_SYSDATE_BUFFER];
//   strftime( sysdate_buffer, NCHAR_SYSDATE_BUFFER, "%a %b %d %H:%M:%S %Z %Y", sys_time ); 
//THIS MUST BE PARALLEL ALSO IF YOU LIKE
//   std::cout << "***** Nicely formatted time: " << sysdate_buffer << std::endl;


//STRING BROADCAST TO ALL THE OTHER PROCESSORS WITH BOOST... TODO
// boost::mpi::communicator world;  //how can I be sure that this is the same communicator as MPI_COMM_WORLD?
// boost::mpi::broadcast(world,new_out,0);

//There is something I dont understand.  I allocate a vector of chars of a certain length,
//smaller than it should be, but when i do std::cout << it prints ALL OF IT,
//even what is OUTSIDE of the VECTOR OF CHARS!

// CHAR BROADCASTING
//we have to convert the string into char*
//we have to know the length of this array
//we only know it after we compose the vector

//so first I have to broadcast the length of the array from processor 0 
//then I can do the allocation
//if proc==0 i fill it with the string again only in processor 0
//then I broadcast the content of the array
//then the string new_out can be filled with the char*
// add a ne to the map

void Files::ComposeOutdirName() {
  
   struct tm* sys_time = NULL;
   std::string new_out;
   char* out_char;
   int outchar_size;
   std::ostringstream outname;

    
   if (paral::get_rank() == 0) {
  
  time_t      curtime = 0; 
 
     curtime = time(NULL);                         //C function
//    time_t curtime_mpi = (time_t)  MPI_Wtime();  //MPI Function for time
   
    sys_time = localtime(&curtime);  //converts time into broken-down representation

  const uint year   = 1900 + sys_time->tm_year;
  const uint month  = 1 +sys_time->tm_mon;
  const uint day    = sys_time->tm_mday;
  const uint hour   = sys_time->tm_hour;
  const uint minute = sys_time->tm_min;
  const uint second = sys_time->tm_sec;

   
        outname  << year
        << "-" << std::setw(2) << std::setfill('0') << month
        << "-" << std::setw(2) << std::setfill('0') << day << "_"
               << std::setw(2) << std::setfill('0') << hour
        << "-" << std::setw(2) << std::setfill('0') << minute
        << "-" << std::setw(2) << std::setfill('0') << second << "/";

  new_out = outname.str(); // declare it outside the block as all the processes will see it

  outchar_size = outname.str().size();
  
    }

#ifdef HAVE_MPI
  MPI_Bcast(&outchar_size,1,MPI_INT,0,MPI_COMM_WORLD);
#endif

out_char = new char[outchar_size +1]; //for the null

  if (paral::get_rank() == 0)   std::strcpy(out_char,new_out.c_str());


#ifdef HAVE_MPI
MPI_Bcast(out_char,outchar_size,MPI_CHAR,0,MPI_COMM_WORLD);
#endif

 out_char[outchar_size] = NULL; //HEY MAYBE NOT 100% ORTHODOX BUT OTHERWISE IT DOESN'T WORK!!!

// printf("^^^^^^^^^^^^^^^^^^^ %s\n",out_char);
//ok at this point the error already comes out: one processor has one more char. when do they stop printing?
//i guess when they find a null... so let us add a null where it is supposed to be!

  // new_out.resize(outchar_size);
  new_out = /*outname.str().c_str()*/out_char;  //uguagliare una stringa ad un char*

 delete [] out_char;

 this->get_frtmap().set("OUTTIME_DIR",new_out);

 std::cout << "iproc = " << paral::get_rank() << " ***** The output dir of this run will be: " << get_frtmap().get("OUTPUT_DIR") << get_frtmap().get("OUTTIME_DIR") << std::endl;

 return; 

}


/// This function copies a file from an input to an output names
//it works by passing the ABSOLUTE PATHS of the source and destination files

void Files::CopyFile(std::string  f_in,std::string  f_out) const {
  
  if (paral::get_rank() == 0) {
  
  std::ifstream f1(f_in.c_str(), std::fstream::binary);
  std::ofstream f2(f_out.c_str(), std::fstream::trunc| std::fstream::binary);
  f2 << f1.rdbuf();
    
  }
  
  return;
}

/// ========== Copy Gencase Files ==========
//it copies the mesh files to the outtime dir
//the multigrid files are not copied here
    /*TODO fileIO*/
  void Files::CopyGencaseFiles() { 

   std::string app_basepath = get_basepath() + "/"; 

    // >>>>>>> outtime dir: COPY FILES   //needs the BASEPATH of the APPLICATION
  //copy mesh.h5
   std::string  mesh_in = app_basepath + "/" + get_frtmap().get("INPUT_DIR") + "/" + get_frtmap().get("BASEMESH") + get_frtmap().get("EXT_H5");
   std::string  mesh_out =  app_basepath + "/" + get_frtmap().get("OUTPUT_DIR") + "/" + get_frtmap().get("OUTTIME_DIR") +  "/" + get_frtmap().get("BASEMESH") + get_frtmap().get("EXT_H5");
/*(iproc==0)*/ CopyFile(mesh_in,mesh_out); 
   
//copy multimesh.xmf
   std::string  mmesh_in = app_basepath + "/" + get_frtmap().get("INPUT_DIR") + "/" + get_frtmap().get("MULTIMESH") + get_frtmap().get("EXT_XDMF");
   std::string  mmesh_out = app_basepath + "/" + get_frtmap().get("OUTPUT_DIR") + "/" + get_frtmap().get("OUTTIME_DIR") +  "/" +  get_frtmap().get("MULTIMESH") + get_frtmap().get("EXT_XDMF");
/*(iproc==0)*/ CopyFile(mmesh_in,mmesh_out); 

//copy MG files
   std::string  op_in = app_basepath + "/" + get_frtmap().get("INPUT_DIR") + "/" + get_frtmap().get("F_MATRIX") + get_frtmap().get("EXT_H5");
   std::string  op_out = app_basepath + "/" + get_frtmap().get("OUTPUT_DIR") + "/" + get_frtmap().get("OUTTIME_DIR") +  "/" +  get_frtmap().get("F_MATRIX") + get_frtmap().get("EXT_H5");
/*(iproc==0)*/ CopyFile(op_in,op_out);

   op_in  = app_basepath + "/" + get_frtmap().get("INPUT_DIR") + "/" + get_frtmap().get("F_REST") + get_frtmap().get("EXT_H5");
   op_out = app_basepath + "/" + get_frtmap().get("OUTPUT_DIR") + "/" + get_frtmap().get("OUTTIME_DIR") +  "/" +  get_frtmap().get("F_REST") + get_frtmap().get("EXT_H5");
/*(iproc==0)*/ CopyFile(op_in,op_out);

   op_in  = app_basepath + "/" + get_frtmap().get("INPUT_DIR") + "/" + get_frtmap().get("F_PROL") + get_frtmap().get("EXT_H5");
   op_out = app_basepath + "/" + get_frtmap().get("OUTPUT_DIR") + "/" + get_frtmap().get("OUTTIME_DIR") +  "/" +  get_frtmap().get("F_PROL") + get_frtmap().get("EXT_H5");
/*(iproc==0)*/ CopyFile(op_in,op_out);

//TODO we should also copy the files in the config/ directory, so we keep track of ALL the input!!!

//barrier so that all the processors will have mesh.h5 and else to read from
//if you copy ALL THE I/O you need into the run, then the restart will read ONLY FROM THAT RUN.
//In fact the restart means RESTARTING FROM A CERTAIN RUN.
//This restricts the possibility of changing processors/levels...
#ifdef HAVE_MPI
MPI_Barrier(MPI_COMM_WORLD);
#endif

    return;
  }




/// This function checks if the needed directories are where they have to
/// This is called by an APP, not by the gencase
/// Beh
void Files::CheckIODirectories() {
  //this does the preparation for a single run
 
 //after reading the file names, we can SET UP the NEW RUN

//the first processor will do it, the others will find it
//these checks are not: check or create, but check or stop
//if these directories arent there you cant do much. I mean, they could also be there
//but be empty, but at least check if they are there.
//These directories have to be THERE and FULL, PREPARED BEFORE BY THE OTHER "gen" programs and by the user

//So its up to the user to PREPARE the INPUT DIRS:
// 'config_dir': set param*.in
               //set the right file names in param_files.in
// INPUT_DIR:  run the gencase
// OUTPUT_DIR: set the run_to_restart_from (could rewrite _utilsmap.in...)

//All the processors can do this check
// I think it is better to program with as few "if proc==0" as possible.
//if you let a check be done by only one processor, then all the others must wait for that check
//so a thing that is useful to everyone but with "if proc==0" must be followed by a MPI_Barrier
 
 
 
 //========= CHECK  NEEDED I/O DIRECTORIES
//INPUT: fem/,config/,input/
                    std::string abs_app = get_basepath() + "/";
/*all procs*/   CheckDirOrAbort(abs_app,DEFAULT_CONFIGDIR);
/*all procs*/   CheckDirOrAbort(abs_app,get_frtmap().get("INPUT_DIR"));


//OUTPUT: output/
/*all procs*/   CheckDir(abs_app,get_frtmap().get("OUTPUT_DIR"));

  //if you reached this point, then output/ exists, so you can do the makedir for outtime!
   //the executable has to make a dir which is  in its ./ directory...
   //we could have made that in a run script, doing those operations with the shell.
   //for instance, you may not want a program (or also a script) to write somewhere else
   //other than its own directory, because you dont know if you CAN write on those other paths
   //but you know you can write on the place where you launched the executable
   //that's why it would be nicer to give just RELATIVE PATHS. But on the other hand,
   //if all the absolute paths you give are within 'femus_dir', it means that you are always 
   //inside your "internal zone".
   //'femus_dir' is a sort of "fence" that tells you: stay inside here.
   //it's better to give relative paths also because the TEXT .xmd files that you write 
   //will only contain RELATIVE PATHS, so they will be readable anywhere


   std::string abs_outputdir = abs_app + "/" + get_frtmap().get("OUTPUT_DIR");
/*(iproc==0)*/  ComposeOutdirName();  //this adds an element to the Files map, so the SetupAll function cannot be constant
/*(iproc==0)*/  CheckDirOrMake(abs_outputdir, get_frtmap().get("OUTTIME_DIR"));

//this is better a makedir! because it never exists
//if you let it do by all processors, then the fastest one will create the file, the others will find it
//AAA there can be one case when the directory already exists! Suppose there is a script running in background 
//and in the meantime i do some runs by hand... then i may start a run in the very precise moment when 
//the scripts starts a run... very unlikely but TODO check that
//if in that very moment the background script prints that dir, then it will be already there for ALL PROCESSORS,
//even for the fastest one, so there will be two CONFLICTING RUNS in the SAME FOLDER!
//so the OUTTIME DIR can be already there for two reasons: either for a background script of another run,
// or for the fastest process of my runs. That is why it is better to let this check be done only by processor zero.
//If the dir is already there, it is because of a contemporary background run

//========= END CHECK  NEEDED I/O DIRECTORIES



 return; 
}


//==========================================
//Here, i can put some alternatives.
//either print on a single file
//or print on terminal with single file
//or print with multiple files

//if you want to DUPLICATE the output both to file and to screen
// you have to do some TEE streams, like for shell... see:
//http://wordaligned.org/articles/cpp-streambufs
//or you could do everywhere "myout << std::cout << ..."
//but, while cout is a GLOBAL variable, it only belongs to the std namespace,
//i should define another variable "myout" which should be 
//global as well, or better belonging to the FEMuS:: namespace.
//in this way i do not have to associate myout to a class, like utils,
//and give it to all the other classes as an instantiation in the constructor.
//A global variable is still an instantiation but not associated to any class
//to be passed to the constructor

  void Files::RedirectCout(std::streambuf* /*sbuf*/,  std::ofstream& file_in) {
//this redirects the standard output to file, only for processor 0
//I dont like the fact that all these functions must be modified
//whenever I change the absolute path,
//I'd prefer passing the file name explicitly

    std::string abs_runlog = get_basepath() + "/" + get_frtmap().get("OUTPUT_DIR") 
    + "/" + get_frtmap().get("OUTTIME_DIR") +  "/" + get_frtmap().get("RUN_LOG");

//  std::ofstream file;  //if a filestream dies, then also its stream-buffer dies ?!? 
//                       //So I have to declare it outside? Yes. This seems to work.
//    std::streambuf* sbuf = std::cout.rdbuf();  //get the current buffer for cout

    

//////// MULTIPRINT // // if (paral::get_rank() == 0) {  //only processor 0 has to open the file stream
 
    std::stringstream runlogproc;
    runlogproc << abs_runlog << "_p" << paral::get_rank() << get_frtmap().get("EXT_LOG"); //multiprint
    
    abs_runlog = runlogproc.str();
    
     file_in.open(abs_runlog.c_str());
    std::streambuf* fbuf = file_in.rdbuf();  //get the current buffer for the file
      std::cout << "Printing to file run_log ..." << std::endl;
      std::cout.rdbuf(fbuf);   //set the buffer of cout to the buffer of the file 
// 
// inside here there is a problem, it gives an error "below main"
//for now,I'll just redirect 

////////// MULTIPRINT  // //   }
//////////  MULTIPRINT // //   else { std::cout.rdbuf(NULL);}

// std::cout.rdbuf() must be 
// NULL for non-zero processors;
// the file for processor 0

std::cout << "The number of processors is " << paral::get_size() << std::endl;
std::cout << "The mode in which this executable was compiled is " << getenv("FM_FEMUS_METHOD") << std::endl;
 
     
    
   return; 
  }
