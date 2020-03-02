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

#include "FemusConfig.hpp" //for log_petsc
#ifdef HAVE_PETSC
#include "petsc.h"
#endif
#include "FemusDefault.hpp"
#include "paral.hpp"//to get iproc HAVE_MPI is inside here


namespace femus {


  //static data ===================================
  std::ofstream Files::file_sbuf;

  
  

  Files::Files()  { }

	
  Files::~Files() { }

  
  
  void Files::ConfigureRestart() {

//         if (paral::get_rank() == 0) { //QUESTA LETTURA LA POSSONO FARE TUTTI I PROCESSORI!
            std::cout << " Reading the  run_to_restart_from file to determine restart status or not" << std::endl;

    std::string app_path = "./";
    std::string lastrun_str;
    lastrun_str = app_path + DEFAULT_OUTPUTDIR + "/" + DEFAULT_LAST_RUN;

    //check if last_run is there, if it's not there go ahead and set restart = FALSE
            std::string lastone;
            std::ifstream last_run;
            last_run.open(lastrun_str.c_str());
            if (!last_run.is_open()) {
                std::cout << "There is no last_run file, it means that someone deleted it or it is the first run. No problem, go ahead. " << std::endl;
   // If the user wants to start he knows that he has to check that file so this procedure is SAFE
	    }

            last_run >> lastone >> lastone;  //"run_to_restart_from" is the STRING i use in the file, so there is one intermediate "buffer"...

            std::string restart_flag_string;
            last_run >> restart_flag_string >> restart_flag_string;
	    _restart_flag = atoi(restart_flag_string.c_str());   //TODO CONVERT STRING TO BOOLEAN
	    std::cout << "The restart flag is " << _restart_flag << std::endl;
	    
// 	}

	    if (_restart_flag) {
	      std::cout << "*** RESTART is activated *****" << std::endl; 
	      //we must set the basepath accordingly
	    
	    _input_path = app_path + DEFAULT_OUTPUTDIR + "/" + lastone + "/";
	    
	      std::cout << "*** The new input path is *****" << _input_path << std::endl; 
	    
	    }
	    else { std::cout << "Normal simulation without restart" << std::endl; 
	    //Notice that the basepath in this case was set by the ARGUMENT of the CONSTRUCTOR
	    
	    _input_path = app_path;
	    
	    }

	    
// maybe excess of safety... but it never hurts
#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif      
     
     return; 
    }

  
  

/// This function just checks if the directory is there
/// The first argument reads the absolute path of the parent directory,
/// the second argument reads the name (no absolute path) of the directory to be checked
///The sum of the first + the second is absolute, as it should be in order to do the check
// TODO this function works if the second argument is EMPTY, I want to BLOCK also this case!!!
void Files::CheckDirOrAbort(const std::string& dir_name_in, const std::string& my_name_in) {
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

void Files::CheckDirOrMake(const std::string& dir_name_in, const std::string& my_name_in) {
  
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
void Files::CheckDir(const std::string& dir_name_in, const std::string& my_name_in) {
  
  if (paral::get_rank() == 0 )   {
  
  // input directory ---------------------------------
  std::ostringstream dirname;
  dirname <<  dir_name_in << "/" << my_name_in ;

  std::ifstream in;  in.open(dirname.str().c_str()); 

  if (!in.is_open()) {
    std::cout << std::endl << " No " << dirname.str()  << " directory: I'll create it. "  << std::endl;

   int status = mkdir(dirname.str().c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
      if (status != 0) {std::cout << "MKDIR error: " << status << std::endl; abort();}

  }
  //else {std::cout << std::endl <<" That's alright, " << my_name_in << " is already there." << std::endl;}

in.close();

  } //end proc==0
  
//after proc==0 put a barrier, so that all the threads will have that directory from this point on
#ifdef HAVE_MPI
MPI_Barrier(MPI_COMM_WORLD);
#endif

  
return;
}




  //the goal of this routine is to print a file with the current OUTTIMEDIR in it
  //it should be called "print outtime dir"
    //here you have to pass the RELATIVE_NAME of the file to be printed, 
  //it composes it wrt the OUTPUTDIR
  //call this function when the OUTPUT DIR is already created
  //and also after OUTTIMEDIR is composed
  
   //print this string to /base_output/new_run
  //Files does not have Utils yet

void Files::PrintRunForRestart(const std::string run_name_in) const {


   if (paral::get_rank() == 0 )   {

   std::string app_path = "./";
   std::string run("");
   run = app_path + DEFAULT_OUTPUTDIR + "/" + run_name_in; //AAA BASE OUTPUT
   std::cout << "Print the run " << run << "to file" << std::endl;

   std::ofstream run_file; run_file.open(run.c_str());
   run_file << DEFAULT_LAST_RUN << " ";
   run_file << _output_time;
  
      run_file << std::endl;
      run_file << "flag_for_restart " << 0;  //TODO if you printed 1 here it would turn out to be a NONSTOPPING CHAIN OF SIMULATIONS!!!!
      
   
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

//****** Proc0 generates the string    
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

//****** Proc0 broadcasts the string size and the string content to every proc    
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

 _output_time = new_out;
 
 std::cout << "iproc = " << paral::get_rank() << " ***** The output dir of this run will be: " << DEFAULT_OUTPUTDIR << "/" << _output_time << std::endl;

 //************************
 //set the input and output_path variables
    std::string app_path = "./";
   _output_path = app_path + DEFAULT_OUTPUTDIR + "/" + _output_time + "/";
 
 
 
 
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


/// ========== Copy Input Files ==========
//it copies the mesh files to the outtime dir
// TODO this function is performed by ALL PROCESSORS but we should do it only for the FIRST ONE...
  void Files::CopyInputFiles() const { 

 std::cout << "TODO: MUST FIND A WAY TO COPY A WHOLE DIRECTORY AND NOT THE SINGLE FILES" << std::endl;
    
CheckDirOrMake(_output_path,DEFAULT_INPUTDIR);

//copy configuration file
   std::string op_in  =   _input_path + "/" + DEFAULT_INPUTDIR + "/" + DEFAULT_RUNTIMECONF;
   std::string op_out =  _output_path + "/" + DEFAULT_INPUTDIR + "/" + DEFAULT_RUNTIMECONF;
/*(iproc==0)*/ CopyFile(op_in,op_out);

//TODO here we should also copy the mesh file from the mesh generator... but we need to know the filename...

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
// OUTPUTDIR: set the run_to_restart_from (could rewrite _utilsmap.in...)

//All the processors can do this check
// I think it is better to program with as few "if proc==0" as possible.
//if you let a check be done by only one processor, then all the others must wait for that check
//so a thing that is useful to everyone but with "if proc==0" must be followed by a MPI_Barrier
 
//============== OUTPUT
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

//==============   
// CheckDirOrMake is better a makedir! because it never exists
//if you let it do by all processors, then the fastest one will create the file, the others will find it
//AAA there can be one case when the directory already exists! Suppose there is a script running in background 
//and in the meantime i do some runs by hand... then i may start a run in the very precise moment when 
//the scripts starts a run... very unlikely but TODO check that
//if in that very moment the background script prints that dir, then it will be already there for ALL PROCESSORS,
//even for the fastest one, so there will be two CONFLICTING RUNS in the SAME FOLDER!
//so the OUTTIME DIR can be already there for two reasons: either for a background script of another run,
// or for the fastest process of my runs. That is why it is better to let this check be done only by processor zero.
//If the dir is already there, it is because of a contemporary background run
 
void Files::CheckIODirectories() {
 
//INPUT
                    std::string abs_app = "./";
/*all procs*/   CheckDirOrAbort(abs_app,DEFAULT_INPUTDIR); //it must be there only to be COPIED (and we don't even need the check in restart case)

/*all procs*/   CheckDir(abs_app,DEFAULT_OUTPUTDIR);

   std::string abs_outputdir = abs_app + "/" + DEFAULT_OUTPUTDIR;
/*(iproc==0)*/  ComposeOutdirName();  //this adds an element to the Files map, so the SetupAll function cannot be constant

/*(iproc==0)*/  CheckDirOrMake(abs_outputdir,_output_time);

// at this point we should copy the input files in the outtime directory

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

  void Files::RedirectCout() const {
//this redirects the standard output to file, only for processor 0
//I dont like the fact that all these functions must be modified
//whenever I change the absolute path,
//I'd prefer passing the file name explicitly

    std::string app_path = "./";
    std::string abs_runlog = app_path + DEFAULT_OUTPUTDIR 
    + "/" + _output_time +  "/" + DEFAULT_RUN_LOG;

//  std::ofstream file;  //if a filestream dies, then also its stream-buffer dies ?!? 
//                       //So I have to declare it outside? Yes. This seems to work.
//    std::streambuf* sbuf = std::cout.rdbuf();  //get the current buffer for cout

    

//////// MULTIPRINT // // if (paral::get_rank() == 0) {  //only processor 0 has to open the file stream
 
    std::stringstream runlogproc;
    runlogproc << abs_runlog << "_p" << paral::get_rank() << DEFAULT_EXT_LOG; //multiprint
    
    abs_runlog = runlogproc.str();
    
     Files::file_sbuf.open(abs_runlog.c_str());
    std::streambuf* fbuf = Files::file_sbuf.rdbuf();  //get the current buffer for the file
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
 
     
    
   return; 
  }


 //it seems like you have to give the stream buffer back to cout !!!
 // http://wordaligned.org/articles/cpp-streambufs
 //well, now I tried to comment this function and EVERYTHING WORKS FINE... so maybe now it is not needed anymore...
 void Files::RedirectCoutFinalize(std::streambuf* sbuf) {
      
     std::cout.rdbuf(sbuf);  //this is a SET function that acts on sbuf

   return; 
     
  }
  
  
  
//==========================================
  void Files::log_petsc() const {
    
  std::string petsc_femus_log = "petsc.log";
  std::ostringstream petsc_log;
  petsc_log << _output_path << "/" << petsc_femus_log;

#ifdef HAVE_PETSC
   PetscViewer my_viewer;
   PetscViewerCreate(MPI_COMM_WORLD, &my_viewer);
   PetscViewerSetType(my_viewer, PETSCVIEWERASCII);
   PetscViewerFileSetName(my_viewer, petsc_log.str().c_str());   
   PetscLogView(my_viewer);
#endif
   
   return;
  }


} //end namespace femus


