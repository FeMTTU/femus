 
//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

#include "Writer_one_level.hpp"
// #include "VTKWriter_one_level.hpp"
// #include "GMVWriter_one_level.hpp"
// #include "XDMFWriter_one_level.hpp"

#include "Solution.hpp"

#include "ElemType.hpp"
#include "FElemTypeEnum_list.hpp"

#include <cstring>


namespace femus {



  
    const std::string Writer_one_level::_name_bdc = "Bdc";
    const std::string Writer_one_level::_name_res = "Res";
    const std::string Writer_one_level::_name_eps = "Eps";

  const unsigned Writer_one_level::FemusToVTKorToXDMFConn[27] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 23, 21, 20, 22, 24, 25, 26};



    Writer_one_level::Writer_one_level() {
        
      initialize_flags();
      
    }
    
    
//     /** Constructor. */
//     Writer_one_level(const Solution * ml_sol) {
//         
//         
//     }



  
  void Writer_one_level::initialize_flags() {
    
    _moving_mesh = false;
    _graph = false;
    _surface = false;
    
    _debugOutput = false;
  }
  



  void Writer_one_level::SetMovingMesh (std::vector<std::string>& movvars_in) {
    _moving_mesh = true;
    _moving_vars = movvars_in;
  }

  void Writer_one_level::SetGraphVariable (const std::string &graphVaraible) {
    _graph = true;
    _surface = false;
    _graphVariable = graphVaraible;
  }

  void Writer_one_level::SetSurfaceVariables (std::vector < std::string > &surfaceVariable) {
    _surface = true;
    _graph = false;
    _surfaceVariables = surfaceVariable;
  }



    
   std::string Writer_one_level::print_sol_bdc_res_eps_name(const std::string solName, const unsigned name) const {
       
            std::string printName;

            if( name == _index_sol )      printName = solName;
            else if( name == _index_bdc ) printName = solName + "_" + _name_bdc;
            else if( name == _index_res ) printName = solName + "_" + _name_res;
            else if( name == _index_eps ) printName = solName + "_" + _name_eps;
            else { abort(); }

       return printName;     
   }
   
   
  
   unsigned Writer_one_level::compute_sol_bdc_res_eps_size(const Solution * solution, const unsigned i) const {
       
       const bool is_unknown = solution->is_unknown_of_system(i);

       const bool is_debug   = _debugOutput;       
       
       const unsigned print_sol_size = 1 + 3 * is_debug * is_unknown;
       return  print_sol_size;
       
   }
   
   
    unsigned Writer_one_level::fe_index(const std::string & order_str) const {
        
        unsigned index = 0;
        
    if( !strcmp( order_str.c_str(), fe_fams_for_files[ FILES_CONTINUOUS_LINEAR ].c_str() ) )           {  index = FILES_CONTINUOUS_LINEAR;  }
    else if( !strcmp( order_str.c_str(), fe_fams_for_files[ FILES_CONTINUOUS_QUADRATIC ].c_str() ) )   {  index = FILES_CONTINUOUS_QUADRATIC;  }
    else if( !strcmp( order_str.c_str(), fe_fams_for_files[ FILES_CONTINUOUS_BIQUADRATIC ].c_str() ) ) {  index = FILES_CONTINUOUS_BIQUADRATIC;  }
    
        return index;
    }

    
    std::string Writer_one_level::get_filename_prefix(const Solution * solution) const {
      
      std::string filename_prefix;
      if( solution != NULL ) filename_prefix = "sol";
      else filename_prefix = "mesh";
    
    return filename_prefix;
    }





} //end namespace femus

