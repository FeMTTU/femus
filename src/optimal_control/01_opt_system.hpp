#ifndef __01_OPT_SYSTEM_HPP__
#define __01_OPT_SYSTEM_HPP__

#include <string>

namespace femus {


class pure_boundary {

public:



    //***** How to identify boundary of boundary in the 3D case - BEGIN ******************
  static const std::string  _node_based_bdry_bdry ;
    //***** How to identify boundary of boundary in the 3D case - END ******************


};

   const std::string  pure_boundary::_node_based_bdry_bdry = "node_based_bdry_bdry_flag";



class lifting_internal {

public:

protected:

    static constexpr double _lifting_internal_penalty_outside_control_domain = 1.e20; // penalty for zero control outside in equation





};


}

#endif
