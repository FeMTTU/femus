#ifndef __01_OPT_SYSTEM_HPP__
#define __01_OPT_SYSTEM_HPP__

#include <string>

namespace femus {


class pure_boundary {

public:



    //***** How to identify boundary of boundary in the 3D case - BEGIN ******************
  static const std::string  _node_based_bdry_bdry;
    //***** How to identify boundary of boundary in the 3D case - END ******************


};

   const std::string  pure_boundary::_node_based_bdry_bdry = "node_based_bdry_bdry_flag";



class lifting_internal {

public:

protected:

    static constexpr double _lifting_internal_penalty_outside_control_domain = 1.e20; // penalty for zero control outside in equation





};




class lifting_external {

public:

protected:

//Flags for enforcement of Adjoint Neumann continuity on \Gamma_c

static constexpr int _neumann_adjoint_explicit = 1;
static constexpr int _neumann_adjoint_continuity_sign = -1;
static constexpr int _u_minus_q_strong = 1;




};

}

#endif
