#include "SlepcSVD.hpp"

using std::cout;
using std::endl;
using namespace femus;

int main(int argc, char** args) {

  SlepcSVD slepcsvd;
  
  double a = slepcsvd.get_2norm_condition_number();
  
  cout << "Estimated condition number: " << a << endl;
  
  return 0;
}


