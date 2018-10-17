
#include <iostream>
#include "adept.h"
#include <math.h>
#include "FemusInit.hpp"


using namespace std;
using namespace adept;
using namespace femus;

adouble f(const adouble x[2]);

int main(int argc, char** args)
{
    // init Petsc-MPI communicator
    FemusInit mpinit(argc, args, MPI_COMM_WORLD);
    adept::Stack& s = FemusInit::_adeptStack;
    adouble x[2]={2.0,3.0};
    s.new_recording();
    adouble y=f(x);
    y.set_gradient(1.0);
    s.compute_adjoint();
    cout << y.value() << endl;
    return 0;
}

adouble f(const adouble x[2]){
    
    adouble y=x[0]*x[0]+x[1]*x[1];
    return y;
}

/** Questions
 * 1. What is the type of "weight" in Gauss loop? 
 * 2. What is the type of phi ???
 * 3. Do I need to know where those 
 * system.CopySolutionToOldSolution();
 * system.MGsolve();
 * are living?
 */

