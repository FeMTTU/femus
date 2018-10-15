
<<<<<<< HEAD
// example: class constructor
#include <iostream>
using namespace std;

class Rectangle {
    int width, height;
  public:
    Rectangle();         // This is default constructor which is called when an object is declared but is not initialized with any arguments.
    Rectangle (int,int); // This is a constructor which initialize the class.
    int area () {return (width*height);}
};
Rectangle::Rectangle () {
  width = 5;
  height = 5;
}
Rectangle::Rectangle (int a, int b) { // This is how we define a constructor, no return type!
  width = a;
  height = b;
}

int main () {
  Rectangle rect; // Notice that this object has no argument.
  Rectangle rectb (5,6);
  cout << "rect area: " << rect.area() << endl;
  cout << "rectb area: " << rectb.area() << endl;
  return 0;
}
=======
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
    adouble x[2]={1.0,2.0};
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
>>>>>>> d4d9dbd30d2ff22cbdb425df655e70954e509ed1
