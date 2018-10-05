
#include <iostream>
#include <adept.h>
#include <vector>
using namespace std;
using namespace adept;

adouble f(const adouble x[2]);
//double alg_and_grad(const double x[2]);

int main() {
    Stack s;
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

// double alg_and_grad(const double x_val[2], double dy_dx[2]){
//     
//     Stack s;
//     adouble x[2]={x_val[0],x_val[1]};
//     s.new_recording();
//     adouble y=f(x);
//     y.set_gradient(1.0);
//     s.compute_adjoint();
//     dy_dx[0] = x[0].get_gradient();
//     dy_dx[1] = x[1].get_gradient();
//     return y.value();
//     
// }
