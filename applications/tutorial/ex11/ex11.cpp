
#include <iostream>
#include "adept.h"
#include <math.h>
#include "FemusInit.hpp"


using namespace std;
using namespace adept;
using namespace femus;

adouble f(const adouble x[2]);


double * AllocateMemory (const unsigned &n){
    
 double *p = new double [n];
 
 return p;
    
}

 void swp(unsigned &i, unsigned &j){
   unsigned temp = i;
   i = j;
   j = temp;
   
   std::cout << i << " " << j <<std::endl;
};


 void swpp(unsigned *i, unsigned *j){
   unsigned temp = *i;
   *i = *j;
   *j = temp;
   
   std::cout << *i << " " << *j <<std::endl;
};
void f(const double &x, double *f_val, double *fder_val);
int main(int argc, char** args)
{
//     // init Petsc-MPI communicator
//     
//     unsigned n1;
//     unsigned n2;
//     
//     std::cout << " enter n1 "<<std::endl;
//     std::cin >> n1;
//     
//     if(n1<1) abort();
//     
//     std::cout << " enter n2 "<<std::endl;
//     std::cin >> n2;
//     
//     if(n2<1) abort();
//     
//     
//     
//     double **A;
//     double *Amemory;
//     
//     Amemory = AllocateMemory(n1*n2);
//     
//     
//     
//     A = new double * [n1];
//     
//     for(unsigned i = 0; i < n1; i++){
//       A[i] = &Amemory[0] + i * n2;  
//     }
//     
//     for(unsigned i = 0; i<n1; i++){
//       for(unsigned j = 0; j<n2; j++){
//         A[i][j] = i * n2 + j;    
//       }     
//     }
//     
//     double *p = Amemory;
//     for(unsigned i = 0; i<n1; i++){
//       for(unsigned j = 0; j<n2; j++){
//         std::cout << *(p + i * n2 +j ) << " ";
//       }
//       std::cout << std::endl;
//     }
//     
//     
//     p = Amemory;
//     for(unsigned i = 0; i<n1; i++){
//       for(unsigned j = 0; j<n2; j++){
//         std::cout << *(p++) << " ";
//       }
//       std::cout << std::endl;
//     }
//     
//     
//     
//     
//     
//     delete [] A;
//     delete [] Amemory;
//     
//     
//     unsigned i = 5, j=6;
//     swp(i,j);
//     
//     std::cout << i << " " << j <<std::endl;
//     
//     swpp(&i,&j);
//     
//     std::cout << i << " " << j <<std::endl;
  
  double x=1.2;
  double f_val,fder_val;
  f(x, &f_val,&fder_val);
  
  std::cout << f_val << " " << fder_val <<std::endl;
  
  
  double * f1_val;
  double *f1der_val;
  
  f1_val = &f_val;
  f1der_val = &fder_val;
  
  x = 2;
  
  f(x, f1_val, f1der_val);
  
  std::cout << f_val << " " << fder_val <<std::endl;
  
    
}



void f(const double &x, double *f_val, double *fder_val){
    
    *f_val= cos(x)-x;
    *fder_val=-sin(x)-1;
}

/** Questions
 * 1. What is the type of "weight" in Gauss loop? 
 * 2. What is the type of phi ???
 * 3. Do I need to know where those 
 * system.CopySolutionToOldSolution();
 * system.MGsolve();
 * are living?
 */

