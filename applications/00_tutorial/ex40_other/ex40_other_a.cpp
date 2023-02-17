

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <cstdlib>
#include <time.h>

using namespace std;

double *create_matrix(int rows, int cols);
void fill_matrix(double *pt_matrix, int rows, int cols);
double &return_element(double *pt_matrix, int i, int j, int coloums);
void multiply_matrices(double *pt_matrixA, double *pt_matrixB, double *pt_matrixC, int rowsA, int n_AB, int colsB);
void print_matrix(double *pt_matrix, int rows, int cols);

int main(){
  int m=2,n=3,r=2;
  double *pt_A=create_matrix(m,n);
  fill_matrix(pt_A,m,n);
  print_matrix(pt_A,m,n);
  double *pt_B=create_matrix(n,r);
  fill_matrix(pt_B,n,r);
  print_matrix(pt_B,n,r);
  double *pt_C=create_matrix(m,r);
  multiply_matrices(pt_A, pt_B, pt_C, m, n, r);
  print_matrix(pt_C,m,r);
}

double *create_matrix(int rows, int cols){
  
  double *pt= new double[rows*cols];
  return pt;
}

double &return_element(double *pt_matrix, int i, int j, int coloums){
  return *(pt_matrix+coloums*i+j);
}

void fill_matrix(double *pt_matrix, int rows, int cols){
  srand(time(0));
  double *pt_end=pt_matrix+rows*cols;
  while (pt_matrix < pt_end){
    *pt_matrix++ = rand()%10;
  }
}

void multiply_matrices(double *pt_matrixA, double *pt_matrixB, double *pt_matrixC, int rowsA, int n_AB, int colsB){
  double *pta;
  double *ptb;
  double *ptc = pt_matrixC;
  for(int k=0; k < rowsA * colsB; k++){
    pta = pt_matrixA + k /colsB * n_AB;
    ptb = pt_matrixB + k % colsB;
    *ptc= 0.;
    for(int k=0;k<n_AB;k++){
      *ptc += (*pta++) * (*ptb);
       ptb += colsB;
    }
    ptc++;
  }
}
void print_matrix(double *pt_matrix, int rows, int cols){
  for (int i=0;i<rows;i++){
    cout << endl;
    for (int j=0;j<cols;j++){
      cout << return_element(pt_matrix,i,j,cols) << " ";
    }
    cout << "\n" ;
  }
  cout << "-------- \n" ;
}






// #include <iostream>
// #include "adept.h"
// #include <math.h>
// #include "FemusInit.hpp"
// 
// 
// using namespace std;
// using namespace adept;
// using namespace femus;
// 
// adouble f(const adouble x[2]);
// 
// 
// double * AllocateMemory (const unsigned &n){
//     
//  double *p = new double [n];
//  
//  return p;
//     
// }
// 
//  void swp(unsigned &i, unsigned &j){
//    unsigned temp = i;
//    i = j;
//    j = temp;
//    
//    std::cout << i << " " << j <<std::endl;
// };
// 
// 
//  void swpp(unsigned *i, unsigned *j){
//    unsigned temp = *i;
//    *i = *j;
//    *j = temp;
//    
//    std::cout << *i << " " << *j <<std::endl;
// };
// void f(const double &x, double *f_val, double *fder_val);
// int main(int argc, char** args)
// {
// //     // init Petsc-MPI communicator
// //     
// //     unsigned n1;
// //     unsigned n2;
// //     
// //     std::cout << " enter n1 "<<std::endl;
// //     std::cin >> n1;
// //     
// //     if(n1<1) abort();
// //     
// //     std::cout << " enter n2 "<<std::endl;
// //     std::cin >> n2;
// //     
// //     if(n2<1) abort();
// //     
// //     
// //     
// //     double **A;
// //     double *Amemory;
// //     
// //     Amemory = AllocateMemory(n1*n2);
// //     
// //     
// //     
// //     A = new double * [n1];
// //     
// //     for(unsigned i = 0; i < n1; i++){
// //       A[i] = &Amemory[0] + i * n2;  
// //     }
// //     
// //     for(unsigned i = 0; i<n1; i++){
// //       for(unsigned j = 0; j<n2; j++){
// //         A[i][j] = i * n2 + j;    
// //       }     
// //     }
// //     
// //     double *p = Amemory;
// //     for(unsigned i = 0; i<n1; i++){
// //       for(unsigned j = 0; j<n2; j++){
// //         std::cout << *(p + i * n2 +j ) << " ";
// //       }
// //       std::cout << std::endl;
// //     }
// //     
// //     
// //     p = Amemory;
// //     for(unsigned i = 0; i<n1; i++){
// //       for(unsigned j = 0; j<n2; j++){
// //         std::cout << *(p++) << " ";
// //       }
// //       std::cout << std::endl;
// //     }
// //     
// //     
// //     
// //     
// //     
// //     delete [] A;
// //     delete [] Amemory;
// //     
// //     
// //     unsigned i = 5, j=6;
// //     swp(i,j);
// //     
// //     std::cout << i << " " << j <<std::endl;
// //     
// //     swpp(&i,&j);
// //     
// //     std::cout << i << " " << j <<std::endl;
//   
//   double x=1.2;
//   double f_val,fder_val;
//   f(x, &f_val,&fder_val);
//   
//   std::cout << f_val << " " << fder_val <<std::endl;
//   
//   
//   double * f1_val;
//   double *f1der_val;
//   
//   f1_val = &f_val;
//   f1der_val = &fder_val;
//   
//   x = 2;
//   
//   f(x, f1_val, f1der_val);
//   
//   std::cout << f_val << " " << fder_val <<std::endl;
//   
//     
// }
// 
// 
// 
// void f(const double &x, double *f_val, double *fder_val){
//     
//     *f_val= cos(x)-x;
//     *fder_val=-sin(x)-1;
// }
// 
// /** Questions
//  * 1. What is the type of "weight" in Gauss loop? 
//  * 2. What is the type of phi ???
//  * 3. Do I need to know where those 
//  * system.CopySolutionToOldSolution();
//  * system.MGsolve();
//  * are living?
//  */
// 
