#include "SlepcSVD.hpp"
#include "iostream"

using std::cout;
using std::endl;
using namespace femus;

/*
  This example, taken from slepc examples, computes the singular values 
  of an nxn Grcar matrix, which is a nonsymmetric Toeplitz matrix:

               |  1  1  1  1               |
               | -1  1  1  1  1            |
               |    -1  1  1  1  1         |
               |       .  .  .  .  .       |
           A = |          .  .  .  .  .    |
               |            -1  1  1  1  1 |
               |               -1  1  1  1 |
               |                  -1  1  1 |
               |                     -1  1 |
*/

int main(int argc, char** args) {

  SlepcSVD slepcsvd;
  
  PetscInt       N=30, Istart, Iend, i, col[5];
  PetscScalar    value[] = { -1, 1, 1, 1, 1 };

  // Generate the matrix using petsc code
  Mat A;
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N, N);
  MatSetFromOptions(A);
  MatSetUp(A);

  MatGetOwnershipRange(A,&Istart,&Iend);
  for (i=Istart;i<Iend;i++) {
    col[0]=i-1; col[1]=i; col[2]=i+1; col[3]=i+2; col[4]=i+3;
    if (i==0) {
      MatSetValues(A,1,&i,PetscMin(4,N-i),col+1,value+1,INSERT_VALUES);
    } else {
      MatSetValues(A,1,&i,PetscMin(5,N-i+1),col,value,INSERT_VALUES);
    }
  }

  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  
  // set the linear operator
  slepcsvd.set_operator(A);

 // initialize the context
  slepcsvd.init();
  
  // compute the 2nd norm condition number
  double cond_numb = slepcsvd.compute_2norm_condition_number();
  
  // print the result
  cout << "Estimated condition number: " << cond_numb << endl;
  
  return 0;
}


