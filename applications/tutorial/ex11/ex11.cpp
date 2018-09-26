
#include<iostream>

int main(int argc, char** args) {


  std::cout<<"Hello Erdi, why is this time consuming?"<<std::endl;

  return 0;
}

// Ask the Big Boss!
// 1. const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim))); conservative: based on line3, quad9, hex27
// Guess: always lagrange so in 1d 3nodes, in  2d 9nodes, in 3d 27 nodes... If so I am a genius!!!
//Answ: 
// 2.   unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
//      unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs
// Guess: These are same? 
// 3. x[k].resize(nDofx); 
// we already resize x above why again?
// 4. How to create new examples?
// Answer:
// 5.Ex6 117--unsigned solVType = mlSol->GetSolutionType(solVIndex[0]); 
// Should we also get for solVIndex[1] which is "v"?
// Answer:
// 6. Ex6 220   sysDof.reserve((dim + 1) * maxSize); that 1 comes from pressure??
// ans:



// Dont be lazy look at:
// You are not good at finite element types, you dont understand. 
