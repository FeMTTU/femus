
#include<iostream>
#include<vector>

int main(int argc, char** args) {


  std::cout<<"Hello Erdi, why is this time consuming?"<<std::endl;
  
  std::vector < double > a(10,1);
  
  for (unsigned i=0;i<a.size();i++) {
    std::cout << a[i]<<" ";
  }
  std::cout<<std::endl;
  
  a.resize(5,3);
  
  for (unsigned i=0;i<a.size();i++) {
    std::cout << a[i]<<" ";
  }
  std::cout<<std::endl;
  
  a.resize(8,3.);
  
  for (unsigned i=0;i<a.size();i++) {
    std::cout << a[i]<<" ";
  }
  std::cout<<std::endl;
  
  a.assign(5,3.);
  
  for (unsigned i=0;i<a.size();i++) {
    std::cout << a[i]<<" ";
  }
  std::cout<<std::endl;
  
  

  return 0;
}
