// pointer to functions
#include <iostream>
using namespace std;
double addition (double a, double b)
{return (a+b); }
double subtraction (double a, double b)
{return (a-b);}
double operation (double a , double b, double (*func) (double, double))
{ double g;
  g = (*func)(a,b);
  return g;
}

int main()
{
  double (*mypoint) (double, double), n, m; 
  mypoint = &subtraction;
  n = operation(5.0, 10.0, addition);
  m = operation(19., n, mypoint);
  cout << n <<" "<< m <<endl;
  return 0;
}










/*
int addition (int a, int b)
{ return (a+b); }
int subtraction(int a, int b)
{return (a-b); }

int operation (int x, int y, int (*functocall)(int,int))
{
  int g;
  g = (*functocall)(x,y);
  return (g);
}

int main ()
{
  int m, n;
  int (*minus)(int, int) = subtraction;
  m = operation (7, 5, addition);
  n = operation (20, m, minus);
  cout <<n;
  return 0;
}
// int main(int argc, char **argv)
// {
// double prandtl;
// double rayleigh;
// unsigned i;
// if (argc >= 2) prandtl = strtod (argv[1],NULL);
// if (argc >= 3) rayleigh = strtod (argv[2],NULL);
// 
// std::vector <unsigned> x[2];
// x[0].reserve(3);
// x[1].reserve(4);
// x[0].resize(2);
// x[1].resize(2);
// for (i=0;i<x[0].size();i++){
// 	x[0][i]=i*2;
// }
// x[0].push_back(10);
// x[0].push_back(11);
// for (i=0;i<x[0].size();i++){
// 	std::cout << x[0][i] << std::endl;
// }
// std::cout <<"111" << std::endl;
// 
// std::vector <unsigned> y(10);
// //y[0] = IT(x);
// cout <<"BBBBB" <<y.size()<<endl;
// y=IT(x[0]);
// cout <<"AAAAA" <<y.size()<<endl;
// for (i=0; i<y.size(); i++){
// 	std::cout << y[i]<< std::endl;	
// }
// 
// char *stdOutfile = new  char [100] ;
// char *infile = new char [100];
// 
// sprintf(stdOutfile, "%s is the most inportant thing","2");
// sprintf(infile, "%s is the number %d","2",2);
// cout << stdOutfile <<"    "<< infile << endl;
//    
// //stdOutfile =&"";
// printf(stdOutfile, "%s is the most inportant thing","3");
// cout << stdOutfile <<"    "<< infile << endl;
// 
// 
// cout << "prandtl="  << prandtl << endl;
// cout << "rayleigh=" << rayleigh << endl; 
// return 0;
// }*/
