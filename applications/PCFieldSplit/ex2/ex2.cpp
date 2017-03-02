#include <vector>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
using namespace std;

vector <unsigned> IT(vector <unsigned>& x)
{
 int i;
// cout << x[0] <<endl;
/*
 for (i=0;i < x[0].size();i++){
    cout << "AAAA"<<"  "<<x[0][i]<<endl;
 }
*/
 x[0]=x[0]+1;
 x[1]=x[1]+2;
 return x;
}

int main(int argc, char **argv)
{
double prandtl;
double rayleigh;
unsigned i;
if (argc >= 2) prandtl = strtod (argv[1],NULL);
if (argc >= 3) rayleigh = strtod (argv[2],NULL);

std::vector <unsigned> x[2];
x[0].reserve(3);
x[1].reserve(4);
x[0].resize(2);
x[1].resize(2);
for (i=0;i<x[0].size();i++){
	x[0][i]=i*2;
}
x[0].push_back(10);
x[0].push_back(11);
for (i=0;i<x[0].size();i++){
	std::cout << x[0][i] << std::endl;
}
std::cout <<"111" << std::endl;

std::vector <unsigned> y(10);
//y[0] = IT(x);
cout <<"BBBBB" <<y.size()<<endl;
y=IT(x[0]);
cout <<"AAAAA" <<y.size()<<endl;
for (i=0; i<y.size(); i++){
	std::cout << y[i]<< std::endl;	
}

char *stdOutfile = new  char [100] ;
char *infile = new char [100];

sprintf(stdOutfile, "%s is the most inportant thing","2");
sprintf(infile, "%s is the number %d","2",2);
cout << stdOutfile <<"    "<< infile << endl;
   
//stdOutfile =&"";
printf(stdOutfile, "%s is the most inportant thing","3");
cout << stdOutfile <<"    "<< infile << endl;


cout << "prandtl="  << prandtl << endl;
cout << "rayleigh=" << rayleigh << endl; 
return 0;
}

