#include <vector>
#include <iostream>

/*
std::vector <unsigned> IT(const std::vector <unsigned> &x)
{
return (x[0]);
}
*/
int main(int argc, char **argv)
{
unsigned i;
std::vector <unsigned> x[2];
x[0].reserve(3);
x[1].reserve(4);
x[0].resize(2);
x[1].resize(2);
for (i=0;i<x[0].size();i++){
	x[0][i]=i;
}
x[0].push_back(10);
x[0].push_back(11);
for (i=0;i<x[0].size();i++){
	std::cout << x[0][i] << std::endl;
}
std::cout <<"111" << std::endl;

std::vector <unsigned> y[2];
//y[0] = IT(x);
y[0]=x[0];
for (i=0; i<y[0].size(); i++){
	std::cout << y[0][i]<< std::endl;	
}
return 0;
}

