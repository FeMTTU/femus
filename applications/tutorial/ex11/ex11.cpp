#include "sphere.hpp"
#include <iostream>
#include <math.h>

using namespace std;


sphere::sphere() // Constructor goes here.
{
    
}
void sphere::assign_centre(double x, double y, double z){
    x_centre=x;
    y_centre=y;
    z_centre=z;
}
double sphere::get_x_coordinate(void){
    return x_centre;
}
double sphere::get_y_coordinate(void){
    return y_centre;
}
double sphere::get_z_coordinate(void){
    return z_centre;
}

 void sphere::assign_radius(double r){
    radius=r;
}
 double sphere::get_radius(){
    return radius;
}

double sphere::volume(void){
    return 4*M_PI*pow(radius,3)/3.0;
}

void sphere::translated_sphere(sphere *pt, double d_x, double d_y, double d_z){
    pt->x_centre=x_centre+d_x; //(*pt).x_centre=x_centre+d_x
    pt->y_centre=y_centre+d_y;
    pt->z_centre=z_centre+d_z;
    pt->radius=radius;
}

int main(){
    
    sphere s1;
    s1.assign_centre(1,2,4);
    s1.assign_radius(3);
    cout << s1.volume() << "\n" ;
}


