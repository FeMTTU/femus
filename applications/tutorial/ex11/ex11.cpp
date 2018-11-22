#include "sphere.hpp"
#include "array.hpp"
#include <iostream>
#include <math.h>

using namespace std;



int main(){
  
}



// Sphere::Sphere():
// _xCenter(0.),
// _yCenter(0.),
// _zCenter(0.),
// _radius(1.){ 
// };
// 
// Sphere::Sphere(const double & x, const double & y, const double & z, const double & radius): // Constructor goes here.
// _xCenter(x),
// _yCenter(y),
// _zCenter(z),
// _radius(radius){
// }
// 
// void Sphere::SetCenter(const double & x, const double & y, const double & z){
//     _xCenter=x;
//     _yCenter=y;
//     _zCenter=z;
// }
// 
// double Sphere::GetXCoordinate(void) const {
//   return _xCenter;
// }
// 
// double Sphere::GetYCoordinate(void) const {
//   return _yCenter;
// }
// 
// double Sphere::GetZCoordinate(void) const {
//   return _zCenter;
// }
// 
// void Sphere::SetRadius(const double & r){
//   _radius=r;
// }
// 
// double Sphere::GetRadius() const {
//   return _radius;
// }
// 
// 
// double Sphere::GetVolume(void) const {
//   return 4 * M_PI * pow(_radius,3)/3.0;
// }
// 
// void Sphere::TranslatedSphere(Sphere *pt, const double & d_x, const double & d_y, const double & d_z) const {
//   pt->_xCenter =_xCenter + d_x; //(*pt).x_centre=x_centre+d_x
//   pt->_yCenter =_yCenter + d_y;
//   pt->_zCenter =_zCenter + d_z;
//   pt->_radius =_radius;
// }
// 
// void Sphere::TranslatedSphere(Sphere &pt, const double & d_x, const double & d_y, const double & d_z) const {
//   pt._xCenter =_xCenter + d_x; //(*pt).x_centre=x_centre+d_x
//   pt._yCenter =_yCenter + d_y;
//   pt._zCenter =_zCenter + d_z;
//   pt._radius =_radius;
// }
// 
// 
// std::ostream& operator<<(std::ostream& os, Sphere& s) {
//   std::cout << "h = "<< s.GetXCoordinate() << ", k = " << s.GetYCoordinate()<<", l = " << s.GetZCoordinate() << ", r = " << s.GetRadius();
// }
// 
// int main(){
//     
//     Sphere s1;
//     
//     Sphere s2(1.,2.,3.,4.);
//     s1.SetCenter(1,2,4);
//     s1.SetRadius(3);
//     cout << s1.GetVolume() << "\n" ;
//     
//     std::cout << s2 << std::endl;
//     
//     s2.TranslatedSphere(s1,1.,1.,1.);
//       
//     std::cout << s1 << std::endl;
//     
//     
// }


