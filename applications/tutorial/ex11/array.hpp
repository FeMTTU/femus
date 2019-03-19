#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>

class array{
public:
  array(const int &size){
    _n = size;
    _pt = new double [_n];
  }
  
  ~array(){
    delete [] _pt;
  }
  
//   Sphere();
//   Sphere(const double & x, const double & y, const double & z, const double & radius);
//   void SetCenter(const double &x, const double &y, const double &z);
//   double GetXCoordinate(void) const;
//   double GetYCoordinate(void) const;
//   double GetZCoordinate(void) const;
//   void SetRadius(const double &r);
//   double GetRadius(void) const;
//   double GetVolume(void) const;
//   void TranslatedSphere(Sphere *pt, const double & d_x, const double & d_y, const double & d_z) const;
//   void TranslatedSphere(Sphere &pt, const double & d_x, const double & d_y, const double & d_z) const;
//    
//   friend std::ostream& operator<<(std::ostream& os, Sphere& s);
//   
//   //static void increment_spheres(int);
//   //static int total_spheres(void);

protected:
  int _n;
  double *_pt;
//   double _xCenter, _yCenter, _zCenter, _radius;
};



#endif // SPHERE_H
