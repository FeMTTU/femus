#ifndef SPHERE_H
#define SPHERE_H

#include <iostream>

class Sphere{
public:
  Sphere();
  Sphere(const double & x, const double & y, const double & z, const double & radius);
  void SetCenter(const double &x, const double &y, const double &z);
  double GetXCoordinate(void) const;
  double GetYCoordinate(void) const;
  double GetZCoordinate(void) const;
  void SetRadius(const double &r);
  double GetRadius(void) const;
  double GetVolume(void) const;
  void TranslatedSphere(Sphere *pt, const double & d_x, const double & d_y, const double & d_z) const;
  void TranslatedSphere(Sphere &pt, const double & d_x, const double & d_y, const double & d_z) const;
   
  friend std::ostream& operator<<(std::ostream& os, Sphere& s);
  
  //static void increment_spheres(int);
  //static int total_spheres(void);

private:
  double _xCenter, _yCenter, _zCenter, _radius;
};



#endif // SPHERE_H
