#ifndef SPHERE_H
#define SPHERE_H

class sphere{
public:
    sphere();
    void assign_centre(double x, double y, double z);
    double get_x_coordinate(void);
    double get_y_coordinate(void);
    double get_z_coordinate(void);
    void assign_radius(double r);
    double get_radius(void);
    double volume(void);
    void translated_sphere(sphere *,double,double,double);
    //static void increment_spheres(int);
    //static int total_spheres(void);

private:
    double x_centre,y_centre,z_centre, radius;
};

#endif // SPHERE_H
