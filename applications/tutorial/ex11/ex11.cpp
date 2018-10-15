
// example: class constructor
#include <iostream>
using namespace std;

class Rectangle {
    int width, height;
  public:
    Rectangle();         // This is default constructor which is called when an object is declared but is not initialized with any arguments.
    Rectangle (int,int); // This is a constructor which initialize the class.
    int area () {return (width*height);}
};
Rectangle::Rectangle () {
  width = 5;
  height = 5;
}
Rectangle::Rectangle (int a, int b) { // This is how we define a constructor, no return type!
  width = a;
  height = b;
}

int main () {
  Rectangle rect; // Notice that this object has no argument.
  Rectangle rectb (5,6);
  cout << "rect area: " << rect.area() << endl;
  cout << "rectb area: " << rectb.area() << endl;
  return 0;
}
