#include <iostream>
#include <vector>
#include <cmath>

void GetRefinedElement(const std::vector< std::vector <double > > &xi,
                       const std::vector <double > &x0, std::vector <double > &x1,
                       const std::vector <double > &x2, std::vector <double > &x3,
                       std::vector< std::vector <double > > &y);

void GetNewPoint(const std::vector <double > &xi,
                 const std::vector <double > &x0, std::vector <double > &x1,
                 const std::vector <double > &x2, std::vector <double > &x3,
                 std::vector <double > &y);


int main(int argc, char** args) {

  std::vector < std::vector< double > > xi;

  xi.resize(10);

  for(int j = 0; j < 10; j++) {
    xi[ j ].resize(3);
  }

  xi[0][0] = 0.;
  xi[0][1] = 0.;
  xi[0][2] = 0.;
  xi[1][0] = 1.;
  xi[1][1] = 0.;
  xi[1][2] = 0.;
  xi[2][0] = 0.;
  xi[2][1] = 1.;
  xi[2][2] = 0.;
  xi[3][0] = 0.;
  xi[3][1] = 0.;
  xi[3][2] = 1.;
  xi[4][0] = 0.5;
  xi[4][1] = 0.;
  xi[4][2] = 0.;
  xi[5][0] = 0.5;
  xi[5][1] = 0.5;
  xi[5][2] = 0.;
  xi[6][0] = 0.;
  xi[6][1] = 0.5;
  xi[6][2] = 0.;
  xi[7][0] = 0.;
  xi[7][1] = 0.;
  xi[7][2] = 0.5;
  xi[8][0] = 0.5;
  xi[8][1] = 0.;
  xi[8][2] = 0.5;
  xi[9][0] = 0.;
  xi[9][1] = 0.5;
  xi[9][2] = 0.5;

  std::vector < std::vector < std::vector< double > > > x;
  x.reserve(1000);

  unsigned n = x.size();
  x.resize(n + 1);
  x[n] = xi;



  std::vector < std::vector< double > > y;
  y.resize(10);

  for(int j = 0; j < 10; j++) {
    y[ j ].resize(3);
  }


  unsigned levelElements[10];
  levelElements[0] = 0;
  levelElements[1] = 1;
  unsigned numberOfLevels = 2;

  for(unsigned level = 1; level < numberOfLevels; level++) {
    std::cout << "level= "<< level << std::endl;
    levelElements[level+1] = levelElements[level];
    std::cout << levelElements[level]<< " " << levelElements[level+1] <<std::endl;
    for(unsigned counter = levelElements[level - 1]; counter < levelElements[level]; counter++){
      std::cout << counter << " ";
      for(unsigned k = 5; k <= 8; k++) {
        if(k == 5) {
          GetRefinedElement(xi, x[counter][5], x[counter][6], x[counter][4], x[counter][7], y);
        }
        else if(k == 6) {
          GetRefinedElement(xi, x[counter][8], x[counter][7], x[counter][5], x[counter][4], y);
        }
        else if(k == 7) {
          GetRefinedElement(xi, x[counter][7], x[counter][9], x[counter][8], x[counter][5], y);
        }
        else if(k == 8) {
          GetRefinedElement(xi, x[counter][9], x[counter][5], x[counter][7], x[counter][6], y);
        }

        // check new element
        bool newElement = true;
        for(unsigned i = 0; i < x.size(); i++) {
          bool sameVertex[4] = { true, false, false, false};
          for(int j = 1; j < 4; j++) {
            if((y[j][0] - x[i][j][0]) * (y[j][0] - x[i][j][0]) +
                (y[j][1] - x[i][j][1]) * (y[j][1] - x[i][j][1]) +
                (y[j][2] - x[i][j][2]) * (y[j][2] - x[i][j][2]) < 1.0e-12) {
              sameVertex[j] = true;
            }
          }

          if(sameVertex[1]*sameVertex[2]*sameVertex[3] == true) {
            newElement = false;
            break;
          }
        }

        if(newElement) {
          n = x.size();
          x.resize(n + 1);
          x[n] = y;
          levelElements[level+1]++;
        }
      }
    }
    std::cout << "\n" << levelElements[level]<< " " << levelElements[level+1] <<std::endl;
    if( levelElements[level+1] > levelElements[level] ) numberOfLevels++;
  }
/*
  std::cout << std::endl;


  for(int i = 0; i < x.size(); i++) {
    std::cout << i << std::endl;

    for(int j = 0; j < 10; j++) {
      std::cout << x[i][j][0] << " " << x[i][j][1] << " " << x[i][j][2] << std::endl;
    }
  }*/


  return 0;
}


void GetRefinedElement(const std::vector< std::vector <double > > &xi,
                       const std::vector <double > &x0, std::vector <double > &x1,
                       const std::vector <double > &x2, std::vector <double > &x3,
                       std::vector< std::vector <double > > &y) {
  for(unsigned i = 0; i < 10; i++) {
    GetNewPoint(xi[i], x0, x1, x2, x3, y[i]);
  }

  // Translate;
  for(unsigned i = 1; i < 10; i++) {
    for(unsigned k = 0; k < 3; k++) {
      y[i][k] -= y[0][k] ;
    }
  }

  y[0][0] = y[0][1] = y[0][2] = 0.;

  // Build rotation Matrix
  double M[3][3];
  double norm = 0;

  for(unsigned k = 0; k < 3; k++) {
    M[0][k] = y[1][k];
    norm += M[0][k] * M[0][k];
    M[1][k] = y[2][k];
  }

  norm = sqrt(norm);

  for(unsigned k = 0; k < 3; k++) M[0][k] /= norm;

  M[2][0] = M[0][1] * M[1][2] - M[0][2] * M[1][1];
  M[2][1] = M[0][2] * M[1][0] - M[0][0] * M[1][2];
  M[2][2] = M[0][0] * M[1][1] - M[0][1] * M[1][0];
  norm = 0;

  for(unsigned k = 0; k < 3; k++) {
    norm += M[2][k] * M[2][k];
  }

  norm = sqrt(norm);

  for(int k = 0; k < 3; k++) M[2][k] /= norm;

  M[1][0] = M[2][1] * M[0][2] - M[2][2] * M[0][1];
  M[1][1] = M[2][2] * M[0][0] - M[2][0] * M[0][2];
  M[1][2] = M[2][0] * M[0][1] - M[2][1] * M[0][0];

  // Rotate and scale
  for(unsigned i = 0; i < 10; i++) {
    double z[3] = {0., 0., 0.};

    for(unsigned k = 0; k < 3; k++) {
      for(unsigned l = 0; l < 3 ; l++) {
        z[k] += 2. * M[k][l] * y[i][l];
      }
    }

    for(unsigned k = 0; k < 3; k++) {
      y[i][k] = z[k];
    }
  }


}

void GetNewPoint(const std::vector <double > &xi,
                 const std::vector <double > &x0, std::vector <double > &x1,
                 const std::vector <double > &x2, std::vector <double > &x3,
                 std::vector <double > &y) {
  for(unsigned k = 0; k < 3; k++) {
    y[k] = (1 - xi[0] - xi[1] - xi[2]) * x0[k] + xi[0] * x1[k] + xi[1] * x2[k] + xi[2] * x3[k];
  }
}
