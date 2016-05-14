#include <iostream>
#include <vector>
#include <cmath>

void GetRefinedElement(const std::vector <double > &x0, std::vector <double > &x1,
                       const std::vector <double > &x2, std::vector <double > &x3,
                       std::vector< std::vector <double > > &y);

void GetNewPoint(const std::vector <double > &xi,
                 const std::vector <double > &x0, std::vector <double > &x1,
                 const std::vector <double > &x2, std::vector <double > &x3,
                 std::vector <double > &y);

void GetRotationMatrix(const std::vector< std::vector <double > > &y, double M[3][3]);

int main(int argc, char** args) {

  std::vector < std::vector< double > > y;
  y.resize(10);

  for(int j = 0; j < 10; j++) {
    y[j].resize(3);
  }

  y[0][0] = 0.;
  y[0][1] = 0.;
  y[0][2] = 0.;
  y[1][0] = 1.2;
  y[1][1] = 0.;
  y[1][2] = 0.;
  y[2][0] = -0.2;
  y[2][1] = 1.;
  y[2][2] = 0.;
  y[3][0] = -0.3;
  y[3][1] = -0.1;
  y[3][2] = 1.3;

  for(int k = 0; k < 3; k++) {
    y[4][k] = 0.5 * (y[0][k] + y[1][k]);
    y[5][k] = 0.5 * (y[1][k] + y[2][k]);
    y[6][k] = 0.5 * (y[0][k] + y[2][k]);
    y[7][k] = 0.5 * (y[0][k] + y[3][k]);
    y[8][k] = 0.5 * (y[1][k] + y[3][k]);
    y[9][k] = 0.5 * (y[2][k] + y[3][k]);
  }


  std::vector < std::vector < std::vector< double > > > x;
  x.reserve(1000);

  unsigned n = x.size();
  x.resize(n + 1);
  x[n] = y;

  unsigned levelElements[10];
  levelElements[0] = 0;
  levelElements[1] = 1;
  unsigned numberOfLevels = 2;

  for(unsigned level = 1; level < numberOfLevels; level++) {
    std::cout << "level: " << level << std::endl;
    levelElements[level + 1] = levelElements[level];

    std::cout << "Father Elements: ";

    for(unsigned iel = levelElements[level - 1]; iel < levelElements[level]; iel++) {
      std::cout << iel << " ";

      for(unsigned k = 5; k <= 8; k++) {
        // our refinement
        if(k == 5) {
          GetRefinedElement(x[iel][5], x[iel][6], x[iel][4], x[iel][7], y);
        }
        else if(k == 6) {
          GetRefinedElement(x[iel][8], x[iel][7], x[iel][5], x[iel][4], y);
        }
        else if(k == 7) {
          GetRefinedElement(x[iel][7], x[iel][9], x[iel][8], x[iel][5], y);
        }
        else if(k == 8) {
          GetRefinedElement(x[iel][9], x[iel][5], x[iel][7], x[iel][6], y);
        }

//         // red refinement
//         if(k == 5) {
//           GetRefinedElement(x[iel][4], x[iel][6], x[iel][7], x[iel][8], y);
//         }
//         else if(k == 6) {
//           GetRefinedElement(x[iel][4], x[iel][6], x[iel][5], x[iel][8], y);
//         }
//         else if(k == 7) {
//           GetRefinedElement(x[iel][6], x[iel][7], x[iel][8], x[iel][9], y);
//         }
//         else if(k == 8) {
//           GetRefinedElement(x[iel][6], x[iel][5], x[iel][8], x[iel][9], y);
//         }


        // check new element
        bool newElement = true;

        for(unsigned i = 0; i < x.size(); i++) {
          bool sameVertex[4] = { true, false, false, false};

          for(int j = 1; j < 4; j++) {
            if((y[j][0] - x[i][j][0]) * (y[j][0] - x[i][j][0]) +
                (y[j][1] - x[i][j][1]) * (y[j][1] - x[i][j][1]) +
                (y[j][2] - x[i][j][2]) * (y[j][2] - x[i][j][2]) < 1.0e-14) {
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
          levelElements[level + 1]++;
        }
      }
    }

    std::cout << "\nNumber of new elements: " << levelElements[level + 1] - levelElements[level] << std::endl << std::endl;

    if(levelElements[level + 1] > levelElements[level]) numberOfLevels++;
  }

  unsigned ind[12][4] = {
    {0, 1, 2, 3},
    {0, 2, 3, 1},
    {0, 3, 1, 2},
    {1, 0, 3, 2},
    {1, 2, 0, 3},
    {1, 3, 2, 0},
    {2, 0, 1, 3},
    {2, 1, 3, 0},
    {2, 3, 0, 1},
    {3, 0, 2, 1},
    {3, 1, 0, 2},
    {3, 2, 1, 0}
  };


  //Check for congruence by translation and rotation
  for(unsigned iel = 0; iel < x.size() - 1; iel++) {
    for(unsigned jel = iel + 1; jel < x.size(); jel++) {
      for(unsigned i = 0; i < 12; i++) {

        for(unsigned j = 0; j < 4; j++) {
          for(int k = 0; k < 3; k++) {
            y[j][k] = x[jel][ind[i][j]][k];
          }
        }

        // Translate;
        for(unsigned j = 1; j < 4; j++) {
          for(unsigned k = 0; k < 3; k++) {
            y[j][k] -= y[0][k] ;
          }
        }
        y[0][0] = y[0][1] = y[0][2] = 0.;

        // Build rotation Matrix
        double M[3][3];
        GetRotationMatrix(y, M);

        // Rotate and scale
        for(unsigned j = 0; j < 4; j++) {
          double z[3] = {0., 0., 0.};

          for(unsigned k = 0; k < 3; k++) {
            for(unsigned l = 0; l < 3 ; l++) {
              z[k] += M[k][l] * y[j][l];
            }
          }

          for(unsigned k = 0; k < 3; k++) {
            y[j][k] = z[k];
          }
        }

        // Compare
        bool sameVertex[4] = { true, false, false, false};

        for(int j = 1; j < 4; j++) {
          if((y[j][0] - x[iel][j][0]) * (y[j][0] - x[iel][j][0]) +
              (y[j][1] - x[iel][j][1]) * (y[j][1] - x[iel][j][1]) +
              (y[j][2] - x[iel][j][2]) * (y[j][2] - x[iel][j][2]) < 1.0e-14) {
            sameVertex[j] = true;
          }
        }

        if(sameVertex[1]*sameVertex[2]*sameVertex[3] == true) {
          std::cout << "Element " << jel << " is element " << iel << std::endl;
          x.erase(x.begin() + jel);
          break;
        }
      }
    }
  }


  //Check for congruence by mirroring translation and rotation
  for(unsigned iel = 0; iel < x.size() - 1; iel++) {
    for(unsigned jel = iel + 1; jel < x.size(); jel++) {
      for(unsigned i = 0; i < 12; i++) {

        for(unsigned j = 0; j < 4; j++) {
          y[j][0] = -x[jel][ind[i][j]][0];
          for(int k = 1; k < 3; k++) {
            y[j][k] = x[jel][ind[i][j]][k];
          }
        }

        // Translate;
        for(unsigned j = 1; j < 4; j++) {
          for(unsigned k = 0; k < 3; k++) {
            y[j][k] -= y[0][k] ;
          }
        }
        y[0][0] = y[0][1] = y[0][2] = 0.;

        // Build rotation Matrix
        double M[3][3];
        GetRotationMatrix(y, M);

        // Rotate and scale
        for(unsigned j = 0; j < 4; j++) {
          double z[3] = {0., 0., 0.};

          for(unsigned k = 0; k < 3; k++) {
            for(unsigned l = 0; l < 3 ; l++) {
              z[k] += M[k][l] * y[j][l];
            }
          }

          for(unsigned k = 0; k < 3; k++) {
            y[j][k] = z[k];
          }
        }

        // Compare
        bool sameVertex[4] = { true, false, false, false};

        for(int j = 1; j < 4; j++) {
          if((y[j][0] - x[iel][j][0]) * (y[j][0] - x[iel][j][0]) +
              (y[j][1] - x[iel][j][1]) * (y[j][1] - x[iel][j][1]) +
              (y[j][2] - x[iel][j][2]) * (y[j][2] - x[iel][j][2]) < 1.0e-14) {
            sameVertex[j] = true;
          }
        }

        if(sameVertex[1]*sameVertex[2]*sameVertex[3] == true) {
          std::cout << "Element " << jel << " is element " << iel << std::endl;
          x.erase(x.begin() + jel);
          break;
        }
      }
    }
  }


  std::cout << "The number of congruent tetrahedron families is " << x.size() << std::endl;
  std::cout << std::endl;


//   for(int i = 0; i < x.size(); i++) {
//     std::cout  << std::endl;
//
//     for(int j = 0; j < 3; j++) {
//       std::cout << 2 * i + x[i][j][0] << " " << x[i][j][1] << " " << x[i][j][2] << std::endl;
//     }
//
//     for(int j = 0; j < 3; j++) {
//       std::cout << x[i][j][0] << " " << x[i][j][1] << " " << x[i][j][2] << std::endl;
//       std::cout << x[i][3][0] << " " << x[i][3][1] << " " << x[i][3][2] << std::endl;
//     }
//   }


  return 0;
}


void GetRefinedElement(const std::vector <double > &x0, std::vector <double > &x1,
                       const std::vector <double > &x2, std::vector <double > &x3,
                       std::vector< std::vector <double > > &y) {

  std::vector < std::vector< double > > xi;

  xi.resize(10);

  for(int j = 0; j < 10; j++) {
    xi[j].resize(3);
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

  for(int k = 0; k < 3; k++) {
    xi[4][k] = 0.5 * (xi[0][k] + xi[1][k]);
    xi[5][k] = 0.5 * (xi[1][k] + xi[2][k]);
    xi[6][k] = 0.5 * (xi[0][k] + xi[2][k]);
    xi[7][k] = 0.5 * (xi[0][k] + xi[3][k]);
    xi[8][k] = 0.5 * (xi[1][k] + xi[3][k]);
    xi[9][k] = 0.5 * (xi[2][k] + xi[3][k]);
  }

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

void GetRotationMatrix(const std::vector< std::vector <double > > &y, double M[3][3]) {
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
}
