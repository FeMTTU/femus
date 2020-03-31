



double hypergeometric(double a, double b, double c, double x)
{
  const double TOLERANCE = 1.0e-10;
  double term = a * b * x / c;
  double value = 1.0 + term;
  int n = 1;

  while(abs(term) > TOLERANCE) {
    a++, b++, c++, n++;
    term *= a * b * x / c / n;
    value += term;
  }

  return value;
}

double Antiderivative1(const double &theta, const double &s, const double &y)
{
  return -(1. / tan(theta) * hypergeometric(0.5, 0.5 - s, 1.5, pow(cos(theta), 2.)) * pow(pow(sin(theta), 2.), 0.5 - s)) /
         (2.* s * pow(y * 1 / sin(theta), 2. * s));
}

double Antiderivative2(const double &theta, const double &s, const double &x)
{
  return (pow(pow(cos(theta), 2.), 0.5 - s) * hypergeometric(0.5, 0.5 - s, 1.5, pow(sin(theta), 2)) * tan(theta)) /
         (2.* s * pow(x * 1. / cos(theta), 2. * s));
}



void GetElementPartition1D(const std::vector <double >  & xg1, const std::vector < std::vector <double > > & x1, const unsigned &split, const unsigned int Nsplit, std::vector < std::vector < std::vector<double>>> &x)
{
  unsigned dim = 1;
  unsigned left = 0;
  unsigned right = 1;

  if(split == 0) { //init
    x.resize(2);
    x[left].resize(dim);
    x[right].resize(dim);
    for(unsigned k = 0; k < dim; k++) {
      x[left][k].resize(x1[0].size());
      x[right][k].resize(x1[0].size());
//       for(unsigned k = 0; k < dim; k++) {
      x[left][k][0] = x1[k][0];
      x[left][k][1] = 0.5 * (x[left][k][0] + xg1[k]);
      x[left][k][2] = 0.5 * (x[left][k][0] + x[left][k][1]);
      x[right][k][1] = x1[k][1];
      x[right][k][0] = 0.5 * (x[right][k][1] + xg1[k]);
      x[right][k][2] = 0.5 * (x[right][k][0] + x[right][k][1]);
//       }
    }
  }
  else if(split == Nsplit) {
    for(unsigned k = 0; k < dim; k++) {
      x[left][k][0] = x[left][k][1];
      x[left][k][1] = xg1[k];
      x[left][k][2] = 0.5 * (x[left][k][0] + x[left][k][1]);

      x[right][k][1] = x[right][k][0];
      x[right][k][0] = xg1[k];
      x[right][k][2] = 0.5 * (x[right][k][0] + x[right][k][1]);
    }
  }
  else {
    for(unsigned k = 0; k < dim; k++) {
      x[left][k][0] = x[left][k][1];
      x[left][k][1] = 0.5 * (x[left][k][0] + xg1[k]);
      x[left][k][2] = 0.5 * (x[left][k][0] + x[left][k][1]);

      x[right][k][1] = x[right][k][0];
      x[right][k][0] = 0.5 * (x[right][k][1] + xg1[k]);
      x[right][k][2] = 0.5 * (x[right][k][0] + x[right][k][1]);
    }
  }
}







void GetElementPartition2D(const std::vector <double >  & xg1, const std::vector < std::vector <double > > & x1, const unsigned &split, const unsigned int Nsplit,  std::vector < std::vector < std::vector<double>>> &x)
{
  unsigned dim = 2;
  unsigned bl = 0; // bottom left
  unsigned br = 1; // bottom right
  unsigned tr = 2; // top left
  unsigned tl = 3; // top right
  unsigned bl1 = 4; // bottom left
  unsigned bl2 = 5; // bottom left
  unsigned br1 = 6; // bottom right
  unsigned br2 = 7; // bottom right
  unsigned tr1 = 8; // top left
  unsigned tr2 = 9; // top left
  unsigned tl1 = 10; // top right
  unsigned tl2 = 11; // top right

  double ex_x_1;
  double ex_x_2;
  double ex_y_1;
  double ex_y_2;

  unsigned size_part = (split != Nsplit) ? 12 : 4;

  if(split == 0) { //init
    x.resize(size_part);
    for(unsigned j = 0; j < size_part; j++) {
      x[j].resize(dim);
      for(unsigned k = 0; k < dim; k++) {
        x[j][k].resize(x1[0].size());
      }
    }
    ex_x_1 = x1[0][0];
    ex_x_2 = x1[0][1];
    ex_y_1 = x1[1][0];
    ex_y_2 = x1[1][3];
  }
  else {
    ex_x_1 = x[bl][0][1];
    ex_x_2 = x[br][0][0];
    ex_y_1 = x[bl][1][2];
    ex_y_2 = x[tl][1][1];
  }


  //     Prototipo: x[quadrante][dim][numero_nodo]
  x[bl][1][0] = x[bl][1][1] = x[br][1][0] = x[br][1][1] = ex_y_1;
  x[tl][1][2] = x[tl][1][3] = x[tr][1][2] = x[tr][1][3] = ex_y_2;
  x[bl][0][0] = x[bl][0][3] = x[tl][0][0] = x[tl][0][3] = ex_x_1;
  x[br][0][1] = x[br][0][2] = x[tr][0][1] = x[tr][0][2] = ex_x_2;


  if(split == Nsplit) {
    x[bl][1][2] = x[bl][1][3] = x[br][1][2] = x[br][1][3] = x[tl][1][0] = x[tl][1][1] = x[tr][1][0] = x[tr][1][1] = xg1[1];
    x[bl][0][1] = x[bl][0][2] = x[tl][0][1] = x[tl][0][2] = x[br][0][0] = x[br][0][3] = x[tr][0][0] = x[tr][0][3] = xg1[0];
  }
  else {
    x[bl][1][2] = x[bl][1][3] = x[br][1][2] = x[br][1][3] = 0.5 * (ex_y_1 + xg1[1]);
    x[tl][1][0] = x[tl][1][1] = x[tr][1][0] = x[tr][1][1] = 0.5 * (ex_y_2 + xg1[1]);
    x[bl][0][1] = x[bl][0][2] = x[tl][0][1] = x[tl][0][2] = 0.5 * (ex_x_1 + xg1[0]);
    x[br][0][0] = x[br][0][3] = x[tr][0][0] = x[tr][0][3] = 0.5 * (ex_x_2 + xg1[0]);
  }

  if(split != Nsplit) {
    x[bl1][1][0] = x[bl1][1][1] = x[br1][1][0] = x[br1][1][1] = ex_y_1;
    x[tl1][1][2] = x[tl1][1][3] = x[tr1][1][2] = x[tr1][1][3] = ex_y_2;
    x[bl1][1][2] = x[bl1][1][3] = x[br1][1][2] = x[br1][1][3] = 0.5 * (ex_y_1 + xg1[1]);
    x[bl2][1][0] = x[bl2][1][1] = x[br2][1][0] = x[br2][1][1] = 0.5 * (ex_y_1 + xg1[1]);
    x[tl1][1][0] = x[tl1][1][1] = x[tr1][1][0] = x[tr1][1][1] = 0.5 * (ex_y_2 + xg1[1]);
    x[tl2][1][2] = x[tl2][1][3] = x[tr2][1][2] = x[tr2][1][3] = 0.5 * (ex_y_2 + xg1[1]);
    x[bl2][1][2] = x[bl2][1][3] = x[br2][1][2] = x[br2][1][3] = xg1[1];
    x[tl2][1][0] = x[tl2][1][1] = x[tr2][1][0] = x[tr2][1][1] = xg1[1];
    x[bl2][0][0] = x[bl2][0][3] = x[tl2][0][0] = x[tl2][0][3] = ex_x_1;
    x[br2][0][1] = x[br2][0][2] = x[tr2][0][1] = x[tr2][0][2] = ex_x_2;
    x[bl2][0][1] = x[bl2][0][2] = x[tl2][0][1] = x[tl2][0][2] = 0.5 * (ex_x_1 + xg1[0]);
    x[bl1][0][0] = x[bl1][0][3] = x[tl1][0][0] = x[tl1][0][3] = 0.5 * (ex_x_1 + xg1[0]);
    x[br2][0][0] = x[br2][0][3] = x[tr2][0][0] = x[tr2][0][3] = 0.5 * (ex_x_2 + xg1[0]);
    x[br1][0][1] = x[br1][0][2] = x[tr1][0][1] = x[tr1][0][2] = 0.5 * (ex_x_2 + xg1[0]);
    x[bl1][0][1] = x[bl1][0][2] = x[tl1][0][1] = x[tl1][0][2] = xg1[0];
    x[br1][0][0] = x[br1][0][3] = x[tr1][0][0] = x[tr1][0][3] = xg1[0];
  }


  for(unsigned qq = 0; qq < size_part; qq++) {
    for(unsigned k = 0; k < dim; k++) { //middle point formula
      x[qq][k][4] = 0.5 * (x[qq][k][0] + x[qq][k][1]);
      x[qq][k][5] = 0.5 * (x[qq][k][1] + x[qq][k][2]);
      x[qq][k][6] = 0.5 * (x[qq][k][2] + x[qq][k][3]);
      x[qq][k][7] = 0.5 * (x[qq][k][3] + x[qq][k][0]);
      x[qq][k][8] = 0.5 * (x[qq][k][0] + x[qq][k][2]);
    }
  }

}


const unsigned ijndex[2][12][2] = {
  { {0, 0}, {3, 0}, {3, 3}, {0, 3},
    {1, 0}, {0, 1},
    {2, 0}, {3, 1},
    {2, 3}, {3, 2},
    {1, 3}, {0, 2}
  },
  {{0, 0}, {1, 0}, {1, 1}, {0, 1}}
};

void GetElementPartitionQuad(const std::vector <double >  & xg1, const std::vector < std::vector <double > > & xNodes, const unsigned & split, const unsigned & totalNumberofSplits,  std::vector < std::vector < std::vector<double>>> &x)
{
  unsigned dim = 2;

  unsigned solType;
  unsigned size = xNodes[0].size();

  if(size == 4) {
    solType = 0; //lagrange linear
  }
  else if(size == 8) {
    solType = 1; //lagrange serendipity
  }
  else if(size == 9) {
    solType = 2; //lagrange quadratic
  }
  else {
    std::cout << "abort in GetElementPartitionQuad" << std::endl;
    abort();
  }


  unsigned bl = 0; // bottom left
  unsigned br = 1; // bottom right
  unsigned tr = 2; // top left
  unsigned tl = 3; // top right

  std::vector < double > XX;
  std::vector < double > YY;

  unsigned size_part = 12;
  unsigned splitType = 0;

  if(split < totalNumberofSplits) { //init && update

    XX.resize(5);
    YY.resize(5);

    if(split == 0) { //init

      x.resize(size_part);
      for(unsigned j = 0; j < size_part; j++) {
        x[j].resize(dim);
        for(unsigned k = 0; k < dim; k++) {
          x[j][k].resize(size);
        }
      }

      XX[0] = xNodes[0][0];
      XX[4] = xNodes[0][1];
      YY[0] = xNodes[1][0];
      YY[4] = xNodes[1][3];

    }
    else { //update
      XX[0] = x[bl][0][1];
      XX[4] = x[br][0][0];
      YY[0] = x[bl][1][2];
      YY[4] = x[tl][1][1];
    }
    XX[2] = xg1[0];
    XX[1] = 0.5 * (XX[0] + XX[2]);
    XX[3] = 0.5 * (XX[2] + XX[4]);

    YY[2] = xg1[1];
    YY[1] = 0.5 * (YY[0] + YY[2]);
    YY[3] = 0.5 * (YY[2] + YY[4]);
  }
  else { //close

    XX.resize(3);
    YY.resize(3);

    XX[0] = x[bl][0][1];
    XX[1] = xg1[0];
    XX[2] = x[br][0][0];
    YY[0] = x[bl][1][2];
    YY[1] = xg1[1];
    YY[2] = x[tl][1][1];

    size_part = 4;
    splitType = 1;
    x.resize(size_part);
    for(unsigned j = 0; j < size_part; j++) {
      x[j].resize(dim);
      for(unsigned k = 0; k < dim; k++) {
        x[j][k].resize(size);
      }
    }
  }

  for(unsigned qq = 0; qq < size_part; qq++) {
    unsigned i = ijndex[splitType][qq][0];
    x[qq][0][0] = x[qq][0][3] = XX[i];
    x[qq][0][1] = x[qq][0][2] = XX[i + 1];

    unsigned j = ijndex[splitType][qq][1];
    x[qq][1][0] = x[qq][1][1] = YY[j];
    x[qq][1][2] = x[qq][1][3] = YY[j + 1];
  }
  if(solType > 0) {
    for(unsigned qq = 0; qq < size_part; qq++) {
      for(unsigned k = 0; k < dim; k++) { //middle point formula
        x[qq][k][4] = 0.5 * (x[qq][k][0] + x[qq][k][1]);
        x[qq][k][5] = 0.5 * (x[qq][k][1] + x[qq][k][2]);
        x[qq][k][6] = 0.5 * (x[qq][k][2] + x[qq][k][3]);
        x[qq][k][7] = 0.5 * (x[qq][k][3] + x[qq][k][0]);
      }
    }
  }

  if(solType > 1) {
    for(unsigned qq = 0; qq < size_part; qq++) {
      for(unsigned k = 0; k < dim; k++) { //middle point formula
        x[qq][k][8] = 0.5 * (x[qq][k][0] + x[qq][k][2]);
      }
    }
  }

}

