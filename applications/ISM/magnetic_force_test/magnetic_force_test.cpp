#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "FElemTypeEnum.hpp"
#include "Files.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"
#include "Line.hpp"
#include "adept.h"
#include "../../FSI/vascular/include/FSITimeDependentAssemblySupg.hpp"
#include <cmath>


#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

using boost::math::ellint_1;
using boost::math::ellint_2;

double scale = 1000.;

using namespace std;
using namespace femus;

void MagneticForce(const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned &material, unsigned forceType);
//------------------------------------------------------------------------------------------------------------------

int main(int argc, char **args) {



  unsigned nx= 1;
  unsigned ny = 20;
  
  std::vector < std::vector < double > >  x(nx*ny);
  std::vector < double > Fm(3,0);

 
  for(unsigned i = 0; i < nx; i++) {
    for(unsigned j = 0; j < ny ;j++) {
      x[i * ny + j].resize(3); 
//       double r = i * 0.01;
//       double theta = 2 * (acos(-1) / ny) * j; 
//       x[i * ny + j][0] = 1. + r*cos(theta);
//       x[i * ny + j][1] = 2. + r*sin(theta);
//       x[i * ny + j][2] = 0.1;
      x[i * ny + j][0] = 0.;
      x[i * ny + j][1] = -0.1 + 0.01*j;
      x[i * ny + j][2] = 0.;
      
      MagneticForce(x[i * ny + j], Fm, 2, 1);
      double samsing = ( 0.04 * 0.04 + x[i * ny + j][2] * x[i * ny + j][2] );
      double Bz = (1.e12 * 0.04 * 0.04) / ( 2 * sqrt(samsing * samsing * samsing) );
      double dBz = (-3 * x[i * ny + j][2] * 1.e12 * 0.04 * 0.04 ) / (2 *  sqrt(samsing * samsing * samsing * samsing * samsing) );
      std::cout <<   x[i * ny + j][0] << " " << x[i * ny + j][1] << " " << x[i * ny + j][2] << std::endl;
      std::cout <<   x[i * ny + j][0] + Fm[0] << " " << x[i * ny + j][1] + Fm[1] << " " << x[i * ny + j][2] + Fm[2] << std::endl;
      std::cout << std::endl << std::endl;
    }
     std::cout << std::endl;
 }







  return 0;
}






void MagneticForce(const std::vector <double> & xMarker, std::vector <double> &Fm, const unsigned &material, unsigned forceType) {

  // std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

  //case 0: infinitely long wire with a current I flowing modelled by the line identified by x and v
  //case 1: current loop of radius a, center x and for z-axis the line identified by x and v

  //BEGIN magnetic and electric parameters

  double PI = acos(-1.);
  double I = 1.857e5; // electric current intensity
  double Msat = 1.e6;  //  magnetic saturation
  double  chi = 3.; //magnetic susceptibility
  double mu0 = 4 * PI * 1.e-7;  //magnetic permeability of the vacuum
  double H;
  std::vector <double> gradH(3);
  std::vector <double> gradHSquared(3);
  std::vector<double> vectorH(3);
  double H0 = Msat / chi;

  //END


  //BEGIN fluid viscosity

  double muf = (material == 2) ? 3.5 * 1.0e-3 : 1.0e100; // fluid viscosity

  //END


  //BEGIN geometric parameters

  double D = 500 * 1.e-9;       //diameter of the particle

  double a = 0.04; //radius of the circular current loop in m

  std::vector <double> v(3);   //case 0: direction vector of the line that identifies the infinite wire
  //case 1: z axis of the current loop (symmetry axis)


//     aortic bifurcation wire
//     v[0] = 0.;
//     v[1] = 0.;
//     v[2] = 1.;



  // aortic bifurcation current loop
//     v[0] = 0.;
//     v[1] = 0.;
//     v[2] = 1.;


//     tube 3D
    v[0] = -1.;
    v[1] = 0.;
    v[2] = 0.;


  //bent tube no FSI wire
//     v[0] = 0.;
//     v[1] = 1.;
//     v[2] = 0.;

  std::vector <double> x(3);   //case 0: point that with v identifies the line of the wire
  //case 1: center of the current loop

  //aeurysm bifurcation wire
//   x[0] = -0.01;
//   x[1] = 0.008;
//   x[2] = 0.;


  //aortic bifurcation wire
//     x[0] = 0.015;
//     x[1] = 0.;
//     x[2] = 0.;



//    //aortic bifurcation current loop
//     x[0] = 0.055;
//     x[1] = 0.;
//     x[2] = 0.;

  //tube 3D
    x[0] = 0.425 - 0.007;
    x[1] = 0.;
    x[2] = 0.;


//    //bent tube no FSI wire
//     x[0] = 9.;
//     x[1] = 0.;
//     x[2] = 3.;


  //END


  //BEGIN extraction of the coordinates of the particle

  std::vector <double> xM(3);
  xM[0] = xMarker[0];
  xM[1] = xMarker[1];
  xM[2] = (xMarker.size() == 3) ? xMarker[2] : 0. ;

  //END


  //BEGIN evaluate H

  switch(forceType) {

    case 0: {  // infinite wire

      double Gamma;
      double Omega;
      std::vector<double> gradOmega(3);

      Gamma = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      Omega = (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) * (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) +
              (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2])) * (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2])) +
              (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1])) * (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1])) ;

      gradOmega[0] = 2 * v[2] * (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2])) + 2 * v[1] * (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1]));

      gradOmega[1] = 2 * v[2] * (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) - 2 * v[0] * (v[1] * (x[0] - xM[0]) - v[0] * (x[1] - xM[1]));

      gradOmega[2] = - 2 * v[1] * (v[2] * (x[1] - xM[1]) - v[1] * (x[2] - xM[2])) - 2 * v[0] * (v[2] * (x[0] - xM[0]) - v[0] * (x[2] - xM[2]));

      H = (I / (2 * PI)) * Gamma * (1 / sqrt(Omega));

      gradHSquared[0] = (I / (2 * PI)) * (I / (2 * PI)) * ((-1) * Gamma * Gamma * (1 / (Omega * Omega))) * gradOmega[0];
      gradHSquared[1] = (I / (2 * PI)) * (I / (2 * PI)) * ((-1) * Gamma * Gamma * (1 / (Omega * Omega))) * gradOmega[1];
      gradHSquared[2] = (I / (2 * PI)) * (I / (2 * PI)) * ((-1) * Gamma * Gamma * (1 / (Omega * Omega))) * gradOmega[2];

      gradH[0] = (I / (2 * PI)) * ((-0.5) * Gamma * gradOmega[0]) / sqrt(Omega * Omega * Omega);
      gradH[1] = (I / (2 * PI)) * ((-0.5) * Gamma * gradOmega[1]) / sqrt(Omega * Omega * Omega);
      gradH[2] = (I / (2 * PI)) * ((-0.5) * Gamma * gradOmega[2]) / sqrt(Omega * Omega * Omega);

    }
    break;

    case 1: {  //current loop

      //BEGIN bulid the rotation Matrix;
      double v2 = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      v[0] /= v2;
      v[1] /= v2;
      v[2] /= v2;

      std::vector <double> u(3); // will be the x'-axis
      std::vector <double> w(3); // will be the y'-axis

      unsigned imax = 0;
      double vmax = fabs(v[0]);
      for(unsigned i = 1; i < 3; i++) {
        if(fabs(v[i]) > vmax) {
          imax = i;
          vmax = fabs(v[i]);
        }
      }

      u[0] = (imax == 0) ? - (v[1] + v[2]) / v[0]  : 1;
      u[1] = (imax == 1) ? - (v[2] + v[0]) / v[1]  : 1;
      u[2] = (imax == 2) ? - (v[0] + v[1]) / v[2]  : 1;
      
//       u[0]=1;
//       u[1]=0;
//       u[2]=0;

      double u2 = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
      u[0] /= u2;
      u[1] /= u2;
      u[2] /= u2;

      w[0] = v[1] * u[2] - v[2] * u[1];
      w[1] = v[2] * u[0] - v[0] * u[2];
      w[2] = v[0] * u[1] - v[1] * u[0];

      double w2 = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
      if(fabs(w2 - 1.) > 1.0e-12) {
        std::cout << "ERRORE AAAAAAAAAAAAAAAAAAA" << std::endl;
      }

      std::vector< std::vector <double> > R(3);
      for(unsigned i = 0; i < 3; i++) {
        R[i].resize(3);
        R[i][0] = u[i];
        R[i][1] = w[i];
        R[i][2] = v[i];
      }
      //END bulid the rotation Matrix;

//           std::cout<<std::endl;
//           for(unsigned i=0;i<3;i++){
// 	    for(unsigned j=0;j<3;j++){
// 	      std::cout << R[i][j] <<  " ";
// 	    }
// 	    std::cout<<std::endl;
// 	  }


      //BEGIN find marker local coordinates
      for(unsigned i = 0; i < 3 ; i++) {
        xM[i] -= x[i];
      }

      
      
      
      std::vector<double> xM1(3, 0.);
      for(unsigned i = 0; i < 3; i++) {
        for(unsigned j = 0; j < 3; j++) {
          xM1[i] += R[j][i] * xM[j];
        }
      }
      xM = xM1;


      //END find marker local coordinates

      double rhoSquared = xM[0] * xM[0] + xM[1] * xM[1];
      double rSquared = rhoSquared + xM[2] * xM[2];
      double alphaSquared = a * a + rSquared - 2 * a * sqrt(rhoSquared);
      double betaSquared = a * a + rSquared + 2 * a * sqrt(rhoSquared);
      double kSquared = 1 - (alphaSquared / betaSquared);
      double gamma = xM[0] * xM[0] - xM[1] * xM[1];
      double C = (I * mu0) / PI;



      std::vector < double > vectorHp(3);
      std::vector < std::vector < double > > jacobianVectorHp(3);
      for(unsigned i = 0; i < 3; i++) {
        jacobianVectorHp[i].resize(3);
      }

      bool nearZAxis = ( rhoSquared < 1e-5 * a * a ) ? true : false;
      
      vectorHp[0] = ( nearZAxis == true) ?  0 : (1 / mu0) * ((C * xM[0] * xM[2]) / (2 * alphaSquared * sqrt(betaSquared) * rhoSquared)) * ((a * a + rSquared) * ellint_2(kSquared) - alphaSquared * ellint_1(kSquared));
      vectorHp[1] = ( nearZAxis == true) ?  0 : (1 / mu0) * ((C * xM[1] * xM[2]) / (2 * alphaSquared * sqrt(betaSquared) * rhoSquared)) * ((a * a + rSquared) * ellint_2(kSquared) - alphaSquared * ellint_1(kSquared));;
      vectorHp[2] = (1 / mu0) * (C / (2 * alphaSquared * sqrt(betaSquared))) * ((a * a - rSquared) * ellint_2(kSquared) + alphaSquared * ellint_1(kSquared));


      jacobianVectorHp[0][0] = ( nearZAxis == true) ?  0 : ((C * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * rhoSquared * mu0)) * (
                                 (a * a * a * a * (- gamma * (3 * xM[2] * xM[2] + a * a) + rhoSquared * (8 * xM[0] * xM[0] - xM[1] * xM[1])) -
                                  a * a * (rhoSquared * rhoSquared * (5 * xM[0] * xM[0] + xM[1] * xM[1]) -
                                           2 * rhoSquared * xM[2] * xM[2] * (2 * xM[0] * xM[0] + xM[1] * xM[1]) + 3 * xM[2] * xM[2] * xM[2] * xM[2] * gamma) -
                                  rSquared * rSquared * (2 * xM[0] * xM[0] * xM[0] * xM[0] + gamma * (xM[1] * xM[1] + xM[2] * xM[2])))  * ellint_2(kSquared) +
                                 (a * a * (gamma * (a * a + 2 * xM[2] * xM[2]) - rhoSquared * (3 * xM[0] * xM[0] - 2 * xM[1] * xM[1])) +
                                  rSquared * (2 * xM[0] * xM[0] * xM[0] * xM[0] + gamma * (xM[1] * xM[1] + xM[2] * xM[2]))) * alphaSquared * ellint_1(kSquared));


      jacobianVectorHp[0][1] = ( nearZAxis == true) ?  0 : ((C * xM[0] * xM[1] * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * rhoSquared * mu0)) * (
                                 (3 * a * a * a * a * (3 * rhoSquared - 2 * xM[2] * xM[2]) - rSquared * rSquared * (2 * rSquared + rhoSquared) - 2 * a * a * a * a * a * a  -
                                  2 * a * a * (2 * rhoSquared * rhoSquared - rhoSquared * xM[2] * xM[2] + 3 * xM[2] * xM[2] * xM[2] * xM[2])) * ellint_2(kSquared) +
                                 (rSquared * (2 * rSquared + rhoSquared) - a * a * (5 * rhoSquared - 4 * xM[2] * xM[2]) + 2 * a * a * a * a) * alphaSquared * ellint_1(kSquared));


      jacobianVectorHp[0][2] = ( nearZAxis == true) ?  0 : ((C * xM[0]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * mu0)) * (
                                 ((rhoSquared - a * a) * (rhoSquared - a * a) * (rhoSquared + a * a) + 2 * xM[2] * xM[2] * (a * a * a * a - 6 * a * a * rhoSquared + rhoSquared * rhoSquared) +
                                  xM[2] * xM[2] * xM[2] * xM[2] * (a * a + rhoSquared)) *  ellint_2(kSquared) -
                                 ((rhoSquared - a * a) * (rhoSquared - a * a) + xM[2] * xM[2] * (rhoSquared + a * a)) * alphaSquared *  ellint_1(kSquared));


      jacobianVectorHp[1][0] = jacobianVectorHp[0][1];


      jacobianVectorHp[1][1] = ( nearZAxis == true) ?  0 : ((C * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * rhoSquared * mu0)) * (
                                 (a * a * a * a * (gamma * (3 * xM[2] * xM[2] + a * a) + rhoSquared * (8 * xM[1] * xM[1] - xM[0] * xM[0])) -
                                  a * a * (rhoSquared * rhoSquared * (5 * xM[1] * xM[1] + xM[0] * xM[0]) - 2 * rhoSquared * xM[2] * xM[2] * (2 * xM[1] * xM[1] + xM[0] * xM[0]) -
                                           3 * xM[2] * xM[2] * xM[2] * xM[2] * gamma)  - rSquared * rSquared * (2 * xM[1] * xM[1] * xM[1] * xM[1] - gamma * (xM[0] * xM[0] + xM[2] * xM[2]))) * ellint_2(kSquared) +
                                 (a * a * (- gamma * (a * a + 2 * xM[2] * xM[2])  - rhoSquared * (3 * xM[1] * xM[1] - 2 * xM[0] * xM[0])) + rSquared * (2 * xM[1] * xM[1] * xM[1] * xM[1] -
                                     gamma * (xM[0] * xM[0] + xM[2] * xM[2]))) * alphaSquared * ellint_1(kSquared));


      jacobianVectorHp[1][2] = ( nearZAxis == true) ?  0 : ((C * xM[1]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * rhoSquared * mu0)) * (
                                 ((rhoSquared - a * a) * (rhoSquared - a * a) * (rhoSquared + a * a) + 2 * xM[2] * xM[2] * (a * a * a * a - 6 * a * a * rhoSquared + rhoSquared * rhoSquared) +
                                  xM[2] * xM[2] * xM[2] * xM[2] * (a * a + rhoSquared)) *  ellint_2(kSquared) -
                                 ((rhoSquared - a * a) * (rhoSquared - a * a) + xM[2] * xM[2] * (rhoSquared + a * a)) * alphaSquared *  ellint_1(kSquared));


      jacobianVectorHp[2][0] = ( nearZAxis == true) ?  0 : jacobianVectorHp[0][2];


      jacobianVectorHp[2][1] = ( nearZAxis == true) ?  0 : jacobianVectorHp[1][2];


      jacobianVectorHp[2][2] = ((C * xM[2]) / (2 * alphaSquared * alphaSquared * sqrt(betaSquared * betaSquared * betaSquared) * mu0)) * (
                                 (6 * a * a * (rhoSquared - xM[2] * xM[2]) - 7 * a * a * a * a + rSquared * rSquared) * ellint_2(kSquared) +
                                 (a * a - rSquared) * alphaSquared * ellint_1(kSquared));


// // 	  std::cout<<std::endl;
// //           for(unsigned i=0;i<3;i++){
// // 	    for(unsigned j=0;j<3;j++){
// // 	      std::cout << jacobianVectorHp[i][j] <<  " ";
// // 	    }
// // 	    std::cout<<std::endl;
// // 	  }
// // 	  std::cout<<std::endl;
// 	  for(unsigned i=0;i<3;i++){
// 	    std::cout << vectorHp[i] <<  " ";
// 	  }
//           std::cout<<std::endl;  std::cout<<std::endl;

      std::vector < std::vector <double> > jacobianVectorH(3);
      for(unsigned i = 0; i < 3; i++) {
        jacobianVectorH[i].resize(3);
      }

      for(unsigned i = 0; i < 3; i++) {
        vectorH[i] = 0.;
        for(unsigned j = 0; j < 3; j++) {
          vectorH[i] += R[i][j] * vectorHp[j];
          jacobianVectorH[i][j] = 0.;
          // std::cout << R[i][j] <<" ";
          // std::cout << jacobianVectorHp[i][j] <<" ";
          for(unsigned k = 0; k < 3; k++) {
            jacobianVectorH[i][j] += R[i][k] * jacobianVectorHp[k][j];
          }
        }
      }


      for(unsigned i = 0; i < 3; i++) {
        for(unsigned j = 0; j < 3; j++) {
          jacobianVectorHp[i][j] = 0. ;
          for(unsigned k = 0; k < 3; k++) {
            jacobianVectorHp[i][j] += jacobianVectorH[i][k] * R[j][k];
          }
        }
      }

      jacobianVectorH = jacobianVectorHp;
      
      
// 	  std::cout<<std::endl;
//           for(unsigned i=0;i<3;i++){
// 	    for(unsigned j=0;j<3;j++){
// 	      std::cout << jacobianVectorH[i][j] <<  " ";
// 	    }
// 	    std::cout<<std::endl;
// 	  }
// 	  std::cout<<std::endl;
// 	  for(unsigned i=0;i<3;i++){
// 	    std::cout << vectorH[i] <<  " ";
// 	  }
//           std::cout<<std::endl;  std::cout<<std::endl;

// 	  std::cout<<std::endl;
//           for(unsigned i=0;i<3;i++){
// 	    for(unsigned j=0;j<3;j++){
// 	      std::cout << jacobianVectorH[i][j] <<  " ";
// 	    }
// 	    std::cout<<std::endl;
// 	  }
// 	  std::cout<<std::endl;
// 	  for(unsigned i=0;i<3;i++){
// 	    std::cout << vectorH[i] <<  " ";
// 	  }
//           std::cout<<std::endl;  std::cout<<std::endl;


      H = 0.;
      for(unsigned  i = 0; i < 3; i++) {
        H += vectorH[i] * vectorH[i];
      }
      H = sqrt(H);

      for(unsigned  i = 0; i < 3; i++) {
        gradHSquared[i] = 0.;
        for(unsigned  j = 0; j < 3; j++) {
          gradHSquared[i] += 2 * vectorH[j] * jacobianVectorH[j][i];
        }
      }


//           gradHSquared[0] = 2 * vectorH[0] * jacobianVectorH[0][0] + 2 * vectorH[1] * jacobianVectorH[1][0] + 2 * vectorH[2] * jacobianVectorH[2][0];
//           gradHSquared[1] = 2 * vectorH[0] * jacobianVectorH[0][1] + 2 * vectorH[1] * jacobianVectorH[1][1] + 2 * vectorH[2] * jacobianVectorH[2][1];
//           gradHSquared[2] = 2 * vectorH[0] * jacobianVectorH[0][2] + 2 * vectorH[1] * jacobianVectorH[1][2] + 2 * vectorH[2] * jacobianVectorH[2][2];

      gradH[0] = 0.5 * (1 / H) * gradHSquared[0];
      gradH[1] = 0.5 * (1 / H) * gradHSquared[1];
      gradH[2] = 0.5 * (1 / H) * gradHSquared[2];

//           for(unsigned i=0;i<3;i++){
// 	    std::cout <<  gradHSquared[i]<<" "<< gradH[i] <<  " ";
// 	  }
//           std::cout<<std::endl;  std::cout<<std::endl;

    }

    break;

  }
  //END valuate H


  //BEGIN evaluate Fm

  if(H < H0) {


    double C1 = (PI * D * D * D * mu0 * chi) / 12;
    // printf("%g\n",C1);
    Fm[0] = C1 * gradHSquared[0];
    Fm[1] = C1 * gradHSquared[1];
    Fm[2] = C1 * gradHSquared[2];
  }

  else {

    double C2 = (PI * D * D * D * mu0 * Msat) / 6;

    //printf("%g\n",C2);
    Fm[0] = C2 * gradH[0];
    Fm[1] = C2 * gradH[1];
    Fm[2] = C2 * gradH[2];

  }


  for(unsigned i = 0 ; i < Fm.size(); i++) {
    Fm[i] = Fm[i] / (3 * PI * D * muf) ;
  }

  //END


  //BEGIN cheating to have attractive force

  for(unsigned i = 0 ; i < Fm.size(); i++) {
    Fm[i] = - Fm[i] ;
    //printf("%g ",Fm[i]);
  }


  //END cheating

}
