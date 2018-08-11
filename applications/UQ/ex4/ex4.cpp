#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

#include "slepceps.h"

#include "../include/sgfem_assembly.hpp"

using namespace femus;

void GetEigenPair(MultiLevelProblem& ml_prob, const int& numberOfEigPairs, std::vector < std::pair<double, double> >& eigenvalues);
//
void GetCoefficientsForQuantityOfInterest(MultiLevelProblem& ml_prob, std::vector <double > &  alphas, const double& domainMeasure);
//
void GetStochasticData(std::vector <double>& alphas);
//
void PlotStochasticData();

//BEGIN stochastic data

unsigned totMoments = 7;
std::vector <double> moments(totMoments, 0.); //initialization
std::vector <double> momentsStandardized(totMoments, 0.); //initialization
std::vector <double> cumulants(totMoments, 0.); //initialization
std::vector <double> cumulantsStandardized(totMoments, 0.); //initialization
double meanQoI = 0.; //initialization
double varianceQoI = 0.; //initialization
double stdDeviationQoI = 0.; //initialization
unsigned M = 1000000;
double startPoint = 0.3678 + 0.051253826666667 * 0.5;
double endPoint = 2.7;

//FOR NORMAL DISTRIBUTION
boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)
boost::normal_distribution<> nd(0., 0.3);
boost::variate_generator < boost::mt19937&,
      boost::normal_distribution<> > var_nor(rng, nd);

//FOR UNIFORM DISTRIBUTION
boost::mt19937 rng1; // I don't seed it on purpouse (it's not relevant)
boost::random::uniform_real_distribution<> un(- 1., 1.);
boost::variate_generator < boost::mt19937&,
      boost::random::uniform_real_distribution<> > var_unif(rng1, un);

//FOR LAPLACE DISTRIBUTION
boost::mt19937 rng2; // I don't seed it on purpouse (it's not relevant)
boost::random::uniform_real_distribution<> un1(- 0.5, 0.49999999999);
boost::variate_generator < boost::mt19937&,
      boost::random::uniform_real_distribution<> > var_unif1(rng2, un1);
double b = 2.;
//END

int main(int argc, char** argv) {

  std::vector<double> QoI(M, 0.);

  for(unsigned m = 0; m < M; m++) {

    double var = var_nor();
    double varunif = var_unif();
    double U = var_unif1();
//     QoI[m] = var * var * var;
//     QoI[m] = exp(var);
    QoI[m] = exp(varunif);

    //exp of truncated gaussian
//     if(fabs(var) <= 1.) {
//       QoI[m] = var / (0.5 * ((1. + erf((1. / 0.3) / sqrt(2))) - (1. + erf((- 1. / 0.3) / sqrt(2))))) ;    //truncated Gaussian
//     }
//     else QoI[m] = 0.;

    //laplace distribution
//     double signU = 0.;
//     if(U < 0) signU = - 1.;
//     else if(U > 0) signU = 1.;
//     QoI[m] = 0. - b * signU * log(1. - 2. * fabs(U)) ;

//     std::cout << "QoI[" << m << "]=" << QoI[m] << std::endl;

  }

  GetStochasticData(QoI);

  PlotStochasticData();

  return 0;

} //end main



//
void GetStochasticData(std::vector <double>& QoI) {

  //let's standardize the quantity of interest after finding moments and standard deviation

  if(totMoments <= 0) {

    std::cout << "ERROR: total number of moments has to be a positive integer" << std::endl;

  }

  else {

    //BEGIN computation of the mean of QoI
    meanQoI = 0.;
    unsigned meanCounter = 0;
    for(unsigned m = 0; m < M; m++) {
      if(QoI[m] < endPoint && QoI[m] > startPoint) {
        meanQoI += QoI[m];
        meanCounter++;
      }
    }
    meanQoI /= meanCounter;
    //END


    //BEGIN computation of the variance and standard deviation of QoI
    varianceQoI = 0;
    unsigned varianceCounter = 0;
    for(unsigned m = 0; m < M; m++) {
      if(QoI[m] < endPoint && QoI[m] > startPoint) {
        varianceQoI += (QoI[m] - meanQoI) * (QoI[m] - meanQoI);
        varianceCounter++;
      }
    }
    varianceQoI /= varianceCounter;

    stdDeviationQoI = sqrt(varianceQoI);
    //END

    int pdfHistogramSize = static_cast <int>(1. + 3.3 * log(M));
    std::vector <double> pdfHistogram(pdfHistogramSize, 0.);
    double lengthOfTheInterval = fabs(endPoint - startPoint);
    double deltat = lengthOfTheInterval / (pdfHistogramSize - 1);

    std::vector < double > QoIStandardized(M, 0.);
    //BEGIN standardization of QoI before computing the moments
    for(unsigned m = 0; m < M; m++) {
//       QoIStandardized[m] = (QoI[m] - meanQoI) / stdDeviationQoI ;
         QoIStandardized[m] = QoI[m] ;
//       std::cout << "standardized QoI " << QoIStandardized[m] << std::endl;

      //BEGIN estimation of the PDF
      bool sampleCaptured = false;
      for(unsigned i = 0; i < pdfHistogramSize; i++) {
        double leftBound = startPoint + i * deltat - deltat * 0.5;
        double rightBound = startPoint + i * deltat + deltat * 0.5;
        if(leftBound <=  QoIStandardized[m] && QoIStandardized[m] < rightBound) {
          pdfHistogram[i]++;
//           std::cout << "leftBound = " << leftBound << " " << "rightBound = " << rightBound << " " << " standardized QoI = " << QoIStandardized[m] << std::endl;
          sampleCaptured = true;
          break;
        }
      }
      if(sampleCaptured == false) {
        std::cout << "WARNING: sample " << QoIStandardized[m] << "is not in any interval" << std::endl;
      }
      //END
    }
    //END
    
    double pdfIntegral = 0;
    for(unsigned i = 0; i < pdfHistogramSize; i++) {
      pdfIntegral += pdfHistogram[i] * deltat;
    }


    //BEGIN histogram check
    double checkHistogram = 0;
    for(unsigned i = 0; i < pdfHistogramSize; i++) {
      double point = startPoint + i * deltat;
      double pdfCheck = pdfHistogram[i]/M;
      pdfHistogram[i] /= pdfIntegral;
      std::cout << point << "  " << pdfHistogram[i]  << std::endl;
//       std::cout << "{" << point << "," << pdfHistogram[i]  << "}," << std::endl;
      checkHistogram += pdfCheck;
    }
    std::cout << "checkHistogram = " << checkHistogram << std::endl;
    //END



    //BEGIN computation of the raw moments
    for(unsigned p = 0; p < totMoments; p++) {
      moments[p] = 0.;
      momentsStandardized[p] = 0.;
      unsigned momentsCounter = 0;
      for(unsigned m = 0; m < M; m++) {
        if(QoI[m] < endPoint && QoI[m] > startPoint) {
          moments[p] += pow(QoI[m], p + 1);
          momentsStandardized[p] += pow(QoIStandardized[m], p + 1);
          momentsCounter++;
        }
      }
      moments[p] /= momentsCounter;
      momentsStandardized[p] /= momentsCounter;
    }
    //END

    cumulants[0] = moments[0];
    cumulantsStandardized[0] = momentsStandardized[0];

    if(totMoments > 1) {
      cumulants[1] = moments[1] - moments[0] * moments[0];

      cumulantsStandardized[1] = momentsStandardized[1] - momentsStandardized[0] * momentsStandardized[0];
//       std::cout.precision(14);
//       std::cout << "AAAAAAAAAAAAAAA" << cumulants[1] << std::endl;
      if(totMoments > 2) {
        cumulants[2] = moments[2] - 3. * moments[1] * moments[0] + 2. * pow(moments[0], 3);

        cumulantsStandardized[2] = momentsStandardized[2] - 3. * momentsStandardized[1] * momentsStandardized[0] + 2. * pow(momentsStandardized[0], 3);
        if(totMoments > 3) {
          cumulants[3] = moments[3] - 4. * moments[2] * moments[0] - 3. * moments[1] * moments[1] + 12. * moments[1] * moments[0] * moments[0] - 6. * pow(moments[0], 4);

          cumulantsStandardized[3] = momentsStandardized[3] - 4. * momentsStandardized[2] * momentsStandardized[0] - 3. * momentsStandardized[1] * momentsStandardized[1]
                                     + 12. * momentsStandardized[1] * momentsStandardized[0] * momentsStandardized[0] - 6. * pow(momentsStandardized[0], 4);
          if(totMoments > 4) {
            cumulants[4] = moments[4] - 5. * moments[3] * moments[0] - 10. * moments[2] * moments[1] + 20. * moments[2] * moments[0] * moments[0]
                           + 30. * moments[1] * moments[1] * moments[0] - 60. * moments[1] * pow(moments[0], 3) + 24. * pow(moments[0], 5);

            cumulantsStandardized[4] = momentsStandardized[4] - 5. * momentsStandardized[3] * momentsStandardized[0] - 10. * momentsStandardized[2] * momentsStandardized[1]
                                       + 20. * momentsStandardized[2] * momentsStandardized[0] * momentsStandardized[0]
                                       + 30. * momentsStandardized[1] * momentsStandardized[1] * momentsStandardized[0]
                                       - 60. * momentsStandardized[1] * pow(momentsStandardized[0], 3) + 24. * pow(momentsStandardized[0], 5);
            if(totMoments > 5) {
              cumulants[5] = moments[5] - 6. * moments[4] * moments[0] - 15. * moments[3] * moments[1] + 30. * moments[3] * moments[0] * moments[0]
                             - 10. * moments[2] * moments[2] + 120. * moments[2] * moments[1] * moments[0] - 120. * moments[2] * pow(moments[0], 3)
                             + 30. * pow(moments[1], 3) - 270. * pow(moments[1], 2) * pow(moments[0], 2) + 360. * moments[1] * pow(moments[0], 4) - 120. * pow(moments[0], 6);

              cumulantsStandardized[5] = momentsStandardized[5] - 6. * momentsStandardized[4] * momentsStandardized[0] - 15. * momentsStandardized[3] * momentsStandardized[1]
                                         + 30. * momentsStandardized[3] * momentsStandardized[0] * momentsStandardized[0]
                                         - 10. * momentsStandardized[2] * momentsStandardized[2] + 120. * momentsStandardized[2] * momentsStandardized[1] * momentsStandardized[0]
                                         - 120. * momentsStandardized[2] * pow(momentsStandardized[0], 3) + 30. * pow(momentsStandardized[1], 3)
                                         - 270. * pow(momentsStandardized[1], 2) * pow(momentsStandardized[0], 2) + 360. * momentsStandardized[1] * pow(momentsStandardized[0], 4)
                                         - 120. * pow(momentsStandardized[0], 6);

              if(totMoments > 6) {
                cumulants[6] = 720. * pow(moments[0], 7) - 2520. * moments[1] * pow(moments[0], 5) + 840. * moments[2] * pow(moments[0], 4)
                               + 2520. * pow(moments[1], 2) * pow(moments[0], 3) - 210. * moments[3] * pow(moments[0], 3) - 1260. * moments[1] * moments[2] * pow(moments[0], 2)
                               + 42. * moments[4] * pow(moments[0], 2) - 630. * pow(moments[1], 3) * moments[0] + 140. * pow(moments[2], 2) * moments[0]
                               + 210. * moments[1] * moments[3] * moments[0] - 7. * moments[5] * moments[0] + 210. * pow(moments[1], 2) * moments[2]
                               - 35. * moments[2] * moments[3] - 21. * moments[1] * moments[4] + moments[6];

                cumulantsStandardized[6] = 720. * pow(momentsStandardized[0], 7) - 2520. * momentsStandardized[1] * pow(momentsStandardized[0], 5) + 840. * momentsStandardized[2] * pow(momentsStandardized[0], 4)
                                           + 2520. * pow(momentsStandardized[1], 2) * pow(momentsStandardized[0], 3) - 210. * momentsStandardized[3] * pow(momentsStandardized[0], 3) - 1260. * momentsStandardized[1] * momentsStandardized[2] * pow(momentsStandardized[0], 2)
                                           + 42. * momentsStandardized[4] * pow(momentsStandardized[0], 2) - 630. * pow(momentsStandardized[1], 3) * momentsStandardized[0] + 140. * pow(momentsStandardized[2], 2) * momentsStandardized[0]
                                           + 210. * momentsStandardized[1] * momentsStandardized[3] * momentsStandardized[0] - 7. * momentsStandardized[5] * momentsStandardized[0] + 210. * pow(momentsStandardized[1], 2) * momentsStandardized[2]
                                           - 35. * momentsStandardized[2] * momentsStandardized[3] - 21. * momentsStandardized[1] * momentsStandardized[4] + momentsStandardized[6];
              }
            }
          }
        }
      }
    }
  }
}
//
//
void PlotStochasticData() {

  std::cout.precision(10);
  std::cout << " the mean is " << meanQoI << std::endl;
  std::cout << " the standard deviation is " << stdDeviationQoI << std::endl;

  std::cout << "Standardized Moments" << std::endl;
  for(unsigned p = 0; p < totMoments; p++) {
    std::cout << " & " << momentsStandardized[p] << "  ";
  }
  std::cout << std::endl;
  std::cout << "Standardized Cumulants" << std::endl;
  for(unsigned p = 0; p < totMoments; p++) {
    std::cout << " & " << cumulantsStandardized[p] << "  ";
  }
  std::cout << std::endl;
  std::cout << " Moments " << std::endl;
  for(unsigned p = 0; p < totMoments; p++) {
    std::cout << " & " << moments[p] << "  ";
  }
  std::cout << std::endl;
  std::cout << " Cumulants " << std::endl;
  for(unsigned p = 0; p < totMoments; p++) {
    std::cout << " & " << cumulants[p] << "  ";
  }

  std::cout << std::endl;
  std::cout << " --------------------------------------------------------------------------------------------- " << std::endl;

  double edgeworth1Term = 0.;
  double edgeworth2Terms = 0.;
  double edgeworth3Terms = 0.;
  double edgeworth4Terms = 0.;
  double edgeworth5Terms = 0.;
  double edgeworth6Terms = 0.;

  double generalizedGC1Term = 0.;
  double generalizedGC2Terms = 0.;
  double generalizedGC3Terms = 0.;
  double generalizedGC4Terms = 0.;
  double generalizedGC5Terms = 0.;
  double generalizedGC6Terms = 0.;
  double generalizedGC7Terms = 0.;

  double lambda3 = 0.;
  double lambda4 = 0.;
  double lambda5 = 0.;
  double lambda6 = 0.;

  double d1gaussian;
  double d2gaussian;
  double d3gaussian;
  double d4gaussian;
  double d5gaussian;
  double d6gaussian;
  double d7gaussian;
  double d8gaussian;
  double d9gaussian;
  double d10gaussian;
  double d12gaussian;

  double t = startPoint;
  double dt = (fabs(endPoint - startPoint)) / 500.;

//   cumulants[0] = 0; //decomment for nonStdGaussian

  //BEGIN GRAM CHARLIER PRINT
  std::cout << " ------------------------- GRAM CHARLIER ------------------------- " << std::endl;
  for(unsigned i = 0; i <= 500; i++) {
    std::cout << t << " ";
//     double t = x - meanQoI; //decomment for nonStdGaussian
    double gaussian = 1. / (sqrt(2 * acos(- 1))) * exp(- 0.5 * (t * t)) ;
    std::cout << gaussian << " ";

    d1gaussian = (- 1.) * gaussian * t ;

    generalizedGC1Term = gaussian - cumulantsStandardized[0] * d1gaussian;

    std::cout << generalizedGC1Term << " ";

    if(totMoments > 1) {

      d2gaussian = (1.) * gaussian * (t * t - 1.) ;
      d3gaussian = (- 1.) * gaussian * (t * t * t - 3. * t) ;

      generalizedGC2Terms = generalizedGC1Term + 0.5 * ((cumulantsStandardized[1] - 1.) + pow(cumulantsStandardized[0], 2)) * d2gaussian ;

      std::cout << generalizedGC2Terms << " ";

      if(totMoments > 2) {

        d4gaussian = (1.) * gaussian * (t * t * t * t - 6. * t * t + 3.) ;
        d6gaussian = (1.) * gaussian * (pow(t, 6) - 15 * pow(t, 4) + 45 * t * t - 15);

        generalizedGC3Terms = generalizedGC2Terms - 1. / 6. * (cumulantsStandardized[2] + 3 * (cumulantsStandardized[1] - 1.) * cumulantsStandardized[0]
                              + pow(cumulantsStandardized[0], 3)) * d3gaussian;

        std::cout << generalizedGC3Terms << " ";

        if(totMoments > 3) {

          d5gaussian = (- 1.) * gaussian * (pow(t, 5) - 10. * t * t * t + 15. * t);
          d7gaussian = (- 1.) * gaussian * (pow(t, 7) - 21. * pow(t, 5) + 105. * t * t * t -  105. * t) ;
          d9gaussian = (- 1.) * gaussian * (pow(t, 9) - 36. * pow(t, 7) + 378. * pow(t, 5) - 1260. * t * t * t + 945. * t) ;

          generalizedGC4Terms = generalizedGC3Terms + 1. / 24. * (cumulantsStandardized[3] + 4. * cumulantsStandardized[2] * cumulantsStandardized[0]
                                + 3. * pow((cumulantsStandardized[1] - 1.), 2) + 6. * (cumulantsStandardized[1] - 1.) * pow(cumulantsStandardized[0], 2)
                                + pow(cumulantsStandardized[0], 4)) * d4gaussian;

          std::cout << generalizedGC4Terms << " ";

          if(totMoments > 4) {

            generalizedGC5Terms = generalizedGC4Terms - 1. / 120. * (cumulantsStandardized[4] + 5. * cumulantsStandardized[3] * cumulantsStandardized[0]
                                  + 10. * cumulantsStandardized[2] * (cumulantsStandardized[1] - 1.) + 10. * cumulantsStandardized[2] * pow(cumulantsStandardized[0], 2)
                                  + 15. * pow((cumulantsStandardized[1] - 1.), 2) * cumulantsStandardized[0] + 10. * (cumulantsStandardized[1] - 1.) * pow(cumulantsStandardized[0], 3)
                                  + pow(cumulantsStandardized[0], 5)) * d5gaussian;

            std::cout << generalizedGC5Terms << " ";

            if(totMoments > 5) {

              generalizedGC6Terms = generalizedGC5Terms + 1. / 720. * (cumulantsStandardized[5] + 6. * cumulantsStandardized[4] * cumulantsStandardized[0]
                                    + 15. * cumulantsStandardized[3] * (cumulantsStandardized[1] - 1.) + 15. * cumulantsStandardized[3] * pow(cumulantsStandardized[0], 2)
                                    + 10. * pow(cumulantsStandardized[2], 2) + 60. * cumulantsStandardized[2] * (cumulantsStandardized[1] - 1.) * cumulantsStandardized[0]
                                    + 20. * cumulantsStandardized[2] * pow(cumulantsStandardized[0], 3) + 15. * pow((cumulantsStandardized[1] - 1.), 3)
                                    + 45. * pow((cumulantsStandardized[1] - 1.), 2) * pow(cumulantsStandardized[0], 2) + 15. * (cumulantsStandardized[1] - 1.) * pow(cumulantsStandardized[0], 4)
                                    +  pow(cumulantsStandardized[0], 6)) * d6gaussian;

              std::cout << generalizedGC6Terms << "  ";

              if(totMoments > 6) {

                generalizedGC7Terms = generalizedGC6Terms - 1. / 5040. * (pow(cumulantsStandardized[0], 7) + 21. * (cumulantsStandardized[1] - 1.) * pow(cumulantsStandardized[0], 5)
                                      + 35. * cumulantsStandardized[2] * pow(cumulantsStandardized[0], 4) + 105. * (cumulantsStandardized[1] - 1.) * (cumulantsStandardized[1] - 1.) * pow(cumulantsStandardized[0], 3)
                                      + 35. * cumulantsStandardized[3] * pow(cumulantsStandardized[0], 3) + 210. * (cumulantsStandardized[1] - 1.) * cumulantsStandardized[2] * pow(cumulantsStandardized[0], 2)
                                      + 21. * cumulantsStandardized[4] * pow(cumulantsStandardized[0], 2) + 105. * pow(cumulantsStandardized[1] - 1., 3) * cumulantsStandardized[0]
                                      + 70. * pow(cumulantsStandardized[2], 2) * cumulantsStandardized[0] + 105. * (cumulantsStandardized[1] - 1.) * cumulantsStandardized[3] * cumulantsStandardized[0]
                                      + 7. * cumulantsStandardized[5] * cumulantsStandardized[0] + 105. * pow(cumulantsStandardized[1] - 1., 2) * cumulantsStandardized[2]
                                      + 35. * cumulantsStandardized[2] * cumulantsStandardized[3] + 21. * (cumulantsStandardized[1] - 1.) * cumulantsStandardized[4]
                                      + cumulantsStandardized[6]) * d7gaussian;

                std::cout << generalizedGC7Terms << " \n ";

              }
            }
          }
        }
      }
    }

    t += dt;
  }

  t = startPoint;
  dt = (fabs(endPoint - startPoint)) / 500.;

  //BEGIN EDGEWORTH PRINT
  std::cout << " ------------------------- EDGEWORTH ------------------------- " << std::endl;
  for(unsigned i = 0; i <= 500; i++) {
    std::cout << t << " ";
//     double t = x - meanQoI; //decomment for nonStdGaussian
    double gaussian = 1. / (sqrt(2 * acos(- 1))) * exp(- 0.5 * (t * t)) ;
    std::cout << gaussian << " ";

    d3gaussian = (- 1.) * gaussian * (t * t * t - 3. * t) ;

    lambda3 = cumulants[2] / pow(stdDeviationQoI, 3);

    edgeworth1Term = gaussian - lambda3 / 6. * d3gaussian;

    std::cout << edgeworth1Term << " ";

    if(totMoments > 1) {

      d4gaussian = (1.) * gaussian * (t * t * t * t - 6. * t * t + 3.) ;
      d6gaussian = (1.) * gaussian * (pow(t, 6) - 15 * pow(t, 4) + 45 * t * t - 15);

      lambda4 = cumulants[3] / pow(stdDeviationQoI, 4);

      edgeworth2Terms = edgeworth1Term + lambda4 / 24. * d4gaussian + lambda3 * lambda3 / 72. * d6gaussian;

      std::cout << edgeworth2Terms << " ";

      if(totMoments > 2) {

        d5gaussian = (- 1.) * gaussian * (pow(t, 5) - 10. * t * t * t + 15. * t);
        d7gaussian = (- 1.) * gaussian * (pow(t, 7) - 21. * pow(t, 5) + 105. * t * t * t -  105. * t) ;
        d9gaussian = (- 1.) * gaussian * (pow(t, 9) - 36. * pow(t, 7) + 378. * pow(t, 5) - 1260. * t * t * t + 945. * t) ;

        lambda5 = cumulants[4] / pow(stdDeviationQoI, 5);

        edgeworth3Terms = edgeworth2Terms - lambda5 / 120. * d5gaussian + lambda3 * lambda4 / 144. * d7gaussian + pow(lambda3, 3) / 1296. * d9gaussian;

        std::cout << edgeworth3Terms << " ";

      }
      if(totMoments > 3) {

        d8gaussian = (1.) * gaussian * (1. / 16. * (1680. - 6720. * t * t + 3360. * pow(t, 4) - 448. * pow(t, 6) + 16. * pow(t, 8)));
        d10gaussian = (1.) * gaussian * (1. / 32. * (- 30240. + 151200. *  t * t - 100800. * pow(t, 4) + 20160. * pow(t, 6) - 1440. * pow(t, 8)
                                         + 32. * pow(t, 10)));
        d12gaussian = (1.) * gaussian * (1. / 64. * (665280. - 3991680. * t * t + 3326400. * pow(t, 4) - 887040. * pow(t, 6) + 95040. * pow(t, 8)
                                         - 4224. * pow(t, 10) + 64. * pow(t, 12)));

        lambda6 = cumulants[5] / pow(stdDeviationQoI, 6);

        edgeworth4Terms = edgeworth3Terms + 1. / 720. * lambda6 * d6gaussian + (1. / 1152. * lambda4 * lambda4 + 1. / 720. * lambda3 * lambda5) * d8gaussian
                          + 1. / 1728. * lambda3 * lambda3 * lambda4 * d10gaussian + 1. / 31104. * pow(lambda3, 4) * d12gaussian;

        std::cout << edgeworth4Terms << " \n ";

      }
    }

    t += dt;
  }
  //END


}
//










