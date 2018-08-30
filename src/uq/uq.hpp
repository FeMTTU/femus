



#ifndef __uq_hpp__
#define __uq_hpp__

#include <vector>
#include <map>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/binomial.hpp>

namespace femus {

  class uq {

    public:

      /// Get Hermite quadrature point coordinates
      const double* GetHermiteQuadraturePoints (const unsigned &numberOfQuadraturePoints);

      /// Get Hermite quadrature weights
      const double* GetHermiteQuadratureWeights (const unsigned &numberOfQuadraturePoints);

      ////////////////////////////////////////////

      /// Compute the tensor product at the key < numberOfQuadraturePoints,  numberOfEigPairs >
      void ComputeTensorProductSet (std::vector < std::vector <unsigned> > & Tp,
                                    const unsigned & numberOfQuadraturePoints,
                                    const unsigned & numberOfEigPairs);

      /// Return the tensor product set at the key < numberOfQuadraturePoints , numberOfEigPairs >
      const std::vector < std::vector <unsigned> > & GetTensorProductSet (const unsigned & numberOfQuadraturePoints,
                                                                          const unsigned & numberOfEigPairs);

      /// Erase the tensor product set at the key < numberOfQuadraturePoints , numberOfEigPairs >, if stored
      void EraseTensorProductSet (const unsigned & numberOfQuadraturePoints, const unsigned & numberOfEigPairs);

      /// Clear all stored tensor product sets
      void ClearTensorProductSet();

      ////////////////////////////////////////////

      /// Compute the Hermite Polynomial at the key < numberOfQuadraturePoints, maxPolyOrder >
      void ComputeHermitePoly (std::vector < std::vector < double > >  & hermitePoly,
                               const unsigned &numberOfQuadraturePoints, const unsigned & maxPolyOrder);


      /// Return the Hermite polynomial at the key < numberOfQuadraturePoints , maxPolyOrder >
      const std::vector < std::vector <double> > & GetHermitePolynomial (const unsigned & numberOfQuadraturePoints,
                                                                         const unsigned & maxPolyOrder);

      /// Erase the Hermite polynomial at the key < numberOfQuadraturePoints , maxPolyOrder >, if stored
      void EraseHermitePolynomial (const unsigned & numberOfQuadraturePoints, const unsigned & maxPolyOrder);

      /// Clear all stored Hermite polynomials
      void ClearHermitePolynomial();

      /////////////////////////////////////////////////

      /// Compute the Hermite Polynomial Histogram at the key < pIndex, samplePoints, numberOfEigPairs >
      const std::vector < std::vector < double > > & GetHermitePolyHistogram (const unsigned & pIndex, const std::vector<double> & samplePoints, const unsigned & numberOfEigPairs);


      void ClearHermitePolynomialHistogram() {
        _hermitePolyHistogram.clear();
      }
      /////////////////////////////////////////////////

      /// Compute the Index Set Jp at the key < p, numberOfEigPairs>
      void ComputeIndexSet (std::vector < std::vector <unsigned> > & Jp, const unsigned & p,
                            const unsigned & numberOfEigPairs);

      /// Return the Index Set at the key < p, numberOfEigPairs>
      const std::vector < std::vector <unsigned> > & GetIndexSet (const unsigned & p, const unsigned & numberOfEigPairs);

      /// Erase the Index Set at the key < p, numberOfEigPairs>, if stored
      void EraseIndexSet (const unsigned & numberOfQuadraturePoints, const unsigned & maxPolyOrder);

      /// Clear all Index Sets Jp
      void ClearIndexSet();

      /////////////////////////////////////////////////

      /// Compute the Integral Matrix at the key < q0, p0>
      void ComputeIntegralMatrix (std::vector < std::vector < std::vector < double > > > &integralMatrix,
                                  const unsigned & q0, const unsigned & p0);

      /// Get the Integral Matrix at the key < q0, p0>
      const std::vector < std::vector < std::vector < double > > > &GetIntegralMatrix (const unsigned & q0, const unsigned & p0);

      /// Erase the Integral Matrix at the key < q0, p0>
      void EraseIntegralMatrix (const unsigned & q0, const unsigned & p0);;

      /// Clear all Integral Matrices
      void ClearIntegralMatrix();

      /////////////////////////////////////////////////

      /// Compute the Integral Matrix at the key < q0, p0, numberOfEigPairs>
      void ComputeStochasticMassMatrix (std::vector < std::vector < std::vector < double > > > & G,
                                        const unsigned & q0, const unsigned & p0, const unsigned & numberOfEigPairs);

      /// Return the Stochastic Mass Matrix at the key < q0, p0, numberOfEigPairs>
      std::vector < std::vector < std::vector < double > > > & GetStochasticMassMatrix (const unsigned & q0, const unsigned & p0,
                                                                                        const unsigned & numberOfEigPairs);


      /// Erase the Stochastic Mass Matrix at the key < q0, p0, numberOfEigPairs>
      void EraseStochasticMassMatrix (const unsigned & q0, const unsigned & p0,
                                      const unsigned & numberOfEigPairs);

      /// Clear all stored Stochastic Mass Matrices
      void ClearStochasticMassMatrix();

      /////////////////////////////////////////////////

      /// Compute the Multivariate Hermite polynomials and weights at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
      void ComputeMultivariateHermite (std::vector < std::vector < double > >  & multivariateHermitePoly,
                                       std::vector < double > & multivariateHermiteQuadratureWeights,
                                       const unsigned & numberOfQuadraturePoints, const unsigned & p,
                                       const unsigned & numberOfEigPairs);

      /// Return the Multivariate Hermite polynomials at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
      const std::vector < std::vector < double > >  & GetMultivariateHermitePolynomial (
        const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs);

      /// Return the Multivariate Hermite weight at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
      const std::vector < double > & GetMultivariateHermiteWeights (
        const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs);

      /// Erase the Multivariate Hermite polynomials and weights at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
      void EraseMultivariateHermite (const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs);

      /// Clear all Multivariate Hermite polynomials and weights
      void ClearMultivariateHermite();

      /////////////////////////////////////////////////

      void Clear() {
        ClearIndexSet();
        ClearIntegralMatrix();
        ClearStochasticMassMatrix();
        ClearMultivariateHermite();
        ClearHermitePolynomialHistogram();
        ClearHermitePolynomial();
        ClearTensorProductSet();
      }

    private:
      static const double _hermiteQuadrature[16][2][16];
      std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <unsigned> > > _Jp;
      std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <unsigned> > > _Tp;
      std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <double> > > _hermitePoly;
      std::vector < std::vector < double > >  _hermitePolyHistogram;
      std::map<std::pair<unsigned, unsigned>, std::vector < std::vector < std::vector < double > > > > _integralMatrix;
      std::map<std::pair < std::pair<unsigned, unsigned>, unsigned >, std::vector < std::vector < std::vector < double > > > > _stochasticMassMatrix;
      std::map<std::pair < std::pair<unsigned, unsigned>, unsigned >, std::vector < std::vector < double > > >  _multivariateHermitePoly;
      std::map<std::pair < std::pair<unsigned, unsigned>, unsigned >, std::vector < double > > _multivariateHermiteWeight;
  };

}

#endif
