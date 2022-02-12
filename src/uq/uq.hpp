
#ifndef __uq_hpp__
#define __uq_hpp__

#include <vector>
#include <map>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "UqQuadratureTypeEnum.hpp"

namespace femus {

  class uq {

    public:
      uq (UqQuadratureType quadratureType) :
        _quadratureType (quadratureType)
      {};

      ~uq() {
        Clear();
      };

      ////////////////////////////////////////////

      /// Get quadrature point coordinates
      const double* GetQuadraturePoints (const unsigned &numberOfQuadraturePoints);

      /// Get quadrature weights
      const double* GetQuadratureWeights (const unsigned &numberOfQuadraturePoints);

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

      /// Return the prescribed polynomial at the key < numberOfQuadraturePoints , maxPolyOrder >
      const std::vector < std::vector <double> > & GetPolynomial (const unsigned & numberOfQuadraturePoints,
                                                                  const unsigned & maxPolyOrder);

      /// Erase the prescribed polynomial at the key < numberOfQuadraturePoints , maxPolyOrder >, if stored
      void ErasePolynomial (const unsigned & numberOfQuadraturePoints, const unsigned & maxPolyOrder);

      /// Clear all stored  polynomials
      void ClearPolynomial ();

      /////////////////////////////////////////////////

      /// Compute the prescribed Polynomial Histogram at the key < pIndex, samplePoints, numberOfEigPairs >
      const std::vector < std::vector < double > >  & GetPolynomialHistogram (
        const unsigned & pIndex, const std::vector<double> & samplePoints, const unsigned & numberOfEigPairs);

      void ClearPolynomialHistogram() {
        _polynomialHistogram.clear();
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

      /// Compute the Multivariate prescribed polynomials and weights at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
      void ComputeMultivariate (
        std::vector < std::vector < double > >  & multivariatePoly,
        std::vector < double > & multivariateQuadratureWeights,
        const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs);

      /// Return the Multivariate prescribed polynomials at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
      const std::vector < std::vector < double > >  & GetMultivariatePolynomial (
        const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs);

      /// Return the Multivariate prescribed weight at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
      const std::vector < double > & GetMultivariateWeights (
        const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs);

      /// Erase the Multivariate prescribed polynomials and weights at the key < numberOfQuadraturePoints, p,numberOfEigPairs>
      void EraseMultivariate (const unsigned & numberOfQuadraturePoints, const unsigned & p, const unsigned & numberOfEigPairs);

      /// Clear all Multivariate polynomials and weights
      void ClearMultivariate ();

      /////////////////////////////////////////////////

      /// Set the standard output for some of the objects;
      void SetOutput (const bool & output);

      /////////////////////////////////////////////////

      void Clear() {
        ClearIndexSet();
        ClearIntegralMatrix();
        ClearStochasticMassMatrix();
        ClearMultivariate();
        ClearPolynomialHistogram();
        ClearPolynomial();
        ClearTensorProductSet();
      }

      //////////////////////////////////////////////////

      const UqQuadratureType & GetQuadratureType() {
        return _quadratureType;
      }

      ////////////////////////////////////////////

    private:
        
      /// Get Hermite quadrature point coordinates
      const double* GetHermiteQuadraturePoints (const unsigned &numberOfQuadraturePoints);

      /// Get Hermite quadrature weights
      const double* GetHermiteQuadratureWeights (const unsigned &numberOfQuadraturePoints);  
      
      /// Get Legenndre quadrature point coordinates
      const double* GetLegendreQuadraturePoints (const unsigned &numberOfQuadraturePoints);

      /// Get Legendrew quadrature weights
      const double* GetLegendreQuadratureWeights (const unsigned &numberOfQuadraturePoints);  
      
        
      /// Compute the Hermite Polynomial at the key < numberOfQuadraturePoints, maxPolyOrder >
      void ComputeHermitePolynomial (std::vector < std::vector < double > >  & hermitePoly,
                                     const unsigned &numberOfQuadraturePoints, const unsigned & maxPolyOrder);

      /// Compute the Legendre Polynomial at the key < numberOfQuadraturePoints, maxPolyOrder >
      void ComputeLegendrePolynomial (std::vector < std::vector < double > >  & legendrePoly,
                                      const unsigned &numberOfQuadraturePoints, const unsigned & maxPolyOrder);
      /// Compute the Hermite Polynomial Histogram at the key < pIndex, samplePoints, numberOfEigPairs >
      const std::vector < std::vector < double > > & GetHermitePolynomialHistogram (const unsigned & pIndex, const std::vector<double> & samplePoints, const unsigned & numberOfEigPairs);

      /// Compute the Legendre Polynomial Histogram at the key < pIndex, samplePoints, numberOfEigPairs >
      const std::vector < std::vector < double > > & GetLegendrePolynomialHistogram (const unsigned & pIndex, const std::vector<double> & samplePoints, const unsigned & numberOfEigPairs);

      ////////////////////////////////////////////

    private:
      const UqQuadratureType _quadratureType;
      static const double _hermiteQuadrature[16][2][16];
      static const double _legendreQuadrature[16][2][16];
      std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <unsigned> > > _Jp;
      std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <unsigned> > > _Tp;
      std::map<std::pair<unsigned, unsigned>, std::vector < std::vector <double> > > _polynomial;
      std::vector < std::vector < double > >  _polynomialHistogram;
      std::map<std::pair<unsigned, unsigned>, std::vector < std::vector < std::vector < double > > > > _integralMatrix;
      std::map<std::pair < std::pair<unsigned, unsigned>, unsigned >, std::vector < std::vector < std::vector < double > > > > _stochasticMassMatrix;
      std::map<std::pair < std::pair<unsigned, unsigned>, unsigned >, std::vector < std::vector < double > > >  _multivariatePolynomial;
      std::map<std::pair < std::pair<unsigned, unsigned>, unsigned >, std::vector < double > > _multivariateWeight;
      bool _output = false;
  };

}

#endif
