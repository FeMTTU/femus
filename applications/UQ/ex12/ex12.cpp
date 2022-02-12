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
#include "sparseGrid.hpp"

//THIS EXTENDS WHAT IS IN EX10 TO SPARSE GRIDS (so stochastic dimension grater than 3)

using namespace femus;



//BEGIN stochastic data

unsigned alpha = 6;
unsigned M = pow ( 10, alpha ); //number of samples
unsigned N = 2; //dimension of the parameter space (each of the M samples has N entries)
unsigned L = 6; //max refinement level
bool output = false; //for debugging
bool matlabView = true;
bool histoView = false;

double xmin = - 5.5;   //-1.5 for uniform // -5.5 for Gaussian
double xmax = 5.5;     //1.5 for uniform // 5.5 for Gaussian

//FOR NORMAL DISTRIBUTION
boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)
boost::normal_distribution<> nd ( 0., 1. );
boost::variate_generator < boost::mt19937&,
      boost::normal_distribution<> > var_nor ( rng, nd );

//FOR UNIFORM DISTRIBUTION
boost::mt19937 rng1; // I don't seed it on purpouse (it's not relevant)
boost::random::uniform_real_distribution<> un ( - 1., 1. );
boost::variate_generator < boost::mt19937&,
      boost::random::uniform_real_distribution<> > var_unif ( rng1, un );

//FOR LAPLACE DISTRIBUTION
boost::mt19937 rng2; // I don't seed it on purpouse (it's not relevant)
boost::random::uniform_real_distribution<> un1 ( - 0.5, 0.49999999999 );
boost::variate_generator < boost::mt19937&,
      boost::random::uniform_real_distribution<> > var_unif1 ( rng2, un1 );
double b = 2.;
//END

int main ( int argc, char** argv )
{

    //BEGIN construction of the sample set
    std::vector < std::vector < double > >  samples;
    samples.resize ( M );

    for ( unsigned m = 0; m < M; m++ ) {
        samples[m].resize ( N );
    }

    for ( unsigned m = 0; m < M; m++ ) {
        for ( unsigned n = 0; n < N; n++ ) {

            double var = var_nor();
            double varunif = var_unif();
            double U = var_unif1();
//     samples[m][n] = var * var * var;
//     samples[m][n] = exp(var);
//         samples[m][n] = exp (varunif);
            samples[m][n] = var;

            //exp of truncated gaussian
//     if(fabs(var) <= 1.) {
//       samples[m][n] = var / (0.5 * ((1. + erf((1. / 0.3) / sqrt(2))) - (1. + erf((- 1. / 0.3) / sqrt(2))))) ;    //truncated Gaussian
//     }
//     else samples[m][n] = 0.;

            //laplace distribution
//     double signU = 0.;
//     if(U < 0) signU = - 1.;
//     else if(U > 0) signU = 1.;
//     samples[m][n] = 0. - b * signU * log(1. - 2. * fabs(U)) ;

//     std::cout << "samples[" << m << "][" << n << "]=" << samples[m][n] << std::endl;

        }
    }

    //END



    //BEGIN initialize grid and compute nodal values
    clock_t total_time = clock();
    clock_t grid_time = clock();
    sparseGrid spg ( samples, xmin, xmax, output );

    std::cout << std::endl << " Builds sparse grid in: " << std::setw ( 11 ) << std::setprecision ( 6 ) << std::fixed
              << static_cast<double> ( ( clock() - grid_time ) ) / CLOCKS_PER_SEC << " s" << std::endl;

    clock_t nodal_time = clock();
    spg.EvaluateNodalValuesPDF ( samples );
//     spg.PrintNodalValuesPDF();

    std::cout << std::endl << " Builds nodal values in: " << std::setw ( 11 ) << std::setprecision ( 6 ) << std::fixed
              << static_cast<double> ( ( clock() - nodal_time ) ) / CLOCKS_PER_SEC << " s" << std::endl;
    //END



    //BEGIN  create grid for plot in 2D

    std::vector < unsigned > refinementLevel ( N );

    refinementLevel[0] = L; //refinement level in x

    if ( N > 1 )  refinementLevel[1] = L; //refinement level in y

    if ( N > 2 )  refinementLevel[2] = L; //refinement level in x

    std::vector < unsigned > gridPoints ( N );
    std::vector < std::vector < double> > gridBounds ( N );

    for ( unsigned n = 0; n < N; n++ ) {
        gridPoints[n] = static_cast<unsigned> ( pow ( 2, refinementLevel[n] ) + 1 );
        gridBounds[n].resize ( 2 );
    }

    unsigned gridSize = 1;

    for ( unsigned n = 0; n < N; n++ ) {
        gridSize *= gridPoints[n];
    }

    gridBounds[0][0] = xmin;
    gridBounds[0][1] = xmax;

    if ( N > 1 ) {

        gridBounds[1][0] = xmin;
        gridBounds[1][1] = xmax;
    }

    if ( N > 2 ) {

        gridBounds[2][0] = xmin;
        gridBounds[2][1] = xmax;
    }

    std::vector < double > h ( N );

    for ( unsigned n = 0; n < N; n++ ) {
        h[n] = ( gridBounds[n][1] - gridBounds[n][0] ) / pow ( 2, refinementLevel[n] );
    }

    std::vector < std::vector < double > > grid;

    unsigned counterGrid = 0;

    if ( N == 1 ) {
        for ( unsigned i = 0; i < gridPoints[0]; i++ ) {
            grid.resize ( counterGrid + 1 );
            grid[counterGrid].resize ( N );
            grid[counterGrid][0] = gridBounds[0][0] + i * h[0];
            counterGrid++;
        }
    }

    else if ( N > 1 ) {
        for ( unsigned j = 0; j < gridPoints[1]; j++ ) {
            for ( unsigned i = 0; i < gridPoints[0]; i++ ) {
                grid.resize ( counterGrid + 1 );
                grid[counterGrid].resize ( N );
                grid[counterGrid][0] = gridBounds[0][0] + i * h[0];
                grid[counterGrid][1] = gridBounds[1][0] + j * h[1];
                counterGrid++;
            }
        }
    }

    //END create grid


    //BEGIN create histogram on finest grid for comparison

    unsigned histoSize1D =  static_cast<unsigned> ( pow ( 2, refinementLevel[0] ) ) ;
    unsigned totNumberOfBins = static_cast<unsigned> ( pow ( histoSize1D, N ) ); //tot number of bins
    std::vector <double> histogram ( totNumberOfBins );

    std::vector < std::vector < double > > histoGrid ( totNumberOfBins );

    if ( histoView ) {

        for ( unsigned i = 0; i < histoGrid.size(); i++ ) {
            histoGrid[i].resize ( N );
        }

        unsigned counterHisto = 0;

        if ( N == 1 ) {
            for ( unsigned i = 0; i < histoSize1D; i++ ) {
                histoGrid[counterHisto][0] = ( gridBounds[0][0] + h[0] * 0.5 ) + i * h[0];
                counterHisto++;
            }
        }

        else if ( N > 1 ) {
            for ( unsigned j = 0; j < histoSize1D; j++ ) {
                for ( unsigned i = 0; i < histoSize1D; i++ ) {
                    histoGrid[counterHisto][0] = ( xmin + h[0] * 0.5 ) + i * h[0];
                    histoGrid[counterHisto][1] = ( xmin + h[1] * 0.5 ) + j * h[1];
                    counterHisto++;
                }
            }
        }

        for ( unsigned m = 0; m < M; m++ ) {
            for ( unsigned i = 0; i < histoGrid.size(); i++ ) {

                std::vector<unsigned> inSupport ( N, 0 );

                for ( unsigned n = 0; n < N; n++ ) {
                    if ( samples[m][n] > ( histoGrid[i][n] - 0.5 * h[n] ) && samples[m][n] <= ( histoGrid[i][n] + 0.5 * h[n] ) ) {
                        inSupport[n] = 1;
                    }
                }

                unsigned sumCheck = 0;

                for ( unsigned n = 0; n < N; n++ ) {
                    sumCheck += inSupport[n];
                }

                if ( sumCheck == N ) {
                    histogram[i]++;
                    break;
                }

            }
        }

        //evaluation of histogram integral

        double histoIntegral = 0.;
        double supportMeasure = pow ( h[0] , N );

        for ( unsigned i = 0; i < histoGrid.size(); i++ ) {
            histoIntegral += supportMeasure * histogram[i];
        }

        for ( unsigned i = 0; i < histoGrid.size(); i++ ) {
            histogram[i] /= histoIntegral;
        }
    }

//END



//BEGIN grid plot
    if ( matlabView ) {
        std::cout << "x=[" << std::endl;

        for ( unsigned i = 0; i < grid.size(); i++ ) {
            std::cout << grid[i][0] << std::endl;
        }

        std::cout << "];" << std::endl;

        if ( N > 1 ) {
            std::cout << "y=[" << std::endl;

            for ( unsigned i = 0; i < grid.size(); i++ ) {
                std::cout << grid[i][1] << std::endl;
            }

            std::cout << "];" << std::endl;
        }

        clock_t pdf_time = clock();

        std::cout << "PDF=[" << std::endl;

        for ( unsigned i = 0; i < grid.size(); i++ ) {
            double pdfValue;
            spg.EvaluatePDF ( pdfValue, grid[i], true );
        }

        std::cout << "];" << std::endl;

        std::cout << std::endl << " Builds PDF in: " << std::setw ( 11 ) << std::setprecision ( 6 ) << std::fixed
                  << static_cast<double> ( ( clock() - pdf_time ) ) / CLOCKS_PER_SEC << " s" << std::endl;



        if ( histoView ) {
            std::cout << "HISTO=[" << std::endl;

            for ( unsigned i = 0; i < grid.size(); i++ ) {

                bool sampleCought  = false;

                for ( unsigned j = 0; j < histoGrid.size(); j++ ) {

                    std::vector <unsigned> inBin ( N, 0 );

                    for ( unsigned n = 0; n < N; n++ ) {

                        double leftBound = histoGrid[j][n] - h[n] * 0.5;

                        double rightBound = histoGrid[j][n] + h[n] * 0.5;

                        if ( grid[i][n] >  leftBound && grid[i][n] <= rightBound ) {

                            inBin[n] = 1;

                        }

                    }

                    unsigned sumCheck = 0;

                    for ( unsigned n = 0; n < N; n++ ) {
                        sumCheck += inBin[n];
                    }

                    if ( sumCheck == N ) {

                        for ( unsigned n = 0; n < N; n++ ) {
//                     std::cout << "grid[" << i << "][" << n << "] = " << grid[i][n] << " " << histoGrid[j][n] - h[n] * 0.5 << " , " <<  histoGrid[j][n] + h[n] * 0.5  << std::endl;
                        }

                        std::cout << histogram[j] << std::endl;

                        sampleCought = true;

                        break;

                    }

                }

                if ( sampleCought == false ) std::cout << 0 << std::endl;

            }

            std::cout << "];" << std::endl;
        }

    }

//END  plot

    double integral;
    spg.EvaluatePDFIntegral ( integral );

    std::cout << " PDF Integral is = " << integral << std::endl;


//BEGIN compute error
    clock_t error_time = clock();

    double aL2E;
    spg.ComputeAvgL2Error ( aL2E, samples, 1 );

    std::cout << " Averaged L2 error is = " << aL2E << std::endl;

    std::cout << " Computes error in: " << std::setw ( 11 ) << std::setprecision ( 6 ) << std::fixed
              << static_cast<double> ( ( clock() - error_time ) ) / CLOCKS_PER_SEC << " s" << std::endl;
//END
    std::cout << " Total time: " << std::setw ( 11 ) << std::setprecision ( 6 ) << std::fixed
              << static_cast<double> ( ( clock() - total_time ) ) / CLOCKS_PER_SEC << " s" << std::endl;

              
              //BEGIN check sparse grid and histo
//     for ( unsigned i = 0; i < histoSize1D; i++ ) {
//         for ( unsigned j = 0; j < histoSize1D; j++ ) {
// 
//             unsigned k = j * histoSize1D + i;
// 
//             std::cout << histogram[k] << " ";
//         }
// 
//         std::cout << std::endl;
//     }
/*    
    std::cout<<"------------------------------------------------------"<<std::endl;

    for ( unsigned i = 0; i < histoSize1D; i++ ) {
        for ( unsigned j = 0; j < histoSize1D; j++ ) {
            
            unsigned k = j * histoSize1D + i;

            double pdfValue;
            spg.EvaluatePDF ( pdfValue, histoGrid[k], false );
            
            std::cout << pdfValue << " ";

        }
        
           std::cout << std::endl;
    }*/
//END

    return 0;

} //end main



