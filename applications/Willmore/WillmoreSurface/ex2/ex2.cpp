/** tutorial/Ex3
 * This example shows how to set and solve the weak form of the nonlinear problem
 *                     -\Delta^2 u = f(x) \text{ on }\Omega,
 *            u=0 \text{ on } \Gamma,
 *      \Delat u=0 \text{ on } \Gamma,
 * on a box domain $\Omega$ with boundary $\Gamma$,
 * by using a system of second order partial differential equation.
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "adept.h"
#include <cstdlib>


using namespace femus;

double a = sqrt(2);

// Torus

bool SetBoundaryConditionTorus(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  double u = x[0];
  double v = x[1];

  if (!strcmp("X", SolName)) {
    value = (a + cos(u)) * cos(v);
  }
  else if (!strcmp("Y", SolName)) {
    value = (a + cos(u)) * sin(v);
  }
  else if (!strcmp("Z", SolName)) {
    value = sin(u);
  }
  else if (!strcmp("H", SolName)) {
    value = 0.5 * (1. + cos(u) / (a + cos(u)));
  }

  return dirichlet;
}

double InitalValueXTorus(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return (a + cos(u)) * cos(v);

}

double InitalValueYTorus(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return (a + cos(u)) * sin(v);

}

double InitalValueZTorus(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return sin(u);

}

double InitalValueHTorus(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return 0.5 * (1. + cos(u) / (a + cos(u)));

}


bool SetBoundaryConditionSphere(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  double u = x[0];
  double v = x[1];

  if (!strcmp("X", SolName)) {
    value = a * sin(v) * cos(u);
  }
  else if (!strcmp("Y", SolName)) {
    value = a * sin(v) * sin(u);
  }
  else if (!strcmp("Z", SolName)) {
    value = a * cos(v);
  }
  else if (!strcmp("H", SolName)) {
    value = 1. / a;
  }

  return dirichlet;
}

double InitalValueXSphere(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return a * sin(v) * cos(u);

}

double InitalValueYSphere(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return  a * sin(v) * sin(u);

}

double InitalValueZSphere(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return a * cos(v);

}

double InitalValueHSphere(const std::vector < double >& x) {
  double u = x[0];
  double v = x[1];

  return 1. / a;

}




double SetVariableTimeStep(const double time) {
  double dt = 1.;
  return dt;
}


void AssembleWillmoreFlow_AD(MultiLevelProblem& ml_prob);

std::pair < double, double > GetErrorNorm(MultiLevelSolution* mlSol);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);


  // define multilevel mesh


  unsigned maxNumberOfMeshes;
 
  MultiLevelMesh mlMsh;
  // read coarse level mesh and generate finers level meshes
  double scalingFactor = 1.;
  
  mlMsh.ReadCoarseMesh("./input/sphere.neu", "seventh", scalingFactor);
  
  unsigned numberOfUniformLevels = 4;
  unsigned numberOfSelectiveLevels = 0;
  mlMsh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

  // erase all the coarse mesh levels
  // mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);

  // print mesh info
  mlMsh.PrintInfo();

  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);

    // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, SECOND, 0);

  mlSol.Initialize("All");
 
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);


  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
  return 0;
}



