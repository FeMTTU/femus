#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>

static const double um = 1.0;
static const double H = 0.40;
static const double L = 1.5;

// Inlet velocity
// extern "C" double InitialValueU(const std::vector < double >& xyz) {
//   const double x = xyz[0];
//   const double y = xyz[1];
//   const double z = xyz[2];
//   const double k = 9./(H*H*H*H);
//   return k*(z*(H - z))*(H*H - y*y)*um*exp(-10.*L*x);
// }

extern "C" double TimeStepFunction(const double time) {
  double dt;
  if(time < 3.0) {
    dt = 0.1;
  } else 
  {
    dt = 0.01; 
  }
  return dt;
}


//---------------------------------------------------------------------------------------------------------------------

extern "C" bool BdcFunction(const std::vector < double >& xyz,const char name[], double &value, const int facename, const double time) {
  bool test=1; //dirichlet
  value=0.;
  if(!strcmp(name,"U")) {
    // inflow
    if(5==facename){
      test=1; // Dirichlet
      const double k = 9./(H*H*H*H);
      const double x = xyz[0];
      const double y = xyz[1];
      const double z = xyz[2];
      if(time < 2.0) {
        value = k*(z*(H - z))*(H*H - y*y)*um* 0.5 * ( 1. - cos( 0.5 * 3.141592653589793 * time) );
      }
      else {
        value = k*(z*(H - z))*(H*H - y*y)*um;
      }
    }
    // fluid wall
    else if(1==facename){
      test=1; // Dirichlet
      value=0.;
    }
    // solid base
    else if(2==facename ){
      test=1; // Dirichlet
      value=0.;
    }
    // slip - fluid symmetric
    else if(3==facename ){
      test=0;
      value=0.;
    }
    // slip - solid symmetric
    else if(4==facename ){   // beam case zero stress
      test=0; // Neumann
      value=0.;
    }
    // outflow
    else if(6==facename ){
      test=0; // Neumann
      value=0.;
    }
  }
  else if(!strcmp(name,"V")) {
    // inflow
    if(5==facename){
      test=1; // Dirichlet
      value=0;
    }
    // fluid wall
    else if(1==facename){
      test=1; // Dirichlet
      value=0.;
    }
    // // solid base
    else if(2==facename ){
      test=1; // Dirichlet
      value=0.;
    }
    // slip - fluid symmetric
    else if(3==facename ){
      test=1; // Dirichlet
      value=0.;
    }
    // slip - solid symmetric
    else if(4==facename ){   // beam case zero stress
      test=1; // Dirichlet
      value=0.;
    }
    // outflow
    else if(6==facename ){
      test=0; // Neumann
      value=0.;
    }
  }
  else if(!strcmp(name,"W")) {
    // inflow
    if(5==facename){
      test=1; // Dirichlet
      value=0;
    }
    // fluid wall
    else if(1==facename){
      test=1; // Dirichlet
      value=0.;
    }
    // solid base
    else if(2==facename ){
      test=1; // Dirichlet
      value=0.;
    }
    // slip - fluid symmetric
    else if(3==facename ){
      test=0; // Dirichlet
      value=0.;
    }
    // slip - solid symmetric
    else if(4==facename ){   // beam case zero stress
      test=0; // Neumann
      value=0.;
    }
    // outflow
    else if(6==facename ){
      test=0; // Neumann
      value=0.;
    }
  }
  else if(!strcmp(name,"P")) {
    // inflow
    if(5==facename){
      test=0; // Neumann
      value=0;
    }
    // fluid wall
    else if(1==facename){
      test=0; // Neumann
      value=0.;
    }
    // solid base
    else if(2==facename ){
      test=0; // Neumann
      value=0.;
    }
    // slip - fluid symmetric
    else if(3==facename ){
      test=0; // Neumann
      value=0.;
    }
    // slip - solid symmetric
    else if(4==facename ){   // beam case zero stress
      test=0; // Neumann
      value=0.;
    }
    // outflow
    else if(6==facename ){
      test=0; // Neumann
      value=0.;
    }
  }
  else if(!strcmp(name,"DX")) {
    // inflow
    if(5==facename){
      test=1; // Dirichlet
      value=0;
    }
    // fluid wall
    else if(1==facename){
      test=1; // Dirichlet
      value=0.;
    }
    // solid base
    else if(2==facename ){
      test=1; // Dirichlet
      value=0.;
    }
    // slip - fluid symmetric
    else if(3==facename ){
      test=0; // Neumann
      value=0.;
    }
    // slip - solid symmetric
    else if(4==facename ){   // beam case zero stress
      test=0; // Neumann
      value=0.;
    }
    // outflow
    else if(6==facename ){
      test=1; // Dirichlet
      value=0.;
    }
  }
  else if(!strcmp(name,"DY")) {
    // inflow
    if(5==facename){
      test=1; // Dirichlet
      value=0;
    }
    // fluid wall
    else if(1==facename){
      test=1; // Dirichlet
      value=0.;
    }
    // solid base
    else if(2==facename ){
      test=1; // Dirichlet
      value=0.;
    }
    // slip - fluid symmetric
    else if(3==facename ){
      test=1; // Dirichlet
      value=0.;
    }
    // slip - solid symmetric
    else if(4==facename ){   // beam case zero stress
      test=1; // Dirichlet
      value=0.;
    }
    // outflow
    else if(6==facename ){
      test=1; // Dirichlet
      value=0.;
    }
  }
  else if(!strcmp(name,"DZ")) {
    // inflow
    if(5==facename){
      test=1; // Dirichlet
      value=0;
    }
    // fluid wall
    else if(1==facename){
      test=1; // Dirichlet
      value=0.;
    }
    // solid base
    else if(2==facename ){
      test=1; // Dirichlet
      value=0.;
    }
    // slip - fluid symmetric
    else if(3==facename ){
      test=0; // Neumann
      value=0.;
    }
    // slip - solid symmetric
    else if(4==facename ){   // beam case zero stress
      test=0; // Neumann
      value=0.;
    }
    // outflow
    else if(6==facename ){
      test=1; // Dirichlet
      value=0.;
    }
  }
  return test;
}