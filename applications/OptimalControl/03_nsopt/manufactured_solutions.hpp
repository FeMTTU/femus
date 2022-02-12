#ifndef MANUFACTURED_SOLUTIONS_NSOPT
#define MANUFACTURED_SOLUTIONS_NSOPT


namespace mms_lid_driven {
//ns full dir - BEGIN ******************************************

//manufactured solution for lid-driven---------------------------------------------

double value_statePress(const std::vector < double >& x) {
  double pi = acos(-1.);
  return sin(2. * pi * x[0]) * sin(2. * pi * x[1]); //p
 };

 void gradient_statePress(const std::vector < double >& x, std::vector < double >& grad_statePress) {
  double pi = acos(-1.);
  grad_statePress[0]  =   2. * pi * cos(2. * pi * x[0]) * sin(2. * pi * x[1]); 
  grad_statePress[1]  =   2. * pi * sin(2. * pi * x[0]) * cos(2. * pi * x[1]);
 };

//manufactured solution for lid-driven---------------------------------------------



//ns full dir - END ******************************************




//opt bdry ctrl - BEGIN ******************************************

//adj---------------------------------------------
void value_adjVel(const std::vector < double >& x, std::vector < double >& val_adjVel) {
  double pi = acos(-1.);
  val_adjVel[0] =   0.5 * sin(pi* x[0]) * sin(pi* x[0]) *  sin(2. * pi * x[1]); //u
  val_adjVel[1] = - 0.5 * sin(2. * pi * x[0]) * sin(pi* x[1]) * sin(pi* x[1]); //v
 };
 
double value_adjPress(const std::vector < double >& x) { 
  double pi = acos(-1.);
  return sin(2. * pi * x[0]) * sin(2. * pi * x[1]); //p
 };
 
 
void gradient_adjVel(const std::vector < double >& x, std::vector < std::vector < double > >& grad_adjVel) {
  double pi = acos(-1.);
  grad_adjVel[0][0]  =   0.5 * pi * sin(2. * pi * x[0]) * sin(2. * pi * x[1]); 
  grad_adjVel[0][1]  =   pi * sin(pi* x[0]) * sin(pi* x[0]) *  cos(2. * pi * x[1]);
  grad_adjVel[1][0]  = - pi * cos(2. * pi * x[0]) * sin(pi * x[1]) * sin(pi * x[1]); 
  grad_adjVel[1][1]  = - 0.5 * pi * sin(2. * pi * x[0]) * sin(2. * pi * x[1]);
 };

 void gradient_adjPress(const std::vector < double >& x, std::vector < double >& grad_adjPress) {
  double pi = acos(-1.);
  grad_adjPress[0]  =   2. * pi * cos(2. * pi * x[0]) * sin(2. * pi * x[1]); 
  grad_adjPress[1]  =   2. * pi * sin(2. * pi * x[0]) * cos(2. * pi * x[1]);
 };
 
 
void laplace_adjVel(const std::vector < double >& x, std::vector < double >& lap_adjVel) {
  double pi = acos(-1.);
  lap_adjVel[0] = pi * pi * cos(2. * pi * x[0]) * sin(2. * pi * x[1]) - 2. * pi * pi * sin(pi* x[0]) * sin(pi* x[0]) *  sin(2. * pi * x[1]);
  lap_adjVel[1] = 2. * pi * pi * sin(2. * pi * x[0]) * sin(pi* x[1]) * sin(pi* x[1]) - pi * pi * sin(2. * pi * x[0]) * cos(2. * pi * x[1]);
};
//adj---------------------------------------------

//state---------------------------------------------
void value_stateVel(const std::vector < double >& x, std::vector < double >& val_stateVel) {
  double pi = acos(-1.);
  val_stateVel[0] =   sin(pi* x[0]) * sin(pi* x[0]) * cos(pi* x[1]) - sin(pi* x[0]) * sin(pi* x[0]);
  val_stateVel[1] = - sin(2. * pi * x[0]) * sin(pi* x[1]) + pi * x[1] * sin(2. * pi * x[0]);
 };
 
 
void gradient_stateVel(const std::vector < double >& x, std::vector < std::vector < double > >& grad_stateVel) {
  double pi = acos(-1.);
  grad_stateVel[0][0]  =   pi * sin(2. * pi * x[0]) * cos(pi* x[1]) - pi * sin(2. * pi * x[0]);
  grad_stateVel[0][1]  = - pi * sin(pi* x[0]) * sin(pi* x[0]) *  sin(pi * x[1]); 
  grad_stateVel[1][0]  = - 2. * pi * cos(2. * pi * x[0]) * sin(pi* x[1]) + 2. * pi * pi * x[1] * cos(2. * pi * x[0]);   
  grad_stateVel[1][1]  = - pi * sin(2. * pi * x[0]) * cos(pi * x[1]) + pi * sin(2. * pi * x[0]); 
 };

  
void laplace_stateVel(const std::vector < double >& x, std::vector < double >& lap_stateVel) {
  double pi = acos(-1.);
  lap_stateVel[0] = - 2. * pi * pi * cos(2. * pi * x[0]) - 0.5 * pi * pi * cos(pi * x[1]) + 2.5 * pi * pi * cos(2. * pi* x[0]) * cos(pi* x[1]);
  lap_stateVel[1] = - 4. * pi * pi * pi * x[1] * sin(2. * pi * x[0]) + 5. * pi * pi * sin(2. * pi * x[0]) * sin(pi * x[1]);
};
//state---------------------------------------------



//opt bdry ctrl - END ******************************************

} //end lid driven




namespace mms_state_control {
//opt lift - BEGIN ******************************************
    
//state---------------------------------------------
void value_stateVel(const std::vector < double >& x, std::vector < double >& val_stateVel) {
  double pi = acos(-1.);
  val_stateVel[0] =   0.5 * sin(pi* x[0]) * sin(pi* x[0]) *  sin(2. * pi * x[1]); //u
  val_stateVel[1] = - 0.5 * sin(2. * pi * x[0]) * sin(pi* x[1]) * sin(pi* x[1]); //v
 };
 
double value_statePress(const std::vector < double >& x) {
  double pi = acos(-1.);
  return sin(2. * pi * x[0]) * sin(2. * pi * x[1]); //p
 };
 
 
void gradient_stateVel(const std::vector < double >& x, std::vector < std::vector < double > >& grad_stateVel) {
  double pi = acos(-1.);
  grad_stateVel[0][0]  =   0.5 * pi * sin(2. * pi * x[0]) * sin(2. * pi * x[1]); 
  grad_stateVel[0][1]  =   pi * sin(pi* x[0]) * sin(pi* x[0]) *  cos(2. * pi * x[1]);
  grad_stateVel[1][0]  = - pi * cos(2. * pi * x[0]) * sin(pi * x[1]) * sin(pi * x[1]); 
  grad_stateVel[1][1]  = - 0.5 * pi * sin(2. * pi * x[0]) * sin(2. * pi * x[1]);
 };

 void gradient_statePress(const std::vector < double >& x, std::vector < double >& grad_statePress) {
  double pi = acos(-1.);
  grad_statePress[0]  =   2. * pi * cos(2. * pi * x[0]) * sin(2. * pi * x[1]); 
  grad_statePress[1]  =   2. * pi * sin(2. * pi * x[0]) * cos(2. * pi * x[1]);
 };
 
 
void laplace_stateVel(const std::vector < double >& x, std::vector < double >& lap_stateVel) {
  double pi = acos(-1.);
  lap_stateVel[0] = pi * pi * cos(2. * pi * x[0]) * sin(2. * pi * x[1]) - 2. * pi * pi * sin(pi* x[0]) * sin(pi* x[0]) *  sin(2. * pi * x[1]);
  lap_stateVel[1] = 2. * pi * pi * sin(2. * pi * x[0]) * sin(pi* x[1]) * sin(pi* x[1]) - pi * pi * sin(2. * pi * x[0]) * cos(2. * pi * x[1]);
};
//state---------------------------------------------


//control---------------------------------------------
void value_ctrlVel(const std::vector < double >& x, std::vector < double >& val_ctrlVel) {
  double pi = acos(-1.);
  val_ctrlVel[0] =   sin(pi* x[0]) * sin(pi* x[0]) * cos(pi* x[1]) - sin(pi* x[0]) * sin(pi* x[0]);
  val_ctrlVel[1] = - sin(2. * pi * x[0]) * sin(pi* x[1]) + pi * x[1] * sin(2. * pi * x[0]);
 };
 
 
void gradient_ctrlVel(const std::vector < double >& x, std::vector < std::vector < double > >& grad_ctrlVel) {
  double pi = acos(-1.);
  grad_ctrlVel[0][0]  =   pi * sin(2. * pi * x[0]) * cos(pi* x[1]) - pi * sin(2. * pi * x[0]);
  grad_ctrlVel[0][1]  = - pi * sin(pi* x[0]) * sin(pi* x[0]) *  sin(pi * x[1]); 
  grad_ctrlVel[1][0]  = - 2. * pi * cos(2. * pi * x[0]) * sin(pi* x[1]) + 2. * pi * pi * x[1] * cos(2. * pi * x[0]);   
  grad_ctrlVel[1][1]  = - pi * sin(2. * pi * x[0]) * cos(pi * x[1]) + pi * sin(2. * pi * x[0]); 
 };

  
void laplace_ctrlVel(const std::vector < double >& x, std::vector < double >& lap_ctrlVel) {
  double pi = acos(-1.);
  lap_ctrlVel[0] = - 2. * pi * pi * cos(2. * pi * x[0]) - 0.5 * pi * pi * cos(pi * x[1]) + 2.5 * pi * pi * cos(2. * pi* x[0]) * cos(pi* x[1]);
  lap_ctrlVel[1] = - 4. * pi * pi * pi * x[1] * sin(2. * pi * x[0]) + 5. * pi * pi * sin(2. * pi * x[0]) * sin(pi * x[1]);
};
//control---------------------------------------------


//opt lift - END ******************************************
}



#endif
