#ifndef __femus_enums_ImplicitRKEnum_hpp__
#define __femus_enums_ImplicitRKEnum_hpp__

enum ImplicitRKScheme {
  LEGENDRE1 = 0,
  LEGENDRE2,
  LEGENDRE3,
  NORSET3
};


// Gauss Legendre 1 point. Order of convergenge 2
const double cLEGENDRE1[1] = {0.5};
const double bLEGENDRE1[1] =  {1.};
const double aLEGENDRE1[1] = {0.5};
const double aiLEGENDRE1[1] = {2.};

// Gauss Legendre 2 points. Order of convergenge 4
const double cLEGENDRE2[2] = {0.5 - sqrt (3.) / 6., 0.5 + sqrt (3) / 6.};
const double bLEGENDRE2[2] =  {0.5, 0.5};
const double aLEGENDRE2[4] = {
  0.25, 0.25 - sqrt (3.) / 6.,
  0.25 + sqrt (3.) / 6., 0.25
};
const double aiLEGENDRE2[4] = {
  3., 0.4641016151377544,
  -6.464101615137755, 3.
};

// Gauss Legendre 3 points. Order of convergenge 6
const double cLEGENDRE3[3] = {0.5 - sqrt (15.) / 10., 0.5, 0.5 + sqrt (15.) / 10.};
const double bLEGENDRE3[3] =  {5. / 18., 4. / 9., 5. / 18.};
const double aLEGENDRE3[9] = {
  5. / 36., 2. / 9. - sqrt (15.) / 15.,   5. / 36. - sqrt (15.) / 30.,
  5. / 36. + sqrt (15.) / 24.,   2. / 9., 5. / 36. - sqrt (15.) / 24.,
  5. / 36. + sqrt (15.) / 30.,   2. / 9. + sqrt (15.) / 15.,   5. / 36.
};
const double aiLEGENDRE3[9] = {
  5., 1.1639777949432233, -0.16397779494322232,
  -5.727486121839513, 2.,  0.7274861218395138,
  10.163977794943225, -9.163977794943223, 5.
};


// Norsett 3 points. Order of convergenge 4
const double cNORSET3[3] =  {1.0685790213016289, 0.5, -.06857902130162885};
const double bNORSET3[3] =    {0.1288864005157204, 0.7422271989685593, 0.1288864005157204};
const double aNORSET3[9] = {
  1.0685790213016289, 0., 0.,
  -0.5685790213016289, 1.0685790213016289, 0.,
  2.1371580426032577, -3.2743160852065154, 1.0685790213016289
};
const double aiNORSET3[9] = {
  0.9358222275240877, 0., 0.,
  0.497940606760015, 0.9358222275240878, 0.,
  -0.34586591580096865, 2.867525668568206, 0.935822227524088
};

const unsigned nRK[4] = {1,2,3,3};
const double *cIRK[4] ={ cLEGENDRE1,  cLEGENDRE2,  cLEGENDRE3,  cNORSET3};
const double *bIRK[4] ={ bLEGENDRE1,  bLEGENDRE2,  bLEGENDRE3,  bNORSET3};
const double *aIRK[4] ={ aLEGENDRE1,  aLEGENDRE2,  aLEGENDRE3,  aNORSET3};
const double *aiIRK[4]={ aiLEGENDRE1, aiLEGENDRE2, aiLEGENDRE3, aiNORSET3};











#endif
