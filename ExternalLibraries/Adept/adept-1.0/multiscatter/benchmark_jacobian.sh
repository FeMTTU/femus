#!/bin/bash
# This script runs each combination of the algorithm to produce speed benchmarks

FILE="benchmark_profile50_AD.in"
PVC="../bin/multiscatter -algorithms fast none -repeat 75000 -quiet"
TDTS="../bin/multiscatter -algorithms none radar -repeat 1100 -quiet"


echo "*** PVC ORIGINAL ALGORITHM"
time $PVC $FILE > /dev/null
echo "*** PVC ADEPT JACOBIAN"
time $PVC -adjoint -adept -analytic-jacobian $FILE > /dev/null
#echo "*** PVC ADOLC JACOBIAN"
#time $PVC -adjoint -adolc -analytic-jacobian $FILE > /dev/null
#echo "*** PVC CPPAD JACOBIAN"
#time $PVC -adjoint -cppad -analytic-jacobian $FILE > /dev/null


echo "*** TDTS ORIGINAL ALGORITHM"
time $TDTS $FILE > /dev/null
echo "*** TDTS ADEPT JACOBIAN"
time $TDTS -adjoint -adept -analytic-jacobian $FILE > /dev/null
#echo "*** TDTS ADOLC JACOBIAN"
#time $TDTS -adjoint -adolc -analytic-jacobian $FILE > /dev/null
#echo "*** TDTS CPPAD JACOBIAN"
#time $TDTS -adjoint -cppad -analytic-jacobian $FILE > /dev/null



