#!/bin/bash
# This script runs each combination of the algorithm to produce speed benchmarks

FILE="benchmark_profile100_AD.in"
PVC="../bin/multiscatter -algorithms fast none -repeat 160 -quiet"

echo "*** PVC ORIGINAL ALGORITHM"
time $PVC $FILE > /dev/null
echo "*** PVC ADEPT ADJOINT"
time $PVC -adjoint -adept $FILE > /dev/null
echo "*** PVC SACADO ADJOINT"
time $PVC -adjoint -sacado $FILE > /dev/null
echo "*** PVC ADOLC ADJOINT"
time $PVC -adjoint -adolc $FILE > /dev/null
echo "*** PVC CPPAD ADJOINT"
time $PVC -adjoint -cppad $FILE > /dev/null

return

echo "*** PVC ADEPT JACOBIAN"
time $PVC -adjoint -adept -analytic-jacobian $FILE > /dev/null
echo "*** PVC ADOLC JACOBIAN"
time $PVC -adjoint -adolc -analytic-jacobian $FILE > /dev/null
echo "*** PVC CPPAD JACOBIAN"
time $PVC -adjoint -cppad -analytic-jacobian $FILE > /dev/null





