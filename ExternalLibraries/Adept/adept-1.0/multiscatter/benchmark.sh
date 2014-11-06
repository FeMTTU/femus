#!/bin/bash
# This script runs each combination of the algorithm to produce speed benchmarks

FILE="benchmark_profile50_AD.in"
PVC_ARGS="-algorithms fast none -repeat 75000 -quiet"
TDTS_ARGS="-algorithms none radar -repeat 1100 -quiet"

PVC="./multiscatter $PVC_ARGS"
TDTS="./multiscatter $TDTS_ARGS"

PVC_REDUCED="./multiscatter_reduced $PVC_ARGS"
TDTS_REDUCED="./multiscatter_reduced $TDTS_ARGS"

echoerr() { echo "$@" 1>&2; }
TIME=
#TIME=time

echoerr "*** PVC ORIGINAL ALGORITHM"
time $PVC_REDUCED $FILE > /dev/null
echoerr "*** PVC HAND-CODED ADJOINT"
$TIME $PVC_REDUCED -adjoint $FILE > /dev/null
echoerr "*** PVC ADEPT ADJOINT"
$TIME $PVC_REDUCED -adjoint -adept $FILE > /dev/null

echoerr "*** TDTS ORIGINAL ALGORITHM"
time $TDTS_REDUCED $FILE > /dev/null
echoerr "*** TDTS HAND-CODED ADJOINT"
$TIME $TDTS_REDUCED -adjoint $FILE > /dev/null
echoerr "*** TDTS ADEPT ADJOINT"
$TIME $TDTS_REDUCED -adjoint -adept $FILE > /dev/null

echoerr "*** PVC SACADO ADJOINT"
$TIME $PVC -adjoint -sacado $FILE > /dev/null
echoerr "*** PVC ADOLC ADJOINT"
$TIME $PVC -adjoint -adolc $FILE > /dev/null
echoerr "*** PVC CPPAD ADJOINT"
$TIME $PVC -adjoint -cppad $FILE > /dev/null

echoerr "*** TDTS SACADO ADJOINT"
$TIME $TDTS -adjoint -sacado $FILE > /dev/null
echoerr "*** TDTS ADOLC ADJOINT"
$TIME $TDTS -adjoint -adolc $FILE > /dev/null
echoerr "*** TDTS CPPAD ADJOINT"
$TIME $TDTS -adjoint -cppad $FILE > /dev/null

echoerr "*** PVC ADEPT JACOBIAN (DEFAULT)"
$TIME $PVC_REDUCED -adjoint -adept -analytic-jacobian $FILE > /dev/null
echoerr "*** PVC ADEPT JACOBIAN (FORWARD)"
$TIME $PVC_REDUCED -adjoint -adept -analytic-jacobian -force-forward-jacobian $FILE > /dev/null
echoerr "*** PVC ADEPT JACOBIAN (REVERSE)"
$TIME $PVC_REDUCED -adjoint -adept -analytic-jacobian -force-reverse-jacobian $FILE > /dev/null
echoerr "*** PVC ADOLC JACOBIAN (DEFAULT)"
$TIME $PVC -adjoint -adolc -analytic-jacobian $FILE > /dev/null
echoerr "*** PVC ADOLC JACOBIAN (FORWARD)"
$TIME $PVC -adjoint -adolc -analytic-jacobian -force-forward-jacobian $FILE > /dev/null
echoerr "*** PVC ADOLC JACOBIAN (REVERSE)"
$TIME $PVC -adjoint -adolc -analytic-jacobian -force-reverse-jacobian $FILE > /dev/null
echoerr "*** PVC CPPAD JACOBIAN (DEFAULT)"
$TIME $PVC -adjoint -cppad -analytic-jacobian $FILE > /dev/null
echoerr "*** PVC CPPAD JACOBIAN (FORWARD)"
$TIME $PVC -adjoint -cppad -analytic-jacobian -force-forward-jacobian $FILE > /dev/null
echoerr "*** PVC CPPAD JACOBIAN (REVERSE)"
$TIME $PVC -adjoint -cppad -analytic-jacobian -force-reverse-jacobian $FILE > /dev/null
echoerr "*** PVC SACADO JACOBIAN"
$TIME $PVC -adjoint -sacado-fad -analytic-jacobian $FILE > /dev/null


echoerr "*** TDTS ADEPT JACOBIAN (DEFAULT)"
$TIME $TDTS_REDUCED -adjoint -adept -analytic-jacobian $FILE > /dev/null
echoerr "*** TDTS ADEPT JACOBIAN (FORWARD)"
$TIME $TDTS_REDUCED -adjoint -adept -analytic-jacobian -force-forward-jacobian $FILE > /dev/null
echoerr "*** TDTS ADEPT JACOBIAN (REVERSE)"
$TIME $TDTS_REDUCED -adjoint -adept -analytic-jacobian -force-reverse-jacobian $FILE > /dev/null
echoerr "*** TDTS ADOLC JACOBIAN (DEFAULT)"
$TIME $TDTS -adjoint -adolc -analytic-jacobian $FILE > /dev/null
echoerr "*** TDTS ADOLC JACOBIAN (FORWARD)"
$TIME $TDTS -adjoint -adolc -analytic-jacobian -force-forward-jacobian $FILE > /dev/null
echoerr "*** TDTS ADOLC JACOBIAN (REVERSE)"
$TIME $TDTS -adjoint -adolc -analytic-jacobian -force-reverse-jacobian $FILE > /dev/null
echoerr "*** TDTS CPPAD JACOBIAN (DEFAULT)"
$TIME $TDTS -adjoint -cppad -analytic-jacobian $FILE > /dev/null
echoerr "*** TDTS CPPAD JACOBIAN (FORWARD)"
$TIME $TDTS -adjoint -cppad -analytic-jacobian -force-forward-jacobian $FILE > /dev/null
echoerr "*** TDTS CPPAD JACOBIAN (REVERSE)"
$TIME $TDTS -adjoint -cppad -analytic-jacobian -force-reverse-jacobian $FILE > /dev/null
echoerr "*** TDTS SACADO JACOBIAN"
$TIME $TDTS -adjoint -sacado-fad -analytic-jacobian $FILE > /dev/null



