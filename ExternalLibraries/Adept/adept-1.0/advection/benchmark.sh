export LD_LIBRARY_PATH=$HOME/lib 
PROGRAM=./run_advection_benchmark
NT=2000
NR=1000
CN=0.125

SCHEME=lax_wendroff
$PROGRAM $SCHEME $NT $NR $CN

SCHEME=toon
$PROGRAM $SCHEME $NT $NR $CN

SCHEME=lax_wendroff
$PROGRAM $SCHEME $NT $NR $CN forward-jacobian
$PROGRAM $SCHEME $NT $NR $CN reverse-jacobian

SCHEME=toon
$PROGRAM $SCHEME $NT $NR $CN forward-jacobian
$PROGRAM $SCHEME $NT $NR $CN reverse-jacobian


