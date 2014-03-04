#!/bin/bash


CONF_FILE="./config/femus_conf.in"

ALPHAMAX=100000

make -j

for myval in $(seq -2 0.5 +2)    #100000   #$(seq 1000000000 1000000000 ${ALPHAMAX})
do
     echo -e "********** INJSUC = " $myval "\n"
#      sed '/alphaT/ c\ alphaT '${myval}' ' -i  ${CONF_FILE}
     sed '/injsuc/ c\ injsuc '${myval}' ' -i  ${CONF_FILE}
    ./temper-$FM_FEMUS_METHOD
    sleep 1
done


# Commenti
# il comparatore di uguaglianza per la shell ha UN SOLO segno di uguale!
