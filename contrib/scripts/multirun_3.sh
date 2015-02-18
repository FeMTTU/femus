#!/bin/bash


CONF_FILE="parameters.in"
MAIN_OUTPUT="./output/"

for mynlevs in 1 #$(seq 1 1 3) #$(seq 1.e6 1.e6 1.e7)
do
for mykappafrac in 1.e-7 #$(seq 1.e-8 1.e-8 1.e-7)
do
# for mypressout in $(seq 1.e6 1.e6 1.e7)
# do
for mykappawell in 1.e-10  #$(seq 1.e-11 1.e-11 1.e-10)
do
 for myyoungfrac in $(seq 1.e5 0.5e5 2.e5) #4e5 # #1.e10 #$(seq 1.e4 5.e4 1.e5)
  do
   for myyoungwell in 1.e10 #$(seq 1.e9 5.e9 1.e10)
    do
# prepare the folder for the single RUN
     echo -e "**********     nlevs = " $mynlevs     #"\n"
#      echo -e "**********     p_out = " $mypressout  #"\n"
     echo -e "********** kappafrac = " $mykappafrac #$mykappafrac #"\n"
     echo -e "********** kappawell = " $mykappawell #"\n"
     echo -e "********** youngfrac = " $myyoungfrac #"\n"
     echo -e "********** youngwell = " $myyoungwell #"\n"
     OUTDIR_TIME=$(date +%Y_%m_%d_%H_%M_%S)
     export OUTFOLDER=$MAIN_OUTPUT$OUTDIR_TIME
# invece di fare export e getenv, potrei scrivere su un file "incoming_run" e leggere da quel file nel main,
#e' piu' portable
    # check that input is there, or make it
    # check that output is there, or make it
    # I think i can do this with
     mkdir -p input
     mkdir -p output
     mkdir $OUTFOLDER
# Copy the template parameters file in that folder, AND THEN MODIFY IT, and then read from that
# E' giusto fare cosi'! cioe' copiare e leggere DOPO, per evitare il rischio di possibili SOVRAPPOSIZIONI di file!
     cp $CONF_FILE $OUTFOLDER/  #it is better to copy and read from the copied   
     sed '/nlevs/ c\nlevs '$mynlevs' '               -i  $OUTFOLDER/${CONF_FILE}
#      sed '/p_out/ c\p_out '$mypressout' '            -i  $OUTFOLDER/${CONF_FILE}
     sed '/kappa_frac/ c\kappa_frac '$mykappafrac' ' -i  $OUTFOLDER/${CONF_FILE}   #'$mykappafrac'
     sed '/kappa_well/ c\kappa_well '$mykappawell' ' -i  $OUTFOLDER/${CONF_FILE}
     sed '/young_frac/ c\young_frac '$myyoungfrac' ' -i  $OUTFOLDER/${CONF_FILE}
     sed '/young_well/ c\young_well '$myyoungwell' ' -i  $OUTFOLDER/${CONF_FILE}
# Run the executable
# Here I could decide to show BOTH TO SCREEN AND TO FILE
# HOW DO I DO DEBUG WITH THIS SCRIPT?!?
    ./main.out 2> ./${OUTFOLDER}/run_log.txt 1> ./${OUTFOLDER}/run_log.txt  #HOW COME 2>&1 DOES NOT WORK?!?
    
    sleep 1
    done
  done
done

# done # pressure
done

done

#if you want to avoid intrusion in the code for printing the output, you can just take the output 
# at the end of the simulation and COPY IT to the related directory (of course you cannot superimpose but it works)
# Commenti
# il comparatore di uguaglianza per la shell ha UN SOLO segno di uguale!

# So, i think it is clear that the generation of the time string OUTSIDE of the main() function, as a shell script,
# is BETTER, because it collects ALL THE OUTPUTS of ALL THE LINKED LIBRARIES. (sometimes i get an "Aborted" outside, but anyway...)
# therefore the "parameter sweep" is performed OUTSIDE of the main function also.
# It cannot be done inside, which would be nice if you had some sort of "ADAPTIVE Parameter Sweep" that decides 
# what parameters to choose BASED on SOME KIND of ESTIMATOR (like error estimators in adaptive finite element).
# then you would have some routine inside the code that generates the new parameters to be considered, 
# and so you would need to perform the parameter loops INSIDE the MAIN(), and therefore the generation 
# of the output folder would be inside as well.
# But, for us, since the choice of the parameters to sweep is done "a priori", then we do everything OUTSIDE.
# In fact the scheme is:
# DECIDE THE VECTOR OF PARAMETERS (by some criterion)
# GENERATE THE OUTPUT DIRECTORY
# COPY THE PARAMETER FILE IN THAT, THEN MODIFY IT
# LAUNCH THE EXECUTABLE AND READ THE PARAMETERS FROM THAT FILE

# if we want to keep the output directory OUTSIDE but the parameter choice INSIDE, we should instead 
# PREPARE an OUTPUT DIRECTORY for some "INCOMING VECTOR OF PARAMETERS"
# DECIDE THE VECTOR OF PARAMETERS (by some criterion)
# ASSOCIATE this NEWLY CREATED vector of parameters to the NEWLY CREATED OUTPUT DIRECTORY (criterion to be sure it is new: the output directory is EMPTY. We must be sure that all the other ones are filled with something)
# COPY THE PARAMETER FILE, THEN MODIFY IT
# LAUNCH THE EXECUTABLE AND READ THE PARAMETERS FROM THAT FILE

# This second one is more complicated. One thing it would allow is the MANAGEMENT OF THE PARAMETERS LOOP 
# to be ONLY INSIDE THE MAIN(), which is good in the sense that it avoids "uncoupling" between the shell script and the main() loop.
# If it were all inside, the parameters could be treated sort of as a "class", "SweptParameter".

# again, the key point is that should perform the REDIRECT operation INSIDE THE MAIN(). But this time not by redirecting 
# the c++ std::cout, but i would say with a SYSTEM CALL inside the main... 
# but, once an executable is launched, is it possible to redirect all its outputs during the execution?!?
# I dont know any shell command that does that... you should stop the program, and then relaunch it...
# I dont think that's possible. So, I guess the SHELL SCRIPT is the best solution. Still, the point is inconvenient is that
# the SHELL SCRIPT and the MAIN() must be AWARE OF EACH OTHER, in terms of the "SWEPT PARAMETERS"
# Maybe doing things INSIDE the MAIN would also allow AVOIDING some REALLOCATIONS, or new operations, or something.
# So, there are some pros and cons. Still, the fact of collecting all the outputs is essential to me, 
# i want to put it as a priority.
# This implies that the loop and the main must be handled accordingly, because they are not aware of each other.

# Another disadvantage of the external shell loop is the use of the DEBUGGER in kdevelop. Still acceptable, even though i would like to find a solution for that.
# In fact, if you have to change an environment variable, then every time you launch INSIDE kdevelop you would need 
# to CLOSE KDEVELOP and REOPEN IT WITH THE NEW SHELL ENVIRONMENT, which is not good at all.
# That's one disadvantage in the use of environment variables for runtime setup of executions.

# How about the management of the COMMON OUTPUT BigTable file? is it better to manage it OUTSIDE or INSIDE?
# Well, the numbers to be printed are generated INSIDE of course, so it is better to generate INSIDE. 
# since the input numbers are generated OUTSIDE instead, i could set them up outside actually... but that would imply
# a sort of splitting in which the bigtable file is handled partly outside (for the inputs) and partly inside (for the outputs).
# Plus, how do i say in a shell script to APPEND to file? Yes, I use the double arrow ">>" instead of the single

