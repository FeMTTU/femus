#!/bin/bash

# Set up default variable values if not supplied by the user.
# This script file is used to execute the hdfview utility
# ... hdfview.root property is for the install location
# ...... default location is system property user.dir
# ... hdfview.workdir property is for the working location to find files
# ...... default location is system property user.home
#


# Adjust the following two variables to match your environment
export JAVABIN= # bin folder of installed java if not found below
 export INSTALLDIR=../ 

if [ ! -d "$JAVABIN" ]; then
    JAVABIN=$(dirname `which java`)
fi

$JAVABIN/java -Xmx1024M -Djava.library.path="$INSTALLDIR/lib" -Dhdfview.root="$INSTALLDIR" -jar "$INSTALLDIR/lib/jhdfview.jar" $*

# alternate invocation when using modules
#$JAVABIN/java -Xmx1024M -Djava.library.path="$INSTALLDIR/lib:$INSTALLDIR/lib/ext" -Dhdfview.root="$INSTALLDIR" -cp "$INSTALLDIR/lib/jhdfview.jar" ncsa.hdf.view.HDFView $*
