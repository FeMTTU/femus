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

export LIBRARY_FOLDER=library

$JAVABIN/java -Xmx1024M -Djava.library.path="$INSTALLDIR/$LIBRARY_FOLDER" -Dhdfview.root="$INSTALLDIR" -jar "$INSTALLDIR/$LIBRARY_FOLDER/jhdfview.jar" $*

# alternate invocation when using modules
#$JAVABIN/java -Xmx1024M -Djava.library.path="$INSTALLDIR/$LIBRARY_FOLDER:$INSTALLDIR/$LIBRARY_FOLDER/ext" -Dhdfview.root="$INSTALLDIR" -cp "$INSTALLDIR/$LIBRARY_FOLDER/jhdfview.jar" ncsa.hdf.view.HDFView $*
