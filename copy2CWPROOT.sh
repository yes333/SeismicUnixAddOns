#!/bin/bash

# copy extended header files to SU root src 

if [ -f $CWPROOT/include ];
then
mkdir $CWPROOT/include
fi

cp include/segyhdr.h $CWPROOT/include

cp include/tapehdr.h $CWPROOT/src/su/include

cp include/tapebhdr.h $CWPROOT/src/su/include

cp include/segy.h $CWPROOT/src/su/include

cp include/bhdr.h $CWPROOT/src/su/include

cp include/bheader.h $CWPROOT/src/su/include

cp include/cwpcmaps.h $CWPROOT/src/psplot/include

