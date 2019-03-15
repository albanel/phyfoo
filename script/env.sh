#!/bin/bash
#source  /local/env/envpython-2.7.sh
source /local/env/envpython-2.7.15.sh

SRC=/groups/dyliss/orthocis_albane/pipeline/src
CONFFILE=config/latest.conf
if [ -f $CONFFILE ]; then
 #TODO: loop to export all keys 
 while read p; do

    p2="$(echo -e "${p}" | tr -d '[:space:]')"
    if [ ! -z $p2 ]; then
     if [  $p2 != "" ]; then
        echo "export  $p "
        ` export $p `
        export $p

     fi
    fi

 done < $CONFFILE
fi


MODINSTALLDIR=${SRC}/pymodule
MODORTHODIR=${SRC}/orthocismodule
HOMOGEN=${SRC}/script/homology
PYTHONPATH=$PYTHONPATH:$MODINSTALLDIR:$MODORTHODIR:$HOMOGEN
export PYTHONPATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/softs/local/gc/gc-7.2/.libs
export ORTHO_CPP_BIN_DIR=/groups/dyliss/orthocis_albane/pipeline/src/script/homology/
