#!/bin/bash

#$ -o /groups/dyliss/orthocis_albane/pipeline/log//run.log
#$ -e /groups/dyliss/orthocis_albane/pipeline/log//run.err

date

echo " start run "

LDIR=/groups/dyliss/orthocis_albane/pipeline/log/
[ -d $LDIR ] || mkdir $LDIR


SRCSC="/local/env/envpython-2.7.15.sh"
if [ "$SRCSC" != "" ]; then 
  source $SRCSC
fi

SRC=/groups/dyliss/orthocis_albane/pipeline/src
cd $SRC
MODINSTALLDIR=${SRC}/pymodule
MODORTHODIR=${SRC}/orthocismodule
CHECKDIR=${SRC}/check_database
PYTHONPATH=$PYTHONPATH:$MODINSTALLDIR:$MODORTHODIR:$CHECKDIR
export PYTHONPATH
CONF=config/latest.conf


#STEPFILE=config/step1.txt
#STEPFILE=config/step2.txt
#STEPFILE=config/step3.txt
#STEPFILE=config/step4.txt
#STEPFILE=config/step5.txt
#STEPFILE=config/step6.txt
#STEPFILE=config/step7.txt

STEPFILE=config/all_steps.txt

fname=$(basename $STEPFILE)
fbname=${fname%.*}
LOG=$LDIR."".$fbname".log"
LOGERR=$LDIR."".$fbname".err"
echo $PYTHONPATH

python BuildOrthocis_check.py -c $CONF -s $STEPFILE -l $LDIR
cat $CONF >$LDIR."/confile.txt"
echo `hostname`  >> $LDIR."/confile.txt"

echo " end run "

date
