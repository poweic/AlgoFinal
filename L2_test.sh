#!/bin/bash -e
cd htk && make && cd -
decoder=htk/HTKLVRec/HDecode.mod

# INPUTS -- mfcc files to decode
scp=data/dev.htk.scp

# INPUTS -- models
LM=data/lm.arpa.txt
AM=data/final.mmf 
lexicon=data/lexicon.txt
tiedlist=data/tiedlist

# OUTPUTS -- lattices and one-best recognition results
lat_dir=lat
rec_result=dev.rec

LOG=`mktemp`
echo "LOG save to $LOG"

# Options
beam=13.0
ac_scale=0.1
options="-A -T 1 -a $ac_scale -s 1.0 -t $beam -z lat -q tvaldm -o M"
$decoder $options -l $lat_dir -i $rec_result -w $LM -H $AM -S $scp $lexicon $tiedlist >> $LOG

sed -i "s%CPU time.*$%%g" $LOG 
#testing/diff.sh $LOG testing/L2_test.log
#rm $LOG
