#!/bin/bash -e
(cd htk && make && cd -) || (printf "\33[31m[Error]\33[0m failed to make htk\n"; exit -1)
printf "\n\33[34mStart testing...\33[0m\n"
decoder=htk/HTKLVRec/HDecode.mod

# INPUTS -- mfcc files to decode
scp=data/small.scp 

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
for beam in 7.0 9.0 11.0; do
  for ac_scale in 0.1 0.3 0.5; do
    options="-A -T 1 -a $ac_scale -s 1.0 -t $beam -z lat -q tvaldm -o M"
    $decoder $options -l $lat_dir -i $rec_result -w $LM -H $AM -S $scp $lexicon $tiedlist >> $LOG
  done
done

sed -i "s%CPU time.*$%%g" $LOG 
testing/diff.sh $LOG testing/L1_test.log
rm $LOG
