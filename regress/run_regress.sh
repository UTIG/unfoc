#!/bin/bash -e

# Run regression and coverage

COV=coverage


S0=`basename $0`
D0=`dirname $0`

echo "Changing to $D0"
pushd $D0
OUTDIR=./covdata

RCFILE=`pwd`/.coveragerc

export PYTHONPATH=${PYTHONPATH}:`pwd`

export COVERAGE_PROCESS_START=$RCFILE

#$COV erase
rm -rf $OUTDIR

COVFLAGS="-a"
# max number of samples to process
UNFOCFLAGS="--nmax 100"

$COV run ../parse_channels.py
echo "Done parse_channel"
# Outputs channels 1 and 2 unmodified.
HICARS2_CHANNELS="[1,1,1,0,0;2,2,1,0,0]"
# 1-4 are pass through
# 5&6 are low gain sum/diff
# 7&8 are high gain sum/diff
# TODO: work out which side is which for the helo
# TODO: Maybe this comment should be moved to pik1/parse_channels.py?
# chan1 -> low gain,
# chan2 -> high gain,
# chan3 -> low gain,
# chan4 -> high gain,
# from ICP6 when the two-phase-center stuff was new
#CHANNELS="[1,1,1,0,0;2,2,1,0,0;3,3,1,0,0;4,4,1,0,0;5,1,1,3,1;6,1,1,3,-1;7,2,1,4,1;8,2,1,4,-1]"
    # Matches how we do it on melt - there's no need to process the differences,
    # and we want to keep ch1 and ch2 with the same meaning as in HiCARS2.
MARFA_CHANNELS="[1,1,1,3,1;5,1,1,0,0;7,3,1,0,0;2,2,1,4,1;6,2,1,0,0;8,4,1,0,0]"


echo "${S0}: Start unfoc"
###############
#test_HiCARS2
echo $COV run $COVFLAGS ../unfoc.py \
    --infile $WAIS/orig/xlob/MBL/JKB2h/Y90a/RADnh3/bxds \
    --outdir $OUTDIR/MBL/JKB2h/Y90a/RADnh3 \
    --channel_def ${HICARS2_CHANNELS} --output_samples 3200 --output_phases \
    --StackDepth 10 --IncoDepth 5  --Scale 20000  --blanking 200 $UNFOCFLAGS

$COV run $COVFLAGS ../unfoc.py \
    --infile $WAIS/orig/xlob/MBL/JKB2h/Y90a/RADnh3/bxds \
    --outdir $OUTDIR/MBL/JKB2h/Y90a/RADnh3 \
    --channel_def ${HICARS2_CHANNELS} --output_samples 3200 --output_phases \
    --StackDepth 10 --IncoDepth 5  --Scale 20000  --blanking 200 $UNFOCFLAGS

# test_MARFA
$COV run $COVFLAGS ../unfoc.py \
    --infile $WAIS/orig/xlob/MBL/MKB2l/Y76a/RADnh3/bxds \
    --outdir $OUTDIR/MBL/MKB2l/Y76a/RADnh3 \
    --channel_def 2  --output_samples 3200 --output_phases \
    --StackDepth 10 --IncoDepth 5  --Scale 20000  --blanking 200 $UNFOCFLAGS


# Read a RADnh5 from GCX (HiCARS2)
$COV run $COVFLAGS ../unfoc.py \
    --infile $WAIS/orig/xlob/ASB1/GCX0e/R40a/RADnh5/bxds \
    --outdir $OUTDIR/ASB1/GCX0e/R40a/RADnh5 \
    --channel_def ${HICARS2_CHANNELS}  --output_samples 3200 --output_phases \
    --StackDepth 10 --IncoDepth 5  --Scale 20000  --blanking 200 $UNFOCFLAGS


# Read a RADnh5 from JKB (MARFA)
$COV run $COVFLAGS ../unfoc.py \
    --infile $WAIS/orig/xlob/ASB1/JKB2s/R06a/RADnh5/bxds \
    --outdir $OUTDIR/ASB1/JKB2s/R06a/RADnh5 \
    --channel_def  ${MARFA_CHANNELS}   --output_samples 3200 --output_phases \
    --StackDepth 10 --IncoDepth 5  --Scale 20000  --blanking 200 $UNFOCFLAGS




#foreach x ( test_MARFA )
#    echo "${S0}: ./$x"
#    ./$x
#end

echo "${S0}: Run '$COV combine && $COV report' to view coverage"
