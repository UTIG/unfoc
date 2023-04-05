#!/bin/bash -e

# Run regression and coverage

COV=coverage3


S0=`basename $0`
D0=`dirname $0`

echo "Changing to $D0"
pushd $D0

# Quick data testing
#OUTDIR=./covdata_quick
#UNFOCFLAGS="--nmax 1000"
OUTDIR=./covdata
UNFOCFLAGS="--nmax 10000"

RCFILE=`pwd`/.coveragerc

export PYTHONPATH=${PYTHONPATH}:`pwd`

export COVERAGE_PROCESS_START=$RCFILE

#$COV erase
rm -rf $OUTDIR

COVFLAGS="-a"

$COV run ../unfoc/parse_channels.py
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
#MARFA_CHANNELS="[1,1,1,3,1;5,1,1,0,0;7,3,1,0,0;2,2,1,4,1;6,2,1,0,0;8,4,1,0,0]"
MARFA_CHANNELS="LoResInco1,LoResInco5,LoResInco3,LoResInco2,LoResInco6,LoResInco8"


#echo "${S0}: Checksum input data"
# TODO: put this into a makefile
#stdbuf -o 0 sha1sum -b \
#    $WAIS/orig/xlob/MBL/JKB2h/Y90a/RADnh3/bxds \
#    $WAIS/orig/xlob/MBL/MKB2l/Y76a/RADnh3/bxds \
#    $WAIS/orig/xlob/ASB1/GCX0e/R40a/RADnh5/bxds \
#    $WAIS/orig/xlob/ASB1/JKB2s/R06a/RADnh5/bxds \
#    | tee $OUTDIR/input_bxds.sha1




echo "${S0}: Running tests"
###############
#test_HiCARS2

# Run a short section for coverage and no blanking
$COV run $COVFLAGS ../run_unfoc.py \
    -i $WAIS/orig/xlob/MBL/JKB2h/Y90a/RADnh3/bxds \
    --outdir $OUTDIR/MBL/JKB2h/Y90a/RADnh3 \
    --channels LoResInco1,LoResInco2 --output_samples 3200 \
    --stackdepth 10 --incodepth 5  --scale 20000 --nmax 100

# Run a short section for coverage and no blanking
$COV run $COVFLAGS ../run_unfoc.py \
    -i $WAIS/orig/xlob/MBL/JKB2h/Y90a/RADnh3/bxds \
    --outdir $OUTDIR/MBL/JKB2h/Y90a/RADnh3 \
    --channels LoResInco1,LoResInco2 --output_samples 3200 \
    --stackdepth 10 --incodepth 5  --scale 20000 --blanking -200 --nmax 100

# Run the same section with blanking from the end
$COV run $COVFLAGS ../run_unfoc.py \
    -i $WAIS/orig/xlob/MBL/JKB2h/Y90a/RADnh3/bxds \
    --outdir $OUTDIR/MBL/JKB2h/Y90a/RADnh3 \
    --channels LoResInco1,LoResInco2 --output_samples 3200 --output_phases \
    --stackdepth 10 --incodepth 5  --scale 20000  --blanking 200 $UNFOCFLAGS

# test_MARFA
$COV run $COVFLAGS ../run_unfoc.py \
    -i $WAIS/orig/xlob/MBL/MKB2l/Y76a/RADnh3/bxds \
    --outdir $OUTDIR/MBL/MKB2l/Y76a/RADnh3 \
    --channels LoResInco2  --output_samples 3200 --output_phases --bandpass \
    --stackdepth 10 --incodepth 5  --scale 20000  --blanking 200 $UNFOCFLAGS


# Read a RADnh5 from GCX (HiCARS2)
$COV run $COVFLAGS ../run_unfoc.py \
    -i $WAIS/orig/xlob/ASB1/GCX0e/R40a/RADnh5/bxds \
    --outdir $OUTDIR/ASB1/GCX0e/R40a/RADnh5 \
    --channels LoResInco1,LoResInco2  --output_samples 3200 --output_phases --bandpass \
    --stackdepth 10 --incodepth 5  --scale 20000  --blanking 200 $UNFOCFLAGS --nmax 2000


# Read a RADnh5 from JKB (MARFA)
$COV run $COVFLAGS ../run_unfoc.py \
    -i $WAIS/orig/xlob/ASB1/JKB2s/R06a/RADnh5/bxds \
    --outdir $OUTDIR/ASB1/JKB2s/R06a/RADnh5 \
    --channels  ${MARFA_CHANNELS}   --output_samples 3200 --output_phases --bandpass \
    --stackdepth 10 --incodepth 5  --scale 20000  --blanking 200


# TODO: use a shortened bxds file

# Use --channels flag
$COV run $COVFLAGS ../run_unfoc.py \
    -i $WAIS/orig/xlob/THW/PBA0a/X66a/RADnh5/bxds \
    --outdir $OUTDIR/unfoc1_THW_PBA0a_X66a_a \
    --channels LoResInco1,LoResInco2,LoResInco3,LoResInco4,LoResInco5,LoResInco6,LoResInco7,LoResInco8 \
    --output_samples 3200 --output_phases --bandpass --stackdepth 10 \
    --incodepth 10 --scale 20000 --blanking 200

# Use --channels flag with just channels 7,8
$COV run $COVFLAGS ../run_unfoc.py \
    -i $WAIS/orig/xlob/THW/PBA0a/X66a/RADnh5/bxds \
    --outdir $OUTDIR/unfoc1_THW_PBA0a_X66a_b \
    --channels LoResInco7,LoResInco8 \
    --output_samples 3200 --output_phases --bandpass --stackdepth 10 \
    --incodepth 10 --scale 20000 --blanking 200

diff $OUTDIR/unfoc1_THW_PBA0a_X66a_{a,b}/MagLoResInco7 && true
diff $OUTDIR/unfoc1_THW_PBA0a_X66a_{a,b}/MagLoResInco8 && true

# Checksum all output data
find $OUTDIR -type f -print0 | sort -z | xargs -r0 sha1sum -b > ${OUTDIR}.sha1
echo "${S0}: Wrote sha1sum to ${OUTDIR}.sha1"





#foreach x ( test_MARFA )
#    echo "${S0}: ./$x"
#    ./$x
#end

echo "${S0}: Run '$COV report -m' to view coverage"
