#!/bin/bash

S0=`basename $0`
D0=`dirname $0`

cd $D0

OUTDIR=covdata

rm -rf $OUTDIR
mkdir -p $OUTDIR

coverage3 run ../run_unfoc.py -h > /dev/null
#coverage3 run -a ../unfoc.py -h > /dev/null
#coverage3 run -a ../unfocb.py -h > /dev/null
coverage3 run -a ../unfoc/read.py -h > /dev/null
#coverage3 run -a ../read.py -i /disk/kea/WAIS/orig/xlob/NIS4/IBH0e/X84b/RADnh5/bxds --mmap
#coverage3 run -a ../read.py -i /disk/kea/WAIS/orig/xlob/NIS4/IBH0e/X84b/RADnh5/bxds
# Try with a RADnh3.
#coverage3 run -a ../read.py -i /disk/kea/WAIS/orig/xlob/ICP4/JKB2g/F19T03a/RADnh3/bxds --mmap
#coverage3 run -a ../read.py -i /disk/kea/WAIS/orig/xlob/ICP4/JKB2g/F19T03a/RADnh3/bxds

coverage3 run -a ./test_parse_channels.py
coverage3 run -a ./test_read.py
coverage3 run -a ./test_unfoc1.py

coverage3 run -a ../run_unfoc.py -i /disk/kea/WAIS/orig/xlob/NIS4/IBH0e/X84b/RADnh5/bxds \
                              -o $OUTDIR/testcmd --nmax 10 \
                              --stackdepth 10 --incodepth 5 --channels LoResInco1


coverage3 run -a ../run_unfoc_1m.py -i /disk/kea/WAIS/targ/xtra/KRT2/FOC/Best_Versions/S2_FIL/NIS4/IBH0e/X84b/bxds2.i \
                                 --stackdepth 1 --incodepth 1  -o $OUTDIR/testcmd --nmax 100 --c 2

rm -rf $OUTDIR

echo "$S0: coverage3 report -m to show coverage results"
