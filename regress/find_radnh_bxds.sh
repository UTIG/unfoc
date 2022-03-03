#!/bin/bash -e

# Collect a list of all RADnh* bxds files for use in directed random testing.
D0=`dirname $0`

TESTLIST=$D0/test_lists/available_radnh_bxds.txt

mkdir -p `dirname $TESTLIST`

# Find bxds files,
# but not RADnh4
# and not ground transects (e.g., F01G01a)
find $WAIS/orig/xlob/ -iname bxds \
   | fgrep '/RADnh' \
   | grep -v '/RADnh4/' \
   | grep -v -P '/F\d\dG\d\d\w/' \
   | sort > $TESTLIST

echo "Tests found: " `wc -l $TESTLIST`

# Compute total disk size
cat $TESTLIST | tr '\n' '\0' | du -ch --files0-from - | tail -n 1
