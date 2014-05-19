#!/bin/sh
NEWEST_LOG=`ls *.log | sort | tail -n 1`

echo Deleting "$NEWEST_LOG*"
rm -rf $NEWEST_LOG*


# rm -rf `ls *.png | sort | tail -n 3`
# rm -rf `ls *.log | sort | tail -n 1`
