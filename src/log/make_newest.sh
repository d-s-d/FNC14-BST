#!/bin/sh
if [ -z "$1" ]; then
	$1="1"
fi

python ../pueschelplot.py `ls $LOG_PREFIX*.log | sort | tail -n $1 | head -n 1`
