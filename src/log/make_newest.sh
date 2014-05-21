#!/bin/sh
if [ -z "$1" ]; then
	ARG="1"
else
	ARG=$(($1*3))
fi

python ../pueschelplot.py `ls $LOG_PREFIX*.log | sort | tail -n $1 | head -n 1`
