#!/bin/sh
if [ -z "$1" ]; then
	$1="1"
fi

python ../pueschelplot.py `ls *.log | sort | tail -n $1 | head -n 1`
