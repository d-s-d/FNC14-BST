#!/bin/sh
if [ -z "$1" ]; then
	$1="1"
fi

start=$(($1 * 3))
okular `ls *.png | sort | tail -n $start | head -n 3`
