#!/bin/sh
python ../pueschelplot.py `ls *.log | sort | tail -n 1`
