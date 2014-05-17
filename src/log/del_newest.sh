#!/bin/sh
rm -rf `ls *.png | sort | tail -n 3`
rm -rf `ls *.log | sort | tail -n 1`
