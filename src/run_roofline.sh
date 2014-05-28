#!/bin/sh

sudo ./roofline_driver $@
mv *.txt roofline/
cd roofline && scala AdvancedRoofline.scala
