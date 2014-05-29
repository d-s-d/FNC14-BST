#!/usr/bin/python

import sys

if __name__ == "__main__":
	name = sys.argv[1]

	f_trans = open("bytes_transferred_bst"+name+".txt")
	transferred = [int(n) for n in f_trans.read().split(' ')]
	f_flops = open("flop_bst"+name+".txt")
	flops = [int(n) for n in f_flops.read().split(' ')]
	f_cycles = open("tsc_bst"+name+".txt")
	cycles = [int(n) for n in f_cycles.read().split(' ')]

	performance = [float(fl)/float(c) for (fl,c) in zip(flops, cycles)]
	opint = [float(fl)/float(b) for (fl,b) in zip(flops, transferred)]

	f_res = open("roofline_data_"+name+".txt", "w")
	print >> f_res, "OpInt [X] %s;Perf [Y] %s" % (name, name)
	for i in range(len(opint)):
		print >> f_res, "%lf;%lf" % (opint[i], performance[i])
