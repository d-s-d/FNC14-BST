#!/usr/bin/env python
__doc__="""
""";

import sys
import json
import numpy as np

def usage():
    print("Usage: %s logfile") % sys.argv[0]
    #print __doc__

def main():
    # input parsing
    filename = sys.argv[1]
    try:
        logfile = open(filename)
    except (KeyError, IndexError):
        usage()
        sys.exit(1)
    
    try:
        jsonLog = json.load(logfile)
    except Exception as e:
        print( "exception type: %s, message: %s" % (type(e), e))
        sys.exit(1)
    
    runs = jsonLog['runs']
    if len(runs) == 0:
        print("No data in logfile.")
        sys.exit(1)

    # data extraction
    N = range(jsonLog['from'], jsonLog['to']+1, jsonLog['step'])
    # fl_add = [(1/3.0*(n**3 + 3*n**2 + 2*n)) for n in N] 
    # fl_cmp = [(1/6.0*(n**3 - n)) for n in N] 
    # flops  = fl_add + fl_cmp # concat instead of vector add
    # flops = [sum(x) for x in zip(fl_add,fl_cmp)]
    # flops = [(n**3 + 4*n**2 + 3*n)/2.0 + n for n in N] # 110_precompute
    # flops = [(n**3/2.0 + 5*n**2/2.0 + 2*n) + n for n in N] # 109_init_opt

    flops  = {}
    cycles = {}
    c_miss = {}

    performance = {}
    hit_rate    = {}

    for implName, results in runs.iteritems():
        flops[implName]  = [float(x['flops']) for x in results]
        cycles[implName] = [int(x['cycles']) for x in results]
        c_miss[implName] = [int(x['cache-misses']) for x in results]

        performance[implName] = [fl/c for (c,fl) in zip(cycles[implName],flops[implName])]

    # output as csv
    with open(filename + ".csv", "w") as csv:
        keys = performance.keys()
        print >> csv, "N,",
        print >> csv, ",".join(keys)
        for i in range(len(N)):
            print >> csv, N[i], ",",
            print >> csv, ",".join([str(performance[k][i]) for k in keys])

if __name__=="__main__":
    main()
