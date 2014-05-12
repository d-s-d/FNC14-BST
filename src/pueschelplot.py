#!/usr/bin/env python
__doc__="""
""";

import json, sys
import matplotlib.pyplot as plt
import numpy as np

def usage():
    print __doc__

def main():
    try:
        logfile = open(sys.argv[1])
    except KeyError:
        usage()
        sys.exit(1)
    
    try:
        jsonLog = json.load(logfile)
    except Exception as e:
        print( "exception type: %s, message: %s" % (type(e), e))
        sys.exit(1)
    
    cycles = {}
    runs = jsonLog['runs']
    if len(runs) > 0:
        N = [int(x['N']) for x in runs[runs.keys()[0]]]
        for implName, results in runs.iteritems():
             if [int(x['N']) for x in results] == N:
                 cycles[implName] = [int(x['cycles']) for x in results]

    yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

    plt.subplot(111,axisbg='#BBBBBB',alpha=0.1)
    plt.grid(color='white', alpha=0.5, linewidth=2, linestyle='-', axis='y')

    for spine_name in ['top', 'left', 'right']:
        plt.gca().spines[spine_name].set_color('none')
    
    plt.ylabel('Runtime [cycles]', **yprops)
    plt.xlabel('N [doubles]')

    plt.gca().tick_params(direction='out', length=0, color='k')

    for implName, results in cycles.iteritems():
        plt.plot(N, results, 'o-', linewidth=2, label=implName)
    plt.legend(loc=2)
    

    plt.gca().set_axisbelow(True)
    
    plt.savefig('pueschelplot.png', dpi=300)
    
    plt.show()

if __name__=="__main__":
    main()
