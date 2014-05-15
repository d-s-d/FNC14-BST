#!/usr/bin/env python
__doc__="""
""";

import sys
import json
import numpy as np

#import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg') # for pyplot to work without X, e.g. SSH session
from matplotlib import pyplot as plt

def usage():
    print("Usage: %s logfile") % sys.argv[0]
    #print __doc__

def plot(title, xlabel, ylabel, data, log, filename):
    yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

    plt.clf()
    plt.subplot(111,axisbg='#BBBBBB',alpha=0.1)
    plt.grid(color='white', alpha=0.5, linewidth=2, linestyle='-', axis='y')

    for spine_name in ['top', 'left', 'right']:
        plt.gca().spines[spine_name].set_color('none')

    plt.gca().tick_params(direction='out', length=0, color='k')
    plt.gca().set_axisbelow(True)
    
    #plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel, **yprops)

    N = range(log['from'], log['to']+1, log['step'])

    for implName, results in data.iteritems():
        plt.plot(N, results, 'o-', linewidth=2, label=implName)
    plt.legend(loc=2)

    # some space on the bottom for comments
    plt.subplots_adjust(bottom=0.2)
    plt.figtext(0.1, 0.02,
            "N = %s:%s:%s; " % (log['from'], log['step'], log['to']) +
            "Seed = %s; " % log['seed'] +
            "Callibration = %s;\n" % log['callibration'] +
            "CFLAGS: %s\n" % log['userflags'] +
            "Git Revision: %s" % log['git-revision'])
    
    plt.savefig(filename, dpi=75)
    print("Graph generated: %s" % filename)

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
    fl_add = [(1/3.0*(n**3 + 3*n**2 + 2*n)) for n in N] 
    fl_cmp = [(1/6.0*(n**3 - n)) for n in N] 
    flops  = fl_add + fl_cmp

    cycles = {}
    c_miss = {}
    c_read = {}

    performance = {}
    hit_rate    = {}

    for implName, results in runs.iteritems():
         cycles[implName] = [int(x['cycles']) for x in results]
         c_miss[implName] = [int(x['cache-misses']) for x in results]
         c_read[implName] = [int(x['cache-references']) for x in results]

         performance[implName] = [fl/c for (c,fl) in zip(cycles[implName],flops)]
         hit_rate[implName] = [1-float(m)/float(r) if r > 0 else 1 for (m,r) in zip(c_miss[implName], c_read[implName])]

    # data plotting
    plot('Runtime', 'N [doubles]', 'Runtime [cycles]',
            cycles, jsonLog, filename + '.cycles.png')
    plot('Performance', 'N [doubles]', 'Performance [flops/cycles]',
            performance, jsonLog, filename + '.performance.png')
    plot('Cache Hit Rate', 'N [doubles]', 'Cache Hit Rate',
            hit_rate, jsonLog, filename + '.cache.png')

if __name__=="__main__":
    main()
