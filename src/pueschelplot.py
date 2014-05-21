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
from matplotlib.pyplot import *
from matplotlib import gridspec

def usage():
    print("Usage: %s logfile") % sys.argv[0]
    #print __doc__

def plot(title, xlabel, ylabel, data, log, filename):
    yprops = dict(rotation=0, y=1.05, horizontalalignment='left')


    plt.clf()
    fig = plt.figure(figsize=(8, 8))
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(13, 6)
    gs = gridspec.GridSpec(1, 2, width_ratios=[5, 2])

    plt.subplot(gs[0],axisbg='#BBBBBB',alpha=0.1)
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

    #plt.legend(loc=3, prop={'size':10})
    legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    # some space on the bottom for comments
    plt.subplots_adjust(bottom=0.2)
    plt.figtext(0.1, 0.02,
            "N = %s:%s:%s; " % (log['from'], log['step'], log['to']) +
            "Seed = %s; " % log['seed'] +
            "Calibration = %s;\n" % log['calibration'] +
            "CFLAGS: %s\n" % log['userflags'] +
            "Git Revision: %s" % log['git-revision'])

    plt.savefig(filename, dpi=150)

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

    # data plotting
    plot('Runtime', 'N [doubles]', 'Runtime [cycles]',
            cycles, jsonLog, filename + '.cycles.png')
    plot('Performance', 'N [doubles]', 'Performance [flops/cycles]',
            performance, jsonLog, filename + '.performance.png')
    plot('', 'N [doubles]', 'Last-Level Cache Misses',
            c_miss, jsonLog, filename + '.cache.png')

if __name__=="__main__":
    main()
