#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

xvals = (112, 224, 336, 448, 560, 672, 784, 896, 1008, 1120, 1232, 1344, 1456)
yvals1 = (1.918499, 1.580462, 1.813013, 1.808579, 1.729736, 1.729384, 1.715734, 1.386305, 1.729841, 1.696438, 1.616113, 1.705391, 1.720684)
yvals2 =  (0.618573, 0.545513, 0.574069, 0.459962, 0.567306, 0.519556, 0.557664, 0.417062, 0.566603, 0.523055, 0.558473, 0.479710, 0.561821)

yprops = dict(rotation=0, y=1.05, horizontalalignment='left')

plt.subplot(111,axisbg='#BBBBBB',alpha=0.1)
plt.grid(color='white', alpha=0.5, linewidth=2, linestyle='-', axis='y')

for spine_name in ['top', 'left', 'right']:
    plt.gca().spines[spine_name].set_color('none')
    
plt.ylabel('Performance [flops/cycle]', **yprops)
plt.xlabel('N [doubles]')

plt.gca().tick_params(direction='out', length=0, color='k')

plt.plot(xvals, yvals1, 'bo-', linewidth=2)
plt.plot(xvals, yvals2, 'go-', linewidth=2)
plt.gca().set_axisbelow(True)

plt.savefig('test.png', dpi=300)

plt.show()

