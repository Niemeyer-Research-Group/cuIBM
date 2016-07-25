#!/usr/bin/env python

#import csv
#import argparse
import numpy as np
from numpy import genfromtxt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
#import os
#import os.path
#import sys

cuibmFolder = '/scratch/src/cuIBM'
caseFolder     = cuibmFolder + '/validation/osc/'
caseFolder =''

print "-"*80
print "Making validaiton plot for an oscillating cylinder in flow."
print "-"*80

a = genfromtxt('%sab/forcesa' % caseFolder,dtype=float,delimiter='\t')
b = genfromtxt('%sab/forcesb' % caseFolder,dtype=float,delimiter='\t')
c = genfromtxt('%scd/forcesc' % caseFolder,dtype=float,delimiter='\t')
d = genfromtxt('%scd/forcesd' % caseFolder,dtype=float,delimiter='\t')
e = genfromtxt('%sef/forcese' % caseFolder,dtype=float,delimiter='\t')
f = genfromtxt('%sef/forcesf' % caseFolder,dtype=float,delimiter='\t')
g = genfromtxt('%sgh/forcesg' % caseFolder,dtype=float,delimiter='\t')
h = genfromtxt('%sgh/forcesh' % caseFolder,dtype=float,delimiter='\t')
atime = [a[i][0] for i in xrange(1,len(a))]
btime = [b[i][0] for i in xrange(1,len(b))]
ctime = [c[i][0] for i in xrange(1,len(c))]
dtime = [d[i][0] for i in xrange(1,len(d))]
etime = [e[i][0] for i in xrange(1,len(e))]
ftime = [f[i][0] for i in xrange(1,len(f))]
gtime = [g[i][0] for i in xrange(1,len(g))]
htime = [h[i][0] for i in xrange(1,len(h))]
af = [a[i][1] for i in xrange(1,len(a))]
bf = [b[i][1] for i in xrange(1,len(b))]
cf = [c[i][1] for i in xrange(1,len(c))]
df = [d[i][1] for i in xrange(1,len(d))]
ef = [e[i][1] for i in xrange(1,len(e))]
ff = [f[i][1] for i in xrange(1,len(f))]
gf = [g[i][1] for i in xrange(1,len(g))]
hf = [h[i][1] for i in xrange(1,len(h))]

f, ((axa, axb), (axc, axd), (axg,axh)) = plt.subplots(2, 3, sharex='col', sharey='row')
axa.plot(atime,af, '-', color='blue', linewidth=1)
axa.set_title('(a)')
axb.plot(btime,bf, '-', color='blue', linewidth=1)
axa.set_title('(b)')
axc.plot(ctime,cf, '-', color='blue', linewidth=1)
axa.set_title('(c)')
axd.plot(dtime,df, '-', color='blue', linewidth=1)
axa.set_title('(d)')
axg.plot(gtime,gf, '-', color='blue', linewidth=1)
axa.set_title('(g)')
axh.plot(htime,hf, '-', color='blue', linewidth=1)
axa.set_title('(h)')
plt.xlabel('Time')
plt.ylabel('Drag')
plt.xlim([0,10])
plt.ylim([0,2])
plt.show()
#plt.savefig('%s/osc.pdf' % (caseFolder))
plt.clf()
	
print "\nDone plotting! Plot saved to file "+outFile+"\n"

