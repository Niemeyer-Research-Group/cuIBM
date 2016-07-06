#!/usr/bin/env python

import argparse
import os
import os.path
import sys
import csv
import matplotlib
from matplotlib import pyplot as plt
import numpy

cuibmFolder = os.path.expandvars("/scratch/src/cuIBM")

validationData = '/osc_Re100_KC5_Dutsch.txt'

execPath       = cuibmFolder + '/bin/cuIBM'
caseFolder     = cuibmFolder + '/validation/osc/static'
validationData = cuibmFolder + '/validation-data' + validationData

print "\n"+"-"*100
print "Running flow induced by Oscillating Cylinder\n"
print "-"*100+"\n"

runCommand = "%s -caseFolder %s" % (execPath, caseFolder)
print runCommand+"\n"
os.system(runCommand)

print "-"*100
print "Plotting drag\n"
print "-"*100

experiment = numpy.genfromtxt(validationData,delimiter='\t')
cfd = numpy.genfromtxt(caseFolder+'/forces',delimiter='\t')

#mult by 5 because the diameter is 0.2
plt.plot(zip(*cfd)[0],[i*5 for i in zip(*cfd)[1]],'-',color='blue',linewidth=2,label='PresentWork')
plt.plot(zip(*experiment)[0],zip(*experiment)[1],'o', color = 'red', markersize = 8, label = 'Dutsch et al')
plt.title('Drag for inline oscillating cylinder Re 100 KC 5')
plt.legend(loc='lower right',numpoints=1, fancybox=True)
plt.xlabel('t/T')
plt.ylabel('Fd')
plt.ylim([-2,2])
plt.savefig('%s/drag.png' % (caseFolder))
plt.clf()


print "\nDone plotting!\n"
