#!/usr/bin/env python
#This script makes the validation plot for an in-line oscilatting cylinder driven flow with Re100 KC5

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

print "-"*100
print "Plotting drag\n"
print "-"*100

experiment = numpy.genfromtxt(validationData,delimiter='\t')
external = numpy.genfromtxt(caseFolder+'/luo/forces_init',delimiter='\t')
embedded = numpy.genfromtxt(caseFolder+'/osc/forces_init',delimiter=',')


#mult by 5 because the diameter is 0.2
#Initial external
plt.plot(zip(*external)[0],[i*5 for i in zip(*external)[1]],'-',color='blue',linewidth=2,label='PresentWork')
plt.plot(zip(*experiment)[0],zip(*experiment)[1],'o', color = 'red', markersize = 8, label = 'Dutsch et al')
plt.title('First period drag of oscillating cylinder in no flow Re 100, KC 5')
plt.legend(loc='lower right',numpoints=1, fancybox=True)
plt.xlabel('t/T')
plt.ylabel('Fd')
plt.ylim([-2,2])
plt.savefig('%s/staticexinit.pdf' % (caseFolder))
plt.clf()

#emb init
plt.plot(zip(*embedded)[0],[i*5 for i in zip(*embedded)[1]],'-',color='blue',linewidth=2,label='PresentWork')
plt.plot(zip(*experiment)[0],zip(*experiment)[1],'o', color = 'red', markersize = 8, label = 'Dutsch et al')
plt.title('First period drag of oscillating cylinder in no flow Re 100, KC 5')
plt.legend(loc='lower right',numpoints=1, fancybox=True)
plt.xlabel('t/T')
plt.ylabel('Fd')
plt.ylim([-2,2])
plt.savefig('%s/staticeminit.pdf' % (caseFolder))
plt.clf()

external = numpy.genfromtxt(caseFolder+'/luo/forces_ss',delimiter=',')
embedded = numpy.genfromtxt(caseFolder+'/osc/forces_ss',delimiter=',')

#ss external
plt.plot(zip(*external)[0],[i*5 for i in zip(*external)[1]],'-',color='blue',linewidth=2,label='PresentWork')
plt.plot(zip(*experiment)[0],zip(*experiment)[1],'o', color = 'red', markersize = 8, label = 'Dutsch et al')
plt.title('Steady state period of oscillating cylinder in no flow Re 100, KC 5')
plt.legend(loc='lower right',numpoints=1, fancybox=True)
plt.xlabel('t/T')
plt.ylabel('Fd')
plt.ylim([-2,2])
plt.savefig('%s/staticexss.pdf' % (caseFolder))
plt.clf()

#ss emb
plt.plot(zip(*embedded)[0],[i*5 for i in zip(*embedded)[1]],'-',color='blue',linewidth=2,label='PresentWork')
plt.plot(zip(*experiment)[0],zip(*experiment)[1],'o', color = 'red', markersize = 8, label = 'Dutsch et al')
plt.title('Steady state period of oscillating cylinder in no flow Re 100, KC 5')
plt.legend(loc='lower right',numpoints=1, fancybox=True)
plt.xlabel('t/T')
plt.ylabel('Fd')
plt.ylim([-2,2])
plt.savefig('%s/staticemss.pdf' % (caseFolder))
plt.clf()


print "\nDone plotting!\n"
