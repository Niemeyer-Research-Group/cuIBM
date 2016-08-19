#!/usr/bin/env python

import csv
import argparse
import numpy as np
from numpy import genfromtxt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import os
import os.path
import sys

cuibmFolder = os.path.expandvars("/scratch/src/cuIBM")

# Parse command line options
parser = argparse.ArgumentParser(description="Runs the validation case for impulsively started flow over a circular cylinder for a specified Reynolds number", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--Re", dest="Re", help="Reynolds number", default='40')
args = parser.parse_args()
Re = args.Re

if Re=='40':
	validationData = '/cylinderRe40-KL95.txt'
	yMax           = '6'
elif Re=='200' or Re =='100' or Re=='150' or Re == '75':
	validationData = ''
	yMax           = '6' 
elif Re=='550':
	validationData = '/cylinderRe550-KL95.txt'
	yMax           = '2' 
elif Re=='3000':
	validationData = '/cylinderRe3000-KL95.txt'
	yMax           = '2' 
else:
	print "Unavailable option for Reynolds number. Choose 40, 75, 100, 150, 200, 550 or 3000."
	sys.exit()

caseFolder     = cuibmFolder + '/validation/cylinder/Re' + Re
validationData = cuibmFolder + '/validation-data' + validationData

print "-"*80
print "Plotting the drag coefficient for flow over an impulsively started circular cylinder at Reynolds number %s" % Re
print "-"*80

outFile        = caseFolder + '/cylRe' + Re + 'Drag.pdf'

my_data = genfromtxt('%s/forces' % caseFolder,dtype=float,delimiter='\t')
time = [my_data[i][0] for i in xrange(1,len(my_data))]
force = [my_data[i][1]*2 for i in xrange(1,len(my_data))]

validation_time=[]
validation_force=[]
with open(validationData) as f2:
	reader = csv.reader(f2,delimiter='\t',quotechar='|')
	i=0
	for row in reader:
		validation_time.append(float(row[0])*0.5)
		validation_force.append(float(row[1]))
		i+=1

plt.plot(validation_time,validation_force, 'o', color = 'red', markersize = 8, label = 'Koumoutsakos and Leonard, 1995')
plt.plot(time,force, '-', color='blue', linewidth=2, label='Present Work')
plt.title('Flow over an impulsively started cylinder at Reynolds number {0}'.format(Re))
plt.legend(loc='upper right', numpoints=1, fancybox=True)
plt.xlabel('Non-dimensional time')
plt.ylabel('Drag Coefficient')
plt.xlim([0,3])
plt.ylim([0,int(yMax)])
plt.savefig('%s/CylinderRe%s.pdf' % (caseFolder,Re))
plt.clf()
#	f.write("'%s/forces' u 1:(2*$2) w l lw 2 lc rgb '#3232ff' title 'Present Work', \\\n" % caseFolder)#
#	f.write("'%s' u (0.5*$1):2 w p pt 6 ps 2 lc rgb '#ff3232' title 'Koumoutsakos and Leonard, 1995'\n" % validationData)
	
print "\nDone plotting! Plot saved to file "+outFile+"\n"

