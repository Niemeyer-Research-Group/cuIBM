#!/usr/bin/env python
#generates validation plots of a lid driven cavity flow without running it
import csv
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import os
import os.path
cuibmFolder = os.path.expandvars("/scratch/src/cuIBM")

import sys
sys.path.insert(0, cuibmFolder+'/scripts/python')

import readData as rd

# Parse command line options
parser = argparse.ArgumentParser(description="Runs the validation case for flow in a lid-driven cavity for a specified Reynolds number", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--Re", dest="Re", help="Reynolds number", default='100')
args = parser.parse_args()
Re = args.Re

if Re=='100':
	uCol           = '1'
	vCol           = '6'
elif Re=='1000':
	uCol           = '2'
	vCol           = '7'
elif Re=='5000':
	uCol           = '3'
	vCol           = '8' 
elif Re=='10000':
	uCol           = '4'
	vCol           = '9'
else:
	print "Unavailable option for Reynolds number. Choose 100, 1000 or 10000."
	sys.exit()

validationData = '/cavity-GGS82.txt'
caseFolder     = cuibmFolder + '/validation/lidDrivenCavity/Re' + Re
validationData = cuibmFolder + '/validation-data' + validationData
execPath       = cuibmFolder + '/bin/cuIBM'

nt, _, _, _ = rd.readSimulationParameters(caseFolder)
nx, ny, dx, dy, _, yu, xv, _ = rd.readGridData(caseFolder)
u, v = rd.readVelocityData(caseFolder, nt, nx, ny, dx, dy)

u2=[]
x2=[]
v2=[]
y2=[]

x2.append(0)
u2.append(0)
for j in range(len(yu)):
	u2.append(u[j*(nx-1)+nx/2-1])
	x2.append(yu[j])
x2.append(1)
u2.append(1)

y2.append(0)
v2.append(0)
for i in range(len(xv)):
	v2.append(v[(ny/2-1)*nx + i])
	y2.append(xv[i])
y2.append(1)
v2.append(0)
	
validation_x=[]
validation_y=[]
validation_v=[]
validation_u=[]
with open(validationData) as f:
	reader = csv.reader(f,delimiter='\t',quotechar='|')
	i=0
	for row in reader:
		validation_x.append(float(row[0]))
		validation_y.append(float(row[5]))
		validation_u.append(float(row[int(uCol)]))
		validation_v.append(float(row[int(vCol)]))
		i+=1

print "-"*80
print "Plotting the centerline velocities for flow in a lid-driven cavity at Reynolds number %s" % Re
print "-"*80

uOutFile       = caseFolder + '/uRe' + Re + '.pdf'
vOutFile       = caseFolder + '/vRe' + Re + '.pdf'

plt.plot(validation_x, validation_u, 'o', color = 'red', markersize = 8, label = 'Ghia et al, 1982')
plt.plot(x2,u2, '-', color='blue', linewidth=2, label='Present Work')
plt.title('Velocity along the vertical centerline Re = {0}'.format(Re))
plt.legend(loc='upper right', numpoints=1, fancybox=True)
plt.xlabel('y')
plt.ylabel('Centerline u-velocity')
plt.xlim([0,1])
plt.ylim([-0.6,1.2])
plt.savefig('%s/u-velocity.pdf' % (caseFolder))
plt.clf()


plt.plot(validation_y, validation_v, 'o', color = 'red', markersize = 8, label = 'Ghia et al, 1982')
plt.plot(y2,v2, '-', color='blue', linewidth=2, label='Present Work')
plt.title('Velocity along the horizontal centerline Re = {0}'.format(Re))
plt.legend(loc='upper right', numpoints=1, fancybox=True)
plt.xlabel('x')
plt.ylabel('Centerline v-velocity')
plt.xlim([0,1])
plt.ylim([-0.6,1.2])
plt.savefig('%s/v-velocity.pdf' % (caseFolder))
plt.clf()

print "\nDone plotting! Plots saved in folder " + caseFolder + "\n"
