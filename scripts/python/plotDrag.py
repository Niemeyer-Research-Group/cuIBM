#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

val_data = np.genfromtxt('/scratch/src/cuIBM-FSI/validation-data/cylinderRe40-KL95.txt')
force = np.genfromtxt('forces.txt')

plt.plot(0.5*val_data[:,0], val_data[:,1], 'o', color = 'red', markersize = 8, label = 'Koumoutsakos and Leonard, 1995')
plt.plot(force[4:-1,0], 2*force[4:-1,1], '-', color='blue', linewidth=2, label='Present Work')
#plt.title('Flow over impulsively started cyliner re 40')
plt.xlabel('Non-Dimensional Time')
plt.ylabel('Drag Coefficient')
plt.legend(loc='upper right', numpoints=1, fancybox=True)
plt.axis([0,5,0,6])
pp = PdfPages('CylinderDragRe40.pdf')
pp.savefig()
pp.close()