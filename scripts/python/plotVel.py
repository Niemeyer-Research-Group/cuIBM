#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
​
val_data = np.genfromtxt('cavity-GGS82.txt')
u = np.genfromtxt('u.txt')
v = np.genfromtxt('v.txt')
​
plt.plot(val_data[:,0], val_data[:,1], 'o', color='red',
         markersize=8, label='Ghia et al., 1982'
         )
plt.plot(u[:,0], u[:,1], '-', color='blue', linewidth=2, label='Present Work')
​
#plt.title('Velocity along the vertical centerline (Re=100)')
plt.xlabel('y-coordinate')
plt.ylabel('Centerline u-velocity')
plt.legend(loc='upper left', numpoints=1, fancybox=True)
​
pp = PdfPages('uRe100.pdf')
pp.savefig()
pp.close()
​
plt.clf()
​
plt.plot(val_data[:,5], val_data[:,6], 'o', color='red',
         markersize=8, label='Ghia et al., 1982'
         )
plt.plot(v[:,0], v[:,1], '-', color='blue', linewidth=2, label='Present Work')
​
#plt.title('Velocity along the horizontal centerline (Re=100)')
plt.xlabel('x-coordinate')
plt.ylabel('Centerline v-velocity')
plt.legend(loc='lower left', numpoints=1, fancybox=True)
​
pp = PdfPages('vRe100.pdf')
pp.savefig()
pp.close()