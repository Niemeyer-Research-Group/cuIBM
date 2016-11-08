#!/usr/bin/env python

#import csv
#import argparse
import numpy as np
from numpy import genfromtxt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from scipy import signal
#import os
#import os.path
#import sys

def main():
	print "-"*80
	print "Plotting error order"
	print "-"*80
	caseFolder = '/scratch/src/cuIBM/validation/error/'
	name = '/scratch/src/cuIBM/validation/error/cylinder/'
	fileid = '/forces'
	grid_cylinder_750ti = [170,208,450,890,135,208,450,890,162,352,372,462]
	time_cylinder_750ti = [117,236,1017,3820,81,127,364,1487,390,894,1099,2880] #short run
	time_cylinder_k20 = [43,61,151,563,20,26,51,174,122,206,281,762]
	grid_osc_750ti = [732, 544, 352, 228, 732, 544, 352, 228]
	time_osc_750ti = [422, 211, 112, 56,  425, 214, 112,  56]
	
	for i in xrange(len(grid_cylinder_750ti)):
		grid_cylinder_750ti[i] = grid_cylinder_750ti[i]**2	

	plt.loglog(grid_cylinder_750ti[0:4], time_cylinder_750ti[0:4] ,'-.o', label='Fadlun-750ti')
	plt.loglog(grid_cylinder_750ti[4:8], time_cylinder_750ti[4:8] ,'-.s', label='External-750ti')
	plt.loglog(grid_cylinder_750ti[8:12],time_cylinder_750ti[8:12],'-.^', label='Embedded-750ti')
	plt.loglog(grid_cylinder_750ti[0:4], time_cylinder_k20[0:4]   ,'-o' , label='Fadlun-k20')
	plt.loglog(grid_cylinder_750ti[4:8], time_cylinder_k20[4:8]   ,'-s' , label='External-k20')
	plt.loglog(grid_cylinder_750ti[8:12],time_cylinder_k20[8:12]  ,'-^' , label='Embedded-k20')
	plt.xlabel('Grid Size')
	plt.ylabel('Run Time')
	plt.title('IBM performance on impulsively started cylinder with Re=40')
	plt.legend(loc='upper left', numpoints=1, fancybox=True)
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/cylinder_performance.pdf')
	plt.clf()

	plt.loglog(grid_osc_750ti[0:4],time_osc_750ti[0:4],'-s', label='External')
	plt.loglog(grid_osc_750ti[4:8],time_osc_750ti[4:8],'-^', label='Embedded')
	plt.xlabel('Grid Size')
	plt.ylabel('Run Time')
	plt.title('IBM performance on oscillating cylinder in flow')
	plt.legend(loc='lower right', numpoints=1, fancybox=True)
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/osc_performance.pdf')
	plt.clf()

if __name__ == "__main__":
	main()

