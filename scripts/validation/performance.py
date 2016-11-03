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
	#
	gridcy = [170,208,450,890,135,208,450,890,162,352,372,462]
	#timecy = [286,521,2138,10484,327,505,1426,5773,773,1749,2139,5488] #long run
	timecy = [117,236,1017,3820,81,127,364,1487,390,894,1099,2880] #short run
	gridosc = [732, 544, 352, 228, 732, 544, 352, 228]
	timeosc = [422, 211, 112, 56,  425, 214, 112,  56]
	
	for i in xrange(len(gridcy)):
		gridcy[i] = gridcy[i]**2	

	plt.loglog(gridcy[0:4],timecy[0:4],'-o', label='Fadlun')
	plt.loglog(gridcy[4:8],timecy[4:8],'-s',label='External')
	plt.loglog(gridcy[8:12],timecy[8:12],'-^',label='Embedded')
	plt.xlabel('Grid Size')
	plt.ylabel('Run Time')
	plt.title('IBM performance on impulsively started cylinder with Re=40')
	plt.legend(loc='lower right', numpoints=1, fancybox=True)
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/cylinder_performance.pdf')
	plt.clf()

	plt.loglog(gridosc[0:4],timeosc[0:4],'-s', label='External')
	plt.loglog(gridosc[4:8],timeosc[4:8],'-^', label='Embedded')
	plt.xlabel('Grid Size')
	plt.ylabel('Run Time')
	plt.title('IBM performance on oscillating cylinder in flow')
	plt.legend(loc='lower right', numpoints=1, fancybox=True)
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/osc_performance.pdf')
	plt.clf()

if __name__ == "__main__":
	main()

