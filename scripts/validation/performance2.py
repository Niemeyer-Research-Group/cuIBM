#!/usr/bin/env python
#This script plots oscillating cylinder performance
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
	name = '/scratch/src/cuIBM/validation/error/oscflow/'
	fileid = '/forces'

	h=[0.0625, 0.05, 0.03125]
	time_em=[454,581,1130]
	time_ex=[53,61,72]
	plt.plot(h,time_ex,'-s',label='External')
	plt.plot(h,time_em,'-^',label='Embedded')
	plt.xlabel('Uniform section grid spacing')
	plt.ylabel('Computational Time(s)')
	plt.title('Drag for flow over in-line oscillating cylinder')
	plt.legend(loc='upper right', numpoints=1, fancybox=True)
	plt.savefig('/scratch/src/cuIBM/validation/error/oscflow2/performance_oscflow2.pdf')
	plt.clf()

	print "\nDone plotting!\n"



if __name__ == "__main__":
	main()

