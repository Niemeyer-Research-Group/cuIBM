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
	print "Making validaiton plot for an oscillating cylinder in flow."
	print "-"*80
	caseFolder = '/scratch/src/cuIBM/validation/osc/VIV'
	name = '/scratch/src/cuIBM/validation/osc/VIV/Ured'
	fileid = '/midPosition'
	data = genfromtxt(name+'4' + fileid,dtype=float,delimiter='\t',skip_header=1)
	t = [data[i][0]*0.005 for i in xrange(1,len(data))]
	y = [data[i][2] for i in xrange(1,len(data))]
	
	for i in xrange(1,len(data)/10):	
		plt.plot(t[0:i*10],y[0:i*10],'k')
		plt.plot(t[i*10],y[i*10],'ro')
		plt.title('t='+str(t[i*10]))
		plt.xlabel('Time')
		plt.ylabel('Height')
		plt.xlim([0,60])
		plt.ylim([-0.6,0.6])
		#plt.savefig('/scratch/src/cuIBM/validation/osc/VIV/pretty/VIV'+str(i*10)+'.jpg')
		plt.savefig('/scratch/src/cuIBM/validation/osc/VIV/pretty/VIV{0:0=5d}.png'.format(i*10))
		plt.clf()
		print i*10, t[i*10], y[i*10]
	
	print "\nDone plotting!\n"

def plot_viv(x,data, peaks, peaks_y):
	plt.plot(x,data,'-')
	plt.plot(peaks,peaks_y,'o')
	plt.title('VIV Position Plot')
	plt.xlabel('Time Step')
	plt.ylabel('Y Position')
	plt.savefig('/scratch/src/cuIBM/validation/osc/VIV/VIV_position.pdf')
	plt.clf()

	print "\nDone plotting!\n"

#run
main()

