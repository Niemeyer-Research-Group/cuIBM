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
	name = '/scratch/src/cuIBM/validation/error/oscflow/'
	fileid = '/forces'
	typeid = ['external', 'embedded']
	val=[0]*2
	place = 0
	for methodtype in typeid:
		d1 = genfromtxt(name + methodtype + '015625' + fileid,dtype=float,delimiter='\t',skip_header=1)
		d2 = genfromtxt(name + methodtype + '02' + fileid,dtype=float,delimiter='\t',skip_header=1)
		d3 = genfromtxt(name + methodtype + '03125' + fileid,dtype=float,delimiter='\t',skip_header=1)
		d4 = genfromtxt(name + methodtype + '0625' + fileid,dtype=float,delimiter='\t',skip_header=1)

		t1 = [d1[i][0] for i in xrange(1,len(d1))]
		t2 = [d2[i][0] for i in xrange(1,len(d2))]
		t3 = [d3[i][0] for i in xrange(1,len(d3))]
		t4 = [d4[i][0] for i in xrange(1,len(d4))]

		y1 = [d1[i][1] for i in xrange(1,len(d1))]
		y2  = [d2[i][1] for i in xrange(1,len(d2))]
		y3  = [d3[i][1] for i in xrange(1,len(d3))]
		y4 = [d4[i][1] for i in xrange(1,len(d4))]

		error = [0]*3
		error[0] = find_error(t1,t4,y1,y4)
		error[1] = find_error(t1,t3,y1,y3)
		error[2] = find_error(t1,t2,y1,y2)
		
		val[place] = error
		place +=1

	h=[0.0625, 0.03125, 0.02]
	plt.loglog(h,val[0],'-s',label='External')
	plt.loglog(h,val[1],'-^',label='Embedded')
	plt.xlabel('Grid Spacing')
	plt.ylabel('Error')
	plt.title('Flow over in-line oscillating cylinder')
	plt.legend(loc='upper left', numpoints=1, fancybox=True)
	plt.savefig('/scratch/src/cuIBM/validation/error/oscflow/error_oscflow.pdf')
	plt.clf()

	print "\nDone plotting!\n"

def find_error(tfine, tcourse, yfine, ycourse):
	error=[0]*len(tcourse)
	for i in xrange(len(tcourse)-1):
		j = 0
		if tfine[j]>10 or j>=len(tfine):
			pass
		else:
			while tfine[j] <= tcourse[i]:
				j+=1
			if tfine[j] == tcourse[i]:
				error[i]=abs(yfine[j]-ycourse[i])/abs(yfine[j])
			else:
				yf = (yfine[j]-yfine[j-1])/(tfine[j]-tfine[j-1])*(tcourse[i]-tfine[j-1])+yfine[j]
				error[i]=abs(yfine[j]-ycourse[i])/abs(yfine[j])
	errorsum = sum(error)
	return errorsum/len(error)

if __name__ == "__main__":
	main()

