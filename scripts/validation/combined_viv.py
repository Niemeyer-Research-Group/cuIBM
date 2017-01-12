#!/usr/bin/env python
#This script makes the plot for all four VIVs results on the same plot
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

	ys = [0.0714, 0.5616, 0.5234, 0.451, 0.371, 0.0696]
	yc = [0.0714, 0.5286, 0.4894, 0.435, 0.381, 0.0696]
	em_lc = [0.07161494006079072,0.5717792720930622,0.5344207235887237, 0.46769951947874355, 0.38874337573979395, 0.14012940247913452]
	em_sc = [0.07406970331657958, 0.5443322122816592, 0.48226014304124076, 0.3941003346010647, 0.310016601470416, 0.0881615228777407]
	ex_lc = [0.07041599084467354, 0.5675364513644794, 0.5202419621998599, 0.4533280558268404, 0.38804988964561826, 0.12004663709542096]
	ex_sc = [0.0742857142857144, 0.5599999999999999, 0.5069387755102042, 0.4342857142857142, 0.35346938775510206, 0.10448979591836749]
	x = [3, 4, 5, 6, 7, 8]

	plt.plot(x,ys,'sk',label='Ann & Kallindens(2006)')
	plt.plot(x,yc,'ok',label='Borazjani et al.(2008)')
	plt.plot(x,ex_lc,'^r',label='External loose')
	plt.plot(x,ex_sc,'dr',label='External strong')
	plt.plot(x,em_lc,'xb',label='Embedded loose')
	plt.plot(x,em_sc,'+b',label='Embedded strong')
	plt.xlabel('Ured')
	plt.ylabel('Maximum Amplitude')
	plt.xlim([2,9])
	plt.ylim([0,0.6])
	plt.legend(loc='best', fancybox=True)
	plt.savefig('/scratch/src/cuIBM/validation/osc/VIV/VIV_combine.pdf')
	plt.clf()

#run
main()

