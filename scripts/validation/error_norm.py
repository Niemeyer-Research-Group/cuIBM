#!/usr/bin/env python

#import csv
#import argparse
import numpy as np
from numpy import genfromtxt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import matplotlib.markers as mmarkers
import matplotlib.lines as mlines
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal
from math import log
from math import sqrt
#import os
#import os.path
#import sys

def main():
	caseFolder = '/scratch/src/cuIBM/validation/error/'
	name = '/scratch/src/cuIBM/validation/error/cylinder/'

	typeid = ['fadlun', 'external', 'embedded']
	t = '5000'
	for methodtype in typeid:	
		y4 = genfromtxt(name+methodtype+'4/xu',dtype=float,delimiter='\t',skip_header=0)
		x4 = genfromtxt(name+methodtype+'4/xu',dtype=float,delimiter='\t',skip_header=0)
		u4 = genfromtxt(name+methodtype+'4/output/'+t+'u.csv',dtype=float,delimiter='\t',skip_header=1)
		tags4 = genfromtxt(name+methodtype+'4/output/'+t+'ghostu.csv',dtype=int,delimiter='\t',skip_header=1)

		y3 = genfromtxt(name+methodtype+'3/xu',dtype=float,delimiter='\t',skip_header=0)
		x3 = genfromtxt(name+methodtype+'3/xu',dtype=float,delimiter='\t',skip_header=0)
		u3 = genfromtxt(name+methodtype+'3/output/'+t+'u.csv',dtype=float,delimiter='\t',skip_header=1)
		tags3 = genfromtxt(name+methodtype+'3/output/'+t+'ghostu.csv',dtype=int,delimiter='\t',skip_header=1)

		y2 = genfromtxt(name+methodtype+'2/xu',dtype=float,delimiter='\t',skip_header=0)
		x2 = genfromtxt(name+methodtype+'2/xu',dtype=float,delimiter='\t',skip_header=0)
		u2 = genfromtxt(name+methodtype+'2/output/'+t+'u.csv',dtype=float,delimiter='\t',skip_header=1)
		tags2 = genfromtxt(name+methodtype+'2/output/'+t+'ghostu.csv',dtype=int,delimiter='\t',skip_header=1)

		y1 = genfromtxt(name+methodtype+'1/xu',dtype=float,delimiter='\t',skip_header=0)
		x1 = genfromtxt(name+methodtype+'1/xu',dtype=float,delimiter='\t',skip_header=0)
		u1 = genfromtxt(name+methodtype+'1/output/'+t+'u.csv',dtype=float,delimiter='\t',skip_header=1)
		tags1 = genfromtxt(name+methodtype+'1/output/'+t+'ghostu.csv',dtype=int,delimiter='\t',skip_header=1)

		L1 = [0]*3
		L2 = [0]*3
		Linf = [0]*3

		L1[0], L2[0], Linf[0] = error_norm(y4,y1,x4,x1,u4,u1,tags1)
		L1[1], L2[1], Linf[1] = error_norm(y4,y2,x4,x2,u4,u2,tags2)
		L1[2], L2[2], Linf[2] = error_norm(y4,y3,x4,x3,u4,u3,tags3)

		print methodtype
		print methodtype, "coarse L1", L1[0]
		print methodtype, "coarse L2", L2[0]
		print methodtype, "coarse Linf", Linf[0]
		print methodtype, "medium L1", L1[1]
		print methodtype, "medium L2", L2[1]
		print methodtype, "medium Linf", Linf[1]
		print methodtype, "fine L1", L1[2]
		print methodtype, "fine L2", L2[2]
		print methodtype, "fine Linf", Linf[2]

def error_norm(yfine,ycoarse,xfine,xcoarse,ufine,ucoarse,tags):
	L1 = np.zeros((len(xcoarse),len(ycoarse)))
	L2 = np.zeros((len(xcoarse),len(ycoarse)))
	uf = 0.0
	count = 0
	for i in xrange(1,len(xcoarse)-1):
		for j in xrange(1,len(ycoarse)-1):
			#interp fine to coarse location
			m=0
			n=0
			while xfine[m]<=xcoarse[i]:
				m+=1
			try:			
				while yfine[n]<=ycoarse[j]:
					n+=1
			except:
				print n, len(yfine)
				print j, len(ycoarse)
				print yfine[n-1], ycoarse[j]
			uf = 1.0/(xfine[m]-xfine[m-1])/(yfine[n]-yfine[n-1]) * (ufine[m-1][n-1]*(xfine[m]-xcoarse[i])*(yfine[n]-ycoarse[j]) + ufine[m][n-1]*(xcoarse[i]-xfine[m-1])*(yfine[n]-ycoarse[j]) + ufine[m-1][n]*(xfine[m]-xcoarse[i])*(ycoarse[j]-yfine[n-1]) + ufine[m][n]*(xcoarse[i]-xfine[m-1])*(ycoarse[j]-yfine[n-1]))
			if tags[i][j] > -1 or tags[i][j+1] > -1 or tags[i][j-1] > -1 or tags[i+1][j] > -1 or tags[i-1][j] > -1 or tags[i][j] == 0 or uf == 0:
				L1[i][j] = 0
				count += 1
			else:
				L1[i][j]=abs(uf-ucoarse[i][j])
				L2[i][j]=L1[i][j]**2
			if L1[i][j] > 5:
				L1[i][j] = 0
				L2[i][j] = 0
				count +=1
	L1_error = L1.sum()
	L2_error = sqrt(L2.sum())
	Linf_error = L1.max()
	return L1_error, L2_error, Linf_error

def fancyPlot():
	circle_ = mlines.Line2D([],[],color='green', linestyle='none', marker ='o', label='Modified Fadlun')
	square_ = mlines.Line2D([],[],color='green', linestyle='none', marker ='s', label='External')
	triangle_ = mlines.Line2D([],[],color='green', linestyle='none', marker ='^', label='Embedded')
	black_solid = mlines.Line2D([],[],color='black', label='L1 error norm')
	red_dash = mlines.Line2D([],[],color='red', linestyle='dashed', label='L2 error norm')
	blue_dot = mlines.Line2D([],[],color='blue', linestyle='dotted', label='Linf error norm')

	#plot error norm
	fadlun_h = [28900, 43264, 202500]
	external_h = [18225, 43264, 202500]
	embedded_h = [26244, 123904, 138384]
	fadlun_L1 = [111.829,152.065, 225.202]
	external_L1 = [64.013,114.79,154.405]
	embedded_L1 = [83.457,194.41,202.39]
	fadlun_L2 = [1.684,1.664,1.348]
	external_L2 = [1.596,1.3817,0.096345]
	embedded_L2 = [1.6256,1.5437,1.3355]
	fadlun_Linf = [0.1246,0.0784,0.0315]
	external_Linf = [0.1754,0.0732,0.0288]
	embedded_Linf = [0.1909,0.09168,0.0447]

	plt.semilogx(fadlun_h,fadlun_L1,'o-k',label="Modified Fadlun")
	plt.semilogx(external_h,external_L1,'s-k',label="External")
	plt.semilogx(embedded_h,embedded_L1,'^-k',label="Embedded")
	plt.xlabel('Grid size(number of cells)')
	plt.ylabel('L1 error norm')
	plt.title('L1 error norm')
	plt.legend(loc="best")
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/L1_error_norm.pdf', bbox_inches='tight')
	plt.clf()

	plt.semilogx(fadlun_h,fadlun_L2,'o-k',label="Modified Fadlun")
	plt.semilogx(external_h,external_L2,'s-k',label="External")
	plt.semilogx(embedded_h,embedded_L2,'^-k',label="Embedded")
	plt.xlabel('Grid size(number of cells)')
	plt.ylabel('L2 error norm')
	plt.title('L2 error norm')
	plt.legend(loc="best")
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/L2_error_norm.pdf', bbox_inches='tight')
	plt.clf()


	plt.semilogx(fadlun_h,fadlun_Linf,'o-k',label="Modified Fadlun")
	plt.semilogx(external_h,external_Linf,'s-k',label="External")
	plt.semilogx(embedded_h,embedded_Linf,'^-k',label="Embedded")
	plt.xlabel('Grid size(number of cells)')
	plt.ylabel('Linf error norm')
	plt.title('L1inf error norm')
	plt.legend(loc="best")
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/Linf_error_norm.pdf', bbox_inches='tight')
	plt.clf()


if __name__ == "__main__":
	#main()
	fancyPlot()

