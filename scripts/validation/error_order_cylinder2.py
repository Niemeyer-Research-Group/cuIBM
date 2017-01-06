#!/usr/bin/env python
#This script calculates and plots the order of accuracy for an impulsively started cylinder using velocity data.
#import csv
#import argparse
import numpy as np
from numpy import genfromtxt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal
from math import log
#import os
#import os.path
#import sys

def main():
	caseFolder = '/scratch/src/cuIBM/validation/error/'
	name = '/scratch/src/cuIBM/validation/error/cylinder/'
	
	typeid = ['fadlun', 'external', 'embedded']
	timestep = ['100','200','300','400','500','600','700','800','900','1000']
	#typeid = ['external']
	#timestep = ['100']
	ooa_fadlun = []
	ooa_ex = []
	ooa_em = []
	for methodtype in typeid:
		for t in timestep:
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

			error = [0]*3
			eoa = [0]*3
			if methodtype == 'fadlun':
				h=[0.03, 0.02, 0.01]
			elif methodtype == 'external':
				h=[0.05, 0.02, 0.01]
			elif methodtype == 'embedded':
				h = [0.0625,0.03125,0.015625]
			else:
				print 'No solver type found'

			L1[0], L2[0] = error_norm(y4,y1,x4,x1,u4,u1,tags1)
			L1[1] = error_norm(y4,y2,x4,x2,u4,u2,tags2)
			error[2] = error_norm(y4,y3,x4,x3,u4,u3,tags3)
			
			eoa[0] = log(error[1]/error[0])/log(h[1]/h[0])
			eoa[1] = log(error[2]/error[1])/log(h[2]/h[1])
			eoa[2] = log(error[2]/error[0])/log(h[2]/h[0])
			print "\n"+methodtype, t
			print "error", error
			print "Order of Accuracy", eoa
			
			if methodtype == 'fadlun':
				ooa_fadlun.append(eoa[1])
			elif methodtype == 'external':
				ooa_ex.append(eoa[1])
			elif methodtype == 'embedded':
				ooa_em.append(eoa[0])
			else:
				print 'No solver type found'
			plt.loglog(h,error,'-o')
	print "\nfadlun"
	print ooa_fadlun
	print "\nexternal"
	print ooa_ex
	print "\nembedded"
	print ooa_em

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
	L1_error = sum(sum(L1))
	L2_error = sqrt(sum(sum(L2))
	return L1_error, L2_error

if __name__ == "__main__":
	#main()
	t1 = [i*0.1+0.1 for i in xrange(5)]
	t2 = [i*0.05+0.05 for i in xrange(10)]
	fadlun = [1.6140209093498898, 1.6148363116, 1.6176514595, 1.6147082774, 1.6073691433, 1.593866169, 1.5897889254, 1.4269754258, 1.5622941351, 1.6658890443071641]
	external = [1.7053037843603716, 1.6785034208, 1.6584672088, 1.6672553451, 1.6962016987, 1.722117897, 1.6719717865, 1.6801085127, 1.6763200642, 1.7155542537]
	embedded = [1.5184468141, 1.4529358104, 1.3968597912, 1.4376764196, 1.3463391108, 1.548904431, 1.2795229804, 1.1966260321, 1.2556144474, 1.1567078918761309]
	plt.plot(t2,fadlun,'o-',label='Modified Fadlun')
	plt.plot(t2,external,'s-',label='External')
	plt.plot(t1,embedded[0:5],'^-',label='Embedded')
	plt.xlabel('Time')
	plt.ylabel('Order of accuracy')
	plt.title('Order of accuracy for impulsively started cylinder')
	plt.legend(loc='lower left', numpoints=1, fancybox=True)
	plt.axis([0,0.5,1,2])
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/error_order_2_plt.pdf')
	plt.clf
