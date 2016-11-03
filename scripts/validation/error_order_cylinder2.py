#!/usr/bin/env python

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

			error[0] = find_error(y4,y1,x4,x1,u4,u1,tags1)
			error[1] = find_error(y4,y2,x4,x2,u4,u2,tags2)
			error[2] = find_error(y4,y3,x4,x3,u4,u3,tags3)
			
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
	
			


def find_error(yfine,ycoarse,xfine,xcoarse,ufine,ucoarse,tags):
	error = np.zeros((len(xcoarse),len(ycoarse)))
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
				error[i][j] = 0
				count += 1
			else:
				error[i][j]=abs(uf-ucoarse[i][j])/abs(uf)
			if error[i][j] >5:
				error[i][j] = 0
				count +=1
	errorsum = sum(sum(error))
	return errorsum/(len(xcoarse)*len(ycoarse)-count)

if __name__ == "__main__":
	main()




