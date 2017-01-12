#!/usr/bin/env python
#This script calculates and plots the order of accuracy for the impuslively started cylinder
#import csv
#import argparse
import numpy as np
from numpy import genfromtxt
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
	typeid = ['fadlun', 'external', 'embedded']
	#typeid = ['fadlun']
	for methodtype in typeid:
		error = [0]*3
		try:
			d005 = genfromtxt(name + methodtype + '4' + fileid,dtype=float,delimiter='\t',skip_header=1)
			t005 = [d005[i][0] for i in xrange(1,len(d005))]
			y005 = [d005[i][1] for i in xrange(1,len(d005))]
		except:
			print methodtype+'4 failed to run'
		
		try:			
			d01 = genfromtxt(name + methodtype + '3' + fileid,dtype=float,delimiter='\t',skip_header=1)
			t01  = [d01[i][0] for i in xrange(1,len(d01))]
			y01  = [d01[i][1] for i in xrange(1,len(d01))]
		except:
			print methodtype+'3 failed to run'

		try:
			d02 = genfromtxt(name + methodtype + '2' + fileid,dtype=float,delimiter='\t',skip_header=1)
			t02  = [d02[i][0] for i in xrange(1,len(d02))]
			y02  = [d02[i][1] for i in xrange(1,len(d02))]
		except:
			print methodtype+'2 failed to run'

		try:
			d025 = genfromtxt(name + methodtype + '1' + fileid,dtype=float,delimiter='\t',skip_header=1)
			t025 = [d025[i][0] for i in xrange(1,len(d025))]
			y025 = [d025[i][1] for i in xrange(1,len(d025))]
		except:
			print methodtype+'1 failed to run'
		
		error[0],errorout = find_error(t005,t025,y005,y025)
		ploterror(t025, errorout, methodtype+'_error_025')
		error[1],errorout = find_error(t005,t02,y005,y02)
		ploterror(t02,errorout,methodtype+'_error_02')
		error[2],errorout = find_error(t005,t01,y005,y01)
		ploterror(t01,errorout,methodtype+'_error_01')

		print error
		if methodtype == 'fadlun':
			h=[0.03, 0.02, 0.01]
		elif methodtype == 'external':
			h=[0.05, 0.02, 0.01]
		elif methodtype == 'embedded':
			h = [0.0625,0.03125,0.02]
		else:
			print 'No solver type found'
		"""#plot all validation cases		
		h.append(0.005)
		for i in xrange(len(h)):
			plotme(methodtype,str(i+1),str(h[i])[2:])"""

		plt.loglog(h,error,'-o')
		plt.xlabel('Grid Spacing')
		plt.ylabel('Error')
		plt.title('Force based error vs grid spacing')
		#plt.xlim([2,9])
		#plt.ylim([0,0.6])
		#plt.legend()
		plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/error_'+methodtype+'.pdf')
		plt.clf()

def main2():
	caseFolder = '/scratch/src/cuIBM/validation/error/'
	name = '/scratch/src/cuIBM/validation/error/cylinder/'
	
	#typeid = ['fadlun', 'external', 'embedded']
	typeid =['fadlun', 'external']
	for methodtype in typeid:
		y005 = genfromtxt(name+methodtype+'4/xu',dtype=float,delimiter='\t',skip_header=0)
		x005 = genfromtxt(name+methodtype+'4/xu',dtype=float,delimiter='\t',skip_header=0)
		u005 = genfromtxt(name+methodtype+'4/output/5000u.csv',dtype=float,delimiter='\t',skip_header=1)
		tags005 = genfromtxt(name+methodtype+'4/output/5000ghostu.csv',dtype=int,delimiter='\t',skip_header=1)

		y01 = genfromtxt(name+methodtype+'3/xu',dtype=float,delimiter='\t',skip_header=0)
		x01 = genfromtxt(name+methodtype+'3/xu',dtype=float,delimiter='\t',skip_header=0)
		u01 = genfromtxt(name+methodtype+'3/output/5000u.csv',dtype=float,delimiter='\t',skip_header=1)
		tags01 = genfromtxt(name+methodtype+'3/output/5000ghostu.csv',dtype=int,delimiter='\t',skip_header=1)

		y02 = genfromtxt(name+methodtype+'2/xu',dtype=float,delimiter='\t',skip_header=0)
		x02 = genfromtxt(name+methodtype+'2/xu',dtype=float,delimiter='\t',skip_header=0)
		u02 = genfromtxt(name+methodtype+'2/output/5000u.csv',dtype=float,delimiter='\t',skip_header=1)
		tags02 = genfromtxt(name+methodtype+'2/output/5000ghostu.csv',dtype=int,delimiter='\t',skip_header=1)

		y025 = genfromtxt(name+methodtype+'1/xu',dtype=float,delimiter='\t',skip_header=0)
		x025 = genfromtxt(name+methodtype+'1/xu',dtype=float,delimiter='\t',skip_header=0)
		u025 = genfromtxt(name+methodtype+'1/output/5000u.csv',dtype=float,delimiter='\t',skip_header=1)
		tags025 = genfromtxt(name+methodtype+'1/output/5000ghostu.csv',dtype=int,delimiter='\t',skip_header=1)

		error = [0]*3
		if methodtype == 'fadlun':
			h=[0.03, 0.02, 0.01]
		elif methodtype == 'external':
			h=[0.05, 0.02, 0.01]
		elif methodtype == 'embedded':
			h = [0.0625,0.03125,0.02]
		else:
			print 'No solver type found'

		error[0] = error_2(y005,y025,x005,x025,u005,u025,tags025)
		error[1] = error_2(y005,y02,x005,x02,u005,u02,tags02)
		error[2] = error_2(y005,y01,x005,x01,u005,u01,tags01)
		print h, error
		plt.loglog(h,error,'-o')
		plt.xlabel('Grid Spacing')
		plt.ylabel('Error')
		plt.title('Field based error vs grid spacing')
		plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/error2_'+methodtype+'.pdf')
		plt.clf()
	
def find_error(tfine, tcourse, yfine, ycourse):
	error=[0]*len(tcourse)
	for i in xrange(len(tcourse)-1):
		j = 0
		while tfine[j] <= tcourse[i]:
			j+=1
		if tfine[j] == tcourse[i]:
			error[i]=abs(yfine[j]-ycourse[i])/abs(yfine[j])	
		else:
			yf = (yfine[j]-yfine[j-1])/(tfine[j]-tfine[j-1])*(tcourse[i]-tfine[j-1])+yfine[j]
			error[i]=abs(yf-ycourse[i])/abs(yf)	
	#cut off initial
	count = 0
	for i in xrange(len(tcourse)-1):
		if ycourse[i] >= 3.0:
			error[i] = 0
			count += 1
	errorsum = sum(error)
	return errorsum/(len(error)-count), error

def error_2(yfine,ycoarse,xfine,xcoarse,ufine,ucoarse,tags):
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
	"""X, Y = np.meshgrid(xcoarse,ycoarse)
	fig=plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(X,Y,error)
	#plt.show()
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/test.pdf')
	plt.clf()"""
	return errorsum/(len(xcoarse)*len(ycoarse)-count)


def plotme(typeid,num,h):
	yMax           = '6'

	cuibmFolder    = '/scratch/src/cuIBM'
	caseFolder     = cuibmFolder + '/validation/error/cylinder/'+typeid+num
	validationData = cuibmFolder + '/validation-data/cylinderRe40-KL95.txt'

	print caseFolder+'/forces'
	my_data = genfromtxt(caseFolder+'/forces',dtype=float,delimiter='\t')
	time = [my_data[i][0] for i in xrange(1,len(my_data))]
	force = [my_data[i][1]*2 for i in xrange(1,len(my_data))]

	validation_data = genfromtxt(validationData, dtype=float, delimiter='\t')
	validation_time=[validation_data[i][0]*0.5 for i in xrange(1,len(validation_data))]
	validation_force=[validation_data[i][1] for i in xrange(1,len(validation_data))]

	plt.plot(validation_time,validation_force, 'o', color = 'red', markersize = 8, label = 'Koumoutsakos and Leonard, 1995')
	plt.plot(time,force, '-', color='blue', linewidth=2, label='Present Work')
	plt.title('Flow over an impulsively started cylinder at Reynolds number 40')
	plt.legend(loc='upper right', numpoints=1, fancybox=True)
	plt.xlabel('Non-dimensional time')
	plt.ylabel('Drag Coefficient')
	plt.xlim([0,3])
	plt.ylim([0,int(yMax)])
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/'+typeid+h+'.pdf')
	plt.clf()

def ploterror(time, error,name):
	if name[0:8] == 'embedded':
		plt.plot(time[0:2500], error[0:2500])
	else:
		plt.plot(time,error)
	plt.xlabel('time')
	plt.ylabel('error')
	plt.axis([0,2.5,0,0.05])
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/'+name+'.pdf')
	plt.clf()
	print 'printing' + name

if __name__ == "__main__":
	main()
