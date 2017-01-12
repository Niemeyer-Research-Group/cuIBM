#!/usr/bin/env python
#this script claculates and plots transient error norms.
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
	L1_error_norm_coarse = []
	L1_error_norm_medium = []
	L1_error_norm_fine = []
	L2_error_norm_coarse = []
	L2_error_norm_medium = []
	L2_error_norm_fine = []
	Linf_error_norm_coarse = []
	Linf_error_norm_medium = []
	Linf_error_norm_fine = []

	typeid = ['fadlun', 'external', 'embedded']
	timestep = ['100','200','300','400','500','600','700','800','900','1000']
	for methodtype in typeid:
		L1_error_norm_coarse = []
		L1_error_norm_medium = []
		L1_error_norm_fine = []
		L2_error_norm_coarse = []
		L2_error_norm_medium = []
		L2_error_norm_fine = []
		Linf_error_norm_coarse = []
		Linf_error_norm_medium = []
		Linf_error_norm_fine = []
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

			L1 = [0]*3
			L2 = [0]*3
			Linf = [0]*3
			if methodtype == 'fadlun':
				h=[0.03, 0.02, 0.01]
			elif methodtype == 'external':
				h=[0.05, 0.02, 0.01]
			elif methodtype == 'embedded':
				h = [0.0625,0.03125,0.015625]
			else:
				print 'No solver type found'

			L1[0], L2[0], Linf[0] = error_norm(y4,y1,x4,x1,u4,u1,tags1)
			L1[1], L2[1], Linf[1] = error_norm(y4,y2,x4,x2,u4,u2,tags2)
			L1[2], L2[2], Linf[2] = error_norm(y4,y3,x4,x3,u4,u3,tags3)

			print methodtype, t
			L1_error_norm_coarse.append(L1[0])
			L1_error_norm_medium.append(L1[1])
			L1_error_norm_fine.append(L1[2])
			L2_error_norm_coarse.append(L2[0])
			L2_error_norm_medium.append(L2[1])
			L2_error_norm_fine.append(L2[2])
			Linf_error_norm_coarse.append(Linf[0])
			Linf_error_norm_medium.append(Linf[1])
			Linf_error_norm_fine.append(Linf[2])
		print methodtype, "L1 coarse", L1_error_norm_coarse
		print methodtype, "L1 medium", L1_error_norm_medium
		print methodtype, "L1 fine", L1_error_norm_fine
		print methodtype, "L2 coarse", L2_error_norm_coarse
		print methodtype, "L2 medium", L2_error_norm_medium
		print methodtype, "L2 fine", L2_error_norm_fine
		print methodtype, "Linf coarse", Linf_error_norm_coarse
		print methodtype, "Linf medium", Linf_error_norm_medium
		print methodtype, "Linf fine", Linf_error_norm_fine

	"""print "\nFadlun, , , , , , , , ,External, , , , , , , , ,Embedded, , , , , , , ,"
	print "L1, , ,L2, , ,L3, , ,L1, , ,L2, , ,L3, , ,L1, , ,L2, , ,L3, ,"
	print "Coarse,Medium,Fine,Coarse,Medium,Fine,Coarse,Medium,Fine,Coarse,Medium,Fine,Coarse,Medium,Fine,Coarse,Medium,Fine,Coarse,Medium,Fine,Coarse,Medium,Fine,Coarse,Medium,Fine"
	for timestep in xrange(10):
		mystr = ""
		for solvertype in range(3):
			mystr+=str(L1_error_norm_coarse[solvertype*9+timestep])+","
			mystr+=str(L1_error_norm_medium[solvertype*9+timestep])+","
			mystr+=str(L1_error_norm_fine[solvertype*9+timestep])+","
			mystr+=str(L2_error_norm_coarse[solvertype*9+timestep])+","
			mystr+=str(L2_error_norm_medium[solvertype*9+timestep])+","
			mystr+=str(L2_error_norm_fine[solvertype*9+timestep])+","
			mystr+=str(Linf_error_norm_coarse[solvertype*9+timestep])+","
			mystr+=str(Linf_error_norm_medium[solvertype*9+timestep])+","
			mystr+=str(Linf_error_norm_fine[solvertype*9+timestep])+","
		print mystr
				
			 
	print "\ncoarse", L1_error_norm_coarse
	print "\nmedium", L1_error_norm_medium
	print "\nfine", L1_error_norm_fine
	print "\nexternal"
	print "\ncoarse",L2_error_norm_coarse
	print "\nmedium",L2_error_norm_medium
	print "\nfine",L2_error_norm_fine
	print "\nembedded"
	print "\ncoarse",Linf_error_norm_coarse
	print "\nmedium",Linf_error_norm_medium
	print "\nfine",Linf_error_norm_fine"""

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
	tc = [i*0.1+0.1 for i in xrange(5)]
	tf = [i*0.05+0.05 for i in xrange(10)]
	circle_ = mlines.Line2D([],[],color='green', linestyle='none', marker ='o', label='Modified Fadlun')
	square_ = mlines.Line2D([],[],color='green', linestyle='none', marker ='s', label='External')
	triangle_ = mlines.Line2D([],[],color='green', linestyle='none', marker ='^', label='Embedded')
	black_solid = mlines.Line2D([],[],color='black', label='Coarse grid')
	red_dash = mlines.Line2D([],[],color='red', linestyle='dashed', label='Medium grid')
	blue_dot = mlines.Line2D([],[],color='blue', linestyle='dotted', label='Fine grid')

	#L1 error norm
	fadlun_coarse = [67.296265509297882, 67.746397279904613, 68.432526803146644, 69.097917564701234, 69.711852297690569, 70.348016282841925, 70.949452688970766, 71.588392071702529, 72.232552114929518, 72.890912914482016]
	external_coarse = [59.610491494672601, 55.770486043447974, 53.878031143345531, 52.749617419797538, 52.03652745731597, 51.569717308254091, 51.437641014348976, 51.453545428849729, 51.527494462439705, 51.712819740403816]
	embedded_coarse = [38.309077658785398, 39.224435128515388, 40.504568488778609, 40.586144215852698, 41.134678976814961]#, 41.857070948085401, 42.696840372922104, 43.682107602726653, 44.771368892569917, 45.969527666426984]
	fadlun_medium = [90.361201294457899, 91.998398937487565, 93.15273186457236, 94.171906745093253, 95.111663335635157, 96.045564375626526, 96.990617248931613, 97.969690710087008, 98.948758756134879, 99.936781848450636]
	external_medium = [94.753106697259454, 95.617915685716937, 96.363353230311986, 97.015438343491269, 97.578934632459621, 98.100784107230965, 98.650865395601357, 99.174385920931527, 99.718622105275486, 100.25712193365493]
	embedded_medium = [57.862937466614085, 60.036084515766959, 63.761696865388885, 60.926249856838943, 60.787997065367989]#, 65.325507751088708, 72.005943955256981, 79.320902692797347, 86.854699383841904, 94.328808742603115]
	fadlun_fine = [125.48708248085832, 127.54720860217121, 129.42732970041942, 132.45886046404195, 136.43955977510743, 140.8153388227891, 145.49107865419126, 150.2308140857796, 154.89773119117342, 159.61346554509868]
	external_fine = [128.64291253333224, 130.01845231202151, 131.01473022991456, 131.87179431205382, 132.52452036186858, 133.23235892328105, 133.80883323929262, 134.56708218553541, 135.15234029381986, 135.77011323594849]
	embedded_fine = [57.573546078806778, 60.74723264401262, 67.111802043461495, 64.888769702880651, 67.193019875526886]#, 75.614225653490792, 85.448552263329987, 95.224441224794219, 104.64100638899572, 113.46364824440037]
	plt.plot(tf,fadlun_coarse,'o-k')
	plt.plot(tf,external_coarse,'s-k')
	plt.plot(tc,embedded_coarse,'^-k')

	plt.plot(tf,fadlun_medium,'o--r')
	plt.plot(tf,external_medium,'s--r')
	plt.plot(tc,embedded_medium,'^--r')

	plt.plot(tf,fadlun_fine,'o:b')
	plt.plot(tf,external_fine,'s:b')
	plt.plot(tc,embedded_fine,'^:b')

	plt.xlabel('Time')
	plt.ylabel('L1 error norm')
	plt.title('L1 error norm')
	plt.legend([circle_,square_,triangle_,black_solid,red_dash,blue_dot],['Modified Fadlun','External','Embedded','Coarse grid','Medium grid','Fine grid'], bbox_to_anchor = (1.05,1.), loc = 2, borderaxespad=0.)
	plt.show()
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/L1_error_norm.pdf', bbox_inches='tight')
	plt.clf()

	#L2 error norm
	fadlun_coarse = [2.1402121439057655, 1.9854659835127535, 1.8773104329016417, 1.8017565703051994, 1.7476349353785527, 1.7074089032787043, 1.676150307406993, 1.651529465497163, 1.6319641654189918, 1.6163005699522757]
	external_coarse = [1.7568359270640101, 1.8122840547120975, 1.800918373267789, 1.7766945307744693, 1.7522478576575649, 1.730558048982929, 1.711798543740224, 1.6958661913664623, 1.6822835478168008, 1.6708329832998114]
	embedded_coarse = [1.5170903400446, 1.5502522047544325, 1.5397654576266595, 1.5207711249220608, 1.5071619834839827]#, 1.4973207850798982, 1.490887062340127, 1.487425183694179, 1.4862193595899662, 1.4869881215683596]
	fadlun_medium = [2.1213614487338304, 1.9103038455794583, 1.7906045782382596, 1.7127680923389772, 1.6576034879987722, 1.616595772428367, 1.585298537120615, 1.5608976404132031, 1.5416344488040095, 1.5262303098013537]
	external_medium = [2.106472259142766, 1.925896180393097, 1.8146612802144462, 1.7380854021893448, 1.6819139464646335, 1.63901733434517, 1.605295662400328, 1.5781828005131162, 1.5560142626166493, 1.537699358871155]
	embedded_medium = [1.5799794285623812, 1.4727568778782925, 1.4120736035381922, 1.3651087071877932, 1.333829818764931]#, 1.315668252483344, 1.3073637375406697, 1.3064406603511152, 1.3113028060009404, 1.3205591882300831]
	fadlun_fine = [1.5627857148369129, 1.3645676985785875, 1.2686720627408516, 1.2094318777768205, 1.1695563073585649, 1.1416194890598892, 1.1218118651008469, 1.1077806056061756, 1.0983186284665978, 1.0921580130685942]
	external_fine = [1.5605919061264413, 1.3862163033728112, 1.2943674512994232, 1.2343627559277248, 1.1916488058505996, 1.1594548911911087, 1.13431849925827, 1.1142374499992647, 1.097819377211876, 1.0842671571274138]
	embedded_fine = [1.2683134118377133, 1.1457551265250319, 1.0844484590002965, 1.0403985890554233, 1.0183481995072254]#, 1.014086847037496, 1.0219319954518051, 1.0378042433306418, 1.0589532361560694, 1.0828550405266408]
	plt.plot(tf,fadlun_coarse,'o-k',label='Modified Fadlun coarse')
	plt.plot(tf,external_coarse,'s-k',label='External coarse')
	plt.plot(tc,embedded_coarse,'^-k',label='Embedded coarse')

	plt.plot(tf,fadlun_medium,'o--r',label='Modified Fadlun medium')
	plt.plot(tf,external_medium,'s--r',label='External medium')
	plt.plot(tc,embedded_medium,'^--r',label='Embedded Fadlun medium')

	plt.plot(tf,fadlun_fine,'o:b',label='Modified Fadlun fine')
	plt.plot(tf,external_fine,'s:b',label='External fine')
	plt.plot(tc,embedded_fine,'^:b',label='Embedded fine')

	plt.xlabel('Time')
	plt.ylabel('L2 error norm')
	plt.title('L2 error norm')
	plt.legend([circle_,square_,triangle_,black_solid,red_dash,blue_dot],['Modified Fadlun','External','Embedded','Coarse grid','Medium grid','Fine grid'], bbox_to_anchor = (1.05,1.), loc = 2, borderaxespad=0.)
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/L2_error_norm.pdf', bbox_inches='tight')
	plt.clf()

	#Linf error norm
	fadlun_coarse = [0.46553109097918444, 0.33015474171164338, 0.27311813801079532, 0.24174994217424939, 0.22113592906707896, 0.20662756360832801, 0.19565915574402559, 0.18713928296067939, 0.18029367771781124, 0.17462234001542104]
	external_coarse = [0.47551926735449779, 0.41832860139077888, 0.36055278491307674, 0.3185536948526082, 0.28707643164399121, 0.27340350342025732, 0.26127222109599424, 0.25226028544973556, 0.24515773771353, 0.23866717394935788]
	embedded_coarse = [0.4870749999999997, 0.38741874999999959, 0.32206249999999981, 0.28391249999999979, 0.26998749999999971]#, 0.25991249999999977, 0.2521749999999997, 0.24614374999999988, 0.24108749999999968, 0.23688124999999982]
	fadlun_medium = [0.24250000000000016, 0.17491979568234528, 0.14586060138781931, 0.13039999999999996, 0.12059999999999998, 0.1152, 0.11109999999999998, 0.10769999999999999, 0.10500000000000004, 0.10269999999999999]
	external_medium = [0.23416141094834408, 0.17750264070933075, 0.15084344641480463, 0.13470391287586858, 0.12367567463377138, 0.11562735158057164, 0.10931040863531327, 0.10443342328450372, 0.10109999999999997, 0.098699999999999954]
	embedded_medium = [0.24565624999999958, 0.18659999999999954, 0.16298124999999947, 0.14938749999999967, 0.14061874999999974]#, 0.13428749999999967, 0.12947499999999967, 0.12559374999999967, 0.12237499999999973, 0.11968749999999978]
	fadlun_fine = [0.091099999999999959, 0.065099999999999991, 0.0567, 0.051600000000000007, 0.048299999999999982, 0.04579999999999998, 0.043999999999999984, 0.04250000000000001, 0.041300000000000003, 0.040299999999999989]
	external_fine = [0.094200000000000172, 0.070400000000000018, 0.059899999999999953, 0.053599999999999981, 0.049399999999999944, 0.046399999999999941, 0.043999999999999984, 0.042099999999999971, 0.04049999999999998, 0.039200000000000013]
	embedded_fine = [0.12070000000000003, 0.091299999999999992, 0.079499999999999904, 0.072399999999999964, 0.067899999999999905]#, 0.06469999999999998, 0.062199999999999978, 0.060199999999999976, 0.058399999999999924, 0.057100000000000012]
	plt.plot(tf,fadlun_coarse,'o-k',label='Modified Fadlun coarse')
	plt.plot(tf,external_coarse,'s-k',label='External coarse')
	plt.plot(tc,embedded_coarse,'^-k',label='Embedded coarse')

	plt.plot(tf,fadlun_medium,'o--r',label='Modified Fadlun medium')
	plt.plot(tf,external_medium,'s--r',label='External medium')
	plt.plot(tc,embedded_medium,'^--r',label='Embedded Fadlun medium')

	plt.plot(tf,fadlun_fine,'o:b',label='Modified Fadlun fine')
	plt.plot(tf,external_fine,'s:b',label='External fine')
	plt.plot(tc,embedded_fine,'^:b',label='Embedded fine')

	plt.xlabel('Time')
	plt.ylabel('Linf error norm')
	plt.title('Linf error norm')
	plt.legend([circle_,square_,triangle_,black_solid,red_dash,blue_dot],['Modified Fadlun','External','Embedded','Coarse grid','Medium grid','Fine grid'], bbox_to_anchor = (1.05,1.), loc = 2, borderaxespad=0.)
	plt.savefig('/scratch/src/cuIBM/validation/error/cylinder/Linf_error_norm.pdf', bbox_inches='tight')
	plt.clf()

if __name__ == "__main__":
	#main()
	fancyPlot()

