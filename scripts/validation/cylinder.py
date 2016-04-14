#!/usr/bin/env python

import argparse
import os
import os.path
import sys

cuibmFolder = os.path.expandvars("/scratch/src/cuIBM-FSI")

# Parse command line options
parser = argparse.ArgumentParser(description="Runs the validation case for impulsively started flow over a circular cylinder for a specified Reynolds number", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--Re", dest="Re", help="Reynolds number", default='40')
args = parser.parse_args()
Re = args.Re

if Re=='40':
	validationData = '/cylinderRe40-KL95.txt'
	yMax           = '6'
elif Re=='200' or Re =='100' or Re=='150' or Re == '75':
	validationData = ''
	yMax           = '6' 
elif Re=='550':
	validationData = '/cylinderRe550-KL95.txt'
	yMax           = '2' 
elif Re=='3000':
	validationData = '/cylinderRe3000-KL95.txt'
	yMax           = '2' 
else:
	print "Unavailable option for Reynolds number. Choose 40, 550 or 3000."
	sys.exit()

execPath       = cuibmFolder + '/bin/cuIBM'
caseFolder     = cuibmFolder + '/validation/cylinder/Re' + Re
validationData = cuibmFolder + '/validation-data' + validationData

print "\n"+"-"*120
print "Running the validation case for flow over an impulsively started circular cylinder at Reynolds number %s" % Re
print "-"*120+"\n"

runCommand = "%s -caseFolder %s" % (execPath, caseFolder)
print runCommand+"\n"
os.system(runCommand)

print "-"*120
print "Plotting the drag coefficient for flow over an impulsively started circular cylinder at Reynolds number %s" % Re
print "-"*120

gnuplotFile    = caseFolder + '/cylinderRe' + Re + '.plt'
outFile        = caseFolder + '/cylRe' + Re + 'Drag.png'

f = open(gnuplotFile, 'w')
f.write("reset;\nset terminal png enhanced font 'Palatino, 11';\n\n");
f.write("set title 'Flow over an impulsively started cylinder at Reynolds number %s'\n" % Re)
f.write("set xlabel 'Non-dimensional time'\n")
f.write("set ylabel 'Drag Coefficient'\n")
f.write("set output '%s'\n" % outFile)
f.write("plot [0:5] [-2:%s] \\\n" % yMax)
if validationData != cuibmFolder+'/validation-data':
	f.write("'%s/forces' u 1:(2*$2) w l lw 2 lc rgb '#3232ff' title 'Present Work', \\\n" % caseFolder)
	f.write("'%s' u (0.5*$1):2 w p pt 6 ps 2 lc rgb '#ff3232' title 'Koumoutsakos and Leonard, 1995'\n" % validationData)
else:
	f.write("'%s/forces' u 1:(2*$2) w l lw 2 lc rgb '#3232ff' title 'Present Work' \\\n" % caseFolder)
	
f.close()

print "\nCreated gnuplot script "+gnuplotFile
runCommand = "gnuplot "+gnuplotFile
os.system(runCommand)
print "\nDone plotting! Plot saved to file "+outFile+"\n"