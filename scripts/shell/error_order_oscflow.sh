#!/bin/sh
#runs all cases
#plots all cases
CUIBM_DIR=/scratch/src/cuIBM
#Run Lid Driven Cavity Validation
cd /scratch/src/cuIBM/src
for var1 in 01 015625 03125 0625 
do
	for var2 in external
	do
	echo "ERROR_ORDER_OSCFLOW.SH: Running "$var2$var1"\n"
	#/scratch/src/cuIBM/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/error/oscflow/$var2$var1
	done
done
for var1 in 03125 #015625 02 03125 0625 
do
	for var2 in external
	do
	echo "ERROR_ORDER_OSCFLOW.SH: Running "$var2$var1"\n"
	echo /scratch/src/cuIBM/validation/error/oscflow/$var2$var1
	/scratch/src/cuIBM/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/error/oscflow/$var2$var1
	done
done
#/scratch/src/cuIBM/scripts/validation/error_order_oscflow.py
