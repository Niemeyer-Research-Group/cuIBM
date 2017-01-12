#!/bin/sh
#runs all cases needed to make the order of accuracy oscflow plot then plots it
CUIBM_DIR=/scratch/src/cuIBM
#Run Lid Driven Cavity Validation
cd /scratch/src/cuIBM/src
for var1 in 015625 02 #03125 0625 
do
	for var2 in external embedded
	do
	echo "ERROR_ORDER_OSCFLOW.SH: Running "$var2$var1"\n"
	echo /scratch/src/cuIBM/validation/error/oscflow/$var2$var1
	/scratch/src/cuIBM/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/error/oscflow/$var2$var1
	done
done
/scratch/src/cuIBM/scripts/validation/error_order_oscflow.py
