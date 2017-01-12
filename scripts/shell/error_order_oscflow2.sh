#!/bin/sh
#runs all cases
#plots all cases
CUIBM_DIR=/scratch/src/cuIBM
#Run Lid Driven Cavity Validation
cd /scratch/src/cuIBM/src
for var1 in 03125 05 0625 
do
	for var2 in external #embedded
	do
		echo "ERROR_ORDER_OSCFLOW.SH: Running "$var2$var1"\n"
		/scratch/src/cuIBM/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/error/oscflow2/$var2$var1
	done
done
#/scratch/src/cuIBM/scripts/validation/error_order_oscflow2.py
