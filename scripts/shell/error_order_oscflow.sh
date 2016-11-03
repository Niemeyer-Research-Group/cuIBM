#!/bin/sh
#runs all cases
#plots all cases
CUIBM_DIR=/scratch/src/cuIBM
#Run Lid Driven Cavity Validation
cd /scratch/src/cuIBM/src
for var1 in 015625 #01 #0625 03125
do
	for var2 in embedded #external
	do
		echo "ERROR_ORDER_OSCFLOW.SH: Running "$var2$var1"\n"
		/scratch/src/cuIBM/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/error/oscflow/$var2$var1
	done
done
/scratch/src/cuIBM/scripts/validation/error_order_oscflow.py
