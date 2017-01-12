#!/bin/sh
#runs all cases needed to make the cylinder order of accuracy plot, then plots it
CUIBM_DIR=/scratch/src/cuIBM
#Run Lid Driven Cavity Validation
cd /scratch/src/cuIBM/src
for var1 in 1 2 3 4
do
	for var2 in fadlun external embedded 
	do
		echo "Running "$var2$var1" impulsively started cylinder for Re40\n"
		/scratch/src/cuIBM/bin/cuIBM -caseFolder /scratch/src/cuIBM/validation/error/cylinder/$var2$var1
	done
done
/scratch/src/cuIBM/scripts/validation/error_order_cylinder.py
#/scratch/src/cuIBM/scripts/validation/performance.py
