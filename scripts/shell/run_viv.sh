#!/bin/sh
#runs all cases
#plots all cases
CUIBM_DIR=/scratch/src/cuIBM
#Run Lid Driven Cavity Validation
cd /scratch/src/cuIBM/src
for var1 in 4 5 6 7 8
do
	echo "Running VIV for Ured " $var1
	make vivUred$var1
done
/scratch/src/cuIBM/scripts/validation/validate_viv.py
