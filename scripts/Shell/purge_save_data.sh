#!/bin/sh
#deletes all the saved data for all cases
CUIBM_DIR=/scratch/src/cuIBM-FSI

CAVITY_DIR=$CUIBM_DIR/validation/lidDrivenCavity
CYLINDER_DIR=$CUIBM_DIR/validation/cylinder

for filename in 100 1000 10000
do
	echo "Cleaning lid driven cavity save data for Reynolds number " $filename
	cd $CAVITY_DIR/Re$filename
	find . ! -name 'output' -type d -exec rm -r {} +
	cd ..
done

for filename in 40 75 100 150 200 550 3000
do
	echo "Cleaning cylinder save data for Reynolds number " $filename
	cd $CYLINDER_DIR/Re$filename
	find . ! -name 'output' -type d -exec rm -r {} +
	cd ..
done
