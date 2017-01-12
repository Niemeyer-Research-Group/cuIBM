#!/bin/sh
#makes validation plots for lid driven cavity and impulsivly started cylinder (they must already have been run)
CUIBM_DIR=/scratch/src/cuIBM
#Run Lid Driven Cavity Validation
for var1 in 100 1000 10000
do
	echo "Plotting lid driven cavity validation for Reynolds number " $var1
	$CUIBM_DIR/scripts/validation/validateCavity.py --Re $var1 >/dev/null
done

#Run Cylinder validation
for var3 in 40 550 3000
do
	echo "Plotting validationg for over an impulsively started cylinder with Reynolds number "$var3
	$CUIBM_DIR/scripts/validation/validateCylinder.py --Re $var3 >/dev/null
done

