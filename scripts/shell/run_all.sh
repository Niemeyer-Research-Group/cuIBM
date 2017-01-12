#!/bin/sh
#runs lid driven cavity and impulsivly started cylinder and makes validation plots for them
CUIBM_DIR=/scratch/src/cuIBM
#Run Lid Driven Cavity Validation
for var1 in 100 1000
do
	echo "Running lid driven cavity flow for Reynolds number " $var1
	$CUIBM_DIR/scripts/validation/cavity.py --Re $var1 >/dev/null
done

#Run Cylinder validation
for var3 in 40 550 3000
do
	echo "Running flow over an impulsively started cylinder with Reynolds number "$var3
	$CUIBM_DIR/scripts/validation/cylinder.py --Re $var3 >/dev/null
done
