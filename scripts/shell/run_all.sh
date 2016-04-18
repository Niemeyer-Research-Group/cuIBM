#!/bin/sh
#runs all cases
#plots all cases
CUIBM_DIR=/scratch/src/cuIBM
#Run Lid Driven Cavity Validation
for var1 in 100 1000 10000
do
	echo "Running lid driven cavity flow for Reynolds number " $var1
	$CUIBM_DIR/scripts/validation/cavity.py --Re $var1 >/dev/null
done


#Plot Lid Driven Cavity
for var2 in 100 1000 10000
do
	echo "Plotting lid driven cavity flow for Reynolds number "$var2
	$CUIBM_DIR/scripts/python/plotVelocity.py -folder $CUIBM_DIR/validation/lidDrivenCavity/Re$var2 >/dev/null
done


#Run Cylinder validation
for var3 in 40 75 100 150 200 550 3000
do
	echo "Running flow over an impulsively started cylinder with Reynolds number "$var3
	$CUIBM_DIR/scripts/validation/cylinder.py --Re $var3 >/dev/null
done

#Plot Cylinder
for var4 in 40 75 100 150 200 550 3000
do
	echo "Plotting flow over an impulsively started cylinder with Reynolds number "$var4
	$CUIBM_DIR/scripts/python/plotVelocity.py -folder $CUIBM_DIR/validation/cylinder/Re$var4 >/dev/null
done
