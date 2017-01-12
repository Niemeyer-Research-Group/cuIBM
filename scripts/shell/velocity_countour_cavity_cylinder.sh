#!/bin/sh
#runs lid driven cavity and impulsivly started cylinder then makes velocity contour plots of them
CUIBM_DIR=/scratch/src/cuIBM
#Plot Lid Driven Cavity
for var2 in 100 1000 10000
do
	echo "Plotting lid driven cavity flow for Reynolds number "$var2
	$CUIBM_DIR/scripts/python/plotVelocity.py -ulim 1 -folder $CUIBM_DIR/validation/lidDrivenCavity/Re$var2 >/dev/null
done

#Plot Cylinder
for var4 in 40 550 3000
do
	echo "Plotting flow over an impulsively started cylinder with Reynolds number "$var4
	$CUIBM_DIR/scripts/python/plotVelocity.py -folder $CUIBM_DIR/validation/cylinder/Re$var4 >/dev/null
done
