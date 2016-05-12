#!/bin/sh
#Converts u velocity png files made with scripts/plotVelocity.py into an mp4
CUIBM_DIR=/scratch/src/cuIBM

CYLINDER_DIR=$CUIBM_DIR/validation/cylinder

#Plot Cylinder
for var4 in 40 75 100 150 200 550 3000
do
	echo "Plotting vorticity CylinderRe"$var4
	$CUIBM_DIR/scripts/python/plotVorticity.py -folder $CUIBM_DIR/validation/cylinder/Re$var4 >/dev/null -vortlim 10 -xmin -3 -xmax 15 -ymin -4 -ymax 4
done

for filename in 40 75 100 150 200 550 3000
do
	echo "Creating video Re" $filename
	cd $CYLINDER_DIR/Re$filename
	/scratch/src/lib/ffmpeg/ffmpeg -f image2 -r 5 -i o%07d.png -vcodec mpeg4 -y VorticityCylinderRe$filename.mp4
	cd ..
done

