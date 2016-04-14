#!/bin/sh
CUIBM_DIR=/scratch/src/cuIBM-FSI

CAVITY_DIR=$CUIBM_DIR/validation/lidDrivenCavity
CYLINDER_DIR=$CUIBM_DIR/validation/cylinder

for filename in 100 1000 10000
do
	#echo "Creating video for lid driven cavity for Reynolds number " $filename
	cd $CAVITY_DIR/Re$filename
	cd ..
done

for filename in 75 #40 75 100 150 200 550
do
	echo "Creating video for cylinder for Reynolds number " $filename
	cd $CYLINDER_DIR/Re$filename
	/scratch/src/lib/ffmpeg/ffmpeg -f image2 -r 5 -i u%07d.png -vcodec mpeg4 -y CylinderRe$filename.mp4
	cd ..
done
