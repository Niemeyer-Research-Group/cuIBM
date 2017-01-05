#!/bin/sh
#Converts u velocity png files made with scripts/plotVelocity.py into an mp4

CUIBM_DIR=/scratch/src/cuIBM
VIV_DIR=$CUIBM_DIR/validation/osc/VIV/Ured4/
pretty_dIR=$CUIBM_DIR/validation/osc/VIV/pretty/
echo "Creating plot gif"
cd $pretty_DIR
convert -delay 10 -loop 0 *.png plotgif.gif

echo "Creating velocity gif"
cd $VIV_DIR/vel_plot
convert -delay 10 -loop 0 *.png velgif.gif

echo "Creating vorticity gif"
cd $VIV_DIR/vort_plot
convert -delay 10 -loop 0 *.png vortgif.gif

echo "Creating pressure gif"
cd $VIV_DIR/p_plot
convert -delay 10 -loop 0 *.png pgif.gif

