#pragma once

namespace kernels
{
__global__
void initialiseU(double *u, double *xu, double *yu, double uInitial, double uPerturb, double pi, double xmax, double xmin,double ymax,double ymin, int nx, int ny);

__global__
void initialiseV(double *u, double *xv, double *yv, double vInitial, double vPerturb, double pi, double xmax, double xmin,double ymax,double ymin, int nx, int ny);




}
