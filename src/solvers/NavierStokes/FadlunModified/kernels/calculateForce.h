#pragma once

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void calcForceFadlun(double *force, double *L, double *Lhat, double *Nold, double *N, double *u, double *uhat, int *tags, double dt, int nx, int ny);

__global__
void dragLeftRight(double *FxX, double *q, double *lambda, double nu, double *dx, double *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

__global__
void dragBottomTop(double *FxY, double *q, double nu, double *dx, double *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

__global__
void dragUnsteady(double *FxU, double *q, double *qOld, int *tagsIn, double *dx, double *dy, double dt,
                   int nx, int ny, int I, int J, int ncx, int ncy);

__global__
void liftLeftRight(double *FyX, double *q, double nu, double *dx, double *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

__global__
void liftBottomTop(double *FyY, double *q, double *lambda, double nu, double *dx, double *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy);

__global__
void liftUnsteady(double *FyU, double *q, double *qOld, int *tagsIn, double *dx, double *dy, double dt,
                  int nx, int ny, int I, int J, int ncx, int ncy);


} // end of namespace kernels

