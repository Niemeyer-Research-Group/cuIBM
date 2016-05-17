#pragma once

namespace kernels
{
__global__
void LHS_BC_X(int *row, int *col, double *val, double *dx, double *dy, double dt, double nu, int nx, int ny);

__global__
void LHS_BC_Y(int *row, int *col, double *val, double *dx, double *dy, double dt, double nu, int nx, int ny);

__global__
void LHS_mid_X_nobody(int *row, int *col, double  *val, double *dx, double *dy, double dt, double nu, int nx, int ny);

__global__
void LHS_mid_Y_nobody(int *row, int *col, double  *val, double *dx, double *dy, double dt, double nu, int nx, int ny);
}
