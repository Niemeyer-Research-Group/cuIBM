#pragma once

namespace kernels
{
__global__
void LHS1_mid_luo_X(int *row, int *col, double *val, int *ghostTagsUV, double *dx, double *dy, double dt, double nu, int nx, int ny);

__global__
void LHS1_mid_luo_Y(int *row, int *col, double *val, int *ghostTagsUV, double *dx, double *dy, double dt, double nu, int nx, int ny);
}
