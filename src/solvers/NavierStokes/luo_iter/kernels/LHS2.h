#pragma once

namespace kernels
{
__global__
void LHS2_mid_iter(int *row, int *col, double *val, double *dx, double* dy, int nx, int ny, double dt);
}
