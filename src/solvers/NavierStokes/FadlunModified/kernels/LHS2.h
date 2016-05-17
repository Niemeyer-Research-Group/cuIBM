#pragma once

namespace kernels
{
__global__
void LHS2_mid(int *row, int *col, double *val, double *distance_from_u_to_body, double *distance_from_v_to_body, int *tagsP, int *tagsPout, double *dx, double* dy, int nx, int ny, double dt);
}
