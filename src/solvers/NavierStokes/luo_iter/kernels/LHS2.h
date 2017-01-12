#pragma once

namespace kernels
{
__global__
void LHS2_mid_iter(int *row, int *col, double *val, double *dx, double *dy, int nx, int ny, double dt,
					int *count, double *ns_rhs, double *interp_rhs, int *hybridTagsP, int *ghostTagsP,
					double *alpha, double *dpdn,
					int *index1, int *index2, int *index3, int *index4,
					double *q1coef, double *q2coef, double *q3coef, double *q4coef,
					double *q1, double *q2, double *q3, double *q4, int timeStep);
}
