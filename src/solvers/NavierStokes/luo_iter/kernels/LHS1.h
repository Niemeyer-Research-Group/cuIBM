#pragma once

namespace kernels
{
__global__
void LHS1_mid_iter_X(int *row, int *col, double *val, double *dx, double *dy, double dt, double nu, int nx, int ny,
						int *hybridTagsUV, int *ghostTagsUV, double *ns_rhs, double *interp_rhs, int *count,
						int *index1, int *index2, int *index3, int *index4,
						double *xu, double *yu, double *alpha, double *uB,
						double *q1coef, double *q2coef, double *q3coef, double *q4coef,
						double *q1, double *q2, double *q3, double *q4);

__global__
void LHS1_mid_iter_Y(int *row, int *col, double *val, double *dx, double *dy, double dt, double nu, int nx, int ny,
						int *hybridTagsUV, int *ghostTagsUV, double *ns_rhs, double *interp_rhs, int *count,
						int *index1, int *index2, int *index3, int *index4,
						double *xv, double *yv, double *alpha, double *vB,
						double *q1coef, double *q2coef, double *q3coef, double *q4coef,
						double *q1, double *q2, double *q3, double *q4);
}
