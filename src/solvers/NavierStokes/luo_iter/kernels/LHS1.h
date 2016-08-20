#pragma once

namespace kernels
{
__global__
void LHS1_mid_iter_X(int *row, int *col, double *val, double *dx, double *dy, double dt, double nu, int nx, int ny,
		int *hybridTagsUV, int *ghostTagsUV, int *ns_rhs, int *interp_rhs, int *count,
		int *index1, int *index2, int *index3, int *index4,
		double *xu, double *yu, double *detA, double *alpha,
		double *b11, double *b12, double *b13, double *b14,
		double *b21, double *b22, double *b23, double *b24,
		double *b31, double *b32, double *b33, double *b34,
		double *b41, double *b42, double *b43, double *b44,
		double *q1, double *q2, double *q3, double *q4
		);

__global__
void LHS1_mid_iter_Y(int *row, int *col, double *val, int *ghostTagsUV, double *dx, double *dy, double dt, double nu, int nx, int ny);
}
