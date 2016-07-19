#pragma once

namespace kernels
{
__global__
void LHS2_mid_luo(int *row, int *col, double *val, double *dx, double *dy, int nx, int ny, double dt, int *count, double *stencilCoef, double *interpCoef,
					double *detA, int *hybridTagsP, int *ghostTagsP, double *alpha,
					double *xv, double *yu,
					double *b11, double *b12, double *b13, double *b14, double *b21, double *b22, double *b23, double *b24,
					double *b31, double *b32, double *b33, double *b34, double *b41, double *b42, double *b43, double *b44,
					/*double *q1, double *q2, double *q3, double *q4,*/ //not used
					bool *q1flag, bool *q2flag, bool *q3flag, bool *q4flag, //not currently used
					/*double *x1, double *x2, double *x3, double *x4, //not used
					double *y1, double *y2, double *y3, double *y4,*/ //not used
					int *index1, int *index2, int *index3, int *index4);
}
