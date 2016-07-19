/***************************************************************************//**
 * \file intermediatePressure.h
 * \author Chris Minar
 * \brief Declaration of kernels to generate the right hand side of step 2: solve for intermediate pressure
 */

#pragma once

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void intermediatePressure_luo(double *rhs2, double *detA, int *hybridTagsP, double *alpha, double *stencilCoef,
								double *xv, double *yu,
								double *b11, double *b12, double *b13, double *b14, double *b21, double *b22, double *b23, double *b24,
								double *b31, double *b32, double *b33, double *b34, double *b41, double *b42, double *b43, double *b44,
								double *q1, double *q2, double *q3, double *q4,
								bool *q1flag, bool *q2flag, bool *q3flag, bool *q4flag,
								int nx, int ny);

__global__
void interpolate_P_HN_setup(double *detA, int *hybridTagsP, double *bx, double *by,
									double *uB, double *uB0, double *vB, double  *vB0,
									double *yu, double *xv,
									double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y,
									int *i_start, int *j_start, int width, int nx, int ny, double dt, double totalPoints,
									double *b11, double *b12, double *b13, double *b14, double *b21, double *b22, double *b23, double *b24,
									double *b31, double *b32, double *b33, double *b34, double *b41, double *b42, double *b43, double *b44,
									double *q1, double *q2, double *q3, double *q4,
									bool *q1flag, bool *q2flag, bool *q3flag, bool *q4flag,
									int *index1, int *index2, int *index3, int *index4,
									double *x1, double *x2, double *x3, double *x4,
									double *y1, double *y2, double *y3, double *y4,
									double *dudt, double *ududx, double *vdudy, double *dvdt, double *udvdx, double *vdvdy);//test

__global__
void hybridPressureNodeCount(int *countD, int *index1, int *index2, int *index3, int *index4, int *hybridTagsP,
								int *i_start, int *j_start, int width, int height, int nx, int ny);
}//end namespace kernels
