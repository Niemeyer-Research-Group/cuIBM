#pragma once

namespace kernels
{
__global__
void interpolateVelocityToGhostNodeX(double *u, int *ghostTagsUV, double *bx, double *by, double *uB, double *yu, double *xu,
									double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
									int *i_start, int *j_start, int width, int nx, int ny,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *image_point_u);//testing variables
__global__
void interpolateVelocityToGhostNodeY(double *u, int *ghostTagsUV, double *bx, double *by, double *vB, double *yv, double *xv,
									double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
									int *i_start, int *j_start, int width, int nx, int ny,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *image_point_u);//testing variables
__global__
void interpolateVelocityToHybridNodeX(double *u, double *ustar, int *hybridTagsUV, double *detA,
										double *bx, double *by, double *uB, double *yu, double *xu,
										double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
										int *i_start, int *j_start, int width, int nx, int ny,
										int *index1, int *index2, int *index3, int *index4,
										double *b11, double *b12, double *b13, double *b14,
										double *b21, double *b22, double *b23, double *b24,
										double *b31, double *b32, double *b33, double *b34,
										double *b41, double *b42, double *b43, double *b44,
										double *x1, double *x2, double *x3, double *x4,
										double *y1, double *y2, double *y3, double *y4,
										double *q1, double *q2, double *q3, double *q4, double *image_point_u);//test
__global__
void interpolateVelocityToHybridNodeY(double *u, double *ustar, int *hybridTagsUV, double *detA,
										double *bx, double *by, double *vB, double *yv, double *xv,
										double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
										int *i_start, int *j_start, int width, int nx, int ny,
										int *index1, int *index2, int *index3, int *index4,
										double *b11, double *b12, double *b13, double *b14,
										double *b21, double *b22, double *b23, double *b24,
										double *b31, double *b32, double *b33, double *b34,
										double *b41, double *b42, double *b43, double *b44,
										double *x1, double *x2, double *x3, double *x4,
										double *y1, double *y2, double *y3, double *y4,
										double *q1, double *q2, double *q3, double *q4, double *image_point_u);//test
__global__
void interpolatePressureToHybridNode(double *pressure, double *pressureStar, double *u, int *hybridTagsP, double *bx, double *by,
									double *uB, double *uB0, double *vB, double  *vB0, double *yu, double *yv, double *xu, double *xv,
									double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y,
									int *i_start, int *j_start, int width, int nx, int ny, double dt,
									double *a0, double *a1, double *a2, double *a3,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, int timeStep);//test
__global__
void interpolatePressureToGhostNode(double *pressure, double *u, int *ghostTagsP, double *bx, double *by,
									double *uB, double *uB0, double *vB, double  *vB0, double *yu, double *yv, double *xu, double *xv,
									double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y, double *body_intercept_p,
									int *i_start, int *j_start, int width, int nx, int ny, double dt,
									double *a0, double *a1, double *a2, double *a3,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4);//test
}
