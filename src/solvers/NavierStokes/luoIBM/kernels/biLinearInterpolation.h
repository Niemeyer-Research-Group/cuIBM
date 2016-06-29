#pragma once

namespace kernels
{
__global__
void interpolateVelocityToGhostNodeX(double *u, int *ghostTagsUV, double *bx, double *by, double *uB, double *yu, double *xu,
									double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
									int *i_start, int *j_start, int width, int nx, int ny,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *ip_u);//testing variables
__global__
void interpolateVelocityToGhostNodeY(double *u, int *ghostTagsUV, double *bx, double *by, double *vB, double *yv, double *xv,
									double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
									int *i_start, int *j_start, int width, int nx, int ny,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *ip_u);//testing variables
__global__
void interpolateVelocityToHybridNodeX(double *u, double *ustar, int *hybridTagsUV, double *bx, double *by, double *uB, double *yu, double *xu,
									double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
									int *i_start, int *j_start, int width, int nx, int ny,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4);//test
__global__
void interpolateVelocityToHybridNodeY(double *u, double *ustar, int *hybridTagsUV, double *bx, double *by, double *vB, double *yv, double *xv,
									double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
									int *i_start, int *j_start, int width, int nx, int ny,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4);//test
__global__
void interpolatePressureToHybridNode(double *pressure, double *pressureStar, double *u, int *hybridTagsP, double *bx, double *by,
									double *uB, double *uB0, double *vB, double  *vB0, double *yu, double *yv, double *xu, double *xv,
									double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y,
									int *i_start, int *j_start, int width, int nx, int ny,
									double *dudt, double *ududx, double *vdudy, double *dvdt, double *udvdx, double *vdvdy,
									double *a0, double *a1, double *a2, double *a3,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, int timeStep);//test
__global__
void interpolatePressureToGhostNode(double *pressure, double *u, int *ghostTagsP, double *bx, double *by,
									double *uB, double *uB0, double *vB, double  *vB0, double *yu, double *yv, double *xu, double *xv,
									double *body_intercept_p_x, double *body_intercept_p_y, double *image_point_p_x, double *image_point_p_y,
									int *i_start, int *j_start, int width, int nx, int ny,
									double *dudt, double *ududx, double *vdudy, double *dvdt, double *udvdx, double *vdvdy,
									double *a0, double *a1, double *a2, double *a3,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4);//test
}
