#pragma once

namespace kernels
{
__global__
void interpolateVelocityToGhostNodeX(double *u, int *ghostTagsUV, double *bx, double *by, double *uB, double *yu, double *xu,
							double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
							int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints,//flag vB needed?
							double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *ip_u);//testing variables
__global__
void interpolateVelocityToGhostNodeY(double *u, int *ghostTagsUV, double *bx, double *by, double *vB, double *yv, double *xv,
							double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
							int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints,
							double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4, double *ip_u);//testing variables

__global__
void interpolateVelocityToHybridNodeX(double *u, double *ustar, int *hybridTagsUV, double *bx, double *by, double *uB, double *yu, double *xu,
									double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
									int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints,
									double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4);//test
}
