#pragma once

namespace kernels
{
__global__
void interpolateVelocityX(double *u, int *ghostTagsUV, double *bx, double *by, double *uB, double *vB, double *yu, double *xu,
							double *image_point_x, double *image_point_y,
							int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints);
__global__
void interpolateVelocityY(double *u, int *ghostTagsUV, double *bx, double *by, double *uB, double *vB, double *yv, double *xv,
							double *image_point_x, double *image_point_y,
							int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints);
}
