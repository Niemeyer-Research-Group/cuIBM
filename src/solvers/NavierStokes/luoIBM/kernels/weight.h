#pragma once

namespace kernels
{

__global__
void weightX(double *uhat, double *ustar, int *ghostTagsUV, int *hybridTagsUV, double *yu, double *xu,
				double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
				int i_start, int j_start, int i_end, int j_end, int nx, int ny);

__global__
void weightY(double *uhat, double *ustar, int *ghostTagsUV, int *hybridTagsUV, double *yv, double *xv,
				double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
				int i_start, int j_start, int i_end, int j_end, int nx, int ny);

__global__
void weightP(double *pressure, double *pressureStar, int *ghostTagsP, int *hybridTagsP, double *yu, double *xv,
				double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
				int i_start, int j_start, int i_end, int j_end, int nx, int ny);
}
