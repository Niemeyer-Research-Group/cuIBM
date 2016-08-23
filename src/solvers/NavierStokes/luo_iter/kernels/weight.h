#pragma once

namespace kernels
{

__global__
void alpha_u(double *alpha, int *ghostTagsUV, int *hybridTagsUV, double *yu, double *xu,
				double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
				int *i_start, int *j_start, int width, int nx, int ny);

__global__
void alpha_v(double *alpha, int *ghostTagsUV, int *hybridTagsUV, double *yv, double *xv,
				double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
				int *i_start, int *j_start, int width, int nx, int ny);

__global__
void alpha_p(double *alpha, int *ghostTagsP, int *hybridTagsP, double *yu, double *xv,
				double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y,
				int *i_start, int *j_start, int width, int nx, int ny);
}
