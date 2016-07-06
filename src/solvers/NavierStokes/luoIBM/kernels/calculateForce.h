#pragma once

namespace kernels
{
__global__
void force_pressure(double *force_pressure, double *body_intercept_p,
					double *body_intercept_p_x, double *body_intercept_p_y,
					double *bx, double *by, double *xv, double *yu, int *ghostTagsP,
					int *i_start, int *j_start, int width, int height, int totalPoints, int nx, int ny, double midX, double midY);
__global__
void force_velocity_x(double *force_dudx, double *uB, double *u,
						double *bx, double *by, double *xu, double *yu,
						int *i_start, int *j_start, int width, int height, int totalPoints, int nx, int ny, double midX, double midY, double dx);
__global__
void force_velocity_y(double *force_dvdx, double *vB, double *u,
						double *bx, double *by, double *xv, double *yv,
						int *i_start, int *j_start, int width, int height, int totalPoints, int nx, int ny, double midX, double midY, double dx);

__global__
void force(double *force_x, double *force_y,  double *pressure, double *dudn, double *dvdn,
			double *bx, double *by,
			int totalPoints, double midX, double midY, double nu);
}
