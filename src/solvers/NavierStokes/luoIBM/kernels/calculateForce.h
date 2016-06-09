#pragma once

namespace kernels
{
__global__
void pressure_at_BI(double *force_pressure, double *pressure, double *u, int *ghostTagsP, int *hybridTagsP, double *bx, double *by,
						double *uB, double *uB0, double *vB, double  *vB0, double *yu, double *xv,
						double *x1, double *x2, double *x3, double *x4, double *y1, double *y2, double *y3, double *y4, double *q1, double *q2, double *q3, double *q4,
						double *point_x, double *point_y, double *point2_x, double *point2_y, double *point3_x, double *point3_y,
						int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY);
}
