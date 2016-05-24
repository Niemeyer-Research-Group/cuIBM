#pragma once

namespace kernels
{
__global__
void tag_u_luo(int *hybridTagsUV, int *ghostTagsUV, int *hybridTagsUV2, double *bx, double *by, double *uB, double *vB, double *yu, double *xu,
			double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y, double *x1, double *y1, double *x2, double *y2,
			double *a, double *b, double *distance_from_u_to_body, double *distance_from_v_to_body, double *uv,
			int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY);

__global__
void tag_v_luo(int *hybridTagsUV, int *ghostTagsUV, int *hybridTagsUV2, double *bx, double *by, double *uB, double *vB, double *yv, double *xv,
				double *body_intercept_x, double *body_intercept_y, double *image_point_x, double *image_point_y, double *x1, double *y1, double *x2, double *y2,
				double *distance_from_intersection_to_node, double *distance_between_nodes_at_IB, double *distance_from_u_to_body, double *distance_from_v_to_body, double *uv,
				int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY);

__global__
void tag_p_luo(int *ghostTagsP, int *hybridTagsP, double *bx, double *by, double *yu, double *xv,
		   int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY);

__global__
void zero_pressure_luo(int *ghostTagsP,  int i_start, int j_start, int i_end, int j_end, int nx, int ny);

__global__
void zero_x_luo(int *ghostTagsUV,  int i_start, int j_start, int i_end, int j_end, int nx, int ny);

__global__
void zero_y_luo(int *ghostTagsUV,  int i_start, int j_start, int i_end, int j_end, int nx, int ny);



}
