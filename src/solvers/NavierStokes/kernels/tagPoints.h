#pragma once

namespace kernels
{
__global__
void tag_u(int *tags, int *tagsIn, int *tags2, double *bx, double *by, double *uB, double *vB, double *yu, double *xu,
			   double *a, double *b, double *distance_from_u_to_body, double *distance_from_v_to_body, double *uv,
			   int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY);

__global__
void tag_v(int *tags, int *tagsIn, int *tags2, double *bx, double *by, double *uB, double *vB, double *yv, double *xv,
		   double *a, double *b, double *distance_from_u_to_body, double *distance_from_v_to_body, double *uv,
		   int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY);

__global__
void tag_p(int *tagsP, int *tagsPOut, double *bx, double *by, double *yu, double *xv,
		   int i_start, int j_start, int i_end, int j_end, int nx, int ny, int totalPoints, double midX, double midY);

__global__
void zero_pressure(int *tagsP,  int i_start, int j_start, int i_end, int j_end, int nx, int ny);

__global__
void zero_x(int *tagsIn,  int i_start, int j_start, int i_end, int j_end, int nx, int ny);

__global__
void zero_y(int *tagsIn,  int i_start, int j_start, int i_end, int j_end, int nx, int ny);



}
