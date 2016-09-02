/***************************************************************************//**
 * \file
 * \author Anush Krishnan (anush@bu.edu), Christopher Minar (minarc@oregonstate.edu)
 * \brief Implementation of the methods of the class \c DirectForcingSolver to tag
 *        points near the immersed boundary using a ray-tracing algorithm.
 */

/**
 * \brief Tags the forcing nodes among the velocity nodes, i.e. the nodes at
 *        which the velocity interpolation is performed.
 */
#include <solvers/NavierStokes/luo_base.h>

#include <solvers/NavierStokes/luo_base/kernels/tagPoints.h>
#include <cusp/print.h>

void luo_base::tagPoints()
{
	logger.startTimer("tagPoints");
	double	i_start = B.startI[0],
			j_start = B.startJ[0],
			width_i = B.numCellsX[0],
			height_j = B.numCellsY[0],
			i_end = i_start + width_i,
			j_end = j_start + height_j;

	cusp::blas::fill(ghostTagsUV, -1);
	cusp::blas::fill(ghostTagsP, -1);
	cusp::blas::fill(hybridTagsUV, -1);
	cusp::blas::fill(hybridTagsP, -1);
	cusp::blas::fill(hybridTagsUV2, -1);
	cusp::blas::fill(distance_from_intersection_to_node, 1);
	cusp::blas::fill(distance_between_nodes_at_IB, 1);
	cusp::blas::fill(x1_ip,0);
	cusp::blas::fill(y1_ip,0);
	cusp::blas::fill(x2_ip,0);
	cusp::blas::fill(y2_ip,0);
	cusp::blas::fill(body_intercept_x, 0);
	cusp::blas::fill(body_intercept_y, 0);
	cusp::blas::fill(image_point_x, 0);
	cusp::blas::fill(image_point_y, 0);

	const int blocksize = 256;

	dim3 dimGrid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 dimBlock(blocksize, 1);
	dim3 dimGrid0(int( (i_end-i_start-0.5)/blocksize ) +1, 1);
	//tag u direction nodes for tags, tagsout and hybridTagsUV2
	kernels::tag_u_luo<<<dimGrid,dimBlock>>>(hybridTagsUV_r, ghostTagsUV_r, hybridTagsUV2_r, B.x_r, B.y_r, B.uB_r, B.vB_r, yu_r, xu_r,
											body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
											x1_r, y1_r, x2_r, y2_r,
											distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r, distance_from_u_to_body_r, distance_from_v_to_body_r, uv_r,
											i_start, j_start, i_end, j_end, nx, ny, B.totalPoints, B.midX, B.midY);
	//tag v direction nodes for tags, tagsout and tag2
	kernels::tag_v_luo<<<dimGrid,dimBlock>>>(hybridTagsUV_r, ghostTagsUV_r, hybridTagsUV2_r, B.x_r, B.y_r, B.uB_r, B.vB_r, yv_r, xv_r,
											body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
											x1_r, y1_r, x2_r, y2_r,
											distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r,distance_from_u_to_body_r, distance_from_v_to_body_r, uv_r,
											i_start, j_start, i_end, j_end, nx, ny, B.totalPoints, B.midX, B.midY);
	//tag pressure nodes for ghostTagsP and hybridTagsP
	kernels::tag_p_luo<<<dimGrid,dimBlock>>>(ghostTagsP_r, hybridTagsP_r,
											B.x_r, B.y_r, yu_r, xv_r,
											body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r, x1_p_r, y1_p_r, x2_p_r, y2_p_r,
											i_start, j_start, i_end, j_end, nx, ny, B.totalPoints, B.midX, B.midY);
	//zero the inside of ghostTagsP
	kernels::zero_pressure_luo<<<dimGrid0, dimBlock>>>(ghostTagsP_r, i_start, j_start, i_end, j_end, nx, ny);
	//zero the inside of ghostTagsUVx
	kernels::zero_x_luo<<<dimGrid0,dimBlock>>>(ghostTagsUV_r, i_start, j_start, i_end, j_end, nx, ny);
	//zero the inside of ghostTagsUVy
	kernels::zero_y_luo<<<dimGrid0,dimBlock>>>(ghostTagsUV_r, i_start, j_start, i_end, j_end, nx, ny);

	//testOutputX();
	//testOutputY();
	logger.stopTimer("tagPoints");
}
