/***************************************************************************//**
 * \file tagPoints.inl
 * \author Anush Krishnan (anush@bu.edu), Christopher Minar (minarc@oregonstate.edu)
 * \brief Implementation of the methods of the class \c DirectForcingSolver to tag
 *        points near the immersed boundary using a ray-tracing algorithm.
 */

/**
 * \brief Tags the forcing nodes among the velocity nodes, i.e. the nodes at
 *        which the velocity interpolation is performed.
 */
#include <solvers/NavierStokes/fadlunModified.h>

#include <solvers/NavierStokes/FadlunModified/kernels/tagPoints.h>
#include <cusp/print.h>

void fadlunModified::tagPoints()
{
	logger.startTimer("tagPoints");
	int  totalPoints = B.totalPoints,
		 i_start = B.startI[0],
		 j_start = B.startJ[0],
		 width_i = B.numCellsX[0],
		 height_j = B.numCellsY[0],
		 i_end = i_start + width_i,
		 j_end = j_start + height_j;

	//testing
	tagsOld = tags;
	tagsPOld = tagsPOut;
	cusp::blas::fill(tags, -1);
	cusp::blas::fill(tags2, -1);
	cusp::blas::fill(tagsIn, -1);
	cusp::blas::fill(tagsP, -1);
	cusp::blas::fill(tagsPOut, -1);
	cusp::blas::fill(distance_from_intersection_to_node, 1);
	cusp::blas::fill(distance_between_nodes_at_IB, 1);

	const int blocksize = 256;

	dim3 dimGrid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 dimBlock(blocksize, 1);
	dim3 dimGrid0(int( (i_end-i_start-0.5)/blocksize ) +1, 1);

	//tag u direction nodes for tags, tagsout and tags2
	kernels::tag_u<<<dimGrid,dimBlock>>>(tags_r, tagsIn_r, tags2_r,
										   B.x_r, B.y_r, B.uB_r, B.vB_r, yu_r, xu_r, distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r, distance_from_u_to_body_r, distance_from_v_to_body_r, uv_r,
										   i_start, j_start, i_end, j_end, nx, ny, totalPoints, B.midX, B.midY);
	//tag v direction nodes for tags, tagsout and tag2
	kernels::tag_v<<<dimGrid,dimBlock>>>(tags_r, tagsIn_r, tags2_r,
										   B.x_r, B.y_r, B.uB_r, B.vB_r, yv_r, xv_r, distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r, distance_from_u_to_body_r, distance_from_v_to_body_r, uv_r,
										   i_start, j_start, i_end, j_end, nx, ny, totalPoints, B.midX, B.midY);
	//tag pressure nodes for tagsp and tagspout
	kernels::tag_p<<<dimGrid,dimBlock>>>(tagsP_r, tagsPOut_r,
											B.x_r, B.y_r, yu_r, xv_r,
											i_start, j_start, i_end, j_end, nx, ny, totalPoints, B.midX, B.midY);
	//zero the inside of tagsp
	kernels::zero_pressure<<<dimGrid0, dimBlock>>>(tagsP_r, i_start, j_start, i_end, j_end, nx, ny);
	//zero the inside of tagsinx
	kernels::zero_x<<<dimGrid0,dimBlock>>>(tagsIn_r, i_start, j_start, i_end, j_end, nx, ny);
	//zero the inside of tagsiny
	kernels::zero_y<<<dimGrid0,dimBlock>>>(tagsIn_r, i_start, j_start, i_end, j_end, nx, ny);

	logger.stopTimer("tagPoints");
}
