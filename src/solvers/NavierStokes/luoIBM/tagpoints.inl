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
#include <solvers/NavierStokes/luoIBM/kernels/tagPoints.h>
#include <cusp/print.h>

void luoIBM::tagPoints()
{
	logger.startTimer("tagPoints");
	int  nx = NavierStokesSolver::domInfo->nx,
		 ny = NavierStokesSolver::domInfo->ny,
		 totalPoints = B.totalPoints,
		 i_start = B.startI[0],
		 j_start = B.startJ[0],
		 width_i = B.numCellsX[0],
		 height_j = B.numCellsY[0],
		 i_end = i_start + width_i,
		 j_end = j_start + height_j;
	
	double	*pressure_r = thrust::raw_pointer_cast ( &(pressure[0]) ),
			*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),//not sure if these are on the host or not
			*by_r		= thrust::raw_pointer_cast ( &(B.y[0]) ),
			*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) ),
			*yu_r		= thrust::raw_pointer_cast ( &(NavierStokesSolver::domInfo->yu[0]) ),
			*xu_r		= thrust::raw_pointer_cast ( &(NavierStokesSolver::domInfo->xu[0]) ),
			*yv_r		= thrust::raw_pointer_cast ( &(NavierStokesSolver::domInfo->yv[0]) ),
			*xv_r		= thrust::raw_pointer_cast ( &(NavierStokesSolver::domInfo->xv[0]) ),
			*a_r		= thrust::raw_pointer_cast ( &(distance_from_intersection_to_node[0]) ),
			*b_r		= thrust::raw_pointer_cast ( &(distance_between_nodes_at_IB[0]) ),
			*dub_r		= thrust::raw_pointer_cast ( &(distance_from_u_to_body[0]) ),
			*dvb_r		= thrust::raw_pointer_cast ( &(distance_from_v_to_body[0]) ),
			*uv_r		= thrust::raw_pointer_cast ( &(uv[0]) ),
			*body_intercept_x_r = thrust::raw_pointer_cast( &(body_intercept_x[0]) ),
			*body_intercept_y_r = thrust::raw_pointer_cast( &(body_intercept_y[0]) ),
			*image_point_x_r = thrust::raw_pointer_cast( &(image_point_x[0]) ),
			*image_point_y_r = thrust::raw_pointer_cast( &(image_point_y[0]) ),
			*x1_r = thrust::raw_pointer_cast( &(x1_ip[0]) ),
			*x2_r = thrust::raw_pointer_cast( &(x2_ip[0]) ),
			*y1_r = thrust::raw_pointer_cast( &(y1_ip[0]) ),
			*y2_r = thrust::raw_pointer_cast( &(y2_ip[0]) ),
			*body_intercept_p_x_r = thrust::raw_pointer_cast( &(body_intercept_p_x[0]) ),
			*body_intercept_p_y_r = thrust::raw_pointer_cast( &(body_intercept_p_y[0]) ),
			*image_point_p_x_r = thrust::raw_pointer_cast( &(image_point_p_x[0]) ),
			*image_point_p_y_r = thrust::raw_pointer_cast( &(image_point_p_y[0]) ),
			*x1_p_r = thrust::raw_pointer_cast( &(x1_ip_p[0]) ),
			*x2_p_r = thrust::raw_pointer_cast( &(x2_ip_p[0]) ),
			*y1_p_r = thrust::raw_pointer_cast( &(y1_ip_p[0]) ),
			*y2_p_r = thrust::raw_pointer_cast( &(y2_ip_p[0]) );
	
	int 	*ghostTagsUV_r		= thrust::raw_pointer_cast ( &(ghostTagsUV[0]) ),
			*hybridTagsUV2_r	= thrust::raw_pointer_cast ( &(hybridTagsUV2[0]) ),
			*hybridTagsUV_r		= thrust::raw_pointer_cast ( &(hybridTagsUV[0]) ),
			*ghostTagsP_r		= thrust::raw_pointer_cast ( &(ghostTagsP[0]) ),
			*hybridTagsP_r		= thrust::raw_pointer_cast ( &(hybridTagsP[0]) );
	
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
	kernels::tag_u_luo<<<dimGrid,dimBlock>>>(hybridTagsUV_r, ghostTagsUV_r, hybridTagsUV2_r, bx_r, by_r, uB_r, vB_r, yu_r, xu_r,
											body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
											x1_r, y1_r, x2_r, y2_r,
											a_r, b_r, dub_r, dvb_r, uv_r,
											i_start, j_start, i_end, j_end, nx, ny, totalPoints, B.midX, B.midY);
	//tag v direction nodes for tags, tagsout and tag2
	kernels::tag_v_luo<<<dimGrid,dimBlock>>>(hybridTagsUV_r, ghostTagsUV_r, hybridTagsUV2_r, bx_r, by_r, uB_r, vB_r, yv_r, xv_r,
											body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
											x1_r, y1_r, x2_r, y2_r,
											a_r, b_r, dub_r, dvb_r, uv_r, 
											i_start, j_start, i_end, j_end, nx, ny, totalPoints, B.midX, B.midY);
	//tag pressure nodes for ghostTagsP and hybridTagsP
	kernels::tag_p_luo<<<dimGrid,dimBlock>>>(ghostTagsP_r, hybridTagsP_r,
											bx_r, by_r, yu_r, xv_r,
											body_intercept_p_x_r, body_intercept_p_y_r, image_point_p_x_r, image_point_p_y_r, x1_p_r, y1_p_r, x2_p_r, y2_p_r,
											i_start, j_start, i_end, j_end, nx, ny, totalPoints, B.midX, B.midY);
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