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
#include <solvers/NavierStokes/kernels/tagPoints.h>

void NavierStokesSolver::tagPoints()
{
	logger.startTimer("tagPoints");
	int  nx = NavierStokesSolver::domInfo->nx,
		 ny = NavierStokesSolver::domInfo->ny,
		 totalPoints = B.totalPoints,
		 start_index_i = B.startI[0],
		 start_index_j = B.startJ[0],
		 width_i = B.numCellsX[0],
		 height_j = B.numCellsY[0],
		 end_index_i = start_index_i + width_i,
		 end_index_j = start_index_j + height_j;
	
	double	*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),//not sure if these are on the host or not
			*by_r		= thrust::raw_pointer_cast ( &(B.y[0]) ),
			*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) ),
			*yu_r		= thrust::raw_pointer_cast ( &(NavierStokesSolver::domInfo->yuD[0]) ),
			*xu_r		= thrust::raw_pointer_cast ( &(NavierStokesSolver::domInfo->xuD[0]) ),
			*yv_r		= thrust::raw_pointer_cast ( &(NavierStokesSolver::domInfo->yvD[0]) ),
			*xv_r		= thrust::raw_pointer_cast ( &(NavierStokesSolver::domInfo->xvD[0]) ),
			*a_r		= thrust::raw_pointer_cast ( &(aD[0]) ),
			*b_r		= thrust::raw_pointer_cast ( &(bD[0]) ),
			*dub_r		= thrust::raw_pointer_cast ( &(distance_from_u_to_bodyD[0]) ),
			*dvb_r		= thrust::raw_pointer_cast ( &(distance_from_v_to_bodyD[0]) ),
			*uv_r		= thrust::raw_pointer_cast ( &(uvD[0]) );
	
	int 	*tags_r		= thrust::raw_pointer_cast ( &(tagsD[0]) ),
			*tags2_r	= thrust::raw_pointer_cast ( &(tags2D[0]) ),
			*tagsIn_r	= thrust::raw_pointer_cast ( &(tagsInD[0]) ),
			*tagsP_r	= thrust::raw_pointer_cast ( &(tagsPD[0]) ),
			*tagsPOut_r	= thrust::raw_pointer_cast ( &(tagsPOutD[0]) );
	
	cusp::blas::fill(tagsD, -1);
	cusp::blas::fill(tags2D, -1);
	cusp::blas::fill(tagsInD, -1);
	cusp::blas::fill(tagsPD, -1);
	cusp::blas::fill(tagsPOutD, -1);
	cusp::blas::fill(aD, 1);
	cusp::blas::fill(bD, 1);
	
	double midX, midY;
	for (int i=0;i<totalPoints;i++)
	{
		midX += B.x[i];
		midY += B.y[i];
	}
	midX /= totalPoints;
	midY /= totalPoints;
	
	const int blocksize = 256;
	dim3 dimGrid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 dimBlock(blocksize, 1);
	dim3 dimGrid0(int( (end_index_i-start_index_i-0.5)/blocksize ) +1, 1);

	//tag u direction nodes for tags, tagsout and tags2
	kernels::tag_u<<<dimGrid,dimBlock>>>(tags_r, tagsIn_r, tags2_r,
										   bx_r, by_r, uB_r, vB_r, yu_r, xu_r, a_r, b_r, dub_r, dvb_r, uv_r, 
										   start_index_i, start_index_j, end_index_i, end_index_j, nx, ny, totalPoints, midX, midY);
	//tag v direction nodes for tags, tagsout and tag2
	kernels::tag_v<<<dimGrid,dimBlock>>>(tags_r, tagsIn_r, tags2_r,
										   bx_r, by_r, uB_r, vB_r, yv_r, xv_r, a_r, b_r, dub_r, dvb_r, uv_r, 
										   start_index_i, start_index_j, end_index_i, end_index_j, nx, ny, totalPoints, midX, midY);
	//tag pressure nodes for tagsp and tagspout
	kernels::tag_p<<<dimGrid,dimBlock>>>(tagsP_r, tagsPOut_r,
											bx_r, by_r, yu_r, xv_r, 
											start_index_i, start_index_j, end_index_i, end_index_j, nx, ny, totalPoints, midX, midY);
	//zero the inside of tagsp
	kernels::zero_pressure<<<dimGrid0, dimBlock>>>(tagsP_r, start_index_i, start_index_j, end_index_i, end_index_j, nx, ny);
	//zero the inside of tagsinx
	kernels::zero_x<<<dimGrid0,dimBlock>>>(tagsIn_r, start_index_i, start_index_j, end_index_i, end_index_j, nx, ny);
	//zero the inside of tagsiny
	kernels::zero_y<<<dimGrid0,dimBlock>>>(tagsIn_r, start_index_i, start_index_j, end_index_i, end_index_j, nx, ny);
	logger.stopTimer("tagPoints");
}
