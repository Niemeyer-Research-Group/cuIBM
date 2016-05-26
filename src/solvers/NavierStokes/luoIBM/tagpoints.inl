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
		 start_index_i = B.startI[0],
		 start_index_j = B.startJ[0],
		 width_i = B.numCellsX[0],
		 height_j = B.numCellsY[0],
		 end_index_i = start_index_i + width_i,
		 end_index_j = start_index_j + height_j;
	
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
			*y2_r = thrust::raw_pointer_cast( &(y2_ip[0]) );
	
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
	dim3 dimGrid0(int( (end_index_i-start_index_i-0.5)/blocksize ) +1, 1);
	
	//tag u direction nodes for tags, tagsout and hybridTagsUV2
	kernels::tag_u_luo<<<dimGrid,dimBlock>>>(hybridTagsUV_r, ghostTagsUV_r, hybridTagsUV2_r, bx_r, by_r, uB_r, vB_r, yu_r, xu_r,
											body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
											x1_r, y1_r, x2_r, y2_r,
											a_r, b_r, dub_r, dvb_r, uv_r,
											start_index_i, start_index_j, end_index_i, end_index_j, nx, ny, totalPoints, B.midX, B.midY);
	//tag v direction nodes for tags, tagsout and tag2
	kernels::tag_v_luo<<<dimGrid,dimBlock>>>(hybridTagsUV_r, ghostTagsUV_r, hybridTagsUV2_r, bx_r, by_r, uB_r, vB_r, yv_r, xv_r,
											body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
											x1_r, y1_r, x2_r, y2_r,
											a_r, b_r, dub_r, dvb_r, uv_r, 
											start_index_i, start_index_j, end_index_i, end_index_j, nx, ny, totalPoints, B.midX, B.midY);
	//tag pressure nodes for ghostTagsP and hybridTagsP
	kernels::tag_p_luo<<<dimGrid,dimBlock>>>(ghostTagsP_r, hybridTagsP_r,
											bx_r, by_r, yu_r, xv_r, 
											start_index_i, start_index_j, end_index_i, end_index_j, nx, ny, totalPoints, B.midX, B.midY);
	//zero the inside of ghostTagsP
	kernels::zero_pressure_luo<<<dimGrid0, dimBlock>>>(ghostTagsP_r, start_index_i, start_index_j, end_index_i, end_index_j, nx, ny);
	//zero the inside of ghostTagsUVx
	kernels::zero_x_luo<<<dimGrid0,dimBlock>>>(ghostTagsUV_r, start_index_i, start_index_j, end_index_i, end_index_j, nx, ny);
	//zero the inside of ghostTagsUVy
	kernels::zero_y_luo<<<dimGrid0,dimBlock>>>(ghostTagsUV_r, start_index_i, start_index_j, end_index_i, end_index_j, nx, ny);
	
	//testOutputX();
	testOutputY();
	logger.stopTimer("tagPoints");
}

void luoIBM::testOutputX()
{
	int iu;
	int nx = NavierStokesSolver::domInfo->nx,
		start_index_i = B.startI[0],
		start_index_j = B.startJ[0],
		width_i = B.numCellsX[0],
		height_j = B.numCellsY[0],
		end_index_i = start_index_i + width_i,
		end_index_j = start_index_j + height_j;
	std::ofstream body_nodes;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/body_nodesX.csv";
	body_nodes.open(out.str().c_str());
	body_nodes << "BN_X1\tBN_Y1\tBN_X2\tBN_Y2\tGN_X\tGN_Y\tBI_X\tBI_Y\tIP_X\tIP_Y\n";
	for (int J=start_index_j;  J<end_index_j;  J++)
	{
		for (int I=start_index_i;  I<end_index_i;  I++)
		{
			iu = J*(nx-1) + I;
			if (hybridTagsUV[iu] >0) //for testing outside interpolation
			//if (ghostTagsUV[iu] >0) //for testing inside interpolation
			{
				std::cout<<iu<<std::endl;
				body_nodes << x1_ip[iu]<<"\t";
				body_nodes << y1_ip[iu]<<"\t";
				body_nodes << x2_ip[iu]<<"\t";
				body_nodes << y2_ip[iu]<<"\t";
				body_nodes << domInfo->xu[I]<<"\t";
				body_nodes << domInfo->yu[J]<<"\t";
				body_nodes << body_intercept_x[iu] <<"\t";
				body_nodes << body_intercept_y[iu] <<"\t";
				body_nodes << image_point_x[iu] <<"\t";
				body_nodes << image_point_y[iu] <<"\n";
			}
		}
	}
	body_nodes.close();
}

void luoIBM::testOutputY()
{
	//test bi, ip
		int iv;
		int nx = NavierStokesSolver::domInfo->nx,
			ny = NavierStokesSolver::domInfo->ny,
			start_index_i = B.startI[0],
			 start_index_j = B.startJ[0],
			 width_i = B.numCellsX[0],
			 height_j = B.numCellsY[0],
			 end_index_i = start_index_i + width_i,
			 end_index_j = start_index_j + height_j;
		std::ofstream body_nodes;
		parameterDB  &db = *NavierStokesSolver::paramDB;
		std::string folder = db["inputs"]["caseFolder"].get<std::string>();
		std::stringstream out;
		out << folder << "/body_nodesY.csv";
		body_nodes.open(out.str().c_str());
		body_nodes << "x1\ty1\tx2\ty2\tg_x\tg_y\tbi_x\tbi_y\tip_x\tip_y\n";
		arrayprint(ghostTagsUV,"gnuvY","y");
		arrayprint(ghostTagsUV,"gnuvX","x");
		for (int J=start_index_j;  J<end_index_j;  J++)
		{
			for (int I=start_index_i;  I<end_index_i;  I++)
			{
				iv = J*(nx) + I + ny*(nx-1);
				if (hybridTagsUV[iv] >0) //for testing outside interpolation
				//if (ghostTagsUV[iv] >0) //for testing inside interpolation
				{
					std::cout<<iv<<std::endl;
					body_nodes << x1_ip[iv]<<"\t";
					body_nodes << y1_ip[iv]<<"\t";
					body_nodes << x2_ip[iv]<<"\t";
					body_nodes << y2_ip[iv]<<"\t";
					body_nodes << domInfo->xv[I]<<"\t";
					body_nodes << domInfo->yv[J]<<"\t";
					body_nodes << body_intercept_x[iv] <<"\t";
					body_nodes << body_intercept_y[iv] <<"\t";
					body_nodes << image_point_x[iv] <<"\t";
					body_nodes << image_point_y[iv] <<"\n";
				}
			}
		}
		
		body_nodes.close();
}