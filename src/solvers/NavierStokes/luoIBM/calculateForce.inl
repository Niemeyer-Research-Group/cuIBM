/***************************************************************************//**
 * \file calculateForce.inl
 * \author Anush Krishnan (anush@bu.edu),
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke kernels that will calculate the force on the immersed body
 */

#include <solvers/NavierStokes/FadlunModified/kernels/calculateForce.h>
#include <solvers/NavierStokes/luoIBM/kernels/calculateForce.h>

/**
 * \brief Calculates forces acting on an immersed body (on the device).
 *
 * Uses the control volume approach explained by Lai and Peskin (2000).
 * This is a general method that can be used with any immersed boundary method.
 * It uses only the velocity and pressure fields to calculate the forces, and
 * does not involve any body forces on the immersed boundary.
 * Currently works only for one body.
 */
void luoIBM::calculateForce()
{
	int  nx = NavierStokesSolver::domInfo->nx,
	     ny = NavierStokesSolver::domInfo->ny;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	double	dt = db["simulation"]["dt"].get<double>(),
			nu = db["flow"]["nu"].get<double>();

	double	*u_r		= thrust::raw_pointer_cast(&NavierStokesSolver::u[0]),
			*uold_r		= thrust::raw_pointer_cast(&NavierStokesSolver::uold[0]),
			*pressure_r	= thrust::raw_pointer_cast(&NavierStokesSolver::pressure[0]),
			*dx		= thrust::raw_pointer_cast(&(NavierStokesSolver::domInfo->dx[0])),
			*dy		= thrust::raw_pointer_cast(&(NavierStokesSolver::domInfo->dy[0]));

	int		*ghostTagsUV_r	= thrust::raw_pointer_cast( &(ghostTagsUV[0]) );

	// Calculating drag
	cusp::array1d<double, cusp::device_memory>
		FxX(B.numCellsY[0]),
		FxY(B.numCellsX[0]+1),
		FxU((B.numCellsX[0]+1)*B.numCellsY[0]);

	double *FxX_r = thrust::raw_pointer_cast(&FxX[0]),
	     *FxY_r = thrust::raw_pointer_cast(&FxY[0]),
	     *FxU_r = thrust::raw_pointer_cast(&FxU[0]);

	const int blockSize = 256;
	dim3 dimGrid( int((B.numCellsX[0]+B.numCellsY[0]+1-0.5)/blockSize)+1, 1 );
	dim3 dimBlock(blockSize, 1);

	kernels::dragLeftRight <<<dimGrid, dimBlock>>> (FxX_r, u_r, pressure_r, nu, dx, dy, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::dragBottomTop <<<dimGrid, dimBlock>>> (FxY_r, u_r, nu, dx, dy, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridX( int( ( (B.numCellsX[0]+1)*B.numCellsY[0]-0.5 )/blockSize )+1, 1 );
	kernels::dragUnsteady <<<dimGridX, dimBlock>>> (FxU_r, u_r, uold_r, ghostTagsUV_r, dx, dy, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	fxx = thrust::reduce(FxX.begin(), FxX.end());
	fxy = thrust::reduce(FxY.begin(), FxY.end());
	fxu = thrust::reduce(FxU.begin(), FxU.end());
	B.forceX =  fxx + fxy + fxu;

	// Calculating lift
	cusp::array1d<double, cusp::device_memory>
		FyX(B.numCellsY[0]+1),
		FyY(B.numCellsX[0]),
		FyU((B.numCellsX[0]+1)*B.numCellsY[0]);

	double	*FyX_r = thrust::raw_pointer_cast(&FyX[0]),
			*FyY_r = thrust::raw_pointer_cast(&FyY[0]),
			*FyU_r = thrust::raw_pointer_cast(&FyU[0]);

	kernels::liftLeftRight <<<dimGrid, dimBlock>>> (FyX_r, u_r, nu, dx, dy, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::liftBottomTop <<<dimGrid, dimBlock>>> (FyY_r, u_r, pressure_r, nu, dx, dy, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridY( int( ( B.numCellsX[0]*(B.numCellsY[0]+1)-0.5 )/blockSize )+1, 1 );
	kernels::liftUnsteady <<<dimGridY, dimBlock>>> (FyU_r, u_r, uold_r, ghostTagsUV_r, dx, dy, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	B.forceY = thrust::reduce(FyX.begin(), FyX.end()) + thrust::reduce(FyY.begin(), FyY.end()) + thrust::reduce(FyU.begin(), FyU.end());
}

void luoIBM::luoForce()
{
	double	*force_pressure_r = thrust::raw_pointer_cast ( &(B.force_pressure[0])),
			*pressure_r = thrust::raw_pointer_cast ( &(pressure[0]) ),
			*u_r 		= thrust::raw_pointer_cast ( &(u[0]) ),
			*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),
			*by_r		= thrust::raw_pointer_cast ( &(B.y[0]) ),
			*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) ),
			*uB0_r		= thrust::raw_pointer_cast ( &(B.uBk[0]) ),
			*vB0_r		= thrust::raw_pointer_cast ( &(B.vBk[0]) ),
			*yu_r		= thrust::raw_pointer_cast ( &(domInfo->yu[0]) ),
			*xv_r		= thrust::raw_pointer_cast ( &(domInfo->xv[0]) );
	
	//test
	double	*x1_r	= thrust::raw_pointer_cast ( &(B.x1[0]) ),
			*x2_r	= thrust::raw_pointer_cast ( &(B.x2[0]) ),
			*x3_r	= thrust::raw_pointer_cast ( &(B.x3[0]) ),
			*x4_r	= thrust::raw_pointer_cast ( &(B.x4[0]) ),
			*y1_r	= thrust::raw_pointer_cast ( &(B.y1[0]) ),
			*y2_r	= thrust::raw_pointer_cast ( &(B.y2[0]) ),
			*y3_r	= thrust::raw_pointer_cast ( &(B.y3[0]) ),
			*y4_r	= thrust::raw_pointer_cast ( &(B.y4[0]) ),
			*q1_r	= thrust::raw_pointer_cast ( &(B.q1[0]) ),
			*q2_r	= thrust::raw_pointer_cast ( &(B.q2[0]) ),
			*q3_r	= thrust::raw_pointer_cast ( &(B.q3[0]) ),
			*q4_r	= thrust::raw_pointer_cast ( &(B.q4[0]) ),
			*px_r	= thrust::raw_pointer_cast ( &(B.point_x[0]) ),
			*py_r	= thrust::raw_pointer_cast ( &(B.point_y[0]) ),
			*px2_r	= thrust::raw_pointer_cast ( &(B.point2_x[0]) ),
			*py2_r	= thrust::raw_pointer_cast ( &(B.point2_y[0]) ),
			*px3_r	= thrust::raw_pointer_cast ( &(B.point3_x[0]) ),
			*py3_r	= thrust::raw_pointer_cast ( &(B.point3_y[0]) );
	
	int 	*ghostTagsP_r		= thrust::raw_pointer_cast ( &(ghostTagsP[0]) ),
			*hybridTagsP_r		= thrust::raw_pointer_cast ( &(hybridTagsP[0]) );

	int nx = domInfo ->nx,
		ny = domInfo ->ny,
		totalPoints = B.totalPoints,
		i_start = B.startI[0],
		j_start = B.startJ[0],
		width_i = B.numCellsX[0],
		height_j = B.numCellsY[0],
		i_end = i_start + width_i,
		j_end = j_start + height_j;
	
	const int blocksize = 256;
	dim3 grid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	
	kernels::pressure_at_BI<<<grid,block>>> (force_pressure_r, pressure_r, u_r, ghostTagsP_r, hybridTagsP_r, bx_r, by_r,
												uB_r, uB0_r, vB_r, vB0_r, yu_r, xv_r,
												x1_r, x2_r, x3_r, x4_r, y1_r, y2_r, y3_r, y4_r, q1_r, q2_r, q3_r, q4_r,
												px_r, py_r, px2_r, py2_r, px3_r, py3_r,
												i_start, j_start, i_end, j_end, nx, ny, totalPoints, B.midX, B.midY);
	
	
	std::cout<<"Outputing for interpolation of the pressure force values\n";
	std::ofstream body_nodes;
	parameterDB  &db = *NavierStokesSolver::paramDB;
	std::string folder = db["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	out << folder << "/interp_test_force.csv";
	body_nodes.open(out.str().c_str());
	body_nodes <<	"x1\t"
					"x2\t"
					"x3\t"
					"x4\t"
					"y1\t"
					"y2\t"
					"y3\t"
					"y4\t"
					"q1\t"
					"q2\t"
					"q3\t"
					"q4\t"
					"p\t"
					"body_node_1_x\t"
					"body_node_1_y\t"
					"body_node_2_x\t"
					"body_node_2_y\t"
					"px\t"
					"py\n"
					;
	for (int i=0;  i<totalPoints;  i++)
	{
		//std::cout<<I<<"\t"<<J<<"\t"<<iv<<"\n";
		body_nodes << B.x1[i] <<"\t";
		body_nodes << B.x2[i] <<"\t";
		body_nodes << B.x3[i] <<"\t";
		body_nodes << B.x4[i] <<"\t";
		body_nodes << B.y1[i] <<"\t";
		body_nodes << B.y2[i] <<"\t";
		body_nodes << B.y3[i] <<"\t";
		body_nodes << B.y4[i] <<"\t";
		body_nodes << B.q1[i] <<"\t";
		body_nodes << B.q2[i] <<"\t";
		body_nodes << B.q3[i] <<"\t";
		body_nodes << B.q4[i] <<"\t";
		body_nodes << B.force_pressure[i] <<"\t";
		body_nodes << B.x[i]<<"\t";
		body_nodes << B.y[i]<<"\t";
		if (i == totalPoints-1)
		{
			body_nodes << B.x[0]<<"\t";
			body_nodes <<B.y[0]<<"\t";
		}
		else
		{
			body_nodes << B.x[i+1]<<"\t";
			body_nodes <<B.y[i+1] <<"\t";
		}
		body_nodes << B.point_x[i]<<"\t";
		body_nodes << B.point_y[i]<<"\t";
		body_nodes << B.point2_x[i]<<"\t";
		body_nodes << B.point2_y[i]<<"\t";
		body_nodes << B.point3_x[i]<<"\t";
		body_nodes << B.point3_y[i]<<"\n";
	}
	body_nodes.close();
	
	
	
	
}