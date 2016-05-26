/***************************************************************************//**
 * \file intermediateVelocity.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */

#include <solvers/NavierStokes/NavierStokes/kernels/intermediateVelocity.h> //generaterhs1
#include <solvers/NavierStokes/luoIBM/kernels/intermediateVelocity.h>//updaterhs1_luo, zeroInside
#include <solvers/NavierStokes/luoIBM/kernels/LHS1.h> //generatelhs_luo _mid
#include <solvers/NavierStokes/NavierStokes/kernels/LHS1.h> //lhs_bc
#include <solvers/NavierStokes/NavierStokes/kernels/N.h> //genN
#include <solvers/NavierStokes/NavierStokes/kernels/L.h> //genL
#include <solvers/NavierStokes/FadlunModified/kernels/intermediateVelocity.h> //updateboundary
#include <solvers/NavierStokes/luoIBM/kernels/biLinearInterpolation.h> //interpolate

void luoIBM::generateRHS1()
{
	double	*u_r    = thrust::raw_pointer_cast( &(u[0]) ),
			*uold_r = thrust::raw_pointer_cast( &(uold[0]) ),
			*N_r    = thrust::raw_pointer_cast( &(N[0]) ),
			*Nold_r = thrust::raw_pointer_cast( &(Nold[0]) ),
			*L_r	= thrust::raw_pointer_cast( &(L[0]) ),
			*rhs_r	= thrust::raw_pointer_cast( &(rhs1[0]) ),
			*bc1_r	= thrust::raw_pointer_cast( &(bc1[0]) ),
			*x_r	= thrust::raw_pointer_cast( &(domInfo->x[0]) ),
			*y_r	= thrust::raw_pointer_cast( &(domInfo->y[0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(bc[XPLUS][0]) ),
			*distance_from_intersection_to_node_r	= thrust::raw_pointer_cast( &(distance_from_intersection_to_node[0]) ),
			*distance_between_nodes_at_IB_r	= thrust::raw_pointer_cast( &(distance_between_nodes_at_IB[0]) ),
			*uv_r	= thrust::raw_pointer_cast( &(uv[0]) );
	
	int	*hybridTagsUV_r	= thrust::raw_pointer_cast( &(hybridTagsUV[0]) ),
		*ghostTagsUV_r	= thrust::raw_pointer_cast( &(ghostTagsUV[0]) );
	
	parameterDB  &db = *paramDB;
	double	nu = db["flow"]["nu"].get<double>();
	double	dt = db["simulation"]["dt"].get<double>();

	int nx = domInfo ->nx,
		ny = domInfo ->ny;
	
	//set correct grid and block size
	const int blocksize = 256;
	dim3 dimGridUV( int( ((nx-1)*ny+nx*(ny-1)-0.5)/blocksize ) + 1, 1);
	dim3 dimBlockUV(blocksize, 1);
	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	//update right boundary of the domain for the convective boundary condition
	updateRobinBoundary();

	//set Nold to N
	Nold = N;
	
	//interoolate u and v to image points, then again to ghost nodes
	preRHS1Interpolation();
	
	//calculate explicit advection terms
	generateN();

	//calculate explicit diffusion terms
	generateL();

	//calculate boundary terms
	generateBC1();

	//sum rhs components
	kernels::generateRHS<<<dimGridUV,dimBlockUV>>>(rhs_r, L_r, Nold_r, N_r, u_r, bc1_r, dt, nx, ny);
}

void luoIBM::updateRobinBoundary()
{
	NavierStokesSolver::logger.startTimer("update Boundary");
	double	*u_r    = thrust::raw_pointer_cast( &(u[0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(bc[XPLUS][0]) ),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dx[0]) );
	double 	Uinf = 1, //need a better way to enforce these, ie read from yaml file
			Vinf = 1;
	int 	nx = domInfo ->nx,
			ny = domInfo ->ny;
	double	dt = (*paramDB)["simulation"]["dt"].get<double>();

	const int blocksize = 256;
	dim3 dimGridBCX( int(ny/blocksize) + 1, 1);
	dim3 dimGridBCY( int(nx/blocksize) + 1, 1);
	dim3 dimBlockBC(blocksize, 1);

	kernels::updateBoundaryX<<<dimGridBCX,dimBlockBC>>>(u_r, xp_r, dx_r, dt, Uinf, nx, ny);
	kernels::updateBoundaryY<<<dimGridBCY,dimBlockBC>>>(u_r, xp_r, dx_r, dt, Vinf, nx, ny);
	NavierStokesSolver::logger.stopTimer("update Boundary");
}

void luoIBM::generateLHS1()
{
	NavierStokesSolver::logger.startTimer("LHS1");

	double 	*val_r	= thrust::raw_pointer_cast( &(LHS1.values[0])),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dx[0])),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dy[0]));

	int 	*row_r				= thrust::raw_pointer_cast( &(LHS1.row_indices[0])),
			*col_r				= thrust::raw_pointer_cast( &(LHS1.column_indices[0])),
			*ghostTagsUV_r		= thrust::raw_pointer_cast( &(ghostTagsUV[0]) );

	double	nu = (*paramDB)["flow"]["nu"].get<double>();
	double	dt = (*paramDB)["simulation"]["dt"].get<double>();
	int 	nx = domInfo->nx,
			ny = domInfo->ny;

	const int blocksize = 256;
	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);

	kernels::LHS1_mid_luo_X<<<gridU,blockU>>>(row_r, col_r, val_r, ghostTagsUV_r, dx_r, dy_r, dt, nu, nx, ny);//flag lhs1 mid luo is the same as lhs1 mid nobody from NS
	kernels::LHS1_mid_luo_Y<<<gridV,blockV>>>(row_r, col_r, val_r, ghostTagsUV_r, dx_r, dy_r, dt, nu, nx, ny);

	kernels::LHS_BC_X<<<gridU,blockU>>>(row_r, col_r, val_r, dx_r, dy_r, dt, nu, nx, ny);
	kernels::LHS_BC_Y<<<gridV,blockV>>>(row_r, col_r, val_r, dx_r, dy_r, dt, nu, nx, ny);

	NavierStokesSolver::logger.stopTimer("LHS1");
}

void luoIBM::preRHS1Interpolation()
{
	double	*u_r = thrust::raw_pointer_cast ( &(u[0]) ),
			*bx_r		= thrust::raw_pointer_cast ( &(B.x[0]) ),//not sure if these are on the host or not
			*by_r		= thrust::raw_pointer_cast ( &(B.y[0]) ),
			*uB_r		= thrust::raw_pointer_cast ( &(B.uB[0]) ),
			*vB_r		= thrust::raw_pointer_cast ( &(B.vB[0]) ),
			*yu_r		= thrust::raw_pointer_cast ( &(domInfo->yu[0]) ),
			*xu_r		= thrust::raw_pointer_cast ( &(domInfo->xu[0]) ),
			*yv_r		= thrust::raw_pointer_cast ( &(domInfo->yv[0]) ),
			*xv_r		= thrust::raw_pointer_cast ( &(domInfo->xv[0]) ),
			*body_intercept_x_r = thrust::raw_pointer_cast( &(body_intercept_x[0]) ),
			*body_intercept_y_r = thrust::raw_pointer_cast( &(body_intercept_y[0]) ),
			*image_point_x_r = thrust::raw_pointer_cast( &(image_point_x[0]) ),
			*image_point_y_r = thrust::raw_pointer_cast( &(image_point_y[0]) );
	
	double	*x1_r = thrust::raw_pointer_cast ( &(x1[0]) ),
			*x2_r = thrust::raw_pointer_cast ( &(x2[0]) ),
			*x3_r = thrust::raw_pointer_cast ( &(x3[0]) ),
			*x4_r = thrust::raw_pointer_cast ( &(x4[0]) ),
			*y1_r = thrust::raw_pointer_cast ( &(y1[0]) ),
			*y2_r = thrust::raw_pointer_cast ( &(y2[0]) ),
			*y3_r = thrust::raw_pointer_cast ( &(y3[0]) ),
			*y4_r = thrust::raw_pointer_cast ( &(y4[0]) ),
			*q1_r = thrust::raw_pointer_cast ( &(q1[0]) ),
			*q2_r = thrust::raw_pointer_cast ( &(q2[0]) ),
			*q3_r = thrust::raw_pointer_cast ( &(q3[0]) ),
			*q4_r = thrust::raw_pointer_cast ( &(q4[0]) ),
			*ip_u_r = thrust::raw_pointer_cast( &(ip_u[0]) );
	
	int 	*ghostTagsUV_r		= thrust::raw_pointer_cast ( &(ghostTagsUV[0]) );
		
	parameterDB  &db = *paramDB;
	double	nu = db["flow"]["nu"].get<double>();
	double	dt = db["simulation"]["dt"].get<double>();

	int nx = domInfo ->nx,
		ny = domInfo ->ny,
		totalPoints = B.totalPoints,
		start_index_i = B.startI[0],
		start_index_j = B.startJ[0],
		width_i = B.numCellsX[0],
		height_j = B.numCellsY[0],
		end_index_i = start_index_i + width_i,
		end_index_j = start_index_j + height_j;
	
	const int blocksize = 256;
	dim3 grid( int( (width_i*height_j-0.5)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	
	kernels::interpolateVelocityToGhostNodeX<<<grid,block>>>(u_r, ghostTagsUV_r, bx_r, by_r, uB_r, yu_r, xu_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													start_index_i, start_index_j, end_index_i, end_index_j, nx, ny, totalPoints,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,ip_u_r);
	kernels::interpolateVelocityToGhostNodeY<<<grid,block>>>(u_r, ghostTagsUV_r, bx_r, by_r, vB_r, yv_r, xv_r,
													body_intercept_x_r, body_intercept_y_r, image_point_x_r, image_point_y_r,
													start_index_i, start_index_j, end_index_i, end_index_j, nx, ny, totalPoints,
													x1_r,x2_r,x3_r,x4_r,y1_r,y2_r,y3_r,y4_r,q1_r,q2_r,q3_r,q4_r,ip_u_r);
	kernels::interpolateVelocityToHybridNodeX<<<grid,block>>>()
			
	testInterpX();
	//testInterpY();
	
}

void luoIBM::testInterpX()
{
	std::cout<<"Outputing for interpolation of the u values\n";
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
	out << folder << "/interp_testX.csv";
	body_nodes.open(out.str().c_str());
	body_nodes << "BN_X1\t"
					"BN_Y1\t"
					"BN_X2\t"
					"BN_Y2\t"
					"GN_X\t"
					"GN_Y\t"
					"BI_X\t"
					"BI_Y\t"
					"IP_X\t"
					"IP_Y\t"
					"x1\t"
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
					"GN_U\t"
					"ip_u\n"
			;
	for (int J=start_index_j;  J<end_index_j;  J++)
	{
		for (int I=start_index_i;  I<end_index_i;  I++)
		{
			iu = J*(nx-1) + I;
			//if (ghostTagsUV[iu] >0)//for inside
			if (hybridTagsUV[iu] >0)//for outside
			{
				body_nodes << x1_ip[iu]<<"\t";
				body_nodes << y1_ip[iu]<<"\t";
				body_nodes << x2_ip[iu]<<"\t";
				body_nodes << y2_ip[iu]<<"\t";
				body_nodes << domInfo->xu[I]<<"\t";
				body_nodes << domInfo->yu[J]<<"\t";
				body_nodes << body_intercept_x[iu] <<"\t";
				body_nodes << body_intercept_y[iu] <<"\t";
				body_nodes << image_point_x[iu] <<"\t";
				body_nodes << image_point_y[iu] <<"\t";
				body_nodes << x1[iu] <<"\t";
				body_nodes << x2[iu] <<"\t";
				body_nodes << x3[iu] <<"\t";
				body_nodes << x4[iu] <<"\t";
				body_nodes << y1[iu] <<"\t";
				body_nodes << y2[iu] <<"\t";
				body_nodes << y3[iu] <<"\t";
				body_nodes << y4[iu] <<"\t";
				body_nodes << q1[iu] <<"\t";
				body_nodes << q2[iu] <<"\t";
				body_nodes << q3[iu] <<"\t";
				body_nodes << q4[iu] <<"\t";
				body_nodes << u[iu] <<"\t";
				body_nodes << ip_u[iu]<<"\n";
			}
		}
	}
	body_nodes.close();
}

void luoIBM::testInterpY()
{
	std::cout<<"Outputing for interpolation of the v values\n";
	int iv;
	int nx = NavierStokesSolver::domInfo->nx,
		ny = domInfo->ny,
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
	out << folder << "/interp_testY.csv";
	body_nodes.open(out.str().c_str());
	body_nodes <<	"BN_X1\t"
					"BN_Y1\t"
					"BN_X2\t"
					"BN_Y2\t"
					"GN_X\t"
					"GN_Y\t"
					"BI_X\t"
					"BI_Y\t"
					"IP_X\t"
					"IP_Y\t"
					"x1\t"
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
					"GN_U\t"
					"ip_u\n";
	for (int J=start_index_j;  J<end_index_j;  J++)
	{
		for (int I=start_index_i;  I<end_index_i;  I++)
		{
			iv = J*nx + I  +  ny*(nx-1);
			//if (ghostTagsUV[iv] >0)//for inside
			if (hybridTagsUV[iv] >0)//for outside
			{
				//std::cout<<I<<"\t"<<J<<"\t"<<iv<<"\n";
				body_nodes << x1_ip[iv]<<"\t";
				body_nodes << y1_ip[iv]<<"\t";
				body_nodes << x2_ip[iv]<<"\t";
				body_nodes << y2_ip[iv]<<"\t";
				body_nodes << domInfo->xv[I]<<"\t";
				body_nodes << domInfo->yv[J]<<"\t";
				body_nodes << body_intercept_x[iv] <<"\t";
				body_nodes << body_intercept_y[iv] <<"\t";
				body_nodes << image_point_x[iv] <<"\t";
				body_nodes << image_point_y[iv] <<"\t";
				body_nodes << x1[iv] <<"\t";
				body_nodes << x2[iv] <<"\t";
				body_nodes << x3[iv] <<"\t";
				body_nodes << x4[iv] <<"\t";
				body_nodes << y1[iv] <<"\t";
				body_nodes << y2[iv] <<"\t";
				body_nodes << y3[iv] <<"\t";
				body_nodes << y4[iv] <<"\t";
				body_nodes << q1[iv] <<"\t";
				body_nodes << q2[iv] <<"\t";
				body_nodes << q3[iv] <<"\t";
				body_nodes << q4[iv] <<"\t";
				body_nodes << u[iv]  <<"\t";
				body_nodes << ip_u[iv]<<"\n";
			}
		}
	}
	body_nodes.close();
}