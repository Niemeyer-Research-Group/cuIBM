/***************************************************************************//**
 * \file intermediateVelocity.inl
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */

#include <solvers/NavierStokes/NavierStokes/kernels/intermediateVelocity.h>
#include <solvers/NavierStokes/FadlunModified/kernels/intermediateVelocity.h>
#include <solvers/NavierStokes/FadlunModified/kernels/LHS1.h>
#include <solvers/NavierStokes/NavierStokes/kernels/LHS1.h>
#include <solvers/NavierStokes/NavierStokes/kernels/N.h>
#include <solvers/NavierStokes/NavierStokes/kernels/L.h>

void fadlunModified::generateRHS1()
{
	double	*u_r    = thrust::raw_pointer_cast( &(NavierStokesSolver::u[0]) ),
			*uold_r = thrust::raw_pointer_cast( &(NavierStokesSolver::uold[0]) ),
			*N_r    = thrust::raw_pointer_cast( &(NavierStokesSolver::N[0]) ),
			*Nold_r = thrust::raw_pointer_cast( &(NavierStokesSolver::Nold[0]) ),
			*L_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::L[0]) ),
			*rhs_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::rhs1[0]) ),
			*bc1_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::bc1[0]) ),
			*x_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::domInfo->x[0]) ),
			*y_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::domInfo->y[0]) ),
			*xp_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::bc[XPLUS][0]) ),
			*distance_from_intersection_to_node_r	= thrust::raw_pointer_cast( &(distance_from_intersection_to_node[0]) ),
			*distance_between_nodes_at_IB_r	= thrust::raw_pointer_cast( &(distance_between_nodes_at_IB[0]) ),
			*uv_r	= thrust::raw_pointer_cast( &(uv[0]) );
	
	int	*tags_r	= thrust::raw_pointer_cast( &(tags[0]) ),
		*tagsIn_r=thrust::raw_pointer_cast( &(tagsIn[0]) );
	
	parameterDB  &db = *NavierStokesSolver::paramDB;
	double	nu = db["flow"]["nu"].get<double>();
	double	dt = db["simulation"]["dt"].get<double>();

	int nx = NavierStokesSolver::domInfo ->nx,
		ny = NavierStokesSolver::domInfo ->ny;
	
	//set correct grid and block size
	const int blocksize = 256;
	dim3 dimGridUV( int( ((nx-1)*ny+nx*(ny-1)-0.5)/blocksize ) + 1, 1);
	dim3 dimBlockUV(blocksize, 1);
	dim3 dimGridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 dimBlockU(blocksize, 1);
	dim3 dimGridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 dimBlockV(blocksize, 1);

	updateRobinBoundary();

	NavierStokesSolver::Nold = NavierStokesSolver::N;
	NavierStokesSolver::generateN();

	NavierStokesSolver::generateL();

	NavierStokesSolver::generateBC1();

	kernels::generateRHS<<<dimGridUV,dimBlockUV>>>(rhs_r, L_r, Nold_r, N_r, u_r, bc1_r, dt, nx, ny);

	kernels::updateRHS1forIBX<<<dimGridU,dimBlockU>>>(tags_r, tagsIn_r, rhs_r, distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r, uv_r, nx, ny);
	kernels::updateRHS1forIBY<<<dimGridV,dimBlockV>>>(tags_r, tagsIn_r, rhs_r, distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r, uv_r, nx, ny);
}

void fadlunModified::updateRobinBoundary()
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

void fadlunModified::generateLHS1()
{
	NavierStokesSolver::logger.startTimer("LHS1");

	double 	*val_r	= thrust::raw_pointer_cast( &(LHS1.values[0])),
			*dx_r	= thrust::raw_pointer_cast( &(domInfo->dx[0])),
			*dy_r	= thrust::raw_pointer_cast( &(domInfo->dy[0])),
			*distance_from_intersection_to_node_r	= thrust::raw_pointer_cast( &(distance_from_intersection_to_node[0]) ),
			*distance_between_nodes_at_IB_r	= thrust::raw_pointer_cast( &(distance_between_nodes_at_IB[0]) );

	int 	*row_r	= thrust::raw_pointer_cast( &(LHS1.row_indices[0])),
			*col_r	= thrust::raw_pointer_cast( &(LHS1.column_indices[0])),
			*tags_r	= thrust::raw_pointer_cast( &(tags[0]) ),
			*tags2_r= thrust::raw_pointer_cast( &(tags2[0]) ),
			*tagsIn_r= thrust::raw_pointer_cast( &(tagsIn[0]) );

	double	nu = (*paramDB)["flow"]["nu"].get<double>();
	double	dt = (*paramDB)["simulation"]["dt"].get<double>();
	int 	nx = domInfo->nx,
			ny = domInfo->ny;

	const int blocksize = 256;
	dim3 gridU( int( ((nx-1)*ny-0.5)/blocksize ) +1, 1);
	dim3 blockU(blocksize, 1);
	dim3 gridV( int( (nx*(ny-1)-0.5)/blocksize ) +1, 1);
	dim3 blockV(blocksize, 1);

	kernels::LHS_mid_X<<<gridU,blockU>>>(row_r, col_r, val_r, tags_r, tags2_r, tagsIn_r, distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r, dx_r, dy_r, dt, nu, nx, ny);
	kernels::LHS_mid_Y<<<gridV,blockV>>>(row_r, col_r, val_r, tags_r, tags2_r, tagsIn_r, distance_from_intersection_to_node_r, distance_between_nodes_at_IB_r, dx_r, dy_r, dt, nu, nx, ny);

	kernels::LHS_BC_X<<<gridU,blockU>>>(row_r, col_r, val_r, dx_r, dy_r, dt, nu, nx, ny);
	kernels::LHS_BC_Y<<<gridV,blockV>>>(row_r, col_r, val_r, dx_r, dy_r, dt, nu, nx, ny);

	NavierStokesSolver::logger.stopTimer("LHS1");
}