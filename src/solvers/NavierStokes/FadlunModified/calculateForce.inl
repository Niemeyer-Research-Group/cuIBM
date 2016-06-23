/***************************************************************************//**
 * \file calculateForce.inl
 * \author Anush Krishnan (anush@bu.edu),
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke kernels that will calculate the force on the immersed body
 */

#include <solvers/NavierStokes/FadlunModified/kernels/calculateForce.h>

/**
 * \brief Calculates forces acting on an immersed body (on the device).
 *
 * Uses the control volume approach explained by Lai and Peskin (2000).
 * This is a general method that can be used with any immersed boundary method.
 * It uses only the velocity and pressure fields to calculate the forces, and
 * does not involve any body forces on the immersed boundary.
 * Currently works only for one body.
 */
void fadlunModified::calculateForce()
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

	int		*tagsIn_r	= thrust::raw_pointer_cast( &(tagsIn[0]) );

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
	kernels::dragUnsteady <<<dimGridX, dimBlock>>> (FxU_r, u_r, uold_r, tagsIn_r, dx, dy, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	fxx = thrust::reduce(FxX.begin(), FxX.end());
	fxy = thrust::reduce(FxY.begin(), FxY.end());
	fxu = thrust::reduce(FxU.begin(), FxU.end());
	B.forceX =  fxx + fxy + fxu;
	//std::cout<<timeStep<<"\tFx\t"<<B.forceX<<"\tFxX\t"<<fxx<<"\tFxY\t"<<fxy<<"\tFxU\t"<<fxu<<"\n";

	// Calculating lift
	cusp::array1d<double, cusp::device_memory>
		FyX(B.numCellsY[0]+1),
		FyY(B.numCellsX[0]),
		FyU((B.numCellsX[0]+1)*B.numCellsY[0]);

	double *FyX_r = thrust::raw_pointer_cast(&FyX[0]),
	     *FyY_r = thrust::raw_pointer_cast(&FyY[0]),
	     *FyU_r = thrust::raw_pointer_cast(&FyU[0]);

	kernels::liftLeftRight <<<dimGrid, dimBlock>>> (FyX_r, u_r, nu, dx, dy, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	kernels::liftBottomTop <<<dimGrid, dimBlock>>> (FyY_r, u_r, pressure_r, nu, dx, dy, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	dim3 dimGridY( int( ( B.numCellsX[0]*(B.numCellsY[0]+1)-0.5 )/blockSize )+1, 1 );
	kernels::liftUnsteady <<<dimGridY, dimBlock>>> (FyU_r, u_r, uold_r, tagsIn_r, dx, dy, dt, \
	                                                nx, ny, B.startI[0], B.startJ[0], B.numCellsX[0], B.numCellsY[0]);

	//B.forceY = thrust::reduce(FyX.begin(), FyX.end()) + thrust::reduce(FyY.begin(), FyY.end()) + thrust::reduce(FyU.begin(), FyU.end());
	//if (timeStep == 138 || timeStep == 139 || timeStep == 146 || timeStep == 147 ||  timeStep == 159 || timeStep == 160)
		//print_forces(FyX, FyY, FyU);
	//std::cout<<timeStep<<"\t";
	//std::cout<<"Fy: "<< B.forceY<<"\t";
	//std::cout<<"FyX: "<< thrust::reduce(FyX.begin(), FyX.end())<<"\t";
	//std::cout<<"FyY: "<< thrust::reduce(FyY.begin(), FyY.end())<<"\t";
	//std::cout<<"FyU: "<< thrust::reduce(FyU.begin(), FyU.end())<<"\n";
	B.forceY = thrust::reduce(FyX.begin(), FyX.end()) + thrust::reduce(FyY.begin(), FyY.end()) + thrust::reduce(FyU.begin(), FyU.end());

}
/*
void NavierStokesSolver::print_forces(cusp::array1d<double, cusp::device_memory> FyX, cusp::array1d<double, cusp::device_memory> FyY, cusp::array1d<double, cusp::device_memory> FyU)
{
	logger.startTimer("output");

	std::ofstream myfile;
	std::string folder = (*paramDB)["inputs"]["caseFolder"].get<std::string>();
	std::stringstream out;
	std::stringstream convert; convert << "/output/FyX"<<timeStep << ".csv";
	std::string folder_name = convert.str();
	out<<folder<<folder_name;
	myfile.open(out.str().c_str());
	myfile<<"FyX\n";
	for (int i = 0; i < B.numCellsY[0]+1; i++)
	{
		myfile<<FyX[i]<<"\n";
	}
	myfile.close();
	std::cout<<"printed FyX\n";
	
	std::ofstream myfile2;
	std::stringstream out2;
	std::stringstream convert2; convert2 << "/output/FyY"<<timeStep << ".csv";
	std::string folder_name2 = convert2.str();
	out2<<folder<<folder_name2;
	myfile2.open(out2.str().c_str());
	myfile2<<"FyY\n";
	for (int i = 0; i < B.numCellsX[0]; i++)
	{
		myfile2<<FyY[i]<<"\n";
	}
	myfile2.close();
	std::cout<<"printed FyY\n";
	
	std::ofstream myfile3;
	std::stringstream out3;
	std::stringstream convert3; convert3 << "/output/FyU"<<timeStep << ".csv";
	std::string folder_name3 = convert3.str();
	out3<<folder<<folder_name3;
	myfile3.open(out3.str().c_str());
	myfile3<<"FyU\n";
	int i,j,idx;
	for (int J=0; J<domInfo->ny; J++)
	{
		for (int I=0; I<domInfo->nx; I++)
		{
			i = I-B.startI[0];
			j = J-B.startJ[0];
			idx = B.numCellsX[0]*j + i;
			if (i >= 0 && j >= 0 && i < B.numCellsX[0]-1 && j < B.numCellsY[0]-1)
			{
				myfile3<<FyU[idx]<<"\t";	
			}
			else
				myfile3<<"0"<<"\t";
			
		}
		myfile3<<"\n";
	}
	myfile3.close();
	std::cout<<"printed FyU\n";
	logger.stopTimer("output");
}*/