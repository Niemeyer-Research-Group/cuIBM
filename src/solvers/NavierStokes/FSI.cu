/*

void FSI::moveBodySC()
{
	parameterDB  &db = *NavierStokesSolver::paramDB;
	NavierStokesSolver::calculateForce();

	double *y_r	= thrust::raw_pointer_cast( &(NavierStokesSolver::B.y[0]) ),
		   *vB_r= thrust::raw_pointer_cast( &(NavierStokesSolver::B.vB[0]) );
	double	Cy	= NavierStokesSolver::B.forceY*2.0,
			U	= NavierStokesSolver::bc[XMINUS][0],
			Mred= 2.0,
			Ured= db["simulation"]["Ured"].get<double>(),
			dt	= db["simulation"]["dt"].get<double>(),
			totalPoints=NavierStokesSolver::B.totalPoints,
			vold= NavierStokesSolver::B.centerVelocityV0,
			yold= NavierStokesSolver::B.midY0,
			vnew,
			ynew,
			relaxation_coeff = 0.75;

	double a = dt*M_PI*M_PI*4/(Ured*Ured),
		   b = dt*dt*2*M_PI*M_PI/(Ured*Ured);

	//calc updated velocity
	vnew = (vold - a*(yold+ dt/2*vold) + dt*Cy/2/Mred)/(1+b);
	vnew = relaxation_coeff * vnew + (1-relaxation_coeff) * NavierStokesSolver::B.centerVelocityV; //relax velocity
	ynew = yold + dt/2*(vnew + vold);

	NavierStokesSolver::B.centerVelocityV = vnew;
	NavierStokesSolver::B.midY = ynew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::update_body_viv<<<grid,block>>>(y_r, vB_r, ynew-yold, vnew, totalPoints);
}

void FSI::SC()
{
	NavierStokesSolver::B.centerVelocityV0 = NavierStokesSolver::B.centerVelocityV;
	NavierStokesSolver::B.midY0 = NavierStokesSolver::B.midY;
	int count = 0;
	do
	{
		NavierStokesSolver::B.forceYk[0] = NavierStokesSolver::B.forceY;
		NavierStokesSolver::generateRHS1();
		NavierStokesSolver::solveIntermediateVelocity();

		NavierStokesSolver::generateRHS2();
		NavierStokesSolver::solvePoisson();

		NavierStokesSolver::velocityProjection();
		//Release the body after a certain timestep
		if (NavierStokesSolver::timeStep >= (*NavierStokesSolver::paramDB)["simulation"]["startStep"].get<int>())
		{
			moveBodySC();
			updateSolver();
		}
		count += 1;
	}
	while (fabs(NavierStokesSolver::B.forceY- NavierStokesSolver::B.forceYk[0]) > 0.0001);
	if (count > 1)
		std::cout<<count<<"\n";
	std::cout<<NavierStokesSolver::timeStep<<"\n";

	NavierStokesSolver::timeStep++;

}*/
