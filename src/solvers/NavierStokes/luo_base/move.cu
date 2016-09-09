/***************************************************************************//**
 * \file
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief functions to invoke the kernels that setup the intermediate velocity solve
 */

#include <solvers/NavierStokes/luo_base.h>

#include <solvers/NavierStokes/luo_base/kernels/structure.h> //update_body_viv

void luo_base::set_movement()
{
	double	t = dt*timeStep,
			f = B.xfrequency,
			xCoeff = B.xCoeff,
			uCoeff = B.uCoeff,
			xPhase = B.xPhase,
			uPhase = B.uPhase,
			totalPoints=B.totalPoints,
			xold= B.midX,
			unew,
			xnew;

	//xnew = -1/(2*M_PI)*sin(2*M_PI*f*t);
	//unew = -f*cos(2*M_PI*f*t);
	xnew = xCoeff*sin(2*M_PI*f*t + xPhase);
	unew = uCoeff*cos(2*M_PI*f*t + uPhase);

	B.centerVelocityU = unew;
	B.midX = xnew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	B.uBk = B.uB;
	kernels::update_body_viv<<<grid,block>>>(B.x_r, B.uB_r, xnew-xold, unew, totalPoints);
}

void luo_base::viv_movement()
{
	double	Cy	= B.forceY*2.0,
			U	= bc[XMINUS][0],
			Mred= 2.0,
			Ured= (*paramDB)["flow"]["Ured"].get<double>(),
			totalPoints=B.totalPoints,
			vold= B.centerVelocityV,
			yold= B.midY,
			vnew,
			ynew;

	double a = dt*M_PI*M_PI*4/(Ured*Ured),
		   b = dt*dt*2*M_PI*M_PI/(Ured*Ured);

	//calc updated velocity
	vnew = (vold - a*(yold+ dt/2*vold) + dt*Cy/2/Mred)/(1+b);
	ynew = yold + dt/2*(vnew + vold);
	std::cout<<"vnew\t"<<vnew<<"\n";
	std::cout<<"vold\t"<<vold<<"\n";
	std::cout<<"ynew\t"<<ynew<<"\n";
	std::cout<<"yold\t"<<yold<<"\n";
	B.centerVelocityV = vnew;
	B.midY = ynew;

	const int blocksize = 256;
	dim3 grid( int( (totalPoints)/blocksize ) +1, 1);
	dim3 block(blocksize, 1);
	kernels::update_body_viv<<<grid,block>>>(B.y_r, B.vB_r, ynew-yold, vnew, totalPoints);
}
