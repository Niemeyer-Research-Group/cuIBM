/***************************************************************************//**
 * \file intermediatePressure.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the right hand side of the poission equation
 */

#include "intermediatePressure.h"

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void intermediatePressure_luo(double *rhs2, double *uhat, int *ghostTagsP, int *hybridTagsP, int *ghostTagsUV, double *distance_from_u_to_body, double *distance_from_v_to_body, double *ym, double *yp, double *xm, double *xp, double *dx, double *dy, int nx, int ny)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*ny)
		return;
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		iu	= (nx-1)*J + I,
		iv	= (nx-1)*ny  +  nx*J +I;

	double temp = 0;

	//Outside immersed body
	if (hybridTagsP[ip] != -1)
	{
		//EAST
		//check if east pressure node is outside of the body
		if (ghostTagsP[ip+1] == -1)
		{
			if (distance_from_u_to_body[ip] > dx[I]/2 && distance_from_u_to_body[ip] < dx[I])
			{
				temp -= uhat[iu]/distance_from_u_to_body[ip];
			}
			else
				temp -= uhat[iu]/dx[I];
		}

		//WEST
		//check if west pressure node is outside of the body
		if (ghostTagsP[ip-1] == -1)
		{
			if (distance_from_u_to_body[ip] > dx[I]/2 && distance_from_u_to_body[ip] < dx[I])
			{
				temp += uhat[iu-1]/distance_from_u_to_body[ip];
			}
			else
				temp += uhat[iu-1]/dx[I];
		}

		//NORTH
		//check if north pressure node is outside of the body
		if (ghostTagsP[ip+nx] == -1)
		{
			if (distance_from_v_to_body[ip] > dy[J]/2 && distance_from_v_to_body[ip] < dy[J])
			{
				temp -= uhat[iv]/distance_from_v_to_body[ip];
			}
			else
				temp -= uhat[iv]/dy[J];
		}

		//SOUTH
		//check if south velocity node is outside of the body
		if (ghostTagsP[ip-nx] == -1)
		{
			if (distance_from_v_to_body[ip] > dy[J]/2 && distance_from_v_to_body[ip] < dy[J])
			{
				temp += uhat[iv-nx]/distance_from_v_to_body[ip];
			}
			else
				temp += uhat[iv-nx]/dy[J];
		}


	}
	//end outside immersed body
	//if just inside body
	else if (ghostTagsP[ip] > 0)
	{
		//EAST
		if (ghostTagsP[ip+1] == 0)
			temp -= uhat[iu]/dx[I];

		//WEST
		if (ghostTagsP[ip-1] == 0)
			temp += uhat[iu - 1]/dx[I];

		//NORTH
		if (ghostTagsP[ip+nx] == 0)
			temp -= uhat[iv]/dy[J];

		//SOUTH
		if (ghostTagsP[ip-nx] == 0)
			temp += uhat[iv-nx]/dy[J];

	}
	//end just inside body
	//everywhere else
	else
	{
		//EAST
		//if not on the east wall and east is outside the body, add east term
		if (I != nx-1)//not at east boundry
			temp -= uhat[iu]/dx[I];
		else if (I == nx-1)//at east boundry
			temp -= xp[J]/dx[I];

		//WEST
		//if not on west wall and west is outside the body, add west term
		if (I != 0)//not at west boundary
			temp += uhat[iu - 1]/dx[I];
		else if (I == 0)//at the west boundary
			temp += xm[J]/dx[I];

		//NORTH
		//if not on north wall and north is outside the body, add north term
		if (J != ny-1)//not at north boundry
			temp -= uhat[iv]/dy[J];
		else if (J == ny-1)//at north boundry
			temp -= yp[(nx-1)+I]/dy[J];

		//SOUTH
		//if not on south wall and south is outside the body, add south term
		if (J != 0)//not at south boundry
			temp += uhat[iv-nx]/dy[J];
		else if (J == 0)//at south boundry
			temp += ym[(nx-1)+I]/dy[J];

	}//end everywhere else

	rhs2[ip] = temp;
}
}
