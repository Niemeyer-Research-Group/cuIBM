/***************************************************************************//**
 * \file LHS1.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief kernels to generate the left hand side for the intermediate velocity solve
 */

#include "LHS1.h"

namespace kernels
{
__global__
void LHS1_mid_iter_X(int *row, int *col, double *val, double *dx, double *dy, double dt, double nu, int nx, int ny,
					int *hybridTagsUV, int *ghostTagsUV, double *ns_rhs, double *interp_rhs, int *count,
					int *index1, int *index2, int *index3, int *index4,
					double *xu, double *yu, double *alpha, double *uB, //xu, yu not used
					double *q1coef, double *q2coef, double *q3coef, double *q4coef,
					double *q1, double *q2, double *q3, double *q4
					)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= (nx-1)*ny)
		return;
	int iu 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= iu % (nx-1),
		J	= iu / (nx-1);
	if (I == 0 || I == nx-2 || J == 0 || J == ny-1)
		return;

	//int numE = i*5;
	//			top row - corner    mid           sides    current row
	int numE = (nx-1)*4 - 2      + (J-1)*(5*(nx-1)  - 2) + I*5 - 1;

	double temp = 1;

	if (hybridTagsUV[iu]>0)
	{
		int interp_index[4] = {index1[iu], index2[iu], index3[iu], index4[iu]};
		double q[4] = {q1[iu], q2[iu], q3[iu], q4[iu]};
		double CInterp[4];
		double Cns[5];
		Cns[0] = -dt*nu/(dy[J+1]*(dy[J]+dy[J+1]));
		Cns[1] = -dt*nu/(dx[I]  *(dx[I]+dx[I+1]));
		Cns[2] = -dt*nu/(dy[J]  *(dy[J]+dy[J+1]));
		Cns[3] = -dt*nu/(dx[I]  *(dx[I]+dx[I-1]));
		Cns[4] = 1-Cns[0] - Cns[1] - Cns[2] - Cns[3];
		CInterp[0] = q1coef[iu];
		CInterp[1] = q2coef[iu];
		CInterp[2] = q3coef[iu];
		CInterp[3] = q4coef[iu];
		for (int l=0; l<4; l++)
		{
			Cns[l] = Cns[l]*(1-alpha[iu])/Cns[4];
			CInterp[l] = CInterp[l]*alpha[iu];
		}
		/*   0  1  2		NW  N   NE
		 *   3  4  5		W   P   E
		 *   6  7  8		SW  S   SE
		 */
		int stencil_index[9]    = {iu + (nx-1) - 1, iu + (nx-1), iu + (nx-1) + 1,
								   iu - 1         , iu         , iu + 1,
								   iu - (nx-1) - 1, iu - (nx-1), iu - (nx-1) + 1};
		double stencil[9] = {0, Cns[0], 0, Cns[3], 1, Cns[1], 0, Cns[2], 0};
		//combine ns and interp stencils
		bool stencil_used[9] = {false, true, false, true, true, true, false, true, false};
		for (int n=0;n<4;n++)
		{
			for (int m=0;m<9;m++)
			{
				if (stencil_index[m] == interp_index[n] && m != 4)
				{
					stencil[m] -= CInterp[n]; //flag should this be minus?
				}
			}
		}
		//add ns to sparse matrix
		for (int m = 0; m<9; m++)
		{
			if (stencil_used[m])
			{
				row[numE] = iu;
				col[numE] = stencil_index[m];
				val[numE] = stencil[m];
				numE++;
			}
		}
		ns_rhs[iu] = (1-alpha[iu])/Cns[4];
		interp_rhs[iu] = 0;
		//calc new numE
		numE = ny*(nx-1)*5 - ny*2 - (nx-1)*2    +   nx*(ny-1)*5 - nx*2 - (ny-1)*2 + count[iu]-1;
		//add interp corner to sparse matrix
		for (int n=0;n<4;n++)
		{
			for (int m=0;m<9;m++)
			{
				if (stencil_index[m] == interp_index[n] && !stencil_used[m])
				{
					row[numE] = iu;
					col[numE] = interp_index[n];
					val[numE] = -CInterp[n]; //this should be minus?
				}
				//else if(stencil_index[m] == interp_index[n] && stencil_used[m])
				else if(stencil_index[m] == interp_index[n] && interp_index[n] == iu)
					interp_rhs[iu] += CInterp[n]*q[n];
			}
		}
	}
	else if (ghostTagsUV[iu]>0)
	{
		int interp_index[4] = {index1[iu], index2[iu], index3[iu], index4[iu]};
		bool interp_in[4] = {false, false, false, false};
		int ns_index[5] = {iu + (nx-1), iu + 1, iu - (nx-1), iu -1, iu}; //n e s w p
		bool ns_overlap[5] = {false, false, false, false, true};
		double q[4] = {q1[iu], q2[iu], q3[iu], q4[iu]};
		double CInterp[4];
		CInterp[0] = q1coef[iu];
		CInterp[1] = q2coef[iu];
		CInterp[2] = q3coef[iu];
		CInterp[3] = q4coef[iu];
		//count the number of nodes the interp is using
		//find how which ns nodes are occupied
		int counter = 0;
		temp = 0;
		for (int l=0; l<4; l++)
		{
			if (ghostTagsUV[interp_index[l]]>0)
			{
				counter +=1;
				interp_in[l] = true;
			}
			for (int n=0; n<5; n++)
			{
				if (interp_index[l] == ns_index[n])
					ns_overlap[n] = true;
			}
		}
		//add center to matrix
		row[numE] = iu;
		col[numE] = iu;
		val[numE] = 1;
		numE++;
		//add real interp values to matrix
		for (int i=0; i<4; i++)
		{
			if (!interp_in[i] && interp_index[i] != iu)
			{
				row[numE] = iu;
				col[numE] = interp_index[i];
				val[numE] = CInterp[i];
				numE++;
			}
			else
			{
				temp -= CInterp[i] * q[i];
			}
		}
		//fill remainder of values
		int counter2 = 0;
		for (int i=0; i<5; i++)
		{
			if (counter2>=counter)
				break;
			if (ns_overlap[i]==false)
			{
				row[numE] = iu;
				col[numE] = ns_index[i];
				val[numE] = 0;
				numE++;
				counter2++;
			}
		}
		ns_rhs[iu] = 0;
		interp_rhs[iu] = 2*uB[0] + temp;//flag this doesn't account for the interpolation part
	}
	else
	{
	temp = 1 + 0.5*dt*nu*(1/(dx[I+1]*(dx[I+1]+dx[I])*0.5)) + 0.5*dt*nu*(1/(dx[I]*(dx[I+1]+dx[I])*0.5)) + 0.5*dt*nu*(1/(dy[J]*(dy[J+1]+dy[J])*0.5)) + 0.5*dt*nu*(1/(dy[J]*(dy[J-1]+dy[J])*0.5));
	//EAST
	row[numE] = iu;
	col[numE] = iu+1;
	val[numE] = -0.5*dt*nu*(1/(dx[I+1]*(dx[I+1]+dx[I])*0.5))/temp;
	numE++;

	//WEST
	row[numE] = iu;
	col[numE] = iu-1;
	val[numE] = -0.5*dt*nu*(1/(dx[I]*(dx[I+1]+dx[I])*0.5))/temp;
	numE++;

	//NORTH
	row[numE] = iu;
	col[numE] = iu+(nx-1);
	val[numE] = -0.5*dt*nu*(1/(dy[J]*(dy[J+1]+dy[J])*0.5))/temp;
	numE++;

	//SOUTH
	row[numE] = iu;
	col[numE] = iu-(nx-1);
	val[numE] = -0.5*dt*nu*(1/(dy[J]*(dy[J-1]+dy[J])*0.5))/temp;
	numE++;

	//CENTER
	row[numE] = iu;
	col[numE] = iu;
	val[numE] = 1;
	numE++;
	ns_rhs[iu] = 1/temp;
	interp_rhs[iu] = 0;
	}
}

__global__
void LHS1_mid_iter_Y(int *row, int *col, double *val, double *dx, double *dy, double dt, double nu, int nx, int ny,
						int *hybridTagsUV, int *ghostTagsUV, double *ns_rhs, double *interp_rhs, int *count,
						int *index1, int *index2, int *index3, int *index4,
						double *xv, double *yv, double *alpha, double *vB,
						double *q1coef, double *q2coef, double *q3coef, double *q4coef,
						double *q1, double *q2, double *q3, double *q4
						)
{
	if (threadIdx.x + blockDim.x * blockIdx.x >= nx*(ny-1))
		return;
	int ip 	= threadIdx.x + blockDim.x * blockIdx.x,
		I	= ip % nx,
		J	= ip / nx,
		iv = ip + (nx-1)*ny;
	if (I == 0 || I == nx-1 || J == 0 || J == ny-2)
		return;

	int numE = (nx-1)*ny*5 - 2*ny-2*(nx-1)  +  nx*4-2  + (J-1)*(nx*5 - 2) + I*5 - 1;
	double temp = 1;

	if (hybridTagsUV[iv]>0)
	{
		int interp_index[4] = {index1[iv], index2[iv], index3[iv], index4[iv]};
		double q[4] = {q1[iv], q2[iv], q3[iv], q4[iv]};
		double CInterp[4];
		double Cns[5];
		Cns[0] = -dt*nu/(dy[J+1]*(dy[J]+dy[J+1]));
		Cns[1] = -dt*nu/(dx[I]*(dx[I]+dx[I+1]));
		Cns[2] = -dt*nu/(dy[J]*(dy[J]+dy[J+1]));
		Cns[3] = -dt*nu/(dx[I]*(dx[I]+dx[I-1]));
		Cns[4] = 1-Cns[0] - Cns[1] - Cns[2] - Cns[3];
		CInterp[0] = q1coef[iv];
		CInterp[1] = q2coef[iv];
		CInterp[2] = q3coef[iv];
		CInterp[3] = q4coef[iv];
		for (int l=0; l<4; l++)
		{
			Cns[l] = Cns[l]*(1-alpha[iv])/Cns[4];
			CInterp[l] = CInterp[l]*alpha[iv];
		}
		/*   0  1  2		NW  N   NE
		 *   3  4  5		W   P   E
		 *   6  7  8		SW  S   SE
		 */
		int stencil_index[9]    = {iv + nx - 1, iv + nx, iv + nx + 1,
								   iv - 1     , iv     , iv + 1,
								   iv - nx - 1, iv - nx, iv - nx + 1};
		double stencil[9] = {0, Cns[0], 0, Cns[3], 1, Cns[1], 0, Cns[2], 0};
		//combine ns and interp stencils
		bool stencil_used[9] = {false, true, false, true, true, true, false, true, false};
		for (int n=0;n<4;n++)
		{
			for (int m=0;m<9;m++)
			{
				if (stencil_index[m] == interp_index[n] && m != 4)
				{
					stencil[m] -= CInterp[n]; //flag should this be minus?
				}
			}
		}
		//add ns to sparse matrix
		for (int m = 0; m<9; m++)
		{
			if (stencil_used[m])
			{
				row[numE] = iv;
				col[numE] = stencil_index[m];
				val[numE] = stencil[m];
				numE++;
			}
		}
		ns_rhs[iv] = (1-alpha[iv])/Cns[4];
		interp_rhs[iv] = 0;
		//calc new numE
		numE = ny*(nx-1)*5 - ny*2 - (nx-1)*2    +   nx*(ny-1)*5 - nx*2 - (ny-1)*2 + count[iv]-1;
		//add interp corner to sparse matrix
		for (int n=0;n<4;n++)
		{
			for (int m=0;m<9;m++)
			{
				if (stencil_index[m] == interp_index[n] && !stencil_used[m])
				{
					row[numE] = iv;
					col[numE] = interp_index[n];
					val[numE] = -CInterp[n];
				}
				//else if(stencil_index[m] == interp_index[n] && stencil_used[m])
				else if(stencil_index[m] == interp_index[n] && interp_index[n] == iv)
					interp_rhs[iv] += CInterp[n]*q[n];
			}
		}
	}
	else if (ghostTagsUV[iv]>0)
	{
		int interp_index[4] = {index1[iv], index2[iv], index3[iv], index4[iv]};
		bool interp_in[4] = {false, false, false, false};
		int ns_index[5] = {iv+nx, iv+1, iv-nx, iv-1, iv}; //n e s w p
		bool ns_overlap[5] = {false, false, false, false, true};
		double q[4] = {q1[iv], q2[iv], q3[iv], q4[iv]};
		double CInterp[4];
		CInterp[0] = q1coef[iv];
		CInterp[1] = q2coef[iv];
		CInterp[2] = q3coef[iv];
		CInterp[3] = q4coef[iv];
		//count the number of nodes the interp is using
		//find how which ns nodes are occupied
		int counter = 0;
		temp = 0;
		for (int l=0; l<4; l++)
		{
			if (ghostTagsUV[interp_index[l]]>0)
			{
				counter +=1;
				interp_in[l] = true;
			}
			for (int n=0; n<5; n++)
			{
				if (interp_index[l] == ns_index[n])
					ns_overlap[n] = true;
			}
		}
		//add center to matrix
		row[numE] = iv;
		col[numE] = iv;
		val[numE] = 1;
		numE++;
		//add real interp values to matrix
		for (int i=0; i<4; i++)
		{
			if (!interp_in[i] && interp_index[i] != iv)
			{
				row[numE] = iv;
				col[numE] = interp_index[i];
				val[numE] = CInterp[i];
				numE++;
			}
			else
			{
				temp -= CInterp[i] * q[i];
			}
		}
		//fill remainder of values
		int counter2 = 0;
		for (int i=0; i<5; i++)
		{
			if (counter2>=counter)
				break;
			if (ns_overlap[i]==false)
			{
				row[numE] = iv;
				col[numE] = ns_index[i];
				val[numE] = 0;
				numE++;
				counter2++;
			}
		}
		ns_rhs[iv] = 0;
		interp_rhs[iv] = 2*vB[0] + temp;//flag this doesn't account for the interpolation part
	}
	else
	{
	temp = 1 + 0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I+1])*0.5)) + 0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I-1])*0.5)) + 0.5*dt*nu*(1/(dy[J+1]*(dy[J]+dy[J+1])*0.5)) + 0.5*dt*nu*(1/(dy[J]*(dy[J]+dy[J+1])*0.5));

	//EAST
	row[numE] = iv;
	col[numE] = iv+1;
	val[numE] = -0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I+1])*0.5))/temp;
	numE++;

	//WEST
	row[numE] = iv;
	col[numE] = iv-1;
	val[numE] = -0.5*dt*nu*(1/(dx[I]*(dx[I]+dx[I-1])*0.5))/temp;
	numE++;

	//NORTH
	row[numE] = iv;
	col[numE] = iv + nx;
	val[numE] = -0.5*dt*nu*(1/(dy[J+1]*(dy[J]+dy[J+1])*0.5))/temp;
	numE++;

	//SOUTH
	row[numE] = iv;
	col[numE] = iv-nx;
	val[numE] = -0.5*dt*nu*(1/(dy[J]*(dy[J]+dy[J+1])*0.5))/temp;
	numE++;

	//CENTER
	row[numE] = iv;
	col[numE] = iv;
	val[numE] = 1;
	numE++;
	ns_rhs[iv] = 1/temp;
	interp_rhs[iv] = 0;
	}
}

}//end kernel
