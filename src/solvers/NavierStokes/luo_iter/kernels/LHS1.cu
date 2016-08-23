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
					double *xu, double *yu, double *detA, double *alpha,
					double *b11, double *b12, double *b13, double *b14,
					double *b21, double *b22, double *b23, double *b24,
					double *b31, double *b32, double *b33, double *b34,
					double *b41, double *b42, double *b43, double *b44,
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
		int ns_index[5] = {iu + (nx-1), iu + 1, iu - (nx-1), iu -1, iu}; //n e s w p
		double q[4] = {q1[iu], q2[iu], q3[iu], q4[iu]};
		double CInterp[4];
		double Cns[5];
		double	x=xu[I],
				y=yu[J];
		Cns[0] = -dt*nu/(dy[J+1]*(dy[J]+dy[J+1]));
		Cns[1] = -dt*nu/(dx[I]  *(dx[I]+dx[I+1]));
		Cns[2] = -dt*nu/(dy[J]  *(dy[J]+dy[J+1]));
		Cns[3] = -dt*nu/(dx[I]  *(dx[I]+dx[I-1]));
		Cns[4] = 1-Cns[0] - Cns[1] - Cns[2] - Cns[3];
		CInterp[0] = (b11[iu] + b21[iu]*x + b31[iu]*y + b41[iu]*x*y)/detA[iu];
		CInterp[1] = (b12[iu] + b22[iu]*x + b32[iu]*y + b42[iu]*x*y)/detA[iu];
		CInterp[2] = (b13[iu] + b23[iu]*x + b33[iu]*y + b43[iu]*x*y)/detA[iu];
		CInterp[3] = (b14[iu] + b24[iu]*x + b34[iu]*y + b44[iu]*x*y)/detA[iu];
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
					stencil[m] += CInterp[n]; //flag should this be minus?
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
					val[numE] = CInterp[n];
				}
				else if(stencil_index[m] == interp_index[n] && stencil_used[m])
					interp_rhs[iu] += CInterp[n]*q[n];
			}
		}
	}
	/*else if (ghostTagsUV[iu]>0)
	{

	}*/
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
						double *xv, double *yv, double *detA, double *alpha,
						double *b11, double *b12, double *b13, double *b14,
						double *b21, double *b22, double *b23, double *b24,
						double *b31, double *b32, double *b33, double *b34,
						double *b41, double *b42, double *b43, double *b44,
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
		int ns_index[5] = {iv + nx, iv + 1, iv - nx, iv -1, iv}; //n e s w p
		double q[4] = {q1[iv], q2[iv], q3[iv], q4[iv]};
		double CInterp[4];
		double Cns[5];
		double	x=xv[I],
				y=yv[J];
		Cns[0] = -dt*nu/(dy[J+1]*(dy[J]+dy[J+1]));
		Cns[1] = -dt*nu/(dx[I]*(dx[I]+dx[I+1]));
		Cns[2] = -dt*nu/(dy[J]*(dy[J]+dy[J+1]));
		Cns[3] = -dt*nu/(dx[I]*(dx[I]+dx[I-1]));
		Cns[4] = 1-Cns[0] - Cns[1] - Cns[2] - Cns[3];
		CInterp[0] = (b11[iv] + b21[iv]*x + b31[iv]*y + b41[iv]*x*y)/detA[iv];
		CInterp[1] = (b12[iv] + b22[iv]*x + b32[iv]*y + b42[iv]*x*y)/detA[iv];
		CInterp[2] = (b13[iv] + b23[iv]*x + b33[iv]*y + b43[iv]*x*y)/detA[iv];
		CInterp[3] = (b14[iv] + b24[iv]*x + b34[iv]*y + b44[iv]*x*y)/detA[iv];
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
					stencil[m] += CInterp[n]; //flag should this be minus?
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
					val[numE] = CInterp[n];
				}
				else if(stencil_index[m] == interp_index[n] && stencil_used[m])
					interp_rhs[iv] += CInterp[n]*q[n];
			}
		}
	}
	/*else if (ghostTagsUV[iv]>0)
	{

	}*/
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
