/***************************************************************************//**
 * \file calculateForce.cu
 * \author Anush Krishnan (anush@bu.edu),
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \based of original cuIBM
 */


#include "calculateForce.h"

namespace kernels
{

__global__
void calcForceFadlun(double *force, double *L, double *Lnew, double *Nold, double *N, double *u, double *uold, int *tags, double dt, int nx, int ny)
{
	if (threadIdx.x + (blockDim.x * blockIdx.x) >= (ny-1)*nx + (nx-1)*ny)
			return;
	int i 	= threadIdx.x + (blockDim.x * blockIdx.x);

	force[i] = (tags[i]!=-1)  *  (
								   u[i]/dt - 0.5*Lnew[i] - 				//u^l+1
								   uold[i]/dt + 1.5*N[i] - 0.5*L[i] - 	//u^l
								   0.5*Nold[i]	 						//u^l-1
								 );
}

/**
 * \brief Calculates drag using a control-volume approach (left-right).
 *
 * Evaluate the contribution from the left and right parts of the control surface.
 *
 * \param FxX raw pointer to the vector storing the drag in the x-direction
 * \param lambda raw pointer to the vector storing all the pressure and Lagrangian forces
 * \param q raw pointer to the vector storing all the fluxes
 * \param nu viscosity
 * \param dx raw pointer to the vector storing the cell widths in the x-direction
 * \param dy raw pointer to the vector storing the cell widths in the y-direction
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direction
 * \param I x-index of the bottom-left corner cell of the control surface
 * \param J y-index of the top-right corner cell of the control surface
 * \param ncx number of cells in the x-direction in the control volume
 * \param ncy number of cells in the y-direction in the control volume
 */
__global__
void dragLeftRight(double *FxX, double *u, double *p, double nu, double *dx, double *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx >= ncy)
		return;
	int  Ip = (J+idx)*nx + I,
	     Iu = (J+idx)*(nx-1) + I;

	FxX[idx] = -(
	              // multiply the pressure with the surface area to get p dy
	              //(p[e]-p[w])*dy
				  (
				   p[Ip+ncx] - p[Ip]
				  )*dy[J+idx]
	              +
	              // ur^2 - ul^2 * dy
	              (
	                  (u[Iu+ncx]+u[Iu+ncx-1])*(u[Iu+ncx]+u[Iu+ncx-1])/4
	                - (u[Iu-1]+u[Iu])*(u[Iu-1]+u[Iu])/4
	              )*dy[J+idx]
	              -
	              // du/dx * dy
	              // approximate using dudx of the inside cell of the cv instead of the lr average
	              nu*
	              (
	                  (u[Iu+ncx] - u[Iu+ncx-1])/dx[I+ncx]
	                - (u[Iu] - u[Iu-1])/dx[I]
	              )*dy[J+idx]
	            );
}

/**
 * \brief Calculate drag using a control-volume approach (bottom-top).
 *
 * Evaluate the contribution from the bottom and top parts of the control surface.
 *
 * \param FxY raw pointer to the vector storing the drag in the y-direction
 * \param q raw pointer to the vector storing all the fluxes
 * \param nu viscosity
 * \param dx raw pointer to the vector storing the cell widths in the x-direction
 * \param dy raw pointer to the vector storing the cell widths in the y-direction
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direction
 * \param I x-index of the bottom-left corner cell of the control surface
 * \param J y-index of the top-right corner cell of the control surface
 * \param ncx number of cells in the x-direction in the control volume
 * \param ncy number of cells in the y-direction in the control volume
 */
__global__
void dragBottomTop(double *FxY, double *u, double nu, double *dx, double *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx > ncx)
		return;
	int  Iu = J*(nx-1) + (I-1+idx),
	     Iv = (nx-1)*ny + (J-1)*nx + I+idx;
	FxY[idx] = -(
	              // multiply by dS
	              (
	                0.25 * ( u[Iu+ncy*(nx-1)] + u[Iu+(ncy-1)*(nx-1)] )
	                     * ( u[Iv+ncy*nx] + u[Iv+ncy*nx-1] )
	                -
	                0.25 * ( u[Iu] + u[Iu-(nx-1)] )
	                     * ( u[Iv] + u[Iv-1] )
	              )
	              -
	              // multiply by dS (cannot use the leftRight trick in this case)
	              nu*
	              (
	                (
	                  (u[Iu+ncy*(nx-1)] - u[Iu+(ncy-1)*(nx-1)])/2.0/(dy[J+ncy]+dy[J+ncy-1]) +
	                  (u[Iv+ncy*nx] - u[Iv+ncy*nx-1])          /2.0/(dx[I+idx]+dx[I+idx-1])
	                )
	                -
	                (
	                  (u[Iu] - u[Iu-(nx-1)])/2.0/(dy[J]+dy[J-1]) +
	                  (u[Iv] - u[Iv-1])     /2.0/(dx[I+idx]+dx[I+idx-1])
	                )
	              )
	            )*0.5*(dx[I+idx]+dx[I+idx-1]);

}

/**
 * \brief Calculate drag using a control-volume approach (unsteady).
 *
 * Evaluate the unsteady contribution of the control volume.
 *
 * \param FxU raw pointer to the vector storing the unsteady drag components
 * \param q raw pointer to the vector storing all the fluxes
 * \param qOld raw pointer to the vector sotring all the fluxes at the previous time-step
 * \param dx raw pointer to the vector storing the cell widths in the x-direction
 * \param dy raw pointer to the vector storing the cell widths in the y-direction
 * \param dt time increment
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direcyion
 * \param I x-index of the bottom-left cell of the control surface
 * \param J y-index of the top-right cell of the control surface
 * \param ncx number of cells in the x-direction in the control volume
 * \param nyc number of cells in the y-direction in the control volume
 */
__global__
void dragUnsteady(double *FxU, double *u, double *uold, int *tagsIn, double *dx, double *dy, double dt,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;

	if(idx >= (ncx+1)*ncy)
		return;

	int i = idx%(ncx+1),
	    j = idx/(ncx+1);

	int Iu = (J+j)*(nx-1) + (I-1+i);

	FxU[idx] = - (tagsIn[Iu] == -1) * ((u[Iu]*dy[J+j] - uold[Iu]*dy[J+j])/dt * 0.5*(dx[I+i]+dx[I-1+i]));
}

/**
 * \brief Calculate lift using a control-volume approach (left-right).
 *
 * Evaluate the contribution from the left and right parts of the control surface.
 *
 * \param FyX raw pointer to the vector storing the lift components in the x-direction
 * \param q raw pointer to the vector storing all the fluxes
 * \param nu viscosity
 * \param dx raw pointer to the vector storing the cell widths in the x-direction
 * \param dy raw pointer to the vector storing the cell widths in the y-direction
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direcyion
 * \param I x-index of the bottom-left cell of the control surface
 * \param J y-index of the top-right cell of the control surface
 * \param ncx number of cells in the x-direction in the control volume
 * \param nyc number of cells in the y-direction in the control volume
 */
__global__
void liftLeftRight(double *FyX, double *u, double nu, double *dx, double *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx > ncy)
		return;
	int  Iu = (J+idx)*(nx-1) + (I-1),
	     Iv = (nx-1)*ny + (J-1+idx)*nx + I;
	FyX[idx] = -(
	              // multiply by dS
	              (
	                0.25 * ( u[Iu+ncx] + u[Iu+ncx-(nx-1)] )
	                     * ( u[Iv+ncx] + u[Iv+ncx-1] )
	                -
	                0.25 * ( u[Iu] + u[Iu-(nx-1)] )
	                     * ( u[Iv] + u[Iv-1] )
	              )
	              -
	              // multiply by dS (cannot use the leftRight trick in this case)
	              nu*
	              (
	                (
	                  (u[Iu+ncx] - u[Iu+ncx-(nx-1)])/2.0/(dy[J+idx]+dy[J-1+idx]) +
	                  (u[Iv+ncx] - u[Iv+ncx-1])/2.0/(dx[I+ncx]+dx[I+ncx-1])
	                )
	                -
	                (
	                  (u[Iu] - u[Iu-(nx-1)])/2.0/(dy[J+idx]+dy[J-1+idx]) +
	                  (u[Iv] - u[Iv-1])/2.0/(dx[I]+dx[I-1])
	                )
	              )
	            )*0.5*(dy[J+idx]+dy[J-1+idx]);
}

/**
 * \brief Calculate lift using a control-volume approach (bottom-top).
 *
 * Evaluate the contribution from the bottom and top parts of the control surface.
 *
 * \param FyY raw pointer to the vector storing the lift components in the y-direction
 * \param q raw pointer to the vector storing all the fluxes
 * \param lambda raw pointer to the vector storing the pressure and Lagrangian forces
 * \param nu viscosity
 * \param dx raw pointer to the vector storing the cell widths in the x-direction
 * \param dy raw pointer to the vector storing the cell widths in the y-direction
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direcyion
 * \param I x-index of the bottom-left cell of the control surface
 * \param J y-index of the top-right cell of the control surface
 * \param ncx number of cells in the x-direction in the control volume
 * \param nyc number of cells in the y-direction in the control volume
 */
__global__
void liftBottomTop(double *FyY, double *u, double *p, double nu, double *dx, double *dy,
                   int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;
	if(idx >= ncx)
		return;
	int  Ip = J*nx + I+idx,
	     Iv = (nx-1)*ny + (J-1)*nx + I+idx;
	FyY[idx] = -(
	              // multiply the pressure with the surface area to get p dx
	              (p[Ip+ncy*nx]-p[Ip-nx])*dx[I+idx]
	              +
	              // divide q^2 by dx, so that just v^2 dx is obtained
	              (
	                  0.25*(u[Iv+(ncy+1)*nx] + u[Iv+ncy*nx])*(u[Iv+(ncy+1)*nx] + u[Iv+ncy*nx])
	                - 0.25*(u[Iv] + u[Iv-nx])*(u[Iv] + u[Iv-nx])
	              )*dx[I+idx]
	              -
	              // no multiplication or division since dv/dy dx = dq/dy
	              nu*
	              (
	                  (u[Iv+(ncy+1)*nx] - u[Iv+ncy*nx])*dx[I+idx]/dy[J+ncy]
	                - (u[Iv] - u[Iv-nx])*dx[I]/dy[J-1]
	              )
	            );
}

/**
 * \brief Calculate lift using a control-volume approach (unsteady).
 *
 * Evaluate the unsteady contribution of the control volume.
 *
 * \param FyU raw pointer to the vector storing the unsteady lift components
 * \param q raw pointer to the vector storing all the fluxes
 * \param qOld raw pointer to the vector sotring all the fluxes at the previous time-step
 * \param dx raw pointer to the vector storing the cell widths in the x-direction
 * \param dy raw pointer to the vector storing the cell widths in the y-direction
 * \param dt time increment
 * \param nx number of cells in the x-direction
 * \param ny number of cells in the y-direcyion
 * \param I x-index of the bottom-left cell of the control surface
 * \param J y-index of the top-right cell of the control surface
 * \param ncx number of cells in the x-direction in the control volume
 * \param nyc number of cells in the y-direction in the control volume
 */
__global__
void liftUnsteady(double *FyU, double *u, double *uold, double *dx, double *dy, double dt,
                  int nx, int ny, int I, int J, int ncx, int ncy)
{
	int  idx = threadIdx.x + blockIdx.x*blockDim.x;

	if( idx >= ncx*(ncy+1) )
		return;

	int i = idx%ncx,
	    j = idx/ncx;

	int Iv = (J-1+j)*nx + (I+i) + (nx-1)*ny;

	FyU[idx] = - (u[Iv]*dx[I+i] - uold[Iv]*dx[I+i])/dt * 0.5*(dy[J+j]+dy[J-1+j]);
}

/**
* \brief To be documented
*/
}

