/***************************************************************************//**
 * \file intermediatePressure.h
 * \author Chris Minar
 * \brief Declaration of kernels to generate the right hand side of step 2: solve for intermediate pressure
 */

#pragma once

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void intermediatePressure_luo(double *rhs2, double *uhat, double *ym, double *yp, double *xm, double *xp, double *dx, double *dy, int nx, int ny);
}//end namespace kernels
