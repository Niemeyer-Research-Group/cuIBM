/***************************************************************************//**
 * \file structure.h
 * \author Chris Minar
 */

#pragma once

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void update_body_viv(double *x, double *uB, double dx, double unew, int totalPoints);
__global__
void initialise_old(double *uB0, double vnew, int totalPoints);
}
