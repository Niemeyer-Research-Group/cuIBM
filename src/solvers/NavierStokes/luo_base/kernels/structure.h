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
void update_body_viv(double *By, double *vB, double *Bdy, double vnew, double midY, int totalPoints);
__global__
void initialise_old(double *uB0, double vnew, int totalPoints);
}
