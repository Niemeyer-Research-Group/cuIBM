/***************************************************************************//**
 * \file projectVelocity.h
 * \author Chris Minar
 * \brief Declaration of kernels to calculate step 3: project u
 */

#pragma once

/**
 * \namespace kernels
 * \brief Contains all the custom-written CUDA kernels.
 */

namespace kernels
{
__global__
void project_velocity_X(double *u, double *uhat, double *uold, double *pressure, int *tagsP, int *tagsIn, double *dx, double dt, int nx, int ny);

__global__
void project_velocity_Y(double *u, double *uhat, double *uold, double *pressure, int *tagsP, int *tagsIn, double *dy, double dt, int nx, int ny);
}
