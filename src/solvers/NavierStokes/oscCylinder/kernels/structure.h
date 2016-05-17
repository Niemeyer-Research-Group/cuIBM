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
void update_body_viv(double *y, double *vB, double dy, double vnew, int totalPoints);
}
