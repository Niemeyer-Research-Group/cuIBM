/*****************************************************************************
 * \file  calculateForce.h
 * \author Christopher Minar
 * \brief Declaration of the kernels to solve the structure equation
 *        Curvilinear immersed boundary method for simulating fluid structure interaction with complex 3D rigid bodies
		Iman Borazjani, Liang Geb, Fotis Sotiropoulos
 */


//not sure if needed
#pragma once

#include <types.h>

namespace kernels
{
__global__ \
void solveStructureHost();
}
