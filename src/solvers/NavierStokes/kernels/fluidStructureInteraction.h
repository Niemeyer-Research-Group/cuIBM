#pragma once
#include <types.h>

namespace kernels
{
__global__
void vorticityInducedVibrationsSC(real *vBk, real *vB, real *y, real *ykp1, real forcey, real Mred, real Ured, real Cy, real dt, real alpha_);

__global__
void vorticityInducedVibrationsLC();

__global__
void freeXYSC();

__global__
void freeXYLC();

__global__
void checkConvergencePosition(real tol, bool *flag, real *yk, real *ykp1);

}
