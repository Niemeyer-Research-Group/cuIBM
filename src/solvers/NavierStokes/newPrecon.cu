/***************************************************************************//**
 * \file newPrecon.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief a temporary way to move the preconditioner out of the main file so it doesn't take forever to compile
 */
#include "newPrecon.h"
#include <cusp/precond/aggregation/smoothed_aggregation.h>//flag
#include <cusp/linear_operator.h>
/*
 * Initialises the new precond class
 * This class exists to seperate the preconditioners from the main code because the preconditioners take forever to compile (~5min)
 */
newPrecon::newPrecon()
{
}

/*
 * Creates new preconditioners for the velocity and poission solve
 * param coo_matrix LHS1 the left hand side matrix for the velocity solve
 * param coo_matrix LHS1 the left hand side matrix for the poisson solve
 * param type1 the type of preconditioner desired for the velocity solve
 * param type2 the type of preconditioner desired for the poisson solve
 */
void newPrecon::generate(cusp::coo_matrix<int, double, cusp::device_memory>LHS1, cusp::coo_matrix<int, double, cusp::device_memory>LHS2, preconditionerType type1, preconditionerType type2)
{
	PC1 = new preconditioner< cusp::coo_matrix<int, double, cusp::device_memory> >(LHS1, type1);
	PC2 = new preconditioner< cusp::coo_matrix<int, double, cusp::device_memory> >(LHS2, type2);
}

/*
 * Updates the precondtioners with new left hand side matricies
 * param coo_matrix LHS1 the left hand side matrix for the velocity solve
 * param coo_matrix LHS1 the left hand side matrix for the poisson solve
 */
void newPrecon::update(cusp::coo_matrix<int, double, cusp::device_memory>LHS1, cusp::coo_matrix<int, double, cusp::device_memory>LHS2)
{
	PC1->update(LHS1);
	PC2->update(LHS2);
}
