/***************************************************************************//**
 * \file newPrecon.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \brief a temporary way to move the preconditioner out of the main file so it doesn't take forever to compile
 */
#include "newPrecon.h"
#include <cusp/precond/aggregation/smoothed_aggregation.h>//flag
#include <cusp/linear_operator.h>
newPrecon::newPrecon()
{
}

void newPrecon::generate(cusp::coo_matrix<int, double, cusp::device_memory>LHS1, cusp::coo_matrix<int, double, cusp::device_memory>LHS2, preconditionerType type1, preconditionerType type2)
{
	//PC1 = new cusp::precond::diagonal<double, cusp::device_memory>(LHS1);
	//cusp::precond::aggregation::smoothed_aggregation<int, double, cusp::device_memory> PC2;
	//PC2.sa_initialize(LHS2);
	PC1 = new preconditioner< cusp::coo_matrix<int, double, cusp::device_memory> >(LHS1, type1);
	PC2 = new preconditioner< cusp::coo_matrix<int, double, cusp::device_memory> >(LHS2, type2);
}

void newPrecon::update(cusp::coo_matrix<int, double, cusp::device_memory>LHS1, cusp::coo_matrix<int, double, cusp::device_memory>LHS2)
{
	PC1->update(LHS1);
	PC2->update(LHS2);
}
