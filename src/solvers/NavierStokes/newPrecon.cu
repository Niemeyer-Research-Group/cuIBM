/***************************************************************************//**
 * \file newPrecon.cu
 * \author Christopher Minar (minarc@oregonstate.edu)
 * \Compiling the preconditioner class takes a very very long time so it has been moved into its own class so it doesn't get flagged for recompiling every time something in the main files is changed.
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

void newPrecon::generate1(cusp::coo_matrix<int, double, cusp::device_memory>LHS1, preconditionerType type1)
{
	PC1 = new preconditioner< cusp::coo_matrix<int, double, cusp::device_memory> >(LHS1, type1);
}

void newPrecon::generate2(cusp::coo_matrix<int, double, cusp::device_memory>LHS2, preconditionerType type2)
{
	PC2 = new preconditioner< cusp::coo_matrix<int, double, cusp::device_memory> >(LHS2, type2);
}

void newPrecon::update1(cusp::coo_matrix<int, double, cusp::device_memory>LHS1)
{
	PC1->update(LHS1);
}

void newPrecon::update2(cusp::coo_matrix<int, double, cusp::device_memory>LHS2)
{
	PC2->update(LHS2);
}

void newPrecon::generate(cusp::coo_matrix<int, double, cusp::device_memory>LHS1, cusp::coo_matrix<int, double, cusp::device_memory>LHS2, preconditionerType type1, preconditionerType type2)
{
	PC1 = new preconditioner< cusp::coo_matrix<int, double, cusp::device_memory> >(LHS1, type1);
	PC2 = new preconditioner< cusp::coo_matrix<int, double, cusp::device_memory> >(LHS2, type2);
}
