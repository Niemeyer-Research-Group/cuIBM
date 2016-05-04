
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
