#pragma once
#include <preconditioner.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/precond/diagonal.h>
#include <cusp/linear_operator.h>

class newPrecon
{

	public:
	//cusp::precond::aggregation::smoothed_aggregation<int, double, cusp::device_memory> PC1;
	//cusp::linear_operator<double, cusp::device_memory, int>* PC1;
	//cusp::precond::aggregation::smoothed_aggregation<int, double, cusp::device_memory> PC2;
	preconditioner< cusp::coo_matrix<int, double, cusp::device_memory> >
		*PC1,		///< preconditioner for the intermediate flux solver
		*PC2;		///< preconditioner for the Poisson solver

	newPrecon();
	void generate1(cusp::coo_matrix<int, double, cusp::device_memory>LHS1, preconditionerType type1);

	void generate2(cusp::coo_matrix<int, double, cusp::device_memory>LHS2, preconditionerType type2);

	void generate(cusp::coo_matrix<int, double, cusp::device_memory>LHS1, cusp::coo_matrix<int, double, cusp::device_memory>LHS2, preconditionerType type1, preconditionerType type2);

	void update1(cusp::coo_matrix<int, double, cusp::device_memory>LHS1);

	void update2(cusp::coo_matrix<int, double, cusp::device_memory>LHS2);

};
