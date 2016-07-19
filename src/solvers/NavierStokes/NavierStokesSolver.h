#pragma once

#include <domain.h>
#include <bodies.h>
#include <io/io.h>
#include "newPrecon.h"
#include <parameterDB.h>
#include <preconditioner.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/precond/diagonal.h>
#include <ctime>
#include <cusp/print.h>

/**
 * \class NavierStokesSolver
 * \brief Solves the Navier-Stokes equations in a rectangular domain.
 *
 * Methods are defined as virtual when they are redefined 
 * in a derived class with the same name.
 *
 */
class NavierStokesSolver
{
protected:
	parameterDB *paramDB;		///< database that contains all the simulation parameters
	domain      *domInfo;		///< computational grid information

	cusp::array1d<double, cusp::device_memory>
		u,			///< velocity vector (u_l and u_l+1, depending on where)
		uhat,		///< intermediate velocity vector
		uold,		///< old velocity vector (u_l-1)
		pressure,			///< pressure vector (p_l+1)
		pressure_old,		///< old pressure vector (p_l)
		Nold,		///< convection term for N(uold)
		N,			///< convection term for N(u)
		L,			///< Laplacian of u
		Lnew,		///< Laplacian of u_l+1
		bc1, 		///< when you take L(uhat) you get boundary terms that need to go on the right side of the equation, those are here
		rhs1,		///< -G*p -1.5N(u) + 0.5 N(uold) + 0.5 L(u)
		rhs2,		///< rhs for the intermediate pressure
		bc[4];		///< array that contains the boundary conditions of the rectangular

	size_t
		timeStep,			///< time iteration number
		iterationCount1,	///< number of iteration to solve the intermediate velocities
		iterationCount2;	///< number of iteration to solve the Poisson equation

	double
		forceX,		///< force acting on the body in the x direction //flag
		forceY,		///< force acting on the body in the y direction //flag
		fxx,
		fxy,
		fxu;

	cusp::coo_matrix<int, double, cusp::device_memory>
		LHS1,		///< Matrix for the unknown uhat
		LHS2;		///< Matrix for the unknown phi

	newPrecon PC;

	Logger logger;	///< instance of the class \c Logger to track time of different tasks
	
	std::ofstream iterationsFile;	///< file that contains the number of iterations
	
	//////////////////////////
	//NavierStokesSolver.cu
	//////////////////////////
	void initialiseNoBody();
	void solveIntermediateVelocity();
	void solvePoisson();
	void arrayprint(cusp::array1d<double, cusp::device_memory> value, std::string name, std::string type, int time);
	void printLHS();

	//////////////////////////
	//intermediateVelocity.inl
	//////////////////////////
	void generateN();
	void generateL();
	void generateBC1();

public:
	// constructor -- copy the database and information about the computational grid
	NavierStokesSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);

	//////////////////////////
	//NavierStokesSolver.cu
	//////////////////////////
	virtual void initialise();
	virtual void initialiseLHS();
	virtual void stepTime();
	virtual void writeCommon();
	virtual void writeData();
	bool finished();
	virtual void shutDown();

	//////////////////////////
	//intermediatepressure.inl
	//////////////////////////
	virtual void generateRHS2();
	virtual void generateLHS2();

	//////////////////////////
	//intermediateVelocity.inl
	//////////////////////////
	virtual void generateRHS1();
	virtual void generateLHS1();

	//////////////////////////
	//projectVelocity.inl
	//////////////////////////
	virtual void velocityProjection();

	
	std::string name()
	{
		return "Navier-Stokes";
	}
};
