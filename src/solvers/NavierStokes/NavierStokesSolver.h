#pragma once

#include <domain.h>
#include <bodies.h>
#include <io/io.h>
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
		force,		///< force at a a node
		bc[4];		///< array that contains the boundary conditions of the rectangular
		
	cusp::array1d<int, cusp::host_memory>
		tags,		///< indices of velocity nodes near the IB on the host
		tags2,		///< indices of 2nd closest velocity node on the host
		tagsIn,		///< indices of velocity nodes inside the IB
		tagsP,		///< indices of pressure nodes inside the IB
		tagsPOut;	///< indices of pressure nodes outside the IB

	cusp::array1d<int, cusp::device_memory>
		tagsD,		///< indices of velocity nodes near the IB on the device
		tags2D,		///< indices of 2nd closest velocity node on the device
		tagsInD,	///< indices of velocity nodes inside the IB
		tagsPD, 	///< indices of pressure nodes inside the IB
		tagsPOutD;	///< indices of pressure nodes outside the IB

	cusp::array1d<double, cusp::host_memory>
		a,			///< distance between IB and tagged node on the host
		b,			///< distance between tags and tags2 on the host
		distance_from_u_to_body,			///< distance between u velocity node and body
		distance_from_v_to_body,			///< distance between v velocity node and body
		uv;			///< velocity at the IB on the host

	cusp::array1d<double, cusp::device_memory>
		aD,			///< distance between IB and tagged node on the device
		bD,			///< distance between tags and tags2 on the device
		distance_from_u_to_bodyD,
		distance_from_v_to_bodyD,
		uvD;		///< velocity at the IB on the device

	size_t
		timeStep,			///< time iteration number
		iterationCount1,	///< number of iteration to solve the intermediate velocities
		iterationCount2;	///< number of iteration to solve the Poisson equation

	double
		forceX,		///< force acting on the body in the x direction
		forceY;		///< force acting on the body in the y direction

	cusp::coo_matrix<int, double, cusp::device_memory>
		LHS1,		///< Matrix for the unknown uhat
		LHS2;		///< Matrix for the unknown phi
	     
	preconditioner< cusp::coo_matrix<int, double, cusp::device_memory> >
	*PC1,		///< preconditioner for the intermediate flux solver
	*PC2;		///< preconditioner for the Poisson solver

	bodies 	B;		///< bodies in the flow

	Logger logger;	///< instance of the class \c Logger to track time of different tasks
	
	std::ofstream iterationsFile;	///< file that contains the number of iterations
	std::ofstream forceFile;
			
	// assemble the right hand-side of the system for the intermediate flux
	void generateRHS1();
	
	// update the boundaries for the convective boundary condition
	void updateRobinBoundary();

	// calculate -1.5 u grad u
	void generateN();

	// calculate discretized laplacian
	void generateL();
	
	// calculate the force based off the fadlun method
	void calculateForceFadlun();

	// calculate the force based off the lia and peskin control volume method
	void calculateForce();

	// fill the left hand side matrix (based off u_hat -dt*0.5*L(u_hat) )
	void generateLHS1();

	//fill the lhs matrix for uhat with no body
	void generateLHS1NoBody();

	// extra terms that come from -dt*0.5*L(u_hat) at the boundaries
	void generateBC1();
	
	// calculate the fake pressure
	void generateRHS2();
	
	//
	void generateLHS2();

	//velocity projection
	void velocityProjection();

	//pressure projection step
	void pressureProjection();

	// assemble the right hand-side of the Poisson system.
	void assembleRHS2();
	
	// solve for the intermediate flux velocity
	void solveIntermediateVelocity();
	
	// solve the Poisson system for the pressure (and the body forces if immersed body)
	void solvePoisson();
	
	// project the flux onto the divergence-free field
	void projectionStep();

	//solves for the new pressure
	void updatePressure();

public:
	// constructor -- copy the database and information about the computational grid
	NavierStokesSolver(parameterDB *pDB=NULL, domain *dInfo=NULL);

	void arrayprint(cusp::array1d<double, cusp::device_memory> value, std::string name, std::string type);

	// tags points on the immersed boundary
	void tagPoints();

	// tags points on the immersed boundary
	void tagPoints(double *bx, double *by, double *uB, double *vB);

	// initialize parameters, arrays and matrices required for the simulation
	void initialise();
	
	// calculate the variables at the next time step
	void stepTime();
	
	// write numerical solution and number of iterations performed in each solver.
	void writeCommon();
	
	// write data into files
	void writeData();
	
	// evaluate the condition required to stop the simulation
	bool finished();
	
	// print timing information and close the different files
	void shutDown();
	/**
	 * \brief Returns the name of the current solver.
	 */
	std::string name()
	{
		return "Navier-Stokes";
	}
};
