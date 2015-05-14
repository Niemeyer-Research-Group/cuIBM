/***************************************************************************//**
 * \file TCFSISolver.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the methods of the class \c TairaColoniusSolver.
 */

#include "TCFSISolver.h"
#include <sys/stat.h>

/**
 * \brief Initializes the solvers.
 *
 * Initializes bodies, arrays and matrices of the intermediate flux solver
 * and the Poisson solver.
 *
 */
template <typename memoryType>
void TCFSISolver<memoryType>::initialise()
{	
	//gets and sets number of points in the x and y from dominfo
	int nx = NavierStokesSolver<memoryType>::domInfo->nx,
        ny = NavierStokesSolver<memoryType>::domInfo->ny;

	int numUV = (nx-1)*ny + nx*(ny-1);
	int numP  = nx*ny;

	NavierStokesSolver<memoryType>::initialiseCommon();

	//sizes the arrays to match the domain and sets initial variables in those arrays
	NSWithBody<memoryType>::initialiseBodies();
	int totalPoints  = NSWithBody<memoryType>::B.totalPoints;

	NavierStokesSolver<memoryType>::initialiseArrays(numUV, numP+2*totalPoints);
	if(totalPoints > 0)
	{
		// for each Lagragian point, 24 velocity nodes around it
		// are influenced by the delta function at the point
		// A 4x3 grid of u-nodes, and a 3x4 grid of v-nodes
		// which gives 12*2 velocity nodes influence by a boundary point in 2D
		E.resize(2*totalPoints, numUV, 24*totalPoints);
	}

	//sets up various matricies required to solve, A,G,H,D,E,Q etc
	NavierStokesSolver<memoryType>::assembleMatrices();
	TCFSISolver<memoryType>::initialiseFSI(totalPoints);
}

/**
 * \brief initializes fluid strucutre interaction variables
 */
template <typename memoryType>
void TCFSISolver<memoryType>::initialiseFSI(int bodyPoints)
{
	alpha = 1.5;
	NSWithBody<memoryType>::B.xk = NSWithBody<memoryType>::B.x;
	NSWithBody<memoryType>::B.yk = NSWithBody<memoryType>::B.y;
}

/**
 * \brief Calculates and writes forces acting on each immersed body at current time.
 */
template <typename memoryType>
void TCFSISolver<memoryType>::writeData()
{	
	NavierStokesSolver<memoryType>::logger.startTimer("output");

	NSWithBody<memoryType>::writeCommon();

	parameterDB  &db = *NavierStokesSolver<memoryType>::paramDB;
	real         dt  = db["simulation"]["dt"].get<real>();
	int          numBodies  = NSWithBody<memoryType>::B.numBodies;

	// print to file
	NSWithBody<memoryType>::forceFile << NavierStokesSolver<memoryType>::timeStep*dt << '\t';
	for(int l=0; l<numBodies; l++)
	{
		NSWithBody<memoryType>::forceFile << NSWithBody<memoryType>::B.forceX[l] << '\t' << NSWithBody<memoryType>::B.forceY[l] << '\t';
	}
	NSWithBody<memoryType>::forceFile << std::endl;

	// print forces calculated using the CV approach
	//NSWithBody<memoryType>::calculateForce();
	//NSWithBody<memoryType>::forceFile << NSWithBody<memoryType>::forceX << '\t' << NSWithBody<memoryType>::forceY << std::endl;

	NavierStokesSolver<memoryType>::logger.stopTimer("output");
}

template <typename memoryType>
void TCFSISolver<memoryType>::stepTime()
{
	int k = 0;
	isSub = true;
	NavierStokesSolver<memoryType>::qOld = NavierStokesSolver<memoryType>::q;

	//solve for intermediate velocity
	//not needed in the loop because it doesn't change
	/*
	NavierStokesSolver<memoryType>::generateRN();
	NavierStokesSolver<memoryType>::generateBC1();
	NavierStokesSolver<memoryType>::assembleRHS1();
	NavierStokesSolver<memoryType>::solveIntermediateVelocity();
	*/

	//release at timestep x
	if (NavierStokesSolver<memoryType>::timeStep < 0)
	{
		NavierStokesSolver<memoryType>::generateRN();
		NavierStokesSolver<memoryType>::generateBC1();
		NavierStokesSolver<memoryType>::assembleRHS1();
		NavierStokesSolver<memoryType>::solveIntermediateVelocity();

		NSWithBody<memoryType>::B.calculateCellIndices(*NavierStokesSolver<memoryType>::domInfo);
		NSWithBody<memoryType>::B.calculateBoundingBoxes(*NavierStokesSolver<memoryType>::paramDB,  *NavierStokesSolver<memoryType>::domInfo);
		//solve poisson equation
		TCFSISolver<memoryType>::updateQT(isSub); //makes new QT based on substep uB,vB	
		NavierStokesSolver<memoryType>::generateC(); //QT*BN*Q, QT and Q made in update QT, BN doesn't need to be remade.
		NavierStokesSolver<memoryType>::PC2->update(NavierStokesSolver<memoryType>::C); //update preconditioner
		TCFSISolver<memoryType>::generateBC2(); //new bc2 based on ubk to enforce no slip
		NavierStokesSolver<memoryType>::assembleRHS2(); //QT*qStar + bc2
		TCFSISolver<memoryType>::solvePoisson();

		////projection step
		TCFSISolver<memoryType>::projectionStep();

		////compute the forces exerted on the boundary by the fluid(step 2)
		TCFSISolver<memoryType>::calculateForce();
	}
	else
	{
		//std::cout<<NavierStokesSolver<memoryType>::timeStep*0.01<<"\n k\ty\t\tv\t\tforce\n";
		do
		{
			////setup
			TCFSISolver<memoryType>::calculateCellIndices(*NavierStokesSolver<memoryType>::domInfo);
			TCFSISolver<memoryType>::calculateBoundingBoxes(*NavierStokesSolver<memoryType>::paramDB,  *NavierStokesSolver<memoryType>::domInfo);
			TCFSISolver<memoryType>::updateQT(isSub); //makes new QT based on substep uB,vB	
			NavierStokesSolver<memoryType>::generateC(); //QT*BN*Q, QT and Q made in update QT, BN doesn't need to be remade.
			NavierStokesSolver<memoryType>::PC2->update(NavierStokesSolver<memoryType>::C); //update preconditioner

			////intermediate velocity
			NavierStokesSolver<memoryType>::generateRN();
			NavierStokesSolver<memoryType>::generateBC1();
			NavierStokesSolver<memoryType>::assembleRHS1();
			NavierStokesSolver<memoryType>::solveIntermediateVelocity();

			////poisson equation
			TCFSISolver<memoryType>::generateBC2(); //new bc2 based on ubk to enforce no slip
			NavierStokesSolver<memoryType>::assembleRHS2(); //QT*qStar + bc2
			TCFSISolver<memoryType>::solvePoisson();
			
			////projection step
			TCFSISolver<memoryType>::projectionStep();
			
			////compute the forces exerted on the boundary by the fluid(step 2)
			TCFSISolver<memoryType>::calculateForce();
			
			////solve the structure equation for Xk+1 and uBk+1(step 3) and update structure location (step 4)
			TCFSISolver<memoryType>::solveStructure();
			k++;
			TCFSISolver<memoryType>::printFSI();
		//}while(TCFSISolver<memoryType>::checkConvergence(k) == false);
		}while(NSWithBody<memoryType>::B.converged[0] == false);
	}
	//update iterating variables
	NavierStokesSolver<memoryType>::logger.startTimer("TRANSFER");
	NavierStokesSolver<memoryType>::q = NavierStokesSolver<memoryType>::qk;
	NSWithBody<memoryType>::B.vB = NSWithBody<memoryType>::B.vBk;
	//NSWithBody<memoryType>::B.x=NSWithBody<memoryType>::B.xk;
	NSWithBody<memoryType>::B.y=NSWithBody<memoryType>::B.yk;
	NavierStokesSolver<memoryType>::logger.stopTimer("TRANSFER");
	//NSWithBody<memoryType>::B.uB = NSWithBody<memoryType>::B.uBk;

	NavierStokesSolver<memoryType>::timeStep++;
}

/**
 * \brief Updates the location of the bodies and re-generates appropriate matrices.
 *
 * Updates only done if the bodies are moving.
 * It re-generates the interpolation matrix, therefore the Poisson matrix too.
 *
 */
//not sure if this is needed
template <typename memoryType>
void TCFSISolver<memoryType>::updateSolverState()
{
	if (NSWithBody<memoryType>::B.bodiesMove)
	{
		NSWithBody<memoryType>::updateBodies();
		isSub = false;
		updateQT(isSub);
		NavierStokesSolver<memoryType>::generateC();
		
		NavierStokesSolver<memoryType>::logger.startTimer("preconditioner2");
		NavierStokesSolver<memoryType>::PC2->update(NavierStokesSolver<memoryType>::C);
		NavierStokesSolver<memoryType>::logger.stopTimer("preconditioner2");
	}
}

/**
 * \ brief
 *  checks if the substep loop is complete
 */
template <typename memoryType>
void TCFSISolver<memoryType>::printFSI()
{
	std::cout<<NavierStokesSolver<memoryType>::timeStep*0.01<<"\n";
	//std::cout<<NSWithBody<memoryType>::B.yk[0]<<"   \t"<<NSWithBody<memoryType>::B.vBk[0]<<"   \t"<<NSWithBody<memoryType>::B.forceY[0]<<"\n";
	NSWithBody<memoryType>::output<<NavierStokesSolver<memoryType>::timeStep*0.01<<"  \t";
	NSWithBody<memoryType>::output<<NSWithBody<memoryType>::B.ykp1[0]<<"  \t";
	NSWithBody<memoryType>::output<<NSWithBody<memoryType>::B.yk[0]<<"  \t";
	NSWithBody<memoryType>::output<<NSWithBody<memoryType>::B.forceY[0]<<"  \t";
	NSWithBody<memoryType>::output<<NSWithBody<memoryType>::B.vBk[0]<<"\n";
}

/**
 * \brief Constructor. Copies the database and information about the computational grid.
 */
template <typename memoryType>
TCFSISolver<memoryType>::TCFSISolver(parameterDB *pDB, domain *dInfo)
{
	NavierStokesSolver<memoryType>::paramDB = pDB;
	NavierStokesSolver<memoryType>::domInfo = dInfo;
}

// include inline files located in "./TairaColonius/"
#include "TCFSI/generateQT.inl"
#include "TCFSI/generateBC2.inl"
#include "TCFSI/calculateForce.inl"
#include "TCFSI/solvePoisson.inl"
#include "TCFSI/projectionStep.inl"
#include "TCFSI/solveStructure.inl"
#include "TCFSI/bounding.inl"

// specialization of the class TCFSISolver
template class TCFSISolver<host_memory>;
template class TCFSISolver<device_memory>;
