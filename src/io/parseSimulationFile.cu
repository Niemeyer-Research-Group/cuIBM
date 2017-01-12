/***************************************************************************//**
 * \file parseSimulationFile.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Parses the file \a simParams.yaml and stores the numerical 
 *        parameters used in the simulation.
 */


#include <fstream>
#include <yaml-cpp/yaml.h>
#include "io.h"
#include <parameterDB.h>
#include <preconditioner.h>
#include <types.h>


/**
 * \namespace io
 * \brief Contains functions related to I/O tasks.
 */
namespace io
{

using std::string;

/**
 * \brief Converts a string to a prconditioner type.
 *
 * \param s the string that describes the preconditioner
 *
 * \return a preconditioner type
 */
preconditionerType preconditionerTypeFromString(string &s)
{
  if (s == "NONE")
    return NONE;
  else if (s == "DIAGONAL")
    return DIAGONAL;
  else if (s == "SMOOTHED_AGGREGATION")
    return SMOOTHED_AGGREGATION;
  else if (s == "AINV")
    return AINV;
  else
    return NONE;
}

solverType solverTypeFromString(string &s)
{
	if (s =="NAVIER_STOKES_SOLVER")
		return NAVIERSTOKES;
	else if (s == "FADLUN")
		return FADLUN;
	else if (s == "OSC_CYLINDER")
		return OSC;
	else if (s == "FSI")
		return FSI;
	else if (s == "LUO")
		return LUO;
	else
		return NAVIERSTOKES;
}

/**
 * \brief Fills the database with the simulation parameters.
 *
 * \param node the parsed file
 * \param DB database that contains the simulation parameters
 */
void parseSimulation(const YAML::Node &node, parameterDB &DB)
{
	double   dt = 0.02,
	       scaleCV = 2.0;
	int    nt = 100,
	       nsave = 100,
	       startStep = 0,
	       VIV = 0,
	       Ured=3;
	string convSch = "ADAMS_BASHFORTH_2";
	bool   restart = false;

	string SolverType = "NAVIER_STOKES_SOLVER";
	// read simulation parameters
	node["dt"] >> dt;
	node["nsave"] >> nsave;
	node["nt"] >> nt;
	try
	{
		node["SolverType"] >> SolverType;
	}
	catch(...)
	{
	}
	try
	{
		node["restart"] >> restart;
	}
	catch(...)
	{
	}
	try
	{
		node["startStep"] >> startStep;
	}
	catch(...)
	{
	}
	try
	{
		node["scaleCV"] >> scaleCV;
	}
	catch(...)
	{
	}
	try
	{
		node["VIV"] >> VIV;
	}
	catch(...)
	{
	}
	try
	{
		node["Ured"] >> Ured;
	}
	catch(...)
	{
	}

	// write to DB
	string dbKey = "simulation";
	DB[dbKey]["dt"].set<double>(dt);
	DB[dbKey]["scaleCV"].set<double>(scaleCV);
	DB[dbKey]["nsave"].set<int>(nsave);
	DB[dbKey]["nt"].set<int>(nt);
	DB[dbKey]["restart"].set<bool>(restart);
	DB[dbKey]["SolverType"].set<solverType>(solverTypeFromString(SolverType));
	DB[dbKey]["VIV"].set<int>(VIV);
	DB[dbKey]["Ured"].set<int>(Ured);

	string system = "velocity", linearSolver = "CG", preconditioner = "DIAGONAL";
	double tol = 1e-5;
	int maxIter = 10000;

	const YAML::Node &solvers = node["linearSolvers"];
	for (unsigned int i=0; i<solvers.size(); i++)
	{
		// read linear solver options
		solvers[i]["system"] >> system;
		solvers[i]["solver"] >> linearSolver;
		try
		{
			solvers[i]["preconditioner"] >> preconditioner;
		}
		catch(...)
		{
		}
		try
		{
			solvers[i]["tolerance"] >> tol;
		}
		catch(...)
		{
		}
		try
		{
			solvers[i]["maxIterations"] >> maxIter;
		}
		catch(...)
		{
		}

		// write to DB
		string dbKey = system + "Solve";
		DB[dbKey]["solver"].set<string>(linearSolver);
		DB[dbKey]["preconditioner"].set<preconditionerType>(preconditionerTypeFromString(preconditioner));
		DB[dbKey]["tolerance"].set<double>(tol);
		DB[dbKey]["maxIterations"].set<int>(maxIter);
	}
}

/**
 * \brief Parses \a simParams.yaml and stores the simulation parameters.
 *
 * \param simFile the file that contains the simulation parameters
 * \param DB the database that will be filled
 */
void parseSimulationFile(std::string &simFile, parameterDB &DB)
{
	std::ifstream fin(simFile.c_str());
	YAML::Parser parser(fin);
	YAML::Node doc;
	parser.GetNextDocument(doc);

	for(unsigned int i=0; i<doc.size(); i++)
		parseSimulation(doc[i], DB);
}

} // end namespace io
