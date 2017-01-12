/***************************************************************************//**
 * \file parseBodiesFile.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Parses the file \a bodies.yaml to get information about immersed bodies.
 */


#include <fstream>
#include <vector>
#include <yaml-cpp/yaml.h>
#include "io.h"
#include <body.h>


/**
 * \namespace io
 * \brief Contains functions related to I/O tasks.
 */
namespace io
{

using std::string;

/**
 * \brief Overloads the operator >>. Stores information about an immersed body.
 *
 * \param node the parsed file
 * \param instance of the class \c body to be filled
 */
void operator >> (const YAML::Node &node, body &Body)
{
	try
	{
		// velocity
		for (int i=0; i<2; i++)
			node["velocity"][i] >> Body.velocity[i];
	}
	catch(...)
	{
	}
	try
	{
		// oscillation in X
		for (int i=0; i<3; i++)
			node["xOscillation"][i] >> Body.xOscillation[i];
		Body.xOscillation[1] *= 2*M_PI;
	}
	catch(...)
	{
	}
	try
	{
		// oscillation in Y
		for (int i=0; i<3; i++)
			node["yOscillation"][i] >> Body.yOscillation[i];
		Body.yOscillation[1] *= 2*M_PI;
	}
	catch(...)
	{
	}

	try
	{
		node["xCoefficient"] >> Body.xCoefficient;
		node["xPhase"] >> Body.xPhase;
	}
	catch(...)
	{
		std::cout<<"x amplitude/phase not found\n\n";
	}
	try
	{
		node["yCoefficient"] >> Body.yCoefficient;
		node["yPhase"] >> Body.yPhase;
	}
	catch(...)
	{
		std::cout<<"y amplitude/phase not found\n\n";
	}
	try
	{
		node["uCoefficient"] >> Body.uCoefficient;
		node["uPhase"] >> Body.uPhase;
	}
	catch(...)
	{
		std::cout<<"u amplitude/phase not found\n\n";
	}
	try
	{
		node["vCoefficient"] >> Body.vCoefficient;
		node["vPhase"] >> Body.vPhase;
	}
	catch(...)
	{
		std::cout<<"v amplitude/phase not found\n\n";
	}
	try
	{
		node["xfrequency"] >> Body.xfrequency;
	}
	catch(...)
	{
		std::cout<<"xfrequency not found\n\n";
	}
	try
	{
		node["yfrequency"] >> Body.yfrequency;
	}
	catch(...)
	{
		std::cout<<"yfrequency not found\n\n";
	}
	// get the type of body and read in appropriate details
	string type;
	node["type"] >> type;
	if (type == "points")
	{
		string fname;
		node["pointsFile"] >> fname;
		fname = "bodyFiles/" + fname;
		// initialise points
		std::ifstream file(fname.c_str());
		file >> Body.numPoints;
		Body.X.resize(Body.numPoints);
		Body.Y.resize(Body.numPoints);
		for(int i=0; i<Body.numPoints; i++)
		{
			file >> Body.X[i] >> Body.Y[i];
		}
		file.close();
	}
	else if (type == "lineSegment")
	{
		double startX, startY, endX, endY;
		int numPoints;
		node["segmentOptions"][0] >> startX;
		node["segmentOptions"][1] >> endX;
		node["segmentOptions"][2] >> startY;
		node["segmentOptions"][3] >> endY;
		node["segmentOptions"][4] >> numPoints;
		Body.numPoints = numPoints;
		// initialise line segment
	}
	else if (type == "circle")
	{
		double cx, cy, R, mid_h;
		int numPoints;
		bool selfSize;
		node["circleOptions"][0] >> cx;
		node["circleOptions"][1] >> cy;
		node["circleOptions"][2] >> R;
		node["circleOptions"][3] >> numPoints;
		try
		{
			node["circleOptions"][4] >> selfSize;
			node["circleOptions"][5] >> mid_h;
		}
		catch(...)
		{
		}
		if (selfSize == true)
		{
			numPoints = 2.0*R*M_PI/mid_h;
			if (numPoints%2!=0)
				numPoints++;
		}

		Body.numPoints = numPoints;
		// initialise circle
		Body.X.resize(numPoints);
		Body.Y.resize(numPoints);
		for(int i=0; i<Body.numPoints; i++)
		{
			Body.X[i] = cx + R*cos(i*2*M_PI/Body.numPoints);
			Body.Y[i] = cy + R*sin(i*2*M_PI/Body.numPoints);
		}
	}
	else
		printf("[E]: unknown Body type\n");
	printf("\nnumPoints: %d",Body.numPoints);
}

/**
 * \brief Parses the \a bodies.yaml file and stores information about the immersed bodies.
 *
 * \param bodiesFile the file that contains information about the immersed bodies
 * \param DB the database that contains the simulation parameters 
 */
void parseBodiesFile(std::string &bodiesFile, parameterDB &DB)
{
	std::ifstream fin(bodiesFile.c_str());
	YAML::Parser  parser(fin);
	YAML::Node    doc;
	parser.GetNextDocument(doc);
	body Body;
	std::vector<body> *B = DB["flow"]["bodies"].get<std::vector<body> *>();
	for (int i=0; i<doc.size(); i++)
	{
		Body.reset();
		doc[i] >> Body;
		B->push_back(Body);
	}
}

} // end namespace io
