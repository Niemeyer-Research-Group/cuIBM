/***************************************************************************//**
 * \file parseDomainFile.cu
 * \author Anush Krishnan (anush@bu.edu), Christopher Minar (mianrc@oregonstate.edu
 * \brief Parse the input file domain.yaml to obtain information about the
 *        computational grid.
 */


#include <fstream>
#include <yaml-cpp/yaml.h>
#include "io.h"

/**
 * \namespace io
 * \brief Contains functions related to I/O tasks.
 */
namespace io
{

using std::string;

/**
 * \brief Overloads the operator >>. Gets information from the parsed domain file.
 *
 * \param node the parsed file
 * \param D instance of the class \c domain to be filled
 */
void operator >> (const YAML::Node &node, domain &D)
{
	string dir;
	double start;
	int  numCells;
	double mid_h;
	node["direction"] >> dir;
	node["start"] >> start;
	node["mid_h"] >> mid_h;
	D.mid_h = mid_h;

	double end, stretchRatio, a, count = 0;
	const YAML::Node &subDomains = node["subDomains"];
	for (unsigned int i=0; i<subDomains.size(); i++) //first pass, gets nx and ny and resizes the arrays appropriately
	{
		subDomains[i]["end"] >> end;
		subDomains[i]["stretchRatio"] >> stretchRatio;
		if (stretchRatio<1)
		{
			//a = mid_h/((end-start)*(stretchRatio-1));
			//numCells = log(a/(a-1))/log(stretchRatio);
			//count += int(numCells);
			a = mid_h/((start-end)*(stretchRatio-1));
			numCells = log(1/a + 1)/log(stretchRatio);
			count += int(numCells);
		}
		else if(stretchRatio == 1)
		{
			numCells = (end-start)/mid_h;
			count += int(numCells);
		}
		else
		{
			a = mid_h/((end-start)*(stretchRatio-1));
			numCells = log(1/a + 1)/log(stretchRatio);
			count += int(numCells);
		}
		start = end;
	}
	if (dir=="x")
	{
		D.nx = count;
		D.x.resize(D.nx);	//x location of the pressure nodes (cell center)
		D.dx.resize(D.nx);	//x width of the pressure nodes
		D.xv.resize(D.nx);	//x location of where v is stored (same x as pressure node)
	}
	else if(dir=="y")
	{
		D.ny = count;
		D.y.resize(D.ny);	//y location of pressure nodes (cell center)
		D.dy.resize(D.ny);	//y height of pressure nodes
		D.yu.resize(D.ny);	//y location of where u is stored
	}

	// second pass, fills x,y,xv,yu,dx,dy
	node["start"] >> start;
	double h;
	int beg = 0;
	for (unsigned int i=0; i<subDomains.size(); i++)
	{
		subDomains[i]["end"] >> end;
		subDomains[i]["stretchRatio"] >> stretchRatio;
		if (stretchRatio<1)
		{
			//a = mid_h/((end-start)*(stretchRatio-1));
			//numCells = log(a/(a-1))/log(stretchRatio);
			a = mid_h/((start-end)*(stretchRatio-1));
			numCells = log(1/a + 1)/log(stretchRatio);
		}
		else if(stretchRatio == 1)
		{
			numCells = (end-start)/mid_h;
		}
		else
		{
			a = mid_h/((end-start)*(stretchRatio-1));
			numCells = log(1/a + 1)/log(stretchRatio);
		}

		if(fabs(stretchRatio-1.0) < 1.0e-6)  //no cell stretching
		{
			h = (end - start)/numCells;
			for (int j=beg; j<beg+numCells; j++)
			{
				if (dir=="x")
				{
					if (j == beg)
					{
						D.x[j] = start + h/2;
						D.xv[j] = D.x[j];     	//xv is measured in the same x location as the cell center
						D.dx[j] = h;
					}
					else
					{
						D.dx[j] = h;
						D.x[j] = D.x[j-1] + 0.5*D.dx[j] + 0.5*D.dx[j-1];//dx is the width of the cell so to move over to the next cell center you need half of each cells dx
						D.xv[j] = D.x[j];
					}
				}
				else if (dir=="y")
				{
					if (j == beg)
					{
						D.y[j] = start + h/2;
						D.yu[j] = D.y[j];
						D.dy[j] = h;
					}
					else
					{
						D.dy[j] = h;
						D.y[j] = D.y[j-1] + 0.5*D.dy[j] + 0.5*D.dy[j-1];
						D.yu[j] = D.y[j];
					}
				}//end x/y elseif
			}//end for
		}//end no cell stretching
		else if (end < 0)
		{
			h = (end - start)*(stretchRatio-1)/(pow(stretchRatio, numCells)-1);
			int jdx = beg;
			for (int j=beg+numCells-1; j>=beg; j--)
			{
				if (dir=="x")
				{
					if (j == beg+numCells-1)
					{
						D.x[j] = end - h/2.0;
						D.xv[j] = D.x[j];
						D.dx[j] = h;
					}
					else
					{
						D.dx[j] = h*pow(stretchRatio, jdx-beg);
						D.x[j] = D.x[j+1] - 0.5*D.dx[j] - 0.5*D.dx[j+1];
						D.xv[j] = D.x[j];
					}
				}
				if (dir=="y")
				{
					if (j == beg+numCells-1)
					{
						D.y[j] = end - h/2.0;
						D.yu[j] = D.y[j];
						D.dy[j] = h;
					}
					else
					{
						D.dy[j] = h*pow(stretchRatio, jdx-beg);
						D.y[j] = D.y[j+1] - 0.5*D.dy[j] - 0.5*D.dy[j+1];
						D.yu[j] = D.y[j];
					}
				}//end x/y elseif
				jdx++;
			}//end for
		}//end cell shrink
		else //cell growing
		{
			h = (end - start)*(stretchRatio-1)/(pow(stretchRatio, numCells)-1);
			for (int j=beg; j<beg+numCells; j++)
			{
				if (dir=="x")
				{
					if (j == beg)
					{
						D.x[j] = start + h*pow(stretchRatio, j-beg)/2;
						D.xv[j] = D.x[j];
						D.dx[j] = h*pow(stretchRatio, j-beg);
					}
					else
					{
						D.dx[j] = h*pow(stretchRatio, j-beg);
						D.x[j] = D.x[j-1] + 0.5*D.dx[j] + 0.5*D.dx[j-1];
						D.xv[j] = D.x[j];
					}
				}
				if (dir=="y")
				{
					if (j == beg)
					{
						D.y[j] = start + h*pow(stretchRatio, j-beg)/2;
						D.yu[j] = D.y[j];
						D.dy[j] = h*pow(stretchRatio, j-beg);
					}
					else
					{
						D.dy[j] = h*pow(stretchRatio, j-beg);
						D.y[j] = D.y[j-1] + 0.5*D.dy[j] + 0.5*D.dy[j-1];
						D.yu[j] = D.y[j];
					}
				}//end x/y elseif
			}//end for
		}//end cell grow
		beg += numCells;
		start = end;
	}
/*
	string dir;
	double start;
	int  numCells;

	node["direction"] >> dir;
	node["start"] >> start;

	if (dir=="x")
		D.nx = 0;
	else if(dir=="y")
		D.ny = 0;

	const YAML::Node &subDomains = node["subDomains"];
	for (unsigned int i=0; i<subDomains.size(); i++) //first pass, gets nx and ny and resizes the arrays appropriately
	{
		subDomains[i]["cells"] >> numCells;
		if (dir=="x")
			D.nx += numCells;
		else if(dir=="y")
			D.ny += numCells;
	}

	// allocate memory
	int  beg = 0;
	if(dir=="x")//x
	{
		D.x.resize(D.nx);	//x location of the pressure nodes (cell center)
		D.dx.resize(D.nx);	//x width of the pressure nodes
		D.xv.resize(D.nx);	//x location of where v is stored (same x as pressure node)
	}
	if(dir=="y")//y
	{
		D.y.resize(D.ny);	//y location of pressure nodes (cell center)
		D.dy.resize(D.ny);	//y height of pressure nodes
		D.yu.resize(D.ny);	//y location of where u is stored
	}

	// second pass, fills x,y,xv,yu,dx,dy
	double end, stretchRatio, h;
	for (unsigned int i=0; i<subDomains.size(); i++)
	{
		subDomains[i]["end"] >> end;
		subDomains[i]["cells"] >> numCells;
		subDomains[i]["stretchRatio"] >> stretchRatio;
		
		if(fabs(stretchRatio-1.0) < 1.0e-6)  //no cell stretching
		{
			h = (end - start)/numCells;
			for (int j=beg; j<beg+numCells; j++)
			{
				if (dir=="x")
				{
					if (j == beg)
					{
						D.x[j] = start + h/2; 	//start designates the boundary location so the location of the cell center is dx/2 away
						D.xv[j] = D.x[j];     	//xv is measured in the same x location as the cell center
						D.dx[j] = h;			//without stretching dx = h
					}
					else
					{
						D.dx[j] = h;			//without stretching dx = h
						D.x[j] = D.x[j-1] + 0.5*D.dx[j] + 0.5*D.dx[j-1];//dx is the width of the cell so to move over to the next cell center you need half of each cells dx
						D.xv[j] = D.x[j];		//x_xv = x_p
					}
				}
				else if (dir=="y")
				{
					if (j == beg)
					{
						D.y[j] = start + h/2;
						D.yu[j] = D.y[j];
						D.dy[j] = h;
					}
					else
					{
						D.dy[j] = h;
						D.y[j] = D.y[j-1] + 0.5*D.dy[j] + 0.5*D.dy[j-1];
						D.yu[j] = D.y[j];
					}
				}//end x/y elseif
			}//end for
		}//end no cell stretching
		else //cell stretching
		{
			h = (end - start)*(stretchRatio-1)/(pow(stretchRatio, numCells)-1); //the initial dx, will either be the largest or smallest value depending on if stretch is greater or less than 1
			for (int j=beg; j<beg+numCells; j++)
			{
				if (dir=="x")
				{
					if (j == beg)
					{
						D.x[j] = start + h*pow(stretchRatio, j-beg)/2;
						D.xv[j] = D.x[j];
						D.dx[j] = h*pow(stretchRatio, j-beg); //dx = biggest possible cell * stretch^(j-beg)
					}
					else
					{
						D.dx[j] = h*pow(stretchRatio, j-beg);
						D.x[j] = D.x[j-1] + 0.5*D.dx[j] + 0.5*D.dx[j-1];
						D.xv[j] = D.x[j];
					}
				}
				if (dir=="y")
				{
					if (j == beg)
					{
						D.y[j] = start + h*pow(stretchRatio, j-beg)/2;
						D.yu[j] = D.y[j];
						D.dy[j] = h*pow(stretchRatio, j-beg);
					}
					else
					{
						D.dy[j] = h*pow(stretchRatio, j-beg);
						D.y[j] = D.y[j-1] + 0.5*D.dy[j] + 0.5*D.dy[j-1];
						D.yu[j] = D.y[j];//D.y[j-1] + 0.5*D.dy[j] + 0.5*D.dy[j-1];
					}
				}//end x/y elseif
			}//end for
		}//end cell stretching
		beg += numCells;
		start = end;
	}
	*/
}

/**
 * \brief Parses the \a domain file and generates the computational grid.
 *
 * \param domFile the file that contains information about the computational grid
 * \param D instance of the class \c domain that will be filled with information about the computational grid
 */
void parseDomainFile(std::string &domFile, domain &D)
{
	std::ifstream fin(domFile.c_str());			//setup to go through casefolder/domain.yaml
	YAML::Parser  parser(fin);
	YAML::Node    doc;
	parser.GetNextDocument(doc);
	for (unsigned int i=0; i<doc.size(); i++)	//go through each node in domain.yaml
		doc[i] >> D;

	D.xu.resize(D.nx-1);	//xu on device
	D.yu.resize(D.ny);		//yu on device
	D.xv.resize(D.nx);		//xv on device
	D.yv.resize(D.ny-1);	//yv on device

	for(int i=0; i<D.nx-1; i++)	//set xu
	{
		D.xu[i] = D.x[i] + D.dx[i]/2;
	}

	for(int j=0; j<D.ny-1; j++)	//set yv
	{
		D.yv[j] = D.y[j] + D.dy[j]/2;
	}
}
} // end namespace io
