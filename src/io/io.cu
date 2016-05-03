/***************************************************************************//**
 * \file io.cu
 * \author Anush Krishnan (anush@bu.edu)
 * \brief Implementation of the functions of the namespace \c io.
 */


#include <sys/stat.h>
#include "io.h"
#include <boundaryCondition.h>
#include "preconditioner.h"

using std::string;
using std::ios;


/**
 * \brief Converts a string to a number.
 *
 * \param str a string
 *
 * \return a number (\c double or \c integer)
 */
template <typename T>
T toNumber(string str)
{
     T num;
     std::stringstream ss(str); //turn the string into a stream
     ss >> num; //convert
     return num;
}

/**
 * \namespace io
 * \brief Contains functions related to I/O tasks.
 */
namespace io
{

/**
 * \brief Splits a string given a delimiter.
 *
 * \param s the string to split
 * \param delim the delimiter
 * \param elems the vector that contains the different elements of the string
 *
 * \return a vector that contains the different elements of the string
 */
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

/**
 * \brief Splits a string given a delimiter.
 *
 * \param s the string to split
 * \param delim the delimiter
 *
 * \return a vector that contains the different elements of the string
 */
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

/**
 * \brief Creates a directory.
 *
 * If the parent directory does not exist, it will be created.
 *
 * \param folderPath the path of the directory to create
 */
void makeDirectory(const std::string folderPath)
{
	std::vector<std::string> x = split(folderPath, '/');
	int n = x.size();
	int i = 0;
	std::string folder = x[i];
	mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	i++;
	while(i<n)
	{
		folder = folder + '/' + x[i];
		mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		i++;
	}
}

//##############################################################################
//                                 INPUT
//##############################################################################

/**
 * \brief Reads data inputs from the command-line and the simulation files.
 *
 * \param argc number of arguments in the command-line
 * \param argv command-line arguments
 * \param DB database that contains all the simulation parameters
 * \param D object of the class \c domain that contains the computational grid
 */
void readInputs(int argc, char **argv, parameterDB &DB, domain &D)
{
	// get a default database
	initialiseDefaultDB(DB);

	// first pass of command line arguments
	commandLineParse1(argc, argv, DB);

	// case folder
	string folder = DB["inputs"]["caseFolder"].get<std::string>();

	// read the simulation file
	string fname = folder + "/simParams.yaml";
	parseSimulationFile(fname, DB);

	// read the flow file
	fname = folder + "/flow.yaml";
	parseFlowFile(fname, DB);

	// read the domain file
	fname = folder + "/domain.yaml";
	parseDomainFile(fname, D);;

	// read the body file
	fname = folder + "/bodies.yaml";
	parseBodiesFile(fname, DB);

	// second pass of command line -- overwrite values in DB
	commandLineParse2(argc, argv, DB);
}

/**
 * \brief Initializes the database with default values.
 *
 * \param DB database that contains all the simulation parameters
 */
void initialiseDefaultDB(parameterDB &DB)
{
	DB["inputs"] = componentParameter();
	DB["flow"] = componentParameter();
	DB["simulation"] = componentParameter();
	DB["velocitySolve"] = componentParameter();
	DB["PoissonSolve"] = componentParameter();

	// default input files
	string inputs = "inputs";
	DB[inputs]["caseFolder"].set<string>("/scratch/cases/cuIBM/cases/cylinder/Re40");
	DB[inputs]["deviceNumber"].set<int>(0);

	// flow parameters
	string flow = "flow";
	DB[flow]["nu"].set<double>(0.01);
	DB[flow]["uInitial"].set<double>(1.0);
	DB[flow]["vInitial"].set<double>(0.0);
	DB[flow]["numBodies"].set<int>(0);
	std::vector<body> *bodyVec = new std::vector<body>;
	DB[flow]["bodies"].set<std::vector<body> *>(bodyVec);

	// boundary conditions
	boundaryCondition **bc = new boundaryCondition*[4];
	for (int i=0; i<4; i++)
		bc[i] = new boundaryCondition[2];
	DB[flow]["boundaryConditions"].set<boundaryCondition **>(bc);

	// simulation parameters
	string sim = "simulation";
	DB[sim]["dt"].set<double>(0.02);
	DB[sim]["nt"].set<int>(100);
	DB[sim]["nsave"].set<int>(100);
	DB[sim]["restart"].set<bool>(false);
	DB[sim]["startStep"].set<bool>(0);
	DB[sim]["FSI"].set<int>(0);
	DB[sim]["Ured"].set<double>(3);

	// velocity solver
	string solver = "velocitySolve";
	DB[solver]["solver"].set<string>("CG");
	DB[solver]["preconditioner"].set<preconditionerType>(DIAGONAL);
	DB[solver]["tolerance"].set<double>(1e-5);
	DB[solver]["maxIterations"].set<int>(10000);

	// Poisson solver
	solver = "PoissonSolve";
	DB[solver]["solver"].set<string>("CG");
	DB[solver]["preconditioner"].set<preconditionerType>(DIAGONAL);
	DB[solver]["tolerance"].set<double>(1e-5);
	DB[solver]["maxIterations"].set<int>(20000);
}

/**
 * \brief Parses the command-line to get the case folder name 
 *        and the device number.
 *
 * \param argc number of arguments in the command-line
 * \param argv arguments of the command-line
 * \param DB database that contains all the simulation parameters
 */
void commandLineParse1(int argc, char **argv, parameterDB &DB)
{
	for (int i=1; i<argc; i++)
	{
		if (strcmp(argv[i],"-caseFolder")==0)
		{
			i++;
			DB["inputs"]["caseFolder"].set<string>(string(argv[i]));
		}
		else if (strcmp(argv[i],"-deviceNumber")==0)
		{
			i++;
			int devNum = toNumber<int>(string(argv[i]));
			DB["inputs"]["deviceNumber"].set<int>(devNum);
			// sets devNum as the current device for the calling host thread
			cudaSetDevice(devNum);
		}
	}
}

/**
 * \brief Overwrites parameters with additional arguments of the command-line. 
 *
 * \param argc number of arguments in the command-line
 * \param argv arguments of the command-line
 * \param DB database that contains all the simulation parameters
 */
void commandLineParse2(int argc, char **argv, parameterDB &DB)
{
	for (int i=1; i<argc; i++)
	{
		// kinematic viscosity
		if ( strcmp(argv[i],"-nu")==0 )
		{
			i++;
			DB["flow"]["nu"].set<double>(toNumber<double>(string(argv[i])));
		}
		//// angle of attack
		//if ( strcmp(argv[i],"-alpha")==0 )
		//{
		//	i++;
		//	DB["flow"]["nu"].set<double>(toNumber<double>(string(argv[i])));
		//}
		// perturbation in the x-velocity
		if ( strcmp(argv[i],"-uPerturb")==0 )
		{
			i++;
			DB["flow"]["uPerturb"].set<double>(toNumber<double>(string(argv[i])));
		}
		// perturbation in the y-velocity
		if ( strcmp(argv[i],"-vPerturb")==0 )
		{
			i++;
			DB["flow"]["vPerturb"].set<double>(toNumber<double>(string(argv[i])));
		}
		// scale the CV with respect to the body
		if ( strcmp(argv[i],"-scaleCV")==0 )
		{
			i++;
			DB["simulation"]["scaleCV"].set<double>(toNumber<double>(string(argv[i])));
		}
		// frequency of saving the data
		if ( strcmp(argv[i],"-nsave")==0 )
		{
			i++;
			DB["simulation"]["nsave"].set<int>(toNumber<int>(string(argv[i])));
		}
		// total number of time steps
		if ( strcmp(argv[i],"-nt")==0 )
		{
			i++;
			DB["simulation"]["nt"].set<int>(toNumber<int>(string(argv[i])));
		}
		// size of time increment
		if ( strcmp(argv[i],"-dt")==0 )
		{
			i++;
			DB["simulation"]["dt"].set<double>(toNumber<double>(string(argv[i])));
		}
		// tolerance for the velocity solve
		if ( strcmp(argv[i],"-velocityTol")==0 )
		{
			i++;
			DB["velocitySolve"]["tolerance"].set<double>(toNumber<double>(string(argv[i])));
		}
		// tolerance for the Poisson solve
		if ( strcmp(argv[i],"-poissonTol")==0 )
		{
			i++;
			DB["PoissonSolve"]["tolerance"].set<double>(toNumber<double>(string(argv[i])));
		}
	}
}

//##############################################################################
//                                OUTPUT
//##############################################################################

/**
 * \brief Converts a \c preconditionerType to a \c std::string.
 *
 * \param s a preconditioner
 *
 * \return a string
 */
/*
string stringFromPreconditionerType(preconditionerType s)
{
  if (s == NONE)
    return "None";
  else if (s == DIAGONAL)
    return "Diagonal";
  else if (s == SMOOTHED_AGGREGATION)
    return "Smoothed Aggregation";
  else if (s == AINV)
    return "Approximate Inverse";
  else
    return "Unrecognised preconditioner";
}*/

/**
 * \brief Prints the parameters of the simulation.
 *
 * \param DB database that contains all the simulation parameters
 * \param D information about the computational grid
 */
void printSimulationInfo(parameterDB &DB, domain &D)
{
	double dt = DB["simulation"]["dt"].get<double>(),
	     scaleCV = DB["simulation"]["scaleCV"].get<double>();
	int  nt = DB["simulation"]["nt"].get<int>(),
	     nsave = DB["simulation"]["nsave"].get<int>(),
	     startStep = DB["simulation"]["startStep"].get<int>();

    std::cout << '\n';
	
	std::cout << "\nFlow parameters" << '\n';
	std::cout << "---------------" << '\n';
	std::cout << "nu = " << DB["flow"]["nu"].get<double>() << '\n';
	std::cout << "Re(1/nu) = " << 1/DB["flow"]["nu"].get<double>() << '\n';

	std::cout << "\nDomain" << '\n';
	std::cout << "------" << '\n';
	std::cout << D.nx << " x " << D.ny << '\n';
	
	std::cout << "\nSimulation parameters" << '\n';
	std::cout << "---------------------" << '\n';
	std::cout << "dt = " << dt << '\n';
	std::cout << "scaleCV = " << scaleCV << '\n';
	std::cout << "startStep = " << startStep << '\n';
	std::cout << "nt = "    << nt << '\n';
	std::cout << "nsave = " << nsave << '\n';
	
	std::cout << "\nVelocity Solve" << '\n';
	std::cout << "--------------" << '\n';
	std::cout << "Solver = " << DB["velocitySolve"]["solver"].get<std::string>() << '\n';
	//std::cout << "Preconditioner = " << stringFromPreconditionerType(DB["velocitySolve"]["preconditioner"].get<preconditionerType>()) << '\n';
	std::cout << "Tolerance = " << DB["velocitySolve"]["tolerance"].get<double>() << '\n';
	
	std::cout << "\nPoisson Solve" << '\n';
	std::cout << "-------------" << '\n';
	std::cout << "Solver = " << DB["PoissonSolve"]["solver"].get<std::string>() << '\n';
	//std::cout << "Preconditioner = " << stringFromPreconditionerType(DB["PoissonSolve"]["preconditioner"].get<preconditionerType>()) << '\n';
	std::cout << "Tolerance = " << DB["PoissonSolve"]["tolerance"].get<double>() << '\n';
	
	std::cout << "\nOutput parameters" << '\n';
	std::cout << "-----------------" << '\n';
	std::cout << "Output folder = " << DB["inputs"]["caseFolder"].get<std::string>() << '\n';
	std::cout << "nsave = " << DB["simulation"]["nsave"].get<int>() << '\n';
	
	cudaDeviceProp deviceProp;
	int gpu = DB["inputs"]["deviceNumber"].get<int>();
	cudaGetDeviceProperties(&deviceProp, gpu);
	std::cout << "\nDevice Properties" << '\n';
	std::cout << "-----------------" << '\n';
	std::cout << "Name = " << deviceProp.name << '\n';
	std::cout << "Number = " << gpu << '\n';
	std::string ecc = deviceProp.ECCEnabled ? "yes" : "no";
	std::cout << "Compute capability = " << deviceProp.major << "." << deviceProp.minor << '\n';
	std::cout << "ECC Enabled = " << ecc << std::endl;
}

/**
 * \brief Prints the time spent to execute tasks.
 *
 * \param logger object that contains the name and time spent of tasks
 */
void printTimingInfo(Logger &logger)
{
	logger.printAllTime();
	std::cout << std::endl;
}

/**
 * \brief Writes information about the run into the file \a run.info.
 *
 * \param DB database that contains all the simulation parameters
 * \param D information about the computational grid
 */
void writeInfoFile(parameterDB &DB, domain &D)
{
	std::string   folder = DB["inputs"]["caseFolder"].get<std::string>();
	std::ofstream infofile((folder+"/run.info").c_str());
	infofile << std::setw(20) << "--nx"  << "\t" << D.nx << '\n';
	infofile << std::setw(20) << "--ny"  << "\t" << D.ny << '\n';
	infofile << std::setw(20) << "--startStep" << "\t" << DB["simulation"]["startStep"].get<int>() << '\n';
	infofile << std::setw(20) << "--nt"     << "\t" << DB["simulation"]["nt"].get<int>() << '\n';
	infofile << std::setw(20) << "--nsave"  << "\t" << DB["simulation"]["nsave"].get<int>() << '\n';
	infofile << std::setw(20) << "--dt"     << "\t" << DB["simulation"]["dt"].get<double>() << '\n';
	infofile << std::setw(20) << "--vortlim"<< "\t" << 15 << '\n';
	infofile << std::setw(20) << "--folder" << "\t" << folder << '\n';
	infofile << std::setw(20) << "--nu"     << "\t" << DB["flow"]["nu"].get<double>() << '\n';
	infofile.close();
}

/**
 * \brief Writes grid-points coordinates into the file \a grid.
 *
 * \param caseFolder the directory of the simulation
 * \param D information about the computational grid
 */
void writeGrid(std::string &caseFolder, domain &D)
{
	std::stringstream out;
	out << caseFolder << "/grid";
	std::ofstream file(out.str().c_str(), ios::binary);
	file.write((char*)(&D.nx), sizeof(int));
	file.write((char*)(&D.x[0]), (D.nx+1)*sizeof(double));
	file.write((char*)(&D.ny), sizeof(int));
	file.write((char*)(&D.y[0]), (D.ny+1)*sizeof(double));
	file.close();
}

/**
 * \brief Writes numerical data at a given time-step (on the device).
 *
 * It creates a directory whose name is the time-step number
 * and writes the flux, the pressure (and eventually the body forces)
 * into the files \a q, \a lambda, respectively.
 *
 * \param caseFolder directory of the simulation
 * \param n the time-step number
 * \param q array that contains the fluxes
 * \param lambda array that contains the pressures (and eventually the body forces)
 * \param D information about the computational grid
 */
template <>
void writeData< cusp::array1d<double, cusp::host_memory> >(std::string &caseFolder, int n, cusp::array1d<double, cusp::host_memory> &u, cusp::array1d<double, cusp::host_memory> &p, domain &D)//, bodies &B)
{
	std::string path;
	std::stringstream out;
	int N;

	out << caseFolder << '/' << std::setfill('0') << std::setw(7) << n;
	path = out.str();

	mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	out.str("");
	out << path << "/q";
	std::ofstream file(out.str().c_str(), ios::binary);
	N = u.size();
	file.write((char*)(&N), sizeof(int));
	file.write((char*)(&u[0]), N*sizeof(double));
	file.close();

	out.str("");
	out << path << "/lambda";
	file.open(out.str().c_str(), ios::binary);
	N = p.size();
	file.write((char*)(&N), sizeof(int));
	file.write((char*)(&p[0]), N*sizeof(double));
	file.close();

	std::cout << "Data saved to folder " << path << std::endl;
}

/**
 * \brief Writes numerical data at a given time-step (on the device).
 *
 * It creates a directory whose name is the time-step number
 * and writes the flux, the pressure (and eventually the body forces)
 * into the files \a q, \a lambda, respectively.
 *
 * \param caseFolder directory of the simulation
 * \param n the time-step number
 * \param q array that contains the fluxes
 * \param lambda array that contains the pressures (and eventually the body forces)
 * \param D information about the computational grid
 */
template <>
void writeData< cusp::array1d<double, cusp::device_memory> >(std::string &caseFolder, int n, cusp::array1d<double, cusp::device_memory> &u, cusp::array1d<double, cusp::device_memory> &p, domain &D)//, bodies &B)
{
	cusp::array1d<double, cusp::host_memory> 	uH = u,
												pH = p;

	writeData(caseFolder, n, uH, pH, D);
}

/**
 * \brief Prints device memory usage.
 *
 * \param label the label of the device
 */
void printDeviceMemoryUsage()
{
	size_t _free, _total;
	cudaMemGetInfo(&_free, &_total);
	std::cout << '\n' << "Initialisation complete\nFlux capacitors charged" << ": Memory Usage " << std::setprecision(3) << (_total-_free)/(1024.0*1024*1024) \
	          << " / " << std::setprecision(3) << _total/(1024.0*1024*1024) << " GB" << std::setprecision(6) << '\n' << std::endl;
}

} // end namespace io
