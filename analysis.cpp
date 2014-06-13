#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <string>
#include <vector>
#include <getopt.h>

using namespace std;

#define NTR 37
#define NMEAN 21
#define MUF 0.01

//#define MM2 // second version of MTE
#define MM1 // first version of MTE

typedef complex<double> dcomplex;

string optionsname = "options";
string tracesname, outputname;

int Ns=4, Nt=4;
int Ns3 = Ns*Ns*Ns;
int nsite = Ns3*Nt;
int nmu, nconf;
bool freeconf=false;

dcomplex *dZdm, *d2Zdm2, *dZdmu, *d2Zdmu2;
double *edZdm, *ed2Zdm2, *edZdmu, *ed2Zdmu2;

char texthelp[]="Usage: exec [OPTION]... [FILE] [TYPE]\n"
		"Calculates the improved Taylor expansion on gauge configurations listed on FILE (or\n"
		"for a free configuration)\n"
		"\n"
		"Mandatory arguments to long options are mandatory for short options too.\n"
		"  -s, --Ns SSIZE             spatial lattice extent (default = 4)\n"
		"  -t, --Nt TSIZE             temporal lattice extent (default = 4)\n"
		"  -n, --nmu NMU         			number of mu in the range [0,1].\n"
		"  -f, --free									free configuration.\n"
		"                             ignored for staggered (default = 0.1)\n"
		"Report bugs to ydelgado83@gmail.com\n";


#include "analysis_observables.h"

void meanValues( dcomplex *mean, dcomplex *obs );
void Jackknife( int imu, dcomplex *mean, dcomplex *obs );
int init( int &argc, char *argv[] );
int readTraces( dcomplex *traces );
void writeResults( );
void allocateArrays( );
void deallocateArrays( );

//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	// Read input parameters
	if( init(argc, argv) != 0)
	{
		cout << "ERROR: Error initializing!" << endl;
		return 1;
	}

	// Allocate global and local arrays
	allocateArrays( );

	dcomplex *traces = new dcomplex[NTR*nconf];
	dcomplex *obs = new dcomplex[NMEAN*nconf];
	dcomplex mean[NMEAN];

	if( readTraces( traces ) != 0)
	{
		cout << "ERROR: Error reading traces file!" << endl;
		return 1;
	}

	// only for FREE configurations!!!!
	if (freeconf){
		cout << "free configuration: dividing all traces by 3." << endl;

		for (int ic=0; ic<nconf*NTR; ic++) traces[ic] /= 3.;
	}

	double dmu = 0.;
	if (nmu>0) dmu = MUF/double(nmu);
	double mu = 0.;
	for (int imu=0; imu<=nmu; imu++ )
	{
		partialObservables( mu, obs, traces ); 
		meanValues( mean, obs );

		// Compute value of observables
		dlnZdmu( &dZdmu[imu], mean );
		d2lnZdmu2( &d2Zdmu2[imu], mean );
		dlnZdm( &dZdm[imu], mean );
		d2lnZdm2( &d2Zdm2[imu], mean );

		// Compute error bars
		Jackknife( imu, mean, obs );

		mu += dmu;
	}

	writeResults( );

	deallocateArrays( );
	delete[] traces;
	delete[] obs;

	return 0;
}

//----------------------------------------------------------------------------------------------------
void meanValues( dcomplex *mean, dcomplex *obs )
{
	// mean values needed for the final observables
	for (int imean=0; imean<NMEAN; imean++)
	{
		// observables
		mean[imean] = dcomplex(0,0);
		for (int iconf=0; iconf<nconf; iconf++)
		{
			mean[imean] += obs[iconf*NMEAN + imean];
		}
		mean[imean] /= double(nconf);
	}

}

//-------------------------------------------------------------------------------------------------
void Jackknife( int imu, dcomplex *mean, dcomplex *obs )
{
	dcomplex jkmean[NMEAN];
	dcomplex jkdZdmu, jkd2Zdmu2, jkdZdm, jkd2Zdm2;

	edZdmu[imu] = 0.;
	ed2Zdmu2[imu] = 0.;
	edZdm[imu] = 0.;
	ed2Zdm2[imu] = 0.;

	if (nconf<2) return;

	for (int iconf=0; iconf<nconf; iconf++)
	{
		for (int imean=0; imean<NMEAN; imean++){
				jkmean[imean] = mean[imean]*double(nconf) - obs[iconf*NMEAN + imean];
				jkmean[imean] /= double(nconf-1);
		}
		dlnZdmu( &jkdZdmu, jkmean );
		d2lnZdmu2( &jkd2Zdmu2, jkmean );
		dlnZdm( &jkdZdm, jkmean );
		d2lnZdm2( &jkd2Zdm2, jkmean );
		
		edZdmu[imu] += std::norm(jkdZdmu - dZdmu[imu]);
		ed2Zdmu2[imu] += std::norm(jkd2Zdmu2 - d2Zdmu2[imu]);
		edZdm[imu] += std::norm(jkdZdm - dZdm[imu]);
		ed2Zdm2[imu] += std::norm(jkd2Zdm2 - d2Zdm2[imu]);
	}
	edZdmu[imu] *= double(nconf-1)/double(nconf);	
	ed2Zdmu2[imu] *= double(nconf-1)/double(nconf);	
	edZdm[imu] *= double(nconf-1)/double(nconf);	
	ed2Zdm2[imu] *= double(nconf-1)/double(nconf);

  edZdmu[imu] = sqrt(edZdmu[imu]);
	ed2Zdmu2[imu] = sqrt(ed2Zdmu2[imu]);
	edZdm[imu] = sqrt(edZdm[imu]);
	ed2Zdm2[imu] = sqrt(ed2Zdm2[imu]);
}

//--------------------------------------------------------------------------------------------------
int readTraces( dcomplex *traces )
{
	std::ifstream file;
	file.open( tracesname.c_str(), ios::in );

	if ( !file.is_open( ) ) return 1;
	file.seekg(0, std::ios::beg);

	cout << "Reading " << nconf << " configurations from " << tracesname.c_str() << endl;

	int ii;
	for (int iconf=0; iconf<nconf; iconf++)
	{
		file >> ii;
		if (ii==iconf)
		{
			for ( int iobs=0; iobs<NTR; iobs++ )
			 file >> traces[iconf*NTR + iobs];
		}
	}
	file.close( );

	return 0;
}

//----------------------------------------------------------------------------------------------------
void writeResults( )
{
	ofstream outfile;
	outfile.open( outputname.c_str(), ios::out );

	double dmu = 0.;
	if (nmu>0) dmu = MUF/double(nmu);
	double mu=0;
	for (int imu=0; imu<=nmu; imu++ )
	{
		mu = imu*dmu;
		outfile << mu << " " << std::real(dZdm[imu])  << " " << edZdm[imu]
									<< " " << std::real(d2Zdm2[imu])  << " " << ed2Zdm2[imu]
									<< " " << std::real(dZdmu[imu]) << " " << edZdmu[imu]
									<< " " << std::real(d2Zdmu2[imu])  << " " << ed2Zdmu2[imu] << endl;
	}
	outfile.close( );
}
 
//--------------------------------------------------------------------------------------------------
int init( int &argc, char *argv[] )
{
	if(argc<1)
	{
		cout << endl << texthelp << endl;
		return 1;
	}

  int c;
       
	while (1)
	{
		static struct option long_options[] =
			{
			/* These options don't set a flag.
			We distinguish them by their indices. */
			{"Ns", required_argument, 0, 's'},
			{"Nt", required_argument, 0, 't'},
			{"nmu", required_argument, 0, 'n'},
			/* These options set a flag. */
			{"free", no_argument, 0, 'f'},
			{0, 0, 0, 0}
			};
			
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "s:t:n:f",
		long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
			case 0:
				/* If this option set a flag, do nothing else now. */
				if (long_options[option_index].flag != 0)
					break;
				printf ("option %s", long_options[option_index].name);
				if (optarg)
					printf (" with arg %s", optarg);
				printf ("\n");

			case 's':
				Ns = atoi(optarg);
				cout << "Ns : " << Ns << endl;
				break;

			case 't':
				Nt = atoi(optarg);
				cout << "Nt : " << Nt << endl;
				break;

			case 'n':
				nmu = atoi(optarg);
				cout << "nmu : " << nmu << endl;
				break;

			case 'f':
				freeconf = true;
				break;

			default:
				cout << endl << texthelp << endl;
				abort();
		}
	}

	/* Print any remaining command line arguments (not options). */
	/* Get the names of the files with the configurations */
	string tmpstring;
	if ( optind >= argc ) 
	{
		cout << "Init.ERROR: file with output/input paths is missing..." << endl;
		return 1;
	}
	else
	{ 
		tmpstring = argv[optind];
	}

	stringstream tmpname;	
	std::ifstream file( tmpstring.c_str(), ios::in );

	if ( !file.is_open( ) )
	{
		cout << "Init.ERROR: file with output/input files couldn't be opened..." << endl;
		return 1;
	}

	// get input and output paths
	file.seekg(0, std::ios::beg);
	tmpstring = "";
	file >> tmpstring;	// path to output

	tmpname.str("");
	tmpname << tmpstring << "finalObs" << ".dat";
	outputname= tmpname.str();

	tmpname.str("");
	tmpname << tmpstring;

	tmpstring = "";
	file >> tmpstring;	// traces file
	tmpname << tmpstring;
	tracesname= tmpname.str();

	/* Set the sizes of the lattice and of matrices*/
	Ns3 = Ns*Ns*Ns;
	nsite = Ns3*Nt;

	file.open( tracesname.c_str(), ios::in );
	if ( !file.is_open( ) ) return 1;
	nconf = std::count(std::istreambuf_iterator<char>(file), 
          std::istreambuf_iterator<char>(), '\n');
	file.close();

	return 0;
}

//----------------------------------------------------------------------------------------------------
void allocateArrays( )
{
	dZdm = new dcomplex[nmu+1];
	d2Zdm2 = new dcomplex[nmu+1];
	dZdmu = new dcomplex[nmu+1];
	d2Zdmu2 = new dcomplex[nmu+1];
	edZdm = new double[nmu+1];
	ed2Zdm2 = new double[nmu+1];
	edZdmu = new double[nmu+1];
	ed2Zdmu2 = new double[nmu+1];
}

//----------------------------------------------------------------------------------------------------
void deallocateArrays( )
{
	delete[] dZdm;
	delete[] d2Zdm2;
	delete[] dZdmu;
	delete[] d2Zdmu2;
	delete[] edZdm;
	delete[] ed2Zdm2;
	delete[] edZdmu;
	delete[] ed2Zdmu2;
}

