#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <string>
#include <vector>
#include <getopt.h>

using namespace std;

// free case for the first version of MTE

typedef complex<double> dcomplex;

string optionsname = "options";

int Ns=4, Nt=8;
int nsite = Ns*Ns*Ns*Nt;

double mass;

dcomplex *dZdm, *d2Zdm2, *dZdmu, *d2Zdmu2;

char texthelp[]="Usage: exec [OPTION]... [FILE] [TYPE]\n"
		"Calculates the improved Taylor expansion on gauge configurations listed on FILE (or\n"
		"for a free configuration)\n"
		"\n"
		"Mandatory arguments to long options are mandatory for short options too.\n"
		"  -s, --Ns SSIZE             spatial lattice extent (default = 4)\n"
		"  -t, --Nt TSIZE             temporal lattice extent (default = 4)\n"
		"  -m, --mass MASS          mass parameter,\n"
		"                               ignored for staggered (default = 0.1)\n"
		"Based on Hans-Peter's code :-)\n"
		"Report bugs to ydelgado83@gmail.com\n";


void finalObservables( );
int init( int &argc, char *argv[] );

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
  finalObservables( );

	return 0;
}

//----------------------------------------------------------------------------------------------------
#define PI2 6.28318530718
void finalObservables( )
{
/*
	Compute the derivatives of lnZ with respect to
		m  : up to 2nd order (>=3rd order terms are too complicated)
		mu : up to 4th order
*/

	double p, cosx, sinx2, cosy, siny2, cosz, sinz2, cost, sint;
	double ks = PI2/double(Ns);
	double kt = PI2/double(Nt);

	int ix,iy,iz,it;

  dcomplex trD = dcomplex(0,0);
  dcomplex trMD = dcomplex(0,0);
  dcomplex trMbD = dcomplex(0,0);

	for (ix=0; ix<Ns; ix++){
		p = ks*ix;
		cosx = cos(p);
		sinx2 = sin(p); sinx2 *= sinx2;

		for (iy=0; iy<Ns; iy++){
			p = ks*iy;
			cosy = cos(p);
			siny2 = sin(p); siny2 *= siny2;

			for (iz=0; iz<Ns; iz++){
				p = ks*iz;
				cosz = cos(p);
				sinz2 = sin(p); sinz2 *= sinz2;

				double s = mass + 4. - cosx - cosy - cosz;
				double ss = sinx2 + siny2 + sinz2;

				for (it=0; it<Nt; it++){
					p = kt*(it+0.5);
	 				cost = cos(p);
					sint = sin(p);

					double den = (s-cost)*(s-cost) + ss + sint*sint ;
										
					trD += (s - cost) / den;
					trMD += (s - cost + dcomplex(0,1.)*sint)*dcomplex(cost,sint) / den;
					trMbD += (-ss - sint*sint - dcomplex(0,1.)*sint*s)*dcomplex(2.*cost,+2.*sint) / (den*den);
							//(s - cost - dcomplex(0,1.)*sint)*dcomplex(cost,-sint) / den;
				}
			}
		}
	}

	trD *= 4.*3.;///double(nsite);
	trMD *= 2.*3.;///double(nsite);
	trMbD *= 2.*3.;///double(nsite);

	cout << "trD " << trD << endl;
	cout << "trMD " << trMD << endl;
	cout << "trMbD " << trMbD << endl;
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
			{"mass", required_argument, 0, 'm'},
			/* These options set a flag. */
			{0, 0, 0, 0}
			};
			
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "s:t:m:k:",
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

			case 'm':
				mass = atof(optarg);
				cout << "mass: " << mass << endl;
				break;

			default:
				cout << endl << texthelp << endl;
				abort();
		}
	}

	/* Set the sizes of the lattice*/
	nsite = Ns*Ns*Ns*Nt;

	return 0;
}

