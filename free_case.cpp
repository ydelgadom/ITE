#include <iostream>
#include <fstream>
#include <algorithm>
#include <complex>
#include <string>
#include <vector>
#include <getopt.h>

using namespace std;

#define MUF 2.0

// free case for the first version of MTE

typedef complex<double> dcomplex;
typedef dcomplex (*func_array) (double *dpar, dcomplex *cpar);

string optionsname = "options";
string outputname;

int Ns=4, Nt=4;
int nsite = Ns*Ns*Ns*Nt;
int nmu, norder;

double mass;

dcomplex *dZdm, *d2Zdm2, *dZdmu, *d2Zdmu2;

char texthelp[]="Usage: exec [OPTION]... [FILE] [TYPE]\n"
		"Calculates the improved Taylor expansion on gauge configurations listed on FILE (or\n"
		"for a free configuration)\n"
		"\n"
		"Mandatory arguments to long options are mandatory for short options too.\n"
		"  -s, --Ns SSIZE           spatial lattice extent (default = 4)\n"
		"  -t, --Nt TSIZE           temporal lattice extent (default = 4)\n"
		"  -o, --O  ORDER						expansion order\n"
		"  -m, --mass MASS         	mass parameter,\n"
		"  -n, --nmu NMU         		number of mu in the range [0,1].\n"
		"Report bugs to ydelgado83@gmail.com\n";


void finalObservables( int imu, double mu );
int init( int &argc, char *argv[] );
void writeResults( );
void allocateArrays( );
void deallocateArrays( );

dcomplex dZdm_0th(double *dpar, dcomplex *cpar);
dcomplex dZdm_1st(double *dpar, dcomplex *cpar);
dcomplex dZdm_2nd(double *dpar, dcomplex *cpar);
dcomplex dZdm_3rd(double *dpar, dcomplex *cpar);
dcomplex dZdm_4rd(double *dpar, dcomplex *cpar);
dcomplex dZdm_exact(double *dpar, dcomplex *cpar);

dcomplex d2Zdm2_0th(double *dpar, dcomplex *cpar);
dcomplex d2Zdm2_1st(double *dpar, dcomplex *cpar);
dcomplex d2Zdm2_2nd(double *dpar, dcomplex *cpar);
dcomplex d2Zdm2_3rd(double *dpar, dcomplex *cpar);
dcomplex d2Zdm2_4rd(double *dpar, dcomplex *cpar);
dcomplex d2Zdm2_exact(double *dpar, dcomplex *cpar);

dcomplex dZdmu_0th(double *dpar, dcomplex *cpar);
dcomplex dZdmu_1st(double *dpar, dcomplex *cpar);
dcomplex dZdmu_2nd(double *dpar, dcomplex *cpar);
dcomplex dZdmu_3rd(double *dpar, dcomplex *cpar);
dcomplex dZdmu_4rd(double *dpar, dcomplex *cpar);
dcomplex dZdmu_exact(double *dpar, dcomplex *cpar);

dcomplex d2Zdmu2_0th(double *dpar, dcomplex *cpar);
dcomplex d2Zdmu2_1st(double *dpar, dcomplex *cpar);
dcomplex d2Zdmu2_2nd(double *dpar, dcomplex *cpar);
dcomplex d2Zdmu2_3rd(double *dpar, dcomplex *cpar);
dcomplex d2Zdmu2_4rd(double *dpar, dcomplex *cpar);
dcomplex d2Zdmu2_exact(double *dpar, dcomplex *cpar);

func_array contribution_dZdm[6]
{ dZdm_0th, dZdm_1st, dZdm_2nd, dZdm_3rd, dZdm_4rd, dZdm_exact};

func_array contribution_d2Zdm2[6]
{ d2Zdm2_0th, d2Zdm2_1st, d2Zdm2_2nd, d2Zdm2_3rd, d2Zdm2_4rd, d2Zdm2_exact};

func_array contribution_dZdmu[6]
{ dZdmu_0th, dZdmu_1st, dZdmu_2nd, dZdmu_3rd, dZdmu_4rd, dZdmu_exact};

func_array contribution_d2Zdmu2[6]
{ d2Zdmu2_0th, d2Zdmu2_1st, d2Zdmu2_2nd, d2Zdmu2_3rd, d2Zdmu2_4rd, d2Zdmu2_exact};

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
	double dmu = 0.;
	if (nmu>0) dmu = MUF/double(nmu);
	double mu = 0.;
	for (int imu=0; imu<=nmu; imu++ )
	{
		finalObservables( imu, mu );
		mu += dmu;
	}

	writeResults( );
	deallocateArrays( );

	return 0;
}

//----------------------------------------------------------------------------------------------------
#define PI2 6.28318530718
void finalObservables( int imu, double mu )
{
/*
	Compute the derivatives of lnZ with respect to
		m  : up to 2nd order (>=3rd order terms are too complicated)
		mu : up to 4th order
*/
	double p, cosx, sinx2, cosy, siny2, cosz, sinz2, cost, sint;
	double R, cp, dRdm;
	double ks = PI2/double(Ns);
	double kt = PI2/double(Nt);

	double dpar[9];
	dcomplex cpar[9];

  dpar[0] = exp(mu) - 1.;
  dpar[1] = exp(-mu) - 1.;
  dpar[2] = exp(mu);
  dpar[3] = exp(-mu);

  dZdmu[imu]=dcomplex(0,0);
  d2Zdmu2[imu]=dcomplex(0,0);
  dZdm[imu]=dcomplex(0,0);
  d2Zdm2[imu]=dcomplex(0,0);

	int ix,iy,iz,it,iorder;

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

				cp = mass + 4. - cosx - cosy - cosz;
				dpar[4] = cp;

				for (it=0; it<Nt; it++){
					p = kt*(it+0.5);
	 				cost = cos(p);
					sint = sin(p);

					R = cp*(cp - 2.*cost) + 1. + sinx2 + siny2 + sinz2;
					dRdm = 2.*(cp-cost);
					dpar[5] = R;
					dpar[6] = dRdm;
					dpar[7] = 1./R;
					dpar[8] = dpar[7]*dpar[7];

					cpar[0] = (cp*dcomplex(cost, sint) - 1.)/R; 
					cpar[1] = (cp*dcomplex(cost, -sint) - 1.)/R;

					cpar[2] = dcomplex(cost, sint)*(1. - cp*dRdm/R)/R + dRdm/(R*R);
					cpar[3] = dcomplex(cost, -sint)*(1. - cp*dRdm/R)/R + dRdm/(R*R);

					cpar[4] = 2.*dcomplex(cost, sint)*(-dRdm - cp + cp*dRdm*dRdm/R)/(R*R)
                          + 2.*(1. - dRdm*dRdm/R)/(R*R);
					cpar[5] = 2.*dcomplex(cost, -sint)*(-dRdm - cp + cp*dRdm*dRdm/R)/(R*R)
                          + 2.*(1. - dRdm*dRdm/R)/(R*R);

					cpar[6] = dpar[0]*cpar[0] + dpar[1]*cpar[1];
					cpar[7] = dpar[0]*cpar[2] + dpar[1]*cpar[3];
					cpar[8] = dpar[0]*cpar[4] + dpar[1]*cpar[5];

					if (norder<=4){

					for (iorder=0; iorder<=norder; iorder++){
						dZdmu[imu] += contribution_dZdmu[iorder](dpar,cpar);
						d2Zdmu2[imu] += contribution_d2Zdmu2[iorder](dpar,cpar);

						dZdm[imu] += contribution_dZdm[iorder](dpar,cpar);
						d2Zdm2[imu] += contribution_d2Zdm2[iorder](dpar,cpar);
					}

					}else{
					dZdm[imu] += contribution_dZdm[5](dpar,cpar);
					d2Zdm2[imu] += contribution_d2Zdm2[5](dpar,cpar);

					dZdmu[imu] += contribution_dZdmu[5](dpar,cpar);
					d2Zdmu2[imu] += contribution_d2Zdmu2[5](dpar,cpar);

					}

				} // t
			} // z
		} // y
	} // x

	dZdmu[imu] *= 2./double(nsite);
	d2Zdmu2[imu] *= 2./double(nsite);
	dZdm[imu] *= 2./double(nsite);
	d2Zdm2[imu] *= 2./double(nsite);
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

dcomplex dZdm_0th(double *dpar, dcomplex *cpar)
{
return dpar[6]*dpar[7];
}

dcomplex dZdm_1st(double *dpar, dcomplex *cpar)
{
return -cpar[7];
}

dcomplex dZdm_2nd(double *dpar, dcomplex *cpar)
{
return -(cpar[7]*cpar[6] + dpar[0]*dpar[1]*dpar[6]*dpar[8]);
}

dcomplex dZdm_3rd(double *dpar, dcomplex *cpar)
{
return -(cpar[7]*cpar[6]*cpar[6]
					+ dpar[0]*dpar[0]*dpar[1]*(cpar[0]*dpar[6]*dpar[8] - cpar[2]*dpar[7])
					+ dpar[0]*dpar[1]*dpar[1]*(cpar[1]*dpar[6]*dpar[8] - cpar[3]*dpar[7]) );
}

dcomplex dZdm_4rd(double *dpar, dcomplex *cpar)
{
return dcomplex(0.,0.);
}

dcomplex dZdm_exact(double *dpar, dcomplex *cpar)
{
return (dpar[6]*dpar[7] 
				+ (-cpar[7] - dpar[0]*dpar[1]*dpar[6]*dpar[8])
					/(1. - cpar[6] + dpar[0]*dpar[1]*dpar[7]) );
}

//--------------------------------------------------------------------------------------------------
dcomplex d2Zdm2_0th(double *dpar, dcomplex *cpar)
{
return ( 2. - dpar[6]*dpar[6]*dpar[7] )*dpar[7] ;
}

dcomplex d2Zdm2_1st(double *dpar, dcomplex *cpar)
{
return -cpar[8];
}

dcomplex d2Zdm2_2nd(double *dpar, dcomplex *cpar)
{
return -(  cpar[8]*cpar[6] 
				 + cpar[7]*cpar[7] 
				 + 2.*dpar[0]*dpar[1]*(1. - dpar[6]*dpar[6]*dpar[7])*dpar[8] );
}

dcomplex d2Zdm2_3rd(double *dpar, dcomplex *cpar)
{
return dcomplex(0.,0.);
}

dcomplex d2Zdm2_4rd(double *dpar, dcomplex *cpar)
{
return dcomplex(0.,0.);
}

dcomplex d2Zdm2_exact(double *dpar, dcomplex *cpar)
{
return ( 	2.*dpar[7] - dpar[6]*dpar[6]*dpar[8]
	- (cpar[7]+dpar[0]*dpar[1]*dpar[6]*dpar[8])*(cpar[7]+dpar[0]*dpar[1]*dpar[6]*dpar[8])
	/((1. - cpar[6] + dpar[0]*dpar[1]*dpar[7])*(1. - cpar[6] + dpar[0]*dpar[1]*dpar[7]))   );
}


//--------------------------------------------------------------------------------------------------
dcomplex dZdmu_0th(double *dpar, dcomplex *cpar)
{
return dcomplex(0.,0.);
}

dcomplex dZdmu_1st(double *dpar, dcomplex *cpar)
{
return -(dpar[2]*cpar[0]-dpar[3]*cpar[1]);
}

dcomplex dZdmu_2nd(double *dpar, dcomplex *cpar)
{
return -( (dpar[2]*cpar[0]-dpar[3]*cpar[1])*cpar[6] 
					+ (-dpar[2]*dpar[1]+dpar[3]*dpar[0])*dpar[7] );
}

dcomplex dZdmu_3rd(double *dpar, dcomplex *cpar)
{
return -( (dpar[2]*cpar[0]-dpar[3]*cpar[1])*cpar[6]*cpar[6] 
					+ (-2.*dpar[0]*dpar[1]*dpar[2] + dpar[0]*dpar[0]*dpar[3])*dpar[7]*cpar[0]
					+ ( 2.*dpar[0]*dpar[1]*dpar[3] - dpar[1]*dpar[1]*dpar[2])*dpar[7]*cpar[1] );
}

dcomplex dZdmu_4rd(double *dpar, dcomplex *cpar)
{
return -( (dpar[2]*cpar[0]-dpar[3]*cpar[1])*cpar[6]*cpar[6]*cpar[6]
					+ dpar[0]*dpar[1]*(dpar[1]*dpar[2] - dpar[0]*dpar[3])*dpar[8]
					+ cpar[6]*cpar[6]*dpar[7]*(dpar[0]*dpar[3]-dpar[1]*dpar[2])
					- 2.*dpar[0]*dpar[1]*dpar[7]*cpar[6]*(dpar[2]*cpar[0]-dpar[3]*cpar[1]) );
}

dcomplex dZdmu_exact(double *dpar, dcomplex *cpar)
{
return (- dpar[2]*cpar[0] + dpar[3]*cpar[1]
        + ( dpar[1]*dpar[2] + dpar[0]*dpar[3])*dpar[7] )
			/(1. - cpar[6] + dpar[0]*dpar[1]*dpar[7]);
}


//--------------------------------------------------------------------------------------------------
dcomplex d2Zdmu2_0th(double *dpar, dcomplex *cpar)
{
return dcomplex(0.,0.);
}

dcomplex d2Zdmu2_1st(double *dpar, dcomplex *cpar)
{
return -(dpar[2]*cpar[0]+dpar[3]*cpar[1]);
}

dcomplex d2Zdmu2_2nd(double *dpar, dcomplex *cpar)
{
return -(  ((dpar[2]*cpar[0]+dpar[3]*cpar[1]))*cpar[6]
				 + (dpar[2]*cpar[0]-dpar[3]*cpar[1])*(dpar[2]*cpar[0]-dpar[3]*cpar[1])
			   + (2. - dpar[2]*dpar[1] - dpar[3]*dpar[0])*dpar[7]  );
}

dcomplex d2Zdmu2_3rd(double *dpar, dcomplex *cpar)
{
return -(  ((dpar[2]*cpar[0]+dpar[3]*cpar[1]))*cpar[6]*cpar[6]
		+ 2.*(dpar[2]*cpar[0]-dpar[3]*cpar[1])*(dpar[2]*cpar[0]-dpar[3]*cpar[1])*cpar[6]
		+ cpar[0]*dpar[7]*(-2.*dpar[2]*dpar[2]*dpar[1] - 2.*dpar[0]*dpar[1]*dpar[2] 
											 + 4.*dpar[0] - dpar[0]*dpar[0]*dpar[3])
		+ cpar[1]*dpar[7]*(-2.*dpar[3]*dpar[3]*dpar[0] - 2.*dpar[0]*dpar[1]*dpar[3]
											 + 4.*dpar[1] - dpar[1]*dpar[1]*dpar[2])  );
}

dcomplex d2Zdmu2_4rd(double *dpar, dcomplex *cpar)
{
return -(  ((dpar[2]*cpar[0]+dpar[3]*cpar[1]))*cpar[6]*cpar[6]*cpar[6]  
		+ 3.*(dpar[2]*cpar[0]-dpar[3]*cpar[1])*(dpar[2]*cpar[0]-dpar[3]*cpar[1])*cpar[6]*cpar[6]
		+ dpar[8]*(dpar[1]*dpar[1]*dpar[2]*dpar[2] + dpar[0]*dpar[1]*dpar[1]*dpar[2] 
							- 4.*dpar[0]*dpar[1] + dpar[0]*dpar[0]*dpar[3]*dpar[3]
							+ dpar[0]*dpar[0]*dpar[1]*dpar[3])

		+ 4.*cpar[6]*(dpar[2]*cpar[0]-dpar[3]*cpar[1])*dpar[7]*
			(-dpar[1]*dpar[2] + dpar[0]*dpar[3])

		+ cpar[6]*cpar[6]*dpar[7]*(2. - dpar[1]*dpar[2] - dpar[0]*dpar[3])

		- 2.*dpar[0]*dpar[1]*dpar[7]*cpar[6]*(dpar[2]*cpar[0]+dpar[3]*cpar[1])

		- 2.*dpar[0]*dpar[1]*dpar[7]*(dpar[2]*cpar[0]-dpar[3]*cpar[1])
			*(dpar[2]*cpar[0]-dpar[3]*cpar[1])  );
}

dcomplex d2Zdmu2_exact(double *dpar, dcomplex *cpar)
{
return (  ( - dpar[2]*cpar[0] - dpar[3]*cpar[1] 
						+ (-2. + dpar[1]*dpar[2] + dpar[0]*dpar[3])*dpar[7] )
					/(1. - cpar[6] + dpar[0]*dpar[1]*dpar[7])
		- (- dpar[2]*cpar[0] + dpar[3]*cpar[1] + ( dpar[1]*dpar[2] + dpar[0]*dpar[3])*dpar[7])
			*(- dpar[2]*cpar[0] + dpar[3]*cpar[1] + ( dpar[1]*dpar[2] + dpar[0]*dpar[3])*dpar[7])
			/((1. - cpar[6] + dpar[0]*dpar[1]*dpar[7])*(1. - cpar[6] + dpar[0]*dpar[1]*dpar[7])) );
}

//--------------------------------------------------------------------------------------------------
 
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
			{"nmu", required_argument, 0, 'n'},
			{"O", required_argument, 0, 'o'},
			/* These options set a flag. */
			{0, 0, 0, 0}
			};
			
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "s:t:m:o:n:k:",
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

			case 'o':
				norder = atoi(optarg);
				cout << "order : " << norder << endl;
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

	/* Print any remaining command line arguments (not options). */
	/* Get the names of the files with the configurations */
	if ( optind >= argc ) 
	{
		cout << "Init.ERROR: output file is missing..." << endl;
		return 1;
	}
	else
	{ 
		outputname = argv[optind];
	}

	/* Set the sizes of the lattice*/
	nsite = Ns*Ns*Ns*Nt;

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
		outfile << mu << " " << std::real(dZdm[imu])  << " " << std::real(d2Zdm2[imu])
									<< " " << std::real(dZdmu[imu]) << " " << std::real(d2Zdmu2[imu])	<< endl;
	}
	outfile.close( );
}

//----------------------------------------------------------------------------------------------------
void allocateArrays( )
{
	dZdm = new dcomplex[nmu+1];
	d2Zdm2 = new dcomplex[nmu+1];
	dZdmu = new dcomplex[nmu+1];
	d2Zdmu2 = new dcomplex[nmu+1];
}

//----------------------------------------------------------------------------------------------------
void deallocateArrays( )
{
	delete[] dZdm;
	delete[] d2Zdm2;
	delete[] dZdmu;
	delete[] d2Zdmu2;
}

