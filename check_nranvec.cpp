#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <limits>
#include <vector>
#include <getopt.h>
#include <cstdio>
#include <algorithm>

using namespace std;

#define NTR 26
#define nn  49
#define cut 0

typedef complex<double> dcomplex;


static int readTraces( dcomplex *traces, int in, int *n )
{
  stringstream tmpname;	
	tmpname.str("");
	tmpname << "traces_" << n[in] << ".dat";
	string tracesname= tmpname.str();

	std::ifstream file;
	file.open( tracesname.c_str(), ios::in );

	if ( !file.is_open( ) ) return 1;
	int ii;
	file >> ii;
	for ( int iobs=0; iobs<NTR; iobs++ )
			 file >> traces[in*NTR + iobs];
	file.close( );

	return 0;
}

//-----------------------------------------------------------------------------------------------
int main( )
{
	// Read input parameters
  int nranvec[nn];
	dcomplex traces[NTR*nn];

	for (int it=0; it<nn; it++)
	{
		nranvec[it] = 100*(it+1);
		if( readTraces( traces, it, nranvec ) != 0)
		{
			cout << "ERROR: Error reading traces file " << it << endl;
			return 1;
		}
	}

	double remean[NTR], immean[NTR];
	for (int io=0; io<NTR; io++)
	{
		remean[io] = 0.;
		immean[io] = 0.;
		for (int it=cut; it<nn; it++)
		{
			remean[io] += std::real(traces[it*NTR+io]);
			immean[io] += std::imag(traces[it*NTR+io]);
		}
		remean[io] /= double(nn-cut);
		immean[io] /= double(nn-cut);
	}

	std::fstream file;
	file.open( "tracesout_re.dat", ios::out );
	for (int it=cut; it<nn; it++)
	{
		file << nranvec[it] << " " ;
		for (int io=0; io<NTR; io++ )
		{
			double ss = (std::real(traces[it*NTR+io]));//-remean[io]);
			//if (remean[io]!=0.) ss = ss/remean[io];
			file << ss << " ";
		}
		file << endl;
	}
	file.close();

	file.open( "tracesout_im.dat", ios::out );
	for (int it=cut; it<nn; it++)
	{
		file << nranvec[it] << " " ;
		for (int io=0; io<NTR; io++ )
		{
			double ss = (std::imag(traces[it*NTR+io]));//-immean[io]);
			//if (immean[io]!=0.) ss = ss/immean[io];
			file << ss << " ";
		}
		file << endl;
	}
	file.close();

	return 0;
}

