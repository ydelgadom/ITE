#include "includes.h"

int main(int argc, char *argv[])
{
	// MPI
  MPI::Init ( argc, argv );
  nproc = MPI::COMM_WORLD.Get_size( );
  myid = MPI::COMM_WORLD.Get_rank( );

	// Initialize random number generator
	rlxd_init( 2, SEED );

	// Set options
	if( init(argc, argv) != 0)
	{
		if ( myid==0 ) cout << "ERROR: Error initializing!" << endl;
		MPI::COMM_WORLD.Abort(2);
	}

	// Allocate lattice vector
	config = new ConfigData(Ns, Ns, Ns, Nt, MD);
	dcomplex *tr = new dcomplex[NTR*mynconf];

	// Loop over all configuration files
  for (int i=0; i<mynconf; i++ )
	{
		int iconf = i + myid*mynconf;

		logfile.open( logname[iconf].c_str(), ios::out );

		// Print settings 
		logfile << endl << "Settings:" << endl;
		logfile << "Lattice data = " << configname[iconf] << endl;
		logfile << "Lattice size = " << Ns << "x" << Ns << "x" << Ns << "x" << Nt << endl;
		logfile << "Job with = " << nproc << " jobs and " << mynconf << " confs/node" << endl;
		logfile << "Hopping parameter kappa = " << kappa << endl;
#ifdef MKL
		logfile << "Compiled with Intel MKL multithreading support" << endl;
#endif

		if ( readConfig( iconf ) != 0 ) MPI::COMM_WORLD.Abort(2);

		// Getting some other informations
		dcomplex poll=config->calcPoll( );
		dcomplex plaq=config->calcPlaq( );

		// Switch to antiperiodic bcs
		config->antiperbc( );
		//config->dumpConfig( );// print it
		//config->runTests( );
		
		// Compute traces
		//test_nranvec( &tr[i*NTR] );
		traces( &tr[i*NTR] );

		// Output
		logfile << endl << "Results:" << endl;
		logfile << "Lattice data = " << configname[iconf] << endl;
		logfile << "Polyakov loop = " << poll << endl;
		logfile << "Plaquette = " << plaq << endl;
		logfile << endl << "Calculation finished!" << endl;

		logfile.close( );
	}

	// Write to disk the traces (needed to compute observables at finite mu)
	write_traces( tr );

	// Memory deallocations
	delete[] configname;
	delete[] logname;
	delete config;	config = 0;
	delete[] tr;

  //  Terminate MPI.
  MPI::Finalize( );
}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
void write_traces( dcomplex *tr )
{
	dcomplex *allTr = new dcomplex[NTR*nconf];

	int nelements = NTR*mynconf;

	if ( myid==0 ) 
		cout << "\nWriting traces..." << endl;

	// Gather results from all configurations
	MPI::COMM_WORLD.Gather( tr, nelements, MPI::DOUBLE_COMPLEX,
													allTr, nelements, MPI::DOUBLE_COMPLEX, 0);

	if ( myid==0 )
	{
		ofstream outfile;
		outfile.open( trname.c_str(), ios::out );

		// Write results to file traces.dat
		for (int iconf=0; iconf<nconf; iconf++ )
		{
			outfile << iconf ;
			for (	int itr=0; itr<NTR; itr++ )
			{
				outfile << " " << allTr[iconf*NTR + itr];
			}
			outfile << endl;
		}
		outfile.close( );
	}

	delete[] allTr;

	if ( myid==0 ) 
		cout << "Write_traces: Done" << endl;
}

//-------------------------------------------------------------------------------
void traces( dcomplex *tr )
{
	cvector v0(12*nsite), v1(12*nsite);
	cvector v2(12*nsite), v3(12*nsite), v4(12*nsite);

	// TRACES
	dcomplex tD=ZERO, tDD=ZERO, tMDDD=ZERO, tMDMDDD=ZERO, tMbDDD=ZERO;
	dcomplex tMbDMbDDD=ZERO, tMbD=ZERO, tMbDD=ZERO, tMbDMbD=ZERO, tMbDMbDD=ZERO; 
	dcomplex tMDMbDMbD=ZERO, tMbDMbDMbD=ZERO, tMDMbDMbDMbD=ZERO, tMbDMbDMbDMbD=ZERO;
  dcomplex tMD=ZERO, tMDD=ZERO, tMDMbD=ZERO, tMDMbDD=ZERO, tMDMbDDD=ZERO, tMDMDMbD=ZERO;
	dcomplex tMDMDMDMD=ZERO, tMDMDMDMbD=ZERO, tMDMD=ZERO, tMDMDD=ZERO, tMDMDMD=ZERO;
	dcomplex tMDMDMbDMbD=ZERO, tMDMbDMDMbD=ZERO;
	dcomplex tMbDMDD=ZERO, tMbDMDDD=ZERO;
	dcomplex tMDMDMDD=ZERO, tMDMDMbDD=ZERO, tMDMbDMbDD=ZERO, tMbDMbDMbDD=ZERO;
	dcomplex tMDMbDMDD=ZERO, tMbDMDMDD=ZERO, tMbDMDMbDD=ZERO, tMbDMbDMDD=ZERO;

	// Compute the trace of a matrix using random vectors
	for ( int ivec=0; ivec<NRANVEC; ivec++ )
	{
		ranvec( &v0 );

	  if ( vec_invDwvec( &v0, &v2 ) ){
			ivec--;
			cout << "000 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tD += dot( v0, v2 );	//	0

	  if ( vec_invDwvec( &v2, &v1 ) ){
			ivec--;
			cout << "0 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tDD += dot( v0, v1 );	//	8

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "8 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		vec_Mvec_Mbvec( &v1, &v3, &v1 );
		tMDDD += dot( v0, v3 );  	//	16
		tMbDDD += dot( v0, v1 );	//	17

	  if ( vec_invDwvec( &v3, &v3 ) ){
			ivec--;
			cout << "16 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		vec_Mvec_Mbvec( &v3, &v3, &v4 );
		tMDMDDD += dot( v0, v3 );	//	23
		tMbDMDDD += dot( v0, v4 ); //	28

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "17 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		vec_Mvec_Mbvec( &v1, &v3, &v1 );
		tMbDMbDDD += dot( v0, v1 );  // 25
		tMDMbDDD += dot( v0, v3 );   // 24

		// observables with Mbar
		vec_Mvec_Mbvec( &v2, &v2, &v1 );
		tMD += dot( v0, v2 );		//	1
		tMbD += dot( v0, v1 );	//	2

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "2 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDD += dot( v0, v1 );	//	7

		vec_Mvec_Mbvec( &v1, &v3, &v1 );
		tMbDMbD += dot( v0, v1 );	//	4	

		if ( vec_invDwvec( &v3, &v3 ) ){
			ivec--;
			cout << "27 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMbDD += dot( v0, v3 );	//	27

		vec_Mvec_Mbvec( &v3, &v3, &v4 );
		if ( vec_invDwvec( &v3, &v3 ) ){
			ivec--;
			cout << "31 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMDMbDD += dot( v0, v3 ); // 31

		if ( vec_invDwvec( &v4, &v4 ) ){
			ivec--;
			cout << "35 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMDMbDD += dot( v0, v4 ); // 35	

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "4 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMbDD += dot( v0, v1 );	//	15

		vec_Mvec_Mbvec( &v1, &v3, &v1 );
		tMDMbDMbD += dot( v0, v3 );	  //	11
		tMbDMbDMbD += dot( v0, v1 );	//	12

		if ( vec_invDwvec( &v3, &v3 ) ){
			ivec--;
			cout << "36 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMbDMbDD += dot( v0, v3 );	// 36

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "12 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMbDMbDD += dot( v0, v1 );	// 30

		vec_Mvec_Mbvec( &v1, &v3, &v1 );
		tMDMbDMbDMbD += dot( v0, v3 );	//	21
		tMbDMbDMbDMbD += dot( v0, v1 );	//	22*/

		// observables with M
		v1 = v2;
	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "1 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDD += dot( v0, v1 );	//	6

		vec_Mvec_Mbvec( &v1, &v1, &v2 );
		tMDMD += dot( v0, v1 );	  //	3
		tMDMbD += dot( v0, v2 );	//	5

	  if ( vec_invDwvec( &v2, &v2 ) ){
			ivec--;
			cout << "5 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMDD += dot( v0, v2 );	//	14

		vec_Mvec_Mbvec( &v2, &v2, &v3 );
		if ( vec_invDwvec( &v2, &v2 ) ){
			ivec--;
			cout << "26 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMbDMDD += dot( v0, v2 ); // 32

		if ( vec_invDwvec( &v3, &v3 ) ){
			ivec--;
			cout << "34 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMbDMDD += dot ( v0, v3 ); // 34

		vec_Mbvec( &v2, &v2 );
		tMDMbDMDMbD += dot( v0, v2 ); // 26

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "13 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMDD += dotu( v0, v1 );	//	13

		vec_Mvec_Mbvec( &v1, &v1, &v2 );
		tMDMDMD += dotu( v0, v1 );	//	9
		tMDMDMbD += dot( v0, v2 );	//	10

	  if ( vec_invDwvec( &v2, &v2 ) ){
			ivec--;
			cout << "10 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMDMDD += dot( v0, v2 ); // 33

		vec_Mvec_Mbvec( &v2, &v2, &v3 );
		tMDMDMDMbD += dot( v0, v2 );	//	19
		tMDMDMbDMbD += dot( v0, v3 );	//	20

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "9 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMDMDD += dotu( v0, v1 );	// 29

		vec_Mvec( &v1, &v1 );
		tMDMDMDMD += dotu( v0, v1 );	//	18*/

		if ( ivec%10==0) cout << "Traces " << ivec << endl;
	}

	tr[0] = tD/double(NRANVEC);
	tr[1] = tMD/double(NRANVEC);
	tr[2] = tMbD/double(NRANVEC);

	tr[3] = tMDMD/double(NRANVEC);
	tr[4] = tMbDMbD/double(NRANVEC);
	tr[5] = tMDMbD/double(NRANVEC);
	tr[6] = tMDD/double(NRANVEC);
	tr[7] = tMbDD/double(NRANVEC);
	tr[8] = tDD/double(NRANVEC);

	tr[9]  = tMDMDMD/double(NRANVEC);	
	tr[10] = tMDMDMbD/double(NRANVEC);
	tr[11] = tMDMbDMbD/double(NRANVEC);
	tr[12] = tMbDMbDMbD/double(NRANVEC);
	tr[13] = tMDMDD/double(NRANVEC);
	tr[14] = tMbDMDD/double(NRANVEC);
	tr[15] = tMbDMbDD/double(NRANVEC);
	tr[16] = tMDDD/double(NRANVEC);
	tr[17] = tMbDDD/double(NRANVEC);

  tr[18] = tMDMDMDMD/double(NRANVEC); 
	tr[19] = tMDMDMDMbD/double(NRANVEC); 
	tr[20] = tMDMDMbDMbD/double(NRANVEC);
	tr[21] = tMDMbDMbDMbD/double(NRANVEC);
	tr[22] = tMbDMbDMbDMbD/double(NRANVEC);
	tr[23] = tMDMDDD/double(NRANVEC);
	tr[24] = tMDMbDDD/double(NRANVEC);
	tr[25] = tMbDMbDDD/double(NRANVEC);
	tr[26] = tMDMbDMDMbD/double(NRANVEC);
	tr[27] = tMDMbDD/double(NRANVEC);
	tr[28] = tMbDMDDD/double(NRANVEC);

	tr[29] = tMDMDMDD/double(NRANVEC);
	tr[30] = tMbDMbDMbDD/double(NRANVEC);
	tr[31] = tMDMDMbDD/double(NRANVEC);
	tr[32] = tMDMbDMDD/double(NRANVEC);
	tr[33] = tMbDMDMDD/double(NRANVEC);
	tr[34] = tMbDMbDMDD/double(NRANVEC);
	tr[35] = tMbDMDMbDD/double(NRANVEC);
	tr[36] = tMDMbDMbDD/double(NRANVEC);	
}

//---------------------------------------------------------------------------------
int readConfig( int iconf )
{
	/* Preparing the configuration lattice (reading and tests) */
	logfile << "readConf: Preparing the configuration lattice... " << endl;

	int ret;

	if( freeconf == false ) 
	{
		switch( config_type )
		{
		  case 'b':
				ret=config->readBinaryConfig( configname[iconf] );
				if(ret!=0){
					logfile << "readConf: ERROR: Reading binary config file:" 
									<< configname[iconf] << endl;
					return 1;
				}
				break;
			case 's':
				ret=config->readBinaryConfig2( configname[iconf] );
				if(ret!=0){
					logfile << "readConf: ERROR: Reading binary config file (new storage format):" 
									<< configname[iconf] << endl;
					return 1;
				}
				break;
			case 'a':
				ret=config->MILCreadConfig( configname[iconf] );
				if(ret!=0){
					logfile << "readConf: ERROR: Reading MILC config file:" 
									<< configname[iconf] << endl;
					return 1;
				}
				break;
			case 'f':
				ret=config->readFConfig( configname[iconf] );
				if(ret!=0){
					logfile << "readConf: ERROR: Reading Fortran config file:" 
									<< configname[iconf] << endl;
					return 1;
				}
				break;
			case 't':
				ret=config->readConfig( configname[iconf] );
				if(ret!=0){
					logfile << "readConf: ERROR: Reading config file:" 
									<< configname[iconf] << endl;
					return 1;
				}
				break;
			default:
				logfile << "readConf: ERROR: Requested TYPE not understood!" << endl;
				return 1;
		}
	}
	else
	{
		config->freeConfig( );
	}

	return 0;
}

//-------------------------------------------------------------------------------
void test_nranvec( dcomplex *tr )
{
	cvector v0(12*nsite), v1(12*nsite);
	cvector v2(12*nsite), v3(12*nsite), v4(12*nsite);

	// TRACES
	dcomplex tD=ZERO, tDD=ZERO, tMDDD=ZERO, tMDMDDD=ZERO, tMbDDD=ZERO;
	dcomplex tMbDMbDDD=ZERO, tMbD=ZERO, tMbDD=ZERO, tMbDMbD=ZERO, tMbDMbDD=ZERO; 
	dcomplex tMDMbDMbD=ZERO, tMbDMbDMbD=ZERO, tMDMbDMbDMbD=ZERO, tMbDMbDMbDMbD=ZERO;
  dcomplex tMD=ZERO, tMDD=ZERO, tMDMbD=ZERO, tMDMbDD=ZERO, tMDMbDDD=ZERO, tMDMDMbD=ZERO;
	dcomplex tMDMDMDMD=ZERO, tMDMDMDMbD=ZERO, tMDMD=ZERO, tMDMDD=ZERO, tMDMDMD=ZERO;
	dcomplex tMDMDMbDMbD=ZERO, tMDMbDMDMbD=ZERO;
	dcomplex tMbDMDD=ZERO, tMbDMDDD=ZERO;
	dcomplex tMDMDMDD=ZERO, tMDMDMbDD=ZERO, tMDMbDMbDD=ZERO, tMbDMbDMbDD=ZERO;
	dcomplex tMDMbDMDD=ZERO, tMbDMDMDD=ZERO, tMbDMDMbDD=ZERO, tMbDMbDMDD=ZERO;

	// Compute the trace of a matrix using random vectors
	for ( int ivec=0; ivec<NRANVEC; ivec++ )
	{
		ranvec( &v0 );

	  if ( vec_invDwvec( &v0, &v2 ) ){
			ivec--;
			cout << "000 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tD += dot( v0, v2 );	//	0

	  if ( vec_invDwvec( &v2, &v1 ) ){
			ivec--;
			cout << "0 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tDD += dot( v0, v1 );	//	8

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "8 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		vec_Mvec_Mbvec( &v1, &v3, &v1 );
		tMDDD += dot( v0, v3 );  	//	16
		tMbDDD += dot( v0, v1 );	//	17

	  if ( vec_invDwvec( &v3, &v3 ) ){
			ivec--;
			cout << "16 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		vec_Mvec_Mbvec( &v3, &v3, &v4 );
		tMDMDDD += dot( v0, v3 );	//	23
		tMbDMDDD += dot( v0, v4 ); //	28

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "17 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		vec_Mvec_Mbvec( &v1, &v3, &v1 );
		tMbDMbDDD += dot( v0, v1 );  // 25
		tMDMbDDD += dot( v0, v3 );   // 24

		// observables with Mbar
		vec_Mvec_Mbvec( &v2, &v2, &v1 );
		tMD += dot( v0, v2 );		//	1
		tMbD += dot( v0, v1 );	//	2

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "2 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDD += dot( v0, v1 );	//	7

		vec_Mvec_Mbvec( &v1, &v3, &v1 );
		tMbDMbD += dot( v0, v1 );	//	4	

		if ( vec_invDwvec( &v3, &v3 ) ){
			ivec--;
			cout << "27 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMbDD += dot( v0, v3 );	//	27

		vec_Mvec_Mbvec( &v3, &v3, &v4 );
		if ( vec_invDwvec( &v3, &v3 ) ){
			ivec--;
			cout << "31 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMDMbDD += dot( v0, v3 ); // 31

		if ( vec_invDwvec( &v4, &v4 ) ){
			ivec--;
			cout << "35 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMDMbDD += dot( v0, v4 ); // 35	

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "4 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMbDD += dot( v0, v1 );	//	15

		vec_Mvec_Mbvec( &v1, &v3, &v1 );
		tMDMbDMbD += dot( v0, v3 );	  //	11
		tMbDMbDMbD += dot( v0, v1 );	//	12

		if ( vec_invDwvec( &v3, &v3 ) ){
			ivec--;
			cout << "36 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMbDMbDD += dot( v0, v3 );	// 36

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "12 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMbDMbDD += dot( v0, v1 );	// 30

		vec_Mvec_Mbvec( &v1, &v3, &v1 );
		tMDMbDMbDMbD += dot( v0, v3 );	//	21
		tMbDMbDMbDMbD += dot( v0, v1 );	//	22*/

		// observables with M
		v1 = v2;
	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "1 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDD += dot( v0, v1 );	//	6

		vec_Mvec_Mbvec( &v1, &v1, &v2 );
		tMDMD += dot( v0, v1 );	  //	3
		tMDMbD += dot( v0, v2 );	//	5

	  if ( vec_invDwvec( &v2, &v2 ) ){
			ivec--;
			cout << "5 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMDD += dot( v0, v2 );	//	14

		vec_Mvec_Mbvec( &v2, &v2, &v3 );
		if ( vec_invDwvec( &v2, &v2 ) ){
			ivec--;
			cout << "26 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMbDMDD += dot( v0, v2 ); // 32

		if ( vec_invDwvec( &v3, &v3 ) ){
			ivec--;
			cout << "34 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMbDMDD += dot ( v0, v3 ); // 34

		vec_Mbvec( &v2, &v2 );
		tMDMbDMDMbD += dot( v0, v2 ); // 26

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "13 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMDD += dotu( v0, v1 );	//	13

		vec_Mvec_Mbvec( &v1, &v1, &v2 );
		tMDMDMD += dotu( v0, v1 );	//	9
		tMDMDMbD += dot( v0, v2 );	//	10

	  if ( vec_invDwvec( &v2, &v2 ) ){
			ivec--;
			cout << "10 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMbDMDMDD += dot( v0, v2 ); // 33

		vec_Mvec_Mbvec( &v2, &v2, &v3 );
		tMDMDMDMbD += dot( v0, v2 );	//	19
		tMDMDMbDMbD += dot( v0, v3 );	//	20

	  if ( vec_invDwvec( &v1, &v1 ) ){
			ivec--;
			cout << "9 CG ERROR in iranvec " << ivec << endl;
			continue;
		}
		tMDMDMDD += dotu( v0, v1 );	// 29

		vec_Mvec( &v1, &v1 );
		tMDMDMDMD += dotu( v0, v1 );	//	18*/

		if (ivec%100==0)
		{
			tr[0] = tD/double(ivec+1);
			tr[1] = tMD/double(ivec+1);
			tr[2] = tMbD/double(ivec+1);

			tr[3] = tMDMD/double(ivec+1);
			tr[4] = tMbDMbD/double(ivec+1);
			tr[5] = tMDMbD/double(ivec+1);
			tr[6] = tMDD/double(ivec+1);
			tr[7] = tMbDD/double(ivec+1);
			tr[8] = tDD/double(ivec+1);

			tr[9]  = tMDMDMD/double(ivec+1);	
			tr[10] = tMDMDMbD/double(ivec+1);
			tr[11] = tMDMbDMbD/double(ivec+1);
			tr[12] = tMbDMbDMbD/double(ivec+1);
			tr[13] = tMDMDD/double(ivec+1);
			tr[14] = tMDMbDD/double(ivec+1);
			tr[15] = tMbDMbDD/double(ivec+1);
			tr[16] = tMDDD/double(ivec+1);
			tr[17] = tMbDDD/double(ivec+1);

			tr[18] = tMDMDMDMD/double(ivec+1); 
			tr[19] = tMDMDMDMbD/double(ivec+1); 
			tr[20] = tMDMDMbDMbD/double(ivec+1);
			tr[21] = tMDMbDMbDMbD/double(ivec+1);
			tr[22] = tMbDMbDMbDMbD/double(ivec+1);
			tr[23] = tMDMDDD/double(ivec+1);
			tr[24] = tMDMbDDD/double(ivec+1);
			tr[25] = tMbDMbDDD/double(ivec+1);
			tr[26] = tMDMbDMDMbD/double(ivec+1);
			tr[27] = tMDMbDD/double(ivec+1);
			tr[28] = tMbDMDDD/double(ivec+1);

			tr[29] = tMDMDMDD/double(ivec+1);
			tr[30] = tMbDMbDMbDD/double(ivec+1);
			tr[31] = tMDMDMbDD/double(ivec+1);
			tr[32] = tMDMbDMDD/double(ivec+1);
			tr[33] = tMbDMDMDD/double(ivec+1);
			tr[34] = tMbDMbDMDD/double(ivec+1);
			tr[35] = tMbDMDMbDD/double(ivec+1);
			tr[36] = tMDMbDMbDD/double(ivec+1);

			char nname[200];
			sprintf( nname, "traces_%d.dat", ivec ); 
			cout << nname << endl;
			ofstream outfile;
			outfile.open( nname, ios::out );

			// Write results to file traces.dat
			int iconf = 0;
			outfile << iconf ;
			for (	int itr=0; itr<NTR; itr++ )
			{
				outfile << " " << tr[itr];
			}
			outfile << endl;
			outfile.close( );

		}
	}	
}
