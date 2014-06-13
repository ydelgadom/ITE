#define DIM 4 // Lattice dimension
#define MD 3 // SU(MD) Matrix dimension

#define NTR 37 // number of traces
#define SEED 6985
#define NRANVEC 500  // number of vectors to compute the traces

#define MM2 // second version of the MTE
//#define MM1 // first version of the MTE

#define ONE complex<double>( 1.0, 0.0 )
#define mONE complex<double>( -1.0, 0.0 )
#define ZERO complex<double>( 0.0, 0.0 )
#define II complex<double>( 0.0, 1.0 )

typedef vector<complex<double> > cvector;
typedef complex<double> dcomplex;

// Global variables
int Ns=4, Nt=4;
int Ns3 = Ns*Ns*Ns;
int nsite = Ns3*Nt;
int nconf;
bool freeconf;
double kappa=0.1, mass;

string outpath, trname, *logname;
string *configname;
ofstream logfile;
char config_type;

// MPI variables
int nproc, myid, mynconf; 

ConfigData *config;

// Functions in main.cpp
void write_traces( dcomplex *tr );
void traces( dcomplex *tr );
int readConfig( int iconf );
void test_nranvec( dcomplex *tr );
