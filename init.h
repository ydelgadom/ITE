#ifndef INIT_H
#define INIT_H

string optionsname = "options";

char texthelp[]="Usage: exec [OPTION]... [FILE] [TYPE]\n"
		"Calculates the improved Taylor expansion on gauge configurations listed on FILE (or\n"
		"for a free configuration)\n"
		"\n"
		"TYPE: a          MILC data structure\n"
		"      b          binary data structure (old)\n"
		"      s          binary data structure v2 (new)\n"
		"      f          Fortran data structure\n"
		"      if empty:  text data structure\n"
		"\n"
		"Mandatory arguments to long options are mandatory for short options too.\n"
		"  -s, --Ns SSIZE             spatial lattice extent (default = 4)\n"
		"  -t, --Nt TSIZE             temporal lattice extent (default = 4)\n"
		"  -k, --kappa KAPPA          hopping parameter (inv. mass) kappa,\n"
		"  -m, --mass MASS            mass parameter,\n"
		"                               ignored for staggered (default = 0.1)\n"
		"  -f, --free                 use a free gauge configuration instead\n"
		"                               of 'CONFIGURATION'\n"
		"  -h  --help                 display this help and exit\n"
		"\n"
		"Exit status:\n"
		" 0  if OK,\n"
		" 1  if ERROR,\n"
		"\n"
		"Based on Hans-Peter's code :-)\n"
		"Report bugs to ydelgado83@gmail.com\n";


int init(int &argc, char *argv[])
{
	if ( myid==0 ) cout << "Init: Initializing... " << endl;

	if(argc<2)
	{
		if ( myid==0 ) cout << endl << texthelp << endl;
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
			{"kappa", required_argument, 0, 'k'},
			{"mass", required_argument, 0, 'm'},
			/* These options set a flag. */
			{"free", no_argument, 0, 'f'},
			{"help", no_argument, 0, 'h'},
			{0, 0, 0, 0}
		};
			
		/* getopt_long stores the option index here. */
		int option_index = 0;

		c = getopt_long (argc, argv, "s:t:k:m:fh",
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
				break;

			case 't':
				Nt = atoi(optarg);
				break;

			case 'k':
				kappa = atof(optarg);
				mass = 1./(2.*kappa) - 4.;
				break;

			case 'm':
				mass = atof(optarg);
				kappa = 1./(2.*mass+4);
				break;

			case 'f':
				freeconf = true;
				break;

			case 'h':
				if (myid == 0) cout << endl << texthelp << endl;
				abort();

			default:
				if (myid == 0) cout << endl << texthelp << endl;
				abort();
		}
	}

	/* Print any remaining command line arguments (not options). */
	/* Get the names of the files with the configurations */
	string tmpstring;
	if ( optind >= argc ) 
	{
		if ( myid==0 ) 
			cout << "Init.ERROR: file with paths is missing..." << endl;
		return 1;
	}
	else tmpstring = argv[optind];

	stringstream tmpname;	
	std::ifstream file( tmpstring.c_str(), ios::in );

	if ( !file.is_open( ) )
	{
		if ( myid==0 ) 
			cout << "Init.ERROR: file with output/input files couldn't be opened..." << endl;
		return 1;
	}

	nconf = std::count(std::istreambuf_iterator<char>(file), 
          std::istreambuf_iterator<char>(), '\n')-2; 
					// -2 to substract the paths to in(out)put files

	if ( nconf==0 && freeconf==false )
	{
		if ( myid==0 ){
			cout << "Init.Warning: No configuration's names specified, switching to free configuration!" << endl;}
		freeconf = true;
	}

	// get input and output paths
	file.seekg(0, std::ios::beg);
	tmpstring = "";
	std::getline(file, tmpstring);	// path to input
	std::getline(file, outpath);	// path to output

	tmpname.str("");
	tmpname << outpath << "traces" << ".dat";
	trname= tmpname.str();

	if( freeconf == false )
	{
		if ( nconf%nproc!=0 && myid==0)
			cout << "Init.Warning: nconf/nproc is not exact!!!!" << endl;

		mynconf = nconf/nproc;

		configname = new string[nconf];
		logname = new string[nconf];

		for ( int iconf=0; iconf<nconf; iconf++ )
		{
			std::getline(file, configname[iconf]);

			tmpname.str("");
			tmpname << outpath << configname[iconf] << ".log";
			logname[iconf] = tmpname.str();

			tmpname.str("");
			tmpname << tmpstring << configname[iconf];
			configname[iconf] = tmpname.str( );

			if ( myid==0 ) cout << "config " << iconf << ":" << configname[iconf] << endl;
			//if ( myid==0 ) cout << "logfile " << iconf << " :" << logname[iconf] << endl; 
		}
	}
	else
	{
		nconf = nproc;
		mynconf = 1;

		configname = new string[nconf];
		logname = new string[nconf];

		for ( int iconf=0; iconf<nconf; iconf++ )
		{
			tmpname.str("");
			tmpname << "free_configuration_" << iconf;	
			configname[iconf] = tmpname.str();

			tmpname.str("");
			tmpname << outpath << configname[iconf] << ".log";
			logname[iconf] = tmpname.str();

			tmpname.str("");
			tmpname << tmpstring << configname[iconf];
			configname[iconf] = tmpname.str( );

			if ( myid==0 ) cout << "config " << iconf << ":" << configname[iconf] << endl;
			//if ( myid==0 ) cout << "logfile " << iconf << " :" << logname[iconf] << endl; 
		}
	}
	file.close( );

	/* Set the sizes of the lattice and of matrices*/
	Ns3 = Ns*Ns*Ns;
	nsite = Ns3*Nt;

	if( freeconf == false )
	{
		if( argc>optind+1 ) 
			config_type = argv[optind+1][0];
		else
			config_type = 't'; // if empty --> text data structure
	}
	else
	{
		config_type = 't';
	}

	if ( myid==0 ) cout << "Init: done!\n" << endl;
	return 0;
}

#endif // INIT_H
