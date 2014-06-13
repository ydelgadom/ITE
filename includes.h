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

#include "mpi.h"

extern "C" 
{
#include "ranlxd.h"
}

#include "include/ConfigData.hpp"
#include "variables.h"
#include "init.h"
#include "matvec_operations.h"
